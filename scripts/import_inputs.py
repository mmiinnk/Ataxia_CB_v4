"""Import the four reference matrices and the curated workbook, or verify the local copies."""
import argparse
import csv
import hashlib
import os
from pathlib import Path
import shutil
from tempfile import NamedTemporaryFile

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "inputs" / "sha256.tsv"


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path)
    parser.add_argument("--workbook", type=Path)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--link", action="store_true",
                        help="reference sources read-only via symlinks instead of copying bytes")
    args = parser.parse_args()
    with (ROOT / "config/reference.tsv").open() as handle:
        names = [row["filename"] for row in csv.DictReader(handle, delimiter="\t")]
    destinations = [ROOT / "inputs/reference" / name for name in names]
    destinations.append(ROOT / "inputs" / "merged.xlsx")
    if args.verify:
        if args.reference_dir or args.workbook:
            parser.error("--verify takes no input paths")
        with MANIFEST.open() as handle:
            hashes = {row["file"]: row["sha256"] for row in csv.DictReader(handle, delimiter="\t")}
        for dest in destinations:
            key = str(dest.relative_to(ROOT))
            if sha256(dest) != hashes[key]:
                raise SystemExit(f"checksum mismatch: {key}")
            print(f"OK  {key}", flush=True)
        return
    if not args.reference_dir or not args.workbook:
        parser.error("provide --reference-dir and --workbook, or use --verify")
    sources = [args.reference_dir / name for name in names] + [args.workbook]
    for source in sources:
        if not source.is_file():
            parser.error(f"missing source: {source}")
    rows = []
    for source, dest in zip(sources, destinations):
        dest.parent.mkdir(exist_ok=True, parents=True)
        digest = sha256(source)
        if dest.is_symlink():
            if Path(os.readlink(dest)).resolve() != source.resolve():
                raise SystemExit(f"existing symlink points elsewhere: {dest}")
        elif dest.exists():
            if sha256(dest) != digest:
                raise SystemExit(f"existing copy differs from source: {dest}")
        elif args.link:
            dest.symlink_to(source)
        else:
            with NamedTemporaryFile(dir=dest.parent, suffix=".partial", delete=False) as tmp:
                partial = Path(tmp.name)
            try:
                shutil.copyfile(source, partial)
                if sha256(partial) != digest:
                    raise SystemExit(f"copy checksum mismatch: {dest}")
                os.replace(partial, dest)
            finally:
                partial.unlink(missing_ok=True)
        rows.append([str(dest.relative_to(ROOT)), dest.stat().st_size, digest])
        print(f"OK  {dest.name}", flush=True)
    with NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=MANIFEST.parent, delete=False) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["file", "bytes", "sha256"])
        writer.writerows(rows)
    os.replace(handle.name, MANIFEST)


if __name__ == "__main__":
    main()
