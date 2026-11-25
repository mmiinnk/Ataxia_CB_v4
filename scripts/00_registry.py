"""Export genotype–disease rows from the curated G3 workbook, preserving the human review fields."""
import argparse
import csv
import os
from pathlib import Path
from tempfile import NamedTemporaryFile

from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
FIELDS = ("row_id", "row_id_h", "hgnc_symbol", "disease_abbrev_h", "disease_name_h",
          "inheritance_h", "mutation_class_h", "mech_class_h", "mech_class_alt_h",
          "primary_pathology_h", "evidence_mech_h", "evidence_path_h", "source_urls_h",
          "review_flag_h", "core_grade_pending_h", "g3_verdict")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, default=ROOT / "inputs/merged.xlsx")
    parser.add_argument("--out", type=Path, default=ROOT / "outputs/disease_registry.tsv")
    args = parser.parse_args()
    workbook = load_workbook(args.workbook, read_only=True, data_only=True)
    try:
        values = iter(workbook["main"].values)
        header = next(values)
        positions = [header.index(field) for field in FIELDS]
        rows = [[row[i] for i in positions] for row in values if any(v is not None for v in row)]
    finally:
        workbook.close()
    ids = [row[FIELDS.index("row_id_h")] for row in rows]
    if not rows or not all(ids) or len(set(ids)) != len(ids):
        raise ValueError("The main sheet must contain a unique nonempty row_id_h for each disease record")
    studies = {key: [] for key in ids}
    with (ROOT / "config/cohorts.tsv").open(encoding="utf-8") as handle:
        for cohort in csv.DictReader(handle, delimiter="\t"):
            for key in filter(None, cohort["row_ids"].split(";")):
                if key not in studies:
                    raise ValueError(f"Unmapped disease record in cohort registry: {key}")
                studies[key].append(cohort["study_id"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=args.out.parent, delete=False) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([*FIELDS, "registered_studies", "review_status"])
        for key, row in zip(ids, rows):
            linked = studies[key]
            writer.writerow([*row, ";".join(linked), "有登记资源，待逐队列质控" if linked else "待逐型核查"])
    os.replace(handle.name, args.out)
    print(f"{len(rows)} disease records, {len({row[2] for row in rows})} distinct genes")


if __name__ == "__main__":
    main()
