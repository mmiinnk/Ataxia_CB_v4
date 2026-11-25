"""Census of all donor/stage strata and optional whole-transcriptome PC count sums."""
import argparse
import csv
import os
from pathlib import Path
from tempfile import TemporaryDirectory

from anndata.io import read_elem, sparse_dataset
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
KEYS = ["donor", "stage", "disease"]


def read_obs(file, spec):
    obs = file["obs"]
    columns = {"donor": spec["donor_column"], "stage": spec["stage_column"],
               "disease": spec["disease_column"], "cell_type": spec["cell_column"]}
    frame = pd.DataFrame({key: read_elem(obs[col]) for key, col in columns.items()})
    if frame[["donor", "stage", "disease"]].isna().any().any():
        raise ValueError(f"Missing stratification metadata in {spec['cohort']}")
    return frame.astype("string")


def census(frame, spec):
    frame["is_pc"] = frame.cell_type.eq(spec["pc_label"]).fillna(False)
    groups = frame.groupby(KEYS, sort=True, observed=True)
    strata = groups.agg(n_nuclei=("is_pc", "size"), n_pc=("is_pc", "sum")).reset_index()
    strata["reference_included"] = ((spec["include_disease"] == "*")
                                     | strata.disease.eq(spec["include_disease"]))
    strata["pooled"] = strata.donor.str.lower().eq("pooled")
    strata["at_min_pc"] = strata.n_pc.ge(int(spec["min_pc"]))
    strata["eligible_individual"] = strata.reference_included & strata.at_min_pc & ~strata.pooled
    strata.insert(0, "cohort", spec["cohort"])
    kept = strata.loc[strata.reference_included]
    summary = {"cohort": spec["cohort"], "n_nuclei_all": len(frame), "n_pc_all": int(frame.is_pc.sum()),
               "n_pc_reference": int(kept.n_pc.sum()), "n_cell_type_missing": int(frame.cell_type.isna().sum()),
               "n_individuals_reference": kept.loc[~kept.pooled, "donor"].nunique(),
               "n_individuals_with_eligible_stratum": strata.loc[strata.eligible_individual, "donor"].nunique()}
    if summary["n_pc_reference"] == 0:
        raise ValueError(f"No reference PCs for {spec['cohort']}. Check cell and disease labels.")
    return strata, summary


def pseudobulk(file, frame, spec, out):
    # Each row is a biological stratum. Keep thin strata visible for downstream QC.
    keep = frame.is_pc.copy()
    if spec["include_disease"] != "*":
        keep &= frame.disease.eq(spec["include_disease"])
    cells = frame.loc[keep, KEYS].copy()
    detail = []
    for name in ("dev_state", "subtype"):
        if name in file["obs"]:
            cells[name] = pd.array(read_elem(file["obs"][name]), dtype="string")[keep.to_numpy(dtype=bool)].fillna("unannotated")
            detail.append(name)
    keys = KEYS + detail
    grouped = cells.groupby(keys, sort=True, observed=True)
    row_ids = np.full(len(frame), -1, dtype=np.int64)
    row_ids[cells.index] = grouped.ngroup().to_numpy()
    samples = grouped.size().rename("n_pc").reset_index()
    samples.insert(0, "row", range(len(samples)))
    samples["pooled"] = samples.donor.str.lower().eq("pooled")
    samples["at_min_pc"] = samples.n_pc.ge(int(spec["min_pc"]))
    if "size_factor" in file["obs"]:
        factors = np.asarray(read_elem(file["obs/size_factor"]), dtype=float)
        denoms = factors[keep.to_numpy(dtype=bool)]
        if not np.isfinite(denoms).all() or (denoms <= 0).any():
            raise ValueError("invalid size_factor")
        samples["sum_size_factor"] = np.bincount(row_ids[cells.index], weights=denoms, minlength=len(samples))
    layer = sparse_dataset(file[spec["count_layer"]])
    counts = sparse.csr_matrix((len(samples), layer.shape[1]), dtype=np.int64)
    for start in range(0, len(frame), 4096):
        stop = min(start + 4096, len(frame))
        rows = row_ids[start:stop]
        local = np.flatnonzero(rows >= 0)
        if len(local) == 0:
            continue
        block = layer[start:stop][local].tocsr()
        if not np.isfinite(block.data).all() or (block.data < 0).any() or (block.data != np.floor(block.data)).any():
            raise ValueError(f"{spec['count_layer']} has non-count values")
        assign = sparse.csr_matrix((np.ones(len(local), dtype=np.int64), (rows[local], np.arange(len(local)))),
                                   shape=(len(samples), len(local)))
        counts += assign @ block.astype(np.int64)
    var_path = "raw/var" if spec["count_layer"].startswith("raw/") else "var"
    var = read_elem(file[var_path])
    ids = var.index.astype(str)
    if not ids.is_unique or len(ids) != counts.shape[1]:
        raise ValueError("feature identifiers must match columns and be unique")
    features = pd.DataFrame({"column": range(len(ids)), "feature_id": ids})
    if "feature_name" in var:
        features["symbol"] = var["feature_name"].astype("string").to_numpy()
    samples["library_counts"] = np.asarray(counts.sum(axis=1)).ravel()
    prefix = spec["cohort"]
    samples.to_csv(out / f"{prefix}_samples.tsv", index=False, sep="\t")
    features.to_csv(out / f"{prefix}_features.tsv", index=False, sep="\t")
    sparse.save_npz(out / f"{prefix}_counts.npz", counts)
    return int(samples.n_pc.sum())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=ROOT / "inputs/reference")
    parser.add_argument("--out", type=Path, default=ROOT / "outputs")
    with (ROOT / "config/reference.tsv").open() as stream:
        configs = list(csv.DictReader(stream, delimiter="\t"))
    parser.add_argument("--pseudobulk", choices=[cfg["cohort"] for cfg in configs],
                        help="also export one cohort's PC count matrix, sample table and feature table")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    strata_tables, summaries = [], []
    with TemporaryDirectory(dir=args.out) as pending:
        tmp = Path(pending)
        for spec in configs:
            with h5py.File(args.data_root / spec["filename"], "r") as file:
                shape = tuple(file[spec["count_layer"]].attrs["shape"])
                if shape != (int(spec["n_obs"]), int(spec["n_vars"])):
                    raise ValueError(f"Unexpected dimensions for {spec['cohort']}: {shape}")
                frame = read_obs(file, spec)
                strata, summary = census(frame, spec)
                if spec["cohort"] == args.pseudobulk:
                    count = pseudobulk(file, frame, spec, tmp)
                    if count != summary["n_pc_reference"]:
                        raise ValueError("pseudobulk cell count differs from census")
                strata_tables.append(strata)
                summaries.append(summary)
        pd.concat(strata_tables, ignore_index=True).to_csv(tmp / "donor_stage_coverage.tsv", index=False, sep="\t")
        pd.DataFrame(summaries).to_csv(tmp / "reference_summary.tsv", index=False, sep="\t")
        for path in tmp.iterdir():
            os.replace(path, args.out / path.name)
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
