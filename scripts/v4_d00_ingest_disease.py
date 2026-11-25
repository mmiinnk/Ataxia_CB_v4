"""v4_d00 疾病数据摄入与盘点（implementation_spec §C-d00；首轮施工图步骤 7）。

职责（--inventory 模式，只读 /mnt，绝不写入 data_root）：
  1. 读 config/disease_datasets.tsv（格式/路径/计数层/独立单位/阻塞登记）；
  2. 读各队列已计算的 QC 表（inputs/disease/<GSE>/…，经符号链接指向 /mnt）：
     GSE246183 pc_coverage_by_individual.tsv（逐个体 PC/zebrin 计数）、
     SCP1300 sample_mapping_qc.tsv + published_pc_coverage.tsv（含 GSE282079 供体重叠列）、
     GSE282079 sample_mapping_qc.tsv（318 文库/71 供体/细胞型）；
  3. 产出 outputs/disease_inventory.tsv（逐数据集事实+阻塞+角色）；
  4. 产出 inputs/disease_manifest.tsv（统一清单：参考四库 sha256 + /mnt 上游 disease sha256
     + 本地衍生 QC 表 sha256 + 解压 h5ad 的"待哈希"登记）；
  5. 预填工作登记表 inputs/sample_manifest.tsv（templates/ 下同名文件保持表头-only 规范）。

硬约束：/mnt（--data-root）只读；大文件（>1GB）不哈希，登记 status=large_hash_pending；
SCP1300 与 GSE282079 的供体重叠必须显式报告（两 A-T 轴独立性证据）；
GSE269430 受 R/Seurat 转换门控（本机无 R），holdout 冻结前不解盲。
纯标准库实现：omics 与 pd_tremor_cl 两环境均可运行。
"""
import argparse
import csv
import hashlib
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HASH_LIMIT = 1 << 30  # 1 GB：超过只登记字节数，不哈希

INV_FIELDS = ["dataset", "study_id", "species", "tissue", "format", "cell_resolution",
              "independent_unit", "n_individuals", "n_cases", "n_controls", "n_libraries",
              "pc_support", "counts_layer", "key_files", "blockers", "role", "access_status"]
MAN_FIELDS = ["file", "bytes", "sha256", "source", "status"]
SAM_FIELDS = ["study_id", "sample_id", "library_id", "individual_id", "genotype_or_diagnosis",
              "control_group", "species", "repeat_or_variant_class", "brain_region",
              "age_or_model_stage", "sex", "batch", "postmortem_interval",
              "cell_annotation_source", "pc_nuclei", "pooled_or_technical_replicate",
              "counts_layer", "qc_status", "reason"]


def sha256(path: Path) -> str:
    with path.open("rb") as fh:
        return hashlib.file_digest(fh, "sha256").hexdigest()


def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def manifest_row(path: Path, source: str, status: str = "ok"):
    if not path.exists():
        return {"file": str(path), "bytes": "", "sha256": "MISSING", "source": source,
                "status": "missing"}
    size = path.stat().st_size
    if size > HASH_LIMIT:
        return {"file": str(path), "bytes": size, "sha256": "NA", "source": source,
                "status": "large_hash_pending"}
    return {"file": str(path), "bytes": size, "sha256": sha256(path), "source": source,
            "status": status}


def inventory(data_root: Path):
    disease_root = ROOT / "inputs" / "disease"
    inv, manifest, samples = [], [], []

    # ---- 统一清单：参考四库 + 上游疾病清单 ----
    for row in read_tsv(ROOT / "inputs" / "sha256.tsv"):
        manifest.append({"file": row["file"], "bytes": row["bytes"], "sha256": row["sha256"],
                         "source": "v3_reference_manifest", "status": "ok"})
    mnt_manifest = data_root / "disease_input_sha256.tsv"
    if mnt_manifest.exists():
        for row in read_tsv(mnt_manifest):
            manifest.append({"file": f"inputs/disease/{row['file']}", "bytes": row["bytes"],
                             "sha256": row["sha256"], "source": row.get("source", "geo_upstream"),
                             "status": "ok"})

    # ---- GSE246183：逐个体 PC 覆盖（人 discovery / 鼠 timeline）----
    pc246 = read_tsv(disease_root / "GSE246183" / "pc_coverage_by_individual.tsv")
    for species in ("human", "mouse"):
        rows = [r for r in pc246 if r["species"] == species]
        cases = [r for r in rows if r["group_for_analysis"] != "control"]
        ctrl = [r for r in rows if r["group_for_analysis"] == "control"]
        elig = [r for r in rows if r["eligible_sub4_ge30"] == "True"]
        inv.append({
            "dataset": f"GSE246183_{species}", "study_id": "GSE246183", "species": species,
            "tissue": "cerebellar_cortex" if species == "human" else "cerebellum",
            "format": "h5ad(extracted)", "cell_resolution": "snRNA-seq",
            "independent_unit": "donor" if species == "human" else "animal",
            "n_individuals": len(rows), "n_cases": len(cases), "n_controls": len(ctrl),
            "n_libraries": sum(int(r["n_libraries"]) for r in rows),
            "pc_support": f"{len(elig)} individuals eligible(PC>=30)",
            "counts_layer": "raw/X",
            "key_files": "GSE246183_*_*.h5ad;pc_coverage_by_individual.tsv",
            "blockers": "冗余 .gz 与解压 h5ad 并存；解压 h5ad 不在上游 sha 清单（见 manifest 待哈希行）",
            "role": "discovery(SCA1)" if species == "human" else "timeline(SCA1)",
            "access_status": "ingested"})
        for r in rows:
            samples.append({
                "study_id": "GSE246183", "sample_id": r["individual_id"],
                "library_id": r["geo_accessions"], "individual_id": r["individual_id"],
                "genotype_or_diagnosis": r["group_for_analysis"],
                "control_group": "control" if r["group_for_analysis"] == "control" else "",
                "species": species,
                "repeat_or_variant_class": "polyQ_cag" if cases and r["group_for_analysis"] != "control" else "",
                "brain_region": "cerebellar_cortex" if species == "human" else "cerebellum",
                "age_or_model_stage": r["age_or_stage"], "sex": r["sex"], "batch": "",
                "postmortem_interval": r["postmortem_hours"],
                "cell_annotation_source": "author_ctype_h5ad", "pc_nuclei": r["n_pc_sub4"],
                "pooled_or_technical_replicate": "technical_merged" if int(r["n_libraries"]) > 1 else "single",
                "counts_layer": "raw/X",
                "qc_status": "eligible" if r["eligible_sub4_ge30"] == "True" else "low_pc",
                "reason": r["qc_notes"]})
    for name in ("GSE246183_230502_human_ctype.h5ad", "GSE246183_211027_mouse.h5ad"):
        manifest.append(manifest_row(disease_root / "GSE246183" / name, "extracted_local",
                                     "extracted_not_in_upstream_manifest"))

    # ---- SCP1300：供体级 + 与 GSE282079 的供体重叠（独立性证据）----
    sc = read_tsv(disease_root / "SCP1300" / "sample_mapping_qc.tsv")
    pub = read_tsv(disease_root / "SCP1300" / "published_pc_coverage.tsv")
    sex_conf = [r for r in sc if r["sex_portal"] != r["sex_table_s1"]]
    overlap = [r for r in pub if r.get("gse282079_candidate_donor")]
    overlap_pc = [r for r in overlap if r.get("gse282079_candidate_has_pc_library") == "True"]
    cases = [r for r in sc if r["diagnosis"] == "A-T"]
    inv.append({
        "dataset": "SCP1300", "study_id": "SCP1300", "species": "human",
        "tissue": "cerebellar_vermis", "format": "tsv.gz raw counts", "cell_resolution": "snRNA-seq",
        "independent_unit": "donor", "n_individuals": len(sc), "n_cases": len(cases),
        "n_controls": len(sc) - len(cases), "n_libraries": len(sc),
        "pc_support": f"{sum(1 for r in sc if r['pc_at_least_30'] == 'True')} donors PC>=30",
        "counts_layer": "raw_counts",
        "key_files": "20210217_ATCB_expressionmatrix.tsv.gz;sample_mapping_qc.tsv",
        "blockers": f"{len(sex_conf)} 例性别冲突按 sex_for_adjustment 裁定；与 GSE282079 重叠供体 "
                    f"{len(overlap)}（其中 {len(overlap_pc)} 有 PC 文库）→ 两 A-T 轴独立性须显式报告",
        "role": "discovery(A-T 轴一)", "access_status": "ingested"})
    for r in sc:
        samples.append({
            "study_id": "SCP1300", "sample_id": r["biosample_id"], "library_id": r["biosample_id"],
            "individual_id": r["donor_id"], "genotype_or_diagnosis": r["diagnosis"],
            "control_group": "" if r["diagnosis"] == "A-T" else "control", "species": "human",
            "repeat_or_variant_class": "biallelic_lof" if r["diagnosis"] == "A-T" else "",
            "brain_region": r["organ"], "age_or_model_stage": r["age"],
            "sex": r["sex_for_adjustment"], "batch": "", "postmortem_interval": "",
            "cell_annotation_source": "portal_cluster", "pc_nuclei": r["n_pc"],
            "pooled_or_technical_replicate": "single", "counts_layer": "raw_counts",
            "qc_status": "eligible" if r["pc_at_least_30"] == "True" else "low_pc",
            "reason": "sex_portal_vs_table_s1_conflict" if r["sex_portal"] != r["sex_table_s1"] else ""})

    # ---- GSE282079：NETSseq 分选群（第二条 A-T 轴，伪汇总 only）----
    net = read_tsv(disease_root / "GSE282079" / "sample_mapping_qc.tsv")
    donors = {r["donor_id"] for r in net}
    at_donors = {r["donor_id"] for r in net if r["diagnosis"] != "Control"}
    pc_libs = [r for r in net if r["cell_type"] == "Purkinje"]
    inv.append({
        "dataset": "GSE282079", "study_id": "GSE282079", "species": "human", "tissue": "cerebellum",
        "format": "NETSseq sorted-population count tables", "cell_resolution": "sorted populations(非单细胞)",
        "independent_unit": "donor", "n_individuals": len(donors), "n_cases": len(at_donors),
        "n_controls": len(donors - at_donors), "n_libraries": len(net),
        "pc_support": f"{len(pc_libs)} Purkinje libraries / {len({r['donor_id'] for r in pc_libs})} donors",
        "counts_layer": "counts_txt",
        "key_files": "GSE282079_RAW.tar;sample_mapping_qc.tsv",
        "blockers": "非单细胞：仅伪汇总/细胞组成分析，不能聚类或亚型；与 SCP1300 供体重叠见 SCP1300 行",
        "role": "discovery(A-T 轴二，独立复现)", "access_status": "ingested"})
    for r in net:
        samples.append({
            "study_id": "GSE282079", "sample_id": r["geo_accession"], "library_id": r["library_id"],
            "individual_id": r["donor_id"], "genotype_or_diagnosis": r["diagnosis"],
            "control_group": "control" if r["diagnosis"] == "Control" else "", "species": "human",
            "repeat_or_variant_class": "biallelic_lof" if r["diagnosis"] != "Control" else "",
            "brain_region": "cerebellum", "age_or_model_stage": r["age"], "sex": r["sex"],
            "batch": "", "postmortem_interval": r["postmortem_hours"],
            "cell_annotation_source": r["cell_type"],
            "pc_nuclei": r["nuclei"] if r["cell_type"] == "Purkinje" else "",
            "pooled_or_technical_replicate": "single", "counts_layer": "counts_txt",
            "qc_status": "ok", "reason": r["territory_as_published"]})

    # ---- GSE269430：SCA7 留出，R/Seurat 转换门 ----
    inv.append({
        "dataset": "GSE269430", "study_id": "GSE269430", "species": "mouse",
        "tissue": "cerebellum_pc_enriched", "format": "Seurat .rds.gz x2 (RAW.tar)",
        "cell_resolution": "PC-enriched snRNA-seq", "independent_unit": "animal(MULTI-seq 还原)",
        "n_individuals": "", "n_cases": "", "n_controls": "", "n_libraries": 2,
        "pc_support": "PC 富集", "counts_layer": "rds",
        "key_files": "GSE269430_RAW.tar;MULTI_seq_barcodes.csv.gz",
        "blockers": "R/Seurat 转换门：本机无 R，转换或取得 h5ad 前不读；本地 5wk/8wk 与方案 8.5 周须核对；"
                    "holdout：预测锁定前不看表达效应",
        "role": "holdout", "access_status": "ingested(conversion gated)"})
    for gsm, stage in (("GSM8315551", "5wk"), ("GSM8315552", "8wk")):
        samples.append({
            "study_id": "GSE269430", "sample_id": gsm, "library_id": gsm, "individual_id": "",
            "genotype_or_diagnosis": "SCA7-266Q", "control_group": "", "species": "mouse",
            "repeat_or_variant_class": "polyQ_cag", "brain_region": "cerebellum_pc_enriched",
            "age_or_model_stage": stage, "sex": "", "batch": "", "postmortem_interval": "",
            "cell_annotation_source": "seurat_rds", "pc_nuclei": "",
            "pooled_or_technical_replicate": "multi_seq_restore_pending", "counts_layer": "rds",
            "qc_status": "blocked_rds_conversion",
            "reason": "MULTI-seq barcode→animal 还原待转换后执行"})

    # ---- 本地衍生 QC 表入统一清单 ----
    for rel in ("GSE246183/sample_mapping_qc.tsv", "GSE246183/pc_coverage_by_individual.tsv",
                "SCP1300/sample_mapping_qc.tsv", "SCP1300/published_pc_coverage.tsv",
                "GSE282079/sample_mapping_qc.tsv"):
        manifest.append(manifest_row(disease_root / rel, "derived_local_qc"))

    write_tsv(ROOT / "outputs" / "disease_inventory.tsv", inv, INV_FIELDS)
    write_tsv(ROOT / "inputs" / "disease_manifest.tsv", manifest, MAN_FIELDS)
    write_tsv(ROOT / "inputs" / "sample_manifest.tsv", samples, SAM_FIELDS)

    print(f"inventory datasets : {len(inv)}")
    print(f"manifest rows      : {len(manifest)} "
          f"(large_hash_pending={sum(1 for m in manifest if m['status'] == 'large_hash_pending')})")
    print(f"sample register    : {len(samples)} rows "
          f"({sum(1 for s in samples if s['qc_status'] == 'eligible')} eligible, "
          f"{sum(1 for s in samples if s['qc_status'].startswith('blocked'))} blocked)")
    for row in inv:
        print(f"  - {row['dataset']:<18} {row['role']:<32} n_ind={row['n_individuals']} "
              f"pc={row['pc_support']}")
    print("SCA3 提示：本地家系/EGA bulk/GSE327324/GSE107958 见 config/cohorts.tsv（盘点中，未入库）。")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", action="store_true", help="运行盘点/清单/登记预填")
    parser.add_argument("--data-root", type=Path, default=Path("/mnt/data/xyy/Ataxia_CB"),
                        help="只读数据根（仅用于读上游 disease_input_sha256.tsv）")
    args = parser.parse_args()
    if not args.inventory:
        parser.error("目前仅支持 --inventory")
    inventory(args.data_root)


if __name__ == "__main__":
    main()
