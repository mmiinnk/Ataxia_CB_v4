"""n00 — 正常 PC metacell 构建(阶段 1 网络层第一步;config network.metacell)。

Sepp 人 PC 核(author_cell_type=='Purkinje')按 donor × stage × subtype 聚合为 metacell:
  · ≥ min_cells_per_metacell(30) 核才成 metacell;不足者记录为不可行(NA 政策:不聚合不插值);
  · cap_large_donor: 每 metacell 降采样至 cap_cells_per_metacell(种子固定)再求和,
    限制 CS22 超大供体(25,993 核/2 供体)的统计权重;
  · 聚合口径 = 项目统一 pseudobulk: CP10K = Σcounts / Σsize_factor × 1e4(作者 size_factor,
    spec normalization 单一口径),取 log1p;
  · 全部参数走 config(硬约束 6);/mnt 只读(硬约束 1);确定性(seed 固定,幂等覆盖)。

产物(results/v{ver}/network/):
  metacells_log.tsv     metacell × 基因 log1p CP10K(列 = ENSG)
  metacells_var.tsv     ENSG → symbol(feature_name)对照
  metacells_meta.tsv    每 metacell: donor/stage/stage_order/subtype/n_cells_total/n_cells_used/capped/disease
  n00_strata_census.tsv 全部 donor×stage×subtype 层(含不可行层)的核数普查
  logs/n00_report.tsv   自动核查(与 s0x_report 同 schema)
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from atx import config as atxconfig            # noqa: E402
from atx import io as atxio                    # noqa: E402

REPORT: list[dict] = []


def rec(check: str, expected: str, observed, level: str, note: str = "") -> None:
    REPORT.append({"check": check, "expected": expected, "observed": str(observed),
                   "level": level, "note": note})
    mark = {"PASS": "✓", "WARN": "⚠", "FAIL": "✗"}.get(level, "?")
    print(f"[{mark} {level}] {check} | 期望={expected} 实测={observed}" + (f" | {note}" if note else ""))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    root = os.path.abspath(args.root)
    cfg = atxconfig.load(args.config or os.path.join(root, "config", "config.yaml"))
    cfg.check_version_gate()
    t0 = time.time()
    mcfg = cfg.raw["network"]["metacell"]
    ccfg = cfg.raw["network"]["coexpression"]
    min_cells = int(mcfg["min_cells_per_metacell"])
    cap = int(mcfg.get("cap_cells_per_metacell") or 0) or None
    # cap_large_donor 真正生效(finding #4 子集): false → 不降采样截尾
    if not bool(mcfg.get("cap_large_donor", True)):
        cap = None
    min_donors_stage = int(ccfg["min_donors_per_stage"])
    sp = cfg.species("human")
    # 列名走 SpeciesCfg(硬约束 6; finding #4 子集), subtype 为实测列名无 config 字段
    time_col, donor_col = sp.time_col, sp.donor_col

    outdir = cfg.results("network", mkdir=True)
    logs = cfg.results("logs", mkdir=True)

    print(f"== n00 metacells | version={cfg.version} | strata={mcfg['strata']} | "
          f"min_cells={min_cells} cap={cap}")
    ad = atxio.read_full(cfg.h5ad("human_h5ad"))
    obs = ad.obs
    atxio.assert_integer_counts(ad.X.data, "X")          # Sepp X = 原始整数计数(spec normalization)
    mask = (obs["author_cell_type"].astype(str) == "Purkinje").to_numpy()
    pos_all = np.flatnonzero(mask)
    pc = obs.iloc[pos_all].copy()
    pc["_pos"] = pos_all
    rec("PC 核普查", "与图谱线锚点一致(45,341)", len(pc),
        "PASS" if len(pc) == 45341 else "FAIL")

    # NA 安全化(pandas 3 陷阱, 硬约束 7; 同 atx/groups.py 模式): astype(str) 不再把 NA 转 "nan",
    # 直接 replace 会漏 → groupby 静默丢弃全部 subtype=NA 的 PC 核(约 1.5 万, 含 CS18/19/17w/20w)。
    sub_raw = pc["subtype"]
    pc["subtype"] = (sub_raw.astype("object").where(sub_raw.notna(), "unannotated").astype(str)
                     .replace({"nan": "unannotated", "__NA__": "unannotated", "": "unannotated"}))
    assert (pc["subtype"] != "nan").all()
    pc["stage_label"] = pc["development_stage"].astype(str).map(sp.label_of)
    assert pc["stage_label"].notna().all(), "存在 time_order 外的 development_stage"
    sf_all = pd.to_numeric(obs["size_factor"], errors="coerce").to_numpy(dtype=float)
    # disease NA 哨兵化(finding #7): astype(str) 在 pandas3 不把 NA 转 "nan", sorted(set(含NA)) 会 TypeError
    disease_all = (obs["disease"].astype("object").where(obs["disease"].notna(), "__NA__")
                   .astype(str).to_numpy())

    rng = np.random.default_rng(cfg.seed)
    var_index = np.asarray(ad.var.index).astype(str)
    mats, meta, census = [], [], []
    grouped = pc.groupby([time_col, donor_col, "subtype"], observed=True, sort=True)
    for (stage, donor, sub), g in grouped:
        cells = g["_pos"].to_numpy()
        n = int(len(cells))
        feasible = n >= min_cells
        cens_row = {"stage": stage, "stage_label": str(sp.label_of[stage]), "donor": donor,
                    "subtype": sub, "n_cells": n, "feasible": feasible, "used": False}
        census.append(cens_row)
        if not feasible:
            continue
        used = cells
        capped = bool(cap is not None and n > cap)
        if capped:
            used = np.sort(rng.choice(cells, cap, replace=False))
        sub_counts = np.asarray(ad.X[used].sum(axis=0)).ravel()
        sf_sum = float(np.nansum(sf_all[used]))
        if not np.isfinite(sf_sum) or sf_sum <= 0:
            # finding #6: 该层 feasible 但被跳过 → used 保持 False, 覆盖率不计入(此前会虚高)
            rec(f"metacell {donor}|{stage}|{sub}", "size_factor 和 > 0", sf_sum, "FAIL")
            continue
        cp = sub_counts / sf_sum * float(cfg.raw["normalization"]["target_sum"])
        mats.append(np.log1p(cp).astype(np.float32))
        diseases = sorted(set(disease_all[used]))
        meta.append({"metacell_id": f"{donor}|{sp.label_of[stage]}|{sub}",
                     "donor": donor, "stage": stage, "stage_label": str(sp.label_of[stage]),
                     "stage_order": int(sp.bin_order[str(sp.label_of[stage])]),
                     "subtype": sub, "n_cells_total": n, "n_cells_used": int(len(used)),
                     "capped": capped, "disease": ";".join(diseases)})
        cens_row["used"] = True

    M = np.vstack(mats)
    meta_df = pd.DataFrame(meta)
    census_df = pd.DataFrame(census)
    n_dropped = int((~census_df.feasible).sum())
    n_skipped = int((census_df.feasible & ~census_df.used).sum())
    nuclei_used = int(census_df.loc[census_df.used, "n_cells"].sum())   # 仅实际建成 metacell 的核
    # 防护核查: 层普查合计必须等于全部 PC 核(groupby 默认 dropna, NaN 键会静默丢行)
    rec("层普查完整性", f"= 全部 PC 核({len(pc)})", int(census_df.n_cells.sum()),
        "PASS" if int(census_df.n_cells.sum()) == len(pc) else "FAIL",
        note="守卫: 任何 NaN 键静默丢弃在此现形(2026-11-26 曾捕获 subtype NA 丢弃 14,856 核)")

    pd.DataFrame(M, index=meta_df.metacell_id, columns=var_index).to_csv(
        os.path.join(outdir, "metacells_log.tsv"), sep="\t", float_format="%.6g")
    pd.DataFrame({"ensg": var_index,
                  "symbol": np.asarray(ad.var["feature_name"]).astype(str)}).to_csv(
        os.path.join(outdir, "metacells_var.tsv"), sep="\t", index=False)
    meta_df.to_csv(os.path.join(outdir, "metacells_meta.tsv"), sep="\t", index=False)
    census_df.to_csv(os.path.join(outdir, "n00_strata_census.tsv"), sep="\t", index=False)

    # ---- 覆盖与可行性核查 ----
    rec("metacell 数", "≥ 建网最低样本量(约 12)", len(meta_df), "PASS" if len(meta_df) >= 12 else "WARN",
        note=f"不可行层 {n_dropped} 个(核数<{min_cells}, 记 NA 不聚合); size_factor 异常跳过 {n_skipped} 个")
    rec("核覆盖率", "实际建成 metacell 的核 / 全部 PC 核", f"{nuclei_used}/{len(pc)} = {nuclei_used / len(pc):.1%}",
        "PASS", note="未覆盖部分为 <30 核的稀薄层(NB/3y/42y 等), 与 proposal 预期一致")
    n_capped = int(meta_df.capped.sum())
    rec("cap_large_donor", f"降采样 cap={cap}", f"{n_capped} 个 metacell 被截尾",
        "PASS" if n_capped >= 1 else "WARN", note="CS22 大供体应触发 cap")
    dis_rows = meta_df[meta_df.disease != "normal"]
    rec("疾病供体 metacell", "0(发现集应为正常参照)", len(dis_rows),
        "PASS" if len(dis_rows) == 0 else "WARN",
        note=";".join(dis_rows.metacell_id) if len(dis_rows) else "全部 normal")
    stage_donors = meta_df.groupby("stage_label")["donor"].nunique()
    ok_stages = stage_donors[stage_donors >= min_donors_stage].index.tolist()
    rec(f"min_donors_per_stage≥{min_donors_stage} 的分期", "≥1 个分期(建网可行)",
        f"{len(ok_stages)} 个: {ok_stages}", "PASS" if ok_stages else "FAIL",
        note="供体数不足的分期仍进 metacell(全轴相关), 但网络稳定性主要由这些分期支撑")
    donors_total = meta_df.donor.nunique()
    rec("供体数", "≥5(LODO 可执行)", int(donors_total), "PASS" if donors_total >= 5 else "WARN")

    pd.DataFrame(REPORT).to_csv(os.path.join(logs, "n00_report.tsv"), sep="\t", index=False)
    n_fail = sum(1 for r in REPORT if r["level"] == "FAIL")
    print(f"\n汇总: PASS={sum(1 for r in REPORT if r['level'] == 'PASS')} "
          f"WARN={sum(1 for r in REPORT if r['level'] == 'WARN')} FAIL={n_fail}  "
          f"产物: {outdir}  ({time.time() - t0:.0f}s)")
    print(meta_df.groupby(["stage_label", "subtype"], observed=True)
          .agg(n_metacells=("metacell_id", "size"), n_donors=("donor", "nunique"),
               cells_used=("n_cells_used", "sum")).to_string())
    del ad
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
