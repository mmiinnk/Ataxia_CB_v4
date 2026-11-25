#!/usr/bin/env python
"""s01 — 细胞级 → donor cube → bin 表(implementation_spec §5 产物定义;鼠先跑, review 门 1)。

用法:
    python -I scripts/s01_build_cubes.py <root> --species mouse [--config PATH] [--min-cells N]

流程(spec §3.4/§3.5/§5):
  read_full → X 计数整数性抽样 → panel 映射(ENSG/ENSMUSG 主键) → derive_obs → assign →
  check_anchors(§4.1 精确, FAIL 即停) → use_cell 过滤(敏感性) → normalize(size_factor, 1e4) →
  切列到 panel∪sanity 基因 → donor_cube → bin_table(全骨架, NA 政策) → 可行性矩阵 vs §4.2
  (WARN 级) → spot_check(n=3) → marker 校验(spec §9 清单 ERROR 级, 其余 WARN) → 趋势锚点
  (鼠, WARNING 级, 正式统计在 s02) → P2c MT 基因检出率(CL-023 冻结行动项)。

产物(results/v<ver>/intermediate/<sp>/):
  cube_donor.h5ad   donor×bin×group pseudobulk(X=mean log1p CP10K; layers mean_cp10k/pct_pos)
  cube_bin.tsv      全骨架 10 组×bin(log_/cp_/pct_/logsem_<gene>; feasible/single_donor/n_donors/n_cells_total)
  cells_panel.h5ad  细胞级 log 面板矩阵(s04/s09 状态剖面用)
日志(logs/): s01_report_<sp>.tsv, s01_feasibility_<sp>.tsv, s01_markers_<sp>.tsv,
  s01_mt_detection_<sp>.tsv, s01_spotcheck_<sp>.tsv, s01_gene_mapping_<sp>.tsv
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sparse   # 别名避开 sp(=species 变量)遮蔽(2026-11-04 P2 真跑踩中)

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from atx import config as atxconfig
from atx import genes as atxgenes
from atx import groups as atxgroups
from atx import io as atxio
from atx import pseudobulk as atxpb
from atx import sanity as atxsanity

# spec §9 review 门 marker 清单(ERROR 级)。鼠免疫 marker 按 dropout 豁免(§3.12, s00 重算),
# 故 PTPRC 仅在人侧 ERROR;其余 marker(含 ATOH1)WARN 级留档人工 review。
SPEC9_MARKERS = {"ITPR1", "CA8", "GABRA6", "CBLN3", "MBP", "PLP1", "AQP4", "SOX2", "MKI67", "PTPRC"}
# spec §9(4) 趋势方向锚点: 定义在 config species.<sp>.trend_anchors(硬约束 6;
# 2026-11-05 人侧裁决后自本文件常量 MOUSE_TREND_ANCHORS 迁入, 经 cfg.trend_anchors(sp) 读取)

REPORT: list[dict] = []
N_FAIL = 0


def rec(check, expected, observed, level, note=""):
    global N_FAIL
    REPORT.append({"check": check, "expected": str(expected), "observed": str(observed),
                   "level": level, "note": note})
    mark = {"PASS": "✓", "WARN": "⚠", "FAIL": "✗"}.get(level, "?")
    print(f"[{mark} {level}] {check} | 期望={expected} 实测={observed}" + (f" | {note}" if note else ""))
    if level == "FAIL":
        N_FAIL += 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--species", required=True, choices=["human", "mouse"])
    ap.add_argument("--config", default=None)
    ap.add_argument("--min-cells", type=int, default=None, help="覆盖 config bins.min_cells(敏感性)")
    ap.add_argument("--out-subdir", default=None,
                    help="产物子目录 results/v{ver}/<subdir>/…; 敏感性变体隔离用, 防覆盖正式产物")
    ap.add_argument("--sens", action="append", default=None, metavar="KEY=VALUE",
                    help="进程内覆盖 config sensitivity 开关(可多次), 如 --sens collapse_adult_bins=true")
    args = ap.parse_args()
    sp = args.species
    t0 = time.time()
    cfg = atxconfig.load(args.config or os.path.join(os.path.abspath(args.root), "config", "config.yaml"))
    atxconfig.apply_sens_overrides(cfg, args.sens)      # 必须在 effective_species/results 之前
    if args.out_subdir:
        cfg.results_subdir = args.out_subdir.strip("/")
    cfg.check_version_gate()
    cfg.write_draft_marker()
    spcfg = cfg.effective_species(sp)
    min_cells = args.min_cells if args.min_cells is not None else int(cfg.bins_cfg["min_cells"])
    inter = cfg.results("intermediate", sp, mkdir=True)
    logs = cfg.results("logs", mkdir=True)
    # collapse_adult_bins 是人侧专属开关(effective_species 只对 human 生效) → 跑鼠时不计入"敏感性开关",
    # 否则鼠侧日志会误报非默认口径。
    sens_on = {k: v for k, v in cfg.sensitivity.items()
               if v not in (False, "mixed") and k not in ("min_cells_grid", "onset_thresholds_grid")
               and not (k == "collapse_adult_bins" and sp != "human")}
    # 口径分类(2026-11-24 用户裁决: collapse_adult_bins 升为人侧**主口径**, 门 2 整改项):
    #  · cells_default = 细胞选择口径未变 → §4.1「9 组细胞数锚点」(ERROR 级)仍**必须**比对;
    #    该锚点按组计数, 与 bin 怎么切无关, 过去被 default_mode 一并跳过是过宽的。
    #  · bins_default  = bin 轴与 spec §4.2 可行性锚点表一致 → 才可比对;
    #    人侧锚点表按 16 bin 记录, 主口径合并成人 bin 后为 13 bin, 二者不可比 → 跳过并写明原因
    #    (锚点表待 v1.1 按新口径重记, 见 spec §11-7)。
    #  · 任何隔离目录(--out-subdir)运行一律视为非默认口径。
    CELL_KEYS = ("exclude_disease_donors", "drop_glutamatergic_uncertain",
                 "split_granule_progenitor", "include_ubc_in_main")
    cells_default = all(not cfg.sens(k, False) for k in CELL_KEYS) \
        and not getattr(spcfg, "exclude_disease_donors", False)
    collapse_on = (sp == "human") and bool(cfg.sens("collapse_adult_bins", False))
    bins_default = (min_cells == int(cfg.bins_cfg["min_cells"])) and not collapse_on \
        and not cfg.results_subdir

    print(f"== s01 {sp} | version={cfg.version}(draft={cfg.draft}) | min_cells={min_cells}"
          + (f" | 敏感性开关: {sens_on}" if sens_on else "")
          + (f" | 产物隔离: results/v{cfg.version}/{cfg.results_subdir}/" if cfg.results_subdir else ""))

    # ---- 全载(硬约束 4:Sepp 全载;backed 列切片会卡死) ----
    path = cfg.h5ad(spcfg.d["h5ad_key"])
    ad = atxio.read_full(path)
    X, obs_raw, var = ad.X, ad.obs, ad.var
    shape = (X.shape[0], X.shape[1])
    exp_shape = (spcfg.n_cells_expected, spcfg.n_vars_expected)
    if shape != exp_shape:
        raise RuntimeError(f"形状锚点失败: {shape} != {exp_shape}")
    rec("形状锚点", exp_shape, shape, "PASS")
    if not sparse.issparse(X):
        raise RuntimeError("X 非稀疏——与实测结构不符, 拒跑")
    samp = X.data[:1_000_000]
    rec("X 计数整数性(前1M nnz 抽样)", "全整数", "全整数" if np.all(samp == np.floor(samp)) else "非整数",
        "PASS" if np.all(samp == np.floor(samp)) else "FAIL")

    # ---- P2c:MT 基因检出率(全矩阵单遍 bincount, 在切列之前;CL-023 冻结行动项) ----
    n_cells_all = X.shape[0]
    pct_all = np.bincount(X.indices, minlength=X.shape[1]) / n_cells_all * 100.0
    fsym = var["feature_name"].astype(str)
    mt_mask = fsym.str.startswith("MT-") if sp == "human" else fsym.str.startswith("mt-")
    mt_rows = pd.DataFrame({"feature_name": fsym[mt_mask], "pct_pos": np.round(pct_all[mt_mask.values], 4)})
    panel = atxgenes.load_panel(cfg)
    in_panel_keys = set(panel["gene_key"].str.upper())
    mt_rows["in_panel_any_tier"] = mt_rows["feature_name"].str.upper().isin(in_panel_keys)
    mt_rows.to_csv(os.path.join(logs, f"s01_mt_detection_{sp}.tsv"), sep="\t", index=False)
    rec("P2c MT 基因检出(CL-023 行动项)", "留档", f"{len(mt_rows)} 个 MT 基因, pct+ 中位 {mt_rows['pct_pos'].median():.2f}%",
        "PASS", note=f"logs/s01_mt_detection_{sp}.tsv")

    # ---- 基因映射与切列 ----
    aliases = atxgenes.load_aliases(cfg)
    gm = atxgenes.map_to_var(panel, var, species=sp, dataset=sp, aliases=aliases)
    gm.to_csv(os.path.join(logs, f"s01_gene_mapping_{sp}.tsv"), sep="\t", index=False)
    hr = atxgenes.hit_rate(gm, sp)
    thr = float(cfg.gene_set["min_hit_rate"][sp])
    miss = gm.loc[(gm["in_stat_panel"].astype(bool)) & (gm[f"status_{sp}"] != "matched"), "gene_key"].tolist()
    rec("panel 命中率", f"≥{thr}(s00 已门)", f"{hr:.4f}", "PASS" if hr >= thr else "WARN",
        note=f"缺失={miss}")
    matched = gm[gm[f"status_{sp}"] == "matched"].copy()
    if matched[f"row_{sp}"].duplicated().any():
        raise RuntimeError("多基因映射到同一 var 行——映射不唯一, 拒跑")
    order = np.sort(matched[f"row_{sp}"].to_numpy(dtype=int))
    key_by_row = dict(zip(matched[f"row_{sp}"].astype(int), matched["gene_key"]))
    var_sub = var.iloc[order].copy()
    var_sub["gene_key"] = [key_by_row[int(r)] for r in order]
    ann = matched.set_index("gene_key")[["tier", "in_stat_panel", "is_sanity", "mech_class",
                                         "onset_typical_y", "onset_class"]]
    for c in ann.columns:
        var_sub[c] = ann.loc[var_sub["gene_key"], c].to_numpy()
    rec("cube 基因数", "panel∪sanity 命中集", f"{len(order)}", "PASS")

    # ---- obs 派生 + 组指派 + 锚点(ERROR 级) ----
    d = atxgroups.derive_obs(obs_raw, spcfg, cfg)
    use = atxgroups.use_cell_mask(d, None, cfg, spcfg)
    if not use.all():
        rec("use_cell 过滤(敏感性)", "-", f"排除 {int((~use).sum())} 核", "WARN")
        d = d.loc[use].reset_index(drop=True)
        X = X[use]
        obs_raw = obs_raw.loc[use].reset_index(drop=True)
    a, _ = atxgroups.assign(d, cfg)
    if cells_default and use.all():   # §4.1 组细胞数锚点与 bin 切法无关 → 合并成人 bin 后仍须比对
        bad = atxgroups.check_anchors(a, sp, cfg)
        rec("9 组细胞数锚点 §4.1", "精确匹配", "匹配" if not bad else bad, "PASS" if not bad else "FAIL")
        if bad:
            _finish(cfg, sp, t0)
            return 1
    lg = atxgroups.expand_long(a, cfg)

    # ---- 归一化(Sepp 用作者 size_factor;硬约束 4) + 切列 ----
    sf = pd.to_numeric(obs_raw["size_factor"], errors="coerce").to_numpy(dtype=np.float64)
    target = float(cfg.raw["normalization"]["target_sum"])
    X_sub = X[:, order]
    X_log, X_cp = atxpb.normalize(X_sub, sf, target)
    rec("归一化", f"X/size_factor×{target:g} → log1p", f"X_log nnz={X_log.nnz}", "PASS")
    del X, X_sub, samp

    # ---- donor cube ----
    cube = atxpb.donor_cube(X_log, X_cp, d, lg, var_sub)
    cube.obs["period"] = cube.obs["bin"].map(spcfg.period_of)
    cube.write_h5ad(os.path.join(inter, "cube_donor.h5ad"))
    rec("donor_cube", "-", f"{cube.n_obs} combos × {cube.n_vars} genes", "PASS",
        note="intermediate/%s/cube_donor.h5ad" % sp)

    # ---- bin 表 + 可行性矩阵 vs §4.2(WARN 级) ----
    bt = atxpb.bin_table(cube, spcfg, cfg, min_cells=min_cells)
    bt.to_csv(os.path.join(inter, "cube_bin.tsv"), sep="\t", index=False, float_format="%.6g")
    ff = atxpb.feasibility_frame(bt)
    ff.to_csv(os.path.join(logs, f"s01_feasibility_{sp}.tsv"), sep="\t")
    if bins_default and use.all():
        diff = atxsanity.check_feas_matrix(ff, sp)
        rec("组×bin 可行性矩阵 §4.2", "全匹配(WARN 级)", "匹配" if diff.empty else f"{len(diff)} 处差异",
            "PASS" if diff.empty else "WARN",
            note="" if diff.empty else diff.head(8).to_string(index=False))
    else:
        reason = ("人侧主口径已合并成人 bin(collapse_adult_bins=true; 2026-11-24 用户裁决, 门 2 整改项): "
                  "spec §4.2 可行性锚点表按 16 bin 记录, 与合并后 13 bin 不可比 → 跳过比对; "
                  "锚点表待 v1.1 按新主口径重记(§11-7)") \
            if collapse_on and min_cells == int(cfg.bins_cfg["min_cells"]) \
            else "非默认口径(敏感性/min-cells 覆盖), 跳过锚点比对"
        rec("组×bin 可行性矩阵 §4.2", "-", reason, "WARN")
    n_feas = int(bt["feasible"].sum())
    rec("bin 表骨架", f"{len(cfg.plot_groups())} 组 × {len(spcfg.labels)} bin",
        f"{len(bt)} 行, 可行 {n_feas}", "PASS")

    # ---- cells_panel.h5ad(s04/s09 用) ----
    import anndata
    obs_cells = d[["donor", "bin", "bin_order", "pooled_bin"]].copy()
    obs_cells["cell_group"] = a["cell_group"].to_numpy()
    # dev_state: 跨物种对齐轴(spec §4.4 共享分化状态; §4.5A 人鼠时间对齐**只**用它, 不做年龄换算),
    # s04 状态剖面与 s09 保真度(traj_r / peak_state_match)的必要列。位置赋值(防索引重对齐, 见 derive_obs 注)。
    if "dev_state" in d.columns:
        obs_cells["dev_state"] = d["dev_state"].to_numpy()
    else:
        raise RuntimeError("obs 无 dev_state 列 —— s04/s09 跨物种对齐轴缺失(spec §4.4), 不静默继续")
    if "disease_flag" in d.columns:
        obs_cells["disease"] = d["disease"].to_numpy()
        obs_cells["disease_flag"] = d["disease_flag"].to_numpy()
    cells = anndata.AnnData(X=X_log, obs=obs_cells, var=var_sub)
    cells.write_h5ad(os.path.join(inter, "cells_panel.h5ad"))
    rec("cells_panel", "-", f"{cells.n_obs} × {cells.n_vars}", "PASS")

    # ---- spot_check(n=3, numpy 直接重算) ----
    sc = atxpb.spot_check(X_log, d, lg, bt, list(var_sub["gene_key"]), n=3, seed=cfg.seed)
    sc.to_csv(os.path.join(logs, f"s01_spotcheck_{sp}.tsv"), sep="\t", index=False)
    ok = (sc["check"] == "PASS").all()
    rec("spot_check 两级聚合重算", "全 PASS", f"{int((sc['check']=='PASS').sum())}/{len(sc)}",
        "PASS" if ok else "FAIL", note="" if ok else sc.to_string(index=False))

    # ---- marker 校验(spec §9:清单内 ERROR, 其余 WARN 留档) ----
    markers = atxsanity.load_markers(cfg.paths["sanity_table"])
    gene_cols = {g: i for i, g in enumerate(var_sub["gene_key"])}
    mt_tab = atxsanity.marker_table(X_cp, a["cell_group"], gene_cols, markers, sp,
                                    cfg.main_groups(), cfg.descriptive_groups())
    mt_tab.to_csv(os.path.join(logs, f"s01_markers_{sp}.tsv"), sep="\t", index=False)
    err_set = set(SPEC9_MARKERS)
    if sp == "mouse":
        err_set -= {"PTPRC", "CSF1R", "CX3CR1", "P2RY12"}   # 鼠免疫 dropout 豁免(§3.12)
    hard = mt_tab[mt_tab["marker"].isin(err_set) & (mt_tab["verdict"] == "FAIL")]
    soft = mt_tab[~mt_tab["marker"].isin(err_set) & (mt_tab["verdict"] == "FAIL")]
    n_abs = int((mt_tab["verdict"] == "ABSENT").sum())
    rec("marker 校验(spec §9 清单)", "期望组=argmax(或≥50%argmax)",
        f"FAIL={len(hard)}(清单内) / {len(soft)}(清单外), ABSENT={n_abs}",
        "PASS" if len(hard) == 0 and len(soft) == 0 else ("FAIL" if len(hard) else "WARN"),
        note="" if len(hard) == 0 else hard[["marker", "expected_group", "argmax_group",
                                             "mean_expected", "mean_argmax"]].to_string(index=False))

    # ---- 趋势方向锚点(config species.<sp>.trend_anchors;WARNING 级, 两物种同路径) ----
    anchors = cfg.trend_anchors(sp)
    if anchors:
        from scipy.stats import spearmanr
        for a in anchors:
            gene, group, sign = a["gene"], a["group"], a["sign"]
            col = f"log_{gene}"
            if col not in bt.columns:
                rec(f"趋势锚点 {gene}×{group}", f"rho {sign}", "基因不在 cube(缺失/未映射)", "WARN")
                continue
            sub = bt[(bt["group"] == group) & bt["feasible"] & bt[col].notna()]
            if len(sub) < int(cfg.bins_cfg["min_bins_for_trend"]):
                rec(f"趋势锚点 {gene}×{group}", f"rho {sign}", f"可行 bin 仅 {len(sub)}", "WARN")
                continue
            rho = float(spearmanr(sub["bin_order"], sub[col]).statistic)
            good = (rho > 0) if sign == "+" else (rho < 0)
            rec(f"趋势锚点 {gene}×{group}", f"rho {sign}", f"rho={rho:+.3f}(n={len(sub)} bin, 粗算)",
                "PASS" if good else "WARN", note="正式 rho/q 在 s02(levels=bin,donor)")

    print(f"\n== s01 {sp} 完成 ({time.time()-t0:.0f}s)")
    return _finish(cfg, sp, t0)


def _finish(cfg, sp, t0):
    df = pd.DataFrame(REPORT, columns=["check", "expected", "observed", "level", "note"])
    p = cfg.results("logs", f"s01_report_{sp}.tsv", mkdir=True)
    df.to_csv(p, sep="\t", index=False)
    n_pass = int((df["level"] == "PASS").sum())
    n_warn = int((df["level"] == "WARN").sum())
    print(f"汇总: PASS={n_pass} WARN={n_warn} FAIL={N_FAIL}  报告: {p}")
    if N_FAIL:
        print("s01 失败(FAIL 即停);WARN 项须在 review 门逐条看过(spec §9)。")
    return 1 if N_FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
