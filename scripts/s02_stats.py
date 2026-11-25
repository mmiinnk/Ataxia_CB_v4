#!/usr/bin/env python
"""s02 — 六表组装(implementation_spec §5 产物定义;第 1–2 层统计全部落表)。

用法:
    python -I scripts/s02_stats.py <root> --species mouse [--config PATH]

前置: s01 已产出 intermediate/<sp>/{cube_bin.tsv, cube_donor.h5ad}(缺即 RuntimeError)。

流程(spec §5/§3.6–3.9):
  读 cube_bin.tsv + cube_donor.h5ad → 基因全集 = 统计 panel(in_stat_panel, 缺失基因
  各表补 status=missing 行)→ build_tables():
    ① trend_by_group   每 plot 组调 stats_trend.trend_table(q_bin 族 = species×group,
       与 config stats.trend.fdr_family 硬校验);
    ② temporal_patterns 每 plot 组调 stats_pattern.patterns + merge 同组 trend 的
       rho_bin/q_bin → stats_pattern.label(§5 优先级);
    ③ specificity_global stats_specificity.aggregate(全期)+ scores;
    ④ specificity_by_period 每期 aggregate(period=p)+ scores(period 列插入 gene 后);
    ⑤ tau_by_bin + emergence(max_group 来自 ③;emergence = 首个主组 fraction ≥
       tau_by_bin.fraction_threshold 的可算 bin, §3.7 差异③④);
    ⑥ peak_by_group   stats_pattern.peak_table(s07 模块指派/s08 泳道输入)。
  → write_tables(outdir)(float %.6g)→ 报告 logs/s02_report_<sp>.tsv;
  鼠: §9 门清单趋势锚点的**正式** rho/q(s01 为粗算 WARNING 级)——锚点基因含
  sanity marker(MBP/ATOH1 不在统计 panel), 故对锚点基因单独调 trend_table,
  不依赖六表行集;方向不符 = WARN(门 1 人工判读, 脚本不 FAIL 停线)。

产物(results/v<ver>/tables/<sp>/): trend_by_group.tsv, temporal_patterns.tsv,
  specificity_global.tsv, specificity_by_period.tsv, tau_by_bin.tsv, emergence.tsv,
  peak_by_group.tsv
build_tables/write_tables 为纯函数(路径显式传入)——敏感性循环与 dev_preview 复用,
测试不触真实 results/。
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from atx import config as atxconfig
from atx import genes as atxgenes
from atx import stats_pattern as atxpattern
from atx import stats_specificity as atxspec
from atx import stats_trend as atxtrend

# spec §5 L462–472 的六表(tau_by_bin 与 emergence 同一表格行、两个文件 → 7 个 TSV)
TABLE_NAMES = ["trend_by_group", "temporal_patterns", "specificity_global",
               "specificity_by_period", "tau_by_bin", "emergence", "peak_by_group"]

_TP_COLS = ["gene", "group", "n_feasible_bins", "peak_bin", "peak_bin_order", "peak_period",
            "peak_internal", "peak_vs_flank", "period_ratio", "temporal_pattern",
            "peak_log", "peak_cp", "rho_bin", "q_bin", "status"]

REPORT: list[dict] = []


def rec(check, observed, level, note=""):
    REPORT.append({"check": check, "observed": str(observed), "level": level, "note": note})
    mark = {"PASS": "✓", "WARN": "⚠", "FAIL": "✗"}.get(level, "?")
    print(f"[{mark} {level}] {check} | {observed}" + (f" | {note}" if note else ""))


def _periods(spcfg) -> list[str]:
    """period 名按首个 bin 的出现序(确定性;鼠 embryonic→perinatal→postnatal→adult)。"""
    return list(dict.fromkeys(spcfg.period_of[lab] for lab in spcfg.labels))


def build_tables(bt: pd.DataFrame, cube, genes: list[str], spcfg, cfg) -> dict[str, pd.DataFrame]:
    """bin 表(+donor cube, 可 None)→ 六表 dict(键 = 文件名去 .tsv;纯内存, 不写盘)。

    行集 = genes(统计 panel 全集;cube 无列的基因 → 各表 status=missing / 标签 NA)。
    组集 = config plot_groups ∩ bt 实有组(trend/temporal_patterns/peak_by_group);
    特异性分母恒为 main_groups(stats_specificity 内部口径)。
    """
    fam = str(cfg.stats_cfg["trend"].get("fdr_family", ""))
    if fam != "species_group":
        raise NotImplementedError(f"stats.trend.fdr_family={fam!r}: 仅实现 species_group"
                                  "(q_bin = trend_table 单次调用[一物种×一组]内 BH)")
    present = set(bt["group"].unique())
    groups = [g for g in cfg.plot_groups() if g in present]

    # ① trend(FDR 族 = species×group:trend_table 每次调用即一族)
    trend = pd.concat([atxtrend.trend_table(bt, cube, genes, g, spcfg, cfg) for g in groups],
                      ignore_index=True)

    # ② temporal_patterns(patterns + 同组 trend 的 rho_bin/q_bin → label)
    tp_frames = []
    for g in groups:
        p = atxpattern.patterns(bt, genes, g, cfg)
        t = trend.loc[trend["group"].to_numpy() == g, ["gene", "rho_bin", "q_bin"]]
        m = p.merge(t, on="gene", how="left")
        if len(m) != len(p) or m["gene"].duplicated().any():
            raise RuntimeError(f"patterns×trend merge 行数不齐({g})——gene 键重复或缺行")
        m["temporal_pattern"] = atxpattern.label(m, cfg)
        tp_frames.append(m[_TP_COLS])
    temporal_patterns = pd.concat(tp_frames, ignore_index=True)

    # ③④ specificity global / by-period
    specificity_global = atxspec.scores(atxspec.aggregate(bt, genes, cfg), genes, cfg)
    bp = []
    for per in _periods(spcfg):
        s = atxspec.scores(atxspec.aggregate(bt, genes, cfg, period=per), genes, cfg)
        s.insert(1, "period", per)
        bp.append(s)
    specificity_by_period = pd.concat(bp, ignore_index=True)

    # ③b P7 裁决①(2026-11-05): not_expressed 基因的 window_high/transient_peak 伪标签 → NA
    atxpattern.mask_lowexpr_labels(temporal_patterns, specificity_global)

    # ⑤ tau_by_bin + emergence(max_group 来自 global;NA → 全 NA 行, §3.7 差异④)
    tau_bin = atxspec.tau_by_bin(bt, genes, cfg)
    emergence = atxspec.emergence(tau_bin, specificity_global.set_index("gene")["max_group"], cfg)

    # ⑥ peak_by_group(组 = plot_groups ∩ 实有, 含 ALL——s07 broad 基因 P* 取 ALL 行)
    peak = atxpattern.peak_table(bt, genes, cfg)

    return {"trend_by_group": trend, "temporal_patterns": temporal_patterns,
            "specificity_global": specificity_global, "specificity_by_period": specificity_by_period,
            "tau_by_bin": tau_bin, "emergence": emergence, "peak_by_group": peak}


def write_tables(tables: dict[str, pd.DataFrame], outdir: str) -> list[str]:
    """六表 → outdir/*.tsv(float %.6g, 与 s01 一致);返回写出路径列表。"""
    os.makedirs(outdir, exist_ok=True)
    paths = []
    for name in TABLE_NAMES:
        p = os.path.join(outdir, f"{name}.tsv")
        tables[name].to_csv(p, sep="\t", index=False, float_format="%.6g")
        paths.append(p)
    return paths


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--species", required=True, choices=["human", "mouse"])
    ap.add_argument("--config", default=None)
    ap.add_argument("--out-subdir", default=None,
                    help="产物子目录 results/v{ver}/<subdir>/…; 须与 s01 同值(读同一隔离 cube)")
    ap.add_argument("--sens", action="append", default=None, metavar="KEY=VALUE",
                    help="进程内覆盖 config sensitivity 开关(可多次); 须与 s01 同口径, 否则 bin_order 硬校验失败")
    args = ap.parse_args()
    sp = args.species
    t0 = time.time()
    cfg = atxconfig.load(args.config or os.path.join(os.path.abspath(args.root), "config", "config.yaml"))
    atxconfig.apply_sens_overrides(cfg, args.sens)      # 必须在 effective_species/results 之前
    if args.out_subdir:
        cfg.results_subdir = args.out_subdir.strip("/")
    cfg.check_version_gate()
    cfg.write_draft_marker()
    spcfg = cfg.effective_species(sp)          # collapse_adult_bins 敏感性须与 s01 建 cube 时一致
    inter = cfg.results("intermediate", sp)

    print(f"== s02 {sp} | version={cfg.version}(draft={cfg.draft})"
          + (f" | 产物隔离: results/v{cfg.version}/{cfg.results_subdir}/" if cfg.results_subdir else ""))
    bt_path = os.path.join(inter, "cube_bin.tsv")
    cube_path = os.path.join(inter, "cube_donor.h5ad")
    for p in (bt_path, cube_path):
        if not os.path.exists(p):
            raise RuntimeError(f"缺 s01 产物 {p}——先跑 s01(spec §9 顺序)")
    bt = pd.read_csv(bt_path, sep="\t")
    import anndata
    cube = anndata.read_h5ad(cube_path)

    # 基因全集 = 统计 panel(sanity marker 不入六表行集;缺失基因 → status=missing)
    panel = atxgenes.load_panel(cfg)
    genes = panel.loc[panel["in_stat_panel"].astype(bool), "gene_key"].astype(str).tolist()
    in_cube = sum(int(f"log_{g}" in bt.columns) for g in genes)
    rec("基因行集(统计 panel)", f"{len(genes)} 基因, cube 命中 {in_cube}, 缺失 {len(genes)-in_cube}(补 missing 行)", "PASS")

    tables = build_tables(bt, cube, genes, spcfg, cfg)
    outdir = cfg.results("tables", sp, mkdir=True)
    paths = write_tables(tables, outdir)
    for name in TABLE_NAMES:
        df = tables[name]
        note = ""
        if "status" in df.columns:
            vc = df["status"].value_counts().to_dict()
            note = f"status={vc}"
        rec(f"{name}.tsv", f"{df.shape[0]} 行 × {df.shape[1]} 列", "PASS", note=note)

    # ---- 汇总口径检查(NA 政策落点, 门 1 人工判读素材) ----
    tp = tables["temporal_patterns"]
    rec("temporal_pattern 标签分布", tp["temporal_pattern"].value_counts().to_dict(), "PASS")
    n_masked = int(((tp["status"].astype(str) == "ok")
                    & (tp["temporal_pattern"].astype(str) == "NA")).sum())
    rec("低表达伪标签掩膜(P7 裁决①)", f"{n_masked} 行 window/peak→NA(not_expressed 基因)", "PASS")
    sg = tables["specificity_global"]
    rec("specific_call 分布(global)", sg["specific_call"].value_counts().to_dict(), "PASS")
    tr = tables["trend_by_group"]
    n_sig = int(((tr["q_bin"] < 0.05) & tr["q_bin"].notna()).sum())
    rec("trend q_bin<0.05 行数(族=species×group)", f"{n_sig}/{len(tr)}", "PASS")
    em = tables["emergence"]
    n_em = int((em["emergence_bin"] != "NA").sum())
    n_lc = int(em["left_censored"].fillna(False).astype(bool).sum())
    rec("emergence 可判基因", f"{n_em}/{len(em)}(左截断 {n_lc})", "PASS")

    # ---- §9(4) 门清单趋势锚点(正式 rho/q;含 sanity marker 基因 → 单独调 trend_table;
    #      锚点定义 = config species.<sp>.trend_anchors, 两物种同路径[2026-11-05 裁决]) ----
    for anch in cfg.trend_anchors(sp):
        gene, group, sign = anch["gene"], anch["group"], anch["sign"]
        if f"log_{gene}" not in bt.columns:
            rec(f"趋势锚点 {gene}×{group}(正式)", "基因不在 cube", "WARN")
            continue
        row = atxtrend.trend_table(bt, cube, [gene], group, spcfg, cfg).iloc[0]
        rho, q = row["rho_bin"], row["q_bin"]
        good = np.isfinite(rho) and ((rho > 0) if sign == "+" else (rho < 0))
        rec(f"趋势锚点 {gene}×{group}(正式)",
            f"rho {sign} | rho_bin={rho:+.3f} q_bin={q:.3g} status={row['status']}",
            "PASS" if good else "WARN",
            note=f"rho_donor={row['rho_donor']:+.3f}" if np.isfinite(row["rho_donor"]) else "donor 级 NA")

    logs = cfg.results("logs", mkdir=True)
    rp = os.path.join(logs, f"s02_report_{sp}.tsv")
    pd.DataFrame(REPORT, columns=["check", "observed", "level", "note"]).to_csv(rp, sep="\t", index=False)
    print(f"\n== s02 {sp} 完成 ({time.time()-t0:.0f}s);表目录: {outdir}")
    print(f"报告: {rp}")
    print("下一步(实现后): s03_figures_species → 人工 review 门(spec §9)。")
    n_fail = sum(1 for r in REPORT if r["level"] == "FAIL")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
