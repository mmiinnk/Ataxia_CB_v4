#!/usr/bin/env python
"""dev_preview — 阶段性预览诊断图(开发轨辅助工具, **非** 流水线正式产物)。

定位:P5 期间应导师要求,把 P2(cube)/P3(trend+specificity)/P4(pattern+模块)在真实
草稿数据上的结果可视化,供人工直观观察;正式图面口径以 P6 的 s03_figures_species 为准,
本脚本产物不进入 Fig/S 组装、不得引用数值(全程 DRAFT 水印)。

纪律:
- 只读 intermediate/<sp>/ 的 s01 产物;**不写 tables/**(正式六表由 run_mouse 在门 1
  流程内落盘)——统计全部经 s02_stats.build_tables 在内存中计算(顺带真数据冒烟 s02);
- 输出仅 results/v<ver>/preview/*.png;
- 颜色/组序/NA 样式/cmap 全部读 config figures(硬约束 6);图内文字用英文
  (运行环境无 CJK 字体, 中文会渲染为方框)。

用法:
    python -I scripts/dev_preview.py <root> --species mouse [--config PATH]

产物(6 张):
    <sp>_p2_cube.png         可行性/细胞数热图 + 每 bin 可行主组数(k≥4 → tau 可算)
    <sp>_p3_stats.png        rho_bin 热图(rho_desc) + tau×max_fraction + emergence 条带
    <sp>_p4_label_matrix.png temporal_pattern 标签矩阵(行按 G* 分块, 左侧 G* 色带)
    <sp>_p4_modules_traj.png 模块 G*×P* 分布 + §9 三个趋势锚点基因轨迹(峰/半峰窗)
    <sp>_p6_b2_<group>.png   P6 正式口径 b2 轨迹小倍数样张(atx.plots, 首组首页)
    <sp>_p6_c1.png           P6 正式口径 c1 dotplot 样张(atx.plots)
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from atx import config as atxconfig
from atx import genes as atxgenes
from atx import plots as atxplots
from atx import stats_enrich as atxenrich
from atx import stats_pattern as atxpattern
from atx import stats_trend as atxtrend
import s02_stats
import s03_figures_species as s03


def _draft_mark(fig, cfg, y=0.004):
    fig.text(0.995, y, f"DRAFT {cfg.version} · dev_preview diagnostic · not for citation",
             ha="right", va="bottom", fontsize=8, color="#cc0000", alpha=0.85)


def _hatch_na(ax, mask, cfg, xoff=0.5, yoff=0.5):
    """NA 格 → config figures.na_style 斜纹(读图约定:表 NA、图斜纹,绝不插值)。"""
    st = dict(cfg.raw["figures"]["na_style"])
    for i, j in zip(*np.where(mask)):
        ax.add_patch(Rectangle((j - xoff, i - yoff), 1, 1, linewidth=0.4, **st))


def _save(fig, cfg, name, t0):
    p = cfg.results("preview", name, mkdir=True)
    fig.savefig(p, dpi=int(cfg.raw["figures"]["dpi"]), bbox_inches="tight")
    plt.close(fig)
    print(f"  [written] {p}  ({time.time()-t0:.1f}s)")
    return p


def _figcfg(cfg):
    return cfg.raw["figures"]


# ---------------------------------------------------------------- P2: cube 总览
def fig_p2_cube(bt, cfg, sp, spcfg):
    t0 = time.time()
    fc = _figcfg(cfg)
    groups = [g for g in cfg.plot_groups() if g in set(bt["group"])]
    bins = (bt[["bin", "bin_order", "period"]].drop_duplicates()
            .sort_values("bin_order").reset_index(drop=True))
    piv = bt.pivot(index="group", columns="bin", values="n_cells_total").reindex(index=groups, columns=bins["bin"])
    feas = bt.pivot(index="group", columns="bin", values="feasible").reindex(index=groups, columns=bins["bin"]).fillna(False).astype(bool)
    M = piv.to_numpy(dtype=float).copy()
    M[~feas.to_numpy()] = np.nan

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.2), gridspec_kw={"width_ratios": [2.4, 1]})
    lm = np.log10(np.maximum(M, 1))
    cmap = plt.get_cmap(fc["sequential_cmap"]).copy()
    cmap.set_bad(fc["na_style"]["facecolor"])
    im = ax1.imshow(np.ma.masked_invalid(lm), cmap=cmap, aspect="auto")
    _hatch_na(ax1, np.isnan(M), cfg)
    ax1.set_yticks(range(len(groups)), groups, fontsize=9)
    ax1.set_xticks(range(len(bins)), bins["bin"], rotation=45, ha="right", fontsize=8)
    ax1.set_title(f"{sp} cube_bin: cells per group×bin (log10), hatched = infeasible\n"
                  f"min_cells={cfg.bins_cfg['min_cells']} · DRAFT", fontsize=10)
    fig.colorbar(im, ax=ax1, shrink=0.8, label="log10 n_cells")

    main = cfg.main_groups()
    k = feas.reindex(index=main).fillna(False).to_numpy().sum(axis=0)
    min_groups = int(cfg.stats_cfg["tau_by_bin"]["min_groups"])
    pcolors = fc["period_colors"]
    ax2.bar(range(len(bins)), k, color=[pcolors[bins["period"][i]] for i in range(len(bins))],
            edgecolor="black", linewidth=0.4)
    ax2.axhline(min_groups, color="red", ls="--", lw=1)
    ax2.text(len(bins) - 0.4, min_groups + 0.12, f"min_groups={min_groups} (tau computable)",
             ha="right", fontsize=8, color="red")
    ax2.set_xticks(range(len(bins)), bins["bin"], rotation=45, ha="right", fontsize=8)
    ax2.set_ylim(0, len(main) + 0.6)
    ax2.set_ylabel("feasible main groups (k)")
    ax2.set_title("Per-bin group feasibility", fontsize=10)
    handles = [Patch(facecolor=pcolors[p], edgecolor="k", label=p)
               for p in dict.fromkeys(bins["period"])]
    ax2.legend(handles=handles, fontsize=7, loc="lower right")
    fig.suptitle(f"[{sp}] P2 preview — s01 cube structure (draft gene table {cfg.version})", fontsize=11)
    _draft_mark(fig, cfg)
    return _save(fig, cfg, f"{sp}_p2_cube.png", t0)


# ---------------------------------------------------------------- P3: trend + specificity
def fig_p3_stats(tabs, cfg, sp, spcfg, genes_shown):
    t0 = time.time()
    fc = _figcfg(cfg)
    if str(fc["row_order_rule"]) != "rho_desc":
        raise NotImplementedError(f"figures.row_order_rule={fc['row_order_rule']!r}: 预览仅实现 rho_desc")
    tr = tabs["trend_by_group"]
    groups = list(dict.fromkeys(tr["group"]))
    anchor_grp = cfg.main_groups()[0]                     # rho_desc 排序锚组 = 主组序首个
    a = tr[tr["group"].to_numpy() == anchor_grp].set_index("gene")["rho_bin"].reindex(genes_shown)

    def _order_key(g):
        v = a.get(g, np.nan)
        v = float(v) if v is not None and np.isfinite(v) else np.nan
        return (0, -v, g) if np.isfinite(v) else (1, 0.0, g)   # 有限 rho 降序在前, NaN 靠右缘

    order = sorted(genes_shown, key=_order_key)
    R = np.full((len(groups), len(order)), np.nan)
    trm = tr.set_index(["group", "gene"])["rho_bin"]
    for gi, g in enumerate(groups):
        for gj, gene in enumerate(order):
            v = trm.get((g, gene), np.nan)
            R[gi, gj] = v if v is not None else np.nan

    sg = tabs["specificity_global"].set_index("gene").reindex(genes_shown)
    em = tabs["emergence"].set_index("gene").reindex(genes_shown)

    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.32, wspace=0.22)
    ax1 = fig.add_subplot(gs[0, :])
    cmap = plt.get_cmap(fc["diverging_cmap"]).copy()
    cmap.set_bad(fc["na_style"]["facecolor"])
    im = ax1.imshow(np.ma.masked_invalid(R), cmap=cmap, vmin=-1, vmax=1, aspect="auto")
    _hatch_na(ax1, np.isnan(R), cfg)
    ax1.set_yticks(range(len(groups)), groups, fontsize=9)
    ax1.set_xticks([])
    ax1.set_xlabel(f"{len(order)} core genes, ordered by {anchor_grp} rho_bin descending "
                   f"(figures.row_order_rule=rho_desc); NaN rho at right edge", fontsize=8)
    ax1.set_title(f"rho_bin (Spearman, bin-level) — {sp}", fontsize=10)
    fig.colorbar(im, ax=ax1, shrink=0.7, label="rho_bin")

    ax2 = fig.add_subplot(gs[1, 0])
    spcfg_s = cfg.stats_cfg["specificity"]
    calls = list(fc["call_colors"])
    for c in calls:
        m = sg["specific_call"].to_numpy() == c
        if m.any():
            ax2.scatter(sg["max_fraction"].to_numpy(dtype=float)[m], sg["tau"].to_numpy(dtype=float)[m],
                        s=16, color=fc["call_colors"][c], label=f"{c} (n={int(m.sum())})",
                        edgecolors="none", alpha=0.85)
    ax2.axvline(float(spcfg_s["tau_threshold"]), color="grey", ls="--", lw=0.8)
    ax2.axhline(float(spcfg_s["fraction_threshold"]), color="grey", ls="--", lw=0.8)
    ax2.set_xlabel("max_fraction")
    ax2.set_ylabel("tau")
    ax2.set_xlim(0, 1)
    ax2.set_title("Specificity (global, 8 main groups)\ndashed = config thresholds", fontsize=10)
    ax2.legend(fontsize=7, loc="upper left")

    ax3 = fig.add_subplot(gs[1, 1])
    main = cfg.main_groups()
    hit = em[em["emergence_bin"].astype(str) != "NA"]
    ys, xs, cols, opens = [], [], [], []
    yidx = {g: i for i, g in enumerate(main)}
    for mg in main:
        sub = hit[hit["max_group"].to_numpy() == mg].sort_index()
        n = len(sub)
        for rank, (gene, r) in enumerate(sub.iterrows()):
            eo = r["emergence_bin_order"]
            if not np.isfinite(eo):
                continue
            xs.append(float(eo))
            ys.append(yidx[mg] + (np.linspace(-0.3, 0.3, n)[rank] if n > 1 else 0.0))
            cols.append(fc["group_colors"][mg])
            opens.append(bool(r["left_censored"]))
    ax3.scatter([x for x, o in zip(xs, opens) if not o], [y for y, o in zip(ys, opens) if not o],
                s=22, c=[c for c, o in zip(cols, opens) if not o], edgecolors="k", linewidths=0.3, label="emergence")
    ax3.scatter([x for x, o in zip(xs, opens) if o], [y for y, o in zip(ys, opens) if o],
                s=26, facecolors="none",
                edgecolors=[c for c, o in zip(cols, opens) if o], linewidths=0.9,
                label="left-censored")
    ax3.set_yticks(range(len(main)), main, fontsize=8)
    ax3.set_xticks(range(1, len(spcfg.labels) + 1), spcfg.labels, rotation=45, ha="right", fontsize=7)
    ax3.set_xlabel("bin")
    ax3.set_title(f"Emergence: first bin with fraction≥{cfg.stats_cfg['tau_by_bin']['fraction_threshold']} "
                  f"in max_group\n(n={int(len(hit))} genes with emergence)", fontsize=9)
    ax3.legend(fontsize=7, loc="upper left")
    fig.suptitle(f"[{sp}] P3 preview — s02 trend & specificity (draft {cfg.version}, DRAFT numbers)", fontsize=11)
    _draft_mark(fig, cfg)
    return _save(fig, cfg, f"{sp}_p3_stats.png", t0)


# ---------------------------------------------------------------- P4a: 标签矩阵
def fig_p4_label_matrix(tabs, mods, cfg, sp, genes_shown):
    t0 = time.time()
    fc = _figcfg(cfg)
    vocab = list(fc["pattern_colors"])
    tp = tabs["temporal_patterns"]
    groups = list(dict.fromkeys(tp["group"]))
    tpm = tp.set_index(["gene", "group"])["temporal_pattern"]
    gstar = mods.set_index("gene")["g_star"].reindex(genes_shown).astype(str)
    cats = [g for g in cfg.main_groups()] + ["broad", "NA"]
    order = sorted(genes_shown, key=lambda g: (cats.index(gstar[g]) if gstar[g] in cats else len(cats), g))
    M = np.full((len(order), len(groups)), -1, dtype=int)
    for i, gene in enumerate(order):
        for j, grp in enumerate(groups):
            lab = tpm.get((gene, grp), "NA")
            lab = str(lab) if lab is not None else "NA"
            if lab not in vocab:
                raise ValueError(f"temporal_pattern 受控词表外标签: {lab!r}")
            M[i, j] = vocab.index(lab)
    strip_colors = {g: fc["group_colors"][g] for g in cfg.main_groups()}
    strip_colors["broad"] = fc["call_colors"]["broad"]      # 口径注: broad 无专属色, 复用 call_colors.broad
    strip_colors["NA"] = fc["call_colors"]["NA"]
    S = np.array([[cats.index(gstar[g]) if gstar[g] in cats else len(cats)] for g in order], dtype=float)
    scmap = ListedColormap([strip_colors[c] for c in cats])

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12.5, 0.17 * len(order) + 2.2),
                                   gridspec_kw={"width_ratios": [0.035, 1]}, sharey=True)
    ax0.imshow(S, cmap=scmap, aspect="auto", vmin=cats.index(cats[0]) - 0.5, vmax=len(cats) - 0.5)
    ax0.set_xticks([])
    ax0.set_title("G*", fontsize=8)
    lc = ListedColormap([fc["pattern_colors"][v] for v in vocab])
    ax1.imshow(M, cmap=lc, aspect="auto", vmin=-0.5, vmax=len(vocab) - 0.5)
    ax1.set_xticks(range(len(groups)), groups, rotation=45, ha="right", fontsize=9)
    ax1.set_yticks([])
    ax1.tick_params(left=False)
    for i in range(0, len(order), 2):                          # 基因名标在 G* 色带左侧
        ax0.text(-0.6, i, order[i], ha="right", va="center", fontsize=4.6,
                 transform=ax0.get_yaxis_transform())
    ax1.set_title(f"temporal_pattern — {len(order)} core genes in {sp} cube × {len(groups)} groups "
                  f"(rows grouped by G*)", fontsize=10)
    handles = [Patch(facecolor=fc["pattern_colors"][v], edgecolor="k", linewidth=0.3, label=v) for v in vocab]
    handles += [Patch(facecolor=strip_colors[c], edgecolor="k", linewidth=0.3, label=f"G*={c}")
                for c in cats if (gstar == c).any()]
    fig.legend(handles=handles, loc="lower center", ncol=7, fontsize=7, frameon=False)
    fig.suptitle(f"[{sp}] P4 preview — temporal_pattern labels (draft {cfg.version})", fontsize=11)
    _draft_mark(fig, cfg, y=0.985)      # 底部留给图例
    return _save(fig, cfg, f"{sp}_p4_label_matrix.png", t0)


# ---------------------------------------------------------------- P4b: 模块分布 + 锚点轨迹
def fig_p4_modules_traj(bt, tabs, mods, cfg, sp, spcfg, genes_shown):
    t0 = time.time()
    fc = _figcfg(cfg)
    trend_anchors = [(a["gene"], a["group"], a["sign"]) for a in cfg.trend_anchors(sp)]
    # §9(4) 门清单锚点(单一来源 = config species.<sp>.trend_anchors, 硬约束 6)
    cats = [g for g in cfg.main_groups() if g in set(mods["g_star"])] + \
           [c for c in ("broad", "NA") if c in set(mods["g_star"])]
    periods = s02_stats._periods(spcfg)
    msub = mods.set_index("gene").reindex(genes_shown)
    counts = {(c, p): int(((msub["g_star"] == c) & (msub["peak_period"] == p)).sum())
              for c in cats for p in periods + ["NA"]}

    fig = plt.figure(figsize=(15, 8.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.1], hspace=0.42, wspace=0.25)
    axb = fig.add_subplot(gs[0, :])
    bottom = np.zeros(len(cats))
    for p in periods + ["NA"]:
        vals = np.array([counts[(c, p)] for c in cats], dtype=float)
        col = fc["period_colors"].get(p, fc["pattern_colors"]["NA"])
        axb.bar(range(len(cats)), vals, bottom=bottom, color=col, edgecolor="black",
                linewidth=0.4, label=p)
        for xi, (v, b0) in enumerate(zip(vals, bottom)):
            if v > 0:
                axb.text(xi, b0 + v / 2, int(v), ha="center", va="center", fontsize=8)
        bottom += vals
    axb.set_xticks(range(len(cats)), cats, rotation=20, fontsize=9)
    axb.set_ylabel("genes")
    axb.set_title(f"Module G*×P* distribution ({int(len(msub))} core genes in cube; "
                  f"broad = max_fraction<{cfg.stats_cfg['enrichment']['broad_fraction_threshold']})", fontsize=10)
    axb.legend(fontsize=8, ncol=len(periods) + 1, loc="upper right")

    for k, (gene, group, sign) in enumerate(trend_anchors[:3]):
        ax = fig.add_subplot(gs[1, k])
        col = fc["group_colors"][group]
        sub = bt[(bt["group"].to_numpy() == group) & bt["feasible"].to_numpy()]
        sub = sub.sort_values("bin_order")
        lcol = f"log_{gene}"
        if lcol not in bt.columns:
            ax.set_title(f"{gene} × {group}: not in cube", fontsize=9)
            continue
        m = sub[lcol].notna().to_numpy()
        x = sub["bin_order"].to_numpy(dtype=float)[m]
        y = pd.to_numeric(sub.loc[m, lcol]).to_numpy(dtype=float)
        ax.plot(x, y, "-o", color=col, ms=3.5, lw=1.4)
        # 锚点含 sanity marker(MBP/ATOH1 不在核心行集)→ 峰/趋势现算, 不查六表
        pk_g = atxpattern.peak_table(bt, [gene], cfg)
        pk_g = pk_g[pk_g["group"].to_numpy() == group]
        prow = pk_g.iloc[0] if len(pk_g) else None
        pats = atxpattern.patterns(bt, [gene], group, cfg).iloc[0]
        trow = atxtrend.trend_table(bt, None, [gene], group, spcfg, cfg).iloc[0]
        if prow is not None and str(prow["status"]) == "ok":
            lo = spcfg.bin_order[str(prow["half_max_lo"])]
            hi = spcfg.bin_order[str(prow["half_max_hi"])]
            ax.axvspan(lo - 0.5, hi + 0.5, color=col, alpha=0.15)
            ax.plot([pats["peak_bin_order"]], [pats["peak_log"]], marker="*", ms=15,
                    color=col, mec="black", mew=0.5, zorder=5)
        rho, q = trow["rho_bin"], trow["q_bin"]
        ok = np.isfinite(rho) and ((rho > 0) if sign == "+" else (rho < 0))
        ax.set_title(f"{gene} × {group}: rho={rho:+.3f} q={q:.1g} "
                     f"{'✓' if ok else '✗ anchor ' + sign}\npeak={pats['peak_bin']} "
                     f"({pats['peak_period']}), shaded=half-max window", fontsize=8.5)
        ax.set_xticks(range(1, len(spcfg.labels) + 1), spcfg.labels, rotation=45, ha="right", fontsize=6)
        ax.set_xlabel("bin")
        ax.set_ylabel("mean log1p CP10K")
    fig.suptitle(f"[{sp}] P4 preview — modules & §9 trend anchors (draft {cfg.version})", fontsize=11)
    _draft_mark(fig, cfg)
    return _save(fig, cfg, f"{sp}_p4_modules_traj.png", t0)


# ---------------------------------------------------------------- P6: 正式口径样张
def fig_p6_formal_views(tabs, bt, cube, cfg, sp, spcfg, genes_shown):
    """P6 正式口径样张: b2 轨迹小倍数(首个细胞组, rho 前 12 基因) + c1 dotplot。

    直接调 atx.plots / s03 的正式原语(与 s03_figures_species 同一代码路径),
    仅落 preview/*.png 并加 DRAFT 水印——供人工提前目检正式画风。
    """
    t0 = time.time()
    out = []
    tr, sg = tabs["trend_by_group"], tabs["specificity_global"]
    cell_groups = s03._cell_groups(cfg, bt)
    g0 = cell_groups[0]                                       # 主组序首个(Purkinje)
    labels = list(spcfg.labels)
    per = int(cfg.raw["figures"]["traj_panels_per_page"])
    rho_g = (tr.loc[tr["group"].to_numpy() == g0].set_index("gene")["rho_bin"]
             .reindex(genes_shown).to_numpy(dtype=float))
    shown = [x for x in atxplots.order_rho_desc(rho_g, genes_shown)
             if np.isfinite(atxplots.group_bin_matrix(bt, g0, [x], labels, "log_")).any()][:per]
    fig, _infos = atxplots.traj_figure(
        bt, cube, shown, g0, cfg, spcfg,
        title=f"[{sp}] P6 preview · b2 {g0} trajectories — top {len(shown)} genes by rho "
              f"(formal s03 style)")
    _draft_mark(fig, cfg)
    out.append(_save(fig, cfg, f"{sp}_p6_b2_{atxplots.safe_name(g0)}.png", t0))

    dot_groups = [g for g in cell_groups
                  if g in set(cfg.main_groups()) | set(cfg.descriptive_groups())]
    c1_order = [g for g in atxplots.order_main_tau(sg, cfg.main_groups()) if g in set(genes_shown)]
    fig2, info = atxplots.dotplot_figure(bt, sg, c1_order, dot_groups, cfg, species=sp)
    _draft_mark(fig2, cfg)
    print(f"  c1 dotplot: {info['n_dots']} dots ({info['n_genes']} genes × {info['n_groups']} groups)")
    out.append(_save(fig2, cfg, f"{sp}_p6_c1.png", t0))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--species", default="mouse", choices=["human", "mouse"])
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    sp = args.species
    t0 = time.time()
    cfg = atxconfig.load(args.config or os.path.join(os.path.abspath(args.root), "config", "config.yaml"))
    cfg.check_version_gate()
    spcfg = cfg.effective_species(sp)
    inter = cfg.results("intermediate", sp)
    bt_path = os.path.join(inter, "cube_bin.tsv")
    if not os.path.exists(bt_path):
        raise RuntimeError(f"缺 {bt_path}——先跑 s01")
    print(f"== dev_preview {sp} | version={cfg.version}(draft={cfg.draft}) | 只读, 仅写 preview/*.png")
    bt = pd.read_csv(bt_path, sep="\t")

    panel = atxgenes.load_panel(cfg)
    universe = panel.loc[panel["in_stat_panel"].astype(bool), "gene_key"].astype(str).tolist()
    genes_shown = [g for g in universe if f"log_{g}" in bt.columns]
    print(f"  核心基因 {len(universe)},cube 命中 {len(genes_shown)}(预览仅展示命中集)")

    tabs = s02_stats.build_tables(bt, None, universe, spcfg, cfg)   # 内存计算(= s02 真数据冒烟)
    mods = atxenrich.assign_modules(tabs["specificity_global"], tabs["peak_by_group"], cfg)
    mods = mods[mods["gene"].isin(genes_shown)]

    # 控制台 DRAFT 摘要(不可引用)
    tp = tabs["temporal_patterns"]
    print("  标签分布(全组×核心基因):", tp.loc[tp["gene"].isin(genes_shown), "temporal_pattern"].value_counts().to_dict())
    sg = tabs["specificity_global"].set_index("gene").reindex(genes_shown)
    print("  specific_call:", sg["specific_call"].value_counts().to_dict())
    print("  G* 分布:", mods["g_star"].value_counts().to_dict())

    paths = [
        fig_p2_cube(bt, cfg, sp, spcfg),
        fig_p3_stats(tabs, cfg, sp, spcfg, genes_shown),
        fig_p4_label_matrix(tabs, mods, cfg, sp, genes_shown),
        fig_p4_modules_traj(bt, tabs, mods, cfg, sp, spcfg, genes_shown),
    ]
    # P6 样张需要 donor cube(灰点/疾病△);缺 s01 产物则跳过并提示
    cube_path = os.path.join(inter, "cube_donor.h5ad")
    if os.path.exists(cube_path):
        import anndata
        cube = anndata.read_h5ad(cube_path)
        paths += fig_p6_formal_views(tabs, bt, cube, cfg, sp, spcfg, genes_shown)
    else:
        print(f"  [skip] P6 样张: 缺 {cube_path}(先跑 s01)")
    print(f"== dev_preview 完成 ({time.time()-t0:.1f}s);{len(paths)} 张 PNG,全部 DRAFT 水印,数值不得引用。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
