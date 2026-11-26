#!/usr/bin/env python3
"""a7_report_figures.py — 《执行情况报告》正文配图（3 张；论文风格；图内不含说明文字与假设检验结果）。

图注（方法、口径、统计量、判读）单独成文：docs/执行情况报告_图注.md。
图内只保留：面板字母、简短标题、坐标轴名/刻度、图例、数据标签。
术语与 docs/执行情况报告_修改后版本.txt 一致（发育分期 / 细胞类群 / 阶段性高表达 /
早期发育程序群、晚期功能成熟群 / 跨物种发育表达一致性评分 等）。

版面（A4 印刷 1:1，页边距 25mm → 正文宽 160mm）：
  图1 160×105mm  a/b 小鼠与人 细胞类群×发育分期 细胞核数矩阵；c 临床异质性四联条图
  图2 160×118mm  a 机制类别表达高峰分布（按类着色+范围线+中位）；
                 b 各类群峰值基因数；c 具有明确发育变化模式的基因占比
  图3 160×122mm  a–d GRM1 与 ATXN3 的人/鼠轨迹（三色实线）；
                 e 疾病级一致性：四指标符号矩阵 + 综合评分棒棒糖图
输出 results/v{ver}/figures/report/report_fig{1,2,3}.{pdf,png}。

配色（dataviz 六项检验实测；不考虑黑白印刷，但仍须色盲安全）：
  · 机制类别 7 色 = Okabe-Ito 系 #0072B2,#56B4E9,#009E73,#E69F00,#D55E00,#CC79A7,#8C510A
    （硬指标全过；#56B4E9/#E69F00/#CC79A7 对白底对比度<3:1 的 WARN 由"每行直接标注类名+色块"补足）。
    项目原 mech_colors 未过（#ff7f0e↔#2ca02c 红色盲 ΔE=0.7），弃用。
  · 轨迹三色 #0072B2,#D55E00,#009E73 六项全过。
  · 人/鼠双系列 #1f77b4/#d62728 六项全过。
  · 有序量=单色相顺序色（Blues/viridis）；缺失=斜纹；e 面板综合评分=viridis 顺序色。
用法: <omics python> -I scripts/a7_report_figures.py <root> [--fig 1|2|3]
"""
from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import gridspec
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atx import config as atxconfig                      # noqa: E402

def _setup_fonts() -> str:
    """字体：首选**微软雅黑**（用户指定）。

    本机（Linux）通常无该字体 → 自动回退 Noto Sans CJK SC（视觉最接近的黑体）。
    若在项目 fonts/ 目录放入 msyh.ttf / msyh.ttc / Microsoft_YaHei.* 之一，则自动注册并优先使用；
    在装有雅黑的机器（如 Windows/Office 环境）上重跑本脚本亦直接得雅黑版。
    返回实际解析到的字体族名，供运行时打印留痕。
    """
    import glob as _glob
    from matplotlib import font_manager as fm
    fdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fonts")
    for pat in ("*yahei*", "*YaHei*", "*msyh*", "*MSYH*"):
        for p in _glob.glob(os.path.join(fdir, pat)):
            try:
                fm.fontManager.addfont(os.path.abspath(p))
            except Exception as e:                      # 字体文件损坏等不静默吞掉
                print(f"[a7] 注册字体失败 {p}: {e}")
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "Microsoft YaHei UI",
                                       "Noto Sans CJK SC", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["pdf.fonttype"] = 42        # PDF 嵌 TrueType, 矢量可编辑
    plt.rcParams["ps.fonttype"] = 42         # EPS 同（Type 42 嵌入, 避免 CJK 转 Type3 出错）
    try:   # 字体缺失时 findfont 会抛异常 → 视为未解析到雅黑
        path = fm.findfont(fm.FontProperties(family="Microsoft YaHei"), fallback_to_default=False)
    except Exception:
        path = ""
    resolved = os.path.basename(path) if path else "Noto Sans CJK SC（回退）"
    return resolved, bool(path) and ("yahei" in path.lower() or "msyh" in path.lower())


FONT_RESOLVED, FONT_IS_YAHEI = _setup_fonts()

MM = 1 / 25.4
C_HUMAN, C_MOUSE = "#1f77b4", "#d62728"
C_INK = "#2f3a45"
C_BAR = "#4C72B0"
MECH_COLORS = {"dna_repair": "#0072B2", "cilia_congenital_dev": "#56B4E9",
               "mito_metabolic": "#E69F00", "proteostasis_RNA": "#CC79A7",
               "polyQ_repeat": "#D55E00", "channel_synapse": "#009E73",
               "cytoskeleton_transport_glia": "#8C510A"}
TRAJ_COLORS = {"ALL": "#0072B2", "Purkinje": "#D55E00", "GranuleLineage": "#009E73"}
TRAJ_LINES = [("ALL", "全小脑"), ("Purkinje", "浦肯野细胞"), ("GranuleLineage", "颗粒细胞谱系")]

GROUP_CN = {"Purkinje": "浦肯野细胞", "GranuleLineage": "颗粒细胞谱系", "DN": "深部核团神经元",
            "Interneuron": "中间神经元", "OligoLineage": "少突胶质谱系", "Astro": "星形胶质细胞",
            "Progenitor": "祖细胞", "Immune": "免疫细胞", "UBC": "单刷状细胞", "ALL": "全小脑"}
GROUP_ORDER = ["Purkinje", "GranuleLineage", "DN", "Interneuron", "OligoLineage",
               "Astro", "Progenitor", "Immune", "UBC", "ALL"]
MAIN8 = GROUP_ORDER[:8]
MECH_CN = {"mito_metabolic": "线粒体代谢类", "channel_synapse": "离子通道与突触类",
           "cilia_congenital_dev": "纤毛与先天发育类", "proteostasis_RNA": "蛋白稳态与RNA代谢类",
           "cytoskeleton_transport_glia": "细胞骨架与运输类", "dna_repair": "DNA修复类",
           "polyQ_repeat": "polyQ类", "other": "其他"}
MECH_ORDER = ["dna_repair", "cilia_congenital_dev", "mito_metabolic", "proteostasis_RNA",
              "polyQ_repeat", "channel_synapse", "cytoskeleton_transport_glia"]
INH_CN = {"AR": "常染色体隐性", "AD": "常染色体显性", "XLR": "X连锁", "XLD": "X连锁"}
ONSET_CN = {"congenital": "先天起病", "childhood": "儿童起病", "juvenile": "青少年起病",
            "adult": "成年起病"}
COURSE_CN = {"progressive": "进行性", "slowly_progressive": "缓慢进行性",
             "static_nonprogressive": "非进行性", "episodic": "发作性", "unknown": "未知"}
PERIOD_CN = {"embryonic": "胚胎期", "perinatal": "围产期", "postnatal": "出生后",
             "fetal": "胎儿期", "infant": "婴儿期", "adult": "成年"}


def label(ax, text, dy=1.04):
    ax.text(-0.05, dy, text, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="bottom", ha="right", color=C_INK)


def na_patch(ax, x, y):
    ax.add_patch(Rectangle((x, y), 1, 1, facecolor="white", hatch="///",
                           edgecolor="#d7dadd", linewidth=0.4, zorder=2))


def clean(ax, grid="y"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(0.5)
    ax.tick_params(length=1.5, labelsize=5.2)
    if grid:
        ax.grid(axis=grid, color="#e8eaec", linewidth=0.5)
        ax.set_axisbelow(True)


def save(fig, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png", "eps"):       # eps 为用户指定交付格式; 矢量, 字体以 Type42 嵌入
        fig.savefig(os.path.join(outdir, f"{stem}.{ext}"), dpi=300)
    plt.close(fig)
    print(f"[a7] {stem} → {outdir}/ (pdf/png/eps)")


# =============================================================== 图1
def fig1(cfg, outdir):
    genes = pd.read_csv(cfg.paths["genes_table"], sep="\t")
    core = genes[genes["tier"] == "core"]
    panel = core[core["in_panel"].astype(str).str.lower().eq("true")]

    fig = plt.figure(figsize=(160 * MM, 105 * MM))
    gs = gridspec.GridSpec(2, 8, figure=fig, height_ratios=[56, 42],
                           left=0.105, right=0.900, top=0.905, bottom=0.105,
                           hspace=0.52, wspace=0.30)
    axes = [fig.add_subplot(gs[0, 0:4]), fig.add_subplot(gs[0, 4:8])]
    cax = fig.add_axes([0.916, 0.475, 0.014, 0.430])

    vmax, im = 1, None
    mats = {}
    for sp in ("mouse", "human"):
        bt = pd.read_csv(cfg.results("intermediate", sp, "cube_bin.tsv"), sep="\t")
        bt = bt.set_index(["group", "bin"])
        spcfg = cfg.effective_species(sp)
        cols = list(spcfg.labels)
        M = np.full((len(GROUP_ORDER), len(cols)), np.nan)
        F = np.zeros((len(GROUP_ORDER), len(cols)), dtype=bool)
        for i, g in enumerate(GROUP_ORDER):
            for j, b in enumerate(cols):
                if (g, b) in bt.index:
                    r = bt.loc[(g, b)]
                    M[i, j], F[i, j] = float(r["n_cells_total"]), bool(r["feasible"])
                    vmax = max(vmax, float(r["n_cells_total"]))
        mats[sp] = (M, F, cols, spcfg)
    vmax_use = max(3e4, vmax)

    for sp, ax in (("mouse", axes[0]), ("human", axes[1])):
        M, F, cols, spcfg = mats[sp]
        im = ax.imshow(np.where(F, M, np.nan), aspect="auto", cmap="Blues",
                       norm=LogNorm(vmin=30, vmax=vmax_use), interpolation="nearest")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if not F[i, j] and np.isfinite(M[i, j]):
                    na_patch(ax, j - 0.5, i - 0.5)
                if np.isfinite(M[i, j]):
                    v = float(M[i, j])
                    frac = (np.log10(max(v, 1)) - np.log10(30)) / (np.log10(vmax_use) - np.log10(30))
                    ax.text(j, i, f"{int(v)}", ha="center", va="center", fontsize=4.3,
                            color=("white" if frac > 0.60 else C_INK) if F[i, j] else "#a8adb3")
        periods = [spcfg.period_of[b] for b in cols]
        uniq = list(dict.fromkeys(periods))
        cmap_p = matplotlib.colormaps["viridis"]
        for j, p in enumerate(periods):
            ax.add_patch(Rectangle((j - 0.5, -1.35), 1, 0.30,
                                   facecolor=cmap_p(uniq.index(p) / max(1, len(uniq) - 1)),
                                   edgecolor="white", linewidth=0.5, clip_on=False))
        for k, p in enumerate(uniq):
            js = [j for j, q in enumerate(periods) if q == p]
            ax.text(np.mean(js), -1.75, PERIOD_CN.get(p, p), ha="center", va="center",
                    fontsize=5.0, color=C_INK, clip_on=False)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(cols, rotation=90, fontsize=5.2)
        ax.set_yticks(range(len(GROUP_ORDER)))
        ax.set_yticklabels([GROUP_CN[g] for g in GROUP_ORDER] if sp == "mouse" else [], fontsize=5.8)
        ax.set_xlim(-0.5, len(cols) - 0.5)
        ax.set_ylim(len(GROUP_ORDER) - 0.5, -2.1)
        ax.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(GROUP_ORDER), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=0.7)
        ax.tick_params(which="both", length=0)
        for s in ax.spines.values():
            s.set_linewidth(0.5)
        ax.set_title(f"{'小鼠' if sp == 'mouse' else '人'}：{len(cols)}个发育分期 × 10个细胞类群",
                     fontsize=7.2, color=C_INK, pad=14)
    label(axes[0], "a")
    label(axes[1], "b")
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("细胞核数（对数刻度）", fontsize=5.8)
    cb.ax.tick_params(labelsize=5.0, length=1.5)
    cb.outline.set_linewidth(0.4)

    d = pd.read_csv(os.path.join(cfg.root, "genes", "ataxia_gene_disease.tsv"),
                    sep="\t", dtype=str).fillna("NA")
    dp = d[d["hgnc_symbol"].isin(set(panel["hgnc_symbol"])) & (d["is_primary_disease"] == "true")]
    sets = [("分子机制类别", core["mech_class"].map(MECH_CN).value_counts(), 8),
            ("遗传方式", core["inheritance"].map(INH_CN).value_counts(), 3),
            ("发病年龄", core["onset_class"].map(ONSET_CN).value_counts(), 4),
            ("病程特征", core["course"].map(COURSE_CN).value_counts(), 5)]
    for k, (title, vc, n) in enumerate(sets):
        ax = fig.add_subplot(gs[1, k * 2:k * 2 + 2])
        vc = vc.head(n)[::-1]
        ax.barh(np.arange(len(vc)), vc.values, height=0.70, color=C_BAR,
                edgecolor="white", linewidth=0.6)
        ax.set_yticks(np.arange(len(vc)))
        ax.set_yticklabels([str(i) for i in vc.index], fontsize=5.2)
        ax.set_xlim(0, vc.values.max() * 1.12)
        for i, v in enumerate(vc.values):
            if v >= vc.values.max() * 0.30:
                ax.text(v - vc.values.max() * 0.02, i, str(int(v)), va="center", ha="right",
                        fontsize=5.2, color="white")
            else:
                ax.text(v + vc.values.max() * 0.03, i, str(int(v)), va="center", ha="left",
                        fontsize=5.2, color=C_INK)
        ax.set_title(title, fontsize=6.6, color=C_INK, pad=3)
        clean(ax, grid="x")
        ax.tick_params(axis="y", length=0)
        if k == 0:
            label(ax, "c")
    fig.suptitle("图1　共济失调致病基因发育表达图谱的构成", fontsize=9, color=C_INK, y=0.985)
    save(fig, outdir, "report_fig1")


# =============================================================== 图2
def fig2(cfg, outdir):
    pr = pd.read_csv(cfg.results("tables", "cross", "peak_rank_by_class.tsv"), sep="\t")
    pr = pr[pr["in_panel"] & ~pr["mech_class"].isin(["other", "unknown"])
            & pr["peak_rank_mouse_only"].notna()].copy()
    ma = pd.read_csv(cfg.results("tables", "cross", "module_assignment.tsv"), sep="\t")
    spcfg = cfg.species("mouse")
    bins = list(spcfg.labels)
    bord = {b: i + 1 for i, b in enumerate(bins)}
    periods = [spcfg.period_of[b] for b in bins]

    fig = plt.figure(figsize=(160 * MM, 118 * MM))
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[58, 44],
                           left=0.140, right=0.985, top=0.895, bottom=0.105,
                           hspace=0.66, wspace=0.30)
    ax_a = fig.add_subplot(gs[0, :])

    # ---- a: 机制类别的表达高峰分布(按类着色 + 范围线 + 中位) ----
    rng = np.random.default_rng(cfg.seed)
    uniq_p = list(dict.fromkeys(periods))
    for k, p in enumerate(uniq_p):
        js = [i for i, q in enumerate(periods) if q == p]
        ax_a.axvspan(min(js) + 0.5, max(js) + 1.5,
                     color="#f2f4f6" if k % 2 == 0 else "white", zorder=0)
        ax_a.text(np.mean(js) + 1, len(MECH_ORDER) + 0.52, PERIOD_CN.get(p, p),
                  ha="center", va="center", fontsize=5.6, color="#6b7480", zorder=3)
    for i, m in enumerate(MECH_ORDER):
        s = pr[pr["mech_class"] == m]
        if not len(s):
            continue
        y = len(MECH_ORDER) - i
        col = MECH_COLORS[m]
        if i % 2:
            ax_a.axhspan(y - 0.5, y + 0.5, color="#fafbfc", zorder=0)
        x = s["peak_bin_mouse"].map(bord).astype(float).to_numpy()
        ax_a.hlines(y, x.min(), x.max(), color=col, lw=1.1, alpha=0.5, zorder=3)
        ax_a.scatter(x, y + rng.uniform(-0.22, 0.22, len(x)), s=14, color=col,
                     edgecolor="white", linewidth=0.4, zorder=4)
        med = float(np.median(x))
        # 中位用该类专属色的细竖线(不用粗黑线, 避免遮挡散点)
        ax_a.plot([med, med], [y - 0.36, y + 0.36], color=col, lw=1.1, zorder=5,
                  solid_capstyle="butt")
        ax_a.text(len(bins) + 0.8, y, bins[min(len(bins) - 1, max(0, int(round(med)) - 1))],
                  va="center", ha="left", fontsize=5.2, color=C_INK)
        # 行身份由"行标签 + 同色散点/范围线"承担, 不再另画色块(色块会压住行标签)
    ax_a.set_yticks(range(1, len(MECH_ORDER) + 1))
    ax_a.set_yticklabels([MECH_CN[m] for m in reversed(MECH_ORDER)], fontsize=6.0)
    ax_a.set_ylim(0.4, len(MECH_ORDER) + 1.05)
    ax_a.set_xlim(0.4, len(bins) + 2.4)
    ax_a.set_xticks(range(1, len(bins) + 1))
    ax_a.set_xticklabels(bins, fontsize=5.4)
    ax_a.set_xlabel("小鼠发育分期", fontsize=6.2)
    clean(ax_a, grid=None)
    ax_a.set_title("不同分子机制类别基因的表达高峰分布（小鼠）", fontsize=7.4, color=C_INK, pad=5)
    label(ax_a, "a", dy=1.10)

    # ---- b: 各类群峰值基因数 ----
    ax_b = fig.add_subplot(gs[1, 0])
    counts = {sp: pd.read_csv(cfg.results("tables", sp, "specificity_global.tsv"), sep="\t")
              ["max_group"].value_counts() for sp in ("human", "mouse")}
    order = sorted(MAIN8, key=lambda g: -(counts["human"].get(g, 0) + counts["mouse"].get(g, 0)))
    y, h = np.arange(len(order)), 0.36
    ax_b.barh(y + h / 2 + 0.01, [counts["human"].get(g, 0) for g in order], height=h,
              color=C_HUMAN, edgecolor="white", linewidth=0.6, label="人")
    ax_b.barh(y - h / 2 - 0.01, [counts["mouse"].get(g, 0) for g in order], height=h,
              color=C_MOUSE, edgecolor="white", linewidth=0.6, label="鼠")
    for i, g in enumerate(order):
        ax_b.text(counts["human"].get(g, 0) + 0.4, i + h / 2 + 0.01,
                  str(counts["human"].get(g, 0)), va="center", fontsize=5.0, color=C_INK)
        ax_b.text(counts["mouse"].get(g, 0) + 0.4, i - h / 2 - 0.01,
                  str(counts["mouse"].get(g, 0)), va="center", fontsize=5.0, color=C_INK)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels([GROUP_CN[g] for g in order], fontsize=5.6)
    ax_b.invert_yaxis()
    ax_b.set_xlim(0, 31)
    clean(ax_b, grid="x")
    ax_b.tick_params(axis="y", length=0)
    ax_b.legend(fontsize=5.6, frameon=False, loc="lower right", ncol=2, handlelength=1.1)
    # 标题给结论; 计数规则用标题下一句大白话说明; 横轴只写"基因数"
    ax_b.set_xlabel("基因数", fontsize=6.0)
    ax_b.set_title("致病基因不集中于单一细胞类群", fontsize=6.8, color=C_INK, pad=12)
    ax_b.text(0.5, 1.02, "每个基因计入其相对表达最高的细胞类群", transform=ax_b.transAxes,
              ha="center", va="bottom", fontsize=5.4, color="#5a6470")
    label(ax_b, "b", dy=1.26)

    # ---- c: 具有明确发育变化模式的基因占比 ----
    ax_c = fig.add_subplot(gs[1, 1])
    rows = []
    for sp in ("human", "mouse"):
        sub = ma[(ma["species"] == sp) & (ma["status"] == "ok") & ma["in_panel"]
                 & ~ma["mech_class"].isin(["other", "unknown"]) & ma["shape"].notna()
                 & (ma["shape"] != "NA")]
        g = sub.groupby("mech_class")["shape"].apply(lambda s: (s != "flat").mean())
        rows += [{"species": sp, "mech": m, "prop": float(v)} for m, v in g.items()]
    prop = pd.DataFrame(rows)
    cls = [m for m in MECH_ORDER if m in set(prop["mech"])]
    x = np.arange(len(cls))
    w = 0.36
    hv = [float(prop[(prop.species == "human") & (prop.mech == m)]["prop"].iloc[0]) for m in cls]
    mv = [float(prop[(prop.species == "mouse") & (prop.mech == m)]["prop"].iloc[0]) for m in cls]
    ax_c.bar(x - w / 2 - 0.01, hv, width=w, color=C_HUMAN, edgecolor="white", linewidth=0.6, label="人")
    ax_c.bar(x + w / 2 + 0.01, mv, width=w, color=C_MOUSE, edgecolor="white", linewidth=0.6, label="鼠")
    for i, m in enumerate(cls):
        if m != "channel_synapse":
            continue
        ax_c.text(x[i] - w / 2 - 0.01, hv[i] + 0.02, f"{hv[i]:.2f}", ha="center",
                  fontsize=5.0, color=C_INK)
        ax_c.text(x[i] + w / 2 + 0.01, mv[i] + 0.02, f"{mv[i]:.2f}", ha="center",
                  fontsize=5.0, color=C_INK)
    ax_c.set_xticks(x)
    ax_c.set_xticklabels([MECH_CN[m] for m in cls], rotation=34, ha="right", fontsize=5.0)
    ax_c.set_ylim(0, 1.14)
    ax_c.set_ylabel("呈现明确发育变化模式\n的基因占比", fontsize=5.6)
    clean(ax_c, grid="y")
    ax_c.legend(fontsize=5.6, frameon=False, loc="upper left", ncol=2, handlelength=1.1)
    ax_c.set_title("离子通道与突触类基因几乎都随发育变化", fontsize=6.8, color=C_INK, pad=3)
    label(ax_c, "c", dy=1.14)

    fig.suptitle("图2　共性体现在发育时序而非细胞类型", fontsize=9, color=C_INK, y=0.985)
    save(fig, outdir, "report_fig2")


# =============================================================== 图3
def fig3(cfg, outdir):
    fig = plt.figure(figsize=(160 * MM, 122 * MM))
    gs = gridspec.GridSpec(3, 4, figure=fig, height_ratios=[30, 30, 40],
                           left=0.122, right=0.900, top=0.845, bottom=0.100,
                           hspace=0.85, wspace=0.26)
    bts = {sp: pd.read_csv(cfg.results("intermediate", sp, "cube_bin.tsv"), sep="\t")
           for sp in ("mouse", "human")}
    spc = {sp: cfg.effective_species(sp) for sp in ("mouse", "human")}
    subs = {"GRM1": {"mouse": "6个分析类群呈上升趋势；高峰：成年期",
                     "human": "5个分析类群呈上升趋势；高峰：3岁"},
            "ATXN3": {"mouse": "全小脑ρ=−0.76、颗粒谱系ρ=−0.70",
                      "human": "全小脑ρ=−0.70，与小鼠同向"}}
    first_ax = None
    for k, gene in enumerate(("GRM1", "ATXN3")):
        for j, sp in enumerate(("mouse", "human")):
            ax = fig.add_subplot(gs[k, j * 2:j * 2 + 2])
            first_ax = first_ax or ax
            bt, spcfg = bts[sp], spc[sp]
            for g, cn in TRAJ_LINES:
                col = f"log_{gene}"
                if col not in bt.columns:
                    continue
                dd = bt[(bt["group"] == g) & bt["feasible"]].sort_values("bin_order")
                if not len(dd):
                    continue
                xx = dd["bin_order"].to_numpy(float)
                yy = np.where(np.isfinite(dd[col].to_numpy(float)), dd[col].to_numpy(float), np.nan)
                gap = np.diff(xx) > 1
                xs = np.array_split(xx, np.where(gap)[0] + 1) if gap.any() else [xx]
                ys = np.array_split(yy, np.where(gap)[0] + 1) if gap.any() else [yy]
                for sx, sy in zip(xs, ys):
                    ax.plot(sx, sy, linestyle="-", linewidth=1.3, color=TRAJ_COLORS[g], zorder=3,
                            marker="o", markersize=2.4, markerfacecolor=TRAJ_COLORS[g],
                            markeredgecolor="white", markeredgewidth=0.3)
            ax.set_xticks(range(1, len(spcfg.labels) + 1))
            ax.set_xticklabels(spcfg.labels, rotation=90, fontsize=4.8)
            ax.set_xlim(0.4, len(spcfg.labels) + 0.6)
            ax.margins(y=0.10)
            if j == 0:
                ax.set_ylabel("平均表达\n(log1p CP10K)", fontsize=5.8)
            clean(ax)
            ax.set_title(f"{gene}（{'小鼠' if sp == 'mouse' else '人'}）", fontsize=6.8,
                         color=C_INK, pad=12)
            ax.text(0.5, 1.012, subs[gene][sp], transform=ax.transAxes, ha="center",
                    va="bottom", fontsize=5.0, color="#5a6470")
            label(ax, "abcd"[k * 2 + j], dy=1.16)

    # ---- e: 全部疾病的综合评分分布(不做 top/bottom 截选, 避免"挑极端"观感) ----
    fd = pd.read_csv(cfg.results("tables", "cross", "fidelity_by_disease.tsv"), sep="\t")
    fd = fd[fd["fidelity"].notna()].copy()
    fd = fd.sort_values("fidelity", ascending=False).reset_index(drop=True)
    fd["rank"] = np.arange(1, len(fd) + 1)
    ax_e = fig.add_subplot(gs[2, :])
    sc = ax_e.scatter(fd["rank"], fd["fidelity"], c=fd["fidelity"].clip(0, 1),
                      cmap=matplotlib.colormaps["viridis"], s=9, linewidths=0, zorder=3)
    med = float(fd["fidelity"].median())
    ax_e.axhline(med, color="#8a929b", lw=0.7, ls="--", zorder=2)
    ax_e.axhline(0.0, color="#8a929b", lw=0.7, ls=":", zorder=2)
    ax_e.text(len(fd) + 2, med, f"中位 {med:.2f}", fontsize=5.0, color="#5a6470", va="center")
    ax_e.text(len(fd) + 2, 0.0, "0", fontsize=5.0, color="#5a6470", va="center")
    # 只标注正文提及或两端各一的少数疾病作方位参照; 全部可评分疾病(fidelity 非 NA)均以点呈现
    OFF = {"AT": (12, 0.10), "SCA44": (0, 0.14), "SCA3": (0, -0.15),
           "SCA10": (-26, -0.08), "JBTS9": (-14, -0.22)}   # 底部两个标注向左偏移, 避免越界到色条
    for dz, (dx, dy) in OFF.items():
        r = fd[fd["disease_abbrev"] == dz]
        if not len(r):
            continue
        r = r.iloc[0]
        ax_e.annotate(f"{dz}（{str(r['genes']).split(',')[0]}）",
                      xy=(r["rank"], r["fidelity"]), xytext=(r["rank"] + dx, r["fidelity"] + dy),
                      fontsize=5.0, color=C_INK, ha="center", va="center",
                      arrowprops=dict(arrowstyle="-", lw=0.5, color="#8a929b"))
    ax_e.set_xlim(0, len(fd) + 16)
    ax_e.set_ylim(-0.45, 1.20)
    ax_e.set_xlabel(f"疾病按综合评分降序排列的序号（共{len(fd)}个疾病条目）", fontsize=6.2)
    ax_e.set_ylabel("综合评分", fontsize=6.2)
    clean(ax_e, grid="y")
    ax_e.set_title("全部疾病的跨物种发育表达一致性评分分布（颜色深浅同评分高低）",
                   fontsize=6.8, color=C_INK, pad=4)
    label(ax_e, "e", dy=1.12)
    cb = fig.colorbar(sc, ax=ax_e, pad=0.012, fraction=0.018)
    cb.set_label("综合评分（0–1 色标；负值见纵轴）", fontsize=4.8)
    cb.ax.tick_params(labelsize=4.6, length=1.5)
    cb.outline.set_linewidth(0.4)

    handles = [plt.Line2D([], [], color=TRAJ_COLORS[g], linewidth=1.3, label=cn)
               for g, cn in TRAJ_LINES]
    first_ax.legend(handles=handles, loc="upper left", fontsize=5.0, frameon=False,
                    handlelength=2.0, borderaxespad=0.4)
    fig.suptitle("图3　代表性致病基因的发育轨迹与疾病级跨物种一致性", fontsize=9, color=C_INK, y=0.985)
    save(fig, outdir, "report_fig3")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--config", default=None)
    ap.add_argument("--fig", default=None, choices=["1", "2", "3"])
    args = ap.parse_args()
    root = os.path.abspath(args.root)
    cfg = atxconfig.load(args.config or os.path.join(root, "config", "config.yaml"))
    outdir = cfg.results("figures", "report", mkdir=True)
    print(f"== a7 | version={cfg.version} | 字体解析: {FONT_RESOLVED}"
          + ("" if FONT_IS_YAHEI else
             "（本机无微软雅黑 → 回退 Noto Sans CJK SC；将 msyh.ttc/msyh.ttf 放入 fonts/ 后重跑即得雅黑版）"))
    for f, fn in (("1", fig1), ("2", fig2), ("3", fig3)):
        if args.fig in (None, f):
            fn(cfg, outdir)
    print(f"[a7] 完成 → {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
