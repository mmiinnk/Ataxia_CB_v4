"""atx.plots — 绘图原语(implementation_spec §3.12;figure_reading_guide 逐 panel 约定)。

口径(全部从 config figures 读, 硬约束 6;NA 政策硬约束 5):
- NA → na_style 斜纹(白底 ///), 绝不插值、绝不填 0;
- sem 带只画在"连续可行段"(bin_order 相邻且 logsem 有限)内; 缺口断线不桥接;
  单供体/pooled bin 的 logsem = NA(pseudobulk 硬断言)→ 该 bin 自然无带;
- b2 轨迹: donor 灰点、疾病供体开口 △、单供体 bin ▲(均见读图指南 Fig2d);
- b2 y 轴 = mean log1p CP10K: bin_table 存储的 sem 即 logsem_(log 尺度), 趋势统计
  rho_bin 亦基于 log_ 轨迹, 故读图指南"y=CP10K"落地为 log1p 尺度(唯一与 logsem_
  自洽的口径; 线性 CP10K 无对应 sem 列, 不得由 logsem_ 反推);
- 热图分页 heatmap_page_size, 色标全局(vmin/vmax 在分页前按整矩阵计算);
- 行序 row_order_rule=rho_desc: rho 降序, NaN 沉底、并列按基因名(确定性);
- footer 统一入口 add_footer(): draft 版本自动附 DRAFT 探索声明; R5 的
  "部分插值"标注(s07/s08 用 peak_rank_human_equiv 时)复用同一 notes 机制;
- 星号 star_levels(b3): q < 阈值取最严档 */**/***;
- 图内文字一律英文(环境无 CJK 字体)。

本模块只产 Figure/落盘, 不做统计; 所有数值来自 s02 六表 + s01 cube(调用方传入)。
"""
from __future__ import annotations

import math
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle


# ---------------------------------------------------------------- 通用小件
def safe_name(s: str) -> str:
    """文件/图注安全化(R9): 非 [0-9A-Za-z._-] 连续段 → '_'。'adult (9 weeks)'→'adult_9_weeks'。"""
    return re.sub(r"[^0-9A-Za-z._-]+", "_", str(s)).strip("_") or "NA"


def star_of_q(q, cfg) -> str:
    """q 值 → 星号(config figures.star_levels, 按阈值降序书写, 取满足的最严档)。"""
    try:
        q = float(q)
    except (TypeError, ValueError):
        return ""
    if not np.isfinite(q):
        return ""
    sym = ""
    for thr, s in cfg.raw["figures"]["star_levels"]:
        if q < float(thr):
            sym = str(s)
    return sym


def zscore_rows(M: np.ndarray) -> np.ndarray:
    """逐行 z(基因内 z, 读图指南 Fig2a/b); NaN 保留; 恒定行或仅 1 个可行 bin → 0。"""
    M = np.asarray(M, dtype=float)
    mu = np.nanmean(M, axis=1, keepdims=True)
    sd = np.nanstd(M, axis=1, keepdims=True)           # ddof=0(展示口径, 非统计推断)
    with np.errstate(invalid="ignore", divide="ignore"):
        Z = (M - mu) / sd
    degen = (~np.isfinite(sd)) | (sd == 0)             # 全 NaN / 恒定 / 单点行
    return np.where(degen & np.isfinite(M), 0.0, Z)


def contig_runs(orders, ok) -> list[tuple[int, int]]:
    """ok=True 且 orders 相邻(差 1)的极大连续段 → [(i0, i1)] 闭区间索引。

    sem 带/均值线只在段内画(缺口断线不桥接); 段长 1 → 画点不画带。
    """
    orders = np.asarray(orders, dtype=float)
    ok = np.asarray(ok, dtype=bool)
    runs: list[tuple[int, int]] = []
    i, n = 0, len(ok)
    while i < n:
        if not ok[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and ok[j + 1] and orders[j + 1] == orders[j] + 1:
            j += 1
        runs.append((i, j))
        i = j + 1
    return runs


def order_rho_desc(rho, names) -> list[str]:
    """rho 降序; NaN 沉底; 并列/NaN 内按名字典序(确定性; row_order_rule=rho_desc)。"""
    rho = np.asarray(rho, dtype=float)
    names = [str(x) for x in names]

    def key(i):
        r = rho[i]
        return (0, -r, names[i]) if np.isfinite(r) else (1, 0.0, names[i])

    return [names[i] for i in sorted(range(len(names)), key=key)]


def order_main_tau(sg: pd.DataFrame, main_groups: list[str]) -> list[str]:
    """c1/c2 行序(读图指南 Fig1c): 主组(max_group)按主组序 → 组内 tau 降序 → 基因名;
    max_group=NA 或 tau NA 沉底。sg 需含 gene/max_group/tau 列。"""
    pos = {g: i for i, g in enumerate(main_groups)}
    genes = sg["gene"].astype(str).tolist()
    mg = sg["max_group"].astype(str).tolist()
    tau = pd.to_numeric(sg["tau"]).to_numpy(dtype=float)

    def key(i):
        t = tau[i]
        fin = np.isfinite(t)
        return (pos.get(mg[i], len(pos)), 0 if fin else 1, -(t if fin else 0.0), genes[i])

    return [genes[i] for i in sorted(range(len(genes)), key=key)]


# ---------------------------------------------------------------- footer / 落盘
def footer_lines(cfg, notes=None) -> list[str]:
    """底部注文本行: 调用方 notes(如 R5 '部分插值') + draft 版本探索声明。"""
    lines = [str(s) for s in (notes or []) if str(s)]
    if getattr(cfg, "draft", False):
        lines.append(f"DRAFT gene_set v{cfg.version} — exploratory, not for citation")
    return lines


def add_footer(fig, cfg, notes=None) -> list[str]:
    lines = footer_lines(cfg, notes)
    if lines:
        fig.text(0.5, 0.002, " · ".join(lines), ha="center", va="bottom",
                 fontsize=6.5, color="#555555")
    return lines


def save_figures(fig, cfg, stem: str, outdir: str, notes=None) -> list[str]:
    """统一 footer → 按 config figures.formats 落盘(dpi/bbox_inches=tight)→ 关闭。"""
    add_footer(fig, cfg, notes)
    os.makedirs(outdir, exist_ok=True)
    fc = cfg.raw["figures"]
    paths = []
    for ext in fc["formats"]:
        p = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(p, dpi=int(fc["dpi"]), bbox_inches="tight")
        paths.append(p)
    plt.close(fig)
    return paths


# ---------------------------------------------------------------- 热图原语
def hatch_na(ax, mask, cfg, xoff: float = 0.5, yoff: float = 0.5):
    """NA 格 → config figures.na_style 斜纹(表 NA、图斜纹)。mask: bool 二维(与数据同形)。"""
    st = dict(cfg.raw["figures"]["na_style"])
    for i, j in zip(*np.where(np.asarray(mask, dtype=bool))):
        ax.add_patch(Rectangle((j - xoff, i - yoff), 1, 1, linewidth=0.4, **st))


def heatmap_figure(M, row_labels, col_labels, cfg, *, cmap_key="diverging_cmap",
                   vmin=None, vmax=None, title="", xrot=90, ylabel="", xlabel="",
                   cell_text=None, figsize=None, hatch=True):
    """单页 imshow 热图: masked NaN + NA 斜纹 + colorbar; cell_text 同形字符串矩阵
    (b3 星号)。返回 (fig, ax, im)。"""
    M = np.asarray(M, dtype=float)
    fc = cfg.raw["figures"]
    nr, nc = M.shape
    if figsize is None:
        cw = float(np.clip(13.0 / max(nc, 1), 0.10, 0.34))
        chh = float(np.clip(14.0 / max(nr, 1), 0.07, 0.20))
        figsize = (cw * nc + 2.6, chh * nr + 1.9)
    fig, ax = plt.subplots(figsize=figsize)
    cmap = plt.get_cmap(fc[cmap_key]).copy()
    cmap.set_bad(fc["na_style"]["facecolor"])
    if vmin is None or vmax is None:
        if np.isfinite(M).any():
            if cmap_key == "diverging_cmap":
                a = float(np.nanmax(np.abs(M)))
                vmin, vmax = -a, a
            else:
                vmin, vmax = float(np.nanmin(M)), float(np.nanmax(M))
        else:
            vmin, vmax = 0.0, 1.0
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
    if hatch:
        hatch_na(ax, ~np.isfinite(M), cfg)
    if cell_text is not None:
        ct = np.asarray(cell_text, dtype=object)
        for i, j in zip(*np.where(ct != "")):
            ax.text(j, i, str(ct[i, j]), ha="center", va="center", fontsize=5, color="black")
    ax.set_xticks(range(nc))
    ax.set_xticklabels([str(c) for c in col_labels], rotation=xrot,
                       ha="right" if xrot else "center", fontsize=max(4.5, min(7, 90 / max(nc, 1) * 0.7)))
    ax.set_yticks(range(nr))
    ax.set_yticklabels([str(r) for r in row_labels], fontsize=6)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=8)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=8)
    ax.set_title(title, fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.85)
    return fig, ax, im


def paged_heatmap(M, row_labels, col_labels, cfg, *, page_axis="rows", cell_text=None,
                  vmin=None, vmax=None, cmap_key="diverging_cmap", title="", xrot=90,
                  ylabel="", xlabel="", header_draw=None):
    """分页热图(heatmap_page_size; 色标全局 = vmin/vmax 按整矩阵先算)。

    返回 [(fig, ax, page_slice)]。header_draw(ax, sl) 供 c2 期色带等页级装饰。
    """
    M = np.asarray(M, dtype=float)
    fc = cfg.raw["figures"]
    ps = int(fc["heatmap_page_size"])
    n = M.shape[0] if page_axis == "rows" else M.shape[1]
    npages = max(1, math.ceil(n / ps))
    if vmin is None or vmax is None:
        if np.isfinite(M).any():
            if cmap_key == "diverging_cmap":
                a = float(np.nanmax(np.abs(M)))
                vmin, vmax = -a, a
            else:
                vmin, vmax = float(np.nanmin(M)), float(np.nanmax(M))
        else:
            vmin, vmax = 0.0, 1.0
    row_labels = np.asarray([str(x) for x in row_labels], dtype=object)
    col_labels = np.asarray([str(x) for x in col_labels], dtype=object)
    out = []
    for pi in range(npages):
        sl = slice(pi * ps, min((pi + 1) * ps, n))
        if page_axis == "rows":
            Mp, rp, cp = M[sl], row_labels[sl], col_labels
            tp = cell_text[sl] if cell_text is not None else None
        else:
            Mp, rp, cp = M[:, sl], row_labels, col_labels[sl]
            tp = cell_text[:, sl] if cell_text is not None else None
        t = title if npages == 1 else f"{title} — page {pi + 1}/{npages}"
        fig, ax, _ = heatmap_figure(Mp, rp, cp, cfg, cmap_key=cmap_key, vmin=vmin, vmax=vmax,
                                    title=t, xrot=xrot, ylabel=ylabel, xlabel=xlabel, cell_text=tp)
        if header_draw is not None:
            header_draw(ax, sl)
        out.append((fig, ax, sl))
    return out


# ---------------------------------------------------------------- 矩阵取数
def group_bin_matrix(bt: pd.DataFrame, group: str, genes: list[str], labels: list[str],
                     prefix: str) -> np.ndarray:
    """bin_table 某组 → genes × labels 矩阵(prefix=log_/cp_/pct_/logsem_);
    缺列基因 → 全 NaN 行(status=missing 路径); 行序 = labels(发育序)。"""
    sub = bt.loc[bt["group"].to_numpy() == group].set_index("bin").reindex(labels)
    M = np.full((len(genes), len(labels)), np.nan)
    for i, g in enumerate(genes):
        c = f"{prefix}{g}"
        if c in sub.columns:
            M[i] = pd.to_numeric(sub[c]).to_numpy(dtype=float)
    return M


# ---------------------------------------------------------------- b2 轨迹
def traj_panel(ax, bt, cube, gene, group, cfg, spcfg, rng=None):
    """单 panel 轨迹(读图指南 Fig2d)。返回计数 dict(供合成测试断言)。

    灰点 = donor pseudobulk(cube X, 疾病供体开口黑 △); 实线 = 组均值(log_ 列,
    只在连续可行段内); 带 = ±logsem(仅连续可行且 sem 有限段, 段长 ≥2);
    单供体可行 bin ▲; 不可行/NA bin 背景斜纹。
    """
    fc = cfg.raw["figures"]
    color = fc["group_colors"].get(group, "#333333")
    labels = list(spcfg.labels)
    x = np.array([float(spcfg.bin_order[b]) for b in labels])
    sub = bt.loc[bt["group"].to_numpy() == group].set_index("bin").reindex(labels)
    lcol, scol = f"log_{gene}", f"logsem_{gene}"
    y = (pd.to_numeric(sub[lcol]).to_numpy(dtype=float) if lcol in sub.columns
         else np.full(len(labels), np.nan))
    sem = (pd.to_numeric(sub[scol]).to_numpy(dtype=float) if scol in sub.columns
           else np.full(len(labels), np.nan))
    feas = sub["feasible"].fillna(False).astype(bool).to_numpy()
    single = sub["single_donor"].fillna(False).astype(bool).to_numpy()
    ok = feas & np.isfinite(y)

    # donor 灰点 / 疾病 △(cube 可 None → 无点)
    dx = dy = np.array([])
    ddis = np.array([], dtype=bool)
    if cube is not None and gene in set(map(str, cube.var["gene_key"])):
        o = cube.obs
        gi = np.where(o["group"].to_numpy().astype(str) == group)[0]
        if len(gi):
            jgene = list(map(str, cube.var["gene_key"])).index(gene)
            Xg = np.asarray(cube.X[gi][:, jgene].todense()).ravel() if hasattr(cube.X, "todense") \
                else np.asarray(cube.X[gi][:, jgene]).ravel()
            ddis = (o["disease_flag"].to_numpy()[gi].astype(bool)
                    if "disease_flag" in o.columns else np.zeros(len(gi), dtype=bool))
            fin = np.isfinite(Xg)
            if rng is None:
                rng = np.random.default_rng(int(cfg.seed))
            jit = rng.uniform(-0.14, 0.14, size=int(fin.sum()))
            dx = o["bin_order"].to_numpy()[gi][fin].astype(float) + jit
            dy = Xg[fin]
            ddis = ddis[fin]
            ax.scatter(dx[~ddis], dy[~ddis], s=4, color="#9e9e9e", alpha=0.55, lw=0, zorder=2)
            if ddis.any():
                ax.scatter(dx[ddis], dy[ddis], s=16, marker="^", facecolors="none",
                           edgecolors="black", linewidths=0.8, zorder=4)

    # 均值线(连续可行段内; 段长 1 → 只画点)
    n_seg = 0
    for i0, i1 in contig_runs(x, ok):
        ax.plot(x[i0:i1 + 1], y[i0:i1 + 1], "-o", color=color, lw=1.2, ms=2.6, zorder=3)
        n_seg += 1
    # sem 带: 可行且 sem 有限的连续段(段长 ≥2)
    oksem = ok & np.isfinite(sem)
    n_band = 0
    for i0, i1 in contig_runs(x, oksem):
        if i1 > i0:
            ax.fill_between(x[i0:i1 + 1], (y - sem)[i0:i1 + 1], (y + sem)[i0:i1 + 1],
                            color=color, alpha=0.22, lw=0, zorder=1)
            n_band += 1
    # 单供体 ▲(开口三角, 与均值线同色)
    smask = ok & single
    if smask.any():
        ax.plot(x[smask], y[smask], "^", color=color, ms=5.5, mfc="none", mew=1.0, zorder=5)

    # y 范围(数据 + donor 点)→ 不可行 bin 背景斜纹
    allv = np.concatenate([y[ok], dy]) if (ok.any() or len(dy)) else np.array([0.0, 1.0])
    allv = allv[np.isfinite(allv)]
    if len(allv) == 0:
        allv = np.array([0.0, 1.0])
    lo, hi = float(allv.min()), float(allv.max())
    pad = 0.08 * (hi - lo) if hi > lo else 0.5
    ax.set_ylim(lo - pad, hi + pad)
    st = dict(fc["na_style"])
    n_hatch = int((~ok).sum())
    for xi in x[~ok]:
        ax.add_patch(Rectangle((xi - 0.5, lo - pad), 1, 2 * pad + (hi - lo), linewidth=0.4,
                               transform=ax.transData, **st))
    ax.set_xlim(x.min() - 0.6, x.max() + 0.6)
    ax.set_ylabel("mean log1p CP10K", fontsize=6.5)
    return {"n_seg": n_seg, "n_band": n_band, "n_single": int(smask.sum()),
            "n_donor_pts": int(len(dx)), "n_disease_pts": int(ddis.sum()), "n_hatch": n_hatch}


def traj_figure(bt, cube, genes, group, cfg, spcfg, *, title=None, ncols=3):
    """一页轨迹小倍数(len(genes) ≤ traj_panels_per_page)。返回 (fig, [panel 计数])。"""
    fc = cfg.raw["figures"]
    n = len(genes)
    nrows = max(1, math.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 2.35 * nrows + 0.5),
                             squeeze=False)
    rng = np.random.default_rng(int(cfg.seed))
    infos = []
    for k, gene in enumerate(genes):
        ax = axes[k // ncols][k % ncols]
        infos.append(traj_panel(ax, bt, cube, gene, group, cfg, spcfg, rng=rng))
        ax.set_title(str(gene), fontsize=8)
        ax.tick_params(labelsize=6)
    for k in range(n, nrows * ncols):
        axes[k // ncols][k % ncols].axis("off")
    xs = [float(spcfg.bin_order[b]) for b in spcfg.labels]
    for k in range(nrows * ncols):
        ax = axes[k // ncols][k % ncols]
        ax.set_xticks(xs)
        ax.set_xticklabels(list(spcfg.labels) if k // ncols == nrows - 1 else [],
                           rotation=90, fontsize=5.5)
    fig.suptitle(title or f"{group} · trajectories (gray = donor, band = ±sem, "
                          f"▲ = single-donor bin, open △ = disease donor, hatch = infeasible)",
                 fontsize=9.5)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    return fig, infos


# ---------------------------------------------------------------- c1 dotplot
def dotplot_figure(bt, sg, gene_order, groups, cfg, *, title=None, species=""):
    """9 组特异性 dotplot(读图指南 Fig1c/d): 行=基因(主组→tau 序), 列=组;
    点大小 = 可行 bin 的 pct_ 均值(表达核比例, 0–1 分数口径; × dot_size_max),
    色 = 可行 bin 的 log_(mean log1p CP10K)均值;
    NA → 无点(不插值); 右缘列 = global tau(tau_cmap, NA 斜纹)。
    返回 (fig, info)。"""
    fc = cfg.raw["figures"]
    ng, ngrp = len(gene_order), len(groups)
    # 组×(基因) 汇总: pct/size 与 log/color 均只在可行 bin 上求均值(nanmean=bin 等权)
    Size = np.full((ng, ngrp), np.nan)
    Color = np.full((ng, ngrp), np.nan)
    for j, g in enumerate(groups):
        sub = bt.loc[bt["group"].to_numpy() == g]
        feas = sub["feasible"].fillna(False).astype(bool).to_numpy()
        sub = sub.loc[feas]
        for i, gene in enumerate(gene_order):
            pc, lc = f"pct_{gene}", f"log_{gene}"
            if pc in sub.columns:
                v = pd.to_numeric(sub[pc]).to_numpy(dtype=float)
                if np.isfinite(v).any():
                    Size[i, j] = float(np.nanmean(v))
            if lc in sub.columns:
                v = pd.to_numeric(sub[lc]).to_numpy(dtype=float)
                if np.isfinite(v).any():
                    Color[i, j] = float(np.nanmean(v))
    s_max = float(fc.get("dot_size_max", 140))
    if np.isfinite(Color).any():
        vmin, vmax = float(np.nanmin(Color)), float(np.nanmax(Color))
    else:
        vmin, vmax = 0.0, 1.0
    fig = plt.figure(figsize=(0.62 * ngrp + 4.6, max(3.2, 0.135 * ng + 2.0)))
    gs = fig.add_gridspec(1, 2, width_ratios=[max(ngrp, 1), 1], wspace=0.08)
    ax = fig.add_subplot(gs[0])
    axT = fig.add_subplot(gs[1])
    n_dots = 0
    cmap = plt.get_cmap(fc["sequential_cmap"])
    for j in range(ngrp):
        m = np.isfinite(Size[:, j]) & np.isfinite(Color[:, j])
        if m.any():
            ax.scatter(np.full(int(m.sum()), float(j)), np.arange(ng)[m],
                       s=Size[m, j] * s_max, c=Color[m, j], cmap=cmap,
                       vmin=vmin, vmax=vmax, edgecolors="#bbbbbb", linewidths=0.3, zorder=3)
            n_dots += int(m.sum())
    ax.set_xlim(-0.6, ngrp - 0.4)
    ax.set_ylim(ng - 0.5, -0.5)
    ax.set_xticks(range(ngrp))
    ax.set_xticklabels([str(g) for g in groups], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(ng))
    ax.set_yticklabels([str(g) for g in gene_order], fontsize=5)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#e8e8e8", lw=0.5, zorder=0)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin, vmax))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, shrink=0.75, label="mean log1p CP10K (feasible bins)")
    # tau 边列
    tau = pd.to_numeric(sg.set_index("gene").reindex(gene_order)["tau"]).to_numpy(dtype=float)
    tcmap = plt.get_cmap(fc["tau_cmap"]).copy()
    tcmap.set_bad(fc["na_style"]["facecolor"])
    im = axT.imshow(np.ma.masked_invalid(tau.reshape(-1, 1)), cmap=tcmap, aspect="auto",
                    vmin=0.0, vmax=1.0)
    hatch_na(axT, ~np.isfinite(tau).reshape(-1, 1), cfg)
    axT.set_xticks([0], ["tau"], fontsize=7)
    axT.set_yticks([])
    fig.colorbar(im, ax=axT, shrink=0.75, label="tau (global)")
    ax.set_title(title or f"{species} c1 · specificity dotplot "
                          f"(size = % nuclei detected, rows: main group → tau)", fontsize=9)
    return fig, {"n_dots": n_dots, "n_genes": ng, "n_groups": ngrp}


# ---------------------------------------------------------------- b3 rho 热图
def rho_heatmap_pages(tr, gene_order, groups, cfg, *, title=""):
    """组×基因 rho 热图 + 星号(读图指南 Fig2f): 色=rho_bin(±1 固定), 白底斜纹=rho NA,
    星=q_bin 最严档。分页同热图原语。返回 [(fig, ax, sl)]。"""
    tri = tr.set_index(["gene", "group"])
    M = np.full((len(gene_order), len(groups)), np.nan)
    T = np.full(M.shape, "", dtype=object)
    for i, gene in enumerate(gene_order):
        for j, g in enumerate(groups):
            if (gene, g) in tri.index:
                r = tri.loc[(gene, g)]
                M[i, j] = float(pd.to_numeric(r["rho_bin"]))
                T[i, j] = star_of_q(r["q_bin"], cfg)
    return paged_heatmap(M, gene_order, groups, cfg, vmin=-1.0, vmax=1.0,
                         cell_text=T, title=title, xrot=45, ylabel="gene", xlabel="group")


# ---------------------------------------------------------------- b4 pattern 堆叠条
def pattern_bars_figure(tp, groups, cfg, *, title=""):
    """每组 temporal_pattern 计数堆叠条(读图指南 Fig2e 左; 词表/色 = config pattern_colors)。"""
    fc = cfg.raw["figures"]
    vocab = [str(k) for k in fc["pattern_colors"].keys()]
    d = tp.assign(_n=1).pivot_table(index="group", columns="temporal_pattern",
                                    values="_n", aggfunc="sum")
    d = d.reindex(index=list(groups), columns=vocab).fillna(0).astype(int)
    fig, ax = plt.subplots(figsize=(7.2, 0.44 * len(groups) + 2.4))
    ypos = np.arange(len(groups))
    left = np.zeros(len(groups))
    handles = []
    for pat in vocab:
        w = d[pat].to_numpy()
        ax.barh(ypos, w, left=left, color=fc["pattern_colors"][pat], edgecolor="black",
                linewidth=0.3, height=0.72)
        handles.append(Patch(facecolor=fc["pattern_colors"][pat], edgecolor="black", label=pat))
        left += w
    ax.set_yticks(ypos, [str(g) for g in groups], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("genes", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.legend(handles=handles, fontsize=6.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------- c2 期×组占比热图
def period_fraction_pages(bp, gene_order, cfg, *, title="", species=""):
    """发育期×组占比热图(读图指南 Fig1e): 行=基因(c1 同序), 列=期×主组,
    值=fraction_<g>(specificity_by_period), NA 斜纹; 顶部期色带 + 该期 n_groups 列头。"""
    fc = cfg.raw["figures"]
    main = cfg.main_groups()
    periods = list(dict.fromkeys(bp["period"].astype(str)))
    cols = [(p, g) for p in periods for g in main]
    idx = bp.set_index(["gene", "period"])
    M = np.full((len(gene_order), len(cols)), np.nan)
    for i, gene in enumerate(gene_order):
        for j, (p, g) in enumerate(cols):
            fc_col = f"fraction_{g}"
            if (gene, p) in idx.index and fc_col in idx.columns:
                M[i, j] = float(pd.to_numeric(idx.loc[(gene, p), fc_col]))
    # 每期 n_groups(列头): 该期基因中 n_groups_used 最大值(缺失基因为 0/NA)
    n_by_per = {}
    for p in periods:
        sub = bp.loc[bp["period"].astype(str) == p, "n_groups_used"]
        v = pd.to_numeric(sub).to_numpy(dtype=float)
        n_by_per[p] = int(np.nanmax(v)) if np.isfinite(v).any() else "NA"
    col_labels = [g for _, g in cols]
    col_periods = [p for p, _ in cols]

    def strip(ax, sl):
        x_colors = fc["period_colors"]
        for p in periods:
            js = [j for j in range(len(cols)) if col_periods[j] == p]
            if not js:
                continue
            a, b = js[0] - 0.5, js[-1] + 0.5
            c = x_colors.get(p, "#999999")
            ax.add_patch(Rectangle((a, -1.75), b - a, 1.05, color=c, clip_on=False, lw=0))
            rgb = matplotlib.colors.to_rgb(c)
            lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
            ax.text((a + b) / 2, -1.22, f"{p}\nn={n_by_per[p]}", ha="center", va="center",
                    fontsize=6, color="black" if lum > 0.55 else "white", clip_on=False)
        ax.set_xlim(-0.5, len(cols) - 0.5)

    pages = paged_heatmap(M, gene_order, col_labels, cfg, cmap_key="sequential_cmap",
                          vmin=0.0, vmax=1.0, title=title, xrot=90, ylabel="gene",
                          header_draw=strip)
    return pages


# ---------------------------------------------------------------- d1 emergence 时间轴
def emergence_figure(em, cfg, spcfg, *, title="", species=""):
    """特异性建立时间轴(读图指南 Fig1f): 行=基因按 emergence 先后, 点=首个主组
    fraction≥阈值 bin(色=主组), 线回零点; 左截断=开口点; 无点行=从未达到。"""
    fc = cfg.raw["figures"]
    gc = fc["group_colors"]
    e = em.copy()
    e["gene"] = e["gene"].astype(str)
    has = e["emergence_bin"].astype(str) != "NA"
    e1 = e[has].sort_values(["emergence_bin_order", "gene"], kind="stable")
    e0 = e[~has].sort_values("gene", kind="stable")
    order = pd.concat([e1, e0], ignore_index=True)
    n = len(order)
    xs = np.array([float(spcfg.bin_order[b]) for b in spcfg.labels])
    x_left = float(xs.min()) - 0.5
    fig, ax = plt.subplots(figsize=(7.6, max(3.0, 0.155 * n + 1.8)))
    seen: dict[str, int] = {}
    for i, r in order.iterrows():
        mg = str(r["max_group"])
        if str(r["emergence_bin"]) == "NA":
            continue
        c = gc.get(mg, "#999999")
        xo = float(r["emergence_bin_order"])
        ax.hlines(i, x_left, xo, color=c, lw=0.9, alpha=0.75, zorder=2)
        lc = bool(r.get("left_censored", False))
        ax.scatter([xo], [i], s=13, marker="o",
                   facecolors="none" if lc else c, edgecolors=c, linewidths=0.9, zorder=3)
        seen[mg] = seen.get(mg, 0) + 1
    ax.set_yticks(range(n), order["gene"].tolist(), fontsize=5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks(xs)
    ax.set_xticklabels([str(b) for b in spcfg.labels], rotation=90, fontsize=6)
    ax.set_xlim(x_left - 0.15, float(xs.max()) + 0.5)
    ax.grid(axis="x", color="#eeeeee", lw=0.5, zorder=0)
    handles = [Patch(facecolor=gc.get(g, "#999999"), edgecolor="black", label=f"{g} ({k})")
               for g, k in sorted(seen.items(), key=lambda kv: -kv[1])]
    handles.append(plt.Line2D([], [], marker="o", ls="", mfc="none", mec="#555555",
                              label="left-censored", ms=5))
    ax.legend(handles=handles, fontsize=6, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.set_title(title or f"{species} d1 · emergence timeline "
                          f"(dot = first bin with main-group fraction ≥ threshold)", fontsize=9)
    fig.tight_layout()
    return fig, {"n_with_emergence": int(has.sum()), "n_rows": n}
