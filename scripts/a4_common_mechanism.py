#!/usr/bin/env python3
"""a4_common_mechanism.py — 共性致病机制筛选(鼠侧, 一次性分析层)。

定位: **不属 sNN 契约层**(不改 s00–s11 的任何产物与口径), 是 A4 阶段任务的独立分析脚本,
为「异质性小脑共济失调的共性致病机制」提供第一块可写入阶段报告的证据。产物写在
results/v{ver}/analysis/a4_common_mechanism_{sp}/ 下, 与流水线产物目录隔离。

做什么(四块):
  1. 时空共识聚类  → 把致病基因按「细胞类型 × 发育时间」表达轨迹聚成共享时序模块
     (consensus clustering: 反复子采样基因/列 → 共分配频率 → 1-consensus 树上切 k;
      k 由 PAC 最小者选定, 全过程 seed=cfg.seed 可复现)。
  2. 模块刻画      → 每模块的 eigen-profile、主导细胞类型、峰期/方向、hub 基因
     (hub = 模块内平均 |corr| 最高者, 即"共性程序"的代表节点)。
  3. 机制类富集    → mech_class × module / mech_class × shape 置换富集。
     **按 CL-027 排除残差兜底类 other/unknown**(它们不是一种机制; atx.stats_enrich
     另有硬保护, 含残差类即 ValueError), 被排除基因逐条打印并写入报告。
  4. 跨组一致性    → 由 trend_by_group.tsv 统计每基因在 8 个主组中显著趋势的方向一致性
     (pan_up / pan_down / mixed / group_specific / ns), 即"共性 vs 异质性"的直接量化。

硬口径(与流水线一致, 不得违反):
  · NA 政策(硬约束 5): 不可行 组×bin 一律 NA — **不插值、不填 0**; 相关性用成对完整观测
    (pairwise-complete, min_periods), 零方差/未表达基因不参与聚类并单独列出。
  · 只用 group_tiers.main 的 8 个主组(UBC=descriptive、ALL=reference 不入聚类)。
  · 单供体 bin 保留(与 s02 同口径)但计数报告; --exclude-single-donor 可做敏感性。
  · 结论属**假设生成**(无疾病数据、单物种), 不得写成因果结论。

用法: /home/victor/.conda/envs/omics/bin/python -I scripts/a4_common_mechanism.py <root>
      [--species mouse] [--n-reps 200] [--k-max 8] [--n-perm 10000] [--exclude-single-donor]
"""
from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atx import config as atxconfig                     # noqa: E402
from atx import stats_enrich as atxenrich               # noqa: E402

MIN_PERIODS = 20        # 成对相关的最小共同观测列数(不足 → 该对相关记 NA)
PAC_LO, PAC_HI = 0.1, 0.9   # PAC 判定区间(标准 ConsensusClusterPlus 口径)
CONS_SUB = 0.8          # 每次子采样比例(基因与列各 80%)
HUB_TOPN = 5


# --------------------------------------------------------------------------- 工具
def _link(dist: np.ndarray) -> np.ndarray:
    d = np.nan_to_num(dist, nan=1.0)          # 相关缺失 → 距离 1(中性), 计数另行报告
    d = np.clip((d + d.T) / 2.0, 0, 2)
    np.fill_diagonal(d, 0.0)
    return linkage(squareform(d, checks=False), method="average")


def _pac(cons: np.ndarray, labels: np.ndarray) -> float:
    """PAC = 共识值落在 (0.1, 0.9) 的基因对比例(越低越稳定)。"""
    iu = np.triu_indices(len(labels), k=1)
    same = (labels[iu[0]] == labels[iu[1]])
    amb = ~same & (cons[iu] > PAC_LO)          # 不同簇但共识中高 → 模糊
    amb |= same & (cons[iu] < PAC_HI)          # 同簇但共识低 → 模糊
    return float(amb.mean())


def _md(df: pd.DataFrame, floatfmt: str | None = None) -> str:
    """DataFrame → GitHub 管道表(不依赖 tabulate; 环境未装 to_markdown 后端)。"""
    d = df.copy()
    if floatfmt:
        for c in d.columns:
            if pd.api.types.is_numeric_dtype(d[c]):
                d[c] = d[c].map(lambda v: "" if pd.isna(v) else format(float(v), floatfmt))
    cols = [str(c) for c in d.columns]
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    # 不用 iterrows(): 它会把混合 int/float 的行上转 float64(k=3 → "3.0")
    rows = d.astype(object).where(d.notna(), None).to_numpy()
    for r in rows:
        out.append("| " + " | ".join("" if v is None else str(v).replace("|", "/") for v in r) + " |")
    return "\n".join(out)


def consensus(Z: pd.DataFrame, k_cons: int, n_reps: int, seed: int):
    """共识矩阵: 反复子采样基因×列 → 层次聚类 → 共分配频率。"""
    rng = np.random.default_rng(seed)
    genes = list(Z.index)
    n = len(genes)
    ssum = np.zeros((n, n))
    scnt = np.zeros((n, n))
    ncol = Z.shape[1]
    for _ in range(n_reps):
        gi = rng.choice(n, size=max(3, int(round(n * CONS_SUB))), replace=False)
        ci = rng.choice(ncol, size=max(5, int(round(ncol * CONS_SUB))), replace=False)
        sub = Z.iloc[gi, ci]
        c = sub.T.corr(min_periods=max(5, MIN_PERIODS // 2)).to_numpy()   # 基因×基因(成对完整)
        lab = fcluster(_link(1.0 - c), k_cons, criterion="maxclust")
        same = (lab[:, None] == lab[None, :]).astype(float)
        ssum[np.ix_(gi, gi)] += same
        scnt[np.ix_(gi, gi)] += 1.0
    with np.errstate(invalid="ignore", divide="ignore"):
        cons = np.where(scnt > 0, ssum / scnt, np.nan)
    np.fill_diagonal(cons, 1.0)
    return pd.DataFrame(cons, index=genes, columns=genes), int((scnt == 0).sum() // 2)


# --------------------------------------------------------------------------- 主流程
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--config", default=None)
    ap.add_argument("--species", default="mouse")
    ap.add_argument("--n-reps", type=int, default=200)
    ap.add_argument("--k-cons", type=int, default=8, help="共识阶段的过指定簇数")
    ap.add_argument("--k-min", type=int, default=3)
    ap.add_argument("--k-max", type=int, default=8)
    ap.add_argument("--n-perm", type=int, default=None, help="默认取 config stats.enrichment.n_perm")
    ap.add_argument("--exclude-single-donor", action="store_true",
                    help="敏感性: 剔除单供体 bin(默认保留, 与 s02 同口径)")
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    cfg = atxconfig.load(args.config or os.path.join(root, "config", "config.yaml"))
    sp = args.species
    ecfg = cfg.stats_cfg.get("enrichment", {}) or {}
    n_perm = int(args.n_perm or ecfg.get("n_perm", 10000))
    residual = list(ecfg.get("exclude_mech_classes") or atxenrich.RESIDUAL_MECH_CLASSES)
    outdir = cfg.results("analysis", f"a4_common_mechanism_{sp}", mkdir=True)
    main_groups = list(cfg.main_groups())
    t0 = datetime.now()

    # ---- 1. 载入 cube 与注释 ----
    inter = cfg.results("intermediate", sp)
    bt = pd.read_csv(os.path.join(inter, "cube_bin.tsv"), sep="\t")
    genes_tab = pd.read_csv(cfg.paths["genes_table"], sep="\t").set_index("hgnc_symbol")
    # in_panel 在 tsv 里是小写 true/false → pandas 会解析成 bool; 统一转字符串比较(防 dtype 漂移)
    in_panel = genes_tab["in_panel"].astype(str).str.lower().eq("true")
    panel = genes_tab[(genes_tab["tier"] == "core") & in_panel]

    bt = bt[bt["feasible"] & bt["group"].isin(main_groups)].copy()
    if args.exclude_single_donor:
        bt = bt[~bt["single_donor"].astype(bool)]
    bt["bin_order"] = bt["bin_order"].astype(int)
    bt = bt.sort_values(["group", "bin_order"])
    bt["col"] = bt["group"] + "|" + bt["bin"].astype(str)
    gcols = [c for c in bt.columns if c.startswith("log_")]
    M = bt.set_index("col")[[c for c in gcols]].T            # 基因 × (组|bin)
    M.index = [i[4:] for i in M.index]
    M = M.loc[~M.index.duplicated()]

    stat_genes = [g for g in panel.index if g in M.index]
    absent = sorted(set(panel.index) - set(M.index))
    X = M.loc[stat_genes].astype(float)

    # ---- 2. 可聚类基因筛选(NA 政策: 不插值不填 0) ----
    spec = pd.read_csv(cfg.results("tables", sp, "specificity_global.tsv"), sep="\t").set_index("gene")
    not_expr = sorted({g for g in spec.index if spec.loc[g, "specific_call"] == "not_expressed"} & set(X.index))
    n_obs = X.notna().sum(axis=1)
    sd = X.std(axis=1, skipna=True)
    zerovar = sorted(X.index[(sd.fillna(0) <= 1e-12) & ~X.index.isin(not_expr)])
    lowobs = sorted(X.index[(n_obs < MIN_PERIODS) & ~X.index.isin(not_expr) & ~X.index.isin(zerovar)])
    clust_genes = [g for g in X.index if g not in set(not_expr) | set(zerovar) | set(lowobs)]
    Z = X.loc[clust_genes]
    Z = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1).replace(0, np.nan), axis=0)   # 逐基因 z

    n_sd_bins = int(bt["single_donor"].astype(bool).sum())
    n_pool_bins = int(bt["pooled"].astype(bool).sum())
    print(f"[a4] {sp}: 统计 panel {len(panel)} 基因, cube 命中 {len(stat_genes)}, 缺失 {len(absent)}{absent}")
    print(f"[a4] 聚类列 = {X.shape[1]} 个可行 组×bin(8 主组; 单供体 {n_sd_bins}, pooled {n_pool_bins}"
          f"{' — 已剔除单供体' if args.exclude_single_donor else ''})")
    print(f"[a4] 排除: not_expressed {len(not_expr)}{not_expr} / 零方差 {len(zerovar)}{zerovar}"
          f" / 观测不足 {len(lowobs)}{lowobs} → 入聚类 {len(clust_genes)}")

    # ---- 3. 共识聚类 + k 选择(PAC) ----
    C_obs = Z.T.corr(min_periods=MIN_PERIODS)      # 基因×基因 相关(成对完整观测, NA 不填充)
    n_nan_pairs = int(C_obs.isna().to_numpy().sum() // 2)
    cons, n_unsampled = consensus(Z, args.k_cons, args.n_reps, cfg.seed)
    cons_f = cons.fillna(cons.mean().mean()).to_numpy()
    L = _link(1.0 - cons_f)
    ks = list(range(args.k_min, args.k_max + 1))
    pac_rows, labels_by_k = [], {}
    for k in ks:
        lab = fcluster(L, k, criterion="maxclust")
        labels_by_k[k] = lab
        pac_rows.append({"k": k, "PAC": round(_pac(cons_f, lab), 4),
                         "min_module_size": int(pd.Series(lab).value_counts().min()),
                         "mean_within_consensus": round(float(np.mean(
                             [cons_f[np.ix_(np.where(lab == c)[0], np.where(lab == c)[0])].mean()
                              for c in sorted(set(lab)) if (lab == c).sum() > 1])), 4)})
    pac_tab = pd.DataFrame(pac_rows)
    pac_tab["k"] = pac_tab["k"].astype(int)
    pac_tab["min_module_size"] = pac_tab["min_module_size"].astype(int)
    k_best = int(pac_tab.sort_values(["PAC", "k"]).iloc[0]["k"])
    labels = labels_by_k[k_best]
    mod = pd.Series(labels, index=clust_genes).astype(str)
    mod.name = "module"
    print(f"[a4] 共识聚类: n_reps={args.n_reps} k_cons={args.k_cons} seed={cfg.seed}; "
          f"未同时子采样到的基因对 {n_unsampled}; 相关缺失对 {n_nan_pairs}")
    print("[a4] PAC 选 k:\n" + pac_tab.to_string(index=False))
    print(f"[a4] 选定 k={k_best}(PAC 最小)")

    # ---- 4. 模块刻画 + hub ----
    corr_f = C_obs.fillna(0.0).to_numpy()
    gene_ix = {g: i for i, g in enumerate(clust_genes)}
    rows = []
    for g in clust_genes:
        m = mod[g]
        peers = [p for p in clust_genes if mod[p] == m and p != g]
        conn = float(np.mean([abs(corr_f[gene_ix[g], gene_ix[p]]) for p in peers])) if peers else np.nan
        rows.append({"gene": g, "module": m, "connectivity": conn,
                     "mech_class": genes_tab.loc[g, "mech_class"] if g in genes_tab.index else "NA",
                     "mech_class_alt": genes_tab.loc[g, "mech_class_alt"] if g in genes_tab.index else "NA",
                     "tier": genes_tab.loc[g, "tier"] if g in genes_tab.index else "NA",
                     "max_group": spec.loc[g, "max_group"] if g in spec.index else "NA",
                     "tau": spec.loc[g, "tau"] if g in spec.index else np.nan,
                     "specific_call": spec.loc[g, "specific_call"] if g in spec.index else "NA",
                     "n_obs_columns": int(n_obs[g])})
    mg = pd.DataFrame(rows)
    # 模块内按 connectivity 降序定 hub 名次; 单基因模块 connectivity=NaN(无同模块邻居)→ 排最后, rank 仍为 1
    mg = mg.sort_values(["module", "connectivity"], ascending=[True, False], na_position="last")
    mg["hub_rank"] = mg.groupby("module").cumcount() + 1
    mg["is_hub"] = mg["hub_rank"] <= HUB_TOPN
    mg = mg.sort_values(["module", "hub_rank"])

    # eigen-profile(模块内 z 均值; 成对完整, 不填 0)
    eig = Z.groupby(mod.to_numpy()).mean()
    eig.index = [str(i) for i in eig.index]
    colmeta = bt.set_index("col")[["group", "bin", "bin_order", "period", "n_cells_total"]]
    eig = eig[colmeta.index]

    def _char(e: pd.Series) -> dict:
        gg = e.groupby(colmeta["group"].to_numpy(), dropna=False).mean()
        top = gg.sort_values(ascending=False)
        pk = e.idxmax()
        lo = e.idxmin()
        first_last = e.dropna()
        slope = "rising" if first_last.iloc[-1] - first_last.iloc[0] > 0.5 else (
            "falling" if first_last.iloc[0] - first_last.iloc[-1] > 0.5 else "flat_or_peak")
        return {"peak_col": pk, "peak_group": colmeta.loc[pk, "group"],
                "peak_bin": colmeta.loc[pk, "bin"], "peak_period": colmeta.loc[pk, "period"],
                "trough_group": colmeta.loc[lo, "group"], "trough_bin": colmeta.loc[lo, "bin"],
                "direction": slope,
                "top_groups": ",".join(f"{g}({v:+.2f})" for g, v in top.head(3).items())}

    msum = []
    for m in sorted(mod.unique(), key=lambda s: (len(s), s)):
        members = mg[mg["module"] == m]
        ch = _char(eig.loc[m])
        mc = members["mech_class"].value_counts()
        msum.append({"module": m, "n_genes": len(members), **ch,
                     "hubs": ",".join(members.nsmallest(HUB_TOPN, "hub_rank")["gene"]),
                     "mech_composition": ",".join(f"{k}:{v}" for k, v in mc.items()),
                     "n_specific": int((members["specific_call"] == "specific").sum()),
                     "mean_connectivity": round(float(members["connectivity"].mean()), 4),
                     "mean_within_consensus": round(float(
                         cons_f[np.ix_([gene_ix[g] for g in members["gene"]],
                                       [gene_ix[g] for g in members["gene"]])].mean()), 4)})
    msum = pd.DataFrame(msum)

    # ---- 5. 机制类富集(CL-027: 排除残差类) ----
    enr_in = mg[~mg["mech_class"].isin(residual) & mg["mech_class"].notna() & (mg["mech_class"] != "NA")]
    excluded = sorted(set(mg["gene"][mg["mech_class"].isin(residual)]))
    print(f"[a4] CL-027 残差类排除: {residual} → 剔除 {len(excluded)} 基因 {excluded}"
          f"(改做单基因描述, 不入富集网格); 入富集 {len(enr_in)} 基因")
    e1 = atxenrich.perm_enrichment(enr_in["mech_class"].to_numpy(), enr_in["module"].to_numpy(),
                                   n_perm=n_perm, seed=cfg.seed)
    e1.insert(0, "test", "mech_class_x_module")
    pat = pd.read_csv(cfg.results("tables", sp, "temporal_patterns.tsv"), sep="\t")
    pat = pat[pat["gene"].isin(enr_in["gene"]) & pat["group"].isin(main_groups)]
    shp = pat.merge(enr_in[["gene", "module"]], on="gene", how="inner")
    shp["shape"] = shp.apply(lambda r: r["temporal_pattern"]
                             if r["group"] == (spec.loc[r["gene"], "max_group"]
                                               if r["gene"] in spec.index else None) else None, axis=1)
    shp = shp.dropna(subset=["shape"])
    shp = shp[~shp["shape"].isin(["NA"])]
    if len(shp) > 5:
        e2 = atxenrich.perm_enrichment(shp["mech_class"].to_numpy() if "mech_class" in shp else
                                       shp["gene"].map(mg.set_index("gene")["mech_class"]).to_numpy(),
                                       shp["shape"].to_numpy(), n_perm=n_perm, seed=cfg.seed)
        e2.insert(0, "test", "mech_class_x_shape")
    else:
        e2 = pd.DataFrame(columns=["test"] + atxenrich._ENR_COLS)
    enr = pd.concat([e1, e2], ignore_index=True)
    enr_sig = enr[(enr["q"] < 0.05) & (enr["obs"] >= 3)].sort_values("log2fe", ascending=False)

    # ---- 6. 跨组趋势一致性(共性 vs 异质性) ----
    tr = pd.read_csv(cfg.results("tables", sp, "trend_by_group.tsv"), sep="\t")
    tr = tr[tr["group"].isin(main_groups) & tr["gene"].isin(stat_genes)]
    sig = tr[(tr["q_bin"] < 0.05) & (tr["status"] == "ok")].copy()
    sig["dir"] = np.sign(sig["rho_bin"])
    g = sig.groupby("gene")["dir"].agg(n_up=lambda s: int((s > 0).sum()),
                                       n_down=lambda s: int((s < 0).sum()),
                                       n_sig="count")
    conc = pd.DataFrame(index=[x for x in stat_genes])
    conc.index.name = "gene"
    conc = conc.join(g).fillna(0).astype({"n_up": int, "n_down": int, "n_sig": int})
    conc["module"] = mod.reindex(conc.index).fillna("NA")
    conc["mech_class"] = genes_tab["mech_class"].reindex(conc.index).fillna("NA")

    def _lab(r):
        if r.n_sig == 0:
            return "ns"
        if r.n_down == 0 and r.n_up >= 5:
            return "pan_up"
        if r.n_up == 0 and r.n_down >= 5:
            return "pan_down"
        if r.n_up > 0 and r.n_down > 0:
            return "mixed_direction"
        if r.n_sig == 1:
            return "single_group"
        return "partial_concordant"

    conc["concordance"] = conc.apply(_lab, axis=1)
    conc = conc.sort_values(["concordance", "n_sig"], ascending=[True, False]).reset_index()

    # ---- 7. 写盘 ----
    mg.to_csv(os.path.join(outdir, "modules.tsv"), sep="\t", index=False, float_format="%.6g")
    msum.to_csv(os.path.join(outdir, "module_summary.tsv"), sep="\t", index=False)
    pac_tab.to_csv(os.path.join(outdir, "consensus_k_selection.tsv"), sep="\t", index=False)
    conc.to_csv(os.path.join(outdir, "cross_group_concordance.tsv"), sep="\t", index=False)
    enr.to_csv(os.path.join(outdir, "enrichment_class_x_module.tsv"), sep="\t", index=False,
               float_format="%.6g")
    Z.to_csv(os.path.join(outdir, "gene_zprofile_matrix.tsv"), sep="\t", float_format="%.6g")
    eig.to_csv(os.path.join(outdir, "module_eigenprofile.tsv"), sep="\t", float_format="%.6g")

    # ---- 8. 图 ----
    order_g = list(mg["gene"])
    Zp = Z.loc[order_g]
    colgrp = colmeta["group"].to_numpy()
    fig, ax = plt.subplots(figsize=(15, 11))
    im = ax.imshow(Zp.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    bounds = np.cumsum([0] + [int((colgrp == g).sum()) for g in pd.unique(colgrp)])
    for b in bounds[1:-1]:
        ax.axvline(b - 0.5, color="k", lw=0.8)
    mcolors = plt.get_cmap("tab10")
    for i, g in enumerate(order_g):
        ax.add_patch(plt.Rectangle((-3.2, i - 0.5), 2.6, 1, color=mcolors(int(mod[g]) % 10), lw=0))
    ax.set_yticks(range(len(order_g)))
    ax.set_yticklabels(order_g, fontsize=5)
    ax.set_xticks(range(X.shape[1]))
    ax.set_xticklabels([f"{c}" for c in Zp.columns], rotation=90, fontsize=4.5)
    ax.set_title(f"[{sp}] A4-1 shared spatiotemporal modules (consensus clustering, k={k_best}, "
                 f"gene_set v{cfg.version})", fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.5, label="z(mean log1p CP10K)")
    fig.text(0.012, 0.5, "module colour bar", rotation=90, va="center", fontsize=7)
    fig.tight_layout(rect=[0.03, 0, 1, 1])
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"a4_1_module_heatmap_{sp}.{ext}"), dpi=200)
    plt.close(fig)

    ks = sorted(mod.unique(), key=lambda s: (len(s), s))
    ncol = min(3, len(ks))
    nrow = int(np.ceil(len(ks) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.4 * nrow), squeeze=False)
    gcolors = plt.get_cmap("tab10")
    for i, m in enumerate(ks):
        a = axes[i // ncol][i % ncol]
        for j, grp in enumerate(pd.unique(colgrp)):
            sel = colmeta["group"].to_numpy() == grp
            x = colmeta["bin_order"].to_numpy()[sel]
            y = eig.loc[m].to_numpy()[sel]
            o = np.argsort(x)
            a.plot(x[o], y[o], "-o", ms=2.5, lw=1.1, color=gcolors(j % 10), label=grp)
        a.set_title(f"module {m} (n={int((mod == m).sum())})", fontsize=9)
        a.set_xlabel("developmental bin order")
        a.set_ylabel("eigen z")
        a.axhline(0, color="0.7", lw=0.6)
        if i == 0:
            a.legend(fontsize=5.5, ncol=2, framealpha=0.6)
    for i in range(len(ks), nrow * ncol):
        axes[i // ncol][i % ncol].axis("off")
    fig.suptitle(f"[{sp}] A4-2 module eigen-profiles by cell type (feasible bins only; NA not interpolated)",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"a4_2_eigenprofiles_{sp}.{ext}"), dpi=200)
    plt.close(fig)

    if len(e1):
        piv = e1.pivot(index="a", columns="b", values="log2fe")
        pq = e1.pivot(index="a", columns="b", values="q")
        fig, ax = plt.subplots(figsize=(2.0 + 1.1 * piv.shape[1], 1.4 + 0.42 * piv.shape[0]))
        im = ax.imshow(piv.to_numpy(), aspect="auto", cmap="RdBu_r",
                       vmin=-np.nanmax(np.abs(piv.to_numpy())), vmax=np.nanmax(np.abs(piv.to_numpy())))
        ax.set_xticks(range(piv.shape[1]), [f"mod {c}" for c in piv.columns])
        ax.set_yticks(range(piv.shape[0]), piv.index, fontsize=7)
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                q = pq.to_numpy()[i, j]
                if np.isfinite(q) and q < 0.05:
                    ax.text(j, i, "*", ha="center", va="center", fontsize=10, color="k")
        ax.set_title(f"[{sp}] A4-3 mech_class × module log2FE (*q<0.05; residual classes "
                     f"{residual} excluded per CL-027)", fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.8, label="log2FE")
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(os.path.join(outdir, f"a4_3_enrichment_{sp}.{ext}"), dpi=200)
        plt.close(fig)

    cc = conc["concordance"].value_counts()
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    ax.bar(range(len(cc)), cc.values, color="#4C72B0")
    ax.set_xticks(range(len(cc)), [f"{k}\n(n={v})" for k, v in cc.items()], fontsize=8)
    ax.set_ylabel("genes")
    ax.set_title(f"[{sp}] A4-4 cross-group trend concordance (8 main groups, q<0.05; "
                 f"pan_* = ≥5 groups same direction)", fontsize=9)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"a4_4_concordance_{sp}.{ext}"), dpi=200)
    plt.close(fig)

    # ---- 9. 报告 ----
    L_ = []
    L_.append(f"# A4 共性机制筛选报告 — {sp}（gene_set v{cfg.version}）\n")
    L_.append(f"> 一次性分析层（`scripts/a4_common_mechanism.py`，非 sNN 契约），生成 {t0:%Y-%m-%d %H:%M}，"
              f"seed={cfg.seed}，n_reps={args.n_reps}，n_perm={n_perm}，"
              f"单供体 bin {'已剔除' if args.exclude_single_donor else '保留(与 s02 同口径)'}。\n")
    L_.append("## 1. 口径与输入\n")
    L_.append(f"- 输入：`results/v{cfg.version}/intermediate/{sp}/cube_bin.tsv`（可行 组×bin）、"
              f"`tables/{sp}/specificity_global.tsv`、`trend_by_group.tsv`、`temporal_patterns.tsv`、"
              f"`genes/ataxia_genes.tsv`。")
    L_.append(f"- 空间维 = 8 个主细胞类型组（{', '.join(main_groups)}）；时间维 = 发育 bin；"
              f"UBC（descriptive）与 ALL（reference）不入聚类。")
    L_.append(f"- 聚类矩阵：{len(clust_genes)} 基因 × {X.shape[1]} 个可行 组×bin 列"
              f"（单供体 {n_sd_bins}、pooled {n_pool_bins}）。")
    L_.append(f"- **NA 政策（硬约束 5）**：不可行 组×bin 为 NA，不插值、不填 0；相关用成对完整观测"
              f"（min_periods={MIN_PERIODS}），相关缺失对 {n_nan_pairs}（按距离 1 中性处理）。")
    L_.append(f"- 排除：`not_expressed` {len(not_expr)}（{', '.join(not_expr) or '无'}）、零方差 {len(zerovar)}"
              f"（{', '.join(zerovar) or '无'}）、观测不足 {len(lowobs)}（{', '.join(lowobs) or '无'}）；"
              f"cube 未命中 {len(absent)}（{', '.join(absent) or '无'}）。")
    L_.append(f"- **CL-027**：机制类富集排除残差兜底类 {residual} → 剔除 {len(excluded)} 基因"
              f"（{', '.join(excluded) or '无'}），这些基因改做单基因描述，不产生"
              f"「残差类富集」结论。\n")
    L_.append("## 2. 共识聚类与 k 选择\n")
    L_.append(_md(pac_tab))
    L_.append(f"\n选定 **k={k_best}**（PAC 最小）；模块内平均共识见下表 "
              "`mean_within_consensus`（越接近 1 越稳定）。\n")
    L_.append("## 3. 模块刻画（共享时空程序）\n")
    L_.append(_md(msum[["module", "n_genes", "direction", "peak_group", "peak_bin", "peak_period",
                        "top_groups", "hubs", "mech_composition", "mean_within_consensus"]]))
    L_.append("\n`hubs` = 模块内平均 |corr| 最高的前 "
              f"{HUB_TOPN} 个基因（共性程序的代表节点，后续网络层 n00–n04 的候选入口）。\n")
    L_.append("## 4. 机制类富集（mech_class × module / × shape）\n")
    if len(enr_sig):
        L_.append(f"显著格（q<0.05 且 obs≥3）共 {len(enr_sig)} 个，按 log2FE 降序前 20：\n")
        L_.append(_md(enr_sig.head(20)[["test", "a", "b", "obs", "exp", "log2fe", "p", "q"]],
                      floatfmt=".4g"))
    else:
        L_.append("无显著格（q<0.05 且 obs≥3）——模块与机制类在本样本量下无可检出的对应关系；"
                  "这是**信息性结果**（共性可能体现在时序而非既定机制分类上），不得解读为「无共性机制」。")
    L_.append("")
    L_.append("## 5. 跨组趋势一致性（共性 vs 异质性）\n")
    L_.append(_md(conc["concordance"].value_counts().rename_axis("concordance").reset_index()))
    for lab in ("pan_up", "pan_down", "mixed_direction"):
        sub = conc[conc["concordance"] == lab]
        if len(sub):
            L_.append(f"\n- **{lab}**（{len(sub)}）：{', '.join(sub['gene'].head(40))}"
                      + ("…" if len(sub) > 40 else ""))
    L_.append("\n判读：`pan_*` = ≥5 个主组同向显著（跨细胞类型的**共性**发育程序）；"
              "`mixed_direction` = 同时存在升与降（**异质性**的细胞类型分歧）；"
              "`single_group` = 仅一组显著（细胞类型特异）。\n")
    L_.append("## 6. 局限（写报告时必须一并说明）\n")
    L_.append("- 单物种（鼠）、正常发育数据，**无疾病对照** → 只能说「致病基因在正常小脑发育中的共享时空程序」，"
              "不能说「共性致病机制」已被证明；后者需要 s04/s05 跨物种 + 疾病投影层（b00/b01、n00–n04）。")
    L_.append("- 模块为无监督发现，属**假设生成**；k 由 PAC 选定，需做 k±1 与 min_cells 敏感性后才可写入正文。")
    L_.append("- 单供体/pooled bin 参与聚类（与 s02 同口径），早期胚胎段的模块稳定性弱于晚期。")
    L_.append(f"- 富集样本量小（入富集 {len(enr_in)} 基因 × k={k_best} 模块），q<0.05 的格需人工复核计数（obs≥3 已过滤）。")
    L_.append("\n## 7. 产物清单\n")
    for f in sorted(os.listdir(outdir)):
        L_.append(f"- `{f}`")
    with open(os.path.join(outdir, f"a4_common_mechanism_{sp}.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L_) + "\n")

    print(f"[a4] 模块数 {len(ks)}; 显著富集格 {len(enr_sig)}; 一致性分布 "
          f"{dict(conc['concordance'].value_counts())}")
    print(f"[a4] 产物 → {outdir}  (耗时 {(datetime.now() - t0).seconds}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
