#!/usr/bin/env python3
"""a6_cross_species_modules.py — 人鼠模块对照(一次性分析层, A4 的跨物种延伸)。

定位: **描述性对照**, 不是 s04/s05(跨物种层)也不是 s09(保守性层)——那三层尚未实现,
且需要 spec §11-4 的 bin 级 `stage_map_mouse_to_human_rank`(6 锚点 + 6 插值)。本脚本刻意
**不做 bin 级对齐**, 只做三件在现有口径下站得住的事:

  1. 模块对应关系: A4 在两物种各自得到的模块, 用列联表 + Fisher + NMI 看结构是否保守
     (NMI≈1 = 划分几乎一致; NMI≈0 = 两物种的模块互不相干)。
  2. **细胞类型偏好保守性**(不需时间轴对齐): 每基因取 8 主组的平均表达谱, 人 vs 鼠做 Spearman
     → `rho_group`。这是最干净的一条跨物种指标, 因为"组"的定义两物种同源。
  3. **粗 4 阶段时序保守性**: 用一个显式、可质疑但可复核的粗映射把两物种 bin 归到 4 个阶段
     (S1 产前早期 / S2 围产 / S3 出生后 / S4 成年), 聚合成 8 组 × 4 阶段 = 32 格后做 Spearman
     → `rho_stage`。粗映射见 STAGE_MAP, 结论只能表述为"阶段级"而非"时点级"。

硬口径: NA 不插值不填 0(不可行 组×bin 直接不进聚合); 只用 group_tiers.main 的 8 组;
相关一律 pairwise-complete 并报告实际用格子数 n_cells_used。

用法: /home/victor/.conda/envs/omics/bin/python -I scripts/a6_cross_species_modules.py <root>
"""
from __future__ import annotations

import os
import sys
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atx import config as atxconfig                     # noqa: E402

# 粗 4 阶段映射(本脚本的分析假设, 不是 spec §11-4 的 bin 级映射)。
# 依据: 人 CS18–11w 与鼠 E10.5–E15.5 同属小脑神经发生/祖细胞扩增期; 人 17w–NB 与鼠 E17.5–P0
# 同属晚胎/围产; 人 6mo–3y 与鼠 P4–P14 同属出生后成熟; 人 42y+ 与鼠 9wk 同属成年。
STAGE_MAP = {
    "human": {"S1_prenatal_early": ["CS18", "CS19", "CS22", "9w", "11w"],
              "S2_perinatal": ["17w", "20w", "NB"],
              "S3_postnatal": ["6mo", "7mo", "9mo", "3y"],
              "S4_adult": ["AD", "42y", "44y", "46y", "52y"]},
    "mouse": {"S1_prenatal_early": ["E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5"],
              "S2_perinatal": ["E17.5", "P0"],
              "S3_postnatal": ["P4", "P7", "P14"],
              "S4_adult": ["9wk"]},
}
STAGES = ["S1_prenatal_early", "S2_perinatal", "S3_postnatal", "S4_adult"]
MIN_CELLS_FOR_RHO = 8      # pairwise-complete 格子数下限, 低于此不出 rho(记 NA)


def _md(df: pd.DataFrame, floatfmt: str | None = None) -> str:
    d = df.copy()
    if floatfmt:
        for c in d.columns:
            if pd.api.types.is_numeric_dtype(d[c]):
                d[c] = d[c].map(lambda v: "" if pd.isna(v) else format(float(v), floatfmt))
    cols = [str(c) for c in d.columns]
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for r in d.astype(object).where(d.notna(), None).to_numpy():
        out.append("| " + " | ".join("" if v is None else str(v).replace("|", "/") for v in r) + " |")
    return "\n".join(out)


def nmi(a: np.ndarray, b: np.ndarray) -> float:
    """归一化互信息(0=互不相干, 1=划分完全一致); 不依赖 sklearn。"""
    la, lb = sorted(set(a.tolist())), sorted(set(b.tolist()))
    if len(la) < 2 or len(lb) < 2:
        return float("nan")
    ia = np.array([la.index(x) for x in a])
    ib = np.array([lb.index(x) for x in b])
    n = len(a)
    H = lambda v, k: -sum((np.sum(v == j) / n) * np.log((np.sum(v == j) / n))
                          for j in range(k) if np.sum(v == j) > 0)
    mi = sum(np.sum((ia == i) & (ib == j)) / n * np.log(
        (np.sum((ia == i) & (ib == j)) / n) / ((np.sum(ia == i) / n) * (np.sum(ib == j) / n)))
        for i in range(len(la)) for j in range(len(lb)) if np.sum((ia == i) & (ib == j)) > 0)
    ha, hb = H(ia, len(la)), H(ib, len(lb))
    return float(2 * mi / (ha + hb)) if (ha + hb) > 0 else float("nan")


def group_stage_profile(cfg, sp: str, genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """→ (8 组平均谱, 8 组 × 4 阶段谱); 只用可行 组×bin, NA 不填。"""
    bt = pd.read_csv(cfg.results("intermediate", sp, "cube_bin.tsv"), sep="\t")
    main_groups = list(cfg.main_groups())
    bt = bt[bt["feasible"] & bt["group"].isin(main_groups)].copy()
    smap = {}
    for st, bins in STAGE_MAP[sp].items():
        for b in bins:
            smap[b] = st
    bt["stage"] = bt["bin"].astype(str).map(smap)
    unmapped = sorted(set(bt["bin"].astype(str)) - set(smap))
    if unmapped:
        print(f"[a6] 警告({sp}): 未映射到 4 阶段的 bin {unmapped} → 该 bin 不进阶段聚合")
        bt = bt[bt["stage"].notna()]
    gcols = [f"log_{g}" for g in genes if f"log_{g}" in bt.columns]
    if not gcols:
        raise SystemExit(f"[a6] {sp}: cube_bin.tsv 里没有目标基因的 log_ 列(基因名不匹配?)")

    def _agg(keys) -> pd.DataFrame:
        t = bt.groupby(keys)[gcols].mean().T          # 基因(带 log_ 前缀) × 该维度的水平
        t.index = [i[len("log_"):] for i in t.index]  # 去前缀 → 基因符号
        if isinstance(t.columns, pd.MultiIndex):      # (group, stage) → "group|stage"
            t.columns = ["|".join(str(x) for x in c) for c in t.columns]
        return t.reindex([g for g in genes if g in set(t.index)])

    # 返回 (基因 × 8 组谱, 基因 × 8组×4阶段=至多 32 格谱); 后者才对得上 MIN_CELLS_FOR_RHO
    return _agg("group"), _agg(["group", "stage"])


def main() -> int:
    root = os.path.abspath(sys.argv[1])
    cfg = atxconfig.load(os.path.join(root, "config", "config.yaml"))
    ver = cfg.version
    a4 = {sp: cfg.results("analysis", f"a4_common_mechanism_{sp}") for sp in ("human", "mouse")}
    for sp, d in a4.items():
        if not os.path.exists(os.path.join(d, "modules.tsv")):
            raise SystemExit(f"[a6] 缺 A4 产物 {d}/modules.tsv —— 先跑 scripts/a4_common_mechanism.py . --species {sp}")
    mh = pd.read_csv(os.path.join(a4["human"], "modules.tsv"), sep="\t")
    mm = pd.read_csv(os.path.join(a4["mouse"], "modules.tsv"), sep="\t")
    genes = sorted(set(mh["gene"]) & set(mm["gene"]))
    print(f"[a6] A4 入聚类基因: 人 {len(mh)} / 鼠 {len(mm)}; 两物种共有 {len(genes)}")

    gt = pd.read_csv(cfg.paths["genes_table"], sep="\t").set_index("hgnc_symbol")
    outdir = cfg.results("analysis", "a6_cross_species", mkdir=True)

    # ---- 1. 模块对应 + NMI ----
    j = mh[["gene", "module"]].merge(mm[["gene", "module"]], on="gene", suffixes=("_h", "_m"))
    j = j[j["gene"].isin(genes)]
    ct = pd.crosstab(j["module_h"], j["module_m"])
    ct.index = [f"H{ i }" for i in ct.index]
    ct.columns = [f"M{c}" for c in ct.columns]
    n_total = int(j.shape[0])
    recs = []
    for h in ct.index:
        for m in ct.columns:
            a = int(ct.loc[h, m])
            if a == 0:
                continue
            nh, nm = int(ct.loc[h].sum()), int(ct[m].sum())
            odds, p = fisher_exact([[a, nh - a], [nm - a, n_total - nh - nm + a]])
            recs.append({"module_human": h, "module_mouse": m, "n_genes": a,
                         "n_human_module": nh, "n_mouse_module": nm,
                         "overlap_frac_of_human": round(a / nh, 3) if nh else np.nan,
                         "overlap_frac_of_mouse": round(a / nm, 3) if nm else np.nan,
                         "p_fisher": p})
    corr = pd.DataFrame(recs).sort_values("p_fisher")
    nmi_v = nmi(j["module_h"].to_numpy(), j["module_m"].to_numpy())
    best = (corr.sort_values(["p_fisher", "n_genes"], ascending=[True, False])
                .drop_duplicates("module_human").set_index("module_human"))
    print(f"[a6] 模块划分 NMI = {nmi_v:.3f}(1=完全一致, 0=无关)")

    # ---- 2/3. 组谱与阶段谱的跨物种相关 ----
    gg_h, gs_h = group_stage_profile(cfg, "human", genes)
    gg_m, gs_m = group_stage_profile(cfg, "mouse", genes)
    common_groups = [g for g in gg_h.columns if g in gg_m.columns]
    Zh = gg_h[common_groups].sub(gg_h[common_groups].mean(axis=1), axis=0)
    Zm = gg_m[common_groups].sub(gg_m[common_groups].mean(axis=1), axis=0)
    common_cells = sorted(set(gs_h.columns) & set(gs_m.columns))   # "组|阶段" 交集(至多 32 格)
    print(f"[a6] 组×阶段格: 人 {gs_h.shape[1]} / 鼠 {gs_m.shape[1]}; 两物种共有 {len(common_cells)}")
    Sh = gs_h[common_cells].sub(gs_h[common_cells].mean(axis=1), axis=0)
    Sm = gs_m[common_cells].sub(gs_m[common_cells].mean(axis=1), axis=0)

    rows = []
    for g in genes:
        def _rho(A, B):
            x, y = A.loc[g].to_numpy(float), B.loc[g].to_numpy(float)
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 3 or np.nanstd(x[ok]) == 0 or np.nanstd(y[ok]) == 0:
                return np.nan, int(ok.sum())
            r = spearmanr(x[ok], y[ok]).statistic
            return (float(r) if np.isfinite(r) else np.nan), int(ok.sum())
        rg, ng = _rho(Zh, Zm)
        rs, ns = _rho(Sh, Sm)
        rows.append({"gene": g,
                     "rho_group": rg if ng >= MIN_CELLS_FOR_RHO else np.nan,
                     "n_groups_used": ng,
                     "rho_stage": rs if ns >= MIN_CELLS_FOR_RHO else np.nan,
                     "n_stage_cells_used": ns,
                     "module_human": j.set_index("gene")["module_h"].get(g),
                     "module_mouse": j.set_index("gene")["module_m"].get(g),
                     "mech_class": gt.loc[g, "mech_class"] if g in gt.index else "NA",
                     "max_group_human": mh.set_index("gene")["max_group"].get(g),
                     "max_group_mouse": mm.set_index("gene")["max_group"].get(g),
                     "same_max_group": bool(mh.set_index("gene")["max_group"].get(g)
                                            == mm.set_index("gene")["max_group"].get(g))})
    prof = pd.DataFrame(rows)
    rc_path = cfg.results("sensitivity", "robust_core_genes.tsv")
    robust = set()
    if os.path.exists(rc_path):
        rc = pd.read_csv(rc_path, sep="\t")
        both = rc.groupby("gene")["species"].nunique()
        robust = set(both[both > 1].index)
    prof["in_cross_species_robust_core"] = prof["gene"].isin(robust)

    # ---- 写盘 ----
    ct.to_csv(os.path.join(outdir, "module_contingency.tsv"), sep="\t")
    corr.to_csv(os.path.join(outdir, "module_correspondence.tsv"), sep="\t", index=False,
                float_format="%.6g")
    prof.to_csv(os.path.join(outdir, "gene_cross_species_profiles.tsv"), sep="\t", index=False,
                float_format="%.6g")

    # ---- 图 ----
    fig, ax = plt.subplots(figsize=(2.2 + 0.9 * ct.shape[1], 1.6 + 0.55 * ct.shape[0]))
    im = ax.imshow(ct.to_numpy(), aspect="auto", cmap="Blues")
    ax.set_xticks(range(ct.shape[1]), ct.columns, fontsize=7)
    ax.set_yticks(range(ct.shape[0]), ct.index, fontsize=7)
    for i in range(ct.shape[0]):
        for k in range(ct.shape[1]):
            v = ct.to_numpy()[i, k]
            if v:
                ax.text(k, i, int(v), ha="center", va="center", fontsize=7,
                        color="w" if v > ct.to_numpy().max() * 0.6 else "k")
    ax.set_title(f"human × mouse A4 module gene counts (NMI={nmi_v:.3f}, n={n_total} genes)",
                 fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.8, label="genes")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"a6_1_module_correspondence.{ext}"), dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    ok = prof.dropna(subset=["rho_group", "rho_stage"])
    ax.scatter(ok.loc[~ok.in_cross_species_robust_core, "rho_group"],
               ok.loc[~ok.in_cross_species_robust_core, "rho_stage"],
               s=22, c="#9DB1C6", label=f"other genes (n={int((~ok.in_cross_species_robust_core).sum())})",
               edgecolors="none")
    ax.scatter(ok.loc[ok.in_cross_species_robust_core, "rho_group"],
               ok.loc[ok.in_cross_species_robust_core, "rho_stage"],
               s=46, c="#C44E52",
               label=f"cross-species robust core (n={int(ok.in_cross_species_robust_core.sum())})",
               edgecolors="k", linewidths=0.4)
    ax.axhline(0, color="0.7", lw=0.7)
    ax.axvline(0, color="0.7", lw=0.7)
    ax.set_xlabel("rho_group: human vs mouse, Spearman over 8 main cell-type groups")
    ax.set_ylabel("rho_stage: human vs mouse, Spearman over 8 groups x 4 coarse stages")
    ax.set_title("cross-species profile conservation (descriptive; coarse 4-stage map)", fontsize=9)
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"a6_2_profile_scatter.{ext}"), dpi=200)
    plt.close(fig)

    # ---- 报告 ----
    q = prof["rho_group"].dropna()
    qs = prof["rho_stage"].dropna()
    same_mg = int(prof["same_max_group"].sum())
    L = [f"# A6 人鼠模块对照报告（gene_set v{ver}）", "",
         f"> 一次性分析层（`scripts/a6_cross_species_modules.py`），生成 {datetime.now():%Y-%m-%d %H:%M}。",
         "> **定位**：描述性对照，**不是** s04/s05（跨物种层）或 s09（保守性层）——那三层未实现，"
         "且需 spec §11-4 的 bin 级 `stage_map_mouse_to_human_rank`（6 锚点 + 6 插值）。"
         "本脚本刻意不做 bin 级对齐。", "",
         "## 1. 口径", "",
         f"- 共有基因（两物种 A4 都入聚类）：**{len(genes)}**；模块划分 **NMI = {nmi_v:.3f}**"
         f"（1=完全一致，0=无关）。",
         f"- 细胞类型偏好：8 主组平均谱，Spearman；阈值 n≥{MIN_CELLS_FOR_RHO} 才出值。",
         f"- 粗 4 阶段：S1 产前早期 / S2 围产 / S3 出生后 / S4 成年，映射见脚本 `STAGE_MAP`；"
         f"人侧 bin 轴 = 主口径（`collapse_adult_bins`，成人合并为 AD）。",
         "- NA 政策：不可行 组×bin 不进聚合，不插值、不填 0；相关为 pairwise-complete 并报实际格子数。", "",
         "## 2. 主要结果", "",
         f"- `rho_group`（细胞类型偏好保守性）：n={len(q)}，中位 **{q.median():+.3f}**，"
         f"四分位 [{q.quantile(.25):+.3f}, {q.quantile(.75):+.3f}]，"
         f">0.6 的基因 {int((q > 0.6).sum())} 个，<0 的 {int((q < 0).sum())} 个。",
         f"- `rho_stage`（粗阶段时序保守性）：n={len(qs)}，中位 **{qs.median():+.3f}**，"
         f">0.6 的 {int((qs > 0.6).sum())} 个。",
         f"- 主组（τ argmax）两物种相同的基因：**{same_mg}/{len(prof)}**。",
         f"- 跨物种稳健核心集（敏感性全口径显著且同向，见 `sensitivity/robust_core_genes.tsv`）"
         f"命中 {int(prof['in_cross_species_robust_core'].sum())} 个基因；其 `rho_group` 中位 "
         f"**{prof.loc[prof.in_cross_species_robust_core, 'rho_group'].median():+.3f}**，"
         f"其余基因中位 {prof.loc[~prof.in_cross_species_robust_core, 'rho_group'].median():+.3f}。", "",
         "## 3. 模块对应（人模块 → 最佳鼠模块，按 Fisher p）", "",
         _md(best.reset_index()[["module_human", "module_mouse", "n_genes", "n_human_module",
                                 "n_mouse_module", "overlap_frac_of_human", "p_fisher"]],
             floatfmt=".4g") if len(best) else "（无）", "",
         "## 4. 列联表（基因数）", "", _md(ct.reset_index().rename(columns={"index": ""})), "",
         "## 5. 细胞类型偏好最保守 / 最不保守的基因", "",
         _md(prof.dropna(subset=["rho_group"]).nlargest(12, "rho_group")
             [["gene", "rho_group", "rho_stage", "mech_class", "max_group_human",
               "max_group_mouse", "in_cross_species_robust_core"]], floatfmt=".3f"), "",
         _md(prof.dropna(subset=["rho_group"]).nsmallest(12, "rho_group")
             [["gene", "rho_group", "rho_stage", "mech_class", "max_group_human",
               "max_group_mouse", "in_cross_species_robust_core"]], floatfmt=".3f"), "",
         "## 6. 按机制类看保守性（中位 rho_group）", "",
         _md(prof.dropna(subset=["rho_group"]).groupby("mech_class")["rho_group"]
             .agg(n="count", median="median", q25=lambda s: s.quantile(.25),
                  q75=lambda s: s.quantile(.75)).reset_index(), floatfmt=".3f"), "",
         "## 7. 局限（引用前必读）", "",
         "- 4 阶段映射是**本脚本的分析假设**，不是 spec §11-4 的 bin 级映射；"
         "换一套映射 `rho_stage` 会变，故只可作阶段级表述。",
         "- 人侧 Purkinje/DN/UBC/Progenitor 产后无数据（右删失），这几组的跨物种相关实际只比较了产前段；"
         "人侧成人段仅 1 个 AD bin（且 2/4 供体为合并症供体，主口径未剔除）。",
         "- 模块由无监督聚类得到，A4 已报告其共识稳定性中等（PAC≈0.48）；模块对应关系属**假设生成**。",
         "- 两物种数据来自不同平台/注释（Sepp 人 snRNA vs Sepp 鼠），`rho_group` 含技术差异成分，"
         "不能直接当作"
         "生物学保守性强度。", "",
         "## 8. 产物", ""]
    for f in sorted(os.listdir(outdir)):
        L.append(f"- `{f}`")
    with open(os.path.join(outdir, "a6_cross_species.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"[a6] NMI={nmi_v:.3f} | rho_group 中位={q.median():+.3f}(n={len(q)}) | "
          f"rho_stage 中位={qs.median():+.3f}(n={len(qs)}) | 主组相同 {same_mg}/{len(prof)}")
    print(f"[a6] 产物 → {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
