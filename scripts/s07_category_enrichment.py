#!/usr/bin/env python3
"""s07_category_enrichment.py — 第 3 层:机制类 × 模块/shape 置换富集 + 统一峰时轴(spec §6)。

输入(两物种): tables/{sp}/specificity_global.tsv、peak_by_group.tsv、temporal_patterns.tsv;
             intermediate/{sp}/cube_bin.tsv(仅 fig3d 代表轨迹用); genes/ataxia_genes.tsv。
产物: tables/cross/{module_assignment, enrichment_class_x_module, enrichment_class_x_shape,
       fisher_class_x_shape, peak_rank_by_class}.tsv
     figures/cross/{fig3a_enrichment_heatmap, fig3b_shape_mosaic, fig3c_peak_ridge,
       fig3d_representative_traj}.{pdf,png}
     logs/s07_report.tsv

口径(spec §6 + §11-4 + CL-027):
  · 模块 = G* × P*(assign_modules 已实现: G*=global τ argmax, max_fraction<broad_fraction_threshold
    → G*=broad 且峰取 ALL 组; G* 可行 bin < min_bins_for_trend → 模块 NA, **不进检验**);
  · 富集 = O/E 的 log2FE + 标签置换(固定模块、打乱机制类, n_perm 取 config, 恒等排列计入 → p ≥ 1/n_perm),
    BH 族 = 单次调用的全部格(= 物种内全格);
  · **CL-027: 残差兜底类 other/unknown 不进富集网格**(perm_enrichment 另有硬保护), 被排除基因逐条打印;
  · 峰时轴 peak_rank_human_equiv 按 config peak_source_rule: G* ∉ human_censored_groups → 人峰 bin 序号;
    G* ∈ {Purkinje, DN, UBC}(人侧产后右删失) → 鼠峰 bin 经 stage_map_mouse_to_human_rank 映射;
    鼠侧也不可行 → NA。**必须同时报三套轴**(mixed 主 / human_only / mouse_only), 三者方向一致才可写正文,
    不一致时以 mouse_only 为准(spec §6);
  · stage_map 有 6 个插值锚点(§11-4) → 凡用鼠峰映射的结果, 图注与日志均标"部分插值"。
用法: <omics python> -I scripts/s07_category_enrichment.py <root>
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# 图注含中文 → 需要 CJK 字体回退(本机 Noto Sans CJK 可用); 缺字体时 matplotlib 只警告并画方框, 不报错,
# 故显式设置并在下方 rec 里留痕, 避免"图能出但中文全是方框"这种静默降级。
plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC", "Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atx import config as atxconfig                      # noqa: E402
from atx import stats_enrich as atxenrich                # noqa: E402
from atx import genes as atxgenes                        # noqa: E402
from atx import plots as atxplots                        # noqa: E402

REPORT: list[dict] = []
N_FAIL = 0


def rec(check, expected, observed, level, note=""):
    """日志行(列约定同 s01_report: check/expected/observed/level/note)。"""
    global N_FAIL
    REPORT.append({"check": check, "expected": str(expected), "observed": str(observed),
                   "level": level, "note": note})
    mark = {"PASS": "✓", "WARN": "⚠", "FAIL": "✗"}.get(level, "?")
    print(f"[{mark} {level}] {check} | 期望={expected} 实测={observed}" + (f" | {note}" if note else ""))
    N_FAIL += level == "FAIL"


def _finish(cfg, t0):
    logs = cfg.results("logs", mkdir=True)
    pd.DataFrame(REPORT).to_csv(os.path.join(logs, "s07_report.tsv"), sep="\t", index=False)
    n_pass = sum(1 for r in REPORT if r["level"] == "PASS")
    n_warn = sum(1 for r in REPORT if r["level"] == "WARN")
    print(f"\n汇总: PASS={n_pass} WARN={n_warn} FAIL={N_FAIL}  "
          f"报告: {os.path.join(logs, 's07_report.tsv')}  ({time.time() - t0:.0f}s)")


def modules_for_species(cfg, sp: str, genes_tab: pd.DataFrame) -> pd.DataFrame:
    """assign_modules + shape(temporal_pattern at peak_group_used) + mech_class/tier 注记。"""
    td = cfg.results("tables", sp)
    spec_g = pd.read_csv(os.path.join(td, "specificity_global.tsv"), sep="\t")
    peak = pd.read_csv(os.path.join(td, "peak_by_group.tsv"), sep="\t")
    patt = pd.read_csv(os.path.join(td, "temporal_patterns.tsv"), sep="\t")
    ma = atxenrich.assign_modules(spec_g, peak, cfg)
    ma["species"] = sp
    # pandas 3 的 str dtype 下 .astype(str) **不再**把 NaN 变成 "nan"(仍是 NA), 故先经 object 显式填 "NA";
    # 被 s02 掩膜的低表达伪标签行(P7 裁决①) temporal_pattern 即为 NA → 这里统一成字符串 "NA",
    # 后续富集**排除** shape=="NA" 的行(NA 政策: 缺证据不当作一个类别, 也不当 False)。
    tp = patt["temporal_pattern"]
    tp = tp.astype("object").where(tp.notna(), "NA").astype(str)
    pidx = patt.assign(_tp=tp.to_numpy()).set_index(["gene", "group"])["_tp"]
    ma["shape"] = [str(pidx.get((r.gene, r.peak_group_used), "NA")) for r in ma.itertuples()]
    ma = ma.merge(genes_tab[["hgnc_symbol", "mech_class", "tier"]],
                  left_on="gene", right_on="hgnc_symbol", how="left").drop(columns=["hgnc_symbol"])
    ma["mech_class"] = ma["mech_class"].fillna("NA")
    return ma


def peak_ranks(cfg, ma_all: pd.DataFrame) -> pd.DataFrame:
    """peak_rank_human_equiv 三套轴(mixed / human_only / mouse_only), spec §6 + config peak_source_rule。"""
    smap = {str(k): float(v) for k, v in cfg.raw["stage_map_mouse_to_human_rank"].items()}
    rule = cfg.raw["peak_source_rule"]
    censored = set(rule["human_censored_groups"])
    # 鼠 bin 标签 → stage_map 键(config 用原始 author_stage 写法 "adult (9 weeks)", 表内是 time_labels "9wk")
    mouse_labels = cfg.species("mouse").labels
    label2key = {lab: lab for lab in mouse_labels}
    for lab in mouse_labels:
        if lab not in smap:                                  # 例: 9wk → "adult (9 weeks)"
            cand = [k for k in smap if k.startswith(lab) or lab in k.replace(" ", "")]
            if not cand:
                cand = [k for k in smap if k.split()[0].lower() in ("adult",) and lab.lower() in ("9wk", "adult")]
            if cand:
                label2key[lab] = cand[0]
    raw_h = cfg.species("human")
    human_order = {lab: i + 1 for i, lab in enumerate(raw_h.labels)}
    # 合并口径(collapse_adult_bins=true 为人侧主口径)下 peak_by_group 用合并标签(如 AD),
    # 不在原始 16-bin 标签里 → 此前被静默丢弃(human_only/mixed 轴缺基因)。给合并 bin 一个
    # 与 stage_map 同尺度(原始人 bin 序号 1..16)的秩 = 其 constituent 原始秩的质心
    # (AD = mean(13,14,15,16) = 14.5, 与鼠 adult 映射 14.5 一致), 由 config 派生不硬编码。
    eff_h = cfg.effective_species("human")
    if eff_h.labels != raw_h.labels:
        for cl in eff_h.labels:
            if cl in human_order:
                continue
            constituents = [raw_h.label_of[st] for st in raw_h.label_of
                            if eff_h.label_of.get(st) == cl and raw_h.label_of[st] in human_order]
            if constituents:
                human_order[cl] = float(np.mean([human_order[c] for c in constituents]))

    hum = ma_all[ma_all["species"] == "human"].set_index("gene")
    mou = ma_all[ma_all["species"] == "mouse"].set_index("gene")
    recs = []
    for gene in sorted(set(hum.index) | set(mou.index)):
        h = hum.loc[gene] if gene in hum.index else None
        m = mou.loc[gene] if gene in mou.index else None
        g_star = None
        for row in (h, m):
            if row is not None and str(row["status"]) == "ok" and str(row["g_star"]) != "NA":
                g_star = str(row["g_star"])
                break
        # 人峰(仅当 G* 非右删失组且人侧模块 ok)
        h_rank = np.nan
        if h is not None and str(h["status"]) == "ok" and g_star and g_star not in censored \
                and g_star != "NA" and str(h["peak_bin"]) in human_order:
            h_rank = float(human_order[str(h["peak_bin"])])
        # 鼠峰 → 人等价秩(部分插值)
        m_rank = np.nan
        if m is not None and str(m["status"]) == "ok" and str(m["peak_bin"]) != "NA":
            key = label2key.get(str(m["peak_bin"]))
            if key in smap:
                m_rank = smap[key]
        if g_star and g_star in censored:
            mixed, src = m_rank, ("mouse" if np.isfinite(m_rank) else "NA")
        elif np.isfinite(h_rank):
            mixed, src = h_rank, "human"
        else:
            mixed, src = (m_rank, "mouse") if np.isfinite(m_rank) else (np.nan, "NA")
        recs.append({"gene": gene, "g_star": g_star or "NA",
                     "peak_bin_human": str(h["peak_bin"]) if h is not None else "NA",
                     "peak_bin_mouse": str(m["peak_bin"]) if m is not None else "NA",
                     "peak_rank_human_equiv": mixed, "peak_source": src,
                     "peak_rank_human_only": h_rank, "peak_rank_mouse_only": m_rank,
                     "censored_group": bool(g_star in censored)})
    return pd.DataFrame(recs)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--config", default=None)
    ap.add_argument("--n-perm", type=int, default=None)
    ap.add_argument("--out-subdir", default=None)
    ap.add_argument("--sens", action="append", default=None, metavar="KEY=VALUE")
    args = ap.parse_args()
    root = os.path.abspath(args.root)
    t0 = time.time()
    cfg = atxconfig.load(args.config or os.path.join(root, "config", "config.yaml"))
    atxconfig.apply_sens_overrides(cfg, args.sens)
    if args.out_subdir:
        cfg.results_subdir = args.out_subdir.strip("/")
    cfg.check_version_gate()
    ecfg = cfg.stats_cfg["enrichment"]
    n_perm = int(args.n_perm or ecfg["n_perm"])
    residual = list(ecfg.get("exclude_mech_classes") or atxenrich.RESIDUAL_MECH_CLASSES)

    print(f"== s07 cross | version={cfg.version}(draft={cfg.draft}) | n_perm={n_perm}"
          + (f" | 产物隔离: results/v{cfg.version}/{cfg.results_subdir}/" if cfg.results_subdir else ""))

    genes_tab = pd.read_csv(cfg.paths["genes_table"], sep="\t")
    # panel 成员走 load_panel(尊重 tiers_main + include_extended_tier 敏感性; 硬约束 6 单一事实来源),
    # 不再硬编码 tier=="core"(否则 extended 敏感性下物种内/跨物种 panel 静默不一致)。
    panel = atxgenes.load_panel(cfg)
    panel = panel[panel["in_stat_panel"]]
    rec("统计 panel", "tiers_main & in_panel(load_panel)", f"{len(panel)} 基因", "PASS")

    # ---- 1. 模块指派 ----
    ma_all = pd.concat([modules_for_species(cfg, sp, genes_tab) for sp in ("human", "mouse")],
                       ignore_index=True)
    crossdir = cfg.results("tables", "cross", mkdir=True)
    for sp in ("human", "mouse"):
        sub = ma_all[ma_all["species"] == sp]
        ok = int((sub["status"] == "ok").sum())
        rec(f"模块指派 {sp}", "status=ok 者进检验",
            f"{ok}/{len(sub)} ok; 其余 {dict(sub['status'].value_counts())}",
            "PASS" if ok else "WARN")
    ma_all["in_panel"] = ma_all["gene"].isin(set(panel["hgnc_symbol"]))

    # ---- 2. 富集(CL-027 排除残差类) ----
    enr_frames, fish_frames = [], []
    for sp in ("human", "mouse"):
        sub = ma_all[(ma_all["species"] == sp) & (ma_all["status"] == "ok") & ma_all["in_panel"]]
        excl = sorted(set(sub.loc[sub["mech_class"].isin(residual), "gene"]))
        keep = sub[~sub["mech_class"].isin(residual) & sub["mech_class"].ne("NA")]
        print(f"[s07] {sp}: CL-027 残差类 {residual} 排除 {len(excl)} 基因 {excl}; 入富集 {len(keep)}")
        rec(f"CL-027 残差类排除 {sp}", f"排除 {residual}",
            f"{len(excl)} 基因: {excl or '无'}; 入富集 {len(keep)}", "PASS",
            note="残差类不是一种机制, 成员改做单基因描述(CL-027/plan §0-3/spec §6)")
        if len(keep) < 5:
            rec(f"富集 {sp}", "样本量", f"入富集基因仅 {len(keep)} → 跳过", "WARN")
            continue
        keep_m = keep[~keep["module"].astype(str).isin(["NA", "nan", ""])]
        keep_s = keep[~keep["shape"].astype(str).isin(["NA", "nan", ""])]
        rec(f"NA 行排除 {sp}", "模块/shape 为 NA 者不进检验(NA 政策)",
            f"class×module 入 {len(keep_m)}/{len(keep)}; class×shape 入 {len(keep_s)}/{len(keep)}"
            f"(shape=NA {len(keep) - len(keep_s)} 行, 多为 not_expressed 的掩膜伪标签)", "PASS")
        e1 = atxenrich.perm_enrichment(keep_m["mech_class"].to_numpy(), keep_m["module"].to_numpy(),
                                       n_perm=n_perm, seed=cfg.seed)
        e1.insert(0, "species", sp)
        e2 = atxenrich.perm_enrichment(keep_s["mech_class"].to_numpy(), keep_s["shape"].to_numpy(),
                                       n_perm=n_perm, seed=cfg.seed)
        e2.insert(0, "species", sp)
        f2 = atxenrich.fisher_by_column(keep_s["mech_class"].to_numpy(), keep_s["shape"].to_numpy())
        f2.insert(0, "species", sp)
        enr_frames += [e1, e2]
        fish_frames.append(f2)
        n_sig = int(((e1["q"] < 0.05) & (e1["obs"] >= 3)).sum() + ((e2["q"] < 0.05) & (e2["obs"] >= 3)).sum())
        top = pd.concat([e1.assign(test="class_x_module"), e2.assign(test="class_x_shape")])
        top = top[(top["q"] < 0.05) & (top["obs"] >= 3)].sort_values("log2fe", ascending=False)
        rec(f"富集显著格 {sp}", "q<0.05 且 obs≥3",
            f"{n_sig} 格" + (f"; 最强 {top.iloc[0]['test']} {top.iloc[0]['a']}×{top.iloc[0]['b']} "
                            f"log2FE={top.iloc[0]['log2fe']:.2f} q={top.iloc[0]['q']:.2g}" if len(top) else ""),
            "PASS", note="无显著格属信息性结果(共性可能体现在时序而非机制分类), 不得解读为「无共性机制」")
    enr = pd.concat(enr_frames, ignore_index=True)
    enr["test"] = np.where(enr["b"].str.contains(r"\|"), "class_x_module", "class_x_shape")
    fish = pd.concat(fish_frames, ignore_index=True)

    # ---- 3. 峰时轴(三套) + 类别检验 ----
    pr = peak_ranks(cfg, ma_all)
    pr = pr.merge(genes_tab[["hgnc_symbol", "mech_class", "tier"]],
                  left_on="gene", right_on="hgnc_symbol", how="left").drop(columns=["hgnc_symbol"])
    pr["in_panel"] = pr["gene"].isin(set(panel["hgnc_symbol"]))
    prk = pr[pr["in_panel"] & ~pr["mech_class"].isin(residual) & pr["mech_class"].notna()].copy()
    axes = {"mixed": "peak_rank_human_equiv", "human_only": "peak_rank_human_only",
            "mouse_only": "peak_rank_mouse_only"}
    kw_rows, pair_frames = [], []
    for axis, col in axes.items():
        sub = prk[np.isfinite(prk[col])]
        kw_p, pairs = atxenrich.group_compare(sub[col].to_numpy(), sub["mech_class"].to_numpy())
        pairs.insert(0, "axis", axis)
        pair_frames.append(pairs)
        med = sub.groupby("mech_class")[col].median()
        kw_rows.append({"axis": axis, "n_genes": len(sub), "n_classes": int(sub["mech_class"].nunique()),
                        "kw_p": kw_p, **{f"median_{k}": round(float(v), 3) for k, v in med.items()}})
        rec(f"峰时轴 {axis}", "KW + 两两 MW(BH) + Cliff's delta",
            f"n={len(sub)} 基因 / {int(sub['mech_class'].nunique())} 类, KW p={kw_p:.3g}"
            + (f", 显著对(q<0.05) {int((pairs['q_mw'] < 0.05).sum())}" if len(pairs) else ""),
            "PASS" if len(sub) >= 8 else "WARN",
            note="鼠峰经 stage_map 映射, 含 6 个插值锚点(§11-4) → 部分插值" if axis != "human_only" else "")
    peak_tab = pd.DataFrame(kw_rows)
    pairs_all = pd.concat(pair_frames, ignore_index=True)

    # 三套轴方向一致性(spec §6: 三者一致才可写正文; 不一致以 mouse_only 为准)
    def _signs(axis):
        col = axes[axis]
        sub = prk[np.isfinite(prk[col])]
        return sub.groupby("mech_class")[col].median()
    med = {a: _signs(a) for a in axes}
    concord, total = 0, 0
    for a in ("mixed", "human_only"):
        common = sorted(set(med[a].index) & set(med["mouse_only"].index))
        for i in range(len(common)):
            for j in range(i + 1, len(common)):
                g1, g2 = common[i], common[j]
                da = med[a][g1] - med[a][g2]
                dm = med["mouse_only"][g1] - med["mouse_only"][g2]
                if abs(da) < 1e-9 or abs(dm) < 1e-9:
                    continue
                total += 1
                concord += int(np.sign(da) == np.sign(dm))
    rate = concord / total if total else float("nan")
    rec("三套峰轴方向一致性", "两两类别中位差符号 vs mouse_only 一致率",
        f"{concord}/{total} = {rate:.1%}" if total else "无可比对",
        "PASS" if total and rate >= 0.8 else "WARN",
        note="spec §6: 三者方向一致才可写正文; 不一致时**以 mouse_only 为准**(单一物种、无映射混杂)")

    # ---- 4. 写盘 ----
    out_ma = ma_all[["gene", "species", "max_group", "max_fraction", "g_star", "peak_group_used",
                     "peak_bin", "peak_period", "n_feasible_bins", "module", "shape", "status",
                     "mech_class", "tier", "in_panel"]]
    out_ma = out_ma.merge(pr[["gene", "peak_rank_human_equiv", "peak_source", "peak_rank_human_only",
                              "peak_rank_mouse_only", "censored_group"]], on="gene", how="left")
    out_ma.to_csv(os.path.join(crossdir, "module_assignment.tsv"), sep="\t", index=False, float_format="%.6g")
    enr.to_csv(os.path.join(crossdir, "enrichment_class_x_module.tsv"), sep="\t", index=False, float_format="%.6g")
    fish.to_csv(os.path.join(crossdir, "fisher_class_x_shape.tsv"), sep="\t", index=False, float_format="%.6g")
    pr.to_csv(os.path.join(crossdir, "peak_rank_by_class.tsv"), sep="\t", index=False, float_format="%.6g")
    pairs_all.to_csv(os.path.join(crossdir, "peak_rank_pairwise.tsv"), sep="\t", index=False, float_format="%.6g")

    # ---- 5. 图 ----
    figdir = cfg.results("figures", "cross", mkdir=True)
    exts = list(cfg.raw["figures"]["formats"])
    stars = cfg.raw["figures"]["star_levels"]
    mech_colors = cfg.raw["figures"]["mech_colors"]

    def _save(fig, stem):
        for e in exts:
            fig.savefig(os.path.join(figdir, f"{stem}.{e}"), dpi=int(cfg.raw["figures"]["dpi"]))
        plt.close(fig)

    # fig3a: 机制类 × 模块 log2FE 热图(两物种并排)
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), squeeze=False)
    for k, sp in enumerate(("human", "mouse")):
        ax = axs[0][k]
        e = enr[(enr["species"] == sp) & (enr["test"] == "class_x_module")]
        if not len(e):
            ax.set_title(f"{sp}: 无富集结果")
            continue
        piv = e.pivot(index="a", columns="b", values="log2fe")
        pq = e.pivot(index="a", columns="b", values="q")
        vmax = float(np.nanmax(np.abs(piv.to_numpy()))) or 1.0
        im = ax.imshow(piv.to_numpy(), aspect="auto", cmap=cfg.raw["figures"]["diverging_cmap"],
                       vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(piv.shape[1]), [str(c).replace("|", "\n") for c in piv.columns], fontsize=5.5)
        ax.set_yticks(range(piv.shape[0]), piv.index, fontsize=7)
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                q = pq.to_numpy()[i, j]
                if np.isfinite(q):
                    lab = atxplots.star_of_q(q, cfg)   # 取满足的最严档(与 b3 一致; 原 next() 取最宽档)
                    if lab:
                        ax.text(j, i, lab, ha="center", va="center", fontsize=8)
        ax.set_title(f"[{sp}] fig3a mech_class × module log2FE\n(G*×P*; *q<0.05; 残差类已排除 CL-027)",
                     fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.8, label="log2FE")
    fig.tight_layout()
    _save(fig, "fig3a_enrichment_heatmap")

    # fig3b: 机制类 × shape 计数马赛克(以计数热图实现, 图注说明)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4), squeeze=False)
    for k, sp in enumerate(("human", "mouse")):
        ax = axs[0][k]
        sub = ma_all[(ma_all["species"] == sp) & (ma_all["status"] == "ok") & ma_all["in_panel"]
                     & ~ma_all["mech_class"].isin(residual)]
        ct = pd.crosstab(sub["mech_class"], sub["shape"])
        im = ax.imshow(ct.to_numpy(), aspect="auto", cmap=cfg.raw["figures"]["sequential_cmap"])
        ax.set_xticks(range(ct.shape[1]), ct.columns, rotation=30, ha="right", fontsize=6.5)
        ax.set_yticks(range(ct.shape[0]), ct.index, fontsize=7)
        for i in range(ct.shape[0]):
            for j in range(ct.shape[1]):
                v = ct.to_numpy()[i, j]
                if v:
                    ax.text(j, i, int(v), ha="center", va="center", fontsize=6,
                            color="w" if v > ct.to_numpy().max() * 0.6 else "k")
        ax.set_title(f"[{sp}] fig3b mech_class × shape 计数\n(Fisher 逐格见 fisher_class_x_shape.tsv)",
                     fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.8, label="genes")
    fig.tight_layout()
    _save(fig, "fig3b_shape_mosaic")

    # fig3c: 峰时 ridge(点图; 三套轴分面)
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    classes = sorted(prk["mech_class"].dropna().unique())
    for k, (axis, col) in enumerate(axes.items()):
        ax = axs[k]
        sub = prk[np.isfinite(prk[col])]
        rng = np.random.default_rng(cfg.seed)
        for i, c in enumerate(classes):
            v = sub.loc[sub["mech_class"] == c, col].to_numpy(float)
            if not len(v):
                continue
            ax.scatter(v, i + rng.uniform(-0.16, 0.16, len(v)), s=16,
                       c=mech_colors.get(c, "#7f7f7f"), edgecolors="none")
            ax.plot([np.median(v)], [i], marker="_", ms=22, mew=2, c="k")
        kw = float(peak_tab.loc[peak_tab["axis"] == axis, "kw_p"].iloc[0])
        ax.set_yticks(range(len(classes)), classes if k == 0 else [], fontsize=7)
        ax.set_xlabel("peak_rank_human_equiv" + ("" if axis == "human_only" else " (鼠峰映射: 部分插值)"),
                      fontsize=7.5)
        ax.set_title(f"fig3c 峰时轴 = {axis}\nn={len(sub)}, KW p={kw:.3g}"
                     + ("  ← 主分析" if axis == "mixed" else ""), fontsize=8)
        ax.grid(axis="x", lw=0.4, color="0.85")
    fig.suptitle("机制类 × 峰时(统一人等价秩): 三套轴必须方向一致才可写正文, 不一致以 mouse_only 为准(spec §6)",
                 fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, "fig3c_peak_ridge")

    # fig3d: 代表轨迹(取每物种最强显著格; 无显著格则取最大模块)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4), squeeze=False)
    for k, sp in enumerate(("human", "mouse")):
        ax = axs[0][k]
        e = enr[(enr["species"] == sp) & (enr["test"] == "class_x_module") & (enr["q"] < 0.05)
                & (enr["obs"] >= 3)]
        sub = ma_all[(ma_all["species"] == sp) & (ma_all["status"] == "ok") & ma_all["in_panel"]]
        if len(e):
            e = e.sort_values("log2fe", ascending=False).iloc[0]
            sel = sub[(sub["mech_class"] == e["a"]) & (sub["module"] == e["b"])]
            title = f"[{sp}] fig3d {e['a']} × {e['b']} (log2FE={e['log2fe']:.2f}, q={e['q']:.2g})"
        elif len(sub):
            top_mod = sub["module"].value_counts().index[0]
            sel = sub[sub["module"] == top_mod]
            title = f"[{sp}] fig3d 最大模块 {top_mod}(无显著富集格)"
        else:
            ax.set_title(f"[{sp}] 无可用模块")
            continue
        bt = pd.read_csv(cfg.results("intermediate", sp, "cube_bin.tsv"), sep="\t")
        bt = bt[bt["feasible"]].sort_values("bin_order")
        g_star = str(sel.iloc[0]["g_star"])
        grp = bt[bt["group"] == (g_star if g_star != "broad" else "ALL")]
        for gene in sel["gene"].head(12):
            col = f"log_{gene}"
            if col in grp.columns:
                ax.plot(grp["bin_order"], grp[col], "-o", ms=2.5, lw=1, label=gene)
        ax.set_title(title + f"\n组={g_star}(broad→ALL), 前 {min(12, len(sel))} 基因", fontsize=8)
        ax.set_xlabel("bin_order")
        ax.set_ylabel("mean log1p CP10K")
        ax.legend(fontsize=5.5, ncol=2, framealpha=0.6)
        ax.grid(lw=0.4, color="0.9")
    fig.tight_layout()
    _save(fig, "fig3d_representative_traj")

    rec("产物", "tables/cross + figures/cross",
        f"module_assignment({len(out_ma)} 行), enrichment({len(enr)}), fisher({len(fish)}), "
        f"peak_rank({len(pr)}), pairwise({len(pairs_all)}); 图 4 张 ×{len(exts)} 格式", "PASS")
    _finish(cfg, t0)
    return 1 if N_FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
