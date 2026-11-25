#!/usr/bin/env python3
"""s04_figures_crossspecies.py — 第 5 层前半:跨物种 dev_state 状态剖面与对齐图(spec §4.5A/§3.10/§8.1)。

做什么:
  1. 用**共享 dev_state**(config `state_alignment.group_states`, spec §4.4)作人鼠共同横轴——
     这是"按细胞分化进度对齐", **不做任何人鼠年龄换算**(§4.5A 时间对齐第二档);
  2. 每物种算 ALL 与各组的状态剖面(mean log1p CP10K), 状态核数 < `state_alignment.min_cells` → **NA**
     (硬约束 5: 不插值、不填 0);
  3. 共识状态序 = 每物种内按状态平均 bin_order 排秩、两物种秩平均(`order_rule: consensus_mean_time_rank`);
  4. 逐基因×组的跨物种统计: rho_human / rho_mouse(各自 bin 轴趋势)、rho_sign_concordant、traj_r
     (共享可行状态上的 Pearson r, ≥3 个才算) → 供 s09 打保真分, 也直接产出 Fig5a–d(=本脚本 x1–x4);
  5. 两物种**绝不合并表达矩阵**; 热图物种内逐基因 z-score、独立 colorbar、**只比形状不比幅值**。

产物: tables/cross/{state_profiles, cross_species_gene_group}.tsv
     figures/cross/{x1_state_heatmap_ALL, x2_state_heatmap_groups, x3_paired_trajectories,
                    x4_rho_concordance}.{pdf,png}
     logs/s04_report.tsv
用法: <omics python> -I scripts/s04_figures_crossspecies.py <root>
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import anndata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atx import align as atxalign                        # noqa: E402
from atx import config as atxconfig                      # noqa: E402

plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC", "Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

REPORT: list[dict] = []
N_FAIL = 0
SHOW_GROUPS = ["Purkinje", "GranuleLineage", "OligoLineage", "Astro"]   # x2 展示的组(spec: 2–4 组)
ANCHOR_GENES = ["ITPR1", "MBP", "ATOH1"]                                # x3 固定展示的趋势锚点基因


def rec(check, expected, observed, level, note=""):
    global N_FAIL
    REPORT.append({"check": check, "expected": str(expected), "observed": str(observed),
                   "level": level, "note": note})
    mark = {"PASS": "✓", "WARN": "⚠", "FAIL": "✗"}.get(level, "?")
    print(f"[{mark} {level}] {check} | 期望={expected} 实测={observed}" + (f" | {note}" if note else ""))
    N_FAIL += level == "FAIL"


def _finish(cfg, t0):
    logs = cfg.results("logs", mkdir=True)
    pd.DataFrame(REPORT).to_csv(os.path.join(logs, "s04_report.tsv"), sep="\t", index=False)
    n_pass = sum(1 for r in REPORT if r["level"] == "PASS")
    n_warn = sum(1 for r in REPORT if r["level"] == "WARN")
    print(f"\n汇总: PASS={n_pass} WARN={n_warn} FAIL={N_FAIL}  "
          f"报告: {os.path.join(logs, 's04_report.tsv')}  ({time.time() - t0:.0f}s)")


def load_species(cfg, sp: str, genes: list[str] | None):
    """读 cells_panel.h5ad → (X_log, obs, gene_key 列表)。"""
    p = cfg.results("intermediate", sp, "cells_panel.h5ad")
    if not os.path.exists(p):
        raise SystemExit(f"[s04] 缺 {p} —— 先跑 s01 --species {sp}")
    ad = anndata.read_h5ad(p)
    gk = [str(g) for g in ad.var["gene_key"].to_numpy()]
    if "dev_state" not in ad.obs.columns:
        raise SystemExit(f"[s04] {sp} cells_panel.obs 无 dev_state —— 需用当前版 s01 重跑(跨物种对齐轴)")
    return ad.X, ad.obs, gk


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--config", default=None)
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
    print(f"== s04 cross | version={cfg.version}(draft={cfg.draft})"
          + (f" | 产物隔离: results/v{cfg.version}/{cfg.results_subdir}/" if cfg.results_subdir else ""))

    sa = cfg.raw["state_alignment"]
    min_cells_state = int(sa["min_cells"])
    group_states = {k: list(v) for k, v in (sa.get("group_states") or {}).items()}
    all_states = sorted({s for v in group_states.values() for s in v})
    main_groups = [g for g in cfg.main_groups() if g in group_states]

    # ---- 载入两物种 ----
    Xh, oh, gh = load_species(cfg, "human", None)
    Xm, om, gm = load_species(cfg, "mouse", None)
    genes_tab = pd.read_csv(cfg.paths["genes_table"], sep="\t")
    in_panel = genes_tab["in_panel"].astype(str).str.lower().eq("true")
    panel = [g for g in genes_tab.loc[(genes_tab["tier"] == "core") & in_panel, "hgnc_symbol"]]
    genes = [g for g in panel if g in gh and g in gm]
    miss = sorted(set(panel) - set(genes))
    rec("跨物种共有统计 panel", "core & in_panel 且两物种 cube 均命中",
        f"{len(genes)} 基因; 仅单物种命中 {len(miss)}: {miss}", "PASS" if genes else "FAIL",
        note="仅单物种命中的基因不进跨物种比较(记 NA, 不静默丢弃)")
    ih = {g: i for i, g in enumerate(gh)}
    im = {g: i for i, g in enumerate(gm)}
    cols_h, cols_m = [ih[g] for g in genes], [im[g] for g in genes]

    # ---- 状态剖面(ALL + 各组) ----
    prof: dict[str, dict[str, pd.DataFrame]] = {}
    meta: dict[str, pd.DataFrame] = {}
    long_rows = []
    for sp, (X, obs, cols) in (("human", (Xh, oh, cols_h)), ("mouse", (Xm, om, cols_m))):
        Xs = X[:, cols]
        group_mask_all = np.ones(len(obs), dtype=bool)
        p_all = atxalign.state_profiles(Xs, obs, genes, all_states, min_cells_state, group_mask_all)
        meta[sp] = atxalign.state_meta(obs, all_states)
        scopes = {"ALL": p_all}
        cg = obs["cell_group"].astype(str).to_numpy()
        for g in main_groups:
            mask = (cg == g)
            scopes[g] = atxalign.state_profiles(Xs, obs, genes, group_states[g], min_cells_state, mask)
        prof[sp] = scopes
        for scope, df in scopes.items():
            feas = {s: int(meta[sp].loc[s, "n_cells"]) if s in meta[sp].index else 0 for s in df.columns}
            for st in df.columns:
                nc = feas[st]
                for gene in df.index:
                    v = df.at[gene, st]
                    long_rows.append({"species": sp, "scope": scope, "gene": gene, "dev_state": st,
                                      "mean_log": (None if pd.isna(v) else float(v)),
                                      "n_cells": nc, "feasible": bool(nc >= min_cells_state),
                                      "rank_in_species": int(list(df.columns).index(st)) + 1})
    profiles = pd.DataFrame(long_rows)

    order = atxalign.consensus_order(prof["human"]["ALL"], prof["mouse"]["ALL"])
    profiles["consensus_rank"] = profiles["dev_state"].map({s: i + 1 for i, s in enumerate(order)})

    n_feas_h = int((prof["human"]["ALL"].notna()).sum(axis=1).gt(0).sum())
    feas_states = {sp: int(sum(1 for s in prof[sp]["ALL"].columns
                               if prof[sp]["ALL"][s].notna().any())) for sp in ("human", "mouse")}
    rec("ALL 状态剖面可行性", f"每状态 ≥{min_cells_state} 核",
        f"人 {feas_states['human']}/{len(prof['human']['ALL'].columns)} 状态可行, "
        f"鼠 {feas_states['mouse']}/{len(prof['mouse']['ALL'].columns)}; 共识序长度 {len(order)}",
        "PASS" if min(feas_states.values()) >= 3 else "WARN")

    # ---- 逐基因×组跨物种统计 ----
    tr = {sp: pd.read_csv(cfg.results("tables", sp, "trend_by_group.tsv"), sep="\t")
          .set_index(["gene", "group"])["rho_bin"] for sp in ("human", "mouse")}
    st = {sp: pd.read_csv(cfg.results("tables", sp, "trend_by_group.tsv"), sep="\t")
          .set_index(["gene", "group"])["status"] for sp in ("human", "mouse")}
    rows = []
    for gene in genes:
        for g in main_groups:
            states = [s for s in group_states[g] if s in prof["human"][g].columns
                      and s in prof["mouse"][g].columns]
            ph = prof["human"][g].loc[gene, states] if gene in prof["human"][g].index else None
            pm = prof["mouse"][g].loc[gene, states] if gene in prof["mouse"][g].index else None
            traj, n_shared = (atxalign._traj_r(ph, pm) if ph is not None and pm is not None
                              else (float("nan"), 0))
            rh = tr["human"].get((gene, g), np.nan) if st["human"].get((gene, g), "ok") == "ok" else np.nan
            rm = tr["mouse"].get((gene, g), np.nan) if st["mouse"].get((gene, g), "ok") == "ok" else np.nan
            rh = float(rh) if pd.notna(rh) else float("nan")
            rm = float(rm) if pd.notna(rm) else float("nan")
            if np.isfinite(rh) and np.isfinite(rm) and abs(rh) >= atxalign.MIN_ABS_RHO \
                    and abs(rm) >= atxalign.MIN_ABS_RHO:
                rs = float(np.sign(rh) == np.sign(rm))
            else:
                rs = float("nan")
            rows.append({"gene": gene, "group": g, "rho_human": rh, "rho_mouse": rm,
                         "rho_sign_concordant": rs, "traj_r": traj, "n_shared_states": n_shared})
    xg = pd.DataFrame(rows)

    crossdir = cfg.results("tables", "cross", mkdir=True)
    profiles.to_csv(os.path.join(crossdir, "state_profiles.tsv"), sep="\t", index=False, float_format="%.6g")
    xg.to_csv(os.path.join(crossdir, "cross_species_gene_group.tsv"), sep="\t", index=False, float_format="%.6g")

    # 同向率(仅两侧 |rho|≥0.2 的可判行)
    jud = xg[xg["rho_sign_concordant"].notna()]
    rate = float(jud["rho_sign_concordant"].mean()) if len(jud) else float("nan")
    rec("趋势方向同向率(组×基因)", "两侧 |rho_bin| ≥ 0.2 才可判(spec §8.1)",
        f"可判 {len(jud)}/{len(xg)} 行, 同向 {rate:.1%}", "PASS" if len(jud) else "WARN")
    n_traj = int(xg["traj_r"].notna().sum())
    rec("traj_r 可算率", f"≥{atxalign.MIN_SHARED_STATES} 个共享可行状态",
        f"{n_traj}/{len(xg)} 行, 中位 {xg['traj_r'].median():+.3f}" if n_traj else "0 行",
        "PASS" if n_traj else "WARN",
        note="人侧右删失组(Purkinje/DN/UBC/Progenitor 产后无数据)会拉低可算率")

    # ---- 图 ----
    figdir = cfg.results("figures", "cross", mkdir=True)
    exts = list(cfg.raw["figures"]["formats"])
    dpi = int(cfg.raw["figures"]["dpi"])
    gcolors = cfg.raw["figures"]["group_colors"]

    def _save(fig, stem):
        for e in exts:
            fig.savefig(os.path.join(figdir, f"{stem}.{e}"), dpi=dpi)
        plt.close(fig)

    def _heat(ax, dfh, dfm, order_cols, title):
        """人/鼠双联热图: 物种内逐基因 z-score, 独立 colorbar, 只比形状(§3.10)。"""
        cols = [s for s in order_cols if s in dfh.columns or s in dfm.columns]
        for k, (df, sp) in enumerate(((dfh, "human"), (dfm, "mouse"))):
            axx = ax[k]
            sub = df.reindex(columns=[c for c in cols if c in df.columns])
            Z = atxalign.zscore_rows(sub)
            cmap = plt.get_cmap(cfg.raw["figures"]["diverging_cmap"]).copy()
            cmap.set_bad(color=cfg.raw["figures"]["na_style"]["facecolor"])
            im = axx.imshow(np.ma.masked_invalid(Z.to_numpy()), aspect="auto", cmap=cmap, vmin=-2, vmax=2)
            axx.set_xticks(range(len(cols)),
                           [c if c in sub.columns else "·" for c in cols], rotation=90, fontsize=5)
            axx.set_yticks(range(len(Z.index)), Z.index, fontsize=5)
            axx.set_title(f"{sp}(物种内 z, 独立色标)", fontsize=7)
            axx.set_xlabel("共享 dev_state(共识序; · = 该物种无此状态列)", fontsize=6)
            fig.colorbar(im, ax=axx, shrink=0.7, label="z")
        ax[0].set_ylabel("gene", fontsize=7)
        ax[0].set_title(f"human — {title}", fontsize=7.5)

    # x1: ALL 热图
    fig, ax = plt.subplots(1, 2, figsize=(15, 9))
    _heat(ax, prof["human"]["ALL"], prof["mouse"]["ALL"], order, "fig5a/x1 ALL × 共享 dev_state")
    fig.suptitle("x1(=Fig5a) dev_state 对齐 ALL 热图: 用细胞分化状态代替日历年龄做人鼠比较\n"
                 f"白格 = 该状态核数 < {min_cells_state}(NA, 未插值)", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, "x1_state_heatmap_ALL")

    # x2: 分组热图
    fig, axs = plt.subplots(len(SHOW_GROUPS), 2, figsize=(14, 4.0 * len(SHOW_GROUPS)))
    for i, g in enumerate(SHOW_GROUPS):
        _heat(axs[i], prof["human"][g], prof["mouse"][g],
              [s for s in order if s in group_states[g]], f"x2 组={g}")
    fig.suptitle("x2(=Fig5b) 状态对齐分组热图(组内状态子集, config state_alignment.group_states)", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    _save(fig, "x2_state_heatmap_groups")

    # x3: 并排轨迹(锚点基因 + traj_r 最高的 3 个基因)
    top = xg.dropna(subset=["traj_r"]).nlargest(3, "traj_r")
    pick = []
    for g in ANCHOR_GENES:
        if g in set(genes):
            pick.append((g, "ALL"))
    for r in top.itertuples():
        if (r.gene, "ALL") not in pick:
            pick.append((r.gene, r.group))
    pick = pick[:6]
    fig, axs = plt.subplots(len(pick), 2, figsize=(11, 2.5 * len(pick)), squeeze=False)
    for i, (gene, note_g) in enumerate(pick):
        for k, sp in enumerate(("human", "mouse")):
            ax = axs[i][k]
            df = prof[sp]["ALL"]
            if gene not in df.index:
                ax.set_title(f"{gene} — {sp}: 不在 panel", fontsize=7)
                continue
            row = df.loc[gene]
            zz = atxalign.zscore_rows(df).loc[gene]
            xs = [order.index(s) for s in row.index if s in order and pd.notna(row[s])]
            ys = [zz[s] for s in row.index if s in order and pd.notna(row[s])]
            ax.plot(xs, ys, "-o", ms=3, lw=1.2, color="#4C72B0" if sp == "human" else "#C44E52")
            ax.set_xticks(range(len(order)), order, rotation=90, fontsize=4.5)
            ax.axhline(0, color="0.8", lw=0.6)
            ax.set_title(f"{gene} — {sp}" + (f"(traj_r 组={note_g})" if note_g != "ALL" else ""), fontsize=7)
            ax.set_ylabel("z", fontsize=6)
    fig.suptitle("x3(=Fig5c) 并排轨迹: 同一共识状态轴上的人/鼠形状比较(锚点基因 + traj_r 最高 3 基因)",
                 fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _save(fig, "x3_paired_trajectories")

    # x4: rho 一致性散点
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    both = xg.dropna(subset=["rho_human", "rho_mouse"])
    for g in main_groups:
        s = both[both["group"] == g]
        if len(s):
            ax.scatter(s["rho_human"], s["rho_mouse"], s=16, alpha=0.75,
                       c=gcolors.get(g, "#7f7f7f"), label=f"{g} (n={len(s)})", edgecolors="none")
    ax.axhline(0, color="0.6", lw=0.7)
    ax.axvline(0, color="0.6", lw=0.7)
    q = [(int(((both.rho_human > 0) & (both.rho_mouse > 0)).sum()), "++"),
         (int(((both.rho_human < 0) & (both.rho_mouse < 0)).sum()), "--"),
         (int(((both.rho_human > 0) & (both.rho_mouse < 0)).sum()), "+-"),
         (int(((both.rho_human < 0) & (both.rho_mouse > 0)).sum()), "-+")]
    ax.set_xlabel("rho_bin(human, 人 13 bin 轴)")
    ax.set_ylabel("rho_bin(mouse, 鼠 12 bin 轴)")
    ax.set_title(f"x4(=Fig5d) 人鼠趋势方向一致性: 可判 {len(jud)} 行(两侧|rho|≥0.2), 同向率 {rate:.1%}\n"
                 f"象限计数 {dict(q)}; 两物种各用自己的 bin 轴(不做年龄换算)", fontsize=8.5)
    ax.legend(fontsize=6.5, markerscale=1.6, framealpha=0.7)
    fig.tight_layout()
    _save(fig, "x4_rho_concordance")

    rec("产物", "tables/cross + figures/cross",
        f"state_profiles({len(profiles)} 行), cross_species_gene_group({len(xg)} 行); "
        f"图 x1–x4 ×{len(exts)} 格式", "PASS")
    rec("Fig5 对应关系", "spec §8.1: fig5a/b/c 用 s04 x1–x3, d 用 x4",
        "x1→fig5a, x2→fig5b, x3→fig5c, x4→fig5d; fig5e(保真度表热图)由 s09 产出", "PASS")
    _finish(cfg, t0)
    return 1 if N_FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
