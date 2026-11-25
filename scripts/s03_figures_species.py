#!/usr/bin/env python
"""s03 — 物种内图形(implementation_spec §5 图产物;绘图原语 = atx/plots.py)。

用法:
    python -I scripts/s03_figures_species.py <root> --species mouse [--config PATH]

前置: s01(cube_bin.tsv + cube_donor.h5ad)与 s02(tables/<sp>/ 七 TSV)已产出。

视图(spec §5 命名 {species}_{view}[_{group}][_pN].{ext}, ext 集合 = figures.formats):
  a1  ALL 组趋势热图(行=基因按 ALL rho 降序, 列=bin, 色=基因内 z, 斜纹=不可行)
  b1  分组热图(行=bin, 列=基因, 色=组内基因 z;9 细胞组各一套, ALL 不重复出)
  b2  轨迹小倍数(donor 灰点/疾病△/均值线/sem 带限连续可行段/单供体▲/斜纹)
  b3  组×基因 rho 热图(星=q_bin 最严档, 白底斜纹=rho NA)
  b4  temporal_pattern 按组堆叠条(词表/色 = config pattern_colors)
  c1  9 组特异性 dotplot(行=主组→tau 序;大小=pct、色=log 均值、tau 边列)
  c2  发育期×组占比热图(顶部期色带 + n_groups 列头;NA 斜纹)
  d1  emergence 时间轴(行=emergence 先后;点=首个达标 bin, 色=主组;左截断开口点)

行序规则 figures.row_order_rule=rho_desc(其余值 NotImplementedError);
b1/b2/c1 的"9 组" = plot_groups ∩ 实有 − ALL(主 8 + 描述组, S2 补充图口径);
c1 的 pct/色值直接从 bin_table 可行 bin 求均值(读图指南 Fig1c: 汇集均值),
tau 边列取 specificity_global.tau。

build_figures() 为纯函数(outdir 显式传入)——dev_preview/测试复用, 不触真实 results/。
统一 footer(atx.plots.add_footer): draft 版本自动附 DRAFT 探索声明; R5 "部分插值"
标注由 s07/s08 经 notes 参数复用同一机制。
"""
from __future__ import annotations

import argparse
import math
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
from atx import plots as atxplots
import s02_stats

REPORT: list[dict] = []


def rec(check, observed, level, note=""):
    REPORT.append({"check": check, "observed": str(observed), "level": level, "note": note})
    mark = {"PASS": "✓", "WARN": "⚠", "FAIL": "✗"}.get(level, "?")
    print(f"[{mark} {level}] {check} | {observed}" + (f" | {note}" if note else ""))


def read_tables(indir: str) -> dict[str, pd.DataFrame]:
    """s02 七 TSV → DataFrame(keep_default_na=False: 受控 'NA' 串保留;空串 → NaN)。"""
    out = {}
    for name in s02_stats.TABLE_NAMES:
        p = os.path.join(indir, f"{name}.tsv")
        if not os.path.exists(p):
            raise RuntimeError(f"缺 s02 产物 {p}——先跑 s02(spec §9 顺序)")
        out[name] = pd.read_csv(p, sep="\t", keep_default_na=False, na_values=[""])
    return out


def _cell_groups(cfg, bt) -> list[str]:
    """9 细胞组 = plot_groups ∩ 实有 − ALL(b1/b2/c1 列集;ALL 由 a1 承载)。"""
    present = set(bt["group"].astype(str))
    return [g for g in cfg.plot_groups() if g in present and g != "ALL"]


def build_figures(tables: dict[str, pd.DataFrame], bt: pd.DataFrame, cube,
                  genes: list[str], spcfg, cfg, outdir: str) -> dict[str, list[str]]:
    """全部视图 → outdir(命名见模块 docstring)。返回 {view: [路径]}(纯函数, 不写日志)。"""
    fc = cfg.raw["figures"]
    if str(fc["row_order_rule"]) != "rho_desc":
        raise NotImplementedError(f"figures.row_order_rule={fc['row_order_rule']!r}: 仅实现 rho_desc")
    sp = spcfg.name
    labels = list(spcfg.labels)
    tr, tp = tables["trend_by_group"], tables["temporal_patterns"]
    sg, bp, em = tables["specificity_global"], tables["specificity_by_period"], tables["emergence"]
    groups_present = [g for g in cfg.plot_groups() if g in set(bt["group"].astype(str))]
    cell_groups = _cell_groups(cfg, bt)
    paths: dict[str, list[str]] = {}

    def save(view, fig, stem, notes=None):
        paths.setdefault(view, []).extend(atxplots.save_figures(fig, cfg, stem, outdir, notes=notes))

    # ---- a1: ALL 趋势热图(基因内 z;行 = ALL rho 降序, NaN 沉底) ----
    rho_all = (tr.loc[tr["group"].to_numpy() == "ALL"].set_index("gene")["rho_bin"]
               .reindex(genes).to_numpy(dtype=float) if (tr["group"] == "ALL").any()
               else np.full(len(genes), np.nan))
    gene_order = atxplots.order_rho_desc(rho_all, genes)
    M = atxplots.group_bin_matrix(bt, "ALL", gene_order, labels, "log_")
    Z = atxplots.zscore_rows(M)
    vmax = float(np.nanmax(np.abs(Z))) if np.isfinite(Z).any() else 1.0
    pages = atxplots.paged_heatmap(
        Z, gene_order, labels, cfg, vmin=-vmax, vmax=vmax, cmap_key="diverging_cmap",
        title=f"{sp} · a1 ALL-group trend heatmap — within-gene z of mean log1p CP10K, "
              f"rows: rho desc, hatch = infeasible/NA", ylabel="gene")
    for i, (fig, _ax, _sl) in enumerate(pages, 1):
        save("a1", fig, f"{sp}_a1_ALL_p{i}")

    # ---- b1: 分组热图(行=bin, 列=基因;组内基因 z;列序 = ALL rho 序, 组间可比) ----
    for g in cell_groups:
        Zg = atxplots.group_bin_matrix(bt, g, gene_order, labels, "log_")
        Zg = atxplots.zscore_rows(Zg).T                       # genes×bins z → bins×genes
        v = float(np.nanmax(np.abs(Zg))) if np.isfinite(Zg).any() else 1.0
        pages = atxplots.paged_heatmap(
            Zg, labels, gene_order, cfg, page_axis="cols", vmin=-v, vmax=v,
            title=f"{sp} · b1 {g} group heatmap — within-gene z across bins, hatch = infeasible/NA",
            xrot=45, ylabel="bin", xlabel="gene (ALL-rho order)")
        for i, (fig, _ax, _sl) in enumerate(pages, 1):
            save("b1", fig, f"{sp}_b1_{atxplots.safe_name(g)}_p{i}")

    # ---- b2: 轨迹小倍数(组内 rho 降序分页;donor 灰点/疾病△/sem 带/单供体▲) ----
    per = int(fc["traj_panels_per_page"])
    ncols = 3
    for g in cell_groups:
        rho_g = (tr.loc[tr["group"].to_numpy() == g].set_index("gene")["rho_bin"]
                 .reindex(genes).to_numpy(dtype=float))
        shown = [x for x in atxplots.order_rho_desc(rho_g, genes)
                 if f"log_{x}" in bt.columns
                 and np.isfinite(atxplots.group_bin_matrix(bt, g, [x], labels, "log_")).any()]
        npages = max(1, math.ceil(len(shown) / per))
        for pi in range(npages):
            chunk = shown[pi * per:(pi + 1) * per]
            if not chunk:
                continue
            fig, _infos = atxplots.traj_figure(
                bt, cube, chunk, g, cfg, spcfg, ncols=ncols,
                title=f"{sp} · b2 {g} trajectories — page {pi + 1}/{npages} "
                      f"(gray = donor, line = group mean ± sem band, ▲ = single-donor bin, "
                      f"open △ = disease donor, hatch = infeasible)")
            save("b2", fig, f"{sp}_b2_{atxplots.safe_name(g)}_p{pi + 1}")

    # ---- b3: 组×基因 rho 热图 + 星 ----
    pages = atxplots.rho_heatmap_pages(
        tr, gene_order, groups_present, cfg,
        title=f"{sp} · b3 group × gene Spearman rho — * q<0.05, ** q<0.01, *** q<0.001 "
              f"(BH within species×group), hatch = rho NA")
    for i, (fig, _ax, _sl) in enumerate(pages, 1):
        stem = f"{sp}_b3_p{i}" if len(pages) > 1 else f"{sp}_b3"
        save("b3", fig, stem)

    # ---- b4: temporal_pattern 堆叠条 ----
    fig = atxplots.pattern_bars_figure(
        tp, groups_present, cfg,
        title=f"{sp} · b4 temporal_pattern counts by group (controlled vocabulary, config colors)")
    save("b4", fig, f"{sp}_b4")

    # ---- c1: 9 组 dotplot(主组+描述组 ∩ 实有;行 = 主组→tau 序) ----
    dot_groups = [g for g in cell_groups
                  if g in set(cfg.main_groups()) | set(cfg.descriptive_groups())]
    c1_order = atxplots.order_main_tau(sg, cfg.main_groups())
    fig, _info = atxplots.dotplot_figure(
        bt, sg, c1_order, dot_groups, cfg, species=sp)
    save("c1", fig, f"{sp}_c1")

    # ---- c2: 期×组占比热图(行序 = c1 同序) ----
    pages = atxplots.period_fraction_pages(
        bp, c1_order, cfg,
        title=f"{sp} · c2 period × group mean fraction — denominator = feasible main groups "
              f"of period (n in header), hatch = NA (k < min_groups / missing)")
    for i, (fig, _ax, _sl) in enumerate(pages, 1):
        save("c2", fig, f"{sp}_c2_p{i}")

    # ---- d1: emergence 时间轴 ----
    fig, _info = atxplots.emergence_figure(em, cfg, spcfg, species=sp)
    save("d1", fig, f"{sp}_d1")

    return paths


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--species", required=True, choices=["human", "mouse"])
    ap.add_argument("--config", default=None)
    ap.add_argument("--out-subdir", default=None,
                    help="产物子目录 results/v{ver}/<subdir>/…; 须与 s01/s02 同值")
    ap.add_argument("--sens", action="append", default=None, metavar="KEY=VALUE",
                    help="进程内覆盖 config sensitivity 开关(可多次); 须与 s01/s02 同口径")
    args = ap.parse_args()
    sp = args.species
    t0 = time.time()
    cfg = atxconfig.load(args.config or os.path.join(os.path.abspath(args.root), "config", "config.yaml"))
    atxconfig.apply_sens_overrides(cfg, args.sens)      # 必须在 effective_species/results 之前
    if args.out_subdir:
        cfg.results_subdir = args.out_subdir.strip("/")
    cfg.check_version_gate()
    cfg.write_draft_marker()
    spcfg = cfg.effective_species(sp)          # 与 s01/s02 同口径(敏感性开关传导)

    print(f"== s03 {sp} | version={cfg.version}(draft={cfg.draft})"
          + (f" | 产物隔离: results/v{cfg.version}/{cfg.results_subdir}/" if cfg.results_subdir else ""))
    inter = cfg.results("intermediate", sp)
    bt_path = os.path.join(inter, "cube_bin.tsv")
    cube_path = os.path.join(inter, "cube_donor.h5ad")
    for p in (bt_path, cube_path):
        if not os.path.exists(p):
            raise RuntimeError(f"缺 s01 产物 {p}——先跑 s01(spec §9 顺序)")
    bt = pd.read_csv(bt_path, sep="\t")
    import anndata
    cube = anndata.read_h5ad(cube_path)
    tables = read_tables(cfg.results("tables", sp))

    panel = atxgenes.load_panel(cfg)
    genes = panel.loc[panel["in_stat_panel"].astype(bool), "gene_key"].astype(str).tolist()
    rec("基因行集(统计 panel)", f"{len(genes)} 基因(cube 命中 {sum(f'log_{g}' in bt.columns for g in genes)})", "PASS")

    outdir = cfg.results("figures", sp, mkdir=True)
    paths = build_figures(tables, bt, cube, genes, spcfg, cfg, outdir)
    for view in ("a1", "b1", "b2", "b3", "b4", "c1", "c2", "d1"):
        ps = paths.get(view, [])
        n_stem = len({os.path.splitext(p)[0] for p in ps})
        rec(f"{view} 视图", f"{n_stem} 文件 × {len(cfg.raw['figures']['formats'])} 格式", "PASS",
            note=", ".join(sorted({os.path.basename(p).rsplit('.', 1)[0] for p in ps}))[:300])

    logs = cfg.results("logs", mkdir=True)
    rp = os.path.join(logs, f"s03_report_{sp}.tsv")
    pd.DataFrame(REPORT, columns=["check", "observed", "level", "note"]).to_csv(rp, sep="\t", index=False)
    total = sum(len(v) for v in paths.values())
    print(f"\n== s03 {sp} 完成 ({time.time()-t0:.0f}s);{total} 个图形文件 → {outdir}")
    print(f"报告: {rp}")
    print("下一步: s06_report(物种报告+门清单)→ 人工 review 门(spec §9)。")
    return 1 if any(r["level"] == "FAIL" for r in REPORT) else 0


if __name__ == "__main__":
    sys.exit(main())
