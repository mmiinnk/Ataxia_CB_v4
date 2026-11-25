#!/usr/bin/env python
"""s00 — 输入与基因表门(run_all.sh gate 第一环;implementation_spec §3.12/§4/§9)。

用法:
    python -I scripts/s00_verify_inputs.py <root> [--config config/config.yaml] [--skip-dropout]

校验清单(级别):
  1. 版本门: config gene_set.version == curation_log.md 冻结版本(FAIL;硬约束 6)
  2. 基因表装载: 统计 panel(core & in_panel) ∪ sanity(报告)
  3. Sepp 人/鼠: 形状(FAIL) / X 计数整数性抽样(FAIL) / panel 命中率 vs min_hit_rate(FAIL) /
     size_factor ↔ X 行和锚点 r≈0.99997, 比值≈0.99(WARN) / bin×供体普查与
     single_donor_bins、pooled 校验(FAIL) / 9 组细胞数锚点 §4.1(FAIL, 精确)
  4. Aldinger/Castonguay: 形状 §4.3(FAIL) / raw/X 整数性 + X 浮点性(FAIL) /
     命中率(Castonguay 仅 ENSG, §11 未决 5)(FAIL) / 标签映射交叉表(未映射显式列出, §4.5F) /
     Castonguay 供体锚点(WARN)
  5. 鼠 dropout 现场重算 vs §3.12 锚点(WARN;--skip-dropout 跳过): pct+ 用
     bincount(X.indices) 单遍计算, 免疫组内 pct+ 用行子集单遍计算
  6. external/ 三个人工文件缺席 → WARN(§11 未决 6, 不阻断)

产物: results/v<ver>/logs/{s00_report.tsv, s00_gene_mapping_<ds>.tsv, s00_bin_census_<sp>.tsv,
s00_label_mapping_<ds>.tsv, s00_dropout.tsv} + DRAFT 标记(draft 版本时)。
退出码: 有 FAIL → 1。
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
from atx import groups as atxgroups
from atx import io as atxio
from atx import sanity as atxsanity

REPORT: list[dict] = []
N_FAIL = 0


def rec(check, dataset, expected, observed, level, note=""):
    global N_FAIL
    REPORT.append({"check": check, "dataset": dataset, "expected": str(expected),
                   "observed": str(observed), "level": level, "note": note})
    mark = {"PASS": "✓", "WARN": "⚠", "FAIL": "✗"}.get(level, "?")
    print(f"[{mark} {level}] {dataset or '-'}: {check} | 期望={expected} 实测={observed}"
          + (f" | {note}" if note else ""))
    if level == "FAIL":
        N_FAIL += 1


def pct_pos_by_gene(X, n_cells: int) -> np.ndarray:
    """全基因 pct+(%):单遍 bincount(X.indices)(O(nnz), 不做逐列切片)。"""
    pos = np.bincount(X.indices, minlength=X.shape[1])
    return pos / max(n_cells, 1) * 100.0


# ---------------------------------------------------------------- Sepp 人/鼠
def check_sepp(sp: str, cfg, panel, aliases):
    spcfg = cfg.species(sp)
    path = cfg.h5ad(spcfg.d["h5ad_key"])
    if not os.path.exists(path):
        rec("文件存在", sp, path, "缺失", "FAIL")
        return None, None
    rec("文件存在", sp, path, f"{os.path.getsize(path)/1e6:.0f} MB", "PASS")

    shape = atxio.h5_shape(path, "X")
    exp_shape = (spcfg.n_cells_expected, spcfg.n_vars_expected)
    rec("X 形状锚点", sp, exp_shape, shape, "PASS" if shape == exp_shape else "FAIL")

    obs, var = atxio.read_obs_var(path)

    # X = 原始计数(硬约束 4):抽样 200 行整数性
    try:
        rng = np.random.default_rng(cfg.seed)
        rows = rng.choice(shape[0], size=min(200, shape[0]), replace=False)
        m, _ = atxio.read_rows_csr(path, rows, layer="X", check_integer=True)
        rec("X 计数整数性(抽样200行)", sp, "全整数", "全整数", "PASS")
        del m
    except AssertionError as e:
        rec("X 计数整数性(抽样200行)", sp, "全整数", str(e), "FAIL")

    # 命中率(FAIL 级;缺失基因逐个列出——不静默)
    gm = atxgenes.map_to_var(panel, var, species=sp, dataset=sp, aliases=aliases)
    hr = atxgenes.hit_rate(gm, sp)
    thr = float(cfg.gene_set["min_hit_rate"][sp])
    miss = gm.loc[(gm["in_stat_panel"].astype(bool)) & (gm[f"status_{sp}"] != "matched"), "gene_key"].tolist()
    amb = gm.loc[gm[f"status_{sp}"] == "ambiguous", "gene_key"].tolist()
    rec("panel 命中率(统计 panel)", sp, f"≥{thr}", f"{hr:.4f}",
        "PASS" if hr >= thr else "FAIL",
        note=f"缺失={miss}" + (f"; ambiguous={amb}" if amb else ""))
    gm.to_csv(cfg.results("logs", f"s00_gene_mapping_{sp}.tsv", mkdir=True), sep="\t", index=False)

    # size_factor 锚点(WARN 级):r≈0.99997, X 行和比 size_factor 低 ≈1%
    sf = pd.to_numeric(obs["size_factor"], errors="coerce").to_numpy(dtype=np.float64)
    if not np.all(np.isfinite(sf)) or np.any(sf <= 0):
        rec("size_factor 有效性", sp, "全为正有限值", f"非正/非有限 {int(np.sum(~(np.isfinite(sf)&(sf>0))))}", "FAIL")
    else:
        rowsum = atxio.row_sums_csr(path, "X")
        r = float(np.corrcoef(rowsum, sf)[0, 1])
        ratio = float(np.median(rowsum / sf))
        a = atxsanity.SIZE_FACTOR_ANCHOR
        lo, hi = a["rowsum_over_sf_median_range"]
        ok = r >= a["corr_min"] and lo <= ratio <= hi
        rec("size_factor↔X行和锚点", sp,
            f"r≥{a['corr_min']}, 中位比值∈[{lo},{hi}]", f"r={r:.5f}, 比值={ratio:.4f}",
            "PASS" if ok else "WARN", note="锚点: r=0.99997, 行和低≈1%(作者基因过滤前总 UMI)")
        del rowsum

    # bin×供体普查 + single_donor_bins / pooled 校验
    tcol = spcfg.time_col
    raw_t = obs[tcol].astype("object").where(obs[tcol].notna(), atxgroups.NA).astype(str)
    lab = raw_t.map(spcfg.label_of)
    if lab.isna().any():
        rec("时间列词表", sp, "全部落在 time_order", f"未映射值 {sorted(set(raw_t[lab.isna()]))[:5]}", "FAIL")
        lab = lab.fillna("<unmapped>")
    don = obs[spcfg.donor_col].astype(str)
    census = (pd.DataFrame({"bin": lab, "donor": don})
              .groupby("bin")["donor"].agg(n_cells="size", n_donors="nunique",
                                           donors=lambda s: ";".join(sorted(set(s))))
              .reindex(spcfg.labels))
    census["bin_order"] = [spcfg.bin_order[b] for b in census.index]
    census[["n_cells", "n_donors"]] = census[["n_cells", "n_donors"]].fillna(0).astype(int)
    census["donors"] = census["donors"].fillna("")
    census.to_csv(cfg.results("logs", f"s00_bin_census_{sp}.tsv", mkdir=True), sep="\t")
    single = {b for b in spcfg.labels if int(census.loc[b, "n_donors"]) == 1}
    if sp == "human":
        rec("single_donor_bins 校验", sp, sorted(spcfg.single_donor_bins), sorted(single),
            "PASS" if single == spcfg.single_donor_bins else "FAIL")
    if sp == "mouse":
        pooled_ok = True
        notes = []
        for b in spcfg.labels:
            donors_b = set(census.loc[b, "donors"].split(";"))
            if b in spcfg.pooled and donors_b != {"pooled"}:
                pooled_ok = False
                notes.append(f"{b}: 期望全 pooled, 实测 {sorted(donors_b)}")
            if b not in spcfg.pooled and "pooled" in donors_b:
                pooled_ok = False
                notes.append(f"{b}: 非 pooled bin 出现 pooled 供体")
        rec("pooled bin 校验", sp, f"pooled={sorted(spcfg.pooled)} 且仅这些 bin", "一致" if pooled_ok else notes,
            "PASS" if pooled_ok else "FAIL")

    # 9 组细胞数锚点(FAIL 级, 精确;§4.1)——只用 obs, 不触矩阵
    d = atxgroups.derive_obs(obs, spcfg, cfg)
    a, _ = atxgroups.assign(d, cfg)
    counts = atxgroups.group_counts(a)
    bad = atxsanity.check_group_counts(counts, sp)
    rec("9 组细胞数锚点 §4.1", sp, "全部精确匹配", "匹配" if not bad else bad,
        "PASS" if not bad else "FAIL")
    return obs, (d, a)


# ---------------------------------------------------------------- 外部数据集
def _ext_matrix_checks(ds: str, path: str, shape_exp, cfg):
    shape = atxio.h5_shape(path, "X")
    rec("形状锚点 §4.3", ds, shape_exp, shape, "PASS" if shape == tuple(shape_exp) else "FAIL")
    rng = np.random.default_rng(cfg.seed)
    rows = rng.choice(shape[0], size=min(200, shape[0]), replace=False)
    try:
        m, _ = atxio.read_rows_csr(path, rows, layer="raw/X", check_integer=True)
        rec("raw/X 整数计数(抽样200行)", ds, "全整数", "全整数", "PASS")
        del m
    except (AssertionError, KeyError) as e:
        rec("raw/X 整数计数(抽样200行)", ds, "全整数", str(e), "FAIL")
    try:
        mx, _ = atxio.read_rows_csr(path, rows[:50], layer="X", check_integer=False)
        nonint = bool(np.any(mx.data != np.floor(mx.data)))
        rec("X 为 log 浮点层(硬约束 4)", ds, "含非整数值", f"非整数={nonint}",
            "PASS" if nonint else "WARN", note="若全整数则 X 可能不是 log 层, 需人工核")
        del mx
    except KeyError as e:
        rec("X 为 log 浮点层(硬约束 4)", ds, "X 存在", str(e), "FAIL")


def check_aldinger(cfg, panel, aliases):
    ds = "aldinger"
    path = cfg.path_for("aldinger_h5ad")
    if not os.path.exists(path):
        rec("文件存在", ds, path, "缺失", "FAIL")
        return
    rec("文件存在", ds, path, f"{os.path.getsize(path)/1e6:.0f} MB", "PASS")
    _ext_matrix_checks(ds, path, atxsanity.EXTERNAL_SHAPES[ds], cfg)
    obs, var = atxio.read_obs_var(path)
    gm = atxgenes.map_to_var(panel, var, species="human", dataset=ds, aliases=aliases)
    hr = atxgenes.hit_rate(gm, ds)
    thr = float(cfg.gene_set["min_hit_rate"][ds])
    miss = gm.loc[(gm["in_stat_panel"].astype(bool)) & (gm[f"status_{ds}"] != "matched"), "gene_key"].tolist()
    rec("panel 命中率", ds, f"≥{thr}", f"{hr:.4f}", "PASS" if hr >= thr else "FAIL", note=f"缺失={miss}")
    gm.to_csv(cfg.results("logs", f"s00_gene_mapping_{ds}.tsv", mkdir=True), sep="\t", index=False)

    # 标签映射交叉表(§4.5F:未映射标签显式列出, 不静默丢弃)
    ext = cfg.raw["external"]["aldinger"]
    gmap = dict(ext["group_map"])
    all_only = set(ext.get("all_only", []))
    lab = obs["fig_cell_type"].astype(str)
    mapped = lab.map(lambda x: gmap.get(x, "ALL_only" if x in all_only else "unmapped"))
    ct = (pd.DataFrame({"fig_cell_type": lab, "group": mapped,
                        "donor": obs["donor_id"].astype(str)})
          .groupby(["fig_cell_type", "group"])["donor"].agg(n_cells="size", n_donors="nunique")
          .reset_index())
    ct.to_csv(cfg.results("logs", f"s00_label_mapping_{ds}.tsv", mkdir=True), sep="\t", index=False)
    unmapped = sorted(set(ct.loc[ct["group"] == "unmapped", "fig_cell_type"]))
    rec("fig_cell_type→组映射完整性", ds, "unmapped=0", f"unmapped={unmapped}",
        "PASS" if not unmapped else "WARN")
    n_pc = int((lab == "H-PC").sum())
    rec("H-PC 细胞数锚点 §4.3", ds, 25724, n_pc, "PASS" if n_pc == 25724 else "WARN")
    rec("供体数锚点 §4.3", ds, 13, obs["donor_id"].nunique(),
        "PASS" if obs["donor_id"].nunique() == 13 else "WARN")


def check_castonguay(cfg, panel):
    ds = "castonguay"
    path = cfg.path_for("castonguay_h5ad")
    if not os.path.exists(path):
        rec("文件存在", ds, path, "缺失", "FAIL")
        return
    rec("文件存在", ds, path, f"{os.path.getsize(path)/1e9:.1f} GB", "PASS")
    _ext_matrix_checks(ds, path, atxsanity.EXTERNAL_SHAPES[ds], cfg)
    obs, var = atxio.read_obs_var(path)
    # §11 未决 5:只按 ENSG 匹配, 缺失记 missing, 不回退符号
    gm = atxgenes.map_to_var(panel, var, species="human", dataset=ds, ensg_only=True)
    hr = atxgenes.hit_rate(gm, ds)
    thr = float(cfg.gene_set["min_hit_rate"][ds])
    miss = gm.loc[(gm["in_stat_panel"].astype(bool)) & (gm[f"status_{ds}"] != "matched"), "gene_key"].tolist()
    rec("panel 命中率(仅 ENSG)", ds, f"≥{thr}", f"{hr:.4f}", "PASS" if hr >= thr else "FAIL",
        note=f"缺失={miss}(status=missing, 不静默回退符号)")
    gm.to_csv(cfg.results("logs", f"s00_gene_mapping_{ds}.tsv", mkdir=True), sep="\t", index=False)

    ext = cfg.raw["external"]["castonguay"]
    lab = obs[ext["cell_type_column"]].astype(str)
    gmap = dict(ext["other_groups"])
    gmap[ext["purkinje_label"]] = "Purkinje"
    mapped = lab.map(lambda x: gmap.get(x, "ALL_only"))
    ct = (pd.DataFrame({"author_cell_type": lab, "group": mapped,
                        "donor": obs["donor_id"].astype(str), "disease": obs["disease"].astype(str)})
          .groupby(["author_cell_type", "group"])["donor"].agg(n_cells="size", n_donors="nunique")
          .reset_index())
    ct.to_csv(cfg.results("logs", f"s00_label_mapping_{ds}.tsv", mkdir=True), sep="\t", index=False)
    known_all = {"Endocytes", "Pericytes"}   # §4.5C:只进 ALL
    unknown = sorted(set(lab[~lab.isin(set(gmap) | known_all)]))
    rec("author_cell_type→组映射完整性", ds, "除 Endocytes/Pericytes 外全映射", f"未识别={unknown}",
        "PASS" if not unknown else "WARN")

    # 供体锚点(§4.3, WARN 级):normal 93 供体 / 31089 PC / 89 供体 ≥30 PC / 中位 257
    normal = obs["disease"].astype(str).isin(ext["use_disease"])
    pc_norm = normal & (lab == ext["purkinje_label"])
    per_don = pd.Series(pc_norm.to_numpy()).groupby(obs["donor_id"].astype(str).to_numpy()).sum()
    per_don = per_don[per_don > 0]
    ge30 = per_don[per_don >= ext["min_pc_per_donor"]]
    a = atxsanity.CASTONGUAY_ANCHORS
    got = {"n_donors_normal": int(obs.loc[normal, "donor_id"].nunique()),
           "n_pc_normal": int(pc_norm.sum()),
           "n_donors_ge30pc": int(len(ge30)),
           "median_pc_ge30": float(ge30.median()) if len(ge30) else np.nan}
    bad = [f"{k}: 实测 {got[k]} != 锚点 {v}" for k, v in a.items() if got[k] != v]
    rec("normal 供体/PC 锚点 §4.3", ds, a, got, "PASS" if not bad else "WARN", note="; ".join(bad))
    n_et = int(obs.loc[obs["disease"].astype(str) == "essential tremor", "donor_id"].nunique())
    rec("ET 供体数(排除口径)", ds, 16, n_et, "PASS" if n_et == 16 else "WARN")


# ---------------------------------------------------------------- 鼠 dropout 重算
def check_mouse_dropout(cfg, X, cell_group, var):
    """§3.12:2026-10-14 实测锚点现场重算(WARN 级)。pct+ 用 bincount 单遍。"""
    syms = {s.lower(): i for i, s in enumerate(var["feature_name"].astype(str))}
    n = X.shape[0]
    pct_all = pct_pos_by_gene(X, n)
    immune_pos = np.flatnonzero((cell_group == "Immune").to_numpy())
    Xi = X[immune_pos]
    pct_immune = pct_pos_by_gene(Xi, len(immune_pos))
    rows = []
    for sym, anchor in atxsanity.MOUSE_DROPOUT_EXEMPT.items():
        j = syms.get(sym.lower())
        if j is None:
            rows.append({"gene": sym, "scope": "全库", "pct_pos": np.nan, "anchor": anchor,
                         "verdict": "ABSENT"})
            continue
        v = float(pct_all[j])
        rows.append({"gene": sym, "scope": "全库", "pct_pos": round(v, 4), "anchor": anchor,
                     "verdict": "PASS" if abs(v - anchor) <= 1.0 else "WARN",
                     "note": "全库 pct+<1% → 豁免成立" if v < 1.0 else "豁免前提(pct+<1%)不成立!"})
        vi = float(pct_immune[j])
        rows.append({"gene": sym, "scope": "Immune组", "pct_pos": round(vi, 4), "anchor": "≤4",
                     "verdict": "PASS" if vi <= 4.0 else "WARN"})
    for sym, anchor in atxsanity.MOUSE_IMMUNE_ANCHORS.items():
        j = syms.get(sym.lower())
        if j is None:
            rows.append({"gene": sym, "scope": "Immune组", "pct_pos": np.nan, "anchor": anchor,
                         "verdict": "ABSENT"})
            continue
        v = float(pct_immune[j])
        rows.append({"gene": sym, "scope": "Immune组", "pct_pos": round(v, 2), "anchor": anchor,
                     "verdict": "PASS" if abs(v - anchor) <= 8.0 else "WARN"})
    df = pd.DataFrame(rows)
    df.to_csv(cfg.results("logs", "s00_dropout.tsv", mkdir=True), sep="\t", index=False)
    n_warn = int((df["verdict"] == "WARN").sum())
    n_abs = int((df["verdict"] == "ABSENT").sum())
    rec("鼠 dropout 现场重算 §3.12", "mouse", f"锚点±容差内(WARN=0)", f"WARN={n_warn}, ABSENT={n_abs}",
        "PASS" if n_warn == 0 and n_abs == 0 else "WARN",
        note="详见 logs/s00_dropout.tsv")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--config", default=None)
    ap.add_argument("--skip-dropout", action="store_true", help="跳过鼠 dropout 重算(需全载 mouse.h5ad)")
    args = ap.parse_args()
    root = os.path.abspath(args.root)
    cfg_path = args.config or os.path.join(root, "config", "config.yaml")
    t0 = time.time()

    cfg = atxconfig.load(cfg_path)
    if os.path.abspath(cfg.root) != root:
        print(f"提示: config project_root={cfg.root} 与命令行 root={root} 不一致;以 config 为准")
    # 1. 版本门(FAIL 即停)
    try:
        fv = cfg.check_version_gate()
        rec("版本门(硬约束 6)", "genes", "config==curation_log", fv, "PASS",
            note="draft 版本 → 全部产物为探索性口径" if cfg.draft else "")
    except RuntimeError as e:
        rec("版本门(硬约束 6)", "genes", "config==curation_log", str(e), "FAIL")
        _finish(cfg, t0)
        return 1
    cfg.write_draft_marker()

    # 2. 基因表
    try:
        panel = atxgenes.load_panel(cfg)
        aliases = atxgenes.load_aliases(cfg)
        n_stat = int(panel["in_stat_panel"].sum())
        n_san = int((~panel["in_stat_panel"].astype(bool)).sum())
        rec("基因表装载", "genes", "core&in_panel ∪ sanity",
            f"统计panel={n_stat}, sanity_only={n_san}, 合计={len(panel)}", "PASS",
            note=f"version={cfg.version}")
    except Exception as e:
        rec("基因表装载", "genes", "可读且无重复 gene_key", repr(e), "FAIL")
        _finish(cfg, t0)
        return 1

    # 3. Sepp
    mouse_pack = None
    for sp in ("human", "mouse"):
        try:
            _, pack = check_sepp(sp, cfg, panel, aliases)
            if sp == "mouse":
                mouse_pack = pack
        except Exception as e:
            rec("Sepp 校验(异常)", sp, "无异常", repr(e), "FAIL")
            import traceback; traceback.print_exc()

    # 4. 外部数据集
    try:
        check_aldinger(cfg, panel, aliases)
    except Exception as e:
        rec("Aldinger 校验(异常)", "aldinger", "无异常", repr(e), "FAIL")
        import traceback; traceback.print_exc()
    try:
        check_castonguay(cfg, panel)
    except Exception as e:
        rec("Castonguay 校验(异常)", "castonguay", "无异常", repr(e), "FAIL")
        import traceback; traceback.print_exc()

    # 5. 鼠 dropout 重算(全载 mouse.h5ad;~8GB, spec §10)
    if args.skip_dropout:
        rec("鼠 dropout 现场重算", "mouse", "-", "--skip-dropout 跳过", "WARN")
    else:
        try:
            import anndata
            d, a = mouse_pack if mouse_pack else (None, None)
            if d is None:
                rec("鼠 dropout 现场重算", "mouse", "obs 派生可用", "鼠 obs 校验未通过", "FAIL")
            else:
                ad = atxio.read_full(cfg.h5ad(cfg.species("mouse").d["h5ad_key"]))
                check_mouse_dropout(cfg, ad.X, a["cell_group"], ad.var)
                del ad
        except Exception as e:
            rec("鼠 dropout 现场重算(异常)", "mouse", "无异常", repr(e), "FAIL")
            import traceback; traceback.print_exc()

    # 6. external/ 人工文件(§11 未决 6:缺席只 WARN)
    ext = cfg.raw["external"]
    for key in ("brainspan", "gtex", "model_mouse_deg"):
        f = ext[key]["file"] if isinstance(ext[key], dict) else ext[key]
        p = f if os.path.isabs(f) else os.path.join(cfg.root, f)
        rec(f"external 文件 {key}", "external", p, "存在" if os.path.exists(p) else "缺席",
            "PASS" if os.path.exists(p) else "WARN",
            note="" if os.path.exists(p) else "s10 对缺失文件只 WARN 并出占位 panel(§11 未决 6)")

    return _finish(cfg, t0)


def _finish(cfg, t0):
    logs = cfg.results("logs", mkdir=True)
    df = pd.DataFrame(REPORT, columns=["check", "dataset", "expected", "observed", "level", "note"])
    df.to_csv(os.path.join(logs, "s00_report.tsv"), sep="\t", index=False)
    n_pass = int((df["level"] == "PASS").sum())
    n_warn = int((df["level"] == "WARN").sum())
    print(f"\ns00 汇总: PASS={n_pass} WARN={n_warn} FAIL={N_FAIL}  ({time.time()-t0:.0f}s)")
    print(f"报告: {os.path.join(logs, 's00_report.tsv')}")
    if N_FAIL:
        print("s00 门失败(失败即停, spec §9)——修复后重跑;WARN 项须逐条看过再放行。")
    return 1 if N_FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
