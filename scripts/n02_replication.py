"""n02 — 模块独立复核与冻结(Aldinger 胎儿复核 + Castonguay 老年检出 + 鼠同源可映射性)。

冻结规则(module_freeze.tsv 模板契约;config network.coexpression.lodo_min_preservation 与
network.replication.*):
  门1 donor_resample_stability   n01 LODO 均值 Jaccard ≥ lodo_min_preservation(0.6);
  门2 aldinger_fetal_replication 模块评分在 Aldinger H-PC(独立胎儿队列, 9–21 wpc) donor×bin
      轨迹方向与 Sepp eigengene 轨迹方向一致(两侧 |Spearman ρ| ≥ min_abs_rho_traj 且同号);
      任一侧不可判(|ρ|<阈或基因覆盖<50%) → indeterminate, 不冻结;
  门3 castonguay_mature_detectability 成员基因在老年 PC **供体级**伪汇总(normal 供体, ≥30 PC,
      每供体 ≤300 核降采样, raw/X 行切片, 先供体内 CP10K 再跨供体均值, spec §8.2(d) 口径)
      检出(均值 CP10K>0)比例 ≥ castonguay_min_detectability。
  三门全过 → frozen;按 LODO 稳定度排序, 至多 module_cap(3) 个;未过者逐项记录原因。

产物(results/v{ver}/network/): module_freeze.tsv(模板契约), replication_details.tsv,
  replication_summary.md, logs/n02_report.tsv。
口径: 数据集间不合并表达矩阵, 只在模块评分/检出层面比较(spec §4.5 总原则);
外部队列无 size_factor → CP10K 用 raw/X 行和(spec §3.5 差异②)。
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from atx import config as atxconfig            # noqa: E402
from atx import io as atxio                    # noqa: E402

REPORT: list[dict] = []
N_FAIL = 0


def rec(check, expected, observed, level, note=""):
    global N_FAIL
    if level == "FAIL":
        N_FAIL += 1
    REPORT.append({"check": check, "expected": expected, "observed": str(observed),
                   "level": level, "note": note})
    mark = {"PASS": "✓", "WARN": "⚠", "FAIL": "✗"}.get(level, "?")
    print(f"[{mark} {level}] {check} | 期望={expected} 实测={observed}" + (f" | {note}" if note else ""))


def resolve_col(obs, candidates, dataset):
    for c in candidates:
        if c in obs.columns:
            return c
    raise SystemExit(f"[n02] {dataset} obs 缺供体列(尝试过 {candidates}); 实际列: {list(obs.columns)[:24]}")


def zscore_cols(A: np.ndarray) -> np.ndarray:
    mu, sd = A.mean(axis=0), A.std(axis=0)
    return (A - mu) / np.where(sd > 0, sd, 1.0)


def aldinger_pseudobulk(cfg):
    """Aldinger H-PC → donor×bin(≥30 核) log1p CP10K(raw/X 行和为 sf)。返回 (矩阵, meta, var_index)。"""
    acfg = cfg.raw["external"]["aldinger"]
    ad = atxio.read_full(cfg.h5ad("aldinger_h5ad"))
    obs = ad.obs
    donor_col = resolve_col(obs, ["donor_id", "donor", "sample_id"], "aldinger")
    pc = (obs[acfg["cell_type_column"]].astype(str) == "H-PC").to_numpy()
    stage2bin = {s: b for b, sts in acfg["time_bins"].items() for s in sts}
    bin_ser = obs[acfg["time_column"]].astype(str).map(stage2bin)
    keep = pc & bin_ser.notna().to_numpy()
    pos = np.flatnonzero(keep)
    if ad.raw is None:
        raise SystemExit("[n02] Aldinger raw 缺失(raw/X 为计数层, 硬约束 4)")
    X = ad.raw.X
    var_index = np.asarray(ad.raw.var.index).astype(str)
    rowsum = np.asarray(X[pos].sum(axis=1)).ravel()          # sf = raw/X 行和(spec §3.5)
    donor_raw = obs[donor_col].astype("object").where(obs[donor_col].notna(), "__NA__").astype(str)
    frame = pd.DataFrame({"bin": bin_ser.to_numpy()[keep],
                          "donor": donor_raw.to_numpy()[keep],
                          "local": np.arange(len(pos))})
    mats, keys = [], []
    for (b, d), g in frame.groupby(["bin", "donor"], observed=True, sort=True):
        loc = g["local"].to_numpy()
        if len(loc) < 30:
            continue
        sub = X[pos[loc]]
        cp = np.asarray(sub.sum(axis=0)).ravel() / max(float(rowsum[loc].sum()), 1e-9) * 1e4
        mats.append(np.log1p(cp))
        keys.append({"bin": b, "donor": d, "n_cells": len(loc)})
    del ad, X
    bins_order = list(cfg.raw["external"]["aldinger"]["time_bins"].keys())
    return np.vstack(mats), pd.DataFrame(keys), var_index, bins_order


def castonguay_aged_pc(cfg):
    """Castonguay normal PC → 供体级 CP10K → 跨供体均值(spec §8.2(d) 口径)。

    返回 (cp_mean, per_donor_detect, var_index, n_donors, n_rows)。
    per_donor_detect[g] = 该基因在多少比例供体中检出(CP10K>0) —— 比"跨供体均值>0"判别力强
    (均值>0 在 ~2 万核上几乎恒真, finding #2)。
    队列过滤读 config external.castonguay.use_disease(硬约束 6), 不硬编码 "normal"(finding #5)。
    """
    kcfg = cfg.raw["external"]["castonguay"]
    path = cfg.h5ad("castonguay_h5ad")
    obs, _ = atxio.read_obs_var(path)
    donor_col = resolve_col(obs, ["donor_id", "donor", "sample_id"], "castonguay")
    use_dis = [str(x) for x in kcfg.get("use_disease", ["normal"])]
    # donor/disease NA 哨兵化, 防 groupby/比较静默丢行(finding #7)
    dis = obs["disease"].astype("object").where(obs["disease"].notna(), "__NA__").astype(str)
    dser = obs[donor_col].astype("object").where(obs[donor_col].notna(), "__NA__").astype(str)
    sel = (dis.isin(use_dis)
           & (obs[kcfg["cell_type_column"]].astype(str) == kcfg["purkinje_label"])).to_numpy()
    pos = np.flatnonzero(sel)
    donors = dser.to_numpy()[pos]
    rng = np.random.default_rng(cfg.seed)
    picks = []
    for d in pd.unique(donors):
        if d == "__NA__":
            continue
        idx = pos[donors == d]
        if len(idx) >= int(kcfg["min_pc_per_donor"]):
            picks.append(np.sort(rng.choice(idx, min(int(kcfg["subsample_pc_per_donor"]), len(idx)),
                                            replace=False)))
    if not picks:
        raise SystemExit(f"[n02] Castonguay 无满足 min_pc_per_donor 的 {use_dis} 供体")
    rows = np.concatenate(picks)
    Xc, var_index = atxio.read_rows_csr(path, rows, layer=kcfg["counts_layer"])
    var_index = np.asarray(var_index).astype(str)
    if len(var_index) != Xc.shape[1]:                        # finding #15: raw/X 列数须与 var 一致
        raise AssertionError(f"Castonguay var 长度 {len(var_index)} != 矩阵列数 {Xc.shape[1]}")
    donor_cps, off = [], 0
    for p in picks:                                          # picks 顺序 = rows 分段顺序
        sub = Xc[off:off + len(p)]
        off += len(p)
        rs = np.asarray(sub.sum(axis=1)).ravel()
        donor_cps.append(np.asarray(sub.sum(axis=0)).ravel() / max(float(rs.sum()), 1e-9) * 1e4)
    donor_mat = np.vstack(donor_cps)                         # n_donors × n_genes
    cp_mean = donor_mat.mean(axis=0)
    per_donor_detect = (donor_mat > 0).mean(axis=0)          # 每基因供体检出比例
    return cp_mean, per_donor_detect, var_index, len(picks), len(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--config", default=None)
    ap.add_argument("--skip-castonguay", action="store_true", help="仅调试; 正式冻结必须跑门3")
    args = ap.parse_args()
    root = os.path.abspath(args.root)
    cfg = atxconfig.load(args.config or os.path.join(root, "config", "config.yaml"))
    cfg.check_version_gate()
    t0 = time.time()
    netdir = cfg.results("network", mkdir=True)
    logs = cfg.results("logs", mkdir=True)
    cc = cfg.raw["network"]["coexpression"]
    rc = cfg.raw["network"].get("replication", {})
    min_rho = float(rc.get("min_abs_rho_traj", 0.2))
    min_detect = float(rc.get("castonguay_min_detectability", 0.5))
    module_cap = int(cfg.raw["network"]["module_cap"])

    summary = pd.read_csv(os.path.join(netdir, "module_summary.tsv"), sep="\t")
    modules = pd.read_csv(os.path.join(netdir, "modules.tsv"), sep="\t")
    cand = summary[summary.frozen_candidate].sort_values("lodo_jaccard", ascending=False)
    rec("冻结候选(门1 LODO 已过)", f"≤{module_cap} 进复核",
        f"{len(cand)} 个: " + ",".join(f"M{m}" for m in cand.module), "PASS")
    mod_genes = {int(m): modules.loc[modules.module == m, "ensg"].astype(str).tolist()
                 for m in summary.module}

    # ================= 门2: Aldinger 胎儿独立复核 =================
    A_log, a_meta, a_var, bins_order = aldinger_pseudobulk(cfg)
    a_col = {g: j for j, g in enumerate(a_var)}
    n_ald_donors = int(a_meta.donor.nunique())
    rec("Aldinger donor×bin 伪汇总", "≥30 核/组(H-PC)",
        f"{len(a_meta)} 组 / {n_ald_donors} 供体 / 基因 {len(a_var)}",
        "PASS" if len(a_meta) >= 8 else "WARN",
        note="组数少 → 轨迹 Spearman 自由度低, 判读须谨慎(finding #3)")
    # 每 bin 的供体数(轨迹每点的支撑; finding #3 透明化)
    donors_per_bin = a_meta.groupby("bin", observed=True)["donor"].nunique().to_dict()
    ald_res = {}
    for mid in sorted(mod_genes):
        gs = [g for g in mod_genes[mid] if g in a_col]
        cover = len(gs) / max(len(mod_genes[mid]), 1)
        rho_sepp = float(summary.loc[summary.module == mid, "traj_rho_stage"].iloc[0])
        if cover < 0.5 or len(gs) < 5:
            ald_res[mid] = {"coverage": cover, "rho_ald": np.nan, "rho_sepp": rho_sepp,
                            "n_bins": 0, "verdict": "indeterminate(基因覆盖不足)"}
            continue
        score = zscore_cols(A_log[:, [a_col[g] for g in gs]]).mean(axis=1)
        per_bin = (pd.Series(score).groupby(a_meta["bin"].to_numpy(), observed=True).mean()
                   .reindex([b for b in bins_order if b in set(a_meta["bin"])]).dropna())
        n_bins = len(per_bin)
        rho_ald = (spearmanr(np.arange(n_bins, dtype=float), per_bin.to_numpy()).statistic
                   if n_bins >= 4 else np.nan)
        if not np.isfinite(rho_ald) or not np.isfinite(rho_sepp) \
                or abs(rho_ald) < min_rho or abs(rho_sepp) < min_rho:
            verdict = "indeterminate(方向不可判)"
        elif np.sign(rho_ald) == np.sign(rho_sepp):
            verdict = "pass"
        else:
            verdict = "fail(方向相反)"
        ald_res[mid] = {"coverage": cover, "rho_ald": rho_ald, "rho_sepp": rho_sepp,
                        "n_bins": n_bins, "verdict": verdict}
    for mid, r in sorted(ald_res.items()):
        rho_a = r["rho_ald"]
        rec(f"门2 Aldinger 复核 M{mid}", f"方向一致(两侧|ρ|≥{min_rho} 且同号)",
            f"ρ_sepp={r['rho_sepp']:+.2f}, ρ_ald={'NA' if not np.isfinite(rho_a) else f'{rho_a:+.2f}'}, "
            f"基因覆盖={r['coverage']:.0%}, Aldinger轨迹点数={r['n_bins']} → {r['verdict']}",
            "PASS" if r["verdict"] == "pass" else "WARN",
            note=f"每bin供体数={donors_per_bin}")
    del A_log

    # ================= 门3: Castonguay 老年 PC 检出 =================
    cast_res = {}
    if args.skip_castonguay:
        rec("门3 Castonguay", "正式冻结必须跑", "--skip-castonguay(调试)", "FAIL",
            note="调试模式不评估门3, 且**不覆盖**已有 module_freeze.tsv(finding #11)")
    else:
        cp_mean, per_donor_det, c_var, n_donors_c, n_rows_c = castonguay_aged_pc(cfg)
        c_col = {g: j for j, g in enumerate(c_var)}
        rec("Castonguay PC 供体级伪汇总", "≥30 PC/供体, ≤300 核/供体, 先供体再跨供体均值",
            f"{n_donors_c} 供体 / {n_rows_c} 核", "PASS")
        for mid in sorted(mod_genes):
            gs = [g for g in mod_genes[mid] if g in c_col]
            cover = len(gs) / max(len(mod_genes[mid]), 1)
            # 供体级检出率(判别力强): 基因在 ≥50% 供体检出才算"可靠检出"; 模块得分=可靠检出基因占比
            det_per_donor = (float(np.mean([per_donor_det[c_col[g]] >= 0.5 for g in gs]))
                             if gs else 0.0)
            det_meanpos = float(np.mean([cp_mean[c_col[g]] > 0 for g in gs])) if gs else 0.0
            cast_res[mid] = {"coverage": cover, "detectability": det_per_donor,
                             "detectability_meanpos": det_meanpos,
                             "median_cp10k": float(np.median([cp_mean[c_col[g]] for g in gs])) if gs else np.nan}
        rec("门3 Castonguay 老年检出", f"成员基因供体级检出率(≥50%供体) ≥{min_detect}",
            ", ".join(f"M{m}={r['detectability']:.0%}(覆盖{r['coverage']:.0%})"
                      for m, r in sorted(cast_res.items())),
            "PASS" if cast_res and all(r["detectability"] >= min_detect for r in cast_res.values()) else "WARN",
            note="改用供体级检出率(均值>0 在~2万核上恒真, 无判别力; finding #2)")
        del cp_mean, per_donor_det

    # ================= 鼠同源可映射性(冻结表参考列; 符号重叠近似, finding #12) =================
    _, var_m = atxio.read_obs_var(cfg.h5ad("mouse_h5ad"))
    mouse_syms = {str(s).upper() for s in np.asarray(var_m["feature_name"])}
    vartab = pd.read_csv(os.path.join(netdir, "metacells_var.tsv"), sep="\t")
    sym_of = dict(zip(vartab.ensg.astype(str), vartab.symbol.astype(str)))
    # 分母 = 全部成员基因; 无符号或符号不在鼠 var → 记不可映射(保守); 这是**符号重叠近似**,
    # 非 ENSG↔ENSMUSG 同源(后者待 v1.1 BioMart, CL-017)。冻结表 notes 标注此口径。
    mouse_frac = {}
    for mid, gs in mod_genes.items():
        if not gs:
            mouse_frac[mid] = float("nan")
            continue
        mouse_frac[mid] = float(np.mean([sym_of.get(g, "").upper() in mouse_syms for g in gs]))

    # ================= 冻结 =================
    frozen, details = [], []
    for _, r in cand.iterrows():
        mid = int(r.module)
        a = ald_res.get(mid, {"verdict": "missing", "rho_ald": np.nan, "n_bins": 0, "coverage": np.nan})
        c = cast_res.get(mid, {"detectability": np.nan, "coverage": np.nan})
        passed = (a.get("verdict") == "pass") and (float(c.get("detectability", 0) or 0) >= min_detect)
        is_frozen = bool(passed and len(frozen) < module_cap)
        details.append({"module": f"M{mid}", "n_genes": int(r.n_genes),
                        "lodo_jaccard": r.lodo_jaccard,
                        "bootstrap_jaccard": r.get("bootstrap_jaccard", np.nan),
                        "traj_rho_stage": r.traj_rho_stage,
                        "aldinger_verdict": a.get("verdict"), "aldinger_rho": a.get("rho_ald"),
                        "aldinger_n_bins": a.get("n_bins"), "aldinger_coverage": a.get("coverage"),
                        "castonguay_detectability": c.get("detectability"),
                        "castonguay_coverage": c.get("coverage"),
                        "mouse_symbol_overlap_frac": mouse_frac.get(mid),
                        "frozen": is_frozen})
        if is_frozen:
            rho_a = a.get("rho_ald")
            frozen.append({"module_id": f"M{mid}", "version": f"v{cfg.version}",
                           "locked_at": time.strftime("%Y-%m-%d"),
                           "members_file": f"results/v{cfg.version}/network/module_M{mid}_members.tsv",
                           "weights_file_hash": r.members_sha256,
                           "scoring_direction": ("up_with_stage" if r.traj_rho_stage >= 0
                                                 else "down_with_stage"),
                           "missing_gene_rule": "评分=可测成员基因 z 均值; 可测成员<50% → NA",
                           "dropout_min_reliable_coverage": 0.5,
                           "donor_resample_stability": round(float(r.lodo_jaccard), 3),
                           "aldinger_fetal_replication": ("pass(rho=%s,n_bins=%d)"
                                                          % ("NA" if not np.isfinite(rho_a) else f"{rho_a:+.2f}",
                                                             int(a.get("n_bins", 0)))),
                           "castonguay_mature_detectability": round(float(c.get("detectability")), 3),
                           "mouse_ortholog_mappable_fraction": round(mouse_frac.get(mid, np.nan), 3),
                           "n_modules_cap": module_cap,
                           "notes": (f"hubs={str(r.hub_genes)[:90]}; mouse_frac=符号重叠近似(非ENSG同源);"
                                     f" cast_detect=供体级≥50%检出占比; ald_n_bins={a.get('n_bins')}")})
    detail_cols = ["module", "n_genes", "lodo_jaccard", "bootstrap_jaccard", "traj_rho_stage",
                   "aldinger_verdict", "aldinger_rho", "aldinger_n_bins", "aldinger_coverage",
                   "castonguay_detectability", "castonguay_coverage", "mouse_symbol_overlap_frac",
                   "frozen"]
    pd.DataFrame(details, columns=detail_cols).to_csv(
        os.path.join(netdir, "replication_details.tsv"), sep="\t", index=False)
    tmpl_cols = pd.read_csv(os.path.join(root, "templates", "module_freeze.tsv"),
                            sep="\t").columns.tolist()
    # finding #11: --skip-castonguay(门3 未评估)时不得用空表覆盖已有冻结契约
    if not args.skip_castonguay:
        pd.DataFrame(frozen, columns=tmpl_cols).to_csv(
            os.path.join(netdir, "module_freeze.tsv"), sep="\t", index=False)
        rec("模块冻结", f"三门全过者 ≤{module_cap}",
            f"{len(frozen)} 个冻结: {[f['module_id'] for f in frozen]}",
            "PASS" if len(frozen) else "WARN",
            note="未冻结候选逐项原因见 replication_details.tsv")
    else:
        rec("模块冻结", "门3 未评估时不写冻结表", "跳过写入(保留已有 module_freeze.tsv)", "WARN")

    def _f(v, fmt="{:.3f}"):
        try:
            return "NA" if v is None or pd.isna(v) else fmt.format(float(v))
        except (TypeError, ValueError):
            return "NA"

    with open(os.path.join(netdir, "replication_summary.md"), "w", encoding="utf-8") as fh:
        fh.write("# n02 模块复核与冻结摘要\n\n")
        fh.write(f"冻结门: ①LODO≥{cc['lodo_min_preservation']} ②Aldinger 方向一致(两侧|ρ|≥{min_rho}) "
                 f"③Castonguay 检出≥{min_detect}; 上限 {module_cap} 个。\n\n")
        fh.write("| 模块 | 基因数 | LODO | bootstrap | Sepp轨迹ρ | Aldinger | ρ_ald | Ald点数 | Cast供体检出 | Cast覆盖 | 鼠符号重叠 | 冻结 |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|---|---|---|\n")
        for d in details:
            mark = "✅" if d["frozen"] else "—"
            fh.write(f"| {d['module']} | {d['n_genes']} | {_f(d['lodo_jaccard'])} | "
                     f"{_f(d['bootstrap_jaccard'])} | {_f(d['traj_rho_stage'], '{:+.2f}')} | "
                     f"{d['aldinger_verdict']} | {_f(d['aldinger_rho'], '{:+.2f}')} | "
                     f"{d.get('aldinger_n_bins', 'NA')} | "
                     f"{_f(d['castonguay_detectability'], '{:.2f}')} | "
                     f"{_f(d['castonguay_coverage'], '{:.2f}')} | "
                     f"{_f(d['mouse_symbol_overlap_frac'], '{:.2f}')} | {mark} |\n")
    pd.DataFrame(REPORT).to_csv(os.path.join(logs, "n02_report.tsv"), sep="\t", index=False)
    print(f"\n汇总: PASS={sum(1 for r in REPORT if r['level'] == 'PASS')} "
          f"WARN={sum(1 for r in REPORT if r['level'] == 'WARN')} FAIL={N_FAIL}  ({time.time() - t0:.0f}s)")
    return 1 if N_FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
