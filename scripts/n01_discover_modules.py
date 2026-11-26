"""n01 — 正常 PC 全转录组共表达模块发现(config network.coexpression;阶段 1)。

流程: n00 metacells → 基因过滤(检出率≥min_detected_frac + 方差前 n_top_var_genes) →
Spearman(metacell 间) → signed weighted 邻接 a=((1+r)/2)^β(β 按 scale-free 判据选) →
TOM → 1−TOM average 层次聚类 → cut_height 定高切割 → <min_module_size 并入 grey →
模块特征基因(eigengene)与 kME → 稳定性: LODO(留一供体重建, 冻结门) + 供体 bootstrap(补充)。

口径注(硬约束 10):
  · config 原注释 "adjacency=|cor|^β" 与方法名 signed_weighted 矛盾; 本实现按 WGCNA signed
    网络约定 a=((1+cor)/2)^β(保留符号, 生物学标准做法); config 注释已同步标注, 阶段门待 ratify。
  · TOM_ij = (l_ij + a_ij)/(min(k_i,k_j) + 1 − a_ij), l = A², k = ΣA − 1(WGCNA 标准式)。
  · SFT 判据: 连通度量分位分箱, log10(p)~log10(k) 线性拟合, 取 R²≥sft_r2 且斜率<0 的最小 β;
    无满足 → max R² 的 β 并记 WARN。
  · metacell 样本量小(~2 位数), 相关噪声大 → LODO 保留度是冻结门槛而非参考值。

产物(results/v{ver}/network/): modules.tsv, module_summary.tsv, module_eigengenes.tsv,
  sft_beta.tsv, lodo_stability.tsv, bootstrap_stability.tsv, module_M*_members.tsv,
  logs/n01_report.tsv
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata, spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from atx import config as atxconfig            # noqa: E402

REPORT: list[dict] = []
N_FAIL = 0


def rec(check: str, expected: str, observed, level: str, note: str = "") -> None:
    global N_FAIL
    if level == "FAIL":
        N_FAIL += 1
    REPORT.append({"check": check, "expected": expected, "observed": str(observed),
                   "level": level, "note": note})
    mark = {"PASS": "✓", "WARN": "⚠", "FAIL": "✗"}.get(level, "?")
    print(f"[{mark} {level}] {check} | 期望={expected} 实测={observed}" + (f" | {note}" if note else ""))


# ------------------------------------------------------------------ 网络构件
def spearman_matrix(M: np.ndarray) -> np.ndarray:
    """metacell × gene 矩阵 → 基因间 Spearman(对列排秩后 Pearson)。恒定列 → NaN 行。"""
    R = rankdata(M, axis=0)
    R = R - R.mean(axis=0, keepdims=True)
    sd = np.sqrt((R ** 2).sum(axis=0))
    ok = sd > 0
    C = np.full((M.shape[1], M.shape[1]), np.nan)
    Rn = R[:, ok] / sd[ok]                    # 中心化后按 L2 归一 → 点积即 Pearson(秩) = Spearman
    C[np.ix_(ok, ok)] = Rn.T @ Rn
    return np.clip(C, -1.0, 1.0)


def _sft_fit(A: np.ndarray, n_bins: int = 10):
    """scale-free 拟合: 连通度 k 等宽直方图, log10(密度)~log10(k中心) 线性拟合 → (slope, R²)。

    WGCNA pickSoftThreshold 的 Python 近似。注: 不能用分位数分箱+计数——分位箱计数恒相等
    使 y 恒定、R² 无定义(2026-11-26 实测踩坑); 必须用等宽箱的**密度**(计数/箱宽)。
    """
    k = A.sum(axis=1) - 1.0
    k = np.clip(k, 1e-9, None)
    lo, hi = float(k.min()), float(k.max())
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.nan, np.nan
    edges = np.linspace(lo, hi, n_bins + 1)
    counts, _ = np.histogram(k, bins=edges)
    width = edges[1] - edges[0]
    centers = (edges[:-1] + edges[1:]) / 2.0
    m = counts > 0
    if m.sum() < 4:
        return np.nan, np.nan
    x = np.log10(centers[m])
    y = np.log10(counts[m] / (len(k) * width))
    slope, inter = np.polyfit(x, y, 1)
    pred = slope * x + inter
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return float(slope), float(r2)


def pick_beta(corr: np.ndarray, lo: int, hi: int, r2_target: float):
    rows = []
    for beta in range(int(lo), int(hi) + 1):
        A = ((1.0 + corr) / 2.0) ** beta
        np.fill_diagonal(A, 1.0)
        slope, r2 = _sft_fit(A)
        rows.append({"beta": beta, "sft_slope": slope, "sft_r2": r2,
                     "mean_connectivity": float((A.sum(axis=1) - 1).mean())})
    tab = pd.DataFrame(rows)
    if tab.sft_r2.isna().all():
        raise RuntimeError("SFT 拟合全 NA(连通度分布退化)——检查相关矩阵与 β 扫描实现")
    ok = tab[(tab.sft_r2 >= r2_target) & (tab.sft_slope < 0)]
    if len(ok):
        beta = int(ok.iloc[0].beta)
        note = "满足 R²≥%.2f 且负斜率的最小 β" % r2_target
    else:
        beta = int(tab.loc[tab.sft_r2.idxmax(), "beta"])
        note = f"无 β 满足 R²≥{r2_target} → 取 max R²(={tab.sft_r2.max():.2f}) 的 β; 网络非无标度, 谨慎解读"
    return beta, tab, note


def build_tom(corr: np.ndarray, beta: float) -> np.ndarray:
    A = ((1.0 + corr) / 2.0) ** float(beta)
    A = np.nan_to_num(A, nan=0.0)
    np.fill_diagonal(A, 1.0)
    k = A.sum(axis=1) - 1.0
    l = A @ A
    denom = np.minimum.outer(k, k) + 1.0 - A
    with np.errstate(divide="ignore", invalid="ignore"):
        tom = (l + A) / denom
    tom = np.clip(np.nan_to_num(tom, nan=0.0), 0.0, 1.0)
    np.fill_diagonal(tom, 1.0)
    return tom


def cut_modules(tom: np.ndarray, cut_height: float, min_size: int) -> np.ndarray:
    diss = 1.0 - tom
    np.fill_diagonal(diss, 0.0)
    diss = (diss + diss.T) / 2.0                       # 数值对称化
    Z = linkage(squareform(diss, checks=False), method="average")
    labels = fcluster(Z, t=cut_height, criterion="distance")
    out = np.zeros_like(labels)
    nxt = 1
    for cl in np.unique(labels):
        m = labels == cl
        if m.sum() >= min_size:
            out[m] = nxt
            nxt += 1
    return out                                         # 0 = grey(小簇并入)


def eigengene(M: np.ndarray, cols: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """模块内基因标准化 → PC1; 符号取向与成员均值谱正相关。返回 (eigengene, kME)。"""
    X = M[:, cols]
    Xz = (X - X.mean(axis=0)) / np.where(X.std(axis=0) > 0, X.std(axis=0), 1.0)
    U, S, Vt = np.linalg.svd(Xz, full_matrices=False)
    pc1 = U[:, 0] * S[0]
    if np.corrcoef(pc1, Xz.mean(axis=1))[0, 1] < 0:
        pc1 = -pc1
    kme = np.array([np.corrcoef(pc1, Xz[:, j])[0, 1] for j in range(Xz.shape[1])])
    return pc1, np.nan_to_num(kme)


def discover(M: np.ndarray, beta: float, cut_height: float, min_size: int) -> np.ndarray:
    corr = spearman_matrix(M)
    corr = np.nan_to_num(corr, nan=0.0)
    return cut_modules(build_tom(corr, beta), cut_height, min_size)


def jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if (a or b) else 0.0


def preservation(ref_sets: dict[int, set], labels: np.ndarray, gene_idx: np.ndarray) -> dict[int, float]:
    """每个参考模块 vs 重采样划分中最佳匹配簇的 Jaccard。"""
    cl_sets = {}
    for cl in np.unique(labels):
        if cl == 0:
            continue
        cl_sets[cl] = set(gene_idx[labels == cl])
    out = {}
    for mid, gs in ref_sets.items():
        out[mid] = max((jaccard(gs, cs) for cs in cl_sets.values()), default=0.0)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--config", default=None)
    ap.add_argument("--skip-bootstrap", action="store_true")
    args = ap.parse_args()
    root = os.path.abspath(args.root)
    cfg = atxconfig.load(args.config or os.path.join(root, "config", "config.yaml"))
    cfg.check_version_gate()
    t0 = time.time()
    cc = cfg.raw["network"]["coexpression"]
    outdir = cfg.results("network", mkdir=True)
    logs = cfg.results("logs", mkdir=True)

    Mdf = pd.read_csv(os.path.join(outdir, "metacells_log.tsv"), sep="\t", index_col=0)
    meta = pd.read_csv(os.path.join(outdir, "metacells_meta.tsv"), sep="\t")
    vartab = pd.read_csv(os.path.join(outdir, "metacells_var.tsv"), sep="\t")
    M = Mdf.to_numpy(dtype=np.float64)
    genes = np.asarray(Mdf.columns)
    print(f"== n01 modules | metacells={M.shape[0]} genes={M.shape[1]}")

    # ---- 基因过滤 ----
    detected = (M > 0).mean(axis=0)
    keep1 = detected >= float(cc["min_detected_frac"])
    var = M[:, keep1].var(axis=0)
    n_top = int(cc["n_top_var_genes"])
    top = np.argsort(var)[::-1][:n_top]
    gsel = np.flatnonzero(keep1)[top]
    gsel = np.sort(gsel)
    Ms = M[:, gsel]
    rec("基因过滤", f"检出率≥{cc['min_detected_frac']} 后取方差前 {n_top}",
        f"检出过滤 {int(keep1.sum())}/{len(genes)} → 入网 {len(gsel)}",
        "PASS" if len(gsel) >= n_top * 0.8 else "WARN")

    # ---- β 选择(全 metacell) ----
    corr_full = np.nan_to_num(spearman_matrix(Ms), nan=0.0)
    lo, hi = cc["beta_scan"]
    beta, sft_tab, beta_note = pick_beta(corr_full, lo, hi, float(cc["sft_r2"]))
    sft_tab.to_csv(os.path.join(outdir, "sft_beta.tsv"), sep="\t", index=False)
    rec("软阈值 β", f"scale-free R²≥{cc['sft_r2']} 最小 β", f"β={beta}",
        "PASS" if "满足" in beta_note else "WARN", note=beta_note)

    # ---- 全数据模块 ----
    labels_full = cut_modules(build_tom(corr_full, beta), float(cc["cut_height"]),
                              int(cc["min_module_size"]))
    mods = {int(c): np.flatnonzero(labels_full == int(c)) for c in np.unique(labels_full) if c != 0}
    ref_sets = {mid: set(gsel[idx]) for mid, idx in mods.items()}
    if not mods:
        # 空模块防护(2026-11-26 实测踩坑: cut_height 与 TOM 尺度不匹配时全 grey → 下游崩):
        # 记 FAIL、写空 summary(带列)、落日志后退出, 不做无意义的稳定性重建。
        rec("模块数(≥min_module_size)", f"cut_height={cc['cut_height']}, min_size={cc['min_module_size']}",
            f"0 个模块: 全部 {len(gsel)} 基因并入 grey", "FAIL",
            note="检查 cut_height 与 TOM 尺度是否匹配(诊断: 合并高度分位与切割扫描)")
        pd.DataFrame(columns=["module", "n_genes", "traj_rho_stage", "hub_genes", "mean_kme",
                              "lodo_jaccard", "bootstrap_jaccard", "frozen_candidate",
                              "members_file", "members_sha256"]).to_csv(
            os.path.join(outdir, "module_summary.tsv"), sep="\t", index=False)
        pd.DataFrame(REPORT).to_csv(os.path.join(logs, "n01_report.tsv"), sep="\t", index=False)
        return 1
    sym_of = dict(zip(vartab.ensg.astype(str), vartab.symbol.astype(str)))
    stage_order = meta.set_index("metacell_id")["stage_order"].reindex(Mdf.index).to_numpy()

    # 清理上一轮可能残留的 members 文件(模块数变化时防陈旧混杂; finding #10)
    import glob
    for stale in glob.glob(os.path.join(outdir, "module_M*_members.tsv")):
        os.unlink(stale)

    rows, eg_cols = [], {}
    for mid, idx in sorted(mods.items()):
        # idx 是"入选基因空间"(0..n_top-1)下标; genes 是全 27410 列索引; gsel 是两者桥梁。
        # 必须 genes[gsel[idx]] 才是真实 ENSG(2026-11-26 修复 finding #1: 原 genes[idx] 张冠李戴,
        # 致冻结成员/hub/权重哈希全部错位; modules.tsv 与 ref_sets 用对了 gsel, 故各门结果本身有效)。
        full_idx = gsel[idx]
        eg, kme = eigengene(Ms, idx)
        eg_cols[f"ME{mid}"] = eg
        rho_traj = spearmanr(eg, stage_order).statistic if len(np.unique(stage_order)) > 2 else np.nan
        hub_i = gsel[idx[np.argsort(kme)[::-1][:8]]]
        rows.append({"module": mid, "n_genes": len(idx),
                     "traj_rho_stage": rho_traj,
                     "hub_genes": ",".join(sym_of.get(g, str(g)) for g in genes[hub_i]),
                     "mean_kme": float(kme.mean())})
        mem = pd.DataFrame({"ensg": genes[full_idx],
                            "symbol": [sym_of.get(g, "") for g in genes[full_idx]],
                            "kme": kme, "eigengene_weight": kme})
        mem = mem.sort_values("kme", ascending=False)
        mem.to_csv(os.path.join(outdir, f"module_M{mid}_members.tsv"), sep="\t", index=False)
    summary = pd.DataFrame(rows)
    pd.DataFrame(eg_cols, index=Mdf.index).to_csv(
        os.path.join(outdir, "module_eigengenes.tsv"), sep="\t", float_format="%.6g")
    mod_of = np.zeros(len(gsel), dtype=int)
    for mid, idx in mods.items():
        mod_of[idx] = mid
    pd.DataFrame({"ensg": genes[gsel], "symbol": [sym_of.get(g, "") for g in genes[gsel]],
                  "module": mod_of, "module_name": [f"M{m}" if m else "grey" for m in mod_of]}
                 ).to_csv(os.path.join(outdir, "modules.tsv"), sep="\t", index=False)
    # 防回归守卫(finding #1 类): members.tsv 的 ENSG 集合必须 == modules.tsv 同模块集合
    mods_df = pd.read_csv(os.path.join(outdir, "modules.tsv"), sep="\t")
    for mid in mods:
        mem_set = set(pd.read_csv(os.path.join(outdir, f"module_M{mid}_members.tsv"),
                                  sep="\t").ensg.astype(str))
        truth = set(mods_df.loc[mods_df.module == mid, "ensg"].astype(str))
        if mem_set != truth:
            rec(f"members↔modules 一致性 M{mid}", "成员文件 ENSG 集合 == modules.tsv",
                f"不一致(交集 {len(mem_set & truth)}/{len(truth)})", "FAIL",
                note="gsel 桥接错误会在此现形(2026-11-26 finding #1)")
    rec("模块数(≥min_module_size)", f"cut_height={cc['cut_height']}, min_size={cc['min_module_size']}",
        f"{len(mods)} 个: " + ", ".join(f"M{m}={len(i)}" for m, i in sorted(mods.items()))
        + f"; grey={int((labels_full == 0).sum())}", "PASS" if mods else "FAIL")

    # ---- LODO 稳定性(冻结门) ----
    donors = meta.set_index("metacell_id")["donor"].reindex(Mdf.index)
    lodo_rows = []
    for d in sorted(donors.unique()):
        keep = (donors != d).to_numpy()
        if keep.sum() < max(len(Mdf) // 2, 5):
            lodo_rows.append({"donor": d, **{f"M{m}": np.nan for m in mods}, "note": "剩余样本过少, 跳过"})
            continue
        lab = discover(Ms[keep], beta, float(cc["cut_height"]), int(cc["min_module_size"]))
        pres = preservation(ref_sets, lab, gsel)
        lodo_rows.append({"donor": d, **{f"M{m}": pres[m] for m in mods}, "note": ""})
    lodo = pd.DataFrame(lodo_rows)
    lodo.to_csv(os.path.join(outdir, "lodo_stability.tsv"), sep="\t", index=False)
    lodo_mean = {m: float(lodo[f"M{m}"].mean(skipna=True)) for m in mods}
    rec("LODO 保留度(均值 Jaccard)", f"冻结门 ≥{cc['lodo_min_preservation']}",
        ", ".join(f"M{m}={v:.3f}" for m, v in sorted(lodo_mean.items())),
        "PASS" if all(v >= float(cc["lodo_min_preservation"]) for v in lodo_mean.values()) else "WARN",
        note="未达门者不冻结, 记录为候选失败")

    # ---- 供体 bootstrap(补充证据, 不做门) ----
    n_boot = 0 if args.skip_bootstrap else int(cc.get("n_donor_bootstrap", 0))
    if n_boot:
        rng = np.random.default_rng(cfg.seed + 101)
        donor_list = np.asarray(sorted(donors.unique()))
        blocks = {d: np.flatnonzero((donors == d).to_numpy()) for d in donor_list}
        boot_rows = []
        for b in range(n_boot):
            draw = rng.choice(donor_list, size=len(donor_list), replace=True)
            rows_idx = np.concatenate([blocks[d] for d in draw])
            lab = discover(Ms[rows_idx], beta, float(cc["cut_height"]), int(cc["min_module_size"]))
            pres = preservation(ref_sets, lab, gsel)
            boot_rows.append({"bootstrap": b, **{f"M{m}": pres[m] for m in mods}})
        boot = pd.DataFrame(boot_rows)
        boot.to_csv(os.path.join(outdir, "bootstrap_stability.tsv"), sep="\t", index=False)
        boot_mean = {m: float(boot[f"M{m}"].mean()) for m in mods}
        rec("供体 bootstrap 保留度", f"{n_boot} 次重建(补充)",
            ", ".join(f"M{m}={v:.3f}" for m, v in sorted(boot_mean.items())), "PASS")
        summary["bootstrap_jaccard"] = [boot_mean.get(m, np.nan) for m in summary.module]
    summary["lodo_jaccard"] = [lodo_mean.get(m, np.nan) for m in summary.module]
    summary["frozen_candidate"] = [lodo_mean.get(m, 0) >= float(cc["lodo_min_preservation"])
                                   for m in summary.module]
    # 成员文件 sha256(冻结表 weights_file_hash 用)
    summary["members_file"] = [f"module_M{m}_members.tsv" for m in summary.module]
    summary["members_sha256"] = [
        hashlib.file_digest(open(os.path.join(outdir, f"module_M{m}_members.tsv"), "rb"),
                            "sha256").hexdigest() for m in summary.module]
    summary.to_csv(os.path.join(outdir, "module_summary.tsv"), sep="\t", index=False)
    rec("冻结候选(LODO≥门)", f"≤ module_cap={cfg.raw['network']['module_cap']} 个进 n02 复核",
        f"{int(summary.frozen_candidate.sum())} 个候选: "
        + ",".join(f"M{m}" for m in summary.loc[summary.frozen_candidate, "module"]),
        "PASS", note="最终冻结还需 n02 Aldinger/Castonguay 两门")

    pd.DataFrame(REPORT).to_csv(os.path.join(logs, "n01_report.tsv"), sep="\t", index=False)
    print(f"\n汇总: PASS={sum(1 for r in REPORT if r['level'] == 'PASS')} "
          f"WARN={sum(1 for r in REPORT if r['level'] == 'WARN')} FAIL={N_FAIL}  ({time.time() - t0:.0f}s)")
    print(summary.to_string(index=False))
    return 1 if N_FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
