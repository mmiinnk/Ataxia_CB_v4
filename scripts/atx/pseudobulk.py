"""atx.pseudobulk — 细胞×基因 → (donor×bin×group) → (group×bin) 两级聚合(spec §3.5)。

借鉴(Depression_CB 已验证口径):
① 归一化 = X / obs.size_factor × target(Sepp size_factor ≡ 作者基因过滤前总 UMI, s00 复验;
   外部数据集无此列时调用方传 size_factor = raw/X 行和);
② 稀疏指示矩阵 S @ X 一次性聚合, 层 mean_log / mean_cp10k / pct_pos;
③ 全骨架 reindex:不可行 bin 显式存在(行在、值 NA), 不静默丢行;
④ 单供体/pooled bin 的 sem = NA(硬断言;硬约束 5:绝不插值、绝不填 0)。

本项目差异:骨架 = plot_groups(10 组) × bins;feasible 由 min_cells 参数化 →
min_cells_grid 敏感性**只重算 bin 表, 不重算 donor cube**;bin 级值 = 供体等权均值
(bin_equal_weight, 与 config stats.specificity.aggregation 一致), 细胞数只进 feasible/n_cells_total。

pandas 3 注意:groupby 列表键返回 tuple、observed=True 丢空组合 → 本模块用显式全骨架
product + numpy 索引, 不依赖 groupby 的组合补全;NaN 不进 bool 列(feasible/single_donor 是
纯 bool, NA 只出现在 float 值列)。
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
import scipy.sparse as sp

META_COLS = ["group", "bin", "bin_order", "period", "n_donors", "n_cells_total",
             "feasible", "single_donor", "pooled"]


def normalize(X, size_factor, target: float = 1e4):
    """X(原始计数 CSR) → (X_log, X_cp):X_cp = X/sf×target;X_log = log1p(X_cp)。float32。"""
    sf = np.asarray(size_factor, dtype=np.float64).ravel()
    if sf.shape[0] != X.shape[0]:
        raise ValueError(f"size_factor 长度 {sf.shape[0]} != X 行数 {X.shape[0]}")
    if not np.all(np.isfinite(sf)) or np.any(sf <= 0):
        bad = int(np.sum(~(np.isfinite(sf) & (sf > 0))))
        raise ValueError(f"size_factor 有 {bad} 个非正/非有限值——不得填 0 或跳过(硬约束 5)")
    if not sp.issparse(X):
        raise ValueError("normalize 只接受稀疏 X(全稠密会爆内存;spec §10)")
    scale = sp.diags((target / sf).astype(np.float32), format="csr")
    X_cp = (scale @ X.astype(np.float32)).tocsr()
    X_log = X_cp.copy()
    X_log.data = np.log1p(X_log.data, dtype=np.float32)
    return X_log, X_cp


def donor_cube(X_log, X_cp, obs: pd.DataFrame, obs_long: pd.DataFrame, var: pd.DataFrame):
    """(donor × bin × group) pseudobulk → AnnData(spec §3.5 签名)。

    X = mean log1p CP10K;layers: mean_cp10k / pct_pos;obs: group/bin/bin_order/donor/
    n_cells/pooled/disease_flag(有则)/period 由调用方经 spcfg 补(此处不依赖 cfg)。
    输入 X_* 必须已按 var(=panel 基因)切好列;obs 行序 = X 行序 = 细胞位置。
    """
    import anndata

    cell_idx = obs_long["cell_pos"].to_numpy(dtype=np.int64)
    grp = obs_long["cell_group"].to_numpy()
    binv = obs["bin"].to_numpy()[cell_idx]
    don = obs["donor"].to_numpy()[cell_idx]
    bord = obs["bin_order"].to_numpy()[cell_idx]
    pooled = obs["pooled_bin"].to_numpy()[cell_idx] if "pooled_bin" in obs else np.zeros(len(cell_idx), bool)
    disease = (obs["disease_flag"].to_numpy()[cell_idx]
               if "disease_flag" in obs.columns else np.zeros(len(cell_idx), bool))

    key = pd.Series(grp).astype(str) + "|" + pd.Series(binv).astype(str) + "|" + pd.Series(don).astype(str)
    codes, uniq = pd.factorize(key, sort=False)
    n_combo = len(uniq)
    # combo 元信息(按 group, bin_order, donor 排序, 输出可读且确定)
    parts = pd.Series(uniq).str.split("|", n=2, expand=True)
    meta = pd.DataFrame({"group": parts[0], "bin": parts[1], "donor": parts[2]})
    meta["bin_order"] = meta["bin"].map(dict(zip(binv, bord))).astype(int)
    meta = meta.sort_values(["group", "bin_order", "donor"], kind="stable").reset_index(drop=True)
    # 旧 combo 位置 → 排序后新位置(显式字典, 键唯一由 factorize 保证)
    new_key = meta["group"] + "|" + meta["bin"] + "|" + meta["donor"]
    remap = pd.Series(np.arange(n_combo), index=new_key).reindex(uniq).to_numpy()
    if np.any(pd.isna(remap)):
        raise AssertionError("donor_cube: combo 重映射出现未命中键(内部错误)")
    codes_sorted = remap[codes].astype(np.int64)

    n_cells = np.bincount(codes_sorted, minlength=n_combo).astype(np.int64)
    S = sp.csr_matrix((np.ones(len(codes_sorted), dtype=np.float32), (codes_sorted, cell_idx)),
                      shape=(n_combo, X_log.shape[0]))
    inv_n = sp.diags((1.0 / np.maximum(n_cells, 1)).astype(np.float32), format="csr")
    mean_log = (inv_n @ (S @ X_log)).tocsr()
    mean_cp = (inv_n @ (S @ X_cp)).tocsr()
    B = X_cp.copy()
    B.data = (B.data > 0).astype(np.float32)
    pct_pos = (inv_n @ (S @ B)).tocsr()

    # donor 级 pooled/disease 一致性(同一 donor×bin×group 内必须恒定)
    df_tmp = pd.DataFrame({"c": codes_sorted, "p": pooled, "d": disease})
    agg = df_tmp.groupby("c").agg(p_nun=("p", "nunique"), d_nun=("d", "nunique"))
    if (agg["p_nun"] > 1).any() or (agg["d_nun"] > 1).any():
        raise AssertionError("donor cube: pooled/disease_flag 在同一 (group,bin,donor) 内不恒定")
    meta["n_cells"] = n_cells
    meta["pooled"] = df_tmp.groupby("c")["p"].first().reindex(np.arange(n_combo)).to_numpy()
    meta["disease_flag"] = df_tmp.groupby("c")["d"].first().reindex(np.arange(n_combo)).to_numpy()
    meta = meta.set_index(meta["group"] + "|" + meta["bin"] + "|" + meta["donor"])
    meta.index.name = "combo"

    cube = anndata.AnnData(X=mean_log, obs=meta[["group", "bin", "bin_order", "donor",
                                                 "n_cells", "pooled", "disease_flag"]],
                           var=var.copy(),
                           layers={"mean_cp10k": mean_cp, "pct_pos": pct_pos})
    return cube


def bin_table(cube, spcfg, cfg, min_cells: int | None = None) -> pd.DataFrame:
    """(group × bin) 全骨架表(spec §3.5/§5)。行 = plot_groups × spcfg.labels(缺行补 NA 行)。

    值 = 供体等权均值(log_/cp_/pct_<gene>);logsem_<gene> = 供体间 sem(ddof=1),
    单供体或 pooled bin → NA(硬断言);n_cells_total < min_cells → feasible=False 且
    全部基因值列 = NA(硬约束 5:表 NA、图斜纹;绝不插值/填 0)。
    min_cells 参数化 → min_cells_grid 敏感性只重算本表(donor cube 复用)。
    """
    if min_cells is None:
        min_cells = int(cfg.bins_cfg["min_cells"])
    genes = list(cube.var["gene_key"])
    o = cube.obs
    Xd = np.asarray(cube.X.todense(), dtype=np.float64)
    CPd = np.asarray(cube.layers["mean_cp10k"].todense(), dtype=np.float64)
    PCd = np.asarray(cube.layers["pct_pos"].todense(), dtype=np.float64)
    g_arr = o["group"].to_numpy()
    b_arr = o["bin"].to_numpy()
    pooled_arr = o["pooled"].to_numpy(dtype=bool)
    ncells_arr = o["n_cells"].to_numpy(dtype=np.int64)
    idx_by = {}
    for i, (g, b) in enumerate(zip(g_arr, b_arr)):
        idx_by.setdefault((g, b), []).append(i)

    groups = cfg.plot_groups()
    labels = spcfg.labels
    recs = []
    for g, b in itertools.product(groups, labels):
        idx = np.asarray(idx_by.get((g, b), []), dtype=np.int64)
        n_donors = len(idx)
        n_cells_total = int(ncells_arr[idx].sum()) if n_donors else 0
        pooled = bool(pooled_arr[idx].any()) if n_donors else False
        feasible = n_cells_total >= min_cells
        single = n_donors <= 1
        rec = {"group": g, "bin": b, "bin_order": spcfg.bin_order[b],
               "period": spcfg.period_of[b], "n_donors": n_donors,
               "n_cells_total": n_cells_total, "feasible": feasible,
               "single_donor": single, "pooled": pooled}
        if n_donors and feasible:
            logm = Xd[idx].mean(axis=0)
            cpm = CPd[idx].mean(axis=0)
            pcm = PCd[idx].mean(axis=0)
            if n_donors >= 2 and not pooled:
                sem = Xd[idx].std(axis=0, ddof=1) / np.sqrt(n_donors)
            else:
                sem = np.full(len(genes), np.nan)
        else:
            logm = cpm = pcm = sem = np.full(len(genes), np.nan)
        for j, gene in enumerate(genes):
            rec[f"log_{gene}"] = logm[j]
        for j, gene in enumerate(genes):
            rec[f"cp_{gene}"] = cpm[j]
        for j, gene in enumerate(genes):
            rec[f"pct_{gene}"] = pcm[j]
        for j, gene in enumerate(genes):
            rec[f"logsem_{gene}"] = sem[j]
        recs.append(rec)
    bt = pd.DataFrame(recs)

    # ---- 硬断言(spec §3.5 ④/硬约束 5) ----
    log_cols = [f"log_{g}" for g in genes]
    sem_cols = [f"logsem_{g}" for g in genes]
    bad_sem = bt.loc[bt["single_donor"] | bt["pooled"], sem_cols].notna().any(axis=None)
    if bad_sem:
        raise AssertionError("bin_table: 单供体/pooled bin 出现非 NA 的 logsem(硬约束 5)")
    bad_feas = bt.loc[~bt["feasible"], log_cols].notna().any(axis=None)
    if bad_feas:
        raise AssertionError("bin_table: 不可行 bin 出现非 NA 基因值(硬约束 5)")
    if len(bt) != len(groups) * len(labels):
        raise AssertionError("bin_table: 骨架不完整")
    return bt


def feasibility_frame(bt: pd.DataFrame) -> pd.DataFrame:
    """bin_table → group×bin 的 n_cells_total 矩阵(与 §4.2 锚点比对用)。"""
    return bt.pivot(index="group", columns="bin", values="n_cells_total")


def spot_check(X_log, obs: pd.DataFrame, obs_long: pd.DataFrame, bt: pd.DataFrame,
               genes: list[str], n: int = 3, seed: int = 0) -> pd.DataFrame:
    """随机抽 n 个可行 (gene, group, bin), 用 numpy 从细胞级直接重算两级均值比对(spec §3.5)。

    返回比对表(check=PASS/FAIL);任何 FAIL 应由调用方升级为断言失败。
    """
    rng = np.random.default_rng(seed)
    feas = bt[bt["feasible"] & (bt["n_donors"] >= 1)]
    picks = feas.sample(n=min(n, len(feas)), random_state=int(rng.integers(1 << 30)))
    out = []
    for _, r in picks.iterrows():
        g, b = r["group"], r["bin"]
        cand = [x for x in genes if pd.notna(r.get(f"log_{x}"))]
        if not cand:
            continue
        gene = cand[rng.integers(len(cand))]
        j = genes.index(gene)
        cells = obs_long.loc[obs_long["cell_group"] == g, "cell_pos"].to_numpy()
        cells = cells[obs["bin"].to_numpy()[cells] == b]
        col = np.asarray(X_log[cells, j].todense(), dtype=np.float64).ravel()
        don = obs["donor"].to_numpy()[cells]
        donor_means = np.array([col[don == d].mean() for d in pd.unique(don)])
        v = donor_means.mean()
        exp = float(r[f"log_{gene}"])
        ok_n = int(r["n_cells_total"]) == len(cells)
        ok_v = abs(v - exp) <= 1e-5 * max(1.0, abs(exp))
        out.append({"gene": gene, "group": g, "bin": b, "n_cells": len(cells),
                    "n_cells_table": int(r["n_cells_total"]), "recomputed": v, "table": exp,
                    "n_donors_recomputed": len(donor_means), "n_donors_table": int(r["n_donors"]),
                    "check": "PASS" if (ok_n and ok_v and len(donor_means) == int(r["n_donors"])) else "FAIL"})
    return pd.DataFrame(out)
