"""tests/synth.py — 合成数据构建(spec §1 第一步:200 核 × 20 基因 × 3 组 × 5 bin)。

细胞布局完全由 PLAN 决定(确定性, 无随机):act ∈ {Purkinje, GC, immune, meningeal}
(meningeal = ungrouped, 只进 ALL);bin ∈ 鼠 E10.5–E14.5(E10.5/E11.5 = pooled, 来自真实
config);供体设计覆盖:多供体 bin(E12.5/E13.5)、单供体 bin(E14.5, sem 必须 NA)、
pooled bin(sem 必须 NA)、低核数组×bin(immune×E12.5 = 5 核 → 不可行)。
计数为确定性整数伪随机(含结构性 0, 供 pct_pos 校验);size_factor = 行和 × 1.01
(模拟 Sepp 作者过滤前总 UMI 与 X 行和 ~1% 偏差)。
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

GENES = [f"G{i}" for i in range(20)]

# (author_cell_type, bin_raw, donor, n_cells) —— 合计 200
PLAN = [
    ("Purkinje", "E10.5", "pooled", 10),
    ("Purkinje", "E11.5", "pooled", 12),
    ("Purkinje", "E12.5", "d1", 20),
    ("Purkinje", "E12.5", "d2", 16),
    ("Purkinje", "E13.5", "d1", 14),
    ("Purkinje", "E13.5", "d2", 12),
    ("Purkinje", "E13.5", "d3", 10),
    ("Purkinje", "E14.5", "d1", 9),
    ("GC", "E12.5", "d1", 15),
    ("GC", "E12.5", "d2", 15),
    ("GC", "E13.5", "d1", 12),
    ("GC", "E13.5", "d2", 12),
    ("GC", "E13.5", "d3", 12),
    ("GC", "E14.5", "d1", 10),
    ("immune", "E12.5", "d2", 5),
    ("meningeal", "E12.5", "d1", 16),
]
N_CELLS = sum(p[3] for p in PLAN)
ACT_TO_GROUP = {"Purkinje": "Purkinje", "GC": "GranuleLineage", "immune": "Immune",
                "meningeal": "ungrouped"}


def counts_matrix() -> np.ndarray:
    """确定性整数计数(200×20, 含结构性 0)。"""
    n, g = N_CELLS, len(GENES)
    ii = np.repeat(np.arange(n), g).reshape(n, g)
    jj = np.tile(np.arange(g), n).reshape(n, g)
    vals = np.where((ii * 3 + jj) % 5 == 0, 0, ((ii * 7 + jj * 13) % 9) + 1)
    return vals.astype(np.float32)


def build():
    """返回 (obs_raw, X_counts_csr, var)——obs_raw 模拟 Sepp obs 的必要列。"""
    assert N_CELLS == 200, f"合成规模须为 200 核(spec §1), 现 {N_CELLS}"
    act = np.repeat([p[0] for p in PLAN], [p[3] for p in PLAN])
    binv = np.repeat([p[1] for p in PLAN], [p[3] for p in PLAN])
    don = np.repeat([p[2] for p in PLAN], [p[3] for p in PLAN])
    n = len(act)
    obs = pd.DataFrame({
        "author_cell_type": act,
        "subtype": np.full(n, np.nan, dtype=object),
        "dev_state": np.full(n, np.nan, dtype=object),
        "cell_type": np.full(n, np.nan, dtype=object),
        "author_stage": binv,
        "donor_id": don,
        "disease": "normal",
    })
    X = sp.csr_matrix(counts_matrix())
    rowsum = np.asarray(X.sum(axis=1)).ravel()
    obs["size_factor"] = rowsum * 1.01
    var = pd.DataFrame({"feature_name": GENES},
                       index=[f"ENSMUSG_FAKE_{j:04d}" for j in range(len(GENES))])
    var.index.name = "gene_ids"
    return obs, X, var


def cell_slices() -> dict[tuple[str, str, str], slice]:
    """PLAN 顺序下每个 (act, bin, donor) 的细胞位置切片(手算期望值用)。"""
    out, start = {}, 0
    for a, b, d, n in PLAN:
        out[(a, b, d)] = slice(start, start + n)
        start += n
    return out


# ---------------------------------------------------------------------------
# stats_* 测试用:bin_table / donor cube 直建(绕开细胞级链路, 期望值可手算)
# ---------------------------------------------------------------------------

def mk_bin_table(cfg, rows, genes, value_kind="log", species="mouse") -> pd.DataFrame:
    """合成 bin_table(列结构与 pseudobulk.bin_table 一致)。

    rows: (group, bin_label, feasible, n_cells_total, {gene: value});
    value_kind='log' → value 为 log_<gene>, cp_=expm1(value);'cp' → value 为 cp_<gene>,
    log_=log1p(value)。pct_=0.5(表达核比例 0–1, 与真实 pct_pos 层同口径; 可行)/NA;
    logsem_=NA;n_donors=2(可行)/0;
    single_donor/pooled=False(趋势与特异性统计不读这两列, 结构性列补齐即可)。
    """
    spcfg = cfg.species(species)
    recs = []
    for group, b, feas, ncells, vals in rows:
        r = {"group": group, "bin": b, "bin_order": spcfg.bin_order[b],
             "period": spcfg.period_of[b], "n_donors": 2 if feas else 0,
             "n_cells_total": int(ncells), "feasible": bool(feas),
             "single_donor": False, "pooled": False}
        for g in genes:
            v = float(vals.get(g, np.nan)) if feas else np.nan
            fin = np.isfinite(v)
            lg = v if (value_kind == "log" and fin) else (np.log1p(v) if fin else np.nan)
            cp = np.expm1(v) if (value_kind == "log" and fin) else (v if fin else np.nan)
            r[f"log_{g}"], r[f"cp_{g}"] = lg, cp
            r[f"pct_{g}"] = 0.5 if feas else np.nan
            r[f"logsem_{g}"] = np.nan
        recs.append(r)
    return pd.DataFrame(recs)


def mk_cube(cfg, rows, genes, species="mouse"):
    """合成 donor cube AnnData(行 = donor×bin×group 组合;X = mean log1p CP10K)。

    rows: (group, bin_label, donor, {gene: log 值})。obs 列与 pseudobulk.donor_cube 产物
    的关键列一致(group/bin/bin_order/donor/period/n_cells/pooled);var 含 gene_key。
    """
    import anndata
    spcfg = cfg.species(species)
    obs_recs, X = [], []
    for group, b, donor, vals in rows:
        obs_recs.append({"group": group, "bin": b, "bin_order": spcfg.bin_order[b],
                         "donor": donor, "period": spcfg.period_of[b],
                         "n_cells": 20, "pooled": False})
        X.append([float(vals.get(g, np.nan)) for g in genes])
    obs = pd.DataFrame(obs_recs)
    var = pd.DataFrame({"gene_key": list(genes)}, index=list(genes))
    return anndata.AnnData(X=np.asarray(X, dtype=np.float64), obs=obs, var=var)
