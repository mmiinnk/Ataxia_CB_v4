"""pseudobulk.py 测试(spec §3.13):bin 均值 = 手算;单供体/pooled sem NA;不可行 bin NA;
全骨架;min_cells 参数化(同一 donor cube 重算 bin 表);spot_check。全部期望值用独立的
numpy 稠密代码手算, 不复用被测模块的实现路径。"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

import _common
import synth
from atx import groups as atxgroups
from atx import pseudobulk as atxpb


def _pipeline(min_cells=None):
    c = _common.cfg()
    spcfg = c.species("mouse")
    obs_raw, X, var = synth.build()
    d = atxgroups.derive_obs(obs_raw, spcfg, c)
    a, _ = atxgroups.assign(d, c)
    lg = atxgroups.expand_long(a, c)
    sf = obs_raw["size_factor"].to_numpy()
    X_log, X_cp = atxpb.normalize(X, sf, target=1e4)
    var_panel = var.copy()
    var_panel["gene_key"] = synth.GENES
    cube = atxpb.donor_cube(X_log, X_cp, d, lg, var_panel)
    bt = atxpb.bin_table(cube, spcfg, c, min_cells=min_cells)
    # 独立手算基准(稠密)
    vals = X.toarray().astype(np.float64)
    cp = vals / sf[:, None] * 1e4
    logv = np.log1p(cp)
    return dict(c=c, spcfg=spcfg, d=d, a=a, lg=lg, X_log=X_log, X_cp=X_cp,
                cube=cube, bt=bt, vals=vals, cp=cp, logv=logv)


def _bt_row(bt, g, b):
    r = bt[(bt["group"] == g) & (bt["bin"] == b)]
    assert len(r) == 1, f"骨架行缺失/重复: {g}×{b}"
    return r.iloc[0]


def test_normalize_math():
    p = _pipeline()
    obs_raw, X, _ = synth.build()
    sfv = obs_raw["size_factor"].to_numpy()
    exp_cp = X.toarray().astype(np.float32) / sfv[:, None].astype(np.float32) * np.float32(1e4)
    got_cp = p["X_cp"].toarray()
    assert np.allclose(got_cp, exp_cp, rtol=1e-5), "X_cp != X/size_factor×1e4"
    assert np.allclose(p["X_log"].toarray(), np.log1p(exp_cp), rtol=1e-5)
    # 稀疏结构:log1p 只作用于 .data, 结构性 0 保持 0
    assert np.allclose(p["X_log"].data, np.log1p(p["X_cp"].data))


def test_normalize_guards():
    obs_raw, X, _ = synth.build()
    sf = obs_raw["size_factor"].to_numpy().copy()
    try:
        sf[0] = 0.0
        atxpb.normalize(X, sf, 1e4)
        raise AssertionError("size_factor=0 未报错(硬约束 5:不得填 0/跳过)")
    except ValueError:
        pass
    try:
        atxpb.normalize(X.toarray(), obs_raw["size_factor"].to_numpy(), 1e4)
        raise AssertionError("稠密 X 未拒绝")
    except ValueError:
        pass


def test_donor_cube_hand():
    p = _pipeline(min_cells=8)
    cube = p["cube"]
    sl = synth.cell_slices()
    j = 3  # 基因 G3
    # combo (Purkinje, E13.5, d1): n_cells=14, X=donor 内细胞 log 均值(手算)
    key = "Purkinje|E13.5|d1"
    assert key in cube.obs_names
    row = cube.obs.loc[key]
    assert int(row["n_cells"]) == 14 and int(row["bin_order"]) == 4
    vkey = cube.var.index[j]      # cube var 索引 = 原 var 索引(ENSMUSG_FAKE_000j), 非 gene_key
    cells = sl[("Purkinje", "E13.5", "d1")]
    exp_mean_log = p["logv"][cells, j].mean()
    got = cube[key, vkey].X.toarray().ravel()[0]
    assert abs(got - exp_mean_log) < 1e-5, f"donor mean_log: {got} != {exp_mean_log}"
    exp_mean_cp = p["cp"][cells, j].mean()
    got_cp = cube[key, vkey].layers["mean_cp10k"].toarray().ravel()[0]
    assert abs(got_cp - exp_mean_cp) < 1e-3
    exp_pct = float((p["vals"][cells, j] > 0).mean())
    got_pct = cube[key, vkey].layers["pct_pos"].toarray().ravel()[0]
    assert abs(got_pct - exp_pct) < 1e-5
    # pooled 旗标
    assert bool(cube.obs.loc["Purkinje|E10.5|pooled", "pooled"]) is True
    assert bool(cube.obs.loc[key, "pooled"]) is False
    assert not cube.obs["disease_flag"].any()


def test_bin_table_hand():
    p = _pipeline(min_cells=8)
    bt, cube = p["bt"], p["cube"]
    sl = synth.cell_slices()
    # --- 多供体 bin:供体等权均值 + sem(ddof=1)/sqrt(n) 手算 ---
    donor_means = np.array([p["logv"][sl[("Purkinje", "E13.5", dn)], 0].mean()
                            for dn in ("d1", "d2", "d3")])
    r = _bt_row(bt, "Purkinje", "E13.5")
    assert int(r["n_donors"]) == 3 and int(r["n_cells_total"]) == 36
    assert r["feasible"] and not r["single_donor"] and not r["pooled"]
    assert abs(r["log_G0"] - donor_means.mean()) < 1e-5
    exp_sem = donor_means.std(ddof=1) / np.sqrt(3)
    assert abs(r["logsem_G0"] - exp_sem) < 1e-5, "sem != 供体间 ddof=1 sem"
    # --- 单供体 bin:值在、sem 必 NA ---
    r14 = _bt_row(bt, "Purkinje", "E14.5")
    assert r14["single_donor"] and r14["feasible"] and int(r14["n_donors"]) == 1
    assert np.isnan(r14["logsem_G0"]) and not np.isnan(r14["log_G0"])
    # --- pooled bin:sem 必 NA ---
    for b in ("E10.5", "E11.5"):
        rp = _bt_row(bt, "Purkinje", b)
        assert rp["pooled"] and np.isnan(rp["logsem_G0"])
    # --- 不可行 bin(Immune×E12.5 = 5 核 < 8):计数在、基因值全 NA ---
    ri = _bt_row(bt, "Immune", "E12.5")
    assert not ri["feasible"] and int(ri["n_cells_total"]) == 5
    assert all(np.isnan(ri[f"log_{g}"]) for g in synth.GENES)
    assert all(np.isnan(ri[f"cp_{g}"]) for g in synth.GENES)
    # --- 无细胞组合:行存在、n=0、NA ---
    r0 = _bt_row(bt, "GranuleLineage", "E10.5")
    assert int(r0["n_donors"]) == 0 and int(r0["n_cells_total"]) == 0 and np.isnan(r0["log_G0"])
    # --- ALL = 全体细胞(ungrouped 只进 ALL) ---
    assert int(_bt_row(bt, "ALL", "E12.5")["n_cells_total"]) == 36 + 30 + 5 + 16
    assert int(_bt_row(bt, "ALL", "E13.5")["n_cells_total"]) == 36 + 36
    # --- 全骨架:10 组 × 12 bin ---
    assert len(bt) == 10 * 12
    assert list(bt["group"].unique()) == p["c"].plot_groups()
    # --- pct 列 = 供体等权 pct 均值 ---
    dm_pct = np.array([float((p["vals"][sl[("Purkinje", "E13.5", dn)], 0] > 0).mean())
                       for dn in ("d1", "d2", "d3")])
    assert abs(_bt_row(bt, "Purkinje", "E13.5")["pct_G0"] - dm_pct.mean()) < 1e-5


def test_min_cells_parameterized_same_cube():
    """min_cells_grid 敏感性:同一 donor cube, 只重算 bin 表(spec §3.5 差异 ②)。"""
    p = _pipeline(min_cells=30)          # config 默认
    bt30 = p["bt"]
    assert _bt_row(bt30, "Purkinje", "E12.5")["feasible"]        # 36 ≥ 30
    assert not _bt_row(bt30, "Purkinje", "E14.5")["feasible"]    # 9 < 30
    assert _bt_row(bt30, "GranuleLineage", "E13.5")["feasible"]  # 36 ≥ 30
    assert not _bt_row(bt30, "Immune", "E12.5")["feasible"]      # 5 < 30
    sum_before = float(p["cube"].X.sum())
    bt8 = atxpb.bin_table(p["cube"], p["spcfg"], p["c"], min_cells=8)
    assert _bt_row(bt8, "Purkinje", "E14.5")["feasible"]         # 9 ≥ 8
    assert not _bt_row(bt8, "Immune", "E12.5")["feasible"]       # 5 < 8
    bt5 = atxpb.bin_table(p["cube"], p["spcfg"], p["c"], min_cells=5)
    assert _bt_row(bt5, "Immune", "E12.5")["feasible"]           # 5 ≥ 5
    # 重算 bin 表不触 donor cube(内容不变)
    assert float(p["cube"].X.sum()) == sum_before


def test_spot_check_passes():
    p = _pipeline(min_cells=8)
    sc = atxpb.spot_check(p["X_log"], p["d"], p["lg"], p["bt"], synth.GENES, n=3, seed=42)
    assert len(sc) >= 1
    assert (sc["check"] == "PASS").all(), f"spot_check FAIL:\n{sc}"


def test_feasibility_frame():
    p = _pipeline(min_cells=8)
    ff = atxpb.feasibility_frame(p["bt"])
    assert ff.loc["Purkinje", "E13.5"] == 36
    assert ff.loc["Immune", "E12.5"] == 5
    assert ff.loc["ALL", "E12.5"] == 87
