"""stats_trend 测试(spec §3.13):精确置换 p 手算/暴力穷举核对(含 ties)、scipy 路、
BH q、trend_table 门与 low_power、donor 级辅证、bin_order 一致性守卫、确定性(R8)。"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

import _common
import synth
from atx import stats_trend as st


def test_spearman_exact_monotone_n5():
    rho, p = st.spearman_exact([1, 2, 3, 4, 5], [0.1, 0.5, 0.9, 1.2, 2.0], exact_max_n=6)
    assert rho == 1.0
    assert abs(p - 2 / 120) < 1e-12      # 双侧: 5!=120 个排列中仅恒等/全逆 |rho|=1
    assert p > 0                          # scipy 在 rho=±1 时给 p=0——精确置换必须 >0


def test_spearman_exact_antitone_n6_cache_deterministic():
    t = list(range(1, 7))
    x = [3.0, 2.5, 2.0, 1.5, 1.0, 0.5]
    r1 = st.spearman_exact(t, x, exact_max_n=6)
    r2 = st.spearman_exact(t, x, exact_max_n=6)
    assert r1 == r2                        # 确定性(无 RNG)
    assert r1[0] == -1.0 and abs(r1[1] - 2 / 720) < 1e-12
    assert len(st._PERM_CACHE[6]) == 720   # 排列表按 n 缓存一次(§3.6 差异③)


def test_spearman_exact_matches_bruteforce_with_ties():
    t = np.array([1.0, 2.0, 3.0, 4.0])
    x = np.array([1.0, 4.0, 2.0, 2.0])     # ties → 平均秩
    rho, p = st.spearman_exact(t, x, exact_max_n=6)
    rt = rankdata(t)
    rho_ref = float(np.corrcoef(rt, rankdata(x))[0, 1])
    assert abs(rho - rho_ref) < 1e-12
    cnt = tot = 0
    for perm in itertools.permutations(x):  # 测试内独立穷举(不复用模块代码)
        rp = np.corrcoef(rt, rankdata(np.array(perm)))[0, 1]
        tot += 1
        if abs(rp) >= abs(rho_ref) - 1e-12:
            cnt += 1
    assert tot == 24 and abs(p - cnt / tot) < 1e-12


def test_spearman_above_exact_max_uses_scipy():
    cfg = _common.cfg()
    rng = np.random.default_rng(cfg.seed)
    t = np.arange(1, 8, dtype=float)
    x = rng.normal(size=7)
    rho, p = st.spearman_exact(t, x, exact_max_n=6)
    ref = spearmanr(t, x)
    assert abs(rho - ref.statistic) < 1e-12 and abs(p - ref.pvalue) < 1e-12


def test_spearman_degenerate_inputs():
    assert all(np.isnan(v) for v in st.spearman_exact([1.0], [2.0]))
    assert all(np.isnan(v) for v in st.spearman_exact([1, 2, 3], [5.0, 5.0, 5.0]))
    rho, p = st.spearman_exact([1, np.nan, 3, 4], [1.0, 9.0, 3.0, 4.0])
    assert rho == 1.0 and abs(p - 2 / 6) < 1e-12   # NaN 对剔除 → n=3, 2/3!


def test_bh_qvalues():
    q = st.bh_qvalues([0.01, 0.02, 0.03, 0.4, np.nan])
    assert np.isnan(q[4])
    # n=4: q_i = p_i*n/i 再从大向小单调化 → [0.04, 0.04, 0.04, 0.4]
    assert np.allclose(q[:4], [0.04, 0.04, 0.04, 0.4])
    q2 = st.bh_qvalues([0.4, 0.01, 0.03, 0.02])    # 乱序输入按原位还原
    assert np.allclose(q2, [0.4, 0.04, 0.04, 0.04])
    assert np.all(np.isnan(st.bh_qvalues([np.nan, np.nan])))


def _build_case(cfg):
    """Immune(4 bin, override→ok) / Purkinje(4 bin→insufficient_bins) /
    DN(5 bin 但总细胞<门→insufficient_cells) / OligoLineage(5 bin→ok) + 缺失基因 GB。"""
    spcfg = cfg.species("mouse")
    L = spcfg.labels
    genes = ["GA", "GB"]
    min_total = int(cfg.bins_cfg["min_cells_total_for_trend"])
    rows, cube_rows = [], []
    for i, b in enumerate(L[1:5]):                       # Immune: 4 可行 bin 单调升
        rows.append(("Immune", b, True, 40, {"GA": 0.5 + i}))
        for j in (1, 2):
            cube_rows.append(("Immune", b, f"d{j}", {"GA": 0.5 + i}))   # bin 内同值 → rho_donor=1
    for i, b in enumerate(L[1:5]):                       # Purkinje: 同 4 bin, 无 override
        rows.append(("Purkinje", b, True, 40, {"GA": 0.5 + i}))
    nc = min_total // 5 - 1                              # DN: 5 bin 总细胞 < min_total
    for i, b in enumerate(L[:5]):
        rows.append(("DN", b, True, nc, {"GA": 0.5 + i}))
    for i, b in enumerate(L[:5]):                        # Oligo: 5 bin 单调升
        rows.append(("OligoLineage", b, True, 40, {"GA": 0.5 + i}))
        for j in (1, 2):
            cube_rows.append(("OligoLineage", b, f"d{j}", {"GA": 0.5 + i}))
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="log")
    bt = bt.drop(columns=[c for c in bt.columns if c.endswith("_GB")])   # 真实 cube: 缺失基因无列
    cube = synth.mk_cube(cfg, cube_rows, ["GA"])
    return spcfg, genes, bt, cube


def test_trend_table_gates_lowpower_donor():
    cfg = _common.cfg()
    spcfg, genes, bt, cube = _build_case(cfg)
    std_min = int(cfg.bins_cfg["min_bins_for_trend"])

    ti = st.trend_table(bt, cube, genes, "Immune", spcfg, cfg)
    assert list(ti.columns) == ["gene", "group", "n_bins", "n_cells_total", "rho_bin",
                                "p_bin", "q_bin", "rho_donor", "p_donor", "n_obs_donor",
                                "low_power", "status"]
    ga = ti[ti["gene"] == "GA"].iloc[0]
    assert ga["status"] == "ok" and ga["n_bins"] == 4
    assert bool(ga["low_power"]) == (4 < std_min)        # override 才通过 → low_power
    assert ga["rho_bin"] == 1.0 and abs(ga["p_bin"] - 2 / 24) < 1e-12   # n=4 精确置换
    assert np.isclose(ga["q_bin"], ga["p_bin"])          # 族内仅 1 个非 NA p → q=p
    assert ga["rho_donor"] == 1.0 and int(ga["n_obs_donor"]) == 8
    gb = ti[ti["gene"] == "GB"].iloc[0]
    assert gb["status"] == "missing" and np.isnan(gb["rho_bin"]) and np.isnan(gb["q_bin"])

    gp = st.trend_table(bt, cube, genes, "Purkinje", spcfg, cfg)
    gp = gp[gp["gene"] == "GA"].iloc[0]
    assert gp["status"] == "insufficient_bins" and np.isnan(gp["rho_bin"]) and gp["low_power"]

    gd = st.trend_table(bt, cube, genes, "DN", spcfg, cfg)
    gd = gd[gd["gene"] == "GA"].iloc[0]
    assert gd["status"] == "insufficient_cells" and not gd["low_power"]

    to = st.trend_table(bt, cube, genes, "OligoLineage", spcfg, cfg)
    go = to[to["gene"] == "GA"].iloc[0]
    assert go["status"] == "ok" and not go["low_power"]
    assert abs(go["p_bin"] - 2 / 120) < 1e-12            # n=5 精确置换
    assert go["rho_donor"] == 1.0 and int(go["n_obs_donor"]) == 10
    p = to["p_bin"].to_numpy(float)
    m = np.isfinite(p)
    assert np.allclose(st.bh_qvalues(p)[m], to["q_bin"].to_numpy(float)[m])


def test_trend_table_bin_order_guard():
    """bin_order 与 spcfg 不一致(如 collapse 敏感性传错 cfg)→ RuntimeError, 不静默。"""
    cfg = _common.cfg()
    spcfg, genes, bt, cube = _build_case(cfg)
    bt2 = bt.copy()
    bt2["bin_order"] = bt2["bin_order"] + 1
    try:
        st.trend_table(bt2, cube, genes, "Immune", spcfg, cfg)
        raise AssertionError("bin_order 不一致未被拦截")
    except RuntimeError as e:
        assert "bin_order" in str(e)


def test_trend_table_deterministic():
    cfg = _common.cfg()
    spcfg, genes, bt, cube = _build_case(cfg)
    t1 = st.trend_table(bt, cube, genes, "OligoLineage", spcfg, cfg)
    t2 = st.trend_table(bt, cube, genes, "OligoLineage", spcfg, cfg)
    pd.testing.assert_frame_equal(t1, t2)                # R8: 两次重跑逐值一致
