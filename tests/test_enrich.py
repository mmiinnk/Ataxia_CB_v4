"""stats_enrich 测试(spec §3.13):perm_enrichment 手算 O/E/log2FE + 无关标签 p=1/近均匀 +
富集格显著 + 恒等排列计入(双侧计数语义);stratified_perm_spearman 在"y 完全由 strata 决定"
时 p_strat==1.0 而 p_naive 显著(混杂消除)+ 层内信号时保持功效;assign_modules 五种路径
(正常/broad/缺失/fraction_na/bin 不足)+ 组级 override;KW/MW/Cliff's δ/Fisher 手算;
同 seed 确定性、异 seed 可区分(R8)。n_perm 全部走 default_rng(seed)。"""
from __future__ import annotations

import numpy as np
import pandas as pd

import _common
import synth  # noqa: F401  (_common.cfg 之外本文件暂不直建 bin_table;保留统一入口)
from atx import stats_enrich as se
from atx import stats_trend as st


def test_cliffs_delta_handcalc():
    assert se.cliffs_delta([1, 2, 3], [4, 5, 6]) == -1.0            # a 全小于 b
    assert se.cliffs_delta([4, 5, 6], [1, 2, 3]) == 1.0             # 正 = a 倾向更大
    assert se.cliffs_delta([1, 1], [1, 1]) == 0.0                   # 全 ties
    assert np.isclose(se.cliffs_delta([1, 2, 3, 4], [3]), -0.25)    # (1−2)/4
    assert np.isclose(se.cliffs_delta([1, np.nan, 3], [2]), 0.0)    # NaN 各自剔除
    assert np.isnan(se.cliffs_delta([], [1])) and np.isnan(se.cliffs_delta([np.nan], [1]))


def test_group_compare_handcalc():
    """三组完全分离 [1,2,3]/[4,5,6]/[7,8,9]: KW H=7.2 → p=exp(−H/2)(df=2 χ² 解析);
    两两 MW p=0.1(3v3 精确双侧最小值);Cliff's δ 全 −1;q=p(全相等 BH 不变);NaN 剔除。"""
    vals = [1.0, 2, 3, 4, 5, 6, 7, 8, 9]
    groups = ["A"] * 3 + ["B"] * 3 + ["C"] * 3
    kw_p, pw = se.group_compare(vals, groups)
    assert np.isclose(kw_p, np.exp(-7.2 / 2), rtol=1e-4)            # chi2.sf(7.2, 2)=e^{−3.6}
    assert list(pw.columns) == ["group1", "group2", "n1", "n2", "p_mw", "q_mw", "cliffs_delta"]
    assert len(pw) == 3
    assert list(pw["group1"]) == ["A", "A", "B"] and list(pw["group2"]) == ["B", "C", "C"]
    assert np.allclose(pw["p_mw"], 0.1, atol=0.02)                  # scipy 精确路 = 0.1
    assert np.allclose(pw["q_mw"], pw["p_mw"])                      # p 全相等 → q=p
    assert np.allclose(pw["cliffs_delta"], -1.0)
    # NaN 值剔除后组样本量不变
    kw2, pw2 = se.group_compare(vals + [np.nan], groups + ["A"])
    assert np.isclose(kw2, kw_p) and list(pw2["n1"]) == [3, 3, 3]
    # 非空组 < 2 → (nan, 空表带列)
    kw3, pw3 = se.group_compare([1.0, 2.0], ["A", "A"])
    assert np.isnan(kw3) and len(pw3) == 0 and list(pw3.columns) == list(pw.columns)


def test_fisher_by_column_handcalc():
    """a=[C1×3,C2×3], b=[m1,m1,m2,m2,m3,m3]: 逐格 2×2 手算。
    (C1,m1): [[2,1],[0,3]] 双侧 p=0.2(左尾)+0.2(右尾)=0.4, OR=inf;
    (C1,m2): [[1,2],[1,2]] OR=1, p=1;BH: 4 个 p=0.4 → q=0.4×6/4=0.6, p=1 → q=1。"""
    a = ["C1", "C1", "C1", "C2", "C2", "C2"]
    b = ["m1", "m1", "m2", "m2", "m3", "m3"]
    df = se.fisher_by_column(a, b)
    assert list(df.columns) == ["a", "b", "n_ab", "n_a", "n_b", "n_total",
                                "odds_ratio", "p_fisher", "q_fisher"]
    assert len(df) == 6 and set(df["n_total"]) == {6}
    c1m1 = df[(df["a"] == "C1") & (df["b"] == "m1")].iloc[0]
    assert c1m1["n_ab"] == 2 and c1m1["n_a"] == 3 and c1m1["n_b"] == 2
    assert np.isinf(c1m1["odds_ratio"]) and np.isclose(c1m1["p_fisher"], 0.4)
    c1m2 = df[(df["a"] == "C1") & (df["b"] == "m2")].iloc[0]
    assert np.isclose(c1m2["odds_ratio"], 1.0) and np.isclose(c1m2["p_fisher"], 1.0)
    assert np.isclose(c1m1["q_fisher"], 0.6)                        # BH: 0.4×6/4 → 0.6
    assert np.isclose(df[(df["p_fisher"] == 1.0)]["q_fisher"].iloc[0], 1.0)
    pd.testing.assert_frame_equal(df, se.fisher_by_column(a, b))    # 无 RNG, 确定性
    try:
        se.fisher_by_column(["C1", np.nan], ["m1", "m2"])
        raise AssertionError("NA 标签未被拦截")
    except ValueError as e:
        assert "NA" in str(e)


def test_perm_enrichment_handcalc_and_errors():
    """a=[A,A,B,B], b=[m1,m1,m2,m2]: O/E/log2FE 手算;p ≥ 1/n_perm(恒等计入);q ≥ p。"""
    cfg = _common.cfg()
    a = ["A", "A", "B", "B"]
    b = ["m1", "m1", "m2", "m2"]
    df = se.perm_enrichment(a, b, n_perm=99, seed=cfg.seed)
    assert list(df.columns) == ["a", "b", "obs", "exp", "log2fe", "p", "q"]
    assert len(df) == 4 and list(df["a"]) == ["A", "A", "B", "B"]   # 升序全格
    r = df[(df["a"] == "A") & (df["b"] == "m1")].iloc[0]
    assert r["obs"] == 2 and r["exp"] == 1.0
    assert np.isclose(r["log2fe"], np.log2(2.5 / 1.5))
    r2 = df[(df["a"] == "A") & (df["b"] == "m2")].iloc[0]
    assert r2["obs"] == 0 and np.isclose(r2["log2fe"], np.log2(0.5 / 1.5))
    assert (df["p"] >= 1 / 99 - 1e-12).all()
    assert (df["q"] >= df["p"] - 1e-12).all() and (df["q"] <= 1.0).all()
    pd.testing.assert_frame_equal(df, se.perm_enrichment(a, b, n_perm=99, seed=cfg.seed))
    try:
        se.perm_enrichment(["A", "B"], ["m1"], n_perm=9, seed=1)
        raise AssertionError("长度不一致未被拦截")
    except ValueError:
        pass
    try:
        se.perm_enrichment(["A", np.nan], ["m1", "m2"], n_perm=9, seed=1)
        raise AssertionError("NA 标签未被拦截")
    except ValueError as e:
        assert "NA" in str(e)
    try:
        se.perm_enrichment(["A"], ["m1"], n_perm=0, seed=1)
        raise AssertionError("n_perm=0 未被拦截")
    except ValueError:
        pass


def test_perm_enrichment_null_gives_p_one():
    """零关联且 O==E(a=[A,A,B,B], b=[m1,m2,m1,m2]): 每格 log2FE=0 →
    任何排列 |lfe_perm| ≥ 0−1e-12 恒真 → p 恰=1.0(双侧计数语义的精确锚)。"""
    cfg = _common.cfg()
    df = se.perm_enrichment(["A", "A", "B", "B"], ["m1", "m2", "m1", "m2"],
                            n_perm=99, seed=cfg.seed)
    assert np.allclose(df["log2fe"], 0.0) and np.allclose(df["obs"], 1.0)
    assert np.allclose(df["p"], 1.0) and np.allclose(df["q"], 1.0)


def test_perm_enrichment_null_uniform():
    """§3.13: 固定种子下无关标签(400 基因 × 4 类 × 6 模块, 独立随机)p 接近均匀。
    断言: p ≥ 1/n_perm;P(p≤0.05) 远低于名义外扩界;mean(p)≈0.5;两跑逐值一致。"""
    cfg = _common.cfg()
    rng = np.random.default_rng(cfg.seed)
    n = 400
    a = np.array([f"C{v}" for v in rng.integers(0, 4, n)])
    b = np.array([f"M{v}" for v in rng.integers(0, 6, n)])
    df = se.perm_enrichment(a, b, n_perm=999, seed=cfg.seed)
    assert len(df) == 24
    assert (df["p"] >= 1 / 999 - 1e-12).all()
    frac05 = float((df["p"] <= 0.05).mean())
    assert frac05 <= 0.25, f"无关标签 p≤0.05 比例过高: {frac05}"   # 期望 ~0.05, 宽界防单 seed 波动
    assert 0.3 <= float(df["p"].mean()) <= 0.7
    pd.testing.assert_frame_equal(df, se.perm_enrichment(a, b, n_perm=999, seed=cfg.seed))


def test_perm_enrichment_enriched_cells():
    """完美对应(4 类 × 4 模块各 25 基因一一对齐): 对角格 log2FE=log2(25.5/6.75) 手算;
    off-diagonal obs=0 → lfe=log2(0.5/6.75)=−3.76,富集臂无排列可达((O+.5)/6.75 ≥ 13.5
    需 O ≥ 90 > 25),仅恒等与罕见 O=0 排列(P≈8e-4)计数 → p ≤ 5/n_perm;
    对角 p < 0.05(双侧含耗竭臂 P(O≤1)≈1%),q < 0.05(16 格 BH)。"""
    cfg = _common.cfg()
    a = np.array([f"C{i // 25}" for i in range(100)])
    b = np.array([f"M{i // 25}" for i in range(100)])
    df = se.perm_enrichment(a, b, n_perm=999, seed=cfg.seed)
    diag = df[df["a"].str[1] == df["b"].str[1]]
    off = df[df["a"].str[1] != df["b"].str[1]]
    assert len(diag) == 4 and len(off) == 12
    assert np.allclose(diag["obs"], 25) and np.allclose(diag["exp"], 6.25)
    assert np.allclose(diag["log2fe"], np.log2(25.5 / 6.75))
    assert np.allclose(off["obs"], 0)
    assert (off["p"] >= 1 / 999 - 1e-12).all() and (off["p"] <= 5 / 999).all()
    assert (diag["p"] < 0.05).all() and (diag["p"] >= 1 / 999 - 1e-12).all()
    assert (diag["q"] < 0.05).all() and (off["q"] < 0.05).all()


def test_stratified_perm_confounding():
    """§3.13 核心断言: y 完全由 strata 决定 → 层内打乱 y 不变 → 每个置换 rho==rho_obs
    → p_strat 恰=1.0;而 x 与 y 的边缘关联真实存在(rho≈0.943 手算)→ p_naive 显著。
    "显著"完全来自混杂, 分层置换正确识别出无可检验的层内信号。"""
    cfg = _common.cfg()
    strata = np.repeat(["S1", "S2", "S3"], 10)
    y = np.concatenate([np.full(10, v) for v in (0.0, 1.0, 2.0)])
    x = y * 2 + 0.1 * np.arange(30, dtype=float)
    rho, p_strat, p_naive = se.stratified_perm_spearman(x, y, strata, n_perm=999, seed=cfg.seed)
    # 手算: rx=1..30, ry=[5.5×10,15.5×10,25.5×10] → rho = 2000/sqrt(2247.5×2000)
    assert np.isclose(rho, 2000 / np.sqrt(2247.5 * 2000))
    assert p_strat == 1.0                                            # 恒等+全部置换同值
    assert p_naive < 0.01                                            # 未分层 → 假显著
    # 退化: n<3 / 恒定列 → (nan,)*3;NaN 对剔除
    assert all(np.isnan(v) for v in se.stratified_perm_spearman([1, 2], [1, 2], ["a", "b"], 99))
    assert all(np.isnan(v) for v in se.stratified_perm_spearman(x, np.zeros(30), strata, 99))
    xm = x.copy(); xm[0] = np.nan
    rho2, ps2, pn2 = se.stratified_perm_spearman(xm, y, strata, n_perm=99, seed=cfg.seed)
    assert np.isfinite(rho2) and ps2 == 1.0 and np.isfinite(pn2)


def test_stratified_perm_within_signal():
    """层内真信号(x 与 y 在每层内同升, 层间无 x 差异): p_strat 仍显著——
    分层置换消除混杂但不损失层内功效(与上一测试对偶)。"""
    cfg = _common.cfg()
    strata = np.repeat(["S1", "S2", "S3"], 10)
    x = np.tile(np.arange(10, dtype=float), 3)
    y = np.concatenate([np.full(10, v) for v in (0.0, 1.0, 2.0)]) + 0.5 * x
    rho, p_strat, p_naive = se.stratified_perm_spearman(x, y, strata, n_perm=999, seed=cfg.seed)
    assert rho > 0.5
    assert p_strat < 0.05 and p_naive < 0.05


def test_stratified_perm_seed_determinism():
    """R8 + spec §3.9 default_rng: 同 seed 两跑逐值一致;真信号 rho 与 seed 无关且两路
    p 均显著;零信号(x⊥y)时置换计数依赖 seed → p 可区分(异 seed 至少一路不同)。"""
    cfg = _common.cfg()
    rng = np.random.default_rng(cfg.seed)
    n = 200
    x = rng.normal(size=n)
    y = 0.3 * x + rng.normal(size=n)
    strata = np.array([f"S{v}" for v in rng.integers(0, 5, n)])
    r1 = se.stratified_perm_spearman(x, y, strata, n_perm=499, seed=cfg.seed)
    r2 = se.stratified_perm_spearman(x, y, strata, n_perm=499, seed=cfg.seed)
    assert r1 == r2                                                  # 同 seed 逐值一致
    assert r1[1] < 0.05 and r1[2] < 0.05                             # 真信号两路都显著
    xn = rng.normal(size=n)                                          # 零信号: x ⊥ y
    s0 = se.stratified_perm_spearman(xn, y, strata, n_perm=499, seed=cfg.seed)
    s1 = se.stratified_perm_spearman(xn, y, strata, n_perm=499, seed=cfg.seed)
    s2 = se.stratified_perm_spearman(xn, y, strata, n_perm=499, seed=cfg.seed + 1)
    assert s0 == s1 and s0[0] == s2[0]                               # rho 与 seed 无关
    assert (s0[1] != s2[1]) or (s0[2] != s2[2])                      # 置换流依赖 seed
    assert s0[1] > 0.05 and s0[2] > 0.05                             # 零信号不显著


def test_assign_modules_paths():
    """五条路径: 正常(G*×P*)/ broad(max_fraction<阈值 → ALL 峰期)/ 缺失基因 /
    fraction_na(k<4 NA 政策 → 模块 NA 非 broad)/ G* 可行 bin 不足(含组级 override 放行)。"""
    cfg = _common.cfg()
    thr = float(cfg.stats_cfg["enrichment"]["broad_fraction_threshold"])
    std_min = int(cfg.bins_cfg["min_bins_for_trend"])
    imm_min = int(cfg.bins_cfg["per_group_override"]["Immune"]["min_bins_for_trend"])
    assert imm_min < std_min                                         # 场景有效性(4<5)
    hi = min(thr * 2, 1.0)
    spec = pd.DataFrame({
        "gene":   ["G1", "G2", "G3", "G4", "G5", "G6"],
        "status": ["ok", "ok", "missing", "ok", "ok", "ok"],
        "max_group": ["Purkinje", "GranuleLineage", "NA", "Immune", "DN", "Immune"],
        "max_fraction": [hi, thr / 2, np.nan, np.nan, hi, hi],
    })
    pk = pd.DataFrame({
        "gene":   ["G1", "G2", "G4", "G5", "G6"],
        "group":  ["Purkinje", "ALL", "Immune", "DN", "Immune"],
        "peak_bin": ["P0", "P4", "E12.5", "E13.5", "E13.5"],
        "peak_period": ["perinatal", "postnatal", "embryonic", "embryonic", "embryonic"],
        "n_feasible_bins": [std_min + 1, std_min + 3, std_min + 1, std_min - 1, imm_min],
    })
    df = se.assign_modules(spec, pk, cfg).set_index("gene")
    assert list(se.assign_modules(spec, pk, cfg).columns) == [
        "gene", "max_group", "max_fraction", "g_star", "peak_group_used", "peak_bin",
        "peak_period", "n_feasible_bins", "module", "status"]
    assert df.loc["G1", "module"] == "Purkinje|perinatal" and df.loc["G1", "status"] == "ok"
    assert df.loc["G2", "g_star"] == "broad" and df.loc["G2", "peak_group_used"] == "ALL"
    assert df.loc["G2", "module"] == "broad|postnatal"               # P* 取 ALL 峰期(口径决策)
    assert df.loc["G3", "module"] == "NA" and df.loc["G3", "status"] == "missing_gene"
    assert df.loc["G4", "module"] == "NA" and df.loc["G4", "status"] == "fraction_na"
    assert df.loc["G4", "g_star"] == "Immune"                        # NA 不当 broad(NA 政策)
    assert df.loc["G5", "module"] == "NA" and df.loc["G5", "status"] == "insufficient_bins_gstar"
    assert df.loc["G6", "module"] == "Immune|embryonic"              # override=4 放行(与 trend 一致)
    pd.testing.assert_frame_equal(df.reset_index(), se.assign_modules(spec, pk, cfg))
