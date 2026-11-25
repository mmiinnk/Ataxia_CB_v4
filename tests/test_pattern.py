"""stats_pattern 测试(spec §3.13):内部峰/端点峰/中间缺口三条轨迹 + window_high(人 16 bin)
断言标签;优先级与 ≥/< 边界(cfg 突变法, 承 test_specificity);half_max 跨缺口;ties 取最早;
缺失/无可行 bin → NA;label 缺列 ValueError;确定性(R8)。阈值全部从 config 读。"""
from __future__ import annotations

import numpy as np
import pandas as pd

import _common
import synth
from atx import stats_pattern as spm


def _lab(pat, cfg, rho=np.nan, q=np.nan) -> pd.Series:
    """label() 输入 = patterns 产物 merge trend 列(s02 的真实调用形态)。"""
    df = pat.copy()
    df["rho_bin"] = rho
    df["q_bin"] = q
    return spm.label(df, cfg)


def test_internal_peak_transient():
    """轨迹①内部峰: 4 可行 bin 尖峰 → transient_peak;pvf/period_ratio 手算。"""
    cfg = _common.cfg()
    L = cfg.species("mouse").labels
    genes = ["GA"]
    rows = [("Purkinje", L[i], True, 40, {"GA": v})
            for i, v in zip(range(2, 6), (1.0, 1.0, 20.0, 1.0))]   # E12.5–E15.5 同期
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="cp")
    r = spm.patterns(bt, genes, "Purkinje", cfg).iloc[0]
    assert r["status"] == "ok" and r["n_feasible_bins"] == 4
    assert r["peak_bin"] == L[4] and r["peak_bin_order"] == 5.0     # bin_order = idx+1
    assert r["peak_period"] == cfg.species("mouse").period_of[L[4]] == "embryonic"
    assert bool(r["peak_internal"])                                  # 0 < 2 < 3
    pvf_exp = np.log1p(20.0) - np.log1p(1.0)                        # 峰 − max(左均, 右均)=log1p(1)
    assert np.isclose(r["peak_vs_flank"], pvf_exp)
    assert np.isclose(r["period_ratio"], 1.0)                       # 全部 bin 同 period
    assert np.isclose(r["peak_cp"], 20.0)
    tp = cfg.stats_cfg["temporal_pattern"]
    assert pvf_exp >= float(tp["pulse_min_peak_vs_flank"])          # 场景有效性
    assert _lab(spm.patterns(bt, genes, "Purkinje", cfg), cfg).iloc[0] == "transient_peak"


def test_endpoint_high_and_priority():
    """轨迹②端点峰: 峰在末可行 bin → peak_internal False、pvf NA;
    period_ratio ≥ 2 → endpoint_high;优先级/边界: monotonic(q<mono_q 且 |rho|≥mono_rho)压过
    endpoint_high, q=mono_q(不含等号)或 |rho|<mono_rho 回落 endpoint_high。"""
    cfg = _common.cfg()
    tp = cfg.stats_cfg["temporal_pattern"]
    mq, mr = float(tp["mono_q"]), float(tp["mono_rho"])
    L = cfg.species("mouse").labels
    genes = ["GA"]
    rows = [("Purkinje", L[i], True, 40, {"GA": v})
            for i, v in zip(range(4, 9), (1.0, 1.0, 1.0, 1.0, 10.0))]  # E14.5–P4
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="cp")
    p = spm.patterns(bt, genes, "Purkinje", cfg)
    r = p.iloc[0]
    assert r["peak_bin"] == L[8] and not bool(r["peak_internal"])   # 末端峰 = 删失单调
    assert np.isnan(r["peak_vs_flank"])                             # 右侧无可行 bin
    ratio_exp = 10.0 / ((1 + 1 + 1 + 1 + 10) / 5)                   # postnatal 均值/全轴均值
    assert np.isclose(r["period_ratio"], ratio_exp) and ratio_exp >= float(tp["window_min_ratio"])
    assert r["peak_period"] == "postnatal"
    assert _lab(p, cfg).iloc[0] == "endpoint_high"                  # 无 trend 信息
    assert _lab(p, cfg, rho=0.9, q=mq * 0.2).iloc[0] == "monotonic_up"      # 优先级压过
    assert _lab(p, cfg, rho=-0.9, q=mq * 0.2).iloc[0] == "monotonic_down"
    assert _lab(p, cfg, rho=0.9, q=mq).iloc[0] == "endpoint_high"   # q < mono_q 严格小于
    assert _lab(p, cfg, rho=mr - 0.1, q=mq * 0.2).iloc[0] == "endpoint_high"  # |rho| ≥ mono_rho


def test_middle_gap_internal_halfmax():
    """轨迹③中间缺口: 不可行 bin 夹在中间, 内部性按可行 bin 序列判定(§3.8 差异①);
    half_max_window 沿可行序列延展、跨缺口(差异②);peak_vs_flank/period_ratio 手算。"""
    cfg = _common.cfg()
    L = cfg.species("mouse").labels
    genes = ["GA"]
    vals = {2: 1.0, 3: 12.0, 4: np.nan, 5: 20.0, 6: 1.0, 7: 1.0}   # L4=E14.5 不可行
    rows = [("Purkinje", L[i], i != 4, 40 if i != 4 else 5, {"GA": vals[i]})
            for i in range(2, 8)]
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="cp")
    r = spm.patterns(bt, genes, "Purkinje", cfg).iloc[0]
    assert r["n_feasible_bins"] == 5                                # L2,L3,L5,L6,L7
    assert r["peak_bin"] == L[5] and bool(r["peak_internal"])       # 序列位 2/5, 缺口不算侧
    pvf_exp = np.log1p(20.0) - (np.log1p(1.0) + np.log1p(12.0)) / 2  # 左均=(L2,L3), 右均=(L6,L7)
    assert np.isclose(r["peak_vs_flank"], pvf_exp)
    ratio_exp = ((1 + 12 + 20) / 3) / ((1 + 12 + 20 + 1 + 1) / 5)   # embryonic(L2,L3,L5)/全轴
    assert np.isclose(r["period_ratio"], ratio_exp)
    pt = spm.peak_table(bt, genes, cfg)
    pk = pt[(pt["group"] == "Purkinje") & (pt["gene"] == "GA")].iloc[0]
    assert pk["half_max_lo"] == L[3] and pk["half_max_hi"] == L[5]  # 12 ≥ 0.5×20 → 跨缺口 L4
    assert pk["peak_bin"] == L[5] and np.isclose(pk["peak_cp"], 20.0)


def test_window_high_human():
    """window_high(内部峰、pvf<0.5、期富集≥2): 人 16 bin, fetal 平台期 9–10、余期 0.5。
    优先级: window_high 压过 monotonic(即便 q/rho 显著)。全部期望值手算。"""
    cfg = _common.cfg()
    sp = cfg.species("human")
    L = sp.labels
    assert (sp.period_of[L[4]] == "fetal"                           # 场景有效性(真 config)
            and sum(p == "fetal" for p in sp.period_of.values()) == 4
            and sum(p == "embryonic" for p in sp.period_of.values()) == 3
            and sum(p == "infant" for p in sp.period_of.values()) == 5
            and sum(p == "adult" for p in sp.period_of.values()) == 4)
    genes = ["GA"]
    cp = [9.0] * 3 + [9.0, 10.0, 9.0, 9.0] + [0.5] * 5 + [0.5] * 4  # 峰=11w(fetal 第 2 bin)
    rows = [("Purkinje", L[i], True, 40, {"GA": v}) for i, v in enumerate(cp)]
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="cp", species="human")
    p = spm.patterns(bt, genes, "Purkinje", cfg)
    r = p.iloc[0]
    assert r["peak_bin"] == L[4] and bool(r["peak_internal"])
    left = np.log1p(9.0)                                            # 左侧 4 bin 全 9
    right = (2 * np.log1p(9.0) + 9 * np.log1p(0.5)) / 11            # 右 11 bin: 2×9 + 9×0.5
    pvf_exp = np.log1p(10.0) - max(left, right)
    assert np.isclose(r["peak_vs_flank"], pvf_exp)
    assert pvf_exp < float(cfg.stats_cfg["temporal_pattern"]["pulse_min_peak_vs_flank"])
    ratio_exp = ((9 + 10 + 9 + 9) / 4) / (68.5 / 16)                # fetal 均值/全轴均值
    assert np.isclose(r["period_ratio"], ratio_exp)
    assert ratio_exp >= float(cfg.stats_cfg["temporal_pattern"]["window_min_ratio"])
    assert _lab(p, cfg).iloc[0] == "window_high"
    assert _lab(p, cfg, rho=0.9, q=0.001).iloc[0] == "window_high"  # 优先级压过 monotonic


def test_ties_flat_missing_nofeasible():
    """ties → 峰取最早可行 bin;恒定轨迹 → flat;缺失基因 → missing/NA;
    组内无可行 bin → no_feasible_bin/NA(NA 而非 flat, 承 specific_call 语义)。"""
    cfg = _common.cfg()
    L = cfg.species("mouse").labels
    genes = ["GA", "GB"]
    rows = [("Purkinje", L[i], True, 40, {"GA": 5.0, "GB": 5.0}) for i in range(2, 6)]
    rows += [("DN", L[i], False, 5, {"GA": np.nan}) for i in range(2, 6)]
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="cp")
    bt = bt.drop(columns=[c for c in bt.columns if c.endswith("_GB")])   # 真实 cube: 缺失基因无列
    p = spm.patterns(bt, genes, "Purkinje", cfg).set_index("gene")
    assert p.loc["GA", "peak_bin"] == L[2]                          # 全 5.0 ties → 最早
    assert not bool(p.loc["GA", "peak_internal"]) and np.isnan(p.loc["GA", "peak_vs_flank"])
    assert np.isclose(p.loc["GA", "period_ratio"], 1.0)
    assert p.loc["GB", "status"] == "missing" and p.loc["GB", "n_feasible_bins"] == 0
    lab = _lab(spm.patterns(bt, genes, "Purkinje", cfg), cfg).set_axis(p.index)
    assert lab.loc["GA"] == "flat" and lab.loc["GB"] == "NA"
    d = spm.patterns(bt, ["GA"], "DN", cfg).iloc[0]
    assert d["status"] == "no_feasible_bin"
    assert _lab(spm.patterns(bt, ["GA"], "DN", cfg), cfg).iloc[0] == "NA"


def test_label_boundaries():
    """≥/< 边界(cfg 突变法): pvf 恰=pulse → transient, +1e-9 → flat;
    ratio 恰=window_min_ratio → window_high, +1e-9 → flat。"""
    cfg = _common.cfg()
    L = cfg.species("mouse").labels
    genes = ["GA"]
    rows = [("Purkinje", L[i], True, 40, {"GA": v})
            for i, v in zip(range(2, 6), (1.0, 1.0, 20.0, 1.0))]
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="cp")
    p = spm.patterns(bt, genes, "Purkinje", cfg)
    pvf = float(p.iloc[0]["peak_vs_flank"])
    cfg_eq = _common.cfg()
    cfg_eq.stats_cfg["temporal_pattern"]["pulse_min_peak_vs_flank"] = pvf
    assert _lab(p, cfg_eq).iloc[0] == "transient_peak"              # ≥ 含等号
    cfg_gt = _common.cfg()
    cfg_gt.stats_cfg["temporal_pattern"]["pulse_min_peak_vs_flank"] = pvf + 1e-9
    assert _lab(p, cfg_gt).iloc[0] == "flat"                        # ratio=1 <2, 内部峰 → flat

    sp = cfg.species("human")
    Lh = sp.labels
    cp = [9.0] * 3 + [9.0, 10.0, 9.0, 9.0] + [0.5] * 5 + [0.5] * 4
    rows_h = [("Purkinje", Lh[i], True, 40, {"GA": v}) for i, v in enumerate(cp)]
    bth = synth.mk_bin_table(cfg, rows_h, genes, value_kind="cp", species="human")
    ph = spm.patterns(bth, genes, "Purkinje", cfg)
    ratio = float(ph.iloc[0]["period_ratio"])
    cfg_eq2 = _common.cfg()
    cfg_eq2.stats_cfg["temporal_pattern"]["window_min_ratio"] = ratio
    assert _lab(ph, cfg_eq2).iloc[0] == "window_high"               # ≥ 含等号
    cfg_gt2 = _common.cfg()
    cfg_gt2.stats_cfg["temporal_pattern"]["window_min_ratio"] = ratio + 1e-9
    assert _lab(ph, cfg_gt2).iloc[0] == "flat"                      # pvf<pulse → 非 transient


def test_label_requires_columns():
    cfg = _common.cfg()
    L = cfg.species("mouse").labels
    bt = synth.mk_bin_table(cfg, [("Purkinje", L[2], True, 40, {"GA": 1.0})], ["GA"],
                            value_kind="cp")
    p = spm.patterns(bt, ["GA"], "Purkinje", cfg)                   # 无 rho_bin/q_bin
    try:
        spm.label(p, cfg)
        raise AssertionError("缺 rho_bin/q_bin 未被拦截")
    except ValueError as e:
        assert "rho_bin" in str(e)


def test_peak_table_structure_determinism():
    """peak_table 列=spec §3.8 差异③;组序=plot_groups ∩ 实有;两跑逐值一致(R8)。"""
    cfg = _common.cfg()
    L = cfg.species("mouse").labels
    genes = ["GA", "GB"]
    rows = []
    for grp in ("ALL", "GranuleLineage", "Purkinje"):               # 故意乱序入表
        for i in range(2, 7):
            rows.append((grp, L[i], True, 40, {"GA": float(i), "GB": 3.0}))
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="cp")
    pt = spm.peak_table(bt, genes, cfg)
    assert list(pt.columns) == ["gene", "group", "peak_bin", "peak_bin_order", "peak_period",
                                "peak_cp", "n_feasible_bins", "half_max_lo", "half_max_hi",
                                "status"]
    groups_seen = list(dict.fromkeys(pt["group"]))
    expect = [g for g in cfg.plot_groups() if g in {"ALL", "GranuleLineage", "Purkinje"}]
    assert groups_seen == expect                                    # plot_groups 序, 非入表序
    assert len(pt) == 3 * len(genes)
    pk = pt[(pt["group"] == "Purkinje") & (pt["gene"] == "GA")].iloc[0]
    # GA 单调升 cp=2,3,4,5,6(L2..L6): 峰=L6, thr=0.5×6=3 → 左走 L5(5)→L4(4)→L3(3, ≥ 含等号),
    # L2(2)<3 停 → 窗=[L3, L6]
    assert pk["peak_bin"] == L[6] and pk["half_max_lo"] == L[3] and pk["half_max_hi"] == L[6]
    gb = pt[(pt["group"] == "ALL") & (pt["gene"] == "GB")].iloc[0]
    assert gb["peak_bin"] == L[2] and gb["half_max_lo"] == L[2] and gb["half_max_hi"] == L[6]
    pd.testing.assert_frame_equal(pt, spm.peak_table(bt, genes, cfg))
    pd.testing.assert_frame_equal(spm.patterns(bt, genes, "Purkinje", cfg),
                                  spm.patterns(bt, genes, "Purkinje", cfg))
