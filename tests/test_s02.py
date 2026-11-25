"""s02 六表合成端到端(spec §5):结构与行数 / 缺失基因与不可行组的 NA 路径 /
global 特异性手算(tau/max_fraction/specific_call)/ by-period groups_used /
emergence 手算(含 left_censored)/ 标签 merge 一致性 + q 族=species×group /
写盘 roundtrip(临时目录, 不触真实 results/)+ 确定性。期望值全部手算。"""
from __future__ import annotations

import os
import tempfile

import numpy as np
import pandas as pd

import _common
import s02_stats
import synth
from atx import stats_trend as st

GENES_IN_BT = ["GT", "GU"]          # GM 故意不入 bt 列 → 缺失基因路径
UNIVERSE = ["GT", "GU", "GM"]
# L 为 0 基 list(L[0]=E10.5):L[1]..L[5] = E11.5..E15.5 全 embryonic;L[6]=E17.5 perinatal
GT_PKJ = {1: 1.0, 2: 1.0, 3: 20.0, 4: 1.0, 5: 1.0}   # 内部尖峰(L[3]=E13.5)


def _setup():
    cfg = _common.cfg()
    L = cfg.species("mouse").labels
    rows = []
    for grp in ("Purkinje", "GranuleLineage", "Astro", "Immune", "ALL"):
        for i in range(1, 6):
            vals = {"GT": GT_PKJ[i] if grp in ("Purkinje", "ALL") else 1.0, "GU": 1.0}
            rows.append((grp, L[i], True, 40, vals))
    rows.append(("GranuleLineage", L[6], True, 40, {"GT": 1.0, "GU": 1.0}))  # perinatal 仅 GL 可行
    for i in range(1, 6):
        rows.append(("Interneuron", L[i], False, 5, {"GT": np.nan}))          # 组实有但全不可行
    bt = synth.mk_bin_table(cfg, rows, GENES_IN_BT, value_kind="cp")
    return cfg, L, bt


def _groups_present(cfg, bt):
    return [g for g in cfg.plot_groups() if g in set(bt["group"])]


def test_build_tables_structure():
    """7 个键、行数 = 组/期/bin × 基因全集、列结构与 §5/§3 一致。"""
    cfg, L, bt = _setup()
    tabs = s02_stats.build_tables(bt, None, UNIVERSE, cfg.species("mouse"), cfg)
    assert list(tabs) == s02_stats.TABLE_NAMES
    gs = _groups_present(cfg, bt)
    assert gs == ["Purkinje", "GranuleLineage", "Interneuron", "Astro", "Immune", "ALL"]
    n_g, n_per = len(gs), 4                                    # 鼠 4 期全出(空期 → missing 行)
    assert len(tabs["trend_by_group"]) == n_g * 3
    assert len(tabs["temporal_patterns"]) == n_g * 3
    assert list(tabs["temporal_patterns"].columns) == s02_stats._TP_COLS
    assert len(tabs["specificity_global"]) == 3
    assert len(tabs["specificity_by_period"]) == n_per * 3
    assert list(tabs["specificity_by_period"].columns)[0:2] == ["gene", "period"]
    assert len(tabs["tau_by_bin"]) == 6 * 3                    # L[1]..L[6] 共 6 bin
    assert len(tabs["emergence"]) == 3
    assert len(tabs["peak_by_group"]) == n_g * 3
    # cube=None → donor 级全 NA(列仍须在, §5 trend_by_group 列集)
    assert tabs["trend_by_group"]["rho_donor"].isna().all()
    for c in ("rho_bin", "p_bin", "q_bin", "rho_donor", "p_donor", "n_bins", "low_power"):
        assert c in tabs["trend_by_group"].columns


def test_missing_gene_and_infeasible_group():
    """GM(无列)→ 各表 missing/NA;Interneuron(全不可行)→ trend insufficient_bins、
    标签 NA(而非 flat)——NA 不当 False。"""
    cfg, L, bt = _setup()
    tabs = s02_stats.build_tables(bt, None, UNIVERSE, cfg.species("mouse"), cfg)
    tr = tabs["trend_by_group"].set_index(["gene", "group"])
    assert tr.loc[("GM", "Purkinje"), "status"] == "missing"
    assert tr.loc[("GT", "Interneuron"), "status"] == "insufficient_bins"
    assert tr.loc[("GT", "Interneuron"), "n_bins"] == 0
    tp = tabs["temporal_patterns"].set_index(["gene", "group"])
    assert tp.loc[("GM", "Purkinje"), "status"] == "missing"
    assert tp.loc[("GM", "Purkinje"), "temporal_pattern"] == "NA"
    assert tp.loc[("GT", "Interneuron"), "status"] == "no_feasible_bin"
    assert tp.loc[("GT", "Interneuron"), "temporal_pattern"] == "NA"
    sg = tabs["specificity_global"].set_index("gene")
    assert sg.loc["GM", "status"] == "missing" and sg.loc["GM", "specific_call"] == "NA"
    pk = tabs["peak_by_group"].set_index(["gene", "group"])
    assert pk.loc[("GM", "Purkinje"), "status"] == "missing"
    tb = tabs["tau_by_bin"]
    assert (tb.loc[tb["gene"] == "GM", "status"] == "missing").all()


def test_global_specificity_handcalc():
    """GT: agg(Purkinje)=4.8, 余组=1(GL 含 L[6] 仍=1)→ tau/max_fraction/specific 手算;
    GU 恒定 → tau=0、ties argmax 取主组序首个、call=broad;fraction_8group: k=4<8 → 全 NA。"""
    cfg, L, bt = _setup()
    tabs = s02_stats.build_tables(bt, None, UNIVERSE, cfg.species("mouse"), cfg)
    sg = tabs["specificity_global"].set_index("gene")
    pc = float(cfg.stats_cfg["specificity"]["pseudocount"])
    tau_exp = 3 * (1 - (1.0 + pc) / (4.8 + pc)) / 3            # Σ(1−x/xmax)/(k−1), k=4
    frac_exp = (4.8 + pc) / (4.8 + pc + 3 * (1.0 + pc))
    r = sg.loc["GT"]
    assert r["status"] == "ok" and r["n_groups_used"] == 4
    assert np.isclose(r["tau"], tau_exp) and np.isclose(r["max_fraction"], frac_exp)
    assert r["max_group"] == "Purkinje" and r["specific_call"] == "specific"     # tau≥0.6, frac≥0.5
    u = sg.loc["GU"]
    assert np.isclose(u["tau"], 0.0) and np.isclose(u["max_fraction"], 0.25)
    assert u["max_group"] == "Purkinje" and u["specific_call"] == "broad"        # ties → 主组序首个
    f8 = [c for c in sg.columns if c.startswith("fraction_8group_")]
    assert len(f8) == 8 and sg.loc[["GT", "GU"], f8].isna().all().all()          # 严格版 8 组全可行才非 NA


def test_by_period_groups_used():
    """perinatal 仅 GL 可行 → k=1<4: tau/fraction NA + groups_used 字符串(§3.7 差异②);
    embryonic 行 = global 值(L[6] 不改任何组的期均值);postnatal 无数据 → missing。"""
    cfg, L, bt = _setup()
    tabs = s02_stats.build_tables(bt, None, UNIVERSE, cfg.species("mouse"), cfg)
    bp = tabs["specificity_by_period"]
    gt = bp[bp["gene"] == "GT"].set_index("period")
    per = gt.loc["perinatal"]
    assert per["n_groups_used"] == 1 and per["groups_used"] == "GranuleLineage"
    assert np.isnan(per["tau"]) and per["specific_call"] == "NA"
    emb = gt.loc["embryonic"]
    pc = float(cfg.stats_cfg["specificity"]["pseudocount"])
    assert np.isclose(emb["tau"], 3 * (1 - (1.0 + pc) / (4.8 + pc)) / 3)
    assert emb["specific_call"] == "specific"
    assert gt.loc["postnatal", "status"] == "missing" and gt.loc["adult", "status"] == "missing"
    # period 顺序 = 首个 bin 出现序(确定性)
    assert list(dict.fromkeys(bp["period"])) == ["embryonic", "perinatal", "postnatal", "adult"]


def test_emergence_handcalc():
    """GT: frac_Purkinje(L4)=(20+pc)/(20+pc+3(1+pc))≥0.6 → emergence=L4, 左截断 False;
    GU 恒 0.25 → NA;GM max_group=NA → 全 NA 行;L[6](k=1)不入可算 bin。"""
    cfg, L, bt = _setup()
    tabs = s02_stats.build_tables(bt, None, UNIVERSE, cfg.species("mouse"), cfg)
    em = tabs["emergence"].set_index("gene")
    pc = float(cfg.stats_cfg["specificity"]["pseudocount"])
    g = em.loc["GT"]
    assert g["max_group"] == "Purkinje"
    assert g["first_computable_bin"] == L[1]
    assert g["emergence_bin"] == L[3] and g["emergence_bin_order"] == 4.0
    assert np.isclose(g["frac_at_emergence"], (20 + pc) / (20 + pc + 3 * (1 + pc)))
    assert not bool(g["left_censored"])
    u = em.loc["GU"]
    assert u["emergence_bin"] == "NA" and np.isnan(u["emergence_bin_order"])
    m = em.loc["GM"]
    assert m["max_group"] == "NA" and m["first_computable_bin"] == "NA"
    assert m["emergence_bin"] == "NA" and not bool(m["left_censored"])


def test_label_merge_and_qfamily():
    """GT×Purkinje = transient_peak 且 rho_bin/q_bin 与 trend 表逐值一致(merge 不串组);
    q_bin = 组内 BH(族=species×group, 逐组重算比对);fdr_family 改口径 → NotImplementedError。"""
    cfg, L, bt = _setup()
    tabs = s02_stats.build_tables(bt, None, UNIVERSE, cfg.species("mouse"), cfg)
    tp = tabs["temporal_patterns"]
    tr = tabs["trend_by_group"]
    row = tp[(tp["gene"] == "GT") & (tp["group"] == "Purkinje")].iloc[0]
    assert row["temporal_pattern"] == "transient_peak"          # pvf=log1p(20)−log1p(1)≥0.5, 内部峰
    assert row["peak_bin"] == L[3] and bool(row["peak_internal"])
    trow = tr[(tr["gene"] == "GT") & (tr["group"] == "Purkinje")].iloc[0]
    assert trow["status"] == "ok" and trow["n_bins"] == 5
    assert row["rho_bin"] == trow["rho_bin"] and row["q_bin"] == trow["q_bin"]
    for g in tr["group"].unique():
        sub = tr[tr["group"].to_numpy() == g]
        q_exp = st.bh_qvalues(sub["p_bin"].to_numpy(dtype=float))
        np.testing.assert_allclose(sub["q_bin"].to_numpy(dtype=float), q_exp, equal_nan=True)
    cfg2 = _common.cfg()
    cfg2.stats_cfg["trend"]["fdr_family"] = "species"
    try:
        s02_stats.build_tables(bt, None, UNIVERSE, cfg.species("mouse"), cfg2)
        raise AssertionError("fdr_family 非 species_group 未被拦截")
    except NotImplementedError as e:
        assert "species_group" in str(e)


def test_write_tables_roundtrip_determinism():
    """write_tables → 临时目录 7 文件;读回列一致、数值 %.6g 容差内;GM 的 'NA' 串保留;
    build_tables 两跑逐帧一致(R8)。"""
    cfg, L, bt = _setup()
    spcfg = cfg.species("mouse")
    tabs = s02_stats.build_tables(bt, None, UNIVERSE, spcfg, cfg)
    tabs2 = s02_stats.build_tables(bt, None, UNIVERSE, spcfg, cfg)
    for name in s02_stats.TABLE_NAMES:
        pd.testing.assert_frame_equal(tabs[name], tabs2[name])
    with tempfile.TemporaryDirectory() as td:
        paths = s02_stats.write_tables(tabs, td)
        assert len(paths) == 7 and all(os.path.exists(p) for p in paths)
        back = pd.read_csv(os.path.join(td, "trend_by_group.tsv"), sep="\t",
                           keep_default_na=False, na_values=[""])   # 浮点列: 空串→NaN
        assert list(back.columns) == list(tabs["trend_by_group"].columns)
        a = tabs["trend_by_group"]["rho_bin"].to_numpy(dtype=float)
        b = back["rho_bin"].to_numpy(dtype=float)
        np.testing.assert_allclose(b[~np.isnan(a)], a[~np.isnan(b)], rtol=1e-5)
        em = pd.read_csv(os.path.join(td, "emergence.tsv"), sep="\t",
                         keep_default_na=False).set_index("gene")   # 原样字符串: "NA" 保留
        assert em.loc["GT", "emergence_bin"] == L[3]
        assert em.loc["GU", "emergence_bin"] == "NA"            # 受控 NA 串可 roundtrip
        tp = pd.read_csv(os.path.join(td, "temporal_patterns.tsv"), sep="\t",
                         keep_default_na=False).set_index(["gene", "group"])
        assert tp.loc[("GT", "Purkinje"), "temporal_pattern"] == "transient_peak"
        assert tp.loc[("GM", "Purkinje"), "temporal_pattern"] == "NA"


def test_mask_lowexpr_labels():
    """P7 裁决①(2026-11-05): not_expressed 基因的窗口/峰伪标签 → NA(status=ok 可识别);
    同形状正常表达基因标签不受影响;引擎函数幂等。

    夹具: L[1..8](embryonic×5 + perinatal L6/L7 + postnatal L8), GLO 内部峰在 L[6]、
    period_ratio≈3.2≥2 → window_high; 但全部 cp≤0.03<min_expr 0.1 → not_expressed。
    GT = 同形状 ×1000(pvf≥0.5 → transient_peak), 高表达不受掩膜影响。
    (低 cp 下 log1p 压缩使 pvf 路径不可触发——真实草稿的伪标签走的正是 period_ratio 路径。)
    """
    cfg, L, _bt = _setup()
    from atx import stats_pattern as atxpattern
    shape = {1: 2.0, 2: 2.0, 3: 2.0, 4: 2.0, 5: 2.0, 6: 30.0, 7: 20.0, 8: 2.0}
    rows = []
    for grp in ("Purkinje", "GranuleLineage", "Astro", "Immune", "ALL"):  # ≥4 主组(min_groups)
        for i in range(1, 9):
            gt = shape[i] if grp in ("Purkinje", "ALL") else 2.0
            rows.append((grp, L[i], True, 40, {"GT": gt, "GLO": gt / 1000.0}))
    bt2 = synth.mk_bin_table(cfg, rows, ["GT", "GLO"], value_kind="cp")
    tabs = s02_stats.build_tables(bt2, None, ["GT", "GLO"], cfg.species("mouse"), cfg)
    sg = tabs["specificity_global"].set_index("gene")
    assert sg.loc["GLO", "specific_call"] == "not_expressed"      # cp≤0.03 < min_expr 0.1
    assert sg.loc["GT", "specific_call"] != "not_expressed"
    tp = tabs["temporal_patterns"].set_index(["gene", "group"])
    assert tp.loc[("GT", "Purkinje"), "temporal_pattern"] == "transient_peak"   # 高表达同形状保留
    r = tp.loc[("GLO", "Purkinje")]
    assert r["temporal_pattern"] == "NA" and r["status"] == "ok"  # 掩膜行 = ok+NA
    n = atxpattern.mask_lowexpr_labels(tabs["temporal_patterns"], tabs["specificity_global"])
    assert n == 0                                                 # 已掩膜 → 重跑 0 行(幂等)
