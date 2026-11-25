"""stats_specificity 测试(spec §3.13):3 组可行 → tau/fraction NA;8 组 tau=手算;
fraction_8group 严格版;specific_call 边界与 NA 语义;by-period 与 bin 等权;
tau_by_bin+emergence 手算(含左截断);确定性(R8)。阈值全部从 config 读, 不硬编码。"""
from __future__ import annotations

import numpy as np
import pandas as pd

import _common
import synth
from atx import stats_specificity as ss


def _par(cfg):
    sp = cfg.stats_cfg["specificity"]
    tb = cfg.stats_cfg["tau_by_bin"]
    return dict(pc=float(sp["pseudocount"]), min_groups=int(tb["min_groups"]),
                tau_thr=float(sp["tau_threshold"]), frac_thr=float(sp["fraction_threshold"]),
                min_expr=float(sp["min_expr"]), em_thr=float(tb["fraction_threshold"]))


def test_three_feasible_groups_all_na():
    """k=3 < min_groups → tau/fraction/specific_call 全 NA;UBC/ALL 不入分母。"""
    cfg = _common.cfg()
    main = cfg.main_groups()
    genes = ["GA"]
    L = cfg.species("mouse").labels
    rows = [(g, L[2], True, 40, {"GA": v}) for g, v in zip(main[:3], (10.0, 9.0, 8.0))]
    rows += [(g, L[2], False, 5, {"GA": np.nan}) for g in main[3:]]
    rows += [(g, L[2], True, 40, {"GA": 100.0}) for g in ("UBC", "ALL")]   # 必须被排除
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="cp")
    sc = ss.scores(ss.aggregate(bt, genes, cfg), genes, cfg)
    r = sc.iloc[0]
    assert r["status"] == "ok" and r["n_groups_used"] == 3
    assert r["groups_used"] == "|".join(main[:3])        # UBC/ALL 未混入
    assert np.isnan(r["tau"]) and np.isnan(r["max_fraction"])
    assert all(np.isnan(r[f"fraction_{g}"]) for g in main)
    assert all(np.isnan(r[f"fraction_8group_{g}"]) for g in main)
    assert r["specific_call"] == "NA"                    # NA 而非 False(§3-14 语义)
    assert r["max_group"] == main[0]                     # 信息列仍给可行组内 argmax
    assert r[f"cp10k_{main[0]}"] == 10.0                 # cp/pct 照常报告


def test_tau_handcalc_8groups():
    """8 组全可行: tau = Σ(1−x'/max')/(k−1) 手算(含 pseudocount);specific_call=specific。"""
    cfg = _common.cfg()
    P = _par(cfg)
    main = cfg.main_groups()
    genes = ["GA"]
    L = cfg.species("mouse").labels
    vals = {main[0]: 8.0, **{g: 1.0 for g in main[1:]}}
    rows = [(g, L[2], True, 40, {"GA": v}) for g, v in vals.items()]
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="cp")
    r = ss.scores(ss.aggregate(bt, genes, cfg), genes, cfg).iloc[0]
    pc = P["pc"]
    tau_exp = 7.0 / (8.0 + pc)                           # 7*(1−(1+pc)/(8+pc))/7
    frac_exp = (8.0 + pc) / (8.0 + pc + 7 * (1.0 + pc))  # = (8+pc)/(15+8pc)
    assert np.isclose(r["tau"], tau_exp)
    assert np.isclose(r[f"fraction_{main[0]}"], frac_exp)
    assert np.isclose(r[f"fraction_8group_{main[0]}"], frac_exp)   # 8 组全可行 → 严格版=同值
    assert r["max_group"] == main[0] and np.isclose(r["max_fraction"], frac_exp)
    assert tau_exp >= P["tau_thr"] and frac_exp >= P["frac_thr"]   # 场景有效性
    assert r["specific_call"] == "specific"


def test_fraction_8group_strict_requires_all_main():
    """k=7: fraction(可行分母)照算, fraction_8group 全 NA(仅 8 组全可行才非 NA)。"""
    cfg = _common.cfg()
    P = _par(cfg)
    main = cfg.main_groups()
    genes = ["GA"]
    L = cfg.species("mouse").labels
    rows = [(main[0], L[2], False, 5, {"GA": np.nan})]
    rows += [(main[1], L[2], True, 40, {"GA": 8.0})]
    rows += [(g, L[2], True, 40, {"GA": 1.0}) for g in main[2:]]
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="cp")
    r = ss.scores(ss.aggregate(bt, genes, cfg), genes, cfg).iloc[0]
    pc = P["pc"]
    assert r["n_groups_used"] == 7
    frac_exp = (8.0 + pc) / (8.0 + pc + 6 * (1.0 + pc))
    assert np.isclose(r[f"fraction_{main[1]}"], frac_exp)         # 分母=7 可行组
    assert all(np.isnan(r[f"fraction_8group_{g}"]) for g in main)
    assert r["max_group"] == main[1]


def test_specific_call_boundaries():
    cfg = _common.cfg()
    P = _par(cfg)
    main = cfg.main_groups()
    genes = ["GA"]
    L = cfg.species("mouse").labels

    def call_of(cp_by_group, cfg_use):
        rows = [(g, L[2], True, 40, {"GA": v}) for g, v in cp_by_group.items()]
        bt = synth.mk_bin_table(cfg_use, rows, genes, value_kind="cp")
        return ss.scores(ss.aggregate(bt, genes, cfg_use), genes, cfg_use).iloc[0]

    uniform = {g: 5.0 for g in main}
    r = call_of(uniform, _common.cfg())
    assert r["tau"] == 0.0 and r["specific_call"] == "broad"      # tau=0 < tau_thr

    tiny = {g: 0.05 for g in main}
    assert 0.05 < P["min_expr"]                                    # 场景有效性
    r = call_of(tiny, _common.cfg())
    assert r["specific_call"] == "not_expressed"                   # 最强组 CP10K < min_expr

    peaked = {main[0]: 8.0, **{g: 1.0 for g in main[1:]}}
    r = call_of(peaked, _common.cfg())
    t0 = float(r["tau"])
    cfg_eq = _common.cfg()
    cfg_eq.stats_cfg["specificity"]["tau_threshold"] = t0          # >= 边界含等号
    assert call_of(peaked, cfg_eq)["specific_call"] == "specific"
    cfg_gt = _common.cfg()
    cfg_gt.stats_cfg["specificity"]["tau_threshold"] = t0 + 1e-9
    assert call_of(peaked, cfg_gt)["specific_call"] == "broad"


def test_aggregate_period_scope_and_bin_equal_weight():
    """bin 等权: 30 核 bin 与 300 核 bin 同权;period 限定只聚合该期。"""
    cfg = _common.cfg()
    spcfg = cfg.species("mouse")
    L = spcfg.labels
    p_emb, p_late = spcfg.period_of[L[2]], spcfg.period_of[L[7]]   # E12.5(胚胎) vs P0(围产)
    assert p_emb == spcfg.period_of[L[3]] and p_emb != p_late      # 场景有效性
    genes = ["GA"]
    rows = [("Purkinje", L[2], True, 30, {"GA": 2.0}),
            ("Purkinje", L[3], True, 300, {"GA": 4.0}),
            ("Purkinje", L[7], True, 40, {"GA": 10.0})]
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="cp")
    agg = ss.aggregate(bt, genes, cfg)
    r = agg[agg["group"] == "Purkinje"].iloc[0]
    assert np.isclose(r["cp10k"], 16.0 / 3.0)                      # (2+4+10)/3, 与细胞数无关
    assert r["n_bins"] == 3 and r["n_cells_total"] == 370          # 30+300+40
    agg_p = ss.aggregate(bt, genes, cfg, period=p_emb)
    assert np.isclose(agg_p[agg_p["group"] == "Purkinje"].iloc[0]["cp10k"], 3.0)
    assert len(ss.aggregate(bt, ["GC_absent"], cfg)) == 0          # 缺失基因不产行


def _rows_for(cfg):
    """6 bin 场景: b0/b1 仅 3 组可行(左截断) → b2 5 组 → b3–b5 8 组;
    GA=晚升(P4 型), GB=均匀(ties), GC=缺失, GD=b2 即特异(左删失)。"""
    spcfg = cfg.species("mouse")
    L = spcfg.labels
    main = list(cfg.main_groups())
    b = [L[2], L[3], L[4], L[5], L[6], L[7]]
    feas = [[main[0], main[2], main[6]],
            [main[0], main[2], main[6]],
            [main[0], main[1], main[2], main[5], main[6]],
            list(main), list(main), list(main)]
    ga = [{g: 1.0 for g in main} for _ in range(6)]
    ga[4][main[0]] = 11.0
    ga[5][main[0]] = 20.0
    gd = [{g: 1.0 for g in main} for _ in range(6)]
    for i in (2, 3, 4, 5):
        gd[i][main[0]] = 11.0
    rows = []
    for i in range(6):
        for g in main:
            f = g in feas[i]
            rows.append((g, b[i], f, 40 if f else 5,
                         {"GA": ga[i][g], "GB": 1.0, "GD": gd[i][g]}))
    return b, main, rows


def test_tau_by_bin_emergence_handcalc():
    cfg = _common.cfg()
    P = _par(cfg)
    b, main, rows = _rows_for(cfg)
    genes = ["GA", "GB", "GC", "GD"]
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="cp")
    tb = ss.tau_by_bin(bt, genes, cfg)
    pc, em_thr = P["pc"], P["em_thr"]

    ga = tb[tb["gene"] == "GA"].set_index("bin")
    assert ga.loc[b[0], "n_groups_used"] == 3 and np.isnan(ga.loc[b[0], "tau"])
    assert np.isnan(ga.loc[b[0], f"frac_{main[0]}"])               # k<min_groups → frac NA
    assert ga.loc[b[2], "n_groups_used"] == 5
    assert np.isclose(ga.loc[b[2], "tau"], 0.0)                    # 5 组全等 → tau=0
    assert np.isclose(ga.loc[b[2], f"frac_{main[0]}"], 1.0 / 5.0)  # (1+pc)/(5(1+pc))=1/5
    frac_b4 = (11.0 + pc) / (11.0 + pc + 7 * (1.0 + pc))
    assert np.isclose(ga.loc[b[4], f"frac_{main[0]}"], frac_b4)
    assert frac_b4 >= em_thr > 1.0 / 5.0 and em_thr > 1.0 / 8.0    # 场景有效性
    assert set(ga.loc[b[3]].index.isin([f"frac_{g}" for g in main]))

    sc = ss.scores(ss.aggregate(bt, genes, cfg), genes, cfg).set_index("gene")
    tau_ga = 1.0 - (1.0 + pc) / (35.0 / 6.0 + pc)                  # main0 均值 35/6, 余 1
    frac_ga = (35.0 / 6.0 + pc) / (35.0 / 6.0 + pc + 7 * (1.0 + pc))
    assert sc.loc["GA", "max_group"] == main[0]
    assert np.isclose(sc.loc["GA", "tau"], tau_ga)
    assert frac_ga < P["frac_thr"] and tau_ga >= P["tau_thr"]      # tau 过而 fraction 不过
    assert sc.loc["GA", "specific_call"] == "broad"                # AND 语义
    assert sc.loc["GB", "max_group"] == main[0]                    # ties 取主组序首个
    assert sc.loc["GB", "tau"] == 0.0
    assert sc.loc["GC", "status"] == "missing" and sc.loc["GC", "specific_call"] == "NA"
    assert sc.loc["GC", "max_group"] == "NA" and np.isnan(sc.loc["GC", "tau"])
    tau_gd = 1.0 - (1.0 + pc) / (46.0 / 6.0 + pc)                  # GD main0 均值 46/6
    frac_gd = (46.0 / 6.0 + pc) / (46.0 / 6.0 + pc + 7 * (1.0 + pc))
    assert np.isclose(sc.loc["GD", "tau"], tau_gd)
    if tau_gd >= P["tau_thr"] and frac_gd >= P["frac_thr"]:
        assert sc.loc["GD", "specific_call"] == "specific"

    em = ss.emergence(tb, sc["max_group"], cfg).set_index("gene")
    assert em.loc["GA", "first_computable_bin"] == b[2]            # 首个 k≥min_groups 的 bin
    assert em.loc["GA", "emergence_bin"] == b[4]                   # frac 首达 ≥0.6
    assert np.isclose(em.loc["GA", "frac_at_emergence"], frac_b4)
    assert not bool(em.loc["GA", "left_censored"])
    assert em.loc["GB", "emergence_bin"] == "NA" and em.loc["GB", "first_computable_bin"] == b[2]
    assert em.loc["GC", "emergence_bin"] == "NA" and em.loc["GC", "first_computable_bin"] == "NA"
    assert em.loc["GD", "emergence_bin"] == b[2]                   # 首个可算 bin 即达阈
    assert bool(em.loc["GD", "left_censored"])                     # → 左删失("不晚于 b2 已特异")


def test_scores_deterministic():
    cfg = _common.cfg()
    b, main, rows = _rows_for(cfg)
    genes = ["GA", "GB", "GC", "GD"]
    bt = synth.mk_bin_table(cfg, rows, genes, value_kind="cp")
    agg = ss.aggregate(bt, genes, cfg)
    pd.testing.assert_frame_equal(ss.scores(agg, genes, cfg), ss.scores(agg, genes, cfg))
    pd.testing.assert_frame_equal(ss.tau_by_bin(bt, genes, cfg), ss.tau_by_bin(bt, genes, cfg))
    pd.testing.assert_frame_equal(agg, ss.aggregate(bt, genes, cfg))   # R8
