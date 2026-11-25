"""可靠性修复回归测试(2026-11-25 代码审查后):

#1 align.fidelity NA 政策 — 缺失 max_group 记 NA 而非伪造 0/1(硬约束 4);
#2 s07.peak_ranks 合并 AD bin 秩 — 合并口径下 AD 峰不再静默丢失, 秩=constituent 质心 14.5;
#4 sanity.marker_table 描述组(UBC)marker 校验 — 不再永久 ABSENT, 且主组 marker 口径不变。
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

import _common
import s07_category_enrichment as s07
from atx import align as atxalign
from atx import sanity as atxsanity


def _prof(states, val):
    return pd.DataFrame([[val] * len(states)], index=["G1"], columns=states, dtype=float)


def test_fidelity_missing_max_group_is_na_not_zero():
    """鼠缺失基因(max_group NA)→ max_group_match/fidelity 均 NA, 不得伪造 0; 真实匹配仍算 1。"""
    cfg = _common.cfg()
    states = ["immune"]
    trend_h = pd.DataFrame({"gene": ["G1", "G2"], "group": ["Immune", "Immune"],
                            "rho_bin": [0.8, 0.8], "status": ["ok", "ok"]})
    trend_m = pd.DataFrame({"gene": ["G1", "G2"], "group": ["Immune", "Immune"],
                            "rho_bin": [np.nan, 0.8], "status": ["missing", "ok"]})
    prof_h = {"ALL": _prof(states, 1.0), "Immune": _prof(states, 1.0)}
    prof_m = {"ALL": _prof(states, np.nan), "Immune": _prof(states, np.nan)}
    # G1 鼠缺失; G2 两物种均 Immune(真实匹配对照)
    for na_form in ("NA", np.nan):   # 受控串 与 被解析成的 NaN 两种形态都应 → NA
        spec_h = pd.DataFrame({"gene": ["G1", "G2"], "max_group": ["Immune", "Immune"]})
        spec_m = pd.DataFrame({"gene": ["G1", "G2"], "max_group": [na_form, "Immune"]})
        fid = atxalign.fidelity(trend_h, trend_m, spec_h, spec_m, prof_h, prof_m, cfg)
        g1 = fid[fid.gene == "G1"]
        assert g1.max_group_match.isna().all(), f"na_form={na_form}: 缺失 max_group 被伪造 match"
        assert g1.fidelity.isna().all() and (g1.n_metrics == 0).all()
        assert (g1.evidence_grade == "none").all()
        g2 = fid[(fid.gene == "G2") & (fid.group == "Immune")].iloc[0]
        assert g2.max_group_match == 1.0, "真实一致 max_group 应仍算 1"


def test_peak_ranks_collapsed_ad_bin_has_centroid_rank():
    """合并口径 AD 峰 → human_only 秩 14.5(constituent 13/14/15/16 质心), peak_source=human。"""
    cfg = _common.cfg()
    assert cfg.sens("collapse_adult_bins", False), "本测试以人侧合并主口径为前提"
    ma = pd.DataFrame([
        {"species": "human", "gene": "G1", "status": "ok", "g_star": "Immune", "peak_bin": "AD"},
        {"species": "mouse", "gene": "G1", "status": "ok", "g_star": "Immune", "peak_bin": "9wk"},
        {"species": "human", "gene": "G2", "status": "ok", "g_star": "Immune", "peak_bin": "CS22"},
        {"species": "mouse", "gene": "G2", "status": "ok", "g_star": "Immune", "peak_bin": "E13.5"},
    ])
    pk = s07.peak_ranks(cfg, ma).set_index("gene")
    assert abs(pk.loc["G1", "peak_rank_human_only"] - 14.5) < 1e-9, "AD 峰秩丢失或非质心"
    assert pk.loc["G1", "peak_source"] == "human", "非删失组 AD 峰不应回退鼠峰"
    assert pk.loc["G2", "peak_rank_human_only"] == 3.0, "非合并 bin(CS22)秩应保持原始 3"


def test_marker_table_validates_descriptive_ubc_group():
    """UBC marker(EOMES)不再永久 ABSENT; 主组 marker(CALB1)argmax 仍只在主组内判。"""
    cg = pd.Series(np.array(["Purkinje"] * 40 + ["UBC"] * 30 + ["GranuleLineage"] * 30))
    X = np.zeros((100, 2))
    X[40:70, 0] = 10.0   # EOMES 高表达于 UBC
    X[:40, 1] = 10.0     # CALB1 高表达于 Purkinje
    X_cp = sp.csr_matrix(X)
    gene_cols = {"EOMES": 0, "CALB1": 1}
    markers = pd.DataFrame({"gene_symbol": ["EOMES", "CALB1", "NOGENE"],
                            "cell_group": ["UBC", "Purkinje", "UBC"],
                            "species_valid": ["both", "both", "both"]})
    mt = atxsanity.marker_table(X_cp, cg, gene_cols, markers, "human",
                                ["Purkinje", "GranuleLineage"], ["UBC"]).set_index("marker")
    assert mt.loc["EOMES", "verdict"] == "PASS" and mt.loc["EOMES", "argmax_group"] == "UBC"
    assert mt.loc["CALB1", "verdict"] == "PASS" and mt.loc["CALB1", "argmax_group"] == "Purkinje"
    assert mt.loc["NOGENE", "verdict"] == "ABSENT"
    assert mt.loc["NOGENE", "absent_reason"] == "gene_not_in_var"