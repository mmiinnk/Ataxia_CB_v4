"""groups.py 测试:合成 200 核指派 / RL 归属裁定 / NA 哨兵 / 互斥断言 / expand_long / 锚点函数。"""
from __future__ import annotations

import numpy as np
import pandas as pd

import _common
import synth
from atx import groups as atxgroups


def _derive_assign():
    c = _common.cfg()
    spcfg = c.species("mouse")
    obs_raw, X, var = synth.build()
    d = atxgroups.derive_obs(obs_raw, spcfg, c)
    a, masks = atxgroups.assign(d, c)
    return c, spcfg, obs_raw, d, a, masks


def test_derive_obs_basics():
    c, spcfg, obs_raw, d, a, _ = _derive_assign()
    assert len(d) == 200
    # NA 安全化:合成 obs 的 subtype/dev_state/cell_type 全 NA → 哨兵
    assert (d["subtype"] == atxgroups.NA).all()
    assert not d["sub_is_progenitor"].any() and d["sub_na"].all()
    # bin 标签与 bin_order(鼠 E10.5..E14.5 → 1..5)
    assert set(d["bin"]) == {"E10.5", "E11.5", "E12.5", "E13.5", "E14.5"}
    assert d.loc[d["bin"] == "E10.5", "bin_order"].eq(1).all()
    assert d.loc[d["bin"] == "E14.5", "bin_order"].eq(5).all()
    # pooled 旗标(E10.5/E11.5)
    assert d.loc[d["bin"].isin(["E10.5", "E11.5"]), "pooled_bin"].all()
    assert not d.loc[d["bin"].isin(["E12.5", "E13.5", "E14.5"]), "pooled_bin"].any()
    # 疾病旗标:config 只在人侧声明 disease_column(鼠无 → 不派生)
    assert spcfg.disease_col is None and "disease_flag" not in d.columns


def test_derive_obs_human_disease_flag():
    """人侧:disease != normal → disease_flag=True(§11 未决 9;donor cube 保留该旗标)。"""
    c = _common.cfg()
    spcfg = c.species("human")
    obs = pd.DataFrame({
        "author_cell_type": ["Purkinje", "Purkinje"],
        "subtype": [None, None], "dev_state": [None, None], "cell_type": [None, None],
        "development_stage": ["44-year-old stage", "newborn stage (0-28 days)"],
        "donor_id": ["h1", "h2"],
        "disease": ["heart failure", "congenital heart disease"],
    })
    d = atxgroups.derive_obs(obs, spcfg, c)
    assert d["bin"].tolist() == ["44y", "NB"] and d["bin_order"].tolist() == [14, 8]
    assert d["disease_flag"].tolist() == [True, True]
    obs.loc[:, "disease"] = "normal"
    d2 = atxgroups.derive_obs(obs, spcfg, c)
    assert not d2["disease_flag"].any()


def test_assign_counts_and_exclusion():
    c, spcfg, obs_raw, d, a, masks = _derive_assign()
    counts = atxgroups.group_counts(a)
    assert counts["Purkinje"] == 103
    assert counts["GranuleLineage"] == 76   # PLAN 手算: 15+15+12+12+12+10
    assert counts["Immune"] == 5
    assert counts["ungrouped"] == 16          # meningeal 只进 ALL
    assert counts["total"] == 200
    # 互斥:任意两组掩码交集为 0(assign 内已断言, 此处再显式核一遍)
    names = list(masks)
    for i, g1 in enumerate(names):
        for g2 in names[i + 1:]:
            assert not (masks[g1] & masks[g2]).any()
    # 细胞级标签与 PLAN 一致
    exp = np.array([synth.ACT_TO_GROUP[x] for x in
                    np.repeat([p[0] for p in synth.PLAN], [p[3] for p in synth.PLAN])])
    assert (a["cell_group"].to_numpy() == exp).all()


def test_expand_long():
    c, spcfg, obs_raw, d, a, _ = _derive_assign()
    lg = atxgroups.expand_long(a, c)
    # 184 个有组细胞各 1 行本组 + 200 行 ALL = 384
    assert len(lg) == 184 + 200
    assert (lg["cell_group"] == "ALL").sum() == 200
    ung = set(a.loc[a["cell_group"] == "ungrouped", "cell_pos"])
    own = lg[lg["cell_group"] != "ALL"]
    assert not (set(own["cell_pos"]) & ung), "ungrouped 只进 ALL(spec §4.1)"


def test_rl_precedence_ruling():
    """astroglia×progenitor_RL(_early) 归 GranuleLineage(2026-11-04 实测裁定, groups.py docstring)。"""
    c = _common.cfg()
    spcfg = c.species("mouse")
    rows = [
        ("astroglia", "progenitor_RL", "progenitor", "GranuleLineage"),
        ("astroglia", "progenitor_RL_early", "progenitor", "GranuleLineage"),
        ("astroglia", "progenitor_gliogenic", "progenitor", "Progenitor"),
        ("astroglia", "progenitor_VZ_early", "progenitor", "Progenitor"),
        ("astroglia", None, "progenitor", "Progenitor"),        # NA subtype 第三项
        ("astroglia", "astro_Bergmann", "astrocyte", "Astro"),
        ("astroglia", "glioblast_PWM", "glioblast", "Astro"),
        ("VZ_neuroblast", None, "VZ_neuroblast_1", "Progenitor"),
        ("NTZ_mixed", None, "progenitor", "Progenitor"),
        ("GC", None, "GC_diff_1", "GranuleLineage"),
        ("glutamatergic_uncertain", None, "GC_defined", "GranuleLineage"),  # 人独有标签
        ("oligo", "OPC", "oligo_progenitor", "OligoLineage"),
        ("interneuron", "interneuron_ML1", "interneuron_diff", "Interneuron"),
        ("glut_DN", "glut_DN_posterior", "glut_DN_defined", "DN"),
        ("UBC", "UBC_Trpc3", "UBC_diff", "UBC"),
        ("immune", "microglia", "immune", "Immune"),
        ("Purkinje", "Purkinje_defined_EB", "Purkinje_defined", "Purkinje"),
        ("ependymal", None, "ependymal", "ungrouped"),          # 鼠独有 → 只进 ALL
    ]
    obs = pd.DataFrame({
        "author_cell_type": [r[0] for r in rows],
        "subtype": [r[1] for r in rows],
        "dev_state": [r[2] for r in rows],
        "cell_type": [None] * len(rows),
        "author_stage": ["E14.5"] * len(rows),
        "donor_id": ["d1"] * len(rows),
    })
    d = atxgroups.derive_obs(obs, spcfg, c)
    a, _ = atxgroups.assign(d, c)   # 互斥断言在 assign 内;RL 裁定失败会当场 AssertionError
    got = a["cell_group"].tolist()
    exp = [r[3] for r in rows]
    assert got == exp, f"RL/Progenitor 归属裁定失败:\n got={got}\n exp={exp}"


def test_sensitivity_expr_rewrites():
    c = _common.cfg()
    base = atxgroups.effective_group_exprs(c)
    assert "glutamatergic_uncertain" in base["GranuleLineage"]
    c.sensitivity["drop_glutamatergic_uncertain"] = True
    e = atxgroups.effective_group_exprs(c)
    assert "glutamatergic_uncertain" not in e["GranuleLineage"]
    assert "'GC'" in e["GranuleLineage"]            # 其余项保留
    c.sensitivity["drop_glutamatergic_uncertain"] = False
    c.sensitivity["split_granule_progenitor"] = True
    e = atxgroups.effective_group_exprs(c)
    assert "progenitor_RL" not in e["GranuleLineage"]
    assert e["GranuleProgenitor"] == "subtype in ['progenitor_RL','progenitor_RL_early']"
    c.sensitivity["split_granule_progenitor"] = False


def test_derive_obs_non_range_index():
    """回归(2026-11-04 s00 实测):真实 h5ad 的 obs 索引是细胞条码字符串。
    derive_obs 必须按**位置**对齐(d 用 RangeIndex), 不得被 pandas 索引重对齐静默填 NA。"""
    c = _common.cfg()
    spcfg = c.species("mouse")
    obs_raw, X, var = synth.build()
    obs_bc = obs_raw.copy()
    obs_bc.index = [f"cell_barcode_{i}" for i in range(len(obs_bc))]
    d = atxgroups.derive_obs(obs_bc, spcfg, c)
    assert list(d.index) == list(range(200))
    assert d["bin"].notna().all() and d["donor"].notna().all()
    assert not (d["subtype"] == atxgroups.NA).all() or (obs_raw["subtype"].isna()).all()
    # 与 RangeIndex 版本逐列一致(位置对齐等价)
    d_ref = atxgroups.derive_obs(obs_raw, spcfg, c)
    for col in d_ref.columns:
        assert (d[col].to_numpy() == d_ref[col].to_numpy()).all(), col


def test_check_anchors_detects_mismatch():
    """合成数据计数必然 ≠ 真实锚点 → check_anchors 必须报出(函数本身有效性测试)。"""
    c, spcfg, obs_raw, d, a, _ = _derive_assign()
    bad = atxgroups.check_anchors(a, "mouse", c)
    assert bad, "合成数据未被锚点核对捕获 = check_anchors 失效"
    assert any("total" in x for x in bad)
