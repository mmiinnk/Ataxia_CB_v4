"""config.py 口径测试:版本门 / 组三层 / 物种派生 / 硬约束 7(人 9w 属胎儿期) / collapse。"""
from __future__ import annotations

import _common
from atx import config as atxconfig


def test_frozen_version_gate():
    fv = atxconfig.frozen_version(_common.CURATION_LOG)
    assert fv == "1.0", f"curation_log 冻结版本解析异常: {fv}"   # 2026-11-23 冻结 v1.0
    c = _common.cfg()
    assert c.check_version_gate() == fv
    assert c.version == fv and c.draft is False


def test_group_tiers():
    c = _common.cfg()
    assert c.main_groups() == ["Purkinje", "GranuleLineage", "DN", "Interneuron",
                               "OligoLineage", "Astro", "Progenitor", "Immune"]
    assert c.descriptive_groups() == ["UBC"]
    pg = c.plot_groups()
    assert pg == c.main_groups() + ["UBC", "ALL"] and len(pg) == 10
    # 敏感性: include_ubc_in_main → UBC 进 main, descriptive 让位, plot 仍 10 组
    c.sensitivity["include_ubc_in_main"] = True
    assert "UBC" in c.main_groups() and c.descriptive_groups() == []
    assert len(c.plot_groups()) == 10
    c.sensitivity["include_ubc_in_main"] = False


def test_species_cfg():
    c = _common.cfg()
    h, m = c.species("human"), c.species("mouse")
    assert len(h.labels) == 16 and len(m.labels) == 12
    assert h.time_col == "development_stage" and m.time_col == "author_stage"
    assert m.pooled == {"E10.5", "E11.5"}
    assert h.label_of["newborn stage (0-28 days)"] == "NB"
    assert m.label_of["adult (9 weeks)"] == "9wk"        # §11 未决 8:文件名安全标签
    assert h.bin_order["CS18"] == 1 and h.bin_order["52y"] == 16
    # 硬约束 7:人 9 wpc 属胎儿期(示意图把它放胚胎期是已知错误, 不得沿用)
    assert h.period_of["9w"] == "fetal", "人 9w 必须属 fetal(硬约束 7)"
    assert h.period_of["CS22"] == "embryonic" and m.period_of["E17.5"] == "perinatal"
    assert set(h.period_of.values()) == {"embryonic", "fetal", "infant", "adult"}
    assert h.single_donor_bins == {"CS19", "20w", "6mo", "7mo", "9mo", "42y", "44y", "46y", "52y"}


def test_collapse_adult_bins():
    c = _common.cfg()
    h = c.species("human").collapsed({"42y": "AD", "44y": "AD", "46y": "AD", "52y": "AD"})
    assert len(h.labels) == 13 and h.labels[-1] == "AD"
    assert h.bin_order["AD"] == 13 and h.period_of["AD"] == "adult"
    # 合并 bin 4 供体 → 不再是 single_donor
    assert not (h.single_donor_bins & {"AD", "42y", "44y", "46y", "52y"})
    assert h.single_donor_bins == {"CS19", "20w", "6mo", "7mo", "9mo"}
    assert len({v for v in h.label_of.values()}) == 13


def test_results_path_and_draft_marker(tmp=None):
    import copy
    import os
    import tempfile
    c = _common.cfg()
    p = c.results("intermediate", "mouse")
    assert p.endswith(os.path.join("results", f"v{c.version}", "intermediate", "mouse"))
    # 冻结 v1.0: write_draft_marker 为 no-op, 不产生 DRAFT 标记
    c.write_draft_marker()
    assert not os.path.exists(c.results("DRAFT"))
    # draft 分支(version != 1.0)仍写 DRAFT; 临时目录隔离, 不触真实 results/
    raw = copy.deepcopy(c.raw)
    raw["gene_set"]["version"] = "0.9.1-draft"
    with tempfile.TemporaryDirectory() as td:
        raw["paths"]["project_root"] = td
        d = atxconfig.Cfg(raw, c.path)
        assert d.draft is True
        d.write_draft_marker()
        assert os.path.exists(d.results("DRAFT"))


def test_trend_anchors():
    """§9(4) 趋势锚点入 config(2026-11-05 裁决): 人鼠各三条, sign 归一, 每条带来源(硬约束 10)。"""
    c = _common.cfg()
    expect = [("ITPR1", "Purkinje", "+"), ("MBP", "OligoLineage", "+"),
              ("ATOH1", "GranuleLineage", "-")]
    for sp in ("mouse", "human"):
        got = [(a["gene"], a["group"], a["sign"]) for a in c.trend_anchors(sp)]
        assert got == expect, f"{sp}: {got}"
        assert all(a["source"] for a in c.trend_anchors(sp))
    # 非法 sign → ValueError;未设物种 → []
    import copy
    raw = copy.deepcopy(c.raw)
    raw["species"]["mouse"]["trend_anchors"] = [{"gene": "X", "group": "G", "sign": "0"}]
    bad = atxconfig.Cfg(raw, c.path)
    try:
        bad.trend_anchors("mouse")
        raise AssertionError("非法 sign 未 ValueError")
    except ValueError:
        pass
    raw2 = copy.deepcopy(c.raw)
    raw2["species"]["human"].pop("trend_anchors")
    assert atxconfig.Cfg(raw2, c.path).trend_anchors("human") == []
