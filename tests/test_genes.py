"""genes.py 测试(spec §3.13):匹配级联(Ensembl 主键 → 符号 → 别名)、ambiguous 判定、
ensg_only(Castonguay 口径, §11 未决 5:不回退符号)、hit_rate;真实 panel 装载口径。"""
from __future__ import annotations

import numpy as np
import pandas as pd

import _common
from atx import genes as atxgenes

STAT = atxgenes.STAT_COLS


def _panel_row(**kw):
    base = {c: pd.NA for c in STAT}
    base.update({"tier": "core", "in_stat_panel": True, "is_sanity": False,
                 "onset_typical_y": np.nan})
    base.update(kw)
    return base


def _fake_panel():
    return pd.DataFrame([
        _panel_row(gene_key="A", hgnc_symbol="A", ensg_id="ENSG0001"),
        _panel_row(gene_key="B", hgnc_symbol="BHS", ensg_id="ENSG9999"),   # ENSG 缺 → 符号回退
        _panel_row(gene_key="C", hgnc_symbol="CSYM", ensg_id=pd.NA),        # 符号双命中 → ambiguous
        _panel_row(gene_key="D", hgnc_symbol="D", ensg_id="ENSG8888"),      # → 别名回退
        _panel_row(gene_key="E", hgnc_symbol="EHS", ensg_id="ENSG7777"),    # 全缺 → missing
    ])


def _fake_var():
    return pd.DataFrame({"feature_name": ["A", "bhs", "CSYM", "CSYM", "DOLD"]},
                        index=["ENSG0001", "ENSG0002", "ENSG0003", "ENSG0004", "ENSG0005"])


def test_map_cascade_human():
    aliases = pd.DataFrame([{"alias": "DOLD", "hgnc_symbol": "D", "ensg_id": pd.NA,
                             "mgi_symbol": pd.NA, "ensmusg_id": pd.NA}])
    gm = atxgenes.map_to_var(_fake_panel(), _fake_var(), species="human", dataset="human",
                             aliases=aliases)
    st = gm.set_index("gene_key")
    assert st.loc["A", "status_human"] == "matched" and st.loc["A", "match_method_human"] == "ensembl"
    assert int(st.loc["A", "row_human"]) == 0
    assert st.loc["B", "match_method_human"] == "symbol" and int(st.loc["B", "row_human"]) == 1
    assert st.loc["C", "status_human"] == "ambiguous" and pd.isna(st.loc["C", "row_human"])
    assert st.loc["D", "match_method_human"] == "alias" and int(st.loc["D", "row_human"]) == 4
    assert st.loc["E", "status_human"] == "missing"
    hr = atxgenes.hit_rate(gm, "human")
    assert abs(hr - 3 / 5) < 1e-9, f"hit_rate={hr}, 期望 0.6(A/B/D)"


def test_ensg_only_no_symbol_fallback():
    """Castonguay 口径(§11 未决 5):ENSG 缺失不得静默回退符号。"""
    aliases = pd.DataFrame([{"alias": "DOLD", "hgnc_symbol": "D", "ensg_id": pd.NA,
                             "mgi_symbol": pd.NA, "ensmusg_id": pd.NA}])
    gm = atxgenes.map_to_var(_fake_panel(), _fake_var(), species="human", dataset="castonguay",
                             aliases=aliases, ensg_only=True)
    st = gm.set_index("gene_key")
    assert st.loc["A", "status_castonguay"] == "matched"
    for k in ("B", "C", "D", "E"):
        assert st.loc[k, "status_castonguay"] == "missing", f"{k} 在 ensg_only 下不得回退"
    assert abs(atxgenes.hit_rate(gm, "castonguay") - 0.2) < 1e-9


def test_map_mouse_alias():
    """鼠别名 Mre11a(实测:Sepp mouse var 以旧名存在;alias 表登记)。"""
    panel = pd.DataFrame([_panel_row(gene_key="MRE11", hgnc_symbol="MRE11",
                                     mgi_symbol="Mre11", ensmusg_id="ENSMUSG00000031928")])
    var = pd.DataFrame({"feature_name": ["Trp53", "Mre11a"]}, index=["ENSMUSG_X", "ENSMUSG_Y"])
    # ① 主键命中(var 有该 ENSMUSG)
    v2 = pd.DataFrame({"feature_name": ["Mre11a"]}, index=["ENSMUSG00000031928"])
    gm = atxgenes.map_to_var(panel, v2, species="mouse", dataset="mouse")
    assert gm.loc[0, "match_method_mouse"] == "ensembl"
    # ② 主键缺 → 别名命中(大小写原样, 鼠 sym 不大写)
    aliases = pd.DataFrame([{"alias": "Mre11a", "hgnc_symbol": "MRE11", "ensg_id": pd.NA,
                             "mgi_symbol": "Mre11a", "ensmusg_id": "ENSMUSG00000031928"}])
    gm = atxgenes.map_to_var(panel, var, species="mouse", dataset="mouse", aliases=aliases)
    assert gm.loc[0, "status_mouse"] == "matched"
    assert gm.loc[0, "match_method_mouse"] == "alias" and int(gm.loc[0, "row_mouse"]) == 1


def test_duplicate_var_index_rejected():
    var = pd.DataFrame({"feature_name": ["A", "A"]}, index=["ENSG0001", "ENSG0001"])
    try:
        atxgenes.map_to_var(_fake_panel(), var, species="human")
        raise AssertionError("重复 var index 未拒绝")
    except ValueError as e:
        assert "重复" in str(e)


def test_load_panel_real():
    c = _common.cfg()
    panel = atxgenes.load_panel(c)
    keys = set(panel["gene_key"])
    # 双标基因:CA8 既是 core panel 又是 sanity marker(plan §8;sanity_markers.tsv note)
    ca8 = panel[panel["gene_key"] == "CA8"]
    assert len(ca8) == 1 and bool(ca8.iloc[0]["in_stat_panel"]) and bool(ca8.iloc[0]["is_sanity"])
    # sanity-only marker(非致病基因, 仅校验):CALB1/PCP2 在 panel 但不进统计
    for k in ("CALB1", "PCP2"):
        r = panel[panel["gene_key"] == k]
        assert len(r) == 1 and r.iloc[0]["tier"] == "sanity_only" and not bool(r.iloc[0]["in_stat_panel"])
    # in_panel=false 的 MT 基因不进 panel(实测:MT-TK/MT-TI 无数据集含有)
    assert "MT-TK" not in keys and "MT-TI" not in keys
    # tier 口径:默认只有 core + sanity_only;peripheral 恒排除(CL-023)
    assert set(panel["tier"]) <= {"core", "sanity_only"}
    # 统计 panel 行数 == ataxia_genes.tsv 独立重算(core & in_panel=true)
    t = pd.read_csv(c.paths["genes_table"], sep="\t", dtype=str, keep_default_na=False)
    n_core = int(((t["tier"] == "core") & (t["in_panel"] == "true")).sum())
    n_stat = int(panel["in_stat_panel"].sum())
    assert n_stat == n_core, f"统计 panel {n_stat} != 表独立重算 {n_core}"
    assert not panel["gene_key"].duplicated().any()
    # include_extended_tier 敏感性:extended 纳入、peripheral 仍排除
    c.sensitivity["include_extended_tier"] = True
    p2 = atxgenes.load_panel(c)
    assert "extended" in set(p2["tier"]) and "peripheral" not in set(p2["tier"])
    c.sensitivity["include_extended_tier"] = False
