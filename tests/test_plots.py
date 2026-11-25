"""atx/plots + s03 合成测试(spec §3.12/§5 图产物;figure_reading_guide 约定)。

覆盖: 安全化/星号/z-score/连续可行段等纯函数手算;分页页数与全局色标;NA 斜纹
patch 数 == NA 格数;b2 轨迹计数(donor 灰点/疾病△/sem 带限连续段/单供体▲/斜纹);
dotplot 点数与轴结构;b3 星号文本数;footer(draft 声明 + notes 复用机制,R5);
双格式落盘;s03 build_figures 端到端(临时目录, 不触真实 results/)与确定性。
"""
from __future__ import annotations

import math
import os
import tempfile

import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection

import _common
import s02_stats
import s03_figures_species as s03
import synth
from atx import plots as P

# ------------------------------------------------------------------ 纯函数
def test_safe_name_star_runs():
    assert P.safe_name("adult (9 weeks)") == "adult_9_weeks"
    assert P.safe_name("GranuleLineage") == "GranuleLineage"
    cfg = _common.cfg()
    assert P.star_of_q(0.2, cfg) == ""
    assert P.star_of_q(0.04, cfg) == "*"
    assert P.star_of_q(0.005, cfg) == "**"
    assert P.star_of_q(0.0005, cfg) == "***"
    assert P.star_of_q(np.nan, cfg) == "" and P.star_of_q(None, cfg) == ""
    # 连续可行段: 缺口(bin_order 不相邻或 ok=False)必须断段, 不桥接
    assert P.contig_runs([2, 3, 4, 5, 6, 7], [True, True, False, True, True, True]) == [(0, 1), (3, 5)]
    assert P.contig_runs([1, 2, 4], [True, True, True]) == [(0, 1), (2, 2)]   # 序号跳变也断
    assert P.contig_runs([1, 2], [False, False]) == []


def test_zscore_and_order_helpers():
    M = np.array([[0.0, 2.0], [1.0, 1.0], [np.nan, np.nan], [5.0, np.nan]])
    Z = P.zscore_rows(M)
    np.testing.assert_allclose(Z[0], [-1.0, 1.0])          # mean1 sd1
    np.testing.assert_allclose(Z[1], [0.0, 0.0])           # 恒定行 → 0
    assert np.isnan(Z[2]).all()                            # 全 NaN 行保留
    assert Z[3, 0] == 0.0 and np.isnan(Z[3, 1])            # 单可行 bin → 0
    assert P.order_rho_desc([np.nan, 0.5, 0.9, 0.5], ["B", "A", "C", "D"]) == ["C", "A", "D", "B"]
    sg = pd.DataFrame({"gene": ["g1", "g2", "g3", "g4"],
                       "max_group": ["Purkinje", "Immune", "Purkinje", "NA"],
                       "tau": [0.5, 0.9, 0.8, np.nan]})
    cfg = _common.cfg()
    # 主组序(Purkinje 首)→ 组内 tau 降序 → NA 沉底
    assert P.order_main_tau(sg, cfg.main_groups()) == ["g3", "g1", "g2", "g4"]


# ------------------------------------------------------------------ 热图原语
def test_paged_heatmap_and_hatch():
    cfg = _common.cfg()
    rng = np.random.default_rng(cfg.seed)
    M = rng.normal(size=(130, 5))
    M[3, 1] = M[70, 2] = M[129, 4] = np.nan
    ps = int(cfg.raw["figures"]["heatmap_page_size"])
    pages = P.paged_heatmap(M, [f"g{i}" for i in range(130)], list("abcde"), cfg,
                            vmin=-2.5, vmax=2.5, title="t")
    assert len(pages) == math.ceil(130 / ps) == 3
    clim = [ax.images[0].get_clim() for _fig, ax, _sl in pages]
    assert clim[0] == clim[1] == clim[2] == (-2.5, 2.5)     # 色标全局(分页前定)
    n_hatch_p1 = len(pages[0][1].patches)                    # 斜纹 patch 数 == 该页 NA 格数
    assert n_hatch_p1 == int((~np.isfinite(M[0:ps])).sum()) == 1
    assert len(pages[1][1].patches) == 1 and len(pages[2][1].patches) == 1
    assert "page 2/3" in pages[1][0].axes[0].get_title()
    import matplotlib.pyplot as plt
    for fig, _a, _s in pages:
        plt.close(fig)


# ------------------------------------------------------------------ b2 轨迹
def _traj_fixture():
    """Purkinje × L[1..6]: L[3] 不可行;L[4] 单供体(sem=NA);其余双供体 sem=0.1。
    cube: L1/L2/L5/L6 各 2 donor(L5 的 d2 疾病), L4 仅 d1 → 灰点 9、疾病△ 1。"""
    cfg = _common.cfg()
    L = cfg.species("mouse").labels
    rows = []
    for i in (1, 2, 4, 5, 6):
        rows.append(("Purkinje", L[i], True, 60, {"GT": float(i)}))
    rows.append(("Purkinje", L[3], False, 5, {"GT": np.nan}))
    bt = synth.mk_bin_table(cfg, rows, ["GT"], value_kind="log")
    bt = bt.set_index(["group", "bin"])
    bt.loc[("Purkinje", L[4]), "single_donor"] = True
    bt.loc[("Purkinje", L[4]), "n_donors"] = 1
    for i in (1, 2, 5, 6):
        bt.loc[("Purkinje", L[i]), "logsem_GT"] = 0.1
    bt = bt.reset_index()
    crows = []
    for i, donors in ((1, ("d1", "d2")), (2, ("d1", "d2")), (4, ("d1",)),
                      (5, ("d1", "d2")), (6, ("d1", "d2"))):
        for d in donors:
            crows.append(("Purkinje", L[i], d, {"GT": float(i) + 0.01}))
    cube = synth.mk_cube(cfg, crows, ["GT"])
    dis = np.zeros(len(cube), dtype=bool)
    for k, (g, b, d) in enumerate((r[0:3] for r in crows)):
        if b == L[5] and d == "d2":
            dis[k] = True
    cube.obs["disease_flag"] = dis
    return cfg, L, bt, cube


def test_traj_panel_counts():
    import matplotlib.pyplot as plt
    cfg, L, bt, cube = _traj_fixture()
    spcfg = cfg.species("mouse")
    fig, ax = plt.subplots()
    info = P.traj_panel(ax, bt, cube, "GT", "Purkinje", cfg, spcfg,
                        rng=np.random.default_rng(cfg.seed))
    # ok = [T,T,F,T,T,T](bin_order 2..7)→ 均值线段 (0,1),(3,5) = 2 段
    assert info["n_seg"] == 2
    # oksem = [T,T,F,F,T,T] → 带段 (0,1),(4,5) = 2(单供体 L4 与缺口 L3 均断带)
    assert info["n_band"] == 2
    assert info["n_single"] == 1                            # L[4] ▲
    assert info["n_donor_pts"] == 9 and info["n_disease_pts"] == 1
    # x 轴 = 全 12 bin 骨架: 不可行 = L[3] + 夹具未覆盖的 7 个 label(全 NA 行)
    assert info["n_hatch"] == len(spcfg.labels) - 5 == 7
    n_fill = sum(isinstance(c, PolyCollection) for c in ax.collections)
    assert n_fill == info["n_band"]                         # fill_between 落地数一致
    plt.close(fig)
    # cube=None → 无 donor 点, 其余不变(合成测试路径)
    fig, ax = plt.subplots()
    info0 = P.traj_panel(ax, bt, None, "GT", "Purkinje", cfg, spcfg)
    assert info0["n_donor_pts"] == 0 and info0["n_band"] == 2 and info0["n_seg"] == 2
    plt.close(fig)


def test_traj_figure_pagination():
    import matplotlib.pyplot as plt
    cfg, L, bt, cube = _traj_fixture()
    spcfg = cfg.species("mouse")
    genes = ["GT"] * 7                                       # 7 panel → ceil(7/3)=3 行
    fig, infos = P.traj_figure(bt, cube, genes, "Purkinje", cfg, spcfg, ncols=3)
    assert len(infos) == 7 and all(i["n_band"] == 2 for i in infos)
    assert len(fig.axes) == 9                                # 3×3 网格(第 9 格 axis off, 无 colorbar)
    plt.close(fig)


# ------------------------------------------------------------------ c1 / b3 / footer / 落盘
def test_dotplot_structure():
    import matplotlib.pyplot as plt
    cfg, L, bt, cube = _traj_fixture()
    sg = pd.DataFrame({"gene": ["GT", "GU"], "max_group": ["Purkinje", "NA"],
                       "tau": [0.7, np.nan]})
    bt2 = bt.copy()                                          # GU 无列 → 无点(不插值)
    fig, info = P.dotplot_figure(bt2, sg, ["GT", "GU"], ["Purkinje"], cfg, species="mouse")
    assert info["n_dots"] == 1                               # 仅 GT×Purkinje 有限值
    assert info["n_genes"] == 2 and info["n_groups"] == 1
    plt.close(fig)


def test_rho_heatmap_stars():
    import matplotlib.pyplot as plt
    cfg = _common.cfg()
    tr = pd.DataFrame({"gene": ["g1", "g1", "g2", "g2"],
                       "group": ["A", "B", "A", "B"],
                       "rho_bin": [0.9, np.nan, -0.5, 0.2],
                       "q_bin": [0.0005, np.nan, 0.03, 0.4]})
    pages = P.rho_heatmap_pages(tr, ["g1", "g2"], ["A", "B"], cfg, title="t")
    assert len(pages) == 1
    fig, ax, _sl = pages[0]
    texts = [t.get_text() for t in ax.texts]
    assert texts.count("***") == 1 and texts.count("*") == 1 and len(texts) == 2
    assert len(ax.patches) == 1                              # rho NA → 1 格斜纹
    assert ax.images[0].get_clim() == (-1.0, 1.0)            # rho 固定色标 ±1
    plt.close(fig)


def test_footer_and_save_formats():
    cfg = _common.cfg()
    assert not cfg.draft                                     # v1.0 已冻结: 无 DRAFT 声明
    lines = P.footer_lines(cfg, notes=["peak axis: partially interpolated"])
    assert lines[0] == "peak axis: partially interpolated"   # R5 notes 复用机制
    assert not any("DRAFT" in s for s in lines)
    assert P.footer_lines(cfg) == lines[1:]
    # draft 分支(version != 1.0)仍带 DRAFT 声明
    import copy
    from atx import config as atxconfig
    raw = copy.deepcopy(cfg.raw)
    raw["gene_set"]["version"] = "0.9.1-draft"
    d = atxconfig.Cfg(raw, cfg.path)
    assert any("DRAFT" in s and d.version in s for s in P.footer_lines(d))
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    with tempfile.TemporaryDirectory() as td:
        paths = P.save_figures(fig, cfg, "mouse_x1_test_p1", td, notes=["n1"])
        exts = sorted(os.path.splitext(p)[1] for p in paths)
        assert exts == sorted("." + e for e in cfg.raw["figures"]["formats"])
        assert all(os.path.getsize(p) > 0 for p in paths)


# ------------------------------------------------------------------ s03 端到端
def _s03_setup():
    cfg = _common.cfg()
    L = cfg.species("mouse").labels
    rows = []
    for grp in ("Purkinje", "GranuleLineage", "Astro", "Immune", "ALL"):
        for i in range(1, 6):
            vals = {"GT": {1: 1.0, 2: 1.0, 3: 20.0, 4: 1.0, 5: 1.0}[i]
                    if grp in ("Purkinje", "ALL") else 1.0, "GU": 1.0}
            rows.append((grp, L[i], True, 40, vals))
    rows.append(("GranuleLineage", L[6], True, 40, {"GT": 1.0, "GU": 1.0}))
    for i in range(1, 6):
        rows.append(("Interneuron", L[i], False, 5, {"GT": np.nan}))
    bt = synth.mk_bin_table(cfg, rows, ["GT", "GU"], value_kind="cp")
    tabs = s02_stats.build_tables(bt, None, ["GT", "GU", "GM"], cfg.species("mouse"), cfg)
    return cfg, bt, tabs


def test_s03_end_to_end_synth():
    """write_tables → read_tables roundtrip → build_figures:8 视图命名/双格式齐;
    缺失基因 GM 与不可行组不炸(NA 路径);两跑文件集一致(确定性)。"""
    import matplotlib.pyplot as plt
    cfg, bt, tabs = _s03_setup()
    spcfg = cfg.species("mouse")
    with tempfile.TemporaryDirectory() as td:
        tdir, fdir = os.path.join(td, "tables"), os.path.join(td, "figs")
        s02_stats.write_tables(tabs, tdir)
        back = s03.read_tables(tdir)
        assert set(back) == set(s02_stats.TABLE_NAMES)
        assert (back["emergence"].set_index("gene").loc["GT", "emergence_bin"]
                == tabs["emergence"].set_index("gene").loc["GT", "emergence_bin"] != "NA")
        try:
            s03.read_tables(os.path.join(td, "nowhere"))
            raise AssertionError("缺表未 RuntimeError")
        except RuntimeError as e:
            assert "s02" in str(e)
        paths = s03.build_figures(back, bt, None, ["GT", "GU", "GM"], spcfg, cfg, fdir)
        assert set(paths) == {"a1", "b1", "b2", "b3", "b4", "c1", "c2", "d1"}
        stems = {os.path.basename(p).rsplit(".", 1)[0] for ps in paths.values() for p in ps}
        # 命名(spec §5): a1 带 ALL;b1/b2 按细胞组(合成里 5 组 − ALL = 4 组有可行/全组出图);
        # b3/b4/c1/d1 单文件;c2 分页 1 页
        for must in ("mouse_a1_ALL_p1", "mouse_b1_Purkinje_p1", "mouse_b2_Purkinje_p1",
                     "mouse_b3", "mouse_b4", "mouse_c1", "mouse_c2_p1", "mouse_d1"):
            assert must in stems, f"缺 {must}: {sorted(stems)}"
        assert "mouse_b2_GranuleLineage_p1" in stems
        assert not any(s.endswith("_ALL_p1") and "_b2_" in s for s in stems)   # ALL 不出 b2
        for p in (p for ps in paths.values() for p in ps):
            assert os.path.getsize(p) > 0
            assert os.path.splitext(p)[1] in {".pdf", ".png"}
        n1 = len(stems)
        paths2 = s03.build_figures(back, bt, None, ["GT", "GU", "GM"], spcfg, cfg, fdir)
        stems2 = {os.path.basename(p).rsplit(".", 1)[0] for ps in paths2.values() for p in ps}
        assert stems2 == stems and n1 == len(stems2)
        plt.close("all")
