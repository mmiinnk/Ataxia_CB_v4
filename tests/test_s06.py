"""s06 物种报告 + 门清单合成测试(spec §9 四条清单;§5A 内容口径)。

覆盖: load_logs 对三种 schema(s00 多 dataset/expected, s01 多 expected, s02 标准)
的归一; 报告九小节齐全 + 缺 s02 层整节跳过 + 确定性; 门清单自动判级(mouse 四条,
human 第 4 条"未设"注记)、FAIL 传播; write_review 防覆盖(.refresh.md)。
全部临时目录, 不触真实 results/ 与 reviews/。
"""
from __future__ import annotations

import os
import tempfile

import numpy as np
import pandas as pd

import _common
import s02_stats
import s06_report as s06
import synth


# ------------------------------------------------------------------ 夹具
def _tables_setup():
    """合成 bt + s02 七表(与 test_plots._s03_setup 同源布局)。"""
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


def _logs_setup(tmp, extra_fail=False, sp="mouse", write_s02=True):
    """三种 schema 的假日志落盘 → load_logs 归一。"""
    s00 = pd.DataFrame([
        {"check": "9 组细胞数锚点 §4.1", "dataset": sp, "expected": "全部精确匹配",
         "observed": "匹配", "level": "PASS", "note": ""},
        {"check": "panel 命中率", "dataset": sp, "expected": "≥0.9",
         "observed": "0.945", "level": "PASS", "note": ""},
        {"check": "external 文件 brainspan", "dataset": "external", "expected": "存在",
         "observed": "缺席", "level": "WARN", "note": "s10 占位(§11.6)"},
    ])
    s01 = pd.DataFrame([
        {"check": "marker 校验(spec §9 清单)", "expected": "FAIL=0",
         "observed": "FAIL=0(清单内) / 1(清单外), ABSENT=2", "level": "PASS", "note": "免疫 dropout 豁免"},
        {"check": "组×bin 可行性矩阵 §4.2", "expected": "匹配",
         "observed": "匹配", "level": "PASS", "note": ""},
    ])
    s02 = pd.DataFrame([
        {"check": "趋势锚点 ITPR1×Purkinje(正式)", "observed": "rho + | rho_bin=+0.915 q_bin=0.001",
         "level": "PASS", "note": ""},
        {"check": "趋势锚点 MBP×OligoLineage(正式)", "observed": "rho + | rho_bin=+0.943",
         "level": "PASS", "note": ""},
        {"check": "趋势锚点 ATOH1×GranuleLineage(正式)", "observed": "rho − | rho_bin=-0.907",
         "level": "PASS", "note": ""},
    ])
    if extra_fail:
        s02 = pd.concat([s02, pd.DataFrame([
            {"check": "X 计数整数性", "observed": "发现非整数", "level": "FAIL", "note": "阻断"}])],
            ignore_index=True)
    os.makedirs(tmp, exist_ok=True)
    s00.to_csv(os.path.join(tmp, "s00_report.tsv"), sep="\t", index=False)
    s01.to_csv(os.path.join(tmp, f"s01_report_{sp}.tsv"), sep="\t", index=False)
    if write_s02:
        s02.to_csv(os.path.join(tmp, f"s02_report_{sp}.tsv"), sep="\t", index=False)
    return s06.load_logs(tmp, sp)


# ------------------------------------------------------------------ 测试
def test_load_logs_normalization():
    with tempfile.TemporaryDirectory() as td:
        logs = _logs_setup(td)
        assert set(logs) == {"s00_report", "s01_report_mouse", "s02_report_mouse"}
        s00 = logs["s00_report"]
        assert list(s00.columns) == ["source", "check", "observed", "level", "note"]
        # dataset 前置 + expected 并入 observed
        assert s00["check"].iloc[0] == "mouse: 9 组细胞数锚点 §4.1"
        assert s00["observed"].iloc[0] == "全部精确匹配 → 匹配"
        # 缺席文件跳过(s03 未落盘不炸)
        assert "s03_report_mouse" not in logs
        fails, warns = s06.scan_levels(logs)
        assert fails == [] and len(warns) == 1 and "brainspan" in warns[0][1]
        logs_f = _logs_setup(os.path.join(td, "f2"), extra_fail=True)
        fails, _ = s06.scan_levels(logs_f)
        assert len(fails) == 1 and fails[0][1] == "X 计数整数性"


def test_species_report_sections_and_determinism():
    cfg, bt, tabs = _tables_setup()
    spcfg = cfg.species("mouse")
    with tempfile.TemporaryDirectory() as td:
        logs = _logs_setup(td)
        stems = ["mouse_a1_ALL_p1", "mouse_b2_Purkinje_p1", "mouse_c1", "mouse_d1"]
        rep = s06.build_species_report(bt, None, tabs, logs, stems, cfg, "mouse", spcfg)
        for sec in ("## 1. 头部", "## 2. 证据结构", "## 3. 第 1 层", "## 4. 第 2 层",
                    "## 5. 锚点", "## 6. 本物种相关 §11", "## 7. 图清单",
                    "## 8. review 焦点", "## 9. 口径注"):
            assert sec in rep, f"缺小节 {sec}"
        assert "✓" in rep and "Interneuron" in rep          # 可行性矩阵(含不可行组)
        assert "条目编号: 3, 8" in rep                       # §5A⑥ 鼠条目
        assert "b2_Purkinje_p1" in rep                      # review 焦点列出 b2/c1
        # 缺 s02 层 → ③④ 整节跳过, 其余小节仍在(§5A)
        rep0 = s06.build_species_report(bt, None, None, logs, [], cfg, "mouse", spcfg)
        assert "## 3. 第 1 层" not in rep0 and "## 4. 第 2 层" not in rep0
        assert "## 5. 锚点" in rep0 and "figures/ 无产物" in rep0
        # 确定性: 同输入两次逐字节一致(无时间戳)
        assert s06.build_species_report(bt, None, tabs, logs, stems, cfg, "mouse", spcfg) == rep


def test_gate_checklist_mouse_and_human():
    cfg = _common.cfg()
    with tempfile.TemporaryDirectory() as td:
        logs = _logs_setup(td)
        text, ok = s06.build_gate_checklist(logs, cfg, "mouse")
        assert ok and "review 门 1 — mouse" in text
        assert text.count("✓ 自动判定通过") == 4             # 四条全过
        assert "ITPR1" in text and "dropout 豁免" in text    # 趋势锚点值 + 鼠免疫豁免注记
        assert "- [ ] 人工确认" in text and "reviewer:" in text
        # human: 门号 2;人侧锚点已裁决增设(2026-11-05, config trend_anchors)→ 四条全过
        logs_h = _logs_setup(os.path.join(td, "h"), sp="human")
        text_h, ok_h = s06.build_gate_checklist(logs_h, cfg, "human")
        assert ok_h and "review 门 2 — human" in text_h
        assert text_h.count("✓ 自动判定通过") == 4
        assert "锚点来源" in text_h and "仅设鼠趋势锚点" not in text_h
        # FAIL 传播: 条目(1) 不过 → 总状态未通过
        logs_f = _logs_setup(os.path.join(td, "f"), extra_fail=True)
        text_f, ok_f = s06.build_gate_checklist(logs_f, cfg, "mouse")
        assert not ok_f and "✗ 自动判定未过" in text_f and "X 计数整数性" in text_f
        # s02 未跑(无正式趋势锚点行)→ 鼠条目(4) 判不过, 不得因空行集误判通过
        logs_n2 = _logs_setup(os.path.join(td, "n2"), write_s02=False)
        text_n2, ok_n2 = s06.build_gate_checklist(logs_n2, cfg, "mouse")
        assert not ok_n2 and "✗ 无行(层未跑或未设)" in text_n2
        # 三态: 锚点行含 WARN(FAIL=0)→ ⚠ 待人工裁决, 总状态仍算自动项通过(人工裁决定夺)
        logs_w = dict(_logs_setup(os.path.join(td, "w")))
        df = logs_w["s01_report_mouse"].copy()
        df.loc[df["check"].str.contains("marker"), "level"] = "WARN"
        logs_w["s01_report_mouse"] = df
        text_w, ok_w = s06.build_gate_checklist(logs_w, cfg, "mouse")
        assert ok_w and "⚠ 有 WARN, 待人工裁决" in text_w
        # 同场景但 WARN 落在趋势锚点行 → 鼠条目(4) 判不过(方向不符即 WARN, 不放行)
        df2 = logs_w["s02_report_mouse"].copy()
        df2.loc[df2["check"].str.contains("ITPR1"), "level"] = "WARN"
        logs_w2 = dict(logs_w, s02_report_mouse=df2)
        _t, ok_w2 = s06.build_gate_checklist(logs_w2, cfg, "mouse")
        assert not ok_w2
        # config 未设 trend_anchors 的物种 → 条目(4) "— 未设", 不阻断其余三条
        import copy
        from atx import config as atxconfig
        raw = copy.deepcopy(cfg.raw)
        raw["species"]["mouse"].pop("trend_anchors")
        cfg_na = atxconfig.Cfg(raw, cfg.path)
        text_na, ok_na = s06.build_gate_checklist(logs, cfg_na, "mouse")
        assert ok_na and "trend_anchors 未设" in text_na


def test_write_review_preserves_human_file():
    with tempfile.TemporaryDirectory() as td:
        p1 = s06.write_review(td, "mouse", "第一版\n")
        assert p1.endswith(os.path.join("reviews", "gate1_mouse.md"))
        with open(p1, "a") as f:
            f.write("- [x] 人工确认\n")                      # 模拟人工勾选
        p2 = s06.write_review(td, "mouse", "重跑版\n")
        assert p2.endswith("gate1_mouse.refresh.md")         # 已存在 → .refresh
        assert "人工确认" in open(p1).read()                 # 原文件未被破坏
        assert open(p2).read() == "重跑版\n"
