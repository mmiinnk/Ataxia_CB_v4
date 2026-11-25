#!/usr/bin/env python
"""s06 — 物种报告 + review 门清单自动生成(spec §9 四条清单;内容口径 §5A,P7 提取补注)。

用法:
    python -I scripts/s06_report.py <root> --species mouse [--config PATH]

前置: s01(cube_bin.tsv + cube_donor.h5ad)。s02/s03 层缺席时对应小节整节跳过
(§5A: 门 1 先于 s07/s08, 报告只呈现已跑层, 不产占位)。

产物:
  results/v{ver}/reports/<sp>/species_report.md   物种报告(九小节, §5A)
  reviews/gate{1|2}_<sp>.md                       §9 四条清单(自动项已填, 人工项留复选框)
  logs/s06_report_<sp>.tsv                        本脚本自身 REPORT

防覆盖(§5A): reviews/gate{N}_<sp>.md 已存在时改写 gate{N}_<sp>.refresh.md——
人工勾选记录不被重跑破坏; 正式采纳哪份由 review 人决定后手工处理。

设计: build_species_report()/build_gate_checklist() 为纯函数(输入显式传入,
outdir/reviews 路径由 main 解析); 无时间戳、无随机 → 同输入逐字节一致, 可 diff。
logs/*_report*.tsv schema 不统一(s00 多 dataset/expected 列, s01 多 expected 列),
load_logs() 归一为 check/observed/level/note + source。
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from atx import config as atxconfig
from atx import genes as atxgenes
import s03_figures_species as s03

REPORT: list[dict] = []


def rec(check, observed, level, note=""):
    REPORT.append({"check": check, "observed": str(observed), "level": level, "note": note})
    mark = {"PASS": "✓", "WARN": "⚠", "FAIL": "✗"}.get(level, "?")
    print(f"[{mark} {level}] {check} | {observed}" + (f" | {note}" if note else ""))


GATE_OF = {"mouse": 1, "human": 2}
# §5A⑥: 本物种相关 §11 未决条目编号(条文以 spec §11 为准, 代码只列编号防漂移)
S11_FLAGS = {"mouse": [3, 8], "human": [2, 7, 9]}
LOG_SOURCES = ("s00_report", "s01_report_{sp}", "s02_report_{sp}", "s03_report_{sp}")


def load_logs(logdir: str, sp: str) -> dict[str, pd.DataFrame]:
    """logs/*_report*.tsv → {source: 归一 DataFrame(check/observed/level/note)}。

    归一规则: check 前置 dataset(s00 有 dataset 列时); expected 列并入 observed
    (显示为 'expected→observed'), 不丢信息。缺席文件跳过(门 1 前 s02/s03 可能未跑)。
    """
    out: dict[str, pd.DataFrame] = {}
    for tpl in LOG_SOURCES:
        src = tpl.format(sp=sp)
        p = os.path.join(logdir, f"{src}.tsv")
        if not os.path.exists(p):
            continue
        df = pd.read_csv(p, sep="\t", keep_default_na=False, na_values=[""])
        for c in ("check", "observed", "level", "note"):
            if c not in df.columns:
                df[c] = ""
        df = df.copy()
        df["check"] = df["check"].astype(str)
        if "dataset" in df.columns:
            ds = df["dataset"].astype(str)
            df["check"] = np.where(ds.str.len() > 0, ds + ": " + df["check"], df["check"])
        if "expected" in df.columns:
            ex = df["expected"].astype(str)
            df["observed"] = np.where(ex.str.len() > 0,
                                      ex + " → " + df["observed"].astype(str),
                                      df["observed"].astype(str))
        df["source"] = src
        out[src] = df[["source", "check", "observed", "level", "note"]].fillna("")
    return out


def scan_levels(logs: dict[str, pd.DataFrame]) -> tuple[list, list]:
    """(fails, warns): 逐行 (source, check, observed, note) 元组。"""
    fails, warns = [], []
    for df in logs.values():
        for r in df.to_dict("records"):
            tup = (r["source"], str(r["check"]), str(r["observed"]), str(r.get("note", "") or ""))
            if r["level"] == "FAIL":
                fails.append(tup)
            elif r["level"] == "WARN":
                warns.append(tup)
    return fails, warns


def _rows(logs: dict[str, pd.DataFrame], pred) -> list:
    out = []
    for df in logs.values():
        for r in df.to_dict("records"):
            if pred(r):
                out.append(r)
    return out


def _md_table(rows: list[list[str]], header: list[str]) -> str:
    if not rows:
        return "_(无行)_\n"
    lines = ["| " + " | ".join(header) + " |",
             "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        lines.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(lines) + "\n"


def _fmt(v, nd=3):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    return "NA" if not np.isfinite(f) else f"{f:.{nd}f}"


def _num(s: pd.Series) -> pd.Series:
    """数值列强制转 float(read_tables 保留受控 'NA' 串; 混合 dtype 会让比较炸)。"""
    return pd.to_numeric(s, errors="coerce")


def _bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


# ---------------------------------------------------------------------------
# 物种报告(§5A 九小节)
# ---------------------------------------------------------------------------

def build_species_report(bt: pd.DataFrame, cube, tables, logs: dict[str, pd.DataFrame],
                         fig_stems: list[str], cfg, sp: str, spcfg) -> str:
    """纯函数: 已跑层的实据 → markdown 文本。tables=None → ③④⑤中趋势锚点跳过。"""
    ver = cfg.version
    L: list[str] = []
    ap = L.append
    ap(f"# {sp} 物种报告 — gene_set v{ver}" + ("(DRAFT 探索性口径)" if cfg.draft else ""))
    ap("")

    # ① 头部
    ap("## 1. 头部")
    panel = atxgenes.load_panel(cfg)
    genes = panel.loc[panel["in_stat_panel"].astype(bool), "gene_key"].astype(str).tolist()
    in_cube = sum(int(f"log_{g}" in bt.columns) for g in genes)
    n_don = n_dis = 0
    if cube is not None:
        n_don = int(cube.obs["donor"].nunique())
        if "disease_flag" in cube.obs:
            n_dis = int(cube.obs.loc[cube.obs["disease_flag"].astype(bool), "donor"].nunique())
    n_feas = int(bt["feasible"].astype(bool).sum())
    ap(f"- 基因表版本: v{ver}(config gene_set.version; 冻结状态见 genes/curation_log.md)")
    ap(f"- 统计 panel: {len(genes)} 基因, cube 命中 {in_cube}, 缺失 {len(genes) - in_cube}(status=missing)")
    ap(f"- 供体: {n_don}" + (f"(疾病供体 {n_dis}, donor cube disease_flag 口径)" if n_dis else ""))
    ap(f"- cube_bin: {len(bt)} 组×bin 行, 可行 {n_feas}, 不可行 {len(bt) - n_feas}")
    ap("")

    # ② 证据结构
    ap("## 2. 证据结构(可行性 / 单供体 / pooled)")
    groups = [g for g in cfg.plot_groups() if g in set(bt["group"].astype(str))]
    bins = sorted(bt["bin"].astype(str).unique(), key=lambda b: spcfg.bin_order.get(b, 999))
    feas = bt.set_index(["group", "bin"])["feasible"].astype(bool)
    ncell = bt.set_index(["group", "bin"])["n_cells_total"]
    rows = []
    for g in groups:
        cells = []
        for b in bins:
            if (g, b) in feas.index:
                cells.append("✓" if bool(feas.loc[(g, b)]) else f"✗{int(ncell.loc[(g, b)])}")
            else:
                cells.append("·")
        rows.append([g] + cells)
    ap(_md_table(rows, ["group"] + bins))
    ap("✓=可行(≥min_cells);✗n=不可行(n 核);·=骨架无此行。")
    sd = bt.loc[bt["single_donor"].astype(bool), ["group", "bin"]]
    pl = bt.loc[bt["pooled"].astype(bool), ["group", "bin"]]
    ap(f"- 单供体 bin(sem=NA, 硬约束 5): {len(sd)} 行 — "
       + (", ".join(f"{g}@{b}" for g, b in sd.itertuples(index=False)) or "无"))
    ap(f"- pooled bin(sem=NA): {len(pl)} 行 — "
       + (", ".join(f"{g}@{b}" for g, b in pl.itertuples(index=False)) or "无"))
    ap("")

    if tables is not None:
        tr, tp = tables["trend_by_group"], tables["temporal_patterns"]
        sg, em = tables["specificity_global"], tables["emergence"]
        pk = tables["peak_by_group"]

        # ③ 第 1 层趋势
        ap("## 3. 第 1 层: 趋势(Spearman rho_bin, BH q<0.05, 族=species×group)")
        rows = []
        for g in groups:
            sub = tr[tr["group"].astype(str) == g]
            q = _num(sub["q_bin"])
            rho = _num(sub["rho_bin"])
            sigm = (q < 0.05) & q.notna()
            up = int(((rho > 0) & sigm).sum())
            dn = int(((rho < 0) & sigm).sum())
            rows.append([g, len(sub), int(sigm.sum()), up, dn, int(rho.isna().sum())])
        ap(_md_table(rows, ["group", "基因行", "q<0.05", "升", "降", "rho NA"]))
        ap("")

        # ④ 第 2 层
        ap("## 4. 第 2 层: 模式 / 特异性 / emergence / 峰期")
        vc = tp["temporal_pattern"].astype(str).value_counts()
        ap("- temporal_pattern 分布(gene×group 行): "
           + ", ".join(f"{k}={v}" for k, v in vc.items()))
        svc = sg["specific_call"].astype(str).value_counts()
        ap("- specific_call 分布: " + ", ".join(f"{k}={v}" for k, v in svc.items()))
        spec = sg.loc[sg["specific_call"].astype(str) == "specific"]
        if len(spec):
            ap("- specific 基因(τ 与主组): " + ", ".join(
                f"{r.gene}→{r.max_group}(τ={_fmt(r.tau)}, frac={_fmt(r.max_fraction)})"
                for r in spec.assign(_t=_num(spec["tau"])).sort_values("_t", ascending=False).itertuples()))
        ne = sg.loc[sg["specific_call"].astype(str) == "not_expressed", "gene"].astype(str).tolist()
        if ne:
            ap(f"- not_expressed 基因({len(ne)}): " + ", ".join(sorted(ne)))
        # P7 裁决①(2026-11-05): not_expressed 基因的 window/peak 标签已在 s02 掩膜为 NA
        masked = tp[(tp["status"].astype(str) == "ok")
                    & (tp["temporal_pattern"].astype(str) == "NA")]
        ap(f"- 低表达伪标签掩膜(P7 裁决①): {len(masked)} 行 window/peak→NA"
           + ("(" + ", ".join(sorted(masked["gene"].astype(str).unique())) + ")" if len(masked) else ""))
        wp = tp[tp["temporal_pattern"].astype(str).isin(["window_high", "transient_peak"])]
        contam = wp[wp["gene"].astype(str).isin(ne)]
        if len(contam):                      # 守卫: 掩膜生效则必为 0
            ap(f"- **异常**: 窗口/峰标签中仍有 not_expressed 基因 {len(contam)} 行——掩膜未生效, 须查 s02")
        resid = tp[(tp["temporal_pattern"].astype(str) == "endpoint_high")
                   & tp["gene"].astype(str).isin(ne)]
        if len(resid):
            ap(f"- endpoint_high 同类残余(不在裁决①范围, 留档): {len(resid)} 行("
               + ", ".join(sorted(resid["gene"].astype(str).unique())) + ")")
        n_em = int((em["emergence_bin"].astype(str) != "NA").sum())
        n_lc = int(_bool(em["left_censored"]).sum())
        ap(f"- emergence 可判: {n_em}/{len(em)}(左截断 {n_lc})")
        per: list[str] = []          # 期序 = bin 时间序首次出现序(spcfg 派生, 不硬编码)
        for b in sorted(spcfg.bin_order, key=lambda x: spcfg.bin_order[x]):
            p = spcfg.period_of.get(b)
            if p and p not in per:
                per.append(p)
        if "peak_period" in pk.columns:
            piv = (pk.assign(_g=pk["group"].astype(str))
                     .pivot_table(index="_g", columns="peak_period", values="gene",
                                  aggfunc="count", fill_value=0))
            order = [g for g in groups if g in piv.index]
            piv = piv.loc[order]
            if per:
                piv = piv[[c for c in per if c in piv.columns]]
            ap("- 峰期×组分布(基因数, peak_by_group):")
            ap(_md_table([[i] + [int(v) for v in r] for i, r in piv.iterrows()],
                         ["group"] + [str(c) for c in piv.columns]))
        ap("")

    # ⑤ 锚点与门自动核查
    ap("## 5. 锚点与门自动核查(§4 + §9)")
    fails, warns = scan_levels(logs)
    ap(f"- 全部 logs: FAIL {len(fails)} 条, WARN {len(warns)} 条(逐条见 reviews/gate{GATE_OF[sp]}_{sp}.md)")
    anchor_rows = _rows(logs, lambda r: "锚点" in str(r["check"]))
    rows = [[r["source"], r["check"], str(r["observed"])[:120], r["level"]] for r in anchor_rows]
    ap(_md_table(rows, ["来源", "check", "observed", "level"]))
    if warns:
        ap("WARN 清单:")
        rows = [[r[0], r[1], r[2][:120], r[3]] for r in warns]
        ap(_md_table(rows, ["来源", "check", "observed", "note"]))
    ap("")

    # ⑥ §11 未决条目
    ap("## 6. 本物种相关 §11 未决条目")
    ap("- 条目编号: " + ", ".join(str(i) for i in S11_FLAGS[sp])
       + "(条文以 implementation_spec §11 为准; 报告只列编号防漂移)")
    if sp == "mouse":
        ap("- 注: §11.1 鼠免疫 marker 按 dropout 豁免(Cx3cr1/P2ry12/Hexb 口径), 见门清单条目(3)。")
    ap("")

    # ⑦ 图清单
    ap("## 7. 图清单(figures/{})".format(sp))
    if fig_stems:
        views: dict[str, list[str]] = {}
        for stem in fig_stems:
            parts = stem.split("_")
            v = parts[1] if len(parts) > 1 else "?"
            views.setdefault(v, []).append(stem)
        rows = [[v, len(ss), ", ".join(sorted(ss)[:6]) + ("…" if len(ss) > 6 else "")]
                for v, ss in sorted(views.items())]
        ap(_md_table(rows, ["视图", "文件茎数", "茎(前 6)"]))
    else:
        ap("_(figures/ 无产物 — s03 未跑或已跳过)_")
    ap("")

    # ⑧ review 焦点
    ap("## 8. review 焦点(spec §9: b2 轨迹与 c1 dotplot)")
    figdir = cfg.results("figures", sp)
    focus = [s for s in fig_stems if "_b2_" in s and s.endswith("p1")] + \
            [s for s in fig_stems if "_c1" in s]
    for stem in sorted(focus):
        exts = cfg.raw["figures"]["formats"]
        ap(f"- {os.path.join(figdir, stem + '.' + exts[0])}")
    if not focus:
        ap("_(s03 未跑: 门 review 前先补 s03)_")
    ap("")

    # ⑨ 口径注
    ap("## 9. 口径注")
    ap("- NA 政策(硬约束 5): 组×bin < min_cells → 表 NA、图斜纹; 不插值、不填 0; 单供体/pooled bin 的 sem=NA。")
    ap("- b2 轨迹 y 轴 = mean log1p CP10K(与存储的 logsem_/rho_bin 自洽; 读图指南 'y=CP10K' 的落地注记)。")
    ap("- c1 dotplot 点大小 = pct(0–1 表达核比例)× dot_size_max; 色 = 可行 bin log 均值。")
    # DRAFT 声明必须随 gene_set 版本条件生成(atx.plots.add_footer 仅在 cfg.draft 时加 footer):
    # v1.0 冻结后图形不带 DRAFT 标, 报告若仍写"全部图形带 DRAFT footer"即为过时文字, 会误导评审。
    if cfg.draft:
        ap(f"- 全部图形带 DRAFT footer(v{ver} 未冻结): 数值为探索性口径, 不得引用为论文数字。")
    else:
        ap(f"- gene_set v{ver} 已冻结 → 图形**不带** DRAFT footer, 数值可引用; 引用时须附口径注"
           "(min_cells 阈值、可行 组×bin 数、单供体/pooled bin 的 sem=NA、右删失的组与 bin 段)。")
    ap("- 低表达伪标签(P7 裁决①, 2026-11-05): not_expressed 基因的 window_high/transient_peak 已在 s02 掩膜为 NA(图灰色/斜纹, 下游自动排除); endpoint_high 同类残余仅计数留档。")
    ap("")
    return "\n".join(L)


# ---------------------------------------------------------------------------
# 门清单(§9 四条)
# ---------------------------------------------------------------------------

def build_gate_checklist(logs: dict[str, pd.DataFrame], cfg, sp: str) -> tuple[str, bool]:
    """纯函数 → (markdown, overall_ok)。四条中 (1)(2)(4) 自动判级, (3) 自动摘录+人工确认。"""
    gate = GATE_OF[sp]
    fails, warns = scan_levels(logs)
    anchors = _rows(logs, lambda r: "锚点" in str(r["check"]))
    markers = _rows(logs, lambda r: "marker" in str(r["check"]).lower())
    trend_anchors = [r for r in anchors if str(r["check"]).startswith("趋势锚点")
                     and "(正式)" in str(r["check"])]

    def grade(rows) -> str:
        """三态判级: pass / warn(FAIL=0 但有非 PASS, 待人工裁决) / fail / empty(无行)。"""
        if not rows:
            return "empty"
        if any(r["level"] == "FAIL" for r in rows):
            return "fail"
        if any(r["level"] != "PASS" for r in rows):
            return "warn"
        return "pass"

    g1 = "pass" if not fails else "fail"
    g2, g3 = grade(anchors), grade(markers)
    # 条目(4): 锚点定义 = config species.<sp>.trend_anchors(2026-11-05 裁决人鼠均设)。
    # s02 对方向不符记 WARN → warn 也判不过; config 未设 → n/a; 设了但无行 = s02 未跑 → 不过。
    cfg_anchors = cfg.trend_anchors(sp)
    if not cfg_anchors:
        g4 = "n/a"
    else:
        g4 = grade(trend_anchors)
        if g4 == "warn":
            g4 = "fail"
    overall = (g1 == "pass" and g2 in ("pass", "warn")
               and g3 in ("pass", "warn") and g4 in ("pass", "n/a"))

    MARK = {"pass": "✓ 自动判定通过", "warn": "⚠ 有 WARN, 待人工裁决(FAIL 0)",
            "fail": "✗ 自动判定未过", "empty": "✗ 无行(层未跑或未设)", "n/a": "— 未设"}

    def mark(g):
        return MARK[g]

    L: list[str] = []
    ap = L.append
    ap(f"# review 门 {gate} — {sp}(spec §9 四条清单; s06 自动生成自动项, 人工项勾选后 commit)")
    ap("")
    ap(f"自动项总状态: {'通过' if overall else '未通过'}"
       f"(FAIL {len(fails)} / WARN {len(warns)} / 锚点 {len(anchors)} / marker {len(markers)} 条)")
    ap("")

    ap("## (1) logs 无 assert 失败且 WARNING 已逐条看过")
    ap(f"- 自动: {mark(g1)} — FAIL {len(fails)} 条"
       + ("" if g1 == "pass" else "(须清零或逐条豁免后才可过门)"))
    if fails:
        ap(_md_table([[f[0], f[1], f[2][:120]] for f in fails], ["来源", "check", "observed"]))
    ap(f"- WARN {len(warns)} 条(判读后勾选):")
    if warns:
        ap(_md_table([[w[0], w[1], w[2][:120], w[3][:80]] for w in warns],
                     ["来源", "check", "observed", "note"]))
    ap("- [ ] 人工确认: WARNING 已逐条看过且非阻断")
    ap("")

    ap("## (2) §4 锚点全部匹配")
    ap(f"- 自动: {mark(g2)}(锚点行 {len(anchors)} 条; FAIL→不过, WARN→人工裁决)")
    ap(_md_table([[r["source"], r["check"], str(r["observed"])[:120], r["level"]] for r in anchors],
                 ["来源", "check", "observed", "level"]))
    ap("- [ ] 人工确认: 锚点数值与 spec §4 逐项核对")
    ap("")

    ap("## (3) marker 校验通过")
    ap(f"- 自动: {mark(g3)}(marker 行 {len(markers)} 条; FAIL→不过, WARN→人工裁决)")
    ap(_md_table([[r["source"], r["check"], str(r["observed"])[:120], r["level"]] for r in markers],
                 ["来源", "check", "observed", "level"]))
    if sp == "mouse":
        ap("- 注: 鼠免疫 marker 按 dropout 豁免(spec §11.1: Neurod1/Ptprc/Aif1/Csf1r/Tmem119;"
           " 免疫锚点改用 Cx3cr1/P2ry12/Hexb)。")
    ap("- [ ] 人工确认: ITPR1/CA8→Purkinje, GABRA6/CBLN3→GranuleLineage, MBP/PLP1→OligoLineage,"
       " AQP4→Astro, SOX2/MKI67→Progenitor, PTPRC→Immune(鼠免疫豁免如上)")
    ap("")

    ap("## (4) 趋势锚点方向正确")
    if cfg_anchors:
        expect = "; ".join(f"{a['gene']}×{a['group']} ρ{'>' if a['sign'] == '+' else '<'}0"
                           for a in cfg_anchors)
        ap(f"- 自动: {mark(g4)}(正式趋势锚点 {len(trend_anchors)} 条; §9(4) 期望: {expect};"
           " 方向不符记 WARN 亦判不过)")
        if not trend_anchors:
            ap("  (无正式趋势锚点行: s02 未跑, 或 config 锚点基因不在 cube——过门前须补齐)")
        ap(_md_table([[r["check"], str(r["observed"])[:120], r["level"]] for r in trend_anchors],
                     ["check", "observed", "level"]))
        ap("- 锚点来源(config trend_anchors.source):")
        for a in cfg_anchors:
            ap(f"  - {a['gene']}×{a['group']}({a['sign']}): {a['source']}")
    else:
        ap("- 自动: — config species.{}.trend_anchors 未设; 冻结前须补(spec §9(4))。".format(sp))
    ap("- [ ] 人工确认: 方向与量级合理")
    ap("")
    ap("---")
    ap("reviewer: ____________  日期: ____________  结论: [ ] 过门 / [ ] 整改后复检")
    ap("")
    return "\n".join(L), bool(overall)


def write_review(root: str, sp: str, text: str) -> str:
    """reviews/gate{N}_{sp}.md; 已存在 → .refresh.md(防覆盖人工勾选, §5A)。"""
    rdir = os.path.join(root, "reviews")
    os.makedirs(rdir, exist_ok=True)
    path = os.path.join(rdir, f"gate{GATE_OF[sp]}_{sp}.md")
    if os.path.exists(path):
        path = os.path.join(rdir, f"gate{GATE_OF[sp]}_{sp}.refresh.md")
    with open(path, "w") as f:
        f.write(text)
    return path


def main():
    ap_ = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap_.add_argument("root")
    ap_.add_argument("--species", required=True, choices=["human", "mouse"])
    ap_.add_argument("--config", default=None)
    args = ap_.parse_args()
    sp = args.species
    t0 = time.time()
    root = os.path.abspath(args.root)
    cfg = atxconfig.load(args.config or os.path.join(root, "config", "config.yaml"))
    cfg.check_version_gate()
    spcfg = cfg.effective_species(sp)

    print(f"== s06 {sp} | version={cfg.version}(draft={cfg.draft})")
    inter = cfg.results("intermediate", sp)
    bt_path = os.path.join(inter, "cube_bin.tsv")
    if not os.path.exists(bt_path):
        raise RuntimeError(f"缺 s01 产物 {bt_path}——先跑 s01(spec §9 顺序)")
    bt = pd.read_csv(bt_path, sep="\t")

    cube = None
    cube_path = os.path.join(inter, "cube_donor.h5ad")
    if os.path.exists(cube_path):
        import anndata
        cube = anndata.read_h5ad(cube_path)

    tdir = cfg.results("tables", sp)
    tables = None
    if all(os.path.exists(os.path.join(tdir, f"{n}.tsv")) for n in s03.s02_stats.TABLE_NAMES):
        tables = s03.read_tables(tdir)
    else:
        rec("s02 六表", "缺席 — 报告跳过第 3/4 小节(§5A 未跑层整节跳过)", "WARN")

    fdir = cfg.results("figures", sp)
    fig_stems = sorted({os.path.basename(p).rsplit(".", 1)[0]
                        for p in (os.path.join(fdir, f) for f in os.listdir(fdir))}
                       ) if os.path.isdir(fdir) else []
    if not fig_stems:
        rec("figures 目录", "缺席或空 — 报告第 7/8 小节记未跑", "WARN")

    logs = load_logs(cfg.results("logs"), sp)
    rec("logs 装载", ", ".join(f"{k}({len(v)})" for k, v in logs.items()) or "无", "PASS")

    report = build_species_report(bt, cube, tables, logs, fig_stems, cfg, sp, spcfg)
    checklist, overall = build_gate_checklist(logs, cfg, sp)

    rdir = cfg.results("reports", sp, mkdir=True)
    rp = os.path.join(rdir, "species_report.md")
    with open(rp, "w") as f:
        f.write(report)
    gp = write_review(cfg.root, sp, checklist)   # cfg.root = paths.project_root(硬约束 6 单一来源)

    fails, warns = scan_levels(logs)
    rec("物种报告", f"{len(report.splitlines())} 行 → {rp}", "PASS")
    rec("门清单", f"门{GATE_OF[sp]} 自动项 {'通过' if overall else '未通过'} → {gp}",
        "PASS" if overall else "WARN",
        note=f"FAIL {len(fails)} / WARN {len(warns)}; 人工勾选后 commit 留痕")

    logdir = cfg.results("logs", mkdir=True)
    lrp = os.path.join(logdir, f"s06_report_{sp}.tsv")
    pd.DataFrame(REPORT, columns=["check", "observed", "level", "note"]).to_csv(lrp, sep="\t", index=False)
    print(f"\n== s06 {sp} 完成 ({time.time()-t0:.0f}s)")
    print(f"报告: {rp}\n门清单: {gp}\n自身报告: {lrp}")
    print(f"下一步: 人工 review 门 {GATE_OF[sp]}(勾选 {gp} 并 commit)。")
    return 1 if any(r["level"] == "FAIL" for r in REPORT) else 0


if __name__ == "__main__":
    sys.exit(main())
