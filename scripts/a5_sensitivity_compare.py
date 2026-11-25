#!/usr/bin/env python3
"""a5_sensitivity_compare.py — 敏感性变体横向比较(一次性分析层, 方案A 配套)。

对照基准 = results/v{ver}/tables/{sp}/(默认口径, 即正式产物);
变体     = results/v{ver}/sensitivity/<变体>_{sp}/tables/{sp}/(由 run_all.sh sensitivity 产出)。

产物(全部写在 results/v{ver}/sensitivity/ 下, 不动正式产物):
  COMPARISON.md          人读比较报告(指标表 + 逐变体判读 + 稳健性分级)
  comparison_summary.tsv 机读逐变体指标
  sig_flips.tsv          显著性翻转明细(gene × group × 变体)

指标: 可行 组×bin 数 / bin 数 / 显著趋势行(q_bin<0.05)与升降数 / temporal_pattern 分布 /
     specific 名单 / emergence 可判数 / 三个趋势锚点 rho /
     **相对基准的显著性保留率、方向一致率、新增与丢失行数**(稳健性的核心量)。
稳健性分级(仅规则化提示, 最终判读由人做): 保留率 ≥0.9 稳健; 0.7–0.9 部分敏感; <0.7 高度口径依赖。

用法: /home/victor/.conda/envs/omics/bin/python -I scripts/a5_sensitivity_compare.py <root>
"""
from __future__ import annotations

import glob
import os
import re
import sys
from datetime import datetime

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atx import config as atxconfig                     # noqa: E402

ANCHORS = [("ITPR1", "Purkinje", "+"), ("MBP", "OligoLineage", "+"), ("ATOH1", "GranuleLineage", "-")]
SIG_Q = 0.05


def _md(df: pd.DataFrame, floatfmt: str | None = None) -> str:
    d = df.copy()
    if floatfmt:
        for c in d.columns:
            if pd.api.types.is_numeric_dtype(d[c]):
                d[c] = d[c].map(lambda v: "" if pd.isna(v) else format(float(v), floatfmt))
    cols = [str(c) for c in d.columns]
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for r in d.astype(object).where(d.notna(), None).to_numpy():
        out.append("| " + " | ".join("" if v is None else str(v).replace("|", "/") for v in r) + " |")
    return "\n".join(out)


def metrics(tabdir: str, interdir: str, logsdir: str, sp: str) -> dict:
    """一套表的口径指标。缺文件 → None(变体跑挂时不崩)。"""
    def rd(name):
        p = os.path.join(tabdir, name)
        return pd.read_csv(p, sep="\t") if os.path.exists(p) else None

    tr, tp, sp_g, em = rd("trend_by_group.tsv"), rd("temporal_patterns.tsv"), \
        rd("specificity_global.tsv"), rd("emergence.tsv")
    if tr is None:
        return {}
    m: dict = {}
    cb = os.path.join(interdir, "cube_bin.tsv")
    if os.path.exists(cb):
        b = pd.read_csv(cb, sep="\t")
        m["n_group_bin_rows"] = int(len(b))
        m["n_feasible"] = int(b["feasible"].sum())
        m["n_bins"] = int(b["bin"].nunique())
        m["n_single_donor_feasible"] = int((b["feasible"] & b["single_donor"].astype(bool)).sum())
    sig = tr[(tr["q_bin"] < SIG_Q) & (tr["status"] == "ok")]
    m["n_gene_group_rows"] = int(len(tr))
    m["n_sig"] = int(len(sig))
    m["n_sig_up"] = int((sig["rho_bin"] > 0).sum())
    m["n_sig_down"] = int((sig["rho_bin"] < 0).sum())
    m["sig_set"] = {(r.gene, r.group): float(r.rho_bin) for r in sig.itertuples()}
    if tp is not None:
        vc = tp["temporal_pattern"].astype(str).value_counts()
        for k in ("flat", "monotonic_up", "monotonic_down", "window_high", "endpoint_high",
                  "transient_peak", "NA"):
            m[f"pat_{k}"] = int(vc.get(k, 0))
    if sp_g is not None:
        m["n_specific"] = int((sp_g["specific_call"] == "specific").sum())
        m["n_not_expressed"] = int((sp_g["specific_call"] == "not_expressed").sum())
        m["specific_set"] = {r.gene: str(r.max_group) for r in sp_g.itertuples()
                             if r.specific_call == "specific"}
    if em is not None:
        m["n_emergence"] = int(em["emergence_bin"].notna().sum())
    # 锚点: MBP/ATOH1 是 sanity marker, **不在统计 panel** → trend_by_group.tsv 里没有它们的行。
    # 正式锚点值由 s02 单独计算并写入 logs/s02_report_<sp>.tsv(形如 "rho_bin=+0.964 q_bin=..."),
    # 故优先从该日志解析; panel 内的锚点基因(如 ITPR1)在日志缺失时回退 trend 表。
    for gene, group, _ in ANCHORS:
        m[f"anchor_{gene}"] = None
    # 必须按物种取: results/v{ver}/logs/ 下同时有 s02_report_human.tsv 与 s02_report_mouse.tsv,
    # 用 glob 排序取第一个会把鼠基准读成人值(实测踩坑)。
    rep_path = os.path.join(logsdir, f"s02_report_{sp}.tsv")
    if os.path.exists(rep_path):
        rep = pd.read_csv(rep_path, sep="\t")
        for r in rep.itertuples():
            chk = str(getattr(r, "check", "") or "")
            if "趋势锚点" not in chk:
                continue
            obs = str(getattr(r, "observed", "") or "")
            mm = re.search(r"rho_bin=([+-]?[0-9.]+)", obs)
            if not mm:
                continue
            for gene, group, _ in ANCHORS:
                if f"{gene}×{group}" in chk:
                    m[f"anchor_{gene}"] = round(float(mm.group(1)), 3)
    for gene, group, _ in ANCHORS:
        if m[f"anchor_{gene}"] is None:
            row = tr[(tr["gene"] == gene) & (tr["group"] == group)]
            if len(row) and pd.notna(row.iloc[0]["rho_bin"]):
                m[f"anchor_{gene}"] = round(float(row.iloc[0]["rho_bin"]), 3)
    return m


def grade(retention: float) -> str:
    if retention >= 0.9:
        return "稳健"
    if retention >= 0.7:
        return "部分敏感"
    return "高度口径依赖"


def main() -> int:
    root = os.path.abspath(sys.argv[1])
    cfg = atxconfig.load(os.path.join(root, "config", "config.yaml"))
    ver = cfg.version
    sens_root = os.path.join(root, "results", f"v{ver}", "sensitivity")
    base = {sp: metrics(os.path.join(root, "results", f"v{ver}", "tables", sp),
                        os.path.join(root, "results", f"v{ver}", "intermediate", sp),
                        os.path.join(root, "results", f"v{ver}", "logs"), sp)
            for sp in ("human", "mouse")}

    rows, flips, notes = [], [], []
    sig_sets: dict[str, dict[str, dict]] = {"human": {}, "mouse": {}}   # sp → 变体 → sig_set
    for sp in ("human", "mouse"):
        if not base[sp]:
            continue
        b = base[sp]
        sig_sets[sp]["默认口径(基准)"] = b["sig_set"]
        rows.append({"变体": "默认口径(基准)", "物种": sp, **{k: v for k, v in b.items()
                                                             if not k.endswith("_set")}})
        for d in sorted(glob.glob(os.path.join(sens_root, f"*_{sp}"))):
            name = os.path.basename(d)[:-(len(sp) + 1)]
            m = metrics(os.path.join(d, "tables", sp), os.path.join(d, "intermediate", sp),
                        os.path.join(d, "logs"), sp)
            if not m:
                notes.append(f"- ⚠ 变体 `{name}`({sp}) 无表产物 —— 该变体未跑成功, 不参与比较。")
                continue
            sig_sets[sp][name] = m["sig_set"]
            bs, vs = b["sig_set"], m["sig_set"]
            kept = set(bs) & set(vs)
            same_dir = sum(1 for k in kept if (bs[k] > 0) == (vs[k] > 0))
            ret = len(kept) / len(bs) if bs else float("nan")
            m["保留率"] = round(ret, 3)
            m["方向一致率"] = round(same_dir / len(kept), 3) if kept else float("nan")
            m["新增显著"] = len(set(vs) - set(bs))
            m["丢失显著"] = len(set(bs) - set(vs))
            m["稳健性"] = grade(ret)
            rows.append({"变体": name, "物种": sp, **{k: v for k, v in m.items()
                                                      if not k.endswith("_set")}})
            for gene, group in sorted(set(bs) | set(vs)):
                st = "both" if gene in bs and gene in vs else ("lost" if gene in bs else "new")
                if st != "both" or (bs[gene] > 0) != (vs[gene] > 0):
                    flips.append({"species": sp, "variant": name, "gene": gene, "group": group,
                                  "change": st, "rho_baseline": bs.get(gene), "rho_variant": vs.get(gene)})
            # 特异性名单变化
            bspec, vspec = b.get("specific_set", {}), m.get("specific_set", {})
            lost_s = sorted(set(bspec) - set(vspec))
            new_s = sorted(set(vspec) - set(bspec))
            moved = sorted(g for g in set(bspec) & set(vspec) if bspec[g] != vspec[g])
            if lost_s or new_s or moved:
                notes.append(f"- `{name}`({sp}) specific 名单变化: 失去 {lost_s or '无'}; "
                             f"新增 {new_s or '无'}; 主组改变 {moved or '无'}。")
            notes.append(f"- `{name}`({sp}): 显著行 {b['n_sig']}→{m['n_sig']}, 保留率 {ret:.1%}, "
                         f"方向一致率 {m['方向一致率']:.1%}, 新增 {m['新增显著']}/丢失 {m['丢失显著']} "
                         f"→ **{m['稳健性']}**; 锚点 ITPR1 {b['anchor_ITPR1']}→{m['anchor_ITPR1']}, "
                         f"MBP {b['anchor_MBP']}→{m['anchor_MBP']}, "
                         f"ATOH1 {b['anchor_ATOH1']}→{m['anchor_ATOH1']}。")

    summ = pd.DataFrame(rows)
    os.makedirs(sens_root, exist_ok=True)
    summ.to_csv(os.path.join(sens_root, "comparison_summary.tsv"), sep="\t", index=False)
    pd.DataFrame(flips).to_csv(os.path.join(sens_root, "sig_flips.tsv"), sep="\t", index=False)

    # ---- 稳健核心集: 基准 + 全部变体都显著、且方向在所有口径下一致的 gene×group 行 ----
    # 这是唯一可以写成"断言"的集合; 其余行须表述为"口径依赖"或不写。
    gt = pd.read_csv(cfg.paths["genes_table"], sep="\t").set_index("hgnc_symbol")
    core_rows, core_notes = [], []
    for sp, variants in sig_sets.items():
        if len(variants) < 2:
            continue
        n_base = len(variants["默认口径(基准)"])
        inter = set.intersection(*[set(v) for v in variants.values()])
        consistent = {k for k in inter
                      if len({(v[k] > 0) for v in variants.values()}) == 1}
        genes: dict[str, list] = {}
        for gene, group in consistent:
            genes.setdefault(gene, []).append(group)
        core_notes.append(
            f"- **{sp}**：{len(variants)} 套口径（含基准）下都显著的 gene×group 行 = {len(inter)}；"
            f"其中方向也完全一致 = **{len(consistent)} 行 / {len(genes)} 个基因**"
            f"（基准共 {n_base} 行显著 → {len(consistent) / max(1, n_base):.0%} 可作跨口径断言）。")
        for gene, groups in sorted(genes.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            signs = {(v[(gene, g)] > 0) for v in variants.values() for g in groups}
            core_rows.append({
                "species": sp, "gene": gene, "n_groups": len(groups),
                "groups": ",".join(sorted(groups)),
                "direction": "up" if signs == {True} else ("down" if signs == {False} else "mixed"),
                "mech_class": gt.loc[gene, "mech_class"] if gene in gt.index else "NA",
                "tier": gt.loc[gene, "tier"] if gene in gt.index else "NA"})
    core_tab = pd.DataFrame(core_rows)
    if len(core_tab):
        core_tab.to_csv(os.path.join(sens_root, "robust_core_genes.tsv"), sep="\t", index=False)

    show = ["变体", "物种", "n_bins", "n_feasible", "n_single_donor_feasible", "n_sig",
            "n_sig_up", "n_sig_down", "pat_monotonic_up", "pat_monotonic_down", "pat_window_high",
            "pat_endpoint_high", "pat_transient_peak", "pat_flat", "n_specific", "n_emergence",
            "保留率", "方向一致率", "新增显著", "丢失显著", "稳健性"]
    L = [f"# 敏感性变体比较 — gene_set v{ver}",
         f"> 生成 {datetime.now():%Y-%m-%d %H:%M}；基准 = `results/v{ver}/tables/<物种>/`（默认口径，正式产物）；"
         f"变体 = `results/v{ver}/sensitivity/<变体>_<物种>/`。显著性判据 q_bin<{SIG_Q}（BH，族=species×group）。",
         "",
         "## 1. 指标总表",
         _md(summ[[c for c in show if c in summ.columns]]),
         "",
         "## 2. 逐变体判读（相对基准）",
         *notes,
         "",
         "## 3. 稳健核心集（基准 + 全部变体都显著且方向一致 → 可写成断言）",
         *core_notes,
         "",
         _md(core_tab.head(60)) if len(core_tab) else "（无）",
         "",
         f"完整清单见 `robust_core_genes.tsv`（{len(core_tab)} 行）。",
         "",
         "## 4. 三个趋势锚点跨变体对照",
         _md(summ[["变体", "物种", "anchor_ITPR1", "anchor_MBP", "anchor_ATOH1"]]),
         "",
         "## 5. 怎么读这份表",
         "- **保留率** = 基准里显著的 gene×group 行，在该变体下仍显著的比例。≥90% 记「稳健」，"
         "70–90%「部分敏感」，<70%「高度口径依赖」。",
         "- **方向一致率** = 保留下来的行里升降方向不变的比例。方向翻转比"
         "「变得不显著」更严重（前者是结论反了，后者只是功效不够）。",
         "- **新增显著** 通常来自口径变化后 bin 数/细胞数改变，不一定是新生物学；须逐条看 "
         "`sig_flips.tsv`。",
         "- 锚点（ITPR1/MBP/ATOH1）若在某个变体下方向翻转 → 该变体口径不可用（属结构性问题，"
         "不是功效问题）。",
         "- 写报告时的用法：主分析用默认口径；只把在多口径下**方向一致且保留率高**的结论写成断言，"
         "其余降级为「口径依赖」或不写。",
         ""]
    with open(os.path.join(sens_root, "COMPARISON.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    print(f"[a5] 变体数 {len(rows) - sum(1 for r in rows if r['变体'].startswith('默认'))};"
          f" 翻转明细 {len(flips)} 行")
    print(f"[a5] 报告 → {os.path.join(sens_root, 'COMPARISON.md')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
