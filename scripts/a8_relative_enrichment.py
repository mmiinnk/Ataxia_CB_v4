"""a8 — M1 相对富集/剂量阈值假说检验（核心假说 docs/核心假说_细胞选择性易感.md §3 M1）。

问题：PC 相对表达占比在成熟期升高的基因，是否更可能 PC 选择性病理（PC_loss）与较晚起病？
数据：鼠 specificity_by_period（人 PC 仅产前，无出生后 PC，故用鼠）+ 基因临床注释。

口径要点（NA 政策硬约束 4；单一事实来源硬约束 6）：
  · 跨期 PC 占比直接用 fraction_Purkinje 会被分母变化混杂（鼠 embryonic 7 组 / postnatal 8 组，
    后者多 OligoLineage）。fraction_8group 严格版在 embryonic 全 NA（早期不足 8 组）不可用。
    → 用 cp10k 列在**两期共有的 7 组**上重算 PC 占比（共同分母），ΔPC = postnatal − embryonic。
  · 另报 postnatal 8 组 fraction_Purkinje（期内全基因同分母，跨基因可比）作横断面指标。
  · PC_loss 来自临床注释 primary_pathology（G3）；onset 关联剔除 NA（CL-024 未核实，初步）。
  · 检验：PC_loss vs 非 PC_loss 用 group_compare（MW + Cliff's δ）；onset 用 mech_class 分层置换
    Spearman（atx.stats_enrich，复用 s08 口径）。关联非因果，属 M1 假说检验。

产物：results/v{ver}/analysis/a8_relative_enrichment/{a8_per_gene.tsv, a8_tests.tsv, a8_relative_enrichment.md}
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from atx import config as atxconfig            # noqa: E402
from atx import stats_enrich as atxenrich      # noqa: E402


def common_groups(bp: pd.DataFrame) -> list[str]:
    """两期共有可行组（embryonic ∩ postnatal 的 groups_used），含 Purkinje。"""
    def groups_of(period):
        sub = bp[bp.period == period]
        gu = sub.groups_used.dropna()
        return set(gu.iloc[0].split("|")) if len(gu) else set()
    common = groups_of("embryonic") & groups_of("postnatal")
    assert "Purkinje" in common, "共有组须含 Purkinje"
    return sorted(common)


def pc_fraction_common(bp: pd.DataFrame, gene: str, period: str, groups: list[str]) -> float:
    row = bp[(bp.gene == gene) & (bp.period == period)]
    if not len(row):
        return float("nan")
    row = row.iloc[0]
    vals = [row.get(f"cp10k_{g}", np.nan) for g in groups]
    if any(pd.isna(v) for v in vals):
        return float("nan")
    den = float(sum(vals))
    return float(row["cp10k_Purkinje"]) / den if den > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--n-perm", type=int, default=10000)
    args = ap.parse_args()
    root = os.path.abspath(args.root)
    cfg = atxconfig.load(os.path.join(root, "config", "config.yaml"))
    cfg.check_version_gate()
    seed = cfg.seed

    gt = pd.read_csv(cfg.paths["genes_table"], sep="\t")
    panel = gt[(gt.tier == "core") & gt.in_panel.astype(str).str.lower().eq("true")].copy()
    bp = pd.read_csv(cfg.results("tables", "mouse", "specificity_by_period.tsv"), sep="\t")
    groups = common_groups(bp)
    print(f"[a8] 共同分母组({len(groups)}): {groups}")

    post = bp[bp.period == "postnatal"].set_index("gene")
    recs = []
    for r in panel.itertuples():
        g = r.hgnc_symbol
        emb = pc_fraction_common(bp, g, "embryonic", groups)
        pst = pc_fraction_common(bp, g, "postnatal", groups)
        delta = pst - emb if np.isfinite(emb) and np.isfinite(pst) else float("nan")
        post8 = float(post.loc[g, "fraction_Purkinje"]) if g in post.index else float("nan")
        recs.append({"gene": g, "mech_class": r.mech_class,
                     "PC_loss": bool("PC_loss" in str(r.primary_pathology)),
                     "onset_class": r.onset_class,
                     "onset_typical_y": pd.to_numeric(r.onset_typical_y, errors="coerce"),
                     "pcfrac_emb_c7": emb, "pcfrac_post_c7": pst,
                     "delta_pc_c7": delta, "pcfrac_post_8g": post8})
    pg = pd.DataFrame(recs)

    # ---- 检验 1: PC_loss vs 非 PC_loss（MW + Cliff's δ）----
    tests = []
    pcloss = pg.PC_loss.map({True: "PC_loss", False: "non_PC_loss"})
    for metric in ("delta_pc_c7", "pcfrac_post_8g"):
        kw_p, pair = atxenrich.group_compare(pg[metric].to_numpy(), pcloss.to_numpy())
        if len(pair):
            pr = pair.iloc[0]
            tests.append({"test": "PC_loss_vs_non", "metric": metric,
                          "n_PC_loss": int((pcloss == "PC_loss").sum()),
                          "n_non": int((pcloss == "non_PC_loss").sum()),
                          "n_used": int(pg[metric].notna().sum()),
                          "p_mw": pr.p_mw, "cliffs_delta": pr.cliffs_delta})
    # ---- 检验 2: onset 分层置换 Spearman（mech_class 分层）----
    y = np.log1p(pg.onset_typical_y.to_numpy(dtype=float))
    for metric in ("delta_pc_c7", "pcfrac_post_8g"):
        rho, p_strat, p_naive = atxenrich.stratified_perm_spearman(
            pg[metric].to_numpy(dtype=float), y, pg.mech_class.to_numpy(), args.n_perm, seed=seed)
        tests.append({"test": "onset_stratified_spearman", "metric": metric,
                      "n_PC_loss": "", "n_non": "", "n_used": int((np.isfinite(y) & pg[metric].notna()).sum()),
                      "rho": rho, "p_strat": p_strat, "p_naive": p_naive})
    tt = pd.DataFrame(tests)

    outdir = cfg.results("analysis", "a8_relative_enrichment", mkdir=True)
    pg.to_csv(os.path.join(outdir, "a8_per_gene.tsv"), sep="\t", index=False)
    tt.to_csv(os.path.join(outdir, "a8_tests.tsv"), sep="\t", index=False)

    # ---- 摘要 md ----
    def fmt(v, nd=3):
        return "NA" if pd.isna(v) else f"{v:.{nd}f}"
    rep = pg[pg.gene.isin(["CACNA1A", "ATXN1", "ATXN3", "ATXN7", "TBP", "ITPR1"])]
    L = ["# a8 — M1 相对富集/剂量阈值假说检验（鼠）", "",
         f"gene_set v{cfg.version}；共同分母 {len(groups)} 组：{', '.join(groups)}。",
         "PC 占比在两期共有组上重算（去除 postnatal 多 OligoLineage 的分母混杂）；ΔPC=postnatal−embryonic。",
         "另报 postnatal 8 组 fraction（期内同分母，跨基因可比）。PC_loss 来自临床注释；onset 剔除 NA（CL-024）。", "",
         "## 检验结果", "",
         "| 检验 | 指标 | n | 统计量 | p | 效应 |", "|---|---|---|---|---|---|"]
    for r in tt.itertuples():
        if r.test == "PC_loss_vs_non":
            L.append(f"| PC_loss vs 非 | {r.metric} | {r.n_PC_loss}/{r.n_non} | MW | {fmt(r.p_mw)} | δ={fmt(r.cliffs_delta)} |")
        else:
            L.append(f"| onset 分层 Spearman | {r.metric} | {r.n_used} | ρ={fmt(r.rho)} | p_strat={fmt(r.p_strat)} | p_naive={fmt(r.p_naive)} |")
    L += ["", "## 代表基因（共同分母 PC 占比）", "",
          "| 基因 | 机制类 | PC_loss | embryonic | postnatal | ΔPC | postnatal 8组 |",
          "|---|---|---|---|---|---|---|"]
    for r in rep.itertuples():
        L.append(f"| {r.gene} | {r.mech_class} | {'是' if r.PC_loss else '否'} | {fmt(r.pcfrac_emb_c7)} | "
                 f"{fmt(r.pcfrac_post_c7)} | {fmt(r.delta_pc_c7)} | {fmt(r.pcfrac_post_8g)} |")
    L += ["", "## 判读边界", "",
          "- 关联非因果；属 M1 假说检验，须疾病线验证（核心假说文档 §4）。",
          "- 仅鼠（人 PC 产前-only）；PC 占比为相对量，受其他组表达影响。",
          "- onset 关联为初步（CL-024 未核实，NA 剔除）；残差类 other 自成 stratum（low n）。",
          "- ΔPC 用共同 7 组分母校正；fraction_8group 严格版 embryonic 全 NA 不可用。"]
    with open(os.path.join(outdir, "a8_relative_enrichment.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")

    print(tt.to_string(index=False))
    print(f"[a8] 产物 → {outdir}")


if __name__ == "__main__":
    main()