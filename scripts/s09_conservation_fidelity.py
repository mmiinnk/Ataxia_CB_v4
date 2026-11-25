#!/usr/bin/env python3
"""s09_conservation_fidelity.py — 第 5 层后半:人鼠保守性四指标 + 加权保真分(spec §8.1)。

输入: 两物种 tables/{sp}/{trend_by_group, specificity_global, peak_by_group}.tsv;
     s04 的 tables/cross/state_profiles.tsv; genes/ataxia_gene_disease.tsv(疾病保真分行来源)。
指标(每 gene × group, 两物种均可行才判; 缺证据 → **NA 不当 0**):
  rho_sign_concordant  两侧 |rho_bin| ≥ 0.2 时判符号是否一致(1/0)
  traj_r               组内共享 dev_state 子集上的 Pearson r(≥3 个共享可行状态)
  max_group_match      两物种 global max_group(τ argmax)是否一致(gene 级)
  peak_state_match     各自峰 dev_state 在**共识序**中相差 ≤ peak_state_tolerance(1/0)
保真分 = 非 NA 指标加权均值(权重 config stats.conservation.*_weight, NA 项按剩余权重重归一);
n_metrics = 非 NA 指标数; 证据等级 = ≥3 high / 2 medium / ≤1 low(spec §8.1)。
**口径注**: traj_r 按原值(−1..1)入权(不做改标度), 故 fidelity 理论区间 [−0.3, 1.0]。

产物: tables/cross/{conservation_gene_by_group, conservation_gene, fidelity_by_disease}.tsv
     figures/cross/fig5e_fidelity_table.{pdf,png}   (fig5a–d = s04 的 x1–x4)
     logs/s09_report.tsv; reviews/gate3_cross.md(人工 review 门 3 清单, 自动项 + 四条勾选)
用法: <omics python> -I scripts/s09_conservation_fidelity.py <root>
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atx import align as atxalign                        # noqa: E402
from atx import config as atxconfig                      # noqa: E402

plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC", "Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

REPORT: list[dict] = []
N_FAIL = 0


def rec(check, expected, observed, level, note=""):
    global N_FAIL
    REPORT.append({"check": check, "expected": str(expected), "observed": str(observed),
                   "level": level, "note": note})
    mark = {"PASS": "✓", "WARN": "⚠", "FAIL": "✗"}.get(level, "?")
    print(f"[{mark} {level}] {check} | 期望={expected} 实测={observed}" + (f" | {note}" if note else ""))
    N_FAIL += level == "FAIL"


def _finish(cfg, t0):
    logs = cfg.results("logs", mkdir=True)
    pd.DataFrame(REPORT).to_csv(os.path.join(logs, "s09_report.tsv"), sep="\t", index=False)
    n_pass = sum(1 for r in REPORT if r["level"] == "PASS")
    n_warn = sum(1 for r in REPORT if r["level"] == "WARN")
    print(f"\n汇总: PASS={n_pass} WARN={n_warn} FAIL={N_FAIL}  "
          f"报告: {os.path.join(logs, 's09_report.tsv')}  ({time.time() - t0:.0f}s)")


def load_profiles(cfg) -> dict[str, dict[str, pd.DataFrame]]:
    """s04 的 state_profiles.tsv(长表) → {species: {scope: 基因×状态宽表}}; 列序按 rank_in_species。"""
    p = cfg.results("tables", "cross", "state_profiles.tsv")
    if not os.path.exists(p):
        raise SystemExit(f"[s09] 缺 {p} —— 先跑 s04_figures_crossspecies.py")
    df = pd.read_csv(p, sep="\t")
    out: dict[str, dict[str, pd.DataFrame]] = {}
    for sp, dsp in df.groupby("species"):
        out[sp] = {}
        for scope, dsc in dsp.groupby("scope"):
            dsc = dsc.sort_values(["gene", "rank_in_species"])
            wide = dsc.pivot(index="gene", columns="dev_state", values="mean_log")
            order = (dsc.drop_duplicates("dev_state").sort_values("rank_in_species")["dev_state"].tolist())
            out[sp][scope] = wide.reindex(columns=order)
    return out


def write_gate3(cfg, root, fid_gene: pd.DataFrame, by_disease: pd.DataFrame) -> str:
    """生成人工 review 门 3 清单(spec §9: run_cross 之后; 焦点 = fig3a/fig4a 统计)。

    与 s06 同语义: 目标文件已存在(= 人工已签署/在签) → 只写 .refresh, **绝不覆盖人工留痕**(硬约束 9)。
    """
    logs = cfg.results("logs")
    frames = []
    for f in ("s04_report.tsv", "s07_report.tsv", "s09_report.tsv"):
        p = os.path.join(logs, f)
        if os.path.exists(p):
            d = pd.read_csv(p, sep="\t")
            d.insert(0, "source", f.replace(".tsv", ""))
            frames.append(d)
    allr = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    fails = allr[allr["level"] == "FAIL"] if len(allr) else allr
    warns = allr[allr["level"] == "WARN"] if len(allr) else allr

    def _tbl(d, cols):
        if not len(d):
            return "| (无) |  |  |\n|---|---|---|"
        head = "| " + " | ".join(cols) + " |\n|" + "|".join(["---"] * len(cols)) + "|"
        body = ["| " + " | ".join("" if pd.isna(r.get(c)) else str(r.get(c))[:150].replace("|", "/")
                                 for c in cols) + " |" for _, r in d.iterrows()]
        return "\n".join([head] + body)

    enr_p = cfg.results("tables", "cross", "enrichment_class_x_module.tsv")
    enr = pd.read_csv(enr_p, sep="\t") if os.path.exists(enr_p) else pd.DataFrame()
    sig = enr[(enr["q"] < 0.05) & (enr["obs"] >= 3)] if len(enr) else enr
    L = ["# review 门 3 — cross/spec §9(s04+s07+s09 自动项; 人工项勾选后 commit)", "",
         f"自动项总状态: {'通过' if not len(fails) else '**有 FAIL**'}"
         f"(FAIL {len(fails)} / WARN {len(warns)} / 日志行 {len(allr)})", "",
         "## (1) logs 无 assert 失败且 WARNING 已逐条看过",
         f"- 自动: {'✓ 通过' if not len(fails) else '✗ 有 FAIL'} — FAIL {len(fails)} 条",
         f"- WARN {len(warns)} 条(判读后勾选):",
         _tbl(warns, ["source", "check", "observed", "note"]), "",
         "- [ ] 人工确认: WARNING 已逐条看过且非阻断", "",
         "## (2) 跨物种层口径核查(锚点表的跨物种等价项)",
         "- 自动项(逐条见 logs/s04_report.tsv、s07_report.tsv):",
         _tbl(allr[allr["source"].isin(["s04_report", "s07_report"])],
              ["source", "check", "expected", "observed", "level"]), "",
         "- [ ] 人工确认: dev_state 对齐口径(§4.4/§4.5A)、CL-027 残差类排除、NA 政策均无异议", "",
         "## (3) marker 校验",
         "- 自动: 跨物种层不重复 marker 校验(marker 属物种内 check, 已在门 1/门 2 判读)。",
         "- [ ] 人工确认: 认可本层免检(如需, 回门 1/门 2 复核)", "",
         "## (4) fig3a 富集与 fig5 保真度统计(门 3 焦点)",
         f"- fig3a 显著富集格(q<0.05 且 obs≥3): **{len(sig)}** 格"
         + (f" — 最强: {sig.iloc[0]['species']} {sig.iloc[0]['a']}×{sig.iloc[0]['b']} "
            f"log2FE={sig.iloc[0]['log2fe']:.2f} q={sig.iloc[0]['q']:.3g}" if len(sig) else
            "(无显著格 = 信息性结果, 不得解读为「无共性机制」)"),
         f"- 保真度: 基因级 {len(fid_gene)} 行; 证据等级 "
         f"{dict(fid_gene['evidence_grade'].value_counts())}; 保真分中位 "
         f"{fid_gene['fidelity'].median():+.3f}(区间 [−0.3, 1.0])",
         f"- 疾病级 {len(by_disease)} 行; 中位保真分 {by_disease['fidelity'].median():+.3f}"
         if len(by_disease) else "- 疾病级: 无",
         "- [ ] 人工确认: fig3a(富集热图)/fig3c(三套峰轴)/fig5a–e 的统计与图面一致、结论表述与证据等级相符",
         "", "---",
         f"reviewer: ____________  日期: ____________  结论: [ ] 过门 / [ ] 整改后复检",
         f"(自动生成 {datetime.now():%Y-%m-%d %H:%M}; gene_set v{cfg.version}; 产物 results/v{cfg.version}/)",
         ""]
    target = os.path.join(root, "reviews", "gate3_cross.md")
    out = target if not os.path.exists(target) else target.replace(".md", ".refresh.md")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--config", default=None)
    ap.add_argument("--out-subdir", default=None)
    ap.add_argument("--sens", action="append", default=None, metavar="KEY=VALUE")
    args = ap.parse_args()
    root = os.path.abspath(args.root)
    t0 = time.time()
    cfg = atxconfig.load(args.config or os.path.join(root, "config", "config.yaml"))
    atxconfig.apply_sens_overrides(cfg, args.sens)
    if args.out_subdir:
        cfg.results_subdir = args.out_subdir.strip("/")
    cfg.check_version_gate()
    print(f"== s09 cross | version={cfg.version}(draft={cfg.draft})"
          + (f" | 产物隔离: results/v{cfg.version}/{cfg.results_subdir}/" if cfg.results_subdir else ""))

    prof = load_profiles(cfg)
    tabs = {}
    for sp in ("human", "mouse"):
        td = cfg.results("tables", sp)
        tabs[sp] = {
            "trend": pd.read_csv(os.path.join(td, "trend_by_group.tsv"), sep="\t"),
            # 与 s03.read_tables 同口径: 保留受控 "NA" 串(空串→NaN), 否则 max_group "NA"
            # 被解析成 NaN→"nan", 在 align.fidelity 中伪造 max_group_match=0(缺证据当不保守)。
            "spec": pd.read_csv(os.path.join(td, "specificity_global.tsv"), sep="\t",
                                keep_default_na=False, na_values=[""]),
            "peak": pd.read_csv(os.path.join(td, "peak_by_group.tsv"), sep="\t"),
        }
    fid = atxalign.fidelity(tabs["human"]["trend"], tabs["mouse"]["trend"],
                            tabs["human"]["spec"], tabs["mouse"]["spec"],
                            prof["human"], prof["mouse"], cfg)
    rec("保真度计算", "每 gene×group 四指标 + 加权分 + n_metrics + 证据等级",
        f"{len(fid)} 行(gene×group); 证据等级 {dict(fid['evidence_grade'].value_counts())}",
        "PASS" if len(fid) else "FAIL")

    # ---- 基因级: 取该基因的主组行(人 max_group 优先, 不可判时取指标最全的行) ----
    mg = tabs["human"]["spec"].set_index("gene")["max_group"].astype(str)
    rows = []
    for gene, g in fid.groupby("gene"):
        want = mg.get(gene, "NA")
        cand = g[g["group"] == want]
        if not len(cand) or cand.iloc[0]["n_metrics"] == 0:
            cand = g.sort_values(["n_metrics", "fidelity"], ascending=False, na_position="last").head(1)
        r = cand.iloc[0].to_dict()
        r["main_group_used"] = r["group"]
        r["main_group_is_human_max_group"] = bool(r["group"] == want)
        rows.append(r)
    fid_gene = pd.DataFrame(rows).sort_values("fidelity", ascending=False, na_position="last")
    fid_gene = fid_gene.merge(
        pd.read_csv(cfg.paths["genes_table"], sep="\t")[
            ["hgnc_symbol", "mech_class", "tier", "disease_primary", "onset_class"]],
        left_on="gene", right_on="hgnc_symbol", how="left").drop(columns=["hgnc_symbol"])

    # ---- 疾病级(spec §8.1: 疾病保真分 = 该疾病基因的均值) ----
    gd = pd.read_csv(os.path.join(root, "genes", "ataxia_gene_disease.tsv"), sep="\t", dtype=str)
    gd = gd[["hgnc_symbol", "disease_abbrev", "disease_name", "clinical_group_draft", "tier"]]
    merged = gd.merge(fid_gene[["gene", "fidelity", "n_metrics", "evidence_grade", "traj_r",
                                "rho_sign_concordant", "max_group_match", "peak_state_match"]],
                      left_on="hgnc_symbol", right_on="gene", how="left")
    by_dis = merged.groupby(["disease_abbrev"], dropna=False).agg(
        disease_name=("disease_name", "first"),
        clinical_group_draft=("clinical_group_draft", "first"),
        n_genes=("hgnc_symbol", "nunique"),
        n_genes_scored=("fidelity", lambda s: int(s.notna().sum())),
        fidelity=("fidelity", "mean"),
        traj_r=("traj_r", "mean"),
        rho_sign_concordant=("rho_sign_concordant", "mean"),
        max_group_match=("max_group_match", "mean"),
        peak_state_match=("peak_state_match", "mean"),
        n_metrics=("n_metrics", "mean"),
        genes=("hgnc_symbol", lambda s: ",".join(sorted(set(s)))),
    ).reset_index().sort_values("fidelity", ascending=False, na_position="last")

    crossdir = cfg.results("tables", "cross", mkdir=True)
    fid.to_csv(os.path.join(crossdir, "conservation_gene_by_group.tsv"), sep="\t", index=False,
               float_format="%.6g")
    fid_gene.to_csv(os.path.join(crossdir, "conservation_gene.tsv"), sep="\t", index=False,
                    float_format="%.6g")
    by_dis.to_csv(os.path.join(crossdir, "fidelity_by_disease.tsv"), sep="\t", index=False,
                  float_format="%.6g")

    grade = fid_gene["evidence_grade"].value_counts()
    rec("证据等级分布", "n_metrics ≥3 high / 2 medium / ≤1 low(spec §8.1)",
        f"high {int(grade.get('high', 0))} / medium {int(grade.get('medium', 0))} / "
        f"low {int(grade.get('low', 0))} / none {int(grade.get('none', 0))}",
        "PASS" if int(grade.get("high", 0)) else "WARN",
        note="low/none 的基因不得单独作为保守性结论(证据不足, 非阴性结果)")
    hi = fid_gene[fid_gene["evidence_grade"] == "high"]
    if len(hi):
        rec("高证据基因保真分", "描述统计",
            f"n={len(hi)}, 中位 {hi['fidelity'].median():+.3f}, "
            f"top5 {list(hi.nlargest(5, 'fidelity')['gene'])}", "PASS")
    rec("四指标可算率", "NA 不当 0",
        "; ".join(f"{c}={int(fid_gene[c].notna().sum())}/{len(fid_gene)}"
                  for c in ("rho_sign_concordant", "traj_r", "max_group_match", "peak_state_match")),
        "PASS")

    # ---- fig5e: 保真度表热图(行=疾病 top/bottom, 列=四指标+综合分) ----
    figdir = cfg.results("figures", "cross", mkdir=True)
    exts = list(cfg.raw["figures"]["formats"])
    dpi = int(cfg.raw["figures"]["dpi"])
    show = pd.concat([by_dis.head(15), by_dis.tail(10)]).drop_duplicates("disease_abbrev")
    show = show[show["fidelity"].notna()]
    cols = ["rho_sign_concordant", "traj_r", "max_group_match", "peak_state_match", "fidelity"]
    M = show[cols].to_numpy(float)
    fig, ax = plt.subplots(figsize=(7.6, max(3.0, 0.30 * len(show))))
    cmap = plt.get_cmap(cfg.raw["figures"]["diverging_cmap"]).copy()
    cmap.set_bad(color=cfg.raw["figures"]["na_style"]["facecolor"])
    im = ax.imshow(np.ma.masked_invalid(M), aspect="auto", cmap=cmap, vmin=-0.3, vmax=1.0)
    ax.set_xticks(range(len(cols)),
                  ["rho 同向", "轨迹 r", "τ 主组一致", "峰状态±1", "综合保真分"], fontsize=7.5)
    ax.set_yticks(range(len(show)),
                  [f"{r.disease_abbrev} (n={r.n_genes_scored}/{r.n_genes})" for r in show.itertuples()],
                  fontsize=6.5)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center", fontsize=5.4,
                        color="w" if M[i, j] > 0.65 else "k")
    ax.set_title("fig5e 疾病级保守性/保真度表(前 15 + 后 10; 白格 = 证据不足 NA)\n"
                 "综合分区间 [−0.3, 1.0]; 权重 config stats.conservation(0.3/0.3/0.2/0.2)", fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.8, label="value")
    fig.tight_layout()
    for e in exts:
        fig.savefig(os.path.join(figdir, f"fig5e_fidelity_table.{e}"), dpi=dpi)
    plt.close(fig)

    rec("产物", "tables/cross + figures/cross",
        f"conservation_gene_by_group({len(fid)}), conservation_gene({len(fid_gene)}), "
        f"fidelity_by_disease({len(by_dis)}); fig5e ×{len(exts)} 格式"
        "(fig5a–d = s04 的 x1–x4)", "PASS")
    g3 = write_gate3(cfg, root, fid_gene, by_dis)
    rec("门 3 清单", "自动项 + 四条人工勾选(硬约束 9: 只能人工签)", f"→ {g3}", "PASS")
    _finish(cfg, t0)
    return 1 if N_FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
