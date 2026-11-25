#!/usr/bin/env python
"""dev_check_human_anchors — 人趋势锚点的独立数据实测核验(dev 轨, P7 裁决②落地)。

背景: spec §9(4) 原仅设鼠趋势锚点;用户裁决(2026-11-05)人侧现在就设。文献方向
(Gold 2003 Neuron: Itpr1 属 RORα 依赖的浦肯野成熟程序;Zečević 1998/Ozarkar 2025:
MBP 随髓鞘化上升;Haldipur 2022 Annu Rev Neurosci: ATOH1 为 RL/EGL 祖细胞标记、
分化即下调)之外, 按硬约束 10 用**独立于人 Sepp 数据**的 Aldinger 2021 人胎小脑
snRNA-seq(输入四文件之一, w9–w21, 13 供体)实测方向, 避免用被验数据自证。

口径: X = log 归一化值(硬约束 4);每 (区室, 发育周) 求基因均值 → Spearman rho vs
周序;周×区室 < 30 核的点子剔除后再算 rho(NA 政策精神, 不插值)。区室定义 =
s00_label_mapping_aldinger.tsv 的 fig_cell_type→group 同源映射(注释列写明成员)。
ATOH1 附加横向对照: 同周内祖细胞区室(RL+GCP)均值应 > 已分化颗粒神经元(GN)。

用法:  python -I scripts/dev_check_human_anchors.py <root> [--config PATH]
产物:  results/v{ver}/logs/dev_human_anchor_check.tsv(DRAFT dev 轨诊断, 非正式层)
"""
from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from atx import config as atxconfig

# 锚点基因 → (Aldinger fig_cell_type 区室成员, 期望方向, 胎窗适用性) —— 区室 = s00 标签映射同源。
# window="postnatal": 方向事件发生在出生后(髓鞘化), Aldinger 胎窗(w9–21)不适用 →
# 判 N/A_WINDOW 而非 FAIL(实测 rho 仍记录在案; 文献来源见 source 与 config trend_anchors)。
ANCHORS = [
    ("ITPR1", "ENSG00000150995", ["H-PC"], "+", "Purkinje 成熟程序(Gold 2003)", "fetal"),
    ("MBP", "ENSG00000197971", ["H-OPC", "H-Committed OPC"], "+", "髓鞘化上升(Zečević 1998)", "postnatal"),
    ("ATOH1", "ENSG00000172238", ["H-RL", "H-GCP"], "-", "祖细胞标记分化下调(Haldipur 2022)", "fetal"),
]
GN_COMPARTMENT = ["H-GN"]     # ATOH1 横向对照: 已分化颗粒神经元
MIN_CELLS = 30                # 与主分析 min_cells 同值(config bins.min_cells 的 dev 复刻;dev 轨只读诊断)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    root = os.path.abspath(args.root)
    cfg = atxconfig.load(args.config or os.path.join(root, "config", "config.yaml"))

    import anndata
    path = cfg.h5ad("aldinger_h5ad")
    print(f"== 全载 Aldinger: {path}")
    ad = anndata.read_h5ad(path)
    obs = ad.obs
    # 发育周序: '9th week post-fertilization stage' → 9
    wk = obs["development_stage"].astype(str).str.extract(r"(\d+)(?:th|st|nd|rd) week")[0].astype(int)
    fct = obs["fig_cell_type"].astype(str)

    var_idx = [str(v) for v in ad.var.index]
    feat = ad.var["feature_name"].astype(str).tolist() if "feature_name" in ad.var.columns else []
    X = ad.X

    recs = []
    for gene, ensg, members, sign, src, window in ANCHORS:
        if ensg in var_idx:
            row = var_idx.index(ensg)
        elif gene in feat:
            row = feat.index(gene)
        else:
            recs.append({"gene": gene, "compartment": "|".join(members), "rho": "NA",
                         "n_stages_used": 0, "expected_sign": sign, "verdict": "ABSENT",
                         "source": src})
            continue
        sel_all = fct.isin(members)
        pts = []
        for w in sorted(wk.unique()):
            sel = sel_all & (wk == w)
            n = int(sel.sum())
            if n < MIN_CELLS:
                continue                      # 低核数周剔除(不插值)
            v = np.asarray(X[sel.to_numpy(), row].todense() if hasattr(X[sel.to_numpy(), row], "todense")
                           else X[sel.to_numpy(), row]).ravel()
            pts.append((w, n, float(np.mean(v))))
        if len(pts) < 3:
            recs.append({"gene": gene, "compartment": "|".join(members), "rho": "NA",
                         "n_stages_used": len(pts), "expected_sign": sign,
                         "verdict": "LOW_POWER", "source": src})
            continue
        ws = [p[0] for p in pts]
        ms = [p[2] for p in pts]
        rho = float(spearmanr(ws, ms).statistic)
        good = (rho > 0) if sign == "+" else (rho < 0)
        verdict = "PASS" if good else ("N/A_WINDOW" if window == "postnatal" else "FAIL")
        recs.append({"gene": gene, "compartment": "|".join(members), "rho": f"{rho:+.3f}",
                     "n_stages_used": len(pts), "expected_sign": sign,
                     "verdict": verdict,
                     "per_stage_mean": ";".join(f"w{w}(n={n}):{m:.3f}" for w, n, m in pts),
                     "source": src + (f"; {window} 事件, 胎窗{'适用' if window == 'fetal' else '不适用(记 N/A)'}"
                                      if not good or window == "postnatal" else "")})

    # ATOH1 横向对照: 每周 祖(RL+GCP) − 已分化(GN)
    if "ATOH1" in var_idx or "ATOH1" in feat:
        row = var_idx.index("ENSG00000172238") if "ENSG00000172238" in var_idx else feat.index("ATOH1")
        diffs = []
        for w in sorted(wk.unique()):
            a = (fct.isin(ANCHORS[2][2])) & (wk == w)
            b = (fct.isin(GN_COMPARTMENT)) & (wk == w)
            if int(a.sum()) < MIN_CELLS or int(b.sum()) < MIN_CELLS:
                continue
            va = np.asarray(X[a.to_numpy(), row].todense() if hasattr(X[a.to_numpy(), row], "todense")
                            else X[a.to_numpy(), row]).ravel()
            vb = np.asarray(X[b.to_numpy(), row].todense() if hasattr(X[b.to_numpy(), row], "todense")
                            else X[b.to_numpy(), row]).ravel()
            diffs.append((w, float(np.mean(va) - np.mean(vb))))
        if diffs:
            frac_pos = sum(1 for _w, d in diffs if d > 0) / len(diffs)
            recs.append({"gene": "ATOH1", "compartment": "(RL+GCP)−GN 横向对照",
                         "rho": f"{frac_pos:.2f}(阳性周占比)", "n_stages_used": len(diffs),
                         "expected_sign": ">0.5", "verdict": "PASS" if frac_pos > 0.5 else "FAIL",
                         "per_stage_mean": ";".join(f"w{w}:{d:+.3f}" for w, d in diffs),
                         "source": "祖>已分化(横向)"})

    df = pd.DataFrame(recs)
    logdir = cfg.results("logs", mkdir=True)
    out = os.path.join(logdir, "dev_human_anchor_check.tsv")
    df.to_csv(out, sep="\t", index=False)
    print(df.drop(columns=["per_stage_mean", "source"], errors="ignore").to_string(index=False))
    print(f"\n== 写入 {out}")
    n_fail = int((df["verdict"] == "FAIL").sum())
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
