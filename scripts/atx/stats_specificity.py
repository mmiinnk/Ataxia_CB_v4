"""atx.stats_specificity — type_fraction / tau / specific_call / tau×bin / emergence(spec §3.7)。

借鉴(Depression_CB README §7.3/§7.4 已验证口径):线性 CP10K、bin 等权聚合(抗发育组成
变化);tau = Σ(1 − x/x_max)/(k−1)(加 pseudocount);specific_call 规则 = tau≥tau_thr 且
fraction≥frac_thr 且最强组 CP10K≥min_expr;缺组期不隐式重归一(pandas sum skipna 陷阱 →
本实现分母显式只在可行组上求和);argmax NA 安全(只对可行组求, ties 取 config 主组序首个)。

本项目差异(§3.7):
① 分母 = 8 主组(config group_tiers.main;UBC/ALL 不入);
② 可行组数 k < min_groups(config tau_by_bin.min_groups)→ tau/fraction/max_fraction/
   specific_call 一律 NA, 输出 groups_used("|" 连接) 与 n_groups_used;数值只在
   groups_used 相同的 bin/期之间可比(Fig 1e 列头须标 n_groups);max_group 仍给可行组内
   argmax(信息列);
③ fraction_8group_<g> 严格版: 固定 8 组分母, 仅当 8 组全可行才非 NA(补充表用);
④ emergence = 首个"主组(global max_group) fraction ≥ tau_by_bin.fraction_threshold"
   的 bin, 仅在 n_groups_used ≥ min_groups 的 bin 上判定;附 first_computable_bin 与
   left_censored(emergence 落在首个可算 bin = 左截断, 只表示"不晚于该 bin 已特异");
⑤ specific_call 受控词表: specific / broad / not_expressed / NA
   (缺失基因或 k<min_groups → "NA" 而非 False——防"非特异"误读, Depression_CB §3-14)。

确定性: 无 RNG;两次重跑逐值一致(R8)。
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def aggregate(bin_table, genes, cfg, period=None) -> pd.DataFrame:
    """bin_table → 长表 [gene, group, cp10k, pct, n_bins, n_cells_total](spec §3.7 签名)。

    口径(config stats.specificity): scale=mean_cp10k_linear(cp_ 列)、
    aggregation=bin_equal_weight(可行 bin 的算术均值, 与细胞数无关)、
    background=main_groups(只聚合主组行);period 非 None 时限定该期 bin。
    缺失基因(无 cp_ 列)不产行——由 scores 补 status=missing 行。
    """
    if str(cfg.stats_cfg["specificity"].get("scale", "")) != "mean_cp10k_linear":
        raise NotImplementedError("stats.specificity.scale 仅实现 mean_cp10k_linear")
    if str(cfg.stats_cfg["specificity"].get("aggregation", "")) != "bin_equal_weight":
        raise NotImplementedError("stats.specificity.aggregation 仅实现 bin_equal_weight")
    main = set(cfg.main_groups())
    bt = bin_table.loc[bin_table["feasible"].to_numpy() & bin_table["group"].isin(main).to_numpy()]
    if period is not None:
        bt = bt.loc[bt["period"].to_numpy() == period]
    recs = []
    for gene in genes:
        col, pcol = f"cp_{gene}", f"pct_{gene}"
        if col not in bin_table.columns:
            continue
        gb = bt.groupby("group", sort=False)
        means = gb[[col, pcol]].mean()
        nb = gb[col].count()
        nc = gb["n_cells_total"].sum()
        for grp in means.index:
            v = means.at[grp, col]
            if not np.isfinite(v):
                continue
            recs.append({"gene": gene, "group": grp, "cp10k": float(v),
                         "pct": float(means.at[grp, pcol]),
                         "n_bins": int(nb[grp]), "n_cells_total": int(nc[grp])})
    return pd.DataFrame(recs, columns=["gene", "group", "cp10k", "pct", "n_bins", "n_cells_total"])


def scores(agg, genes, cfg) -> pd.DataFrame:
    """aggregate 长表 → gene 级特异性评分(global 或 by-period, 由 agg 的口径决定)。

    列(spec §5 specificity_global + §3.7 差异): gene, status, n_groups_used, groups_used,
    8 组 cp10k_/pct_/fraction_/fraction_8group_, tau, max_group, max_fraction, specific_call。
    NA 政策见模块 docstring ②③⑤。
    """
    sp = cfg.stats_cfg["specificity"]
    tb = cfg.stats_cfg["tau_by_bin"]
    pc = float(sp["pseudocount"])
    min_groups = int(tb["min_groups"])
    tau_thr = float(sp["tau_threshold"])
    frac_thr = float(sp["fraction_threshold"])
    min_expr = float(sp["min_expr"])
    main = list(cfg.main_groups())

    if len(agg):
        wc = agg.pivot(index="gene", columns="group", values="cp10k")
        wp = agg.pivot(index="gene", columns="group", values="pct")
    else:
        wc = wp = pd.DataFrame()
    wc = wc.reindex(index=list(genes), columns=main)
    wp = wp.reindex(index=list(genes), columns=main)

    recs = []
    for gene in genes:
        x = pd.to_numeric(wc.loc[gene]).to_numpy(dtype=float)
        pv = pd.to_numeric(wp.loc[gene]).to_numpy(dtype=float)
        ok = np.isfinite(x)
        k = int(ok.sum())
        main_ok = [g for g, o in zip(main, ok) if o]
        rec = {"gene": gene, "status": "ok" if k else "missing",
               "n_groups_used": k, "groups_used": "|".join(main_ok)}
        for g, xv, pvv in zip(main, x, pv):
            rec[f"cp10k_{g}"] = xv if np.isfinite(xv) else np.nan
            rec[f"pct_{g}"] = pvv if np.isfinite(pvv) else np.nan
        tau = max_frac = np.nan
        call = "NA"
        mg = main_ok[int(np.argmax(x[ok]))] if k else "NA"   # ties 取主组序首个(信息列)
        if k >= min_groups and k >= 2:
            xp = x[ok] + pc
            mx = float(xp.max())
            tau = float((1.0 - xp / mx).sum() / (k - 1))
            frac = xp / xp.sum()
            max_frac = float(frac[int(np.argmax(xp))])
            for gi, g in enumerate(main_ok):
                rec[f"fraction_{g}"] = float(frac[gi])
                if k == len(main):                          # 严格版: 8 组全可行才非 NA
                    rec[f"fraction_8group_{g}"] = float(frac[gi])
            raw_max = float(np.max(x[ok]))
            if raw_max < min_expr:
                call = "not_expressed"
            elif tau >= tau_thr and max_frac >= frac_thr:
                call = "specific"
            else:
                call = "broad"
        rec.update(tau=tau, max_group=mg, max_fraction=max_frac, specific_call=call)
        recs.append(rec)
    df = pd.DataFrame(recs)
    for g in main:
        for pre in ("fraction_", "fraction_8group_"):
            if f"{pre}{g}" not in df.columns:
                df[f"{pre}{g}"] = np.nan
    cols = (["gene", "status", "n_groups_used", "groups_used"]
            + [f"cp10k_{g}" for g in main] + [f"pct_{g}" for g in main]
            + [f"fraction_{g}" for g in main] + [f"fraction_8group_{g}" for g in main]
            + ["tau", "max_group", "max_fraction", "specific_call"])
    return df[cols].reset_index(drop=True)


def tau_by_bin(bin_table, genes, cfg) -> pd.DataFrame:
    """每 bin 的 tau 与各可行组 fraction(spec §3.7 签名;长表 gene×bin)。

    k = 该 bin 可行主组数;k < min_groups → tau/frac 全 NA(左截断 bin 显式存在, 不丢行)。
    frac_<g> 分母 = 该 bin 可行组之和(fraction_in_feasible 口径, 读图指南 Fig1e)。
    """
    sp = cfg.stats_cfg["specificity"]
    tb = cfg.stats_cfg["tau_by_bin"]
    pc = float(sp["pseudocount"])
    min_groups = int(tb["min_groups"])
    main = list(cfg.main_groups())
    bins = (bin_table[["bin", "bin_order"]].drop_duplicates()
            .sort_values("bin_order").reset_index(drop=True))
    bt = bin_table.loc[bin_table["group"].isin(main).to_numpy()]
    recs = []
    for gene in genes:
        col = f"cp_{gene}"
        if col not in bin_table.columns:
            for _, b in bins.iterrows():
                recs.append({"gene": gene, "bin": b["bin"], "bin_order": int(b["bin_order"]),
                             "status": "missing", "n_groups_used": 0, "groups_used": "",
                             "tau": np.nan})
            continue
        blist = bins["bin"].tolist()
        piv = bt.pivot(index="bin", columns="group", values=col).reindex(index=blist, columns=main)
        fpiv = bt.pivot(index="bin", columns="group", values="feasible").reindex(index=blist, columns=main)
        for _, b in bins.iterrows():
            bn = b["bin"]
            xv = pd.to_numeric(piv.loc[bn]).to_numpy(dtype=float)
            feas = fpiv.loc[bn].fillna(False).astype(bool).to_numpy()
            ok = feas & np.isfinite(xv)
            k = int(ok.sum())
            main_ok = [g for g, o in zip(main, ok) if o]
            rec = {"gene": gene, "bin": bn, "bin_order": int(b["bin_order"]),
                   "status": "ok" if k else "no_feasible",
                   "n_groups_used": k, "groups_used": "|".join(main_ok), "tau": np.nan}
            if k >= min_groups and k >= 2:
                xp = xv[ok] + pc
                mx = float(xp.max())
                rec["tau"] = float((1.0 - xp / mx).sum() / (k - 1))
                frac = xp / xp.sum()
                for gi, g in enumerate(main_ok):
                    rec[f"frac_{g}"] = float(frac[gi])
            recs.append(rec)
    df = pd.DataFrame(recs)
    for g in main:
        if f"frac_{g}" not in df.columns:
            df[f"frac_{g}"] = np.nan
    cols = (["gene", "bin", "bin_order", "status", "n_groups_used", "groups_used", "tau"]
            + [f"frac_{g}" for g in main])
    return df[cols].reset_index(drop=True)


def emergence(tau_bin, max_group_map, cfg) -> pd.DataFrame:
    """首个"主组 fraction ≥ tau_by_bin.fraction_threshold"的 bin(spec §3.7 签名/差异③④)。

    max_group_map: dict 或 Series(gene → global max_group, 来自 scores();"NA" → 全 NA 行)。
    仅在 n_groups_used ≥ min_groups 的 bin 上判定;附 first_computable_bin(_order) 与
    left_censored(= emergence 落在首个可算 bin, 左截断语义, Depression_CB §7.4)。
    """
    tb = cfg.stats_cfg["tau_by_bin"]
    thr = float(tb["fraction_threshold"])
    min_groups = int(tb["min_groups"])
    items = (list(max_group_map.items()) if isinstance(max_group_map, dict)
             else list(zip(max_group_map.index, max_group_map)))
    recs = []
    for gene, mg in items:
        mg = str(mg)
        sub = tau_bin.loc[tau_bin["gene"].to_numpy() == gene].sort_values("bin_order")
        comp = sub.loc[sub["n_groups_used"].to_numpy() >= min_groups]
        rec = {"gene": gene, "max_group": mg,
               "first_computable_bin": str(comp.iloc[0]["bin"]) if len(comp) else "NA",
               "first_computable_bin_order": float(comp.iloc[0]["bin_order"]) if len(comp) else np.nan}
        fcol = f"frac_{mg}"
        hit = comp.iloc[0:0]
        if mg != "NA" and fcol in comp.columns and len(comp):
            fv = pd.to_numeric(comp[fcol]).to_numpy(dtype=float)
            hit = comp.loc[np.isfinite(fv) & (fv >= thr)]
        if len(hit):
            r0 = hit.iloc[0]
            rec.update(emergence_bin=str(r0["bin"]), emergence_bin_order=float(r0["bin_order"]),
                       frac_at_emergence=float(r0[fcol]),
                       n_groups_used_at_emergence=int(r0["n_groups_used"]),
                       left_censored=bool(str(r0["bin"]) == rec["first_computable_bin"]))
        else:
            rec.update(emergence_bin="NA", emergence_bin_order=np.nan, frac_at_emergence=np.nan,
                       n_groups_used_at_emergence=np.nan, left_censored=False)
        recs.append(rec)
    cols = ["gene", "max_group", "first_computable_bin", "first_computable_bin_order",
            "emergence_bin", "emergence_bin_order", "frac_at_emergence",
            "n_groups_used_at_emergence", "left_censored"]
    return pd.DataFrame(recs, columns=cols)
