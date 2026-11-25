"""atx.stats_trend — Spearman 趋势(bin 级主、donor 级辅) + BH-FDR(implementation_spec §3.6)。

借鉴(Depression_CB README §7.2 已验证口径):可行 bin ≤ exact_perm_max_bins 用精确置换 p
(scipy 在 rho=±1 时退化出 p=0, 假显著);FDR 族 = species×group(本函数恰好按一物种×一组
调用 → q_bin 即族内 BH);时间轴统一 bin_order;门槛 = 可行 bin ≥ min_bins_for_trend
(组级 override)且可行 bin 总细胞 ≥ min_cells_total_for_trend。

本项目差异(§3.6):
① min_bins_for_trend 按组可覆盖(config bins.per_group_override: UBC/Immune=4);
② low_power 标记 = 可行 bin 数低于**标准**阈值、仅靠组级 override 才被分析
   (承 Depression_CB 鼠 MICROGLIA 4-bin 全链路标注先例);
③ 置换排列表按 n 缓存一次(模块级 _PERM_CACHE, n! 穷举上限 _EXACT_HARD_MAX 防组合爆炸)。

status 受控词表: ok / missing / insufficient_bins / insufficient_cells。
确定性: 全程无 RNG(精确置换=穷举, scipy 路=解析 p)→ 两次重跑逐值一致(R8)。
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

_PERM_CACHE: dict[int, np.ndarray] = {}
_EXACT_HARD_MAX = 8      # n! 穷举硬上限(8!=40320);结构保护, 非统计参数
_MIN_DONOR_OBS = 3       # 数学退化保护: n=2 时 Spearman 恒 ±1 无意义;非可调统计阈值


def bh_qvalues(p) -> np.ndarray:
    """Benjamini-Hochberg step-up q 值。NaN 保留;输入可乱序, 输出按原位还原。"""
    p = np.asarray(p, dtype=float).ravel()
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    po = p[ok]
    n = po.size
    if n == 0:
        return q
    order = np.argsort(po, kind="stable")
    qv = po[order] * n / np.arange(1, n + 1)
    qv = np.minimum.accumulate(qv[::-1])[::-1]      # step-up: 从大 p 向小单调化
    out = np.empty(n)
    out[order] = np.clip(qv, 0.0, 1.0)
    q[ok] = out
    return q


def _perm_matrix(n: int) -> np.ndarray:
    """arange(n) 全排列 (n!, n), 按 n 缓存一次(spec §3.6 差异③)。"""
    if n not in _PERM_CACHE:
        if n > _EXACT_HARD_MAX:
            raise ValueError(f"精确置换穷举仅限 n≤{_EXACT_HARD_MAX}(n! 组合爆炸), 得 n={n}")
        _PERM_CACHE[n] = np.array(list(itertools.permutations(np.arange(n))))
    return _PERM_CACHE[n]


def _rho_ranks(rt: np.ndarray, rx: np.ndarray) -> float:
    """秩向量间 Pearson(= Spearman, 平均秩处理 ties);常数列 → nan。"""
    rtc = rt - rt.mean()
    rxc = rx - rx.mean()
    den = np.sqrt((rtc ** 2).sum() * (rxc ** 2).sum())
    return float((rtc * rxc).sum() / den) if den > 0 else np.nan


def spearman_exact(t, x, exact_max_n: int = 6):
    """Spearman rho + 双侧 p → (rho, p)(spec §3.6 签名)。

    NaN 对先剔除;n<2 或任一列恒定 → (nan, nan)。
    n ≤ exact_max_n: 精确置换 p = #{|rho_perm| ≥ |rho_obs|−1e-12} / n!(穷举;观测排列
    自身计入 → p ≥ 1/n!, 不会像 scipy 那样在 rho=±1 时给 p=0)。观测 rho 与置换分布用
    同一秩-Pearson 公式, 保证计数一致。
    n > exact_max_n: scipy 解析 p。
    """
    t = np.asarray(t, dtype=float).ravel()
    x = np.asarray(x, dtype=float).ravel()
    if t.shape != x.shape:
        raise ValueError(f"t/x 长度不一致: {t.shape} vs {x.shape}")
    m = np.isfinite(t) & np.isfinite(x)
    t, x = t[m], x[m]
    n = int(t.size)
    if n < 2:
        return np.nan, np.nan
    rt = rankdata(t)
    rx = rankdata(x)
    if np.ptp(rt) == 0 or np.ptp(rx) == 0:
        return np.nan, np.nan
    rho = _rho_ranks(rt, rx)
    if n <= exact_max_n:
        P = _perm_matrix(n)
        RX = rankdata(x[P], axis=1)                  # (n!, n) 每行平均秩
        rtc = rt - rt.mean()
        RXc = RX - RX.mean(axis=1, keepdims=True)
        num = RXc @ rtc
        den = np.sqrt((RXc ** 2).sum(axis=1) * (rtc ** 2).sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            rp = np.where(den > 0, num / den, 0.0)
        p = float(np.mean(np.abs(rp) >= abs(rho) - 1e-12))
        return rho, p
    return rho, float(spearmanr(t, x).pvalue)


def trend_table(bin_table, cube, genes, group, spcfg, cfg) -> pd.DataFrame:
    """一物种×一组的 gene 级趋势表(spec §3.6 签名;§5 trend_by_group.tsv 列)。

    列: gene, group, n_bins, n_cells_total, rho_bin, p_bin, q_bin, rho_donor, p_donor,
        n_obs_donor, low_power, status。
    门(config bins, 全读不硬编码): n_bins ≥ min_bins_for_trend(组级 override)且可行 bin
    总细胞 ≥ min_cells_total_for_trend, 否则 NA 行 + status。q_bin = 本表(=species×group 族)
    内 BH。donor 级为辅证: 观测数 ≥ _MIN_DONOR_OBS 才算(退化保护), 无独立细胞数门。
    bin_order 与 spcfg 一致性硬校验(collapse_adult_bins 敏感性必须传 collapsed spcfg)。
    """
    bins_cfg = cfg.bins_cfg
    std_min = int(bins_cfg["min_bins_for_trend"])
    ovr = (bins_cfg.get("per_group_override") or {}).get(group) or {}
    min_bins = int(ovr.get("min_bins_for_trend", std_min))
    min_total = int(bins_cfg["min_cells_total_for_trend"])
    tcfg = cfg.stats_cfg["trend"]
    if str(tcfg.get("method", "spearman")) != "spearman":
        raise NotImplementedError(f"stats.trend.method={tcfg.get('method')!r}: 仅实现 spearman")
    exact_max = int(tcfg["exact_perm_max_bins"])

    bt = bin_table[bin_table["group"].to_numpy() == group]
    if len(bt):
        exp = bt["bin"].map(spcfg.bin_order).to_numpy()
        if not np.array_equal(exp, bt["bin_order"].to_numpy()):
            raise RuntimeError("bin_table.bin_order 与 spcfg.bin_order 不一致——spcfg 传错"
                               "(collapse_adult_bins 敏感性须用 effective_species 的产物)")
    bt = bt.sort_values("bin_order")

    gmask = None
    gene_idx: dict[str, int] = {}
    if cube is not None:
        gmask = cube.obs["group"].to_numpy() == group
        gene_idx = {g: i for i, g in enumerate(cube.var["gene_key"])}

    rows = []
    for gene in genes:
        col = f"log_{gene}"
        rec: dict = {"gene": gene, "group": group}
        if col not in bin_table.columns:
            rec.update(status="missing", n_bins=0, n_cells_total=0, low_power=True)
            rows.append(rec)
            continue
        sub = bt.loc[bt["feasible"].to_numpy() & bt[col].notna().to_numpy()]
        n_bins = len(sub)
        n_tot = int(sub["n_cells_total"].sum())
        rec.update(n_bins=n_bins, n_cells_total=n_tot, low_power=bool(n_bins < std_min))
        if n_bins < min_bins:
            rec["status"] = "insufficient_bins"
        elif n_tot < min_total:
            rec["status"] = "insufficient_cells"
        else:
            rho, p = spearman_exact(sub["bin_order"].to_numpy(float),
                                    sub[col].to_numpy(float), exact_max)
            rec.update(rho_bin=rho, p_bin=p, status="ok")
            if gmask is not None and gene in gene_idx and gmask.any():
                cx = cube.X[:, gene_idx[gene]]
                xd = np.asarray(cx.todense() if hasattr(cx, "todense") else cx).ravel()
                xd = xd.astype(float)[gmask]
                td = cube.obs.loc[gmask, "bin_order"].to_numpy(float)
                mm = np.isfinite(xd) & np.isfinite(td)
                rec["n_obs_donor"] = int(mm.sum())
                if mm.sum() >= _MIN_DONOR_OBS:
                    rd, pdon = spearman_exact(td[mm], xd[mm], exact_max)
                    rec.update(rho_donor=rd, p_donor=pdon)
        rows.append(rec)

    df = pd.DataFrame(rows)
    df["q_bin"] = (bh_qvalues(df["p_bin"].to_numpy(float))
                   if "p_bin" in df.columns else np.nan)
    cols = ["gene", "group", "n_bins", "n_cells_total", "rho_bin", "p_bin", "q_bin",
            "rho_donor", "p_donor", "n_obs_donor", "low_power", "status"]
    for c in cols:
        if c not in df.columns:
            df[c] = np.nan
    return df[cols].reset_index(drop=True)
