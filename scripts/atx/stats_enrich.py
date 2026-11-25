"""atx.stats_enrich — 第 3/4 层统计(implementation_spec §3.9, 全新模块)。

职责: 基因→模块指派(s07 §6)、机制类×模块置换富集(Fig3a/3b)、分层置换 Spearman
(s08 §7 主检验)、Kruskal–Wallis + 两两 Mann–Whitney + Cliff's delta(Fig4b)、
Fisher 列联退化为逐格 2×2(Fig3b/4c)。

口径(spec §6/§7 + Depression_CB README 借鉴):
- 模块 = G* × P*: G* = max_group(global tau argmax);max_fraction <
  broad_fraction_threshold(config enrichment, spec §6 L484)→ G* = broad, P* 改取
  **ALL 组**峰 period(口径决策: spec 未指明 broad 的 P* 来源;ALL 为参考组、全骨架
  恒有行, 以全小脑轨迹定 broad 基因的期归属;max_fraction NA[k<min_groups, NA 政策]
  → 模块 NA 而非 broad——"NA 不当 False"语义, Depression_CB §3-14);
  G*(broad 用 ALL)可行 bin < min_bins_for_trend(组级 override)→ 模块 NA 不进检验;
- 富集(spec §6): O = 类×模块计数, E = n_c × n_m / N, log2FE = log2((O+.5)/(E+.5));
  置换 = 固定模块、打乱类标签 n_perm 次;双侧 p = 比例(|log2FE_perm| ≥ |log2FE_obs|
  − 1e-12);**观测排列(恒等)计入 n_perm → p ≥ 1/n_perm**(与 stats_trend 精确置换
  约定一致, 防 p=0 假显著);BH 族 = 单次调用的全部 (c,m) 格(= 物种内全格);
- 分层置换(spec §7): y 在 strata 内打乱(s08 传 mech_class 列, config
  onset.stratify_by), 消除机制类×发病年龄共线混杂;同时报未分层 p_naive 对照;
  rho = 秩化 Pearson(与置换分布同公式, 计数一致);y 的置换 ≡ 其秩向量置换
  (多重集不变, 平均秩随位置平移)→ 秩只算一次;
- Cliff's delta = (#{a>b} − #{a<b}) / (n_a·n_b), 正 = a 倾向更大, ties 贡献 0;
- Fisher(spec §6 L490 "2×k 退化为逐列 2×2 + BH"): 逐 (类, 列) 格 2×2 =
  [[n_ab, n_a−n_ab], [n_b−n_ab, N−n_a−n_b+n_ab]], fisher_exact 双侧, BH 族=全格;
- 全部置换 np.random.default_rng(seed);**正式调用一律传 cfg.seed**(spec §3.9
  "所有置换用 default_rng(cfg.seed)";签名中 seed 为显式参数, 默认 0 仅兜底);
  同 seed 两次重跑逐值一致(R8)。
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, kruskal, mannwhitneyu, rankdata

from .stats_trend import bh_qvalues

_MIN_PERM_OBS = 3      # 结构保护: n=2 时 Spearman 恒 ±1, 置换检验无意义;非可调统计阈值

# 残差兜底类(CL-027, 2026-11-24 用户裁决; plan §0-3/§4.1; spec §6/§7)。
# mech_class 的 other/unknown 是"小类并入"与"无法归入现有类"的落脚点, **不是一种机制**:
#   · s07 机制类×模块 / 机制类×shape 富集: 必须**排除**(perm_enrichment 内有硬保护, 见下);
#     其成员改做单基因描述。config: stats.enrichment.exclude_mech_classes。
#   · s08 分层置换: **不排除**(排除会改样本集、损失功效), 改为自成 stratum 并标 low n;
#     n=3 时层内仅 6 种排列 → 该层分层 p 无解释力, 结论以 p_naive 与"剔除残差层"敏感性为准。
# v1.0 冻结集下残差类成员(core): WWOX(SCAR12)、CCDC88C(SCA40)、MME(SCA43)。
RESIDUAL_MECH_CLASSES = ("other", "unknown")

_MOD_COLS = ["gene", "max_group", "max_fraction", "g_star", "peak_group_used", "peak_bin",
             "peak_period", "n_feasible_bins", "module", "status"]
_ENR_COLS = ["a", "b", "obs", "exp", "log2fe", "p", "q"]
_PAIR_COLS = ["group1", "group2", "n1", "n2", "p_mw", "q_mw", "cliffs_delta"]
_FISH_COLS = ["a", "b", "n_ab", "n_a", "n_b", "n_total", "odds_ratio", "p_fisher", "q_fisher"]


def assign_modules(spec_global, peak_by_group, cfg) -> pd.DataFrame:
    """基因 → G*×P* 模块指派(spec §3.9 签名;§6 规则;s07 每物种调用一次)。

    输入: spec_global = stats_specificity.scores() 产物(读 gene/status/max_group/
    max_fraction);peak_by_group = stats_pattern.peak_table() 产物。
    输出列: gene, max_group, max_fraction, g_star, peak_group_used, peak_bin,
    peak_period, n_feasible_bins, module("G*|P*"), status。
    status 受控词表: ok / missing_gene / fraction_na / insufficient_bins_gstar。
    shape(=temporal_pattern(G*)) 与 species 列由 s07 merge temporal_patterns 追加
    (§3.9 签名不含该输入, 保持与 spec 一致)。
    """
    ecfg = cfg.stats_cfg["enrichment"]
    broad_thr = float(ecfg["broad_fraction_threshold"])
    bins_cfg = cfg.bins_cfg
    std_min = int(bins_cfg["min_bins_for_trend"])
    ovr = bins_cfg.get("per_group_override") or {}
    pk: dict[tuple[str, str], dict] = {}
    for row in peak_by_group.to_dict("records"):
        pk[(str(row["gene"]), str(row["group"]))] = row
    recs = []
    for row in spec_global.to_dict("records"):
        gene = str(row["gene"])
        mg = str(row.get("max_group", "NA"))
        mf = row.get("max_fraction", np.nan)
        mf = float(mf) if mf is not None and np.isfinite(mf) else np.nan
        rec = {"gene": gene, "max_group": mg, "max_fraction": mf}

        def _na(g_star: str, src: str, status: str, nfeas: int = 0) -> dict:
            rec.update(g_star=g_star, peak_group_used=src, peak_bin="NA", peak_period="NA",
                       n_feasible_bins=nfeas, module="NA", status=status)
            return rec

        if str(row.get("status", "ok")) != "ok" or mg == "NA":
            recs.append(_na("NA", "NA", "missing_gene"))
            continue
        if not np.isfinite(mf):
            recs.append(_na(mg, "NA", "fraction_na"))
            continue
        if mf < broad_thr:
            g_star, src = "broad", "ALL"
            min_bins = std_min                      # ALL 无组级 override
        else:
            g_star, src = mg, mg
            min_bins = int((ovr.get(g_star) or {}).get("min_bins_for_trend", std_min))
        prow = pk.get((gene, src))
        nfeas = int(prow["n_feasible_bins"]) if prow is not None else 0
        if prow is None or nfeas < min_bins or not isinstance(prow.get("peak_period"), str):
            recs.append(_na(g_star, src, "insufficient_bins_gstar", nfeas))
            continue
        rec.update(g_star=g_star, peak_group_used=src, peak_bin=str(prow["peak_bin"]),
                   peak_period=str(prow["peak_period"]), n_feasible_bins=nfeas,
                   module=f"{g_star}|{prow['peak_period']}", status="ok")
        recs.append(rec)
    return pd.DataFrame(recs, columns=_MOD_COLS)


def perm_enrichment(labels_a, labels_b, n_perm, seed=0, allow_residual: bool = False) -> pd.DataFrame:
    """类标签 × 模块标签置换富集(spec §3.9 签名;§6 检验)。

    labels_a/labels_b: 等长类别标签(每基因一个;NA/NaN 须由调用方先行过滤——
    §6 "模块=NA 不进检验")。输出 = a_levels × b_levels 全格(升序),
    列 [a, b, obs, exp, log2fe, p, q];q = 全格 BH(族=单次调用)。
    置换: 固定 b、打乱 a(行置换 one-hot 矩阵);恒等排列计入 → p ≥ 1/n_perm。

    **残差类硬保护(CL-027; spec §6)**: labels_a 为 mech_class 时, 其中**不得**含
    RESIDUAL_MECH_CLASSES(other/unknown)——残差类不是一种机制, "残差类富集于某模块"
    无生物学解释。调用方(s07)必须先过滤这些基因并打印清单; 若确需纳入(仅方法学对照),
    显式传 allow_residual=True 且在报告注明。默认含残差类 → ValueError(防静默算错)。
    """
    n_perm = int(n_perm)
    if n_perm < 1:
        raise ValueError(f"n_perm 须 ≥1, 得 {n_perm}")
    # dtype=object 必写: np.asarray(["A", np.nan]) 会提升为 <U32 字符串数组,
    # nan 被转成 "nan" 字符串 → pd.isna 漏检(实测陷阱)
    a = np.asarray(labels_a, dtype=object).ravel()
    b = np.asarray(labels_b, dtype=object).ravel()
    if a.shape != b.shape:
        raise ValueError(f"labels 长度不一致: {a.shape} vs {b.shape}")
    if a.size == 0:
        return pd.DataFrame(columns=_ENR_COLS)
    if pd.isna(a).any() or pd.isna(b).any():
        raise ValueError("labels 含 NA/NaN——须先过滤模块=NA 的基因(spec §6 不进检验)")
    a = a.astype(str)
    b = b.astype(str)
    if not allow_residual:
        hit = sorted({v for v in a.tolist() if v in RESIDUAL_MECH_CLASSES})
        if hit:
            raise ValueError(
                f"labels_a 含残差兜底类 {hit} —— 按 CL-027 / spec §6, mech_class 的 other/unknown "
                "不是一种机制, 不得进入富集网格。调用方(s07)须先过滤这些基因并打印清单, 改做单基因描述; "
                "确需方法学对照才传 allow_residual=True 且在报告注明。")
    a_levels = sorted(set(a.tolist()))
    b_levels = sorted(set(b.tolist()))
    A = (a[:, None] == np.asarray(a_levels)[None, :]).astype(float)      # (n, na)
    B = (b[:, None] == np.asarray(b_levels)[None, :]).astype(float)      # (n, nb)
    N = int(a.size)
    E = np.outer(A.sum(0), B.sum(0)) / N
    O = A.T @ B
    lfe0 = np.log2((O + 0.5) / (E + 0.5))
    abs0 = np.abs(lfe0)
    rng = np.random.default_rng(int(seed))
    cnt = np.ones_like(lfe0)                        # 恒等排列(观测自身)计入
    for _ in range(n_perm - 1):
        idx = rng.permutation(N)
        lp = np.log2((A[idx].T @ B + 0.5) / (E + 0.5))
        cnt += np.abs(lp) >= abs0 - 1e-12
    p = cnt / float(n_perm)
    q = bh_qvalues(p.ravel()).reshape(p.shape)
    recs = [{"a": av, "b": bv, "obs": int(O[i, j]), "exp": float(E[i, j]),
             "log2fe": float(lfe0[i, j]), "p": float(p[i, j]), "q": float(q[i, j])}
            for i, av in enumerate(a_levels) for j, bv in enumerate(b_levels)]
    return pd.DataFrame(recs, columns=_ENR_COLS)


def stratified_perm_spearman(x, y, strata, n_perm, seed=0):
    """分层置换 Spearman → (rho, p_strat, p_naive)(spec §3.9 签名;§7 主检验)。

    NaN 的 x/y 对先剔除(strata 同步);n < _MIN_PERM_OBS 或 x/y 恒定 → (nan,)*3。
    分层置换: 每 stratum 内部独立重排 y 的位置;p_naive 用全排列对照。
    两路置换均把恒等排列计入 → p ≥ 1/n_perm;双侧 |rho| 比较(容差 1e-12)。
    秩只算一次: rankdata(y[idx]) == rankdata(y)[idx](置换不改变多重集与平均秩)。
    """
    n_perm = int(n_perm)
    if n_perm < 1:
        raise ValueError(f"n_perm 须 ≥1, 得 {n_perm}")
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    s = np.asarray(strata).ravel()
    if not (x.shape == y.shape == s.shape):
        raise ValueError(f"x/y/strata 长度不一致: {x.shape}/{y.shape}/{s.shape}")
    m = np.isfinite(x) & np.isfinite(y)
    x, y, s = x[m], y[m], s[m]
    n = int(x.size)
    if n < _MIN_PERM_OBS:
        return np.nan, np.nan, np.nan
    rx, ry = rankdata(x), rankdata(y)
    rxc, ryc = rx - rx.mean(), ry - ry.mean()
    den = float(np.sqrt((rxc ** 2).sum() * (ryc ** 2).sum()))
    if den == 0:
        return np.nan, np.nan, np.nan
    rho = float((rxc * ryc).sum() / den)
    abs0 = abs(rho)
    rng = np.random.default_rng(int(seed))
    strata_pos = [np.flatnonzero(s == u) for u in np.unique(s)]
    base = np.arange(n)
    c_strat = c_naive = 1.0                         # 恒等排列计入两路
    for _ in range(n_perm - 1):
        idx = base.copy()
        for pos in strata_pos:
            if pos.size > 1:
                idx[pos] = rng.permutation(pos)
        c_strat += abs(float((rxc * ryc[idx]).sum() / den)) >= abs0 - 1e-12
        idx2 = rng.permutation(n)
        c_naive += abs(float((rxc * ryc[idx2]).sum() / den)) >= abs0 - 1e-12
    return rho, c_strat / float(n_perm), c_naive / float(n_perm)


def cliffs_delta(a, b) -> float:
    """Cliff's delta = (#{a>b} − #{a<b})/(n_a·n_b);正 = a 倾向更大;ties 计 0。

    NaN 各自剔除;任一为空 → nan。规模为基因级(≤数百), 外积比较足够。
    """
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return np.nan
    gt = int((a[:, None] > b[None, :]).sum())
    lt = int((a[:, None] < b[None, :]).sum())
    return (gt - lt) / float(a.size * b.size)


def group_compare(values, groups):
    """KW + 两两 MW(BH)+ Cliff's delta → (kw_p, pairwise_df)(spec §3.9 签名;§7 类别检验)。

    NaN 值先剔除(组标签同步);非空组 < 2 → (nan, 空表)。pairwise 列 =
    [group1, group2, n1, n2, p_mw, q_mw, cliffs_delta];q = 单次调用内全对 BH;
    组水平升序 → 对 (g1<g2) 确定性枚举。全同值输入 scipy 会给 NaN/警告——透传不掩盖。
    """
    v = np.asarray(values, dtype=float).ravel()
    g = np.asarray(groups).ravel()
    if v.shape != g.shape:
        raise ValueError(f"values/groups 长度不一致: {v.shape} vs {g.shape}")
    m = np.isfinite(v)
    v, g = v[m], g[m]
    levels = sorted({str(x) for x in g.tolist()})
    samples = {lv: v[g.astype(str) == lv] for lv in levels}
    samples = {lv: s for lv, s in samples.items() if s.size > 0}
    if len(samples) < 2:
        return np.nan, pd.DataFrame(columns=_PAIR_COLS)
    lvls = sorted(samples)
    kw_p = float(kruskal(*[samples[lv] for lv in lvls]).pvalue)
    recs = []
    for i in range(len(lvls)):
        for j in range(i + 1, len(lvls)):
            a_, b_ = samples[lvls[i]], samples[lvls[j]]
            recs.append({"group1": lvls[i], "group2": lvls[j],
                         "n1": int(a_.size), "n2": int(b_.size),
                         "p_mw": float(mannwhitneyu(a_, b_, alternative="two-sided").pvalue),
                         "cliffs_delta": cliffs_delta(a_, b_)})
    df = pd.DataFrame(recs, columns=_PAIR_COLS)
    df["q_mw"] = bh_qvalues(df["p_mw"].to_numpy(dtype=float))
    return kw_p, df[_PAIR_COLS].reset_index(drop=True)


def fisher_by_column(labels_a, labels_b) -> pd.DataFrame:
    """2×k 列联退化为逐格 2×2 Fisher + BH(spec §3.9 职责" Fisher 列联";§6 L490)。

    每 (a 水平, b 水平) 格: 2×2 = [[n_ab, n_a−n_ab], [n_b−n_ab, N−n_a−n_b+n_ab]],
    fisher_exact 双侧;q = 单次调用全格 BH。NA/NaN 标签 → ValueError(调用方先过滤)。
    确定性(无 RNG)。列 = [a, b, n_ab, n_a, n_b, n_total, odds_ratio, p_fisher, q_fisher]。
    """
    a = np.asarray(labels_a, dtype=object).ravel()      # dtype=object: 见 perm_enrichment 注
    b = np.asarray(labels_b, dtype=object).ravel()
    if a.shape != b.shape:
        raise ValueError(f"labels 长度不一致: {a.shape} vs {b.shape}")
    if a.size == 0:
        return pd.DataFrame(columns=_FISH_COLS)
    if pd.isna(a).any() or pd.isna(b).any():
        raise ValueError("labels 含 NA/NaN——须先过滤(spec §6)")
    a = a.astype(str)
    b = b.astype(str)
    N = int(a.size)
    a_levels = sorted(set(a.tolist()))
    b_levels = sorted(set(b.tolist()))
    recs = []
    ps = []
    for av in a_levels:
        am = a == av
        n_a = int(am.sum())
        for bv in b_levels:
            bm = b == bv
            n_b = int(bm.sum())
            n_ab = int((am & bm).sum())
            table = [[n_ab, n_a - n_ab], [n_b - n_ab, N - n_a - n_b + n_ab]]
            orr, pf = fisher_exact(table, alternative="two-sided")
            recs.append({"a": av, "b": bv, "n_ab": n_ab, "n_a": n_a, "n_b": n_b,
                         "n_total": N, "odds_ratio": float(orr), "p_fisher": float(pf)})
            ps.append(float(pf))
    df = pd.DataFrame(recs, columns=_FISH_COLS)
    df["q_fisher"] = bh_qvalues(np.asarray(ps))
    return df[_FISH_COLS].reset_index(drop=True)
