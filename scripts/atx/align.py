"""atx.align — 第 5 层跨物种对齐(implementation_spec §3.10; 口径 §4.4/§4.5A/§8.1)。

职责(spec §3.10): 跨物种 dev_state 剖面、共识状态排序、组内状态子集、保真度评分。

口径(§4.5 A, 硬约束):
  · 两物种**绝不合并表达矩阵**, 只在 pseudobulk 统计量层面比较;
  · 跨物种时间轴 = 共享 `dev_state`(§4.4; config `state_alignment.group_states` 列出各组的状态子集),
    **不做任何人鼠年龄/日历换算**;
  · 每状态取该物种该状态全部核的均值(scope=ALL)或组内核的均值(scope=组);
  · 状态剖面也按 `state_alignment.min_cells` 判可行; 不可行 → **NA**(硬约束 5: 不插值、不填 0);
  · 状态顺序 = 共识序(`order_rule: consensus_mean_time_rank`): 每物种内按状态的**平均 bin_order**
    排秩, 两物种秩取平均后升序;
  · 热图/轨迹作图时物种内逐基因 z-score、两物种独立 colorbar、**只比形状不比幅值**(§3.10 借鉴条)。

签名与 spec 的一致性说明(硬约束 10 留痕):
  `state_profiles(X_log, obs, genes, states, min_cells, group_mask=None)`、
  `consensus_order(profiles_h, profiles_m)`、
  `fidelity(trend_h, trend_m, spec_h, spec_m, prof_h, prof_m, cfg)` 三个签名按 spec 原样保留;
  其中 (1) state_profiles 返回的**列序即该物种的平均 bin_order 秩**, 这样 consensus_order 只需两个
  profile 就能算共识序(与 spec 签名一致); (2) fidelity 的 prof_h/prof_m 为 {scope: 基因×状态宽表} 字典
  (spec 未规定类型; 需要 ALL 与各组两套剖面才能同时算 traj_r 与 peak_state_match)。
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

NA = "__NA__"           # atx.groups 的 NA 哨兵(obs 字符串化后缺失值)
MIN_ABS_RHO = 0.2       # spec §8.1: rho 符号一致性只在两侧 |rho| ≥ 0.2 时判定, 否则 NA
MIN_SHARED_STATES = 3   # spec §8.1: traj_r 需 ≥3 个共享可行状态


# --------------------------------------------------------------------------- 剖面
def state_meta(obs: pd.DataFrame, states: list[str], group_mask: np.ndarray | None = None) -> pd.DataFrame:
    """每个 dev_state 的核数 / 可行性 / 平均 bin_order(共识序的输入)。"""
    ds = obs["dev_state"].astype(str).to_numpy()
    bo = pd.to_numeric(obs["bin_order"], errors="coerce").to_numpy(dtype=float)
    gm = np.ones(len(ds), dtype=bool) if group_mask is None else np.asarray(group_mask, dtype=bool)
    recs = []
    for s in states:
        m = gm & (ds == s)
        n = int(m.sum())
        recs.append({"dev_state": s, "n_cells": n,
                     "mean_bin_order": float(np.nanmean(bo[m])) if n and np.isfinite(bo[m]).any() else np.nan})
    return pd.DataFrame(recs).set_index("dev_state")


def state_profiles(X_log, obs: pd.DataFrame, genes: list[str], states: list[str],
                   min_cells: int, group_mask: np.ndarray | None = None) -> pd.DataFrame:
    """基因 × dev_state 的 mean log1p CP10K 剖面; n_cells < min_cells → 该状态整列 NA。

    X_log: 细胞 × 基因 稀疏矩阵(列序 == genes); obs: 含 dev_state / bin_order;
    **返回列序 = 该物种内平均 bin_order 升序**(即该物种的秩序, 供 consensus_order 直接用)。
    """
    meta = state_meta(obs, states, group_mask)
    order = [s for s in meta.sort_values(["mean_bin_order", "dev_state"], na_position="last").index]
    ds = obs["dev_state"].astype(str).to_numpy()
    gm = np.ones(len(ds), dtype=bool) if group_mask is None else np.asarray(group_mask, dtype=bool)
    out = pd.DataFrame(np.nan, index=list(genes), columns=order, dtype=float)
    for s in order:
        if int(meta.loc[s, "n_cells"]) < int(min_cells):
            continue                                   # 不可行 → 保持 NA(不插值不填 0)
        m = gm & (ds == s)
        out[s] = np.asarray(X_log[m].mean(axis=0)).ravel()
    return out


def consensus_order(profiles_h: pd.DataFrame, profiles_m: pd.DataFrame) -> list[str]:
    """共识状态序: 两物种各自列序(=平均 bin_order 秩)取平均后升序; 单物种独有状态用其单侧秩。"""
    rh = {s: i for i, s in enumerate(profiles_h.columns)}
    rm = {s: i for i, s in enumerate(profiles_m.columns)}
    states = sorted(set(rh) | set(rm))

    def _avg(s: str) -> float:
        vals = [r for r in (rh.get(s), rm.get(s)) if r is not None]
        return sum(vals) / len(vals)

    return sorted(states, key=lambda s: (_avg(s), s))


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    """逐基因 z-score(跨列); 全 NA 或零方差 → 全 NA(不填 0)。作图/形状比较用。"""
    mu = df.mean(axis=1)
    sd = df.std(axis=1)
    return df.sub(mu, axis=0).div(sd.where(sd > 0), axis=0)


# --------------------------------------------------------------------------- 保真度
def _traj_r(a: pd.Series, b: pd.Series) -> tuple[float, int]:
    """共享可行状态上的 Pearson r(spec §8.1 traj_r)。

    逐物种 z-score 后求 Pearson 与直接对原始值求 Pearson **等价**(z 是逐向量仿射变换),
    故这里直接对共享可行状态的原值计算; 共享可行状态 < MIN_SHARED_STATES 或任一侧零方差 → (nan, n)。
    """
    ok = a.notna() & b.notna()
    n = int(ok.sum())
    if n < MIN_SHARED_STATES:
        return float("nan"), n
    x, y = a[ok].to_numpy(float), b[ok].to_numpy(float)
    if np.std(x) == 0 or np.std(y) == 0:
        return float("nan"), n
    r = pearsonr(x, y).statistic
    return (float(r) if np.isfinite(r) else float("nan")), n


def fidelity(trend_h: pd.DataFrame, trend_m: pd.DataFrame, spec_h: pd.DataFrame,
             spec_m: pd.DataFrame, prof_h: dict[str, pd.DataFrame],
             prof_m: dict[str, pd.DataFrame], cfg) -> pd.DataFrame:
    """4 指标 + 加权保真分 + n_metrics + 证据等级(spec §8.1)。行粒度 = gene × group。

    指标(缺证据 → NA, 不当 0):
      rho_sign_concordant  两侧 |rho_bin| ≥ 0.2 时判 sign 是否一致(1/0), 否则 NA
      traj_r               组内共享 dev_state 子集上的 Pearson r(≥3 个共享可行状态)
      max_group_match      两物种 global max_group 是否一致(gene 级, 广播到该基因各组行)
      peak_state_match     各自峰 dev_state 在**共识序**中相差 ≤ peak_state_tolerance(1/0)
    保真分 = 非 NA 指标的加权均值(权重 config `stats.conservation.*_weight`, NA 项按剩余权重重归一);
    n_metrics = 非 NA 指标数; 证据等级 = n_metrics ≥3 high / 2 medium / ≤1 low。
    """
    ccfg = cfg.stats_cfg["conservation"]
    W = {"rho_sign_concordant": float(ccfg["rho_sign_weight"]),
         "traj_r": float(ccfg["traj_r_weight"]),
         "max_group_match": float(ccfg["max_group_weight"]),
         "peak_state_match": float(ccfg["peak_state_weight"])}
    tol = int(ccfg["peak_state_tolerance"])
    sa = cfg.raw["state_alignment"]
    group_states = {k: list(v) for k, v in (sa.get("group_states") or {}).items()}
    main_groups = [g for g in cfg.main_groups() if g in group_states]

    order = consensus_order(prof_h["ALL"], prof_m["ALL"])
    pos = {s: i for i, s in enumerate(order)}

    def _idx(df: pd.DataFrame, keys=("gene", "group")) -> pd.DataFrame:
        return df.set_index(list(keys))

    th, tm = _idx(trend_h), _idx(trend_m)
    mg_h = spec_h.set_index("gene")["max_group"].astype(str)
    mg_m = spec_m.set_index("gene")["max_group"].astype(str)
    genes = sorted(set(trend_h["gene"].astype(str)) & set(trend_m["gene"].astype(str)))

    def _rho(t, gene, group):
        try:
            r = t.loc[(gene, group)]
        except KeyError:
            return float("nan"), "missing"
        if isinstance(r, pd.DataFrame):
            r = r.iloc[0]
        if str(r.get("status", "ok")) != "ok":
            return float("nan"), str(r.get("status"))
        v = r.get("rho_bin")
        return (float(v) if pd.notna(v) else float("nan")), "ok"

    def _peak(prof, gene, states):
        """组内(限定该组的状态子集)峰 dev_state; 无可行状态 → None。"""
        if gene not in prof.index:
            return None
        row = prof.loc[gene]
        cand = [s for s in states if s in row.index and pd.notna(row[s])]
        if not cand:
            return None
        return max(cand, key=lambda s: (float(row[s]), -pos.get(s, 0)))

    recs = []
    for gene in genes:
        a, b = mg_h.get(gene, "NA"), mg_m.get(gene, "NA")
        # NA 政策(硬约束 4): max_group 缺失(NA/nan/非字符串)→ 不可判, 记 NA 而非 0/1。
        # 防御调用方读表把 "NA" 解析成 NaN(.astype(str)→"nan")而静默伪造 match=0 或 1。
        a = a if isinstance(a, str) and a not in ("NA", "nan", "") else "NA"
        b = b if isinstance(b, str) and b not in ("NA", "nan", "") else "NA"
        mgm = float("nan") if (a == "NA" or b == "NA") else float(a == b)
        for g in main_groups:
            states = group_states[g]
            rh, sh = _rho(th, gene, g)
            rm, sm = _rho(tm, gene, g)
            if np.isfinite(rh) and np.isfinite(rm) and abs(rh) >= MIN_ABS_RHO and abs(rm) >= MIN_ABS_RHO:
                rs = float(np.sign(rh) == np.sign(rm))
            else:
                rs = float("nan")
            ph = prof_h.get(g)
            pm = prof_m.get(g)
            if ph is not None and pm is not None and gene in ph.index and gene in pm.index:
                cols = [s for s in states if s in ph.columns and s in pm.columns]
                tr, n_shared = _traj_r(ph.loc[gene, cols], pm.loc[gene, cols])
                pk_h = _peak(ph, gene, states)
                pk_m = _peak(pm, gene, states)
            else:
                tr, n_shared, pk_h, pk_m = float("nan"), 0, None, None
            if pk_h is None or pk_m is None:
                psm = float("nan")
            else:
                psm = float(abs(pos.get(pk_h, -999) - pos.get(pk_m, -999)) <= tol)
            vals = {"rho_sign_concordant": rs, "traj_r": tr,
                    "max_group_match": mgm, "peak_state_match": psm}
            # spec §8.1 原文 = "加权均值(权重见 config, NA 项按剩余权重重归一)"。traj_r 按**原值**(−1..1)
            # 入权, 不做 (r+1)/2 之类改标度(口径不自创, 硬约束 10)。故 fidelity 理论区间 [−0.3, 1.0]:
            # 上限 = 四指标全满分(0.3+0.3+0.2+0.2); 下限 = 仅 traj_r=−1 其余为 0。按"越高越保守"使用, 不当百分比。
            num = sum(W[k] * v for k, v in vals.items() if np.isfinite(v))
            den = sum(W[k] for k, v in vals.items() if np.isfinite(v))
            n_met = int(sum(1 for v in vals.values() if np.isfinite(v)))
            fid = (num / den) if den > 0 else float("nan")
            recs.append({"gene": gene, "group": g, "rho_human": rh, "rho_mouse": rm,
                         "status_human": sh, "status_mouse": sm,
                         "rho_sign_concordant": rs, "traj_r": tr, "n_shared_states": n_shared,
                         "max_group_match": mgm, "max_group_human": a, "max_group_mouse": b,
                         "peak_state_human": pk_h or "NA", "peak_state_mouse": pk_m or "NA",
                         "peak_state_match": psm, "fidelity": fid, "n_metrics": n_met,
                         "evidence_grade": ("high" if n_met >= 3 else "medium" if n_met == 2
                                            else "low" if n_met == 1 else "none")})
    return pd.DataFrame(recs)
