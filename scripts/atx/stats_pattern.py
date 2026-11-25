"""atx.stats_pattern — temporal_pattern 标签与 peak_by_group(implementation_spec §3.8/§5)。

借鉴(Depression_CB README §7.5 已验证口径):峰 = 可行 bin 的 log_ 轨迹 argmax;
peak_internal = 峰两侧均有可行 bin(端点峰 = 删失单调, 与 window_high 单列区分——
如人 Purkinje 产后仅 NB 可行时一切晚升基因都会"富集于 perinatal", 是删失伪象非时间窗);
peak_vs_flank = 峰 − max(左侧可行 bin log 均值, 右侧可行 bin log 均值), 仅内部峰有值;
period_ratio = 峰所在 period 的 CP10K bin 均值 / 全轴 CP10K bin 均值(线性, bin 等权);
标签优先级 transient_peak > window_high > monotonic_up/down > endpoint_high > flat。

本项目差异(§3.8):
① 中间缺口(人 Interneuron 6mo/44y、Astro 11w)时"两侧有可行 bin"按**可行 bin 序列**
   判定而非相邻 bin: 峰在序列第 i 位, internal ⟺ 0 < i < n−1;
② 新增 half_max_window = CP10K ≥ 峰值 × half_max_fraction(config)的含峰最长连续区间,
   **沿可行 bin 序列**延展(可跨缺口), 供 Fig4d 泳道宽度(由数据决定, 非示意宽度);
③ peak_table 输出 peak_by_group.tsv(gene × group: peak_bin/peak_period/peak_cp/
   n_feasible_bins/half_max_lo/half_max_hi), 为 s07 模块指派与 s08 泳道的输入。

口径:
- 峰 ties 取最早 bin(np.argmax 首最大, 确定性);
- 可行行 = feasible 且 log_ 非 NA(与 stats_trend.trend_table 同口径);
- label() 要求 merged 帧含 peak_bin/peak_internal/peak_vs_flank/period_ratio/rho_bin/
  q_bin(后两者来自 trend_table——spec §5: monotonic 判据 q < mono_q 且 |rho| ≥ mono_rho);
  缺列 ValueError, 不静默跳过 monotonic;
- 标签受控词表: transient_peak / window_high / monotonic_up / monotonic_down /
  endpoint_high / flat / NA(无峰行[缺失基因或无可行 bin] → "NA" 而非 flat,
  与 specific_call 的 NA 语义一致);阈值比较均为"≥/< 含边界"(边界测试见 tests)。
确定性: 无 RNG;两次重跑逐值一致(R8)。
"""
from __future__ import annotations

import numpy as np
import pandas as pd

_PAT_COLS = ["gene", "group", "n_feasible_bins", "peak_bin", "peak_bin_order", "peak_period",
             "peak_log", "peak_cp", "peak_internal", "peak_vs_flank", "period_ratio", "status"]
_PEAK_COLS = ["gene", "group", "peak_bin", "peak_bin_order", "peak_period", "peak_cp",
              "n_feasible_bins", "half_max_lo", "half_max_hi", "status"]
_LABEL_REQUIRED = ("peak_bin", "peak_internal", "peak_vs_flank", "period_ratio",
                   "rho_bin", "q_bin")


def _peak_records(bin_table, genes, group, cfg) -> list[dict]:
    """一组×多基因的峰统计记录(patterns 与 peak_table 的共享内核, 单遍计算)。"""
    tp = cfg.stats_cfg["temporal_pattern"]
    half_frac = float(tp["half_max_fraction"])
    bt = bin_table.loc[bin_table["group"].to_numpy() == group]
    if len(bt):
        bt = bt.loc[bt["feasible"].to_numpy()].sort_values("bin_order")
    recs: list[dict] = []
    for gene in genes:
        col, ccol = f"log_{gene}", f"cp_{gene}"
        rec: dict = {"gene": gene, "group": group}
        if col not in bin_table.columns:
            rec.update(status="missing", n_feasible_bins=0)
            recs.append(rec)
            continue
        keep = bt[col].notna().to_numpy() if len(bt) else np.zeros(0, dtype=bool)
        b = bt["bin"].to_numpy()[keep]
        o = bt["bin_order"].to_numpy()[keep]
        per = bt["period"].to_numpy()[keep]
        lg = pd.to_numeric(bt.loc[keep, col]).to_numpy(dtype=float)
        cp = pd.to_numeric(bt.loc[keep, ccol]).to_numpy(dtype=float)
        n = int(b.size)
        if n == 0:
            rec.update(status="no_feasible_bin", n_feasible_bins=0)
            recs.append(rec)
            continue
        i = int(np.argmax(lg))                              # ties → 最早 bin(确定性)
        internal = bool(0 < i < n - 1)                      # 可行 bin 序列口径(§3.8 差异①)
        rec.update(status="ok", n_feasible_bins=n,
                   peak_bin=str(b[i]), peak_bin_order=float(o[i]), peak_period=str(per[i]),
                   peak_log=float(lg[i]), peak_cp=float(cp[i]), peak_internal=internal)
        if internal:
            left = float(np.mean(lg[:i]))
            right = float(np.mean(lg[i + 1:]))
            rec["peak_vs_flank"] = float(lg[i] - max(left, right))
        else:
            rec["peak_vs_flank"] = np.nan                   # 端点峰: 一侧无可行 bin
        pmask = per == per[i]
        tot = np.nanmean(cp) if np.isfinite(cp).any() else np.nan
        pk_mean = np.nanmean(cp[pmask]) if np.isfinite(cp[pmask]).any() else np.nan
        rec["period_ratio"] = (float(pk_mean / tot)
                               if np.isfinite(pk_mean) and np.isfinite(tot) and tot > 0
                               else np.nan)
        thr = half_frac * cp[i]                             # half_max_window(§3.8 差异②)
        lo = hi = i
        while lo - 1 >= 0 and np.isfinite(cp[lo - 1]) and cp[lo - 1] >= thr:
            lo -= 1
        while hi + 1 < n and np.isfinite(cp[hi + 1]) and cp[hi + 1] >= thr:
            hi += 1
        rec.update(half_max_lo=str(b[lo]), half_max_hi=str(b[hi]))
        recs.append(rec)
    return recs


def _frame(recs: list[dict], cols: list[str]) -> pd.DataFrame:
    df = pd.DataFrame(recs)
    for c in cols:
        if c not in df.columns:
            df[c] = np.nan
    return df[cols].reset_index(drop=True)


def patterns(bin_table, genes, group, cfg) -> pd.DataFrame:
    """一物种×一组的 gene 级峰/窗口统计(spec §3.8 签名;§5 temporal_patterns.tsv 列)。

    列 = spec §5(peak_bin/peak_bin_order/peak_period/peak_internal/peak_vs_flank/
    period_ratio)+ 诊断列(gene/group/n_feasible_bins/peak_log/peak_cp/status);
    temporal_pattern 由 label() 在 s02 merge trend_table 后追加。
    status: ok / missing(无 log_ 列) / no_feasible_bin(组内无可行 bin)。
    """
    return _frame(_peak_records(bin_table, genes, group, cfg), _PAT_COLS)


LOWEXPR_MASKED_PATTERNS = ("window_high", "transient_peak")


def mask_lowexpr_labels(tp: pd.DataFrame, spec_global: pd.DataFrame) -> int:
    """not_expressed 基因的窗口/峰伪标签 → temporal_pattern='NA'(原地修改, 返回掩膜行数)。

    P7 裁决①(用户确认 2026-11-05): 全组均低于 min_expr 的基因, 其"窗口高表达/瞬时峰"
    标签是接近零的噪声上的假峰(真实鼠草稿: 61 行窗口/峰标签中 15 行属此类;Purkinje 组
    仅有的 2 个窗口标签 BEAN1/TTBK2 均属此类)。与 NA 政策同构: 不可靠证据不给标签,
    图上按 NA 斜纹/灰色, 下游 s07/s08 自动排除。endpoint_high 同类残余不在本次裁决,
    s06 物种报告计数留档。掩膜行可识别: status=='ok' 且 temporal_pattern=='NA'
    (label() 对 status ok 行必给非 NA 标签)。
    """
    ne = set(spec_global.loc[spec_global["specific_call"].astype(str) == "not_expressed",
                             "gene"].astype(str))
    m = (tp["gene"].astype(str).isin(ne)
         & tp["temporal_pattern"].astype(str).isin(LOWEXPR_MASKED_PATTERNS))
    tp.loc[m, "temporal_pattern"] = "NA"
    return int(m.sum())


def label(df, cfg) -> pd.Series:
    """temporal_pattern 标签(spec §5 优先级;返回与 df 同 index 的 Series)。

    优先级(高→低): transient_peak(pvf ≥ pulse_min_peak_vs_flank 且内部峰)
    > window_high(period_ratio ≥ window_min_ratio 且内部峰)
    > monotonic_up/down(q_bin < mono_q 且 rho_bin ≥/≤ ±mono_rho)
    > endpoint_high(period_ratio ≥ window_min_ratio 但端点峰——与删失单调不可区分, 单列)
    > flat;无峰行 → "NA"。实现按低→高覆盖, 等价于优先级 cascade。
    """
    tp = cfg.stats_cfg["temporal_pattern"]
    pulse = float(tp["pulse_min_peak_vs_flank"])
    wr = float(tp["window_min_ratio"])
    mq = float(tp["mono_q"])
    mr = float(tp["mono_rho"])
    miss = [c for c in _LABEL_REQUIRED if c not in df.columns]
    if miss:
        raise ValueError(f"label() 缺列 {miss}——须先 merge trend_table 的 rho_bin/q_bin"
                         "(spec §5: monotonic 判据 q<mono_q 且 |rho|≥mono_rho)")
    has_peak = df["peak_bin"].notna().to_numpy()
    internal = df["peak_internal"].eq(True).to_numpy()      # NaN.eq(True)=False, 安全
    pvf = pd.to_numeric(df["peak_vs_flank"]).to_numpy(dtype=float)
    pr = pd.to_numeric(df["period_ratio"]).to_numpy(dtype=float)
    rho = pd.to_numeric(df["rho_bin"]).to_numpy(dtype=float)
    q = pd.to_numeric(df["q_bin"]).to_numpy(dtype=float)
    lab = np.full(len(df), "flat", dtype=object)
    lab[~has_peak] = "NA"
    lab[has_peak & ~internal & (pr >= wr)] = "endpoint_high"
    lab[has_peak & (q < mq) & (rho >= mr)] = "monotonic_up"
    lab[has_peak & (q < mq) & (rho <= -mr)] = "monotonic_down"
    lab[has_peak & internal & (pr >= wr)] = "window_high"
    lab[has_peak & internal & (pvf >= pulse)] = "transient_peak"
    return pd.Series(lab, index=df.index, name="temporal_pattern")


def peak_table(bin_table, genes, cfg) -> pd.DataFrame:
    """peak_by_group.tsv(spec §3.8 差异③): gene × group 峰表, s07/s08 输入。

    组 = config plot_groups ∩ bin_table 实有组(含 UBC/ALL——broad 基因的模块指派
    用 ALL 行, s07 §6);列 = gene, group, peak_bin, peak_bin_order, peak_period,
    peak_cp, n_feasible_bins, half_max_lo, half_max_hi, status。
    """
    present = set(bin_table["group"].unique())
    groups = [g for g in cfg.plot_groups() if g in present]
    recs: list[dict] = []
    for g in groups:
        recs.extend(_peak_records(bin_table, genes, g, cfg))
    return _frame(recs, _PEAK_COLS)
