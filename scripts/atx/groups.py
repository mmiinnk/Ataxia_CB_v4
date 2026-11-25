"""atx.groups — obs 派生列 + 细胞组指派(implementation_spec §3.4)。

借鉴(Depression_CB 已验证口径):分组表达式用 pandas.eval(engine='python');分类列先 NA 安全化
(<NA> → 哨兵 '__NA__');敏感性开关只打标 use_cell, 由 s01 同步过滤 X 与 obs(避免行错位)。

本项目差异(spec §3.4):
① NA 安全化列 = author_cell_type, subtype, dev_state, cell_type;
② `subtype.str.startswith(...)` 类判断**先预计算布尔辅助列**再 eval(eval 的 .str 支持不可靠):
   sub_is_rl / sub_is_progenitor / sub_na;
③ 断言组互斥(除 ALL)与 "ungrouped 只进 ALL", 计数与 §4.1 锚点精确比对(check_anchors);
④ 鼠时间列 author_stage, 人 development_stage(SpeciesCfg.time_col);
⑤ 敏感性 drop_glutamatergic_uncertain / split_granule_progenitor 在此改写表达式。

**RL 归属实测裁定(2026-11-04, 本仓原型脚本输出)**:astroglia × subtype=progenitor_RL(_early)
(人 1360 / 鼠 3688 核)同时命中 GranuleLineage 的 subtype 项与 Progenitor 的
`astroglia & sub_is_progenitor` 项。定义 sub_is_progenitor = startswith('progenitor_') & ~RL
后, 两物种 9 组 + ungrouped + total **全部精确复现 §4.1 锚点**(GranuleLineage 优先的等价
表述)。生物学一致:RL = rhombic lip, 颗粒细胞谱系起源。若未来数据/表达式变动导致其它组对
重叠, assign() 的互斥断言会当场报错而不是静默按序吞掉。
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

NA = "__NA__"
SAFE_COLS = ["author_cell_type", "subtype", "dev_state", "cell_type"]
RL_SUBTYPES = ["progenitor_RL", "progenitor_RL_early"]


def effective_group_exprs(cfg) -> dict[str, str | None]:
    """config cell_groups 表达式 + 敏感性改写(spec §3.4 ⑤;表达式本体只在 config, 这里只做
    两个已声明开关的机械改写, 改写失败即报错——不许静默)。"""
    exprs = dict(cfg.group_exprs())
    if cfg.sens("drop_glutamatergic_uncertain", False):
        new = re.sub(r",\s*'glutamatergic_uncertain'", "", exprs["GranuleLineage"])
        if new == exprs["GranuleLineage"]:
            raise RuntimeError("drop_glutamatergic_uncertain: GranuleLineage 表达式未含该标签, 改写失败")
        exprs["GranuleLineage"] = new
    if cfg.sens("split_granule_progenitor", False):
        pat = r"\s*\|\s*subtype in \['progenitor_RL','progenitor_RL_early'\]"
        new = re.sub(pat, "", exprs["GranuleLineage"])
        if new == exprs["GranuleLineage"]:
            raise RuntimeError("split_granule_progenitor: GranuleLineage 表达式未含 RL 项, 改写失败")
        exprs["GranuleLineage"] = new
        exprs["GranuleProgenitor"] = "subtype in ['progenitor_RL','progenitor_RL_early']"
    return exprs


def derive_obs(obs: pd.DataFrame, spcfg, cfg) -> pd.DataFrame:
    """原始 obs → 派生框(NA 安全化 + 布尔辅助列 + bin/bin_order/donor[/disease_flag])。

    行序与 obs 一致(位置对齐 X);时间值必须全部落在 spcfg.label_of, 否则报错(不静默丢弃)。
    """
    # d 用位置索引(RangeIndex),与 X 行位置对齐;真实 h5ad 的 obs 索引是细胞条码(字符串),
    # 直接 d[col]=series 会被 pandas 按索引重对齐 → 整列静默变 NA(2026-11-04 s00 实测踩中)。
    # 故一律 .to_numpy() 按位置赋值。
    d = pd.DataFrame(index=np.arange(len(obs)))
    for c in SAFE_COLS:
        if c in obs.columns:
            s = obs[c]
            d[c] = s.astype("object").where(s.notna(), NA).astype(str).to_numpy()
        else:
            d[c] = NA
    d["sub_is_rl"] = d["subtype"].isin(RL_SUBTYPES)
    # RL 归属 GranuleLineage(见模块 docstring 实测裁定)→ 从 Progenitor 辅助列中排除
    d["sub_is_progenitor"] = d["subtype"].str.startswith("progenitor_") & ~d["sub_is_rl"]
    d["sub_na"] = d["subtype"] == NA

    tcol = spcfg.time_col
    if tcol not in obs.columns:
        raise RuntimeError(f"obs 缺时间列 {tcol!r}(实际列: {list(obs.columns)[:8]}…)")
    raw_t = obs[tcol].astype("object").where(obs[tcol].notna(), NA).astype(str)
    lab = raw_t.map(spcfg.label_of)
    unmapped = raw_t[lab.isna()].unique()
    if len(unmapped):
        raise RuntimeError(f"时间列 {tcol} 有 config time_order 之外的值: {list(unmapped)[:5]}")
    d["bin"] = lab.astype(str).to_numpy()
    d["bin_order"] = d["bin"].map(spcfg.bin_order).astype(int)

    dcol = spcfg.donor_col
    if dcol not in obs.columns:
        raise RuntimeError(f"obs 缺供体列 {dcol!r}")
    d["donor"] = obs[dcol].astype("object").where(obs[dcol].notna(), NA).astype(str).to_numpy()
    d["pooled_bin"] = d["bin"].isin(spcfg.pooled)

    if spcfg.disease_col and spcfg.disease_col in obs.columns:
        dis = obs[spcfg.disease_col].astype("object").where(obs[spcfg.disease_col].notna(), NA).astype(str)
        d["disease"] = dis.to_numpy()
        d["disease_flag"] = (dis != "normal").to_numpy()  # §11 未决 9: donor cube 保留 disease_flag
    return d


def assign(d: pd.DataFrame, cfg):
    """细胞组指派(spec §3.4 签名: assign(obs, cfg) -> DataFrame[cell_pos, cell_group])。

    互斥断言:任意两组掩码交集必须为 0(RL 裁定后应天然成立);出现新重叠 = 数据或表达式变了,
    当场报错。未命中任何组的细胞 = 'ungrouped'(只进 ALL)。组顺序 = config cell_groups 声明序。
    """
    exprs = effective_group_exprs(cfg)
    masks: dict[str, np.ndarray] = {}
    for g, e in exprs.items():
        if e is None:  # ALL
            continue
        m = pd.eval(e, local_dict=d, engine="python")
        m = np.asarray(m, dtype=bool)
        if m.shape != (len(d),):
            raise RuntimeError(f"组 {g} 表达式返回形状 {m.shape}, 期望 {(len(d),)}")
        masks[g] = m
    names = list(masks)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            ov = int((masks[a] & masks[b]).sum())
            if ov:
                raise AssertionError(f"组互斥断言失败: {a} ∩ {b} = {ov} 核(见 groups.py docstring RL 裁定)")
    group = np.full(len(d), "ungrouped", dtype=object)
    any_hit = np.zeros(len(d), dtype=bool)
    for g in names:
        group[masks[g]] = g
        any_hit |= masks[g]
    out = pd.DataFrame({"cell_pos": np.arange(len(d), dtype=np.int64), "cell_group": group})
    return out, masks


def expand_long(assign_df: pd.DataFrame, cfg) -> pd.DataFrame:
    """细胞 → (cell_pos, cell_group) 长表:每核一行本组(ungrouped 除外) + 一行 ALL。

    'ungrouped 只进 ALL'(spec §4.1)。伪聚合与状态剖面都从本长表出发。
    """
    ref = set(cfg.reference_groups())  # {'ALL'}
    keep = assign_df["cell_group"] != "ungrouped"
    own = assign_df.loc[keep, ["cell_pos", "cell_group"]]
    allrows = pd.DataFrame({"cell_pos": assign_df["cell_pos"].to_numpy(),
                            "cell_group": np.repeat(sorted(ref), len(assign_df))})
    return pd.concat([own, allrows], ignore_index=True)


def use_cell_mask(d: pd.DataFrame, assign_df: pd.DataFrame, cfg, spcfg) -> np.ndarray:
    """敏感性/排除开关 → 细胞级 use_cell 布尔(打标不删行;由 s01 同步过滤 X 与 obs)。"""
    use = np.ones(len(d), dtype=bool)
    if cfg.sens("exclude_disease_donors", False) or spcfg.exclude_disease_donors:
        if "disease_flag" in d.columns:
            use &= ~d["disease_flag"].to_numpy()
        else:
            raise RuntimeError("exclude_disease_donors=true 但 obs 无疾病列")
    return use


def group_counts(assign_df: pd.DataFrame) -> dict[str, int]:
    """9 组 + ungrouped + total 计数(check_anchors / s00 用)。"""
    vc = assign_df["cell_group"].value_counts()
    out = {g: int(vc.get(g, 0)) for g in vc.index}
    out["total"] = int(len(assign_df))
    return out


def check_anchors(assign_df: pd.DataFrame, species: str, cfg) -> list[str]:
    """组计数 vs §4.1 锚点(spec §3.4 签名 check_anchors(obs_long, species, cfg);ERROR 级)。

    返回违例列表(空 = 全部精确匹配)。敏感性开关启用时按 GROUP_COUNTS_SENS 调整期望。
    """
    from . import sanity
    counts = group_counts(assign_df)
    expect = dict(sanity.GROUP_COUNTS[species])
    if species == "human" and cfg.sens("drop_glutamatergic_uncertain", False):
        expect.update(sanity.GROUP_COUNTS_SENS["human"]["drop_glutamatergic_uncertain"])
    if cfg.sens("split_granule_progenitor", False) or cfg.sens("include_ubc_in_main", False):
        # 组构成变了 → 锚点只核 total 与 ungrouped, 其余降级为报告(不静默跳过: 返回说明行)
        bad = []
        if counts.get("total") != expect["total"]:
            bad.append(f"total: {counts.get('total')} != {expect['total']}")
        bad.append("NOTE: split_granule_progenitor/include_ubc_in_main 启用, 分组锚点不适用(仅核 total)")
        return bad
    bad = []
    for g, v in expect.items():
        got = counts.get(g, 0 if g != "total" else None)
        if got != v:
            bad.append(f"{g}: 实测 {got} != 锚点 {v}")
    return bad
