"""atx.sanity — 实测锚点字典与校验(implementation_spec §3.12/§4)。

锚点全部来自 2026-10-14/10-28 对原始数据的实测(spec §4.1/§4.2/§4.3/§3.12),不是设计估计。
级别约定:组计数锚点(§4.1)与形状锚点 = ERROR(断言失败即停);组×bin 可行性矩阵(§4.2)与
趋势方向锚点 = WARNING(spec §3.12);dropout 复验 = WARNING(现场重算并与锚点比对)。

marker 期望组从 genes/sanity_markers.tsv 读(不硬编码, 硬约束 6 只约束参数;锚点是实测事实,
按 spec §3.12 固化在本模块)。
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# ---- §4.1 9 组细胞数(互斥;ungrouped 只进 ALL)。ERROR 级。----
GROUP_COUNTS = {
    "human": {
        "Purkinje": 45341, "GranuleLineage": 72112, "DN": 6352, "Interneuron": 6805,
        "UBC": 3847, "OligoLineage": 11064, "Astro": 3150, "Progenitor": 25543,
        "Immune": 1770, "ungrouped": 4972, "total": 180956,
    },
    "mouse": {
        "Purkinje": 18271, "GranuleLineage": 35247, "DN": 6656, "Interneuron": 5286,
        "UBC": 1208, "OligoLineage": 881, "Astro": 1843, "Progenitor": 36419,
        "Immune": 371, "ungrouped": 9100, "total": 115282,
    },
}

# 敏感性 drop_glutamatergic_uncertain=true 时的人 GranuleLineage(§4.1 注)
GROUP_COUNTS_SENS = {"human": {"drop_glutamatergic_uncertain": {"GranuleLineage": 71005}}}

# ---- §4.2 组 × bin 细胞数(min_cells=30 可行性矩阵的原始计数)。WARNING 级。----
HUMAN_BINS = ["CS18", "CS19", "CS22", "9w", "11w", "17w", "20w", "NB",
              "6mo", "7mo", "9mo", "3y", "42y", "44y", "46y", "52y"]
MOUSE_BINS = ["E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5",
              "E17.5", "P0", "P4", "P7", "P14", "9wk"]

FEAS_MATRIX = {
    "human": {  # 列 = HUMAN_BINS
        "Purkinje":       [80, 656, 25993, 9300, 6136, 1887, 1231, 53, 0, 0, 0, 4, 1, 0, 0, 0],
        "GranuleLineage": [354, 206, 875, 762, 2155, 3144, 2198, 13187, 2969, 4606, 4070, 13432, 14102, 496, 6214, 3342],
        "DN":             [58, 463, 3807, 1163, 525, 336, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "Interneuron":    [2, 6, 282, 812, 1220, 1594, 933, 930, 17, 78, 62, 87, 214, 27, 399, 142],
        "UBC":            [0, 0, 51, 95, 1279, 1920, 444, 58, 0, 0, 0, 0, 0, 0, 0, 0],
        "OligoLineage":   [0, 0, 2, 0, 73, 289, 124, 208, 15, 37, 71, 398, 466, 3209, 908, 5264],
        "Astro":          [0, 0, 2, 0, 8, 453, 73, 585, 134, 132, 199, 161, 583, 167, 392, 261],
        "Progenitor":     [4832, 4436, 6851, 1301, 1658, 3251, 2582, 619, 5, 3, 3, 2, 0, 0, 0, 0],
        "Immune":         [6, 9, 166, 45, 47, 110, 46, 142, 3, 11, 8, 46, 115, 675, 49, 292],
    },
    "mouse": {  # 列 = MOUSE_BINS
        "Purkinje":       [21, 38, 2420, 7664, 2272, 2756, 1517, 1111, 231, 190, 30, 21],
        "GranuleLineage": [1477, 876, 688, 1461, 657, 1474, 2727, 4107, 4692, 4176, 4781, 8131],
        "DN":             [16, 101, 1597, 2557, 825, 741, 437, 284, 37, 45, 9, 7],
        "Interneuron":    [38, 60, 211, 833, 311, 552, 698, 759, 560, 498, 400, 366],
        "UBC":            [0, 0, 0, 0, 27, 102, 311, 536, 115, 86, 17, 14],
        "OligoLineage":   [6, 4, 0, 1, 6, 6, 90, 158, 155, 155, 149, 151],
        "Astro":          [0, 0, 0, 13, 70, 76, 245, 458, 403, 307, 125, 146],
        "Progenitor":     [7595, 10634, 4760, 5656, 975, 1574, 1487, 1732, 1013, 896, 67, 30],
        "Immune":         [56, 41, 14, 27, 6, 37, 16, 41, 41, 52, 21, 19],
    },
}
FEAS_BINS = {"human": HUMAN_BINS, "mouse": MOUSE_BINS}

# ---- §3.12 鼠 dropout 锚点(2026-10-14 原始 mouse.h5ad 实测; s00 现场重算比对)。WARNING 级。----
# 全库 pct+(%) 低于 1% → 豁免(不作为鼠侧 marker 失败依据)
MOUSE_DROPOUT_EXEMPT = {"Neurod1": 0.66, "Ptprc": 0.07, "Aif1": 0.06, "Csf1r": 0.08, "Tmem119": 0.05}
# 鼠免疫锚点改用这三个(免疫组内 pct+ %)
MOUSE_IMMUNE_ANCHORS = {"Cx3cr1": 52.0, "P2ry12": 42.0, "Hexb": 50.0}
# 人侧免疫 marker 组内 pct+ 正常范围(%, PTPRC/CSF1R/CX3CR1/P2RY12)
HUMAN_IMMUNE_PCT_RANGE = (59.0, 68.0)

# ---- §4.3 外部数据集形状锚点。ERROR 级。----
EXTERNAL_SHAPES = {"aldinger": (69174, 26326), "castonguay": (1004112, 47789)}

# Castonguay 供体锚点(§4.3, normal 供体)。WARNING 级。
CASTONGUAY_ANCHORS = {"n_donors_normal": 93, "n_pc_normal": 31089,
                      "n_donors_ge30pc": 89, "median_pc_ge30": 257}

# ---- size_factor 锚点(§2 normalization 注释/§3.12):与 X 行和 r=0.99997, 行和约低 1%。WARNING 级。----
SIZE_FACTOR_ANCHOR = {"corr_min": 0.999, "rowsum_over_sf_median_range": (0.95, 1.02)}


def check_group_counts(counts: dict[str, int], species: str) -> list[str]:
    """组计数 vs §4.1 锚点(ERROR 级)。counts 键 = 9 组 + ungrouped + total。返回违例列表(空=通过)。"""
    anchor = GROUP_COUNTS[species]
    bad = []
    for g, v in anchor.items():
        got = counts.get(g)
        if got is None:
            bad.append(f"{g}: 缺失(期望 {v})")
        elif int(got) != int(v):
            bad.append(f"{g}: 实测 {got} != 锚点 {v}")
    return bad


def check_feas_matrix(n_cells: pd.DataFrame, species: str) -> pd.DataFrame:
    """组×bin 细胞数 vs §4.2 锚点(WARNING 级)。

    n_cells: index=group, columns=bin 标签, 值=n_cells_total(来自 bin_table)。
    返回差异表 [group, bin, observed, anchor, diff](空=全匹配)。只比 9 主组(锚点无 UBC/ALL 行
    ——§4.2 矩阵含 UBC 行时才比;本表含 UBC)。
    """
    anchor = FEAS_MATRIX[species]
    bins = FEAS_BINS[species]
    rows = []
    for g, arow in anchor.items():
        if g not in n_cells.index:
            rows.append({"group": g, "bin": "<row>", "observed": np.nan, "anchor": sum(arow), "diff": "缺行"})
            continue
        for b, a in zip(bins, arow):
            obs = n_cells.loc[g, b] if b in n_cells.columns else np.nan
            obs = 0 if pd.isna(obs) else int(obs)
            if obs != int(a):
                rows.append({"group": g, "bin": b, "observed": obs, "anchor": int(a), "diff": obs - int(a)})
    return pd.DataFrame(rows, columns=["group", "bin", "observed", "anchor", "diff"])


def load_markers(path: str) -> pd.DataFrame:
    """genes/sanity_markers.tsv → DataFrame(gene_symbol, cell_group, species_valid, is_also_panel, …)。"""
    m = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=["NA", ""])
    for c in ("gene_symbol", "cell_group", "species_valid"):
        if c not in m.columns:
            raise ValueError(f"sanity_markers.tsv 缺列 {c}")
    return m


def marker_table(X_cp, cell_group: pd.Series, gene_cols: dict[str, int],
                 markers: pd.DataFrame, species: str, main_groups: list[str],
                 descriptive_groups: list[str] | None = None) -> pd.DataFrame:
    """marker 校验(spec §9 review 门第 3 条):每个 marker 的期望组必须是评估组内均值 CP10K 的 argmax。

    X_cp: 细胞×基因 CSR(线性 CP10K); cell_group: 每细胞组标签(含 'ungrouped');
    gene_cols: 基因符号 → X_cp 列号; markers: load_markers 输出; 按 species_valid 过滤。
    返回逐 marker 结果表(verdict=PASS/FAIL/ABSENT);FAIL 定义 = 期望组非 argmax 或期望组均值 < argmax 的一半。
    (鼠缺 Ca8/Sox2 等 → ABSENT, 依 species_valid 与实测 var, 不算失败, CL-018。)

    评估组口径: 期望组属主组的 marker 仍只在主组内判 argmax(保持既有门 verdict 不变);
    期望组属描述组(如 UBC)的 marker 在 主组+描述组 内判 —— 补上此前 UBC marker 因期望组
    非主组而永久 ABSENT、细胞类型校验存在盲区的缺陷。ABSENT 区分 gene_not_in_var / expected_group_empty。
    """
    import scipy.sparse as sp

    main = list(main_groups)
    groups = main + [g for g in (descriptive_groups or []) if g not in main]
    mask = {g: (cell_group == g).to_numpy() for g in groups}
    n = {g: int(m.sum()) for g, m in mask.items()}
    valid = markers[markers["species_valid"].isin(["both", species])]
    out = []
    for _, r in valid.iterrows():
        sym = r["gene_symbol"]
        exp_g = r["cell_group"]
        if sym not in gene_cols or n.get(exp_g, 0) == 0:
            out.append({"marker": sym, "expected_group": exp_g, "argmax_group": "",
                        "mean_expected": np.nan, "mean_argmax": np.nan, "verdict": "ABSENT",
                        "absent_reason": ("gene_not_in_var" if sym not in gene_cols
                                          else "expected_group_empty")})
            continue
        eval_groups = main if exp_g in main else groups
        col = X_cp[:, gene_cols[sym]].toarray().ravel().astype(np.float64)
        means = {g: float(col[mask[g]].mean()) for g in eval_groups if n[g] > 0}
        amax = max(means, key=means.get)
        ok = amax == exp_g or means[exp_g] >= 0.5 * means[amax]
        out.append({"marker": sym, "expected_group": exp_g, "argmax_group": amax,
                    "mean_expected": means[exp_g], "mean_argmax": means[amax],
                    "verdict": "PASS" if ok else "FAIL",
                    **{f"mean_{g}": means.get(g, np.nan) for g in eval_groups}})
    return pd.DataFrame(out)
