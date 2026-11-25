"""atx.genes — panel 装载与数据集 var 映射(implementation_spec §3.2)。

主键 = ENSG/ENSMUSG(四个数据集 var.index 都是 Ensembl ID),符号只是回退;别名表
genes/gene_aliases.tsv 作第三级回退(如鼠 Mre11a→MRE11)。映射结果固化为 GeneMap 表并写明
match_method(借鉴 Depression_CB:映射必须可审计;panel 基因不得落在重复符号上)。

Castonguay 特例(§11 未决 5):只按 ENSG 匹配, ENSG 缺失记 status=missing, **不得**回退符号后
静默使用(ensg_only=True)。

tier 口径(硬约束/CL-023):统计 panel = tier ∈ config tiers_main(默认 [core]) 且 in_panel=true;
sensitivity.include_extended_tier=true 时加 effective tier=extended(peripheral 恒排除);
sanity marker 并入映射(is_sanity=True)但不进统计 panel(除非 is_also_panel)。
"""
from __future__ import annotations

import pandas as pd

STAT_COLS = ["gene_key", "hgnc_symbol", "mgi_symbol", "ensg_id", "ensmusg_id", "tier",
             "mech_class", "onset_typical_y", "onset_class", "in_stat_panel", "is_sanity"]


def _read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=["NA", ""])


def _bool(s: pd.Series) -> pd.Series:
    return s.astype("string").str.strip().str.lower().isin(["true", "1", "yes"]).fillna(False).astype(bool)


def load_panel(cfg) -> pd.DataFrame:
    """ataxia_genes.tsv(+sanity_markers.tsv)→ panel 表(spec §3.2 签名)。

    行 = 统计 panel(tier∈tiers_main & in_panel) ∪ sanity-only marker(tier='sanity_only')。
    列含 STAT_COLS;onset_typical_y 转 float(NA 保留)。
    """
    t = _read_tsv(cfg.paths["genes_table"])
    tiers = list(cfg.gene_set["tiers_main"])
    if cfg.sens("include_extended_tier", False) and "extended" not in tiers:
        tiers = tiers + ["extended"]
    in_panel = _bool(t["in_panel"])
    stat = t["tier"].isin(tiers) & in_panel
    panel = t.loc[stat].copy()
    panel["in_stat_panel"] = True
    panel["is_sanity"] = False

    markers = pd.read_csv(cfg.paths["sanity_table"], sep="\t", dtype=str,
                          keep_default_na=False, na_values=["NA", ""])
    have = set(panel["hgnc_symbol"])
    extra = markers[~markers["gene_symbol"].isin(have)].copy()
    if len(extra):
        rows = pd.DataFrame({
            "hgnc_symbol": extra["gene_symbol"],
            "mgi_symbol": extra["mgi_symbol"],
            "ensg_id": extra["ensg_id"],
            "ensmusg_id": extra["ensmusg_id"],
            "tier": "sanity_only",
            "in_stat_panel": False,
            "is_sanity": True,
        })
        panel = pd.concat([panel, rows], ignore_index=True, sort=False)
    # panel 内已存在的 marker 基因(如 CA8)双标 is_sanity
    panel.loc[panel["hgnc_symbol"].isin(set(markers["gene_symbol"])), "is_sanity"] = True

    panel["gene_key"] = panel["hgnc_symbol"]
    if panel["gene_key"].duplicated().any():
        dup = panel.loc[panel["gene_key"].duplicated(), "gene_key"].tolist()
        raise ValueError(f"panel gene_key 重复(符号冲突, 须先走 aliases/s00a 流程): {dup}")
    panel["onset_typical_y"] = pd.to_numeric(panel["onset_typical_y"], errors="coerce")
    for c in STAT_COLS:
        if c not in panel.columns:
            panel[c] = pd.NA
    return panel[STAT_COLS].reset_index(drop=True)


def load_aliases(cfg) -> pd.DataFrame:
    return _read_tsv(cfg.paths["aliases_table"])


def map_to_var(panel: pd.DataFrame, var: pd.DataFrame, species: str, dataset: str | None = None,
               aliases: pd.DataFrame | None = None, id_col: str = "index",
               symbol_col: str = "feature_name", ensg_only: bool = False) -> pd.DataFrame:
    """panel × 数据集 var → GeneMap 行(spec §3.2 签名;dataset 缺省 = species)。

    匹配级联:① Ensembl 主键(人及外部人类数据 ensg_id;鼠 ensmusg_id);② 符号(人 hgnc_symbol
    大写比对;鼠 mgi_symbol, 兼容首字母大写形态);③ 别名表。ensg_only=True(Castonguay)只走 ①。
    status ∈ matched / ambiguous(符号多命中, 不匹配) / missing;match_method ∈ ensembl / symbol / alias。
    返回列: STAT_COLS + row_<ds>, status_<ds>, match_method_<ds>。
    """
    dataset = dataset or species
    key_col = "ensmusg_id" if species == "mouse" else "ensg_id"
    sym_col = "mgi_symbol" if species == "mouse" else "hgnc_symbol"

    vidx = pd.Index(var.index) if id_col == "index" else pd.Index(var[id_col].astype(str))
    if vidx.duplicated().any():
        dups = vidx[vidx.duplicated()].unique()[:5].tolist()
        raise ValueError(f"{dataset}: var {id_col} 有重复 Ensembl ID(前5): {dups}")
    id_pos = {k: i for i, k in enumerate(vidx)}

    vsym = None
    sym_pos: dict[str, list[int]] = {}
    if symbol_col and symbol_col in getattr(var, "columns", []):
        vsym = var[symbol_col].astype(str)
        up = vsym.str.upper() if species != "mouse" else vsym
        for i, s in enumerate(up):
            sym_pos.setdefault(s, []).append(i)

    # 别名: alias 符号 → 本 panel 基因(按 canonical ensg/ensmusg/symbol 关联)
    alias_map: dict[str, str] = {}
    if aliases is not None and len(aliases):
        id2key = {}
        for _, p in panel.iterrows():
            id2key[str(p["ensg_id"])] = p["gene_key"]
            id2key[str(p["ensmusg_id"])] = p["gene_key"]
            id2key[str(p["hgnc_symbol"]).upper()] = p["gene_key"]
            if pd.notna(p["mgi_symbol"]):
                id2key[str(p["mgi_symbol"])] = p["gene_key"]
        for _, a in aliases.iterrows():
            for ref in (a.get("ensg_id"), a.get("ensmusg_id"), a.get("hgnc_symbol"), a.get("mgi_symbol")):
                if pd.notna(ref) and str(ref) in id2key:
                    # 人 sym_pos 键为大写、鼠为原样(如 Mre11a) → 原形与大写都登记
                    alias_map[str(a["alias"])] = id2key[str(ref)]
                    alias_map[str(a["alias"]).upper()] = id2key[str(ref)]
                    break

    rows = []
    for _, p in panel.iterrows():
        # method 用 None 而非 pd.NA:pandas 3 中 `NA != "alias"` 返回 NA, 进布尔上下文即
        # TypeError: boolean value of NA is ambiguous(实测踩中, Depression_CB 同款陷阱)
        row, status, method = pd.NA, "missing", None
        key = p[key_col]
        if pd.notna(key) and str(key) in id_pos:
            row, status, method = id_pos[str(key)], "matched", "ensembl"
        elif not ensg_only:
            cands = []
            # 候选符号按固定顺序尝试(本列符号 → HGNC 大写);set 会引入不确定序, 禁用
            for s in ([str(p[sym_col])] if pd.notna(p[sym_col]) else []) + \
                     [str(p["hgnc_symbol"]).upper()]:
                if s and s != "nan" and s in sym_pos:
                    cands = sym_pos[s]
                    method = "symbol"
                    break
            if not cands:  # 别名回退
                for al, canon in alias_map.items():
                    if canon == p["gene_key"] and al in sym_pos:
                        cands = sym_pos[al]
                        method = "alias"
                        break
            if len(cands) == 1:
                row, status = cands[0], "matched"
            elif len(cands) > 1:
                status, method = "ambiguous", None
        out = {c: p[c] for c in STAT_COLS}
        out.update({f"row_{dataset}": row, f"status_{dataset}": status,
                    f"match_method_{dataset}": method if method is not None else pd.NA})
        rows.append(out)
    return pd.DataFrame(rows)


def merge_maps(*gms: pd.DataFrame) -> pd.DataFrame:
    """多数据集 GeneMap 按 gene_key 合并(s00 汇总用)。"""
    out = gms[0]
    for gm in gms[1:]:
        out = out.merge(gm.drop(columns=[c for c in STAT_COLS if c != "gene_key"]),
                        on="gene_key", how="outer")
    return out


def hit_rate(gm: pd.DataFrame, dataset: str, stat_panel_only: bool = True) -> float:
    """命中率 = matched / 应命中行数(spec §3.2;分母默认 = 统计 panel, sanity-only 不计入门)。"""
    sub = gm[gm["in_stat_panel"].astype(bool)] if stat_panel_only else gm
    if len(sub) == 0:
        return float("nan")
    return float((sub[f"status_{dataset}"] == "matched").mean())


def matched_rows(gm: pd.DataFrame, dataset: str) -> pd.DataFrame:
    return gm[gm[f"status_{dataset}"] == "matched"].copy()
