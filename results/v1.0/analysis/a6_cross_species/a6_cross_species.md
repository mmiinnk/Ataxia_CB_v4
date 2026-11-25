# A6 人鼠模块对照报告（gene_set v1.0）

> 一次性分析层（`scripts/a6_cross_species_modules.py`），生成 2026-11-24 03:17。
> **定位**：描述性对照，**不是** s04/s05（跨物种层）或 s09（保守性层）——那三层未实现，且需 spec §11-4 的 bin 级 `stage_map_mouse_to_human_rank`（6 锚点 + 6 插值）。本脚本刻意不做 bin 级对齐。

## 1. 口径

- 共有基因（两物种 A4 都入聚类）：**90**；模块划分 **NMI = 0.220**（1=完全一致，0=无关）。
- 细胞类型偏好：8 主组平均谱，Spearman；阈值 n≥8 才出值。
- 粗 4 阶段：S1 产前早期 / S2 围产 / S3 出生后 / S4 成年，映射见脚本 `STAGE_MAP`；人侧 bin 轴 = 主口径（`collapse_adult_bins`，成人合并为 AD）。
- NA 政策：不可行 组×bin 不进聚合，不插值、不填 0；相关为 pairwise-complete 并报实际格子数。

## 2. 主要结果

- `rho_group`（细胞类型偏好保守性）：n=90，中位 **+0.429**，四分位 [+0.101, +0.714]，>0.6 的基因 32 个，<0 的 15 个。
- `rho_stage`（粗阶段时序保守性）：n=90，中位 **+0.429**，>0.6 的 28 个。
- 主组（τ argmax）两物种相同的基因：**34/90**。
- 跨物种稳健核心集（敏感性全口径显著且同向，见 `sensitivity/robust_core_genes.tsv`）命中 21 个基因；其 `rho_group` 中位 **+0.714**，其余基因中位 +0.310。

## 3. 模块对应（人模块 → 最佳鼠模块，按 Fisher p）

| module_human | module_mouse | n_genes | n_human_module | n_mouse_module | overlap_frac_of_human | p_fisher |
|---|---|---|---|---|---|---|
| H7 | M5 | 25 | 29 | 54 | 0.862 | 0.000481 |
| H5 | M2 | 10 | 16 | 24 | 0.625 | 0.0009437 |
| H1 | M4 | 3 | 8 | 4 | 0.375 | 0.001825 |
| H3 | M3 | 2 | 18 | 3 | 0.111 | 0.1007 |
| H4 | M2 | 6 | 15 | 24 | 0.4 | 0.2142 |
| H8 | M2 | 1 | 1 | 24 | 1 | 0.2667 |
| H2 | M2 | 1 | 2 | 24 | 0.5 | 0.4644 |
| H6 | M5 | 1 | 1 | 54 | 1 | 1 |

## 4. 列联表（基因数）

|  | M1 | M2 | M3 | M4 | M5 | M6 | M7 |
|---|---|---|---|---|---|---|---|
| H1 | 0 | 1 | 1 | 3 | 3 | 0 | 0 |
| H2 | 0 | 1 | 0 | 0 | 1 | 0 | 0 |
| H3 | 0 | 3 | 2 | 0 | 12 | 1 | 0 |
| H4 | 1 | 6 | 0 | 0 | 8 | 0 | 0 |
| H5 | 1 | 10 | 0 | 1 | 4 | 0 | 0 |
| H6 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| H7 | 0 | 2 | 0 | 0 | 25 | 1 | 1 |
| H8 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |

## 5. 细胞类型偏好最保守 / 最不保守的基因

| gene | rho_group | rho_stage | mech_class | max_group_human | max_group_mouse | in_cross_species_robust_core |
|---|---|---|---|---|---|---|
| CACNA1G | 0.976 | 0.712 | channel_synapse | Purkinje | Purkinje | 0.000 |
| PPP2R2B | 0.976 | 0.721 | proteostasis_RNA | Astro | OligoLineage | 1.000 |
| SLC9A6 | 0.952 | 0.555 | cytoskeleton_transport_glia | DN | DN | 0.000 |
| RNF216 | 0.929 | 0.786 | proteostasis_RNA | Immune | Immune | 0.000 |
| ATCAY | 0.905 | 0.866 | cytoskeleton_transport_glia | OligoLineage | DN | 0.000 |
| ATP2B3 | 0.905 | 0.778 | channel_synapse | Interneuron | GranuleLineage | 1.000 |
| SYT14 | 0.905 | 0.803 | cytoskeleton_transport_glia | Purkinje | Purkinje | 0.000 |
| FXN | 0.881 | 0.337 | mito_metabolic | Progenitor | Progenitor | 0.000 |
| NPHP1 | 0.881 | 0.578 | cilia_congenital_dev | Astro | Astro | 1.000 |
| WWOX | 0.881 | 0.841 | other | Astro | Astro | 0.000 |
| TRPC3 | 0.833 | 0.386 | channel_synapse | Purkinje | Purkinje | 1.000 |
| CAMTA1 | 0.810 | 0.705 | cilia_congenital_dev | DN | DN | 0.000 |

| gene | rho_group | rho_stage | mech_class | max_group_human | max_group_mouse | in_cross_species_robust_core |
|---|---|---|---|---|---|---|
| SETX | -0.952 | -0.339 | dna_repair | Immune | Interneuron | 0.000 |
| STUB1 | -0.667 | -0.539 | proteostasis_RNA | Interneuron | Immune | 0.000 |
| APTX | -0.619 | -0.111 | dna_repair | Immune | Purkinje | 0.000 |
| TDP2 | -0.262 | -0.292 | dna_repair | Immune | GranuleLineage | 0.000 |
| UBA5 | -0.167 | 0.439 | proteostasis_RNA | Astro | Interneuron | 0.000 |
| PMPCA | -0.143 | 0.199 | mito_metabolic | OligoLineage | Progenitor | 0.000 |
| PNPLA6 | -0.143 | 0.278 | mito_metabolic | Immune | DN | 0.000 |
| CCDC88C | -0.119 | 0.216 | other | Progenitor | Interneuron | 1.000 |
| PCNA | -0.071 | 0.065 | dna_repair | Progenitor | Immune | 0.000 |
| XRCC1 | -0.071 | 0.410 | dna_repair | Immune | GranuleLineage | 0.000 |
| COQ2 | -0.048 | 0.602 | mito_metabolic | Immune | Astro | 0.000 |
| FAT2 | -0.048 | 0.465 | channel_synapse | GranuleLineage | GranuleLineage | 1.000 |

## 6. 按机制类看保守性（中位 rho_group）

| mech_class | n | median | q25 | q75 |
|---|---|---|---|---|
| channel_synapse | 15.000 | 0.667 | 0.429 | 0.786 |
| cilia_congenital_dev | 16.000 | 0.571 | 0.190 | 0.762 |
| cytoskeleton_transport_glia | 10.000 | 0.607 | 0.185 | 0.857 |
| dna_repair | 10.000 | 0.071 | -0.214 | 0.494 |
| mito_metabolic | 19.000 | 0.286 | 0.155 | 0.512 |
| other | 2.000 | 0.381 | 0.131 | 0.631 |
| polyQ_repeat | 6.000 | 0.190 | 0.155 | 0.476 |
| proteostasis_RNA | 12.000 | 0.381 | 0.042 | 0.708 |

## 7. 局限（引用前必读）

- 4 阶段映射是**本脚本的分析假设**，不是 spec §11-4 的 bin 级映射；换一套映射 `rho_stage` 会变，故只可作阶段级表述。
- 人侧 Purkinje/DN/UBC/Progenitor 产后无数据（右删失），这几组的跨物种相关实际只比较了产前段；人侧成人段仅 1 个 AD bin（且 2/4 供体为合并症供体，主口径未剔除）。
- 模块由无监督聚类得到，A4 已报告其共识稳定性中等（PAC≈0.48）；模块对应关系属**假设生成**。
- 两物种数据来自不同平台/注释（Sepp 人 snRNA vs Sepp 鼠），`rho_group` 含技术差异成分，不能直接当作生物学保守性强度。

## 8. 产物

- `a6_1_module_correspondence.pdf`
- `a6_1_module_correspondence.png`
- `a6_2_profile_scatter.pdf`
- `a6_2_profile_scatter.png`
- `a6_cross_species.md`
- `gene_cross_species_profiles.tsv`
- `module_contingency.tsv`
- `module_correspondence.tsv`
