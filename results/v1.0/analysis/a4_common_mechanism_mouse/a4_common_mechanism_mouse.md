# A4 共性机制筛选报告 — mouse（gene_set v1.0）

> 一次性分析层（`scripts/a4_common_mechanism.py`，非 sNN 契约），生成 2026-11-24 01:40，seed=20261013，n_reps=200，n_perm=10000，单供体 bin 保留(与 s02 同口径)。

## 1. 口径与输入

- 输入：`results/v1.0/intermediate/mouse/cube_bin.tsv`（可行 组×bin）、`tables/mouse/specificity_global.tsv`、`trend_by_group.tsv`、`temporal_patterns.tsv`、`genes/ataxia_genes.tsv`。
- 空间维 = 8 个主细胞类型组（Purkinje, GranuleLineage, DN, Interneuron, OligoLineage, Astro, Progenitor, Immune）；时间维 = 发育 bin；UBC（descriptive）与 ALL（reference）不入聚类。
- 聚类矩阵：93 基因 × 75 个可行 组×bin 列（单供体 11、pooled 10）。
- **NA 政策（硬约束 5）**：不可行 组×bin 为 NA，不插值、不填 0；相关用成对完整观测（min_periods=20），相关缺失对 0（按距离 1 中性处理）。
- 排除：`not_expressed` 8（BEAN1, COQ8A, MME, SPG7, TTBK2, TTPA, TWNK, WDR81）、零方差 0（无）、观测不足 0（无）；cube 未命中 6（ATN1, CA8, KIAA0586, PDYN, SAMD9L, TCTN2）。
- **CL-027**：机制类富集排除残差兜底类 ['other', 'unknown'] → 剔除 2 基因（CCDC88C, WWOX），这些基因改做单基因描述，不产生「残差类富集」结论。

## 2. 共识聚类与 k 选择

| k | PAC | min_module_size | mean_within_consensus |
|---|---|---|---|
| 3 | 0.5503 | 1 | 0.5615 |
| 4 | 0.5136 | 1 | 0.6123 |
| 5 | 0.5068 | 1 | 0.6869 |
| 6 | 0.4972 | 1 | 0.6828 |
| 7 | 0.4834 | 1 | 0.6834 |
| 8 | 0.4834 | 1 | 0.7153 |

选定 **k=7**（PAC 最小）；模块内平均共识见下表 `mean_within_consensus`（越接近 1 越稳定）。

## 3. 模块刻画（共享时空程序）

| module | n_genes | direction | peak_group | peak_bin | peak_period | top_groups | hubs | mech_composition | mean_within_consensus |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | falling | Immune | E11.5 | embryonic | Progenitor(+0.53),Interneuron(+0.31),Immune(+0.18) | ANO10,PTF1A | channel_synapse:1,cilia_congenital_dev:1 | 0.8517 |
| 2 | 24 | flat_or_peak | GranuleLineage | E10.5 | embryonic | GranuleLineage(+0.20),Progenitor(+0.10),Purkinje(+0.06) | CSPP1,ATXN10,RFC1,CWF19L1,EEF2 | dna_repair:7,proteostasis_RNA:6,mito_metabolic:5,cilia_congenital_dev:4,other:1,polyQ_repeat:1 | 0.7102 |
| 3 | 3 | flat_or_peak | Immune | E15.5 | embryonic | Immune(+1.29),DN(+0.82),Progenitor(-0.09) | OPHN1,ATM,ZFHX3 | cytoskeleton_transport_glia:1,dna_repair:1,cilia_congenital_dev:1 | 0.5541 |
| 4 | 4 | flat_or_peak | OligoLineage | E17.5 | perinatal | Astro(+1.22),OligoLineage(+0.78),Progenitor(+0.11) | PPP2R2B,GRID2,SLC1A3,POLG | channel_synapse:2,proteostasis_RNA:1,mito_metabolic:1 | 0.6513 |
| 5 | 57 | rising | Purkinje | P14 | postnatal | Purkinje(+0.37),DN(+0.27),Interneuron(+0.05) | GRM1,CAMTA1,CACNB4,KCNC3,VPS13D | channel_synapse:14,mito_metabolic:13,cilia_congenital_dev:9,cytoskeleton_transport_glia:8,polyQ_repeat:5,proteostasis_RNA:5,dna_repair:2,other:1 | 0.6815 |
| 6 | 2 | rising | Immune | P7 | postnatal | Immune(+0.96),Purkinje(+0.16),Astro(+0.16) | KIF1C,STUB1 | cytoskeleton_transport_glia:1,proteostasis_RNA:1 | 0.6515 |
| 7 | 1 | rising | GranuleLineage | E17.5 | perinatal | GranuleLineage(+1.79),Progenitor(-0.11),DN(-0.29) | RELN | cilia_congenital_dev:1 | 1.0 |

`hubs` = 模块内平均 |corr| 最高的前 5 个基因（共性程序的代表节点，后续网络层 n00–n04 的候选入口）。

## 4. 机制类富集（mech_class × module / × shape）

显著格（q<0.05 且 obs≥3）共 1 个，按 log2FE 降序前 20：

| test | a | b | obs | exp | log2fe | p | q |
|---|---|---|---|---|---|---|---|
| mech_class_x_shape | channel_synapse | flat | 4 | 10.84 | -1.333 | 0.0001 | 0.0042 |

## 5. 跨组趋势一致性（共性 vs 异质性）

| concordance | count |
|---|---|
| ns | 45 |
| single_group | 28 |
| partial_concordant | 20 |
| mixed_direction | 5 |
| pan_up | 3 |

- **pan_up**（3）：FAT2, GRM1, SYNE1

- **mixed_direction**（5）：TRPC3, MTTP, ELOVL5, FMR1, NPHP1

判读：`pan_*` = ≥5 个主组同向显著（跨细胞类型的**共性**发育程序）；`mixed_direction` = 同时存在升与降（**异质性**的细胞类型分歧）；`single_group` = 仅一组显著（细胞类型特异）。

## 6. 局限（写报告时必须一并说明）

- 单物种（鼠）、正常发育数据，**无疾病对照** → 只能说「致病基因在正常小脑发育中的共享时空程序」，不能说「共性致病机制」已被证明；后者需要 s04/s05 跨物种 + 疾病投影层（b00/b01、n00–n04）。
- 模块为无监督发现，属**假设生成**；k 由 PAC 选定，需做 k±1 与 min_cells 敏感性后才可写入正文。
- 单供体/pooled bin 参与聚类（与 s02 同口径），早期胚胎段的模块稳定性弱于晚期。
- 富集样本量小（入富集 91 基因 × k=7 模块），q<0.05 的格需人工复核计数（obs≥3 已过滤）。

## 7. 产物清单

- `a4_1_module_heatmap_mouse.pdf`
- `a4_1_module_heatmap_mouse.png`
- `a4_2_eigenprofiles_mouse.pdf`
- `a4_2_eigenprofiles_mouse.png`
- `a4_3_enrichment_mouse.pdf`
- `a4_3_enrichment_mouse.png`
- `a4_4_concordance_mouse.pdf`
- `a4_4_concordance_mouse.png`
- `a4_common_mechanism_mouse.md`
- `consensus_k_selection.tsv`
- `cross_group_concordance.tsv`
- `enrichment_class_x_module.tsv`
- `gene_zprofile_matrix.tsv`
- `module_eigenprofile.tsv`
- `module_summary.tsv`
- `modules.tsv`
