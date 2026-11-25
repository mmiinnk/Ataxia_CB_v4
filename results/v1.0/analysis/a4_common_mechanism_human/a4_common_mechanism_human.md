# A4 共性机制筛选报告 — human（gene_set v1.0）

> 一次性分析层（`scripts/a4_common_mechanism.py`，非 sNN 契约），生成 2026-11-24 03:14，seed=20261013，n_reps=200，n_perm=10000，单供体 bin 保留(与 s02 同口径)。

## 1. 口径与输入

- 输入：`results/v1.0/intermediate/human/cube_bin.tsv`（可行 组×bin）、`tables/human/specificity_global.tsv`、`trend_by_group.tsv`、`temporal_patterns.tsv`、`genes/ataxia_genes.tsv`。
- 空间维 = 8 个主细胞类型组（Purkinje, GranuleLineage, DN, Interneuron, OligoLineage, Astro, Progenitor, Immune）；时间维 = 发育 bin；UBC（descriptive）与 ALL（reference）不入聚类。
- 聚类矩阵：103 基因 × 69 个可行 组×bin 列（单供体 22、pooled 0）。
- **NA 政策（硬约束 5）**：不可行 组×bin 为 NA，不插值、不填 0；相关用成对完整观测（min_periods=20），相关缺失对 0（按距离 1 中性处理）。
- 排除：`not_expressed` 4（KCNA1, PDYN, TMEM240, VAMP1）、零方差 0（无）、观测不足 0（无）；cube 未命中 0（无）。
- **CL-027**：机制类富集排除残差兜底类 ['other', 'unknown'] → 剔除 3 基因（CCDC88C, MME, WWOX），这些基因改做单基因描述，不产生「残差类富集」结论。

## 2. 共识聚类与 k 选择

| k | PAC | min_module_size | mean_within_consensus |
|---|---|---|---|
| 3 | 0.7514 | 1 | 0.4282 |
| 4 | 0.5726 | 1 | 0.5229 |
| 5 | 0.569 | 1 | 0.6422 |
| 6 | 0.4428 | 1 | 0.6715 |
| 7 | 0.4365 | 1 | 0.7215 |
| 8 | 0.4298 | 1 | 0.7388 |

选定 **k=8**（PAC 最小）；模块内平均共识见下表 `mean_within_consensus`（越接近 1 越稳定）。

## 3. 模块刻画（共享时空程序）

| module | n_genes | direction | peak_group | peak_bin | peak_period | top_groups | hubs | mech_composition | mean_within_consensus |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 10 | flat_or_peak | Astro | 7mo | infant | Astro(+1.47),OligoLineage(+0.34),Progenitor(+0.18) | TTPA,NPHP1,CYP27A1,ATXN1,WWOX | mito_metabolic:2,polyQ_repeat:2,channel_synapse:2,cilia_congenital_dev:1,other:1,proteostasis_RNA:1,cytoskeleton_transport_glia:1 | 0.8499 |
| 2 | 3 | flat_or_peak | Progenitor | NB | infant | Progenitor(+0.61),Astro(+0.31),Interneuron(+0.10) | TWNK,COQ6,MTTP | mito_metabolic:3 | 0.7274 |
| 3 | 21 | rising | Immune | AD | adult | Immune(+1.10),DN(+0.17),Purkinje(+0.01) | RNF216,SETX,VPS13D,PNPLA6,ATM | mito_metabolic:7,dna_repair:4,cytoskeleton_transport_glia:4,polyQ_repeat:3,proteostasis_RNA:2,cilia_congenital_dev:1 | 0.6335 |
| 4 | 17 | rising | DN | 17w | fetal | DN(+0.75),Purkinje(+0.61),Progenitor(+0.24) | RPGRIP1L,CAMTA1,SACS,CEP290,CC2D2A | cilia_congenital_dev:9,proteostasis_RNA:3,channel_synapse:2,dna_repair:1,cytoskeleton_transport_glia:1,mito_metabolic:1 | 0.743 |
| 5 | 16 | flat_or_peak | GranuleLineage | CS19 | embryonic | Progenitor(+0.23),OligoLineage(+0.03),GranuleLineage(+0.01) | RFC1,UBA5,ABCB7,OFD1,APTX | dna_repair:5,mito_metabolic:5,proteostasis_RNA:3,cilia_congenital_dev:2,polyQ_repeat:1 | 0.58 |
| 6 | 3 | rising | Astro | 6mo | infant | GranuleLineage(+0.36),Purkinje(+0.23),Astro(+0.20) | MME,BEAN1,PHYH | other:1,proteostasis_RNA:1,mito_metabolic:1 | 0.8791 |
| 7 | 32 | rising | Purkinje | NB | infant | Purkinje(+0.54),Interneuron(+0.46),GranuleLineage(+0.15) | ATP2B3,ELOVL4,GRM1,CACNA1A,KCNC3 | channel_synapse:12,cilia_congenital_dev:6,mito_metabolic:5,cytoskeleton_transport_glia:5,proteostasis_RNA:3,polyQ_repeat:1 | 0.7585 |
| 8 | 1 | flat_or_peak | GranuleLineage | 9w | fetal | Progenitor(+0.65),GranuleLineage(+0.46),Purkinje(+0.11) | CCDC88C | other:1 | 1.0 |

`hubs` = 模块内平均 |corr| 最高的前 5 个基因（共性程序的代表节点，后续网络层 n00–n04 的候选入口）。

## 4. 机制类富集（mech_class × module / × shape）

显著格（q<0.05 且 obs≥3）共 1 个，按 log2FE 降序前 20：

| test | a | b | obs | exp | log2fe | p | q |
|---|---|---|---|---|---|---|---|
| mech_class_x_shape | channel_synapse | flat | 3 | 10.08 | -1.596 | 0.0005 | 0.0175 |

## 5. 跨组趋势一致性（共性 vs 异质性）

| concordance | count |
|---|---|
| partial_concordant | 37 |
| ns | 32 |
| single_group | 25 |
| mixed_direction | 11 |
| pan_up | 2 |

- **pan_up**（2）：FGF14, GRID2

- **mixed_direction**（11）：SLC1A3, TRPC3, CA8, CACNB4, CCDC88C, NPHP1, TCTN2, TMEM67, PDSS2, PPP2R2B, SPTBN2

判读：`pan_*` = ≥5 个主组同向显著（跨细胞类型的**共性**发育程序）；`mixed_direction` = 同时存在升与降（**异质性**的细胞类型分歧）；`single_group` = 仅一组显著（细胞类型特异）。

## 6. 局限（写报告时必须一并说明）

- 单物种（鼠）、正常发育数据，**无疾病对照** → 只能说「致病基因在正常小脑发育中的共享时空程序」，不能说「共性致病机制」已被证明；后者需要 s04/s05 跨物种 + 疾病投影层（b00/b01、n00–n04）。
- 模块为无监督发现，属**假设生成**；k 由 PAC 选定，需做 k±1 与 min_cells 敏感性后才可写入正文。
- 单供体/pooled bin 参与聚类（与 s02 同口径），早期胚胎段的模块稳定性弱于晚期。
- 富集样本量小（入富集 100 基因 × k=8 模块），q<0.05 的格需人工复核计数（obs≥3 已过滤）。

## 7. 产物清单

- `a4_1_module_heatmap_human.pdf`
- `a4_1_module_heatmap_human.png`
- `a4_2_eigenprofiles_human.pdf`
- `a4_2_eigenprofiles_human.png`
- `a4_3_enrichment_human.pdf`
- `a4_3_enrichment_human.png`
- `a4_4_concordance_human.pdf`
- `a4_4_concordance_human.png`
- `consensus_k_selection.tsv`
- `cross_group_concordance.tsv`
- `enrichment_class_x_module.tsv`
- `gene_zprofile_matrix.tsv`
- `module_eigenprofile.tsv`
- `module_summary.tsv`
- `modules.tsv`
