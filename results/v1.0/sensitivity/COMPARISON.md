# 敏感性变体比较 — gene_set v1.0
> 生成 2026-11-24 03:13；基准 = `results/v1.0/tables/<物种>/`（默认口径，正式产物）；变体 = `results/v1.0/sensitivity/<变体>_<物种>/`。显著性判据 q_bin<0.05（BH，族=species×group）。

## 1. 指标总表
| 变体 | 物种 | n_bins | n_feasible | n_single_donor_feasible | n_sig | n_sig_up | n_sig_down | pat_monotonic_up | pat_monotonic_down | pat_window_high | pat_endpoint_high | pat_transient_peak | pat_flat | n_specific | n_emergence | 保留率 | 方向一致率 | 新增显著 | 丢失显著 | 稳健性 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 默认口径(基准) | human | 13 | 88 | 29 | 185 | 146 | 39 | 135 | 37 | 30 | 73 | 33 | 761 | 9 | 32 |  |  |  |  |  |
| adult_16bin+exclude_disease | human | 16 | 93 | 41 | 189 | 148 | 41 | 130 | 34 | 45 | 48 | 39 | 765 | 11 | 34 | 0.838 | 1.0 | 34.0 | 30.0 | 部分敏感 |
| adult_16bin | human | 16 | 105 | 52 | 256 | 207 | 49 | 185 | 36 | 57 | 42 | 35 | 702 | 10 | 36 | 0.946 | 1.0 | 81.0 | 10.0 | 稳健 |
| exclude_disease | human | 13 | 87 | 29 | 166 | 126 | 40 | 116 | 37 | 32 | 68 | 33 | 783 | 10 | 31 | 0.849 | 1.0 | 9.0 | 28.0 | 部分敏感 |
| min_cells_100 | human | 13 | 71 | 22 | 112 | 77 | 35 | 66 | 33 | 16 | 64 | 28 | 862 | 9 | 19 | 0.53 | 1.0 | 14.0 | 87.0 | 高度口径依赖 |
| min_cells_50 | human | 13 | 83 | 27 | 184 | 145 | 39 | 134 | 37 | 25 | 71 | 29 | 773 | 9 | 27 | 0.995 | 1.0 | 0.0 | 1.0 | 稳健 |
| 默认口径(基准) | mouse | 12 | 92 | 13 | 136 | 85 | 51 | 79 | 51 | 35 | 93 | 10 | 727 | 6 | 14 |  |  |  |  |  |
| min_cells_100 | mouse | 12 | 74 | 7 | 107 | 62 | 45 | 54 | 45 | 23 | 76 | 6 | 700 | 8 | 15 | 0.713 | 1.0 | 10.0 | 39.0 | 部分敏感 |
| min_cells_50 | mouse | 12 | 82 | 10 | 124 | 79 | 45 | 71 | 45 | 17 | 97 | 9 | 759 | 6 | 15 | 0.801 | 1.0 | 15.0 | 27.0 | 部分敏感 |

## 2. 逐变体判读（相对基准）
- `adult_16bin+exclude_disease`(human) specific 名单变化: 失去 无; 新增 ['NPTX1', 'TWNK']; 主组改变 无。
- `adult_16bin+exclude_disease`(human): 显著行 185→189, 保留率 83.8%, 方向一致率 100.0%, 新增 34/丢失 30 → **部分敏感**; 锚点 ITPR1 1.0→1.0, MBP 0.929→0.95, ATOH1 -0.564→-0.581。
- `adult_16bin`(human) specific 名单变化: 失去 无; 新增 ['NPTX1']; 主组改变 无。
- `adult_16bin`(human): 显著行 185→256, 保留率 94.6%, 方向一致率 100.0%, 新增 81/丢失 10 → **稳健**; 锚点 ITPR1 1.0→1.0, MBP 0.929→0.964, ATOH1 -0.564→-0.601。
- `exclude_disease`(human) specific 名单变化: 失去 无; 新增 ['TWNK']; 主组改变 无。
- `exclude_disease`(human): 显著行 185→166, 保留率 84.9%, 方向一致率 100.0%, 新增 9/丢失 28 → **部分敏感**; 锚点 ITPR1 1.0→1.0, MBP 0.929→0.929, ATOH1 -0.564→-0.564。
- `min_cells_100`(human): 显著行 185→112, 保留率 53.0%, 方向一致率 100.0%, 新增 14/丢失 87 → **高度口径依赖**; 锚点 ITPR1 1.0→1.0, MBP 0.929→1.0, ATOH1 -0.564→-0.564。
- `min_cells_50`(human): 显著行 185→184, 保留率 99.5%, 方向一致率 100.0%, 新增 0/丢失 1 → **稳健**; 锚点 ITPR1 1.0→1.0, MBP 0.929→1.0, ATOH1 -0.564→-0.564。
- `min_cells_100`(mouse) specific 名单变化: 失去 ['PRKCG']; 新增 ['PTF1A', 'TTPA', 'ZFHX3']; 主组改变 无。
- `min_cells_100`(mouse): 显著行 136→107, 保留率 71.3%, 方向一致率 100.0%, 新增 10/丢失 39 → **部分敏感**; 锚点 ITPR1 0.915→0.976, MBP 0.943→0.9, ATOH1 -0.907→-0.907。
- `min_cells_50`(mouse) specific 名单变化: 失去 ['PRKCG']; 新增 ['ZFHX3']; 主组改变 无。
- `min_cells_50`(mouse): 显著行 136→124, 保留率 80.1%, 方向一致率 100.0%, 新增 15/丢失 27 → **部分敏感**; 锚点 ITPR1 0.915→0.976, MBP 0.943→0.943, ATOH1 -0.907→-0.907。

## 3. 稳健核心集（基准 + 全部变体都显著且方向一致 → 可写成断言）
- **human**：6 套口径（含基准）下都显著的 gene×group 行 = 86；其中方向也完全一致 = **86 行 / 46 个基因**（基准共 185 行显著 → 46% 可作跨口径断言）。
- **mouse**：3 套口径（含基准）下都显著的 gene×group 行 = 96；其中方向也完全一致 = **96 行 / 47 个基因**（基准共 136 行显著 → 71% 可作跨口径断言）。

| species | gene | n_groups | groups | direction | mech_class | tier |
|---|---|---|---|---|---|---|
| human | COQ8A | 4 | ALL,GranuleLineage,Progenitor,Purkinje | up | mito_metabolic | core |
| human | GRID2 | 4 | GranuleLineage,Interneuron,Progenitor,Purkinje | up | channel_synapse | core |
| human | SYNE1 | 4 | ALL,GranuleLineage,Interneuron,Progenitor | up | cytoskeleton_transport_glia | core |
| human | ATP2B3 | 3 | ALL,GranuleLineage,Purkinje | up | channel_synapse | core |
| human | ATXN1 | 3 | ALL,GranuleLineage,Progenitor | up | polyQ_repeat | core |
| human | CACNB4 | 3 | ALL,Interneuron,Progenitor | mixed | channel_synapse | core |
| human | FGF14 | 3 | GranuleLineage,Interneuron,Progenitor | up | channel_synapse | core |
| human | GRM1 | 3 | ALL,GranuleLineage,Interneuron | up | channel_synapse | core |
| human | KCNC3 | 3 | ALL,GranuleLineage,Interneuron | up | channel_synapse | core |
| human | KCND3 | 3 | ALL,GranuleLineage,Interneuron | up | channel_synapse | core |
| human | TMEM67 | 3 | ALL,GranuleLineage,Progenitor | mixed | cilia_congenital_dev | core |
| human | ZFHX3 | 3 | ALL,GranuleLineage,Progenitor | down | cilia_congenital_dev | core |
| human | AHI1 | 2 | Interneuron,Purkinje | up | cilia_congenital_dev | core |
| human | ATN1 | 2 | ALL,GranuleLineage | up | polyQ_repeat | core |
| human | CACNA1A | 2 | ALL,GranuleLineage | up | polyQ_repeat | core |
| human | CCDC88C | 2 | ALL,GranuleLineage | down | other | core |
| human | FAT2 | 2 | ALL,GranuleLineage | up | channel_synapse | core |
| human | KCNA1 | 2 | ALL,GranuleLineage | up | channel_synapse | core |
| human | NPTX1 | 2 | ALL,GranuleLineage | up | proteostasis_RNA | core |
| human | SACS | 2 | ALL,GranuleLineage | down | proteostasis_RNA | core |
| human | SLC1A3 | 2 | GranuleLineage,Progenitor | mixed | channel_synapse | core |
| human | SPTBN2 | 2 | GranuleLineage,Progenitor | mixed | cytoskeleton_transport_glia | core |
| human | STUB1 | 2 | ALL,GranuleLineage | up | proteostasis_RNA | core |
| human | TCTN2 | 2 | ALL,GranuleLineage | down | cilia_congenital_dev | core |
| human | TRPC3 | 2 | GranuleLineage,Interneuron | down | channel_synapse | core |
| human | AFG3L2 | 1 | Interneuron | up | mito_metabolic | core |
| human | CA8 | 1 | Interneuron | down | channel_synapse | core |
| human | CASK | 1 | ALL | down | channel_synapse | core |
| human | CYP27A1 | 1 | Progenitor | up | mito_metabolic | core |
| human | DAB1 | 1 | GranuleLineage | up | cilia_congenital_dev | core |
| human | ELOVL4 | 1 | GranuleLineage | up | mito_metabolic | core |
| human | ITPR1 | 1 | Purkinje | up | channel_synapse | core |
| human | NPHP1 | 1 | GranuleLineage | down | cilia_congenital_dev | core |
| human | PDSS2 | 1 | GranuleLineage | down | mito_metabolic | core |
| human | PHYH | 1 | Purkinje | up | mito_metabolic | core |
| human | PMM2 | 1 | Progenitor | down | mito_metabolic | core |
| human | PPP2R2B | 1 | Interneuron | down | proteostasis_RNA | core |
| human | PTF1A | 1 | ALL | down | cilia_congenital_dev | core |
| human | RFC1 | 1 | Interneuron | up | dna_repair | core |
| human | RPGRIP1L | 1 | ALL | down | cilia_congenital_dev | core |
| human | SCYL1 | 1 | Interneuron | up | cytoskeleton_transport_glia | core |
| human | SETX | 1 | Interneuron | up | dna_repair | core |
| human | SYT14 | 1 | Progenitor | down | cytoskeleton_transport_glia | core |
| human | TMEM240 | 1 | GranuleLineage | up | proteostasis_RNA | core |
| human | TPP1 | 1 | Progenitor | up | mito_metabolic | core |
| human | TTPA | 1 | GranuleLineage | down | mito_metabolic | core |
| mouse | ATP2B3 | 5 | ALL,DN,GranuleLineage,Interneuron,Purkinje | up | channel_synapse | core |
| mouse | FAT2 | 5 | ALL,GranuleLineage,Interneuron,Progenitor,Purkinje | up | channel_synapse | core |
| mouse | SYNE1 | 5 | ALL,DN,GranuleLineage,Progenitor,Purkinje | up | cytoskeleton_transport_glia | core |
| mouse | EEF2 | 4 | ALL,GranuleLineage,Interneuron,Progenitor | down | proteostasis_RNA | core |
| mouse | GRID2 | 4 | DN,GranuleLineage,Interneuron,Purkinje | up | channel_synapse | core |
| mouse | GRM1 | 4 | ALL,DN,GranuleLineage,Purkinje | up | channel_synapse | core |
| mouse | KCNA1 | 4 | ALL,GranuleLineage,Interneuron,Purkinje | up | channel_synapse | core |
| mouse | KCNC3 | 4 | ALL,GranuleLineage,Interneuron,Purkinje | up | channel_synapse | core |
| mouse | NPHP1 | 3 | ALL,DN,GranuleLineage | mixed | cilia_congenital_dev | core |
| mouse | PHYH | 3 | ALL,DN,Purkinje | up | mito_metabolic | core |
| mouse | PRKCG | 3 | ALL,GranuleLineage,Purkinje | up | channel_synapse | core |
| mouse | TRPC3 | 3 | DN,GranuleLineage,Purkinje | mixed | channel_synapse | core |
| mouse | ABCB7 | 2 | ALL,GranuleLineage | down | mito_metabolic | core |
| mouse | APTX | 2 | ALL,GranuleLineage | down | dna_repair | core |

完整清单见 `robust_core_genes.tsv`（93 行）。

## 4. 三个趋势锚点跨变体对照
| 变体 | 物种 | anchor_ITPR1 | anchor_MBP | anchor_ATOH1 |
|---|---|---|---|---|
| 默认口径(基准) | human | 1.0 | 0.929 | -0.564 |
| adult_16bin+exclude_disease | human | 1.0 | 0.95 | -0.581 |
| adult_16bin | human | 1.0 | 0.964 | -0.601 |
| exclude_disease | human | 1.0 | 0.929 | -0.564 |
| min_cells_100 | human | 1.0 | 1.0 | -0.564 |
| min_cells_50 | human | 1.0 | 1.0 | -0.564 |
| 默认口径(基准) | mouse | 0.915 | 0.943 | -0.907 |
| min_cells_100 | mouse | 0.976 | 0.9 | -0.907 |
| min_cells_50 | mouse | 0.976 | 0.943 | -0.907 |

## 5. 怎么读这份表
- **保留率** = 基准里显著的 gene×group 行，在该变体下仍显著的比例。≥90% 记「稳健」，70–90%「部分敏感」，<70%「高度口径依赖」。
- **方向一致率** = 保留下来的行里升降方向不变的比例。方向翻转比「变得不显著」更严重（前者是结论反了，后者只是功效不够）。
- **新增显著** 通常来自口径变化后 bin 数/细胞数改变，不一定是新生物学；须逐条看 `sig_flips.tsv`。
- 锚点（ITPR1/MBP/ATOH1）若在某个变体下方向翻转 → 该变体口径不可用（属结构性问题，不是功效问题）。
- 写报告时的用法：主分析用默认口径；只把在多口径下**方向一致且保留率高**的结论写成断言，其余降级为「口径依赖」或不写。
