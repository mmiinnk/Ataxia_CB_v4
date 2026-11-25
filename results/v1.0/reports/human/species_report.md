# human 物种报告 — gene_set v1.0

## 1. 头部
- 基因表版本: v1.0(config gene_set.version; 冻结状态见 genes/curation_log.md)
- 统计 panel: 107 基因, cube 命中 107, 缺失 0(status=missing)
- 供体: 27(疾病供体 3, donor cube disease_flag 口径)
- cube_bin: 130 组×bin 行, 可行 88, 不可行 42

## 2. 证据结构(可行性 / 单供体 / pooled)
| group | CS18 | CS19 | CS22 | 9w | 11w | 17w | 20w | NB | 6mo | 7mo | 9mo | 3y | AD |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Purkinje | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗0 | ✗0 | ✗0 | ✗4 | ✗1 |
| GranuleLineage | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| DN | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗0 | ✗0 | ✗0 | ✗0 | ✗0 | ✗0 | ✗0 |
| Interneuron | ✗2 | ✗6 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗17 | ✓ | ✓ | ✓ | ✓ |
| OligoLineage | ✗0 | ✗0 | ✗2 | ✗0 | ✓ | ✓ | ✓ | ✓ | ✗15 | ✓ | ✓ | ✓ | ✓ |
| Astro | ✗0 | ✗0 | ✗2 | ✗0 | ✗8 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Progenitor | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗5 | ✗3 | ✗3 | ✗2 | ✗0 |
| Immune | ✗6 | ✗9 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗3 | ✗11 | ✗8 | ✓ | ✓ |
| UBC | ✗0 | ✗0 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗0 | ✗0 | ✗0 | ✗0 | ✗0 |
| ALL | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

✓=可行(≥min_cells);✗n=不可行(n 核);·=骨架无此行。
- 单供体 bin(sem=NA, 硬约束 5): 68 行 — Purkinje@CS19, Purkinje@20w, Purkinje@6mo, Purkinje@7mo, Purkinje@9mo, Purkinje@3y, Purkinje@AD, GranuleLineage@CS19, GranuleLineage@20w, GranuleLineage@6mo, GranuleLineage@7mo, GranuleLineage@9mo, DN@CS19, DN@20w, DN@NB, DN@6mo, DN@7mo, DN@9mo, DN@3y, DN@AD, Interneuron@CS19, Interneuron@20w, Interneuron@6mo, Interneuron@7mo, Interneuron@9mo, OligoLineage@CS18, OligoLineage@CS19, OligoLineage@CS22, OligoLineage@9w, OligoLineage@20w, OligoLineage@6mo, OligoLineage@7mo, OligoLineage@9mo, Astro@CS18, Astro@CS19, Astro@CS22, Astro@9w, Astro@20w, Astro@6mo, Astro@7mo, Astro@9mo, Astro@3y, Progenitor@CS19, Progenitor@20w, Progenitor@6mo, Progenitor@7mo, Progenitor@9mo, Progenitor@3y, Progenitor@AD, Immune@CS19, Immune@20w, Immune@6mo, Immune@7mo, Immune@9mo, UBC@CS18, UBC@CS19, UBC@CS22, UBC@20w, UBC@6mo, UBC@7mo, UBC@9mo, UBC@3y, UBC@AD, ALL@CS19, ALL@20w, ALL@6mo, ALL@7mo, ALL@9mo
- pooled bin(sem=NA): 0 行 — 无

## 3. 第 1 层: 趋势(Spearman rho_bin, BH q<0.05, 族=species×group)
| group | 基因行 | q<0.05 | 升 | 降 | rho NA |
|---|---|---|---|---|---|
| Purkinje | 107 | 47 | 45 | 2 | 3 |
| GranuleLineage | 107 | 33 | 19 | 14 | 0 |
| DN | 107 | 35 | 34 | 1 | 3 |
| Interneuron | 107 | 28 | 23 | 5 | 0 |
| OligoLineage | 107 | 0 | 0 | 0 | 1 |
| Astro | 107 | 0 | 0 | 0 | 1 |
| Progenitor | 107 | 14 | 9 | 5 | 2 |
| Immune | 107 | 1 | 1 | 0 | 2 |
| UBC | 107 | 0 | 0 | 0 | 3 |
| ALL | 107 | 27 | 15 | 12 | 0 |


## 4. 第 2 层: 模式 / 特异性 / emergence / 峰期
- temporal_pattern 分布(gene×group 行): flat=761, monotonic_up=135, endpoint_high=73, monotonic_down=37, transient_peak=33, window_high=30, NA=1
- specific_call 分布: broad=94, specific=9, not_expressed=4
- specific 基因(τ 与主组): ITPR1→Purkinje(τ=0.976, frac=0.856), TRPC3→Purkinje(τ=0.972, frac=0.836), CA8→Purkinje(τ=0.957, frac=0.770), SAMD9L→Immune(τ=0.938, frac=0.698), PTF1A→Progenitor(τ=0.933, frac=0.682), FAT2→GranuleLineage(τ=0.916, frac=0.629), RELN→GranuleLineage(τ=0.911, frac=0.615), TTPA→Astro(τ=0.907, frac=0.606), DAB1→Purkinje(τ=0.894, frac=0.574)
- not_expressed 基因(4): KCNA1, PDYN, TMEM240, VAMP1
- 低表达伪标签掩膜(P7 裁决①): 1 行 window/peak→NA(TMEM240)
- endpoint_high 同类残余(不在裁决①范围, 留档): 15 行(KCNA1, PDYN, TMEM240, VAMP1)
- emergence 可判: 32/107(左截断 2)
- 峰期×组分布(基因数, peak_by_group):
| group | embryonic | fetal | infant | adult |
|---|---|---|---|---|
| Purkinje | 23 | 23 | 61 | 0 |
| GranuleLineage | 45 | 12 | 19 | 31 |
| DN | 32 | 75 | 0 | 0 |
| Interneuron | 13 | 11 | 39 | 44 |
| OligoLineage | 0 | 21 | 59 | 27 |
| Astro | 0 | 27 | 70 | 10 |
| Progenitor | 57 | 30 | 20 | 0 |
| Immune | 22 | 40 | 24 | 21 |
| UBC | 36 | 40 | 31 | 0 |
| ALL | 33 | 18 | 25 | 31 |


## 5. 锚点与门自动核查(§4 + §9)
- 全部 logs: FAIL 0 条, WARN 7 条(逐条见 reviews/gate2_human.md)
| 来源 | check | observed | level |
|---|---|---|---|
| s00_report | human: X 形状锚点 | (180956, 27410) → (180956, 27410) | PASS |
| s00_report | human: size_factor↔X行和锚点 | r≥0.999, 中位比值∈[0.95,1.02] → r=1.00000, 比值=0.9974 | PASS |
| s00_report | human: 9 组细胞数锚点 §4.1 | 全部精确匹配 → 匹配 | PASS |
| s00_report | mouse: X 形状锚点 | (115282, 20409) → (115282, 20409) | PASS |
| s00_report | mouse: size_factor↔X行和锚点 | r≥0.999, 中位比值∈[0.95,1.02] → r=0.99997, 比值=0.9852 | PASS |
| s00_report | mouse: 9 组细胞数锚点 §4.1 | 全部精确匹配 → 匹配 | PASS |
| s00_report | aldinger: 形状锚点 §4.3 | (69174, 26326) → (69174, 26326) | PASS |
| s00_report | aldinger: H-PC 细胞数锚点 §4.3 | 25724 → 25724 | PASS |
| s00_report | aldinger: 供体数锚点 §4.3 | 13 → 13 | PASS |
| s00_report | castonguay: 形状锚点 §4.3 | (1004112, 47789) → (1004112, 47789) | PASS |
| s00_report | castonguay: normal 供体/PC 锚点 §4.3 | {'n_donors_normal': 93, 'n_pc_normal': 31089, 'n_donors_ge30pc': 89, 'median_pc_ge30': 257} → {'n_donors_normal': 93, 'n | WARN |
| s01_report_human | 形状锚点 | (180956, 27410) → (180956, 27410) | PASS |
| s01_report_human | 9 组细胞数锚点 §4.1 | 精确匹配 → 匹配 | PASS |
| s01_report_human | 趋势锚点 ITPR1×Purkinje | rho + → rho=+1.000(n=8 bin, 粗算) | PASS |
| s01_report_human | 趋势锚点 MBP×OligoLineage | rho + → rho=+0.929(n=8 bin, 粗算) | PASS |
| s01_report_human | 趋势锚点 ATOH1×GranuleLineage | rho - → rho=-0.564(n=13 bin, 粗算) | PASS |
| s02_report_human | 趋势锚点 ITPR1×Purkinje(正式) | rho + | rho_bin=+1.000 q_bin=0 status=ok | PASS |
| s02_report_human | 趋势锚点 MBP×OligoLineage(正式) | rho + | rho_bin=+0.929 q_bin=0.000863 status=ok | PASS |
| s02_report_human | 趋势锚点 ATOH1×GranuleLineage(正式) | rho - | rho_bin=-0.564 q_bin=0.0448 status=ok | PASS |

WARN 清单:
| 来源 | check | observed | note |
|---|---|---|---|
| s00_report | castonguay: normal 供体/PC 锚点 §4.3 | {'n_donors_normal': 93, 'n_pc_normal': 31089, 'n_donors_ge30pc': 89, 'median_pc_ge30': 257} → {'n_donors_normal': 93, 'n | median_pc_ge30: 实测 264.0 != 锚点 257 |
| s00_report | mouse: 鼠 dropout 现场重算 §3.12 | 锚点±容差内(WARN=0) → WARN=2, ABSENT=0 | 详见 logs/s00_dropout.tsv |
| s00_report | external: external 文件 brainspan | /home/victor/Projects/Ataxia_CB_v4/external/brainspan_cerebellum_rpkm.tsv → 缺席 | s10 对缺失文件只 WARN 并出占位 panel(§11 未决 6) |
| s00_report | external: external 文件 gtex | /home/victor/Projects/Ataxia_CB_v4/external/gtex_v8_cerebellum_tpm.tsv → 缺席 | s10 对缺失文件只 WARN 并出占位 panel(§11 未决 6) |
| s00_report | external: external 文件 model_mouse_deg | /home/victor/Projects/Ataxia_CB_v4/external/model_mouse_deg_signatures.tsv → 缺席 | s10 对缺失文件只 WARN 并出占位 panel(§11 未决 6) |
| s01_report_human | 组×bin 可行性矩阵 §4.2 | - → 人侧主口径已合并成人 bin(collapse_adult_bins=true; 2026-11-24 用户裁决, 门 2 整改项): spec §4.2 可行性锚点表按 16 bin 记录, 与合并后 13 bin 不可比 → 跳 |  |
| s01_report_human | marker 校验(spec §9 清单) | 期望组=argmax(或≥50%argmax) → FAIL=0(清单内) / 3(清单外), ABSENT=4 |  |


## 6. 本物种相关 §11 未决条目
- 条目编号: 2, 7, 9(条文以 implementation_spec §11 为准; 报告只列编号防漂移)

## 7. 图清单(figures/human)
| 视图 | 文件茎数 | 茎(前 6) |
|---|---|---|
| a1 | 2 | human_a1_ALL_p1, human_a1_ALL_p2 |
| b1 | 18 | human_b1_Astro_p1, human_b1_Astro_p2, human_b1_DN_p1, human_b1_DN_p2, human_b1_GranuleLineage_p1, human_b1_GranuleLineage_p2… |
| b2 | 81 | human_b2_Astro_p1, human_b2_Astro_p2, human_b2_Astro_p3, human_b2_Astro_p4, human_b2_Astro_p5, human_b2_Astro_p6… |
| b3 | 2 | human_b3_p1, human_b3_p2 |
| b4 | 1 | human_b4 |
| c1 | 1 | human_c1 |
| c2 | 2 | human_c2_p1, human_c2_p2 |
| d1 | 1 | human_d1 |


## 8. review 焦点(spec §9: b2 轨迹与 c1 dotplot)
- /home/victor/Projects/Ataxia_CB_v4/results/v1.0/figures/human/human_b2_Astro_p1.pdf
- /home/victor/Projects/Ataxia_CB_v4/results/v1.0/figures/human/human_b2_DN_p1.pdf
- /home/victor/Projects/Ataxia_CB_v4/results/v1.0/figures/human/human_b2_GranuleLineage_p1.pdf
- /home/victor/Projects/Ataxia_CB_v4/results/v1.0/figures/human/human_b2_Immune_p1.pdf
- /home/victor/Projects/Ataxia_CB_v4/results/v1.0/figures/human/human_b2_Interneuron_p1.pdf
- /home/victor/Projects/Ataxia_CB_v4/results/v1.0/figures/human/human_b2_OligoLineage_p1.pdf
- /home/victor/Projects/Ataxia_CB_v4/results/v1.0/figures/human/human_b2_Progenitor_p1.pdf
- /home/victor/Projects/Ataxia_CB_v4/results/v1.0/figures/human/human_b2_Purkinje_p1.pdf
- /home/victor/Projects/Ataxia_CB_v4/results/v1.0/figures/human/human_b2_UBC_p1.pdf
- /home/victor/Projects/Ataxia_CB_v4/results/v1.0/figures/human/human_c1.pdf

## 9. 口径注
- NA 政策(硬约束 5): 组×bin < min_cells → 表 NA、图斜纹; 不插值、不填 0; 单供体/pooled bin 的 sem=NA。
- b2 轨迹 y 轴 = mean log1p CP10K(与存储的 logsem_/rho_bin 自洽; 读图指南 'y=CP10K' 的落地注记)。
- c1 dotplot 点大小 = pct(0–1 表达核比例)× dot_size_max; 色 = 可行 bin log 均值。
- gene_set v1.0 已冻结 → 图形**不带** DRAFT footer, 数值可引用; 引用时须附口径注(min_cells 阈值、可行 组×bin 数、单供体/pooled bin 的 sem=NA、右删失的组与 bin 段)。
- 低表达伪标签(P7 裁决①, 2026-11-05): not_expressed 基因的 window_high/transient_peak 已在 s02 掩膜为 NA(图灰色/斜纹, 下游自动排除); endpoint_high 同类残余仅计数留档。
