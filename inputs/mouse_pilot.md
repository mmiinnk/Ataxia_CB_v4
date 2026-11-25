# mouse 物种报告 — gene_set v0.9.1-draft(DRAFT 探索性口径)

## 1. 头部
- 基因表版本: v0.9.1-draft(config gene_set.version; 冻结状态见 genes/curation_log.md)
- 统计 panel: 109 基因, cube 命中 103, 缺失 6(status=missing)
- 供体: 24
- cube_bin: 120 组×bin 行, 可行 92, 不可行 28

## 2. 证据结构(可行性 / 单供体 / pooled)
| group | E10.5 | E11.5 | E12.5 | E13.5 | E14.5 | E15.5 | E17.5 | P0 | P4 | P7 | P14 | 9wk |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Purkinje | ✗21 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗21 |
| GranuleLineage | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| DN | ✗16 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗9 | ✗7 |
| Interneuron | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| OligoLineage | ✗6 | ✗4 | ✗0 | ✗1 | ✗6 | ✗6 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Astro | ✗0 | ✗0 | ✗0 | ✗13 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Progenitor | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Immune | ✓ | ✓ | ✗14 | ✗27 | ✗6 | ✓ | ✗16 | ✓ | ✓ | ✓ | ✗21 | ✗19 |
| UBC | ✗0 | ✗0 | ✗0 | ✗0 | ✗27 | ✓ | ✓ | ✓ | ✓ | ✓ | ✗17 | ✗14 |
| ALL | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

✓=可行(≥min_cells);✗n=不可行(n 核);·=骨架无此行。
- 单供体 bin(sem=NA, 硬约束 5): 30 行 — Purkinje@E10.5, Purkinje@E11.5, GranuleLineage@E10.5, GranuleLineage@E11.5, DN@E10.5, DN@E11.5, Interneuron@E10.5, Interneuron@E11.5, OligoLineage@E10.5, OligoLineage@E11.5, OligoLineage@E12.5, OligoLineage@E13.5, Astro@E10.5, Astro@E11.5, Astro@E12.5, Astro@E13.5, Astro@E14.5, Progenitor@E10.5, Progenitor@E11.5, Immune@E10.5, Immune@E11.5, UBC@E10.5, UBC@E11.5, UBC@E12.5, UBC@E13.5, UBC@E14.5, UBC@P14, UBC@9wk, ALL@E10.5, ALL@E11.5
- pooled bin(sem=NA): 16 行 — Purkinje@E10.5, Purkinje@E11.5, GranuleLineage@E10.5, GranuleLineage@E11.5, DN@E10.5, DN@E11.5, Interneuron@E10.5, Interneuron@E11.5, OligoLineage@E10.5, OligoLineage@E11.5, Progenitor@E10.5, Progenitor@E11.5, Immune@E10.5, Immune@E11.5, ALL@E10.5, ALL@E11.5

## 3. 第 1 层: 趋势(Spearman rho_bin, BH q<0.05, 族=species×group)
| group | 基因行 | q<0.05 | 升 | 降 | rho NA |
|---|---|---|---|---|---|
| Purkinje | 109 | 28 | 27 | 1 | 7 |
| GranuleLineage | 109 | 39 | 12 | 27 | 7 |
| DN | 109 | 12 | 11 | 1 | 7 |
| Interneuron | 109 | 17 | 15 | 2 | 6 |
| OligoLineage | 109 | 0 | 0 | 0 | 7 |
| Astro | 109 | 3 | 2 | 1 | 7 |
| Progenitor | 109 | 17 | 9 | 8 | 7 |
| Immune | 109 | 0 | 0 | 0 | 13 |
| UBC | 109 | 0 | 0 | 0 | 9 |
| ALL | 109 | 24 | 10 | 14 | 6 |


## 4. 第 2 层: 模式 / 特异性 / emergence / 峰期
- temporal_pattern 分布(gene×group 行): flat=742, endpoint_high=93, monotonic_up=80, NA=75, monotonic_down=54, window_high=36, transient_peak=10
- specific_call 分布: broad=89, not_expressed=8, NA=6, specific=6
- specific 基因(τ 与主组): ITPR1→Purkinje(τ=0.954, frac=0.756), DAB1→Purkinje(τ=0.899, frac=0.587), CYP27A1→OligoLineage(τ=0.899, frac=0.586), SLC1A3→Astro(τ=0.884, frac=0.553), RELN→GranuleLineage(τ=0.879, frac=0.542), PRKCG→Purkinje(τ=0.865, frac=0.515)
- not_expressed 基因(8): BEAN1, COQ8A, MME, SPG7, TTBK2, TTPA, TWNK, WDR81
- 低表达伪标签掩膜(P7 裁决①): 15 行 window/peak→NA(BEAN1, MME, SPG7, TTBK2, TTPA, TWNK, WDR81)
- endpoint_high 同类残余(不在裁决①范围, 留档): 18 行(BEAN1, COQ8A, MME, SPG7, TTPA, WDR81)
- emergence 可判: 14/109(左截断 2)
- 峰期×组分布(基因数, peak_by_group):
| group | embryonic | perinatal | postnatal | adult |
|---|---|---|---|---|
| Purkinje | 21 | 25 | 57 | 0 |
| GranuleLineage | 51 | 13 | 6 | 33 |
| DN | 17 | 40 | 46 | 0 |
| Interneuron | 34 | 6 | 6 | 57 |
| OligoLineage | 0 | 45 | 27 | 31 |
| Astro | 22 | 28 | 21 | 32 |
| Progenitor | 23 | 15 | 35 | 30 |
| Immune | 80 | 13 | 10 | 0 |
| UBC | 17 | 57 | 29 | 0 |
| ALL | 41 | 29 | 4 | 29 |


## 5. 锚点与门自动核查(§4 + §9)
- 全部 logs: FAIL 0 条, WARN 6 条(逐条见 reviews/gate1_mouse.md)
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
| s01_report_mouse | 形状锚点 | (115282, 20409) → (115282, 20409) | PASS |
| s01_report_mouse | 9 组细胞数锚点 §4.1 | 精确匹配 → 匹配 | PASS |
| s01_report_mouse | 趋势锚点 ITPR1×Purkinje | rho + → rho=+0.915(n=10 bin, 粗算) | PASS |
| s01_report_mouse | 趋势锚点 MBP×OligoLineage | rho + → rho=+0.943(n=6 bin, 粗算) | PASS |
| s01_report_mouse | 趋势锚点 ATOH1×GranuleLineage | rho - → rho=-0.907(n=12 bin, 粗算) | PASS |
| s02_report_mouse | 趋势锚点 ITPR1×Purkinje(正式) | rho + | rho_bin=+0.915 q_bin=0.000204 status=ok | PASS |
| s02_report_mouse | 趋势锚点 MBP×OligoLineage(正式) | rho + | rho_bin=+0.943 q_bin=0.0167 status=ok | PASS |
| s02_report_mouse | 趋势锚点 ATOH1×GranuleLineage(正式) | rho - | rho_bin=-0.907 q_bin=4.63e-05 status=ok | PASS |

WARN 清单:
| 来源 | check | observed | note |
|---|---|---|---|
| s00_report | castonguay: normal 供体/PC 锚点 §4.3 | {'n_donors_normal': 93, 'n_pc_normal': 31089, 'n_donors_ge30pc': 89, 'median_pc_ge30': 257} → {'n_donors_normal': 93, 'n | median_pc_ge30: 实测 264.0 != 锚点 257 |
| s00_report | mouse: 鼠 dropout 现场重算 §3.12 | 锚点±容差内(WARN=0) → WARN=2, ABSENT=0 | 详见 logs/s00_dropout.tsv |
| s00_report | external: external 文件 brainspan | /home/victor/Projects/Ataxia_CB/external/brainspan_cerebellum_rpkm.tsv → 缺席 | s10 对缺失文件只 WARN 并出占位 panel(§11 未决 6) |
| s00_report | external: external 文件 gtex | /home/victor/Projects/Ataxia_CB/external/gtex_v8_cerebellum_tpm.tsv → 缺席 | s10 对缺失文件只 WARN 并出占位 panel(§11 未决 6) |
| s00_report | external: external 文件 model_mouse_deg | /home/victor/Projects/Ataxia_CB/external/model_mouse_deg_signatures.tsv → 缺席 | s10 对缺失文件只 WARN 并出占位 panel(§11 未决 6) |
| s01_report_mouse | marker 校验(spec §9 清单) | 期望组=argmax(或≥50%argmax) → FAIL=0(清单内) / 1(清单外), ABSENT=2 |  |


## 6. 本物种相关 §11 未决条目
- 条目编号: 3, 8(条文以 implementation_spec §11 为准; 报告只列编号防漂移)
- 注: §11.1 鼠免疫 marker 按 dropout 豁免(Cx3cr1/P2ry12/Hexb 口径), 见门清单条目(3)。

## 7. 图清单(figures/mouse)
| 视图 | 文件茎数 | 茎(前 6) |
|---|---|---|
| a1 | 2 | mouse_a1_ALL_p1, mouse_a1_ALL_p2 |
| b1 | 18 | mouse_b1_Astro_p1, mouse_b1_Astro_p2, mouse_b1_DN_p1, mouse_b1_DN_p2, mouse_b1_GranuleLineage_p1, mouse_b1_GranuleLineage_p2… |
| b2 | 81 | mouse_b2_Astro_p1, mouse_b2_Astro_p2, mouse_b2_Astro_p3, mouse_b2_Astro_p4, mouse_b2_Astro_p5, mouse_b2_Astro_p6… |
| b3 | 2 | mouse_b3_p1, mouse_b3_p2 |
| b4 | 1 | mouse_b4 |
| c1 | 1 | mouse_c1 |
| c2 | 2 | mouse_c2_p1, mouse_c2_p2 |
| d1 | 1 | mouse_d1 |


## 8. review 焦点(spec §9: b2 轨迹与 c1 dotplot)
- /home/victor/Projects/Ataxia_CB/results/v0.9.1-draft/figures/mouse/mouse_b2_Astro_p1.pdf
- /home/victor/Projects/Ataxia_CB/results/v0.9.1-draft/figures/mouse/mouse_b2_DN_p1.pdf
- /home/victor/Projects/Ataxia_CB/results/v0.9.1-draft/figures/mouse/mouse_b2_GranuleLineage_p1.pdf
- /home/victor/Projects/Ataxia_CB/results/v0.9.1-draft/figures/mouse/mouse_b2_Immune_p1.pdf
- /home/victor/Projects/Ataxia_CB/results/v0.9.1-draft/figures/mouse/mouse_b2_Interneuron_p1.pdf
- /home/victor/Projects/Ataxia_CB/results/v0.9.1-draft/figures/mouse/mouse_b2_OligoLineage_p1.pdf
- /home/victor/Projects/Ataxia_CB/results/v0.9.1-draft/figures/mouse/mouse_b2_Progenitor_p1.pdf
- /home/victor/Projects/Ataxia_CB/results/v0.9.1-draft/figures/mouse/mouse_b2_Purkinje_p1.pdf
- /home/victor/Projects/Ataxia_CB/results/v0.9.1-draft/figures/mouse/mouse_b2_UBC_p1.pdf
- /home/victor/Projects/Ataxia_CB/results/v0.9.1-draft/figures/mouse/mouse_c1.pdf

## 9. 口径注
- NA 政策(硬约束 5): 组×bin < min_cells → 表 NA、图斜纹; 不插值、不填 0; 单供体/pooled bin 的 sem=NA。
- b2 轨迹 y 轴 = mean log1p CP10K(与存储的 logsem_/rho_bin 自洽; 读图指南 'y=CP10K' 的落地注记)。
- c1 dotplot 点大小 = pct(0–1 表达核比例)× dot_size_max; 色 = 可行 bin log 均值。
- 全部图形带 DRAFT footer(v0.9.1-draft 期间): 数值不得引用为论文数字。
- 低表达伪标签(P7 裁决①, 2026-11-05): not_expressed 基因的 window_high/transient_peak 已在 s02 掩膜为 NA(图灰色/斜纹, 下游自动排除); endpoint_high 同类残余仅计数留档。
