# 阶段报告 ↔ 图片/表格 对应清单

> 配套文件：`docs/progress_report_2025H2_atlas.md`（2025.9–2026.2 阶段进展报告，正文不含图）。
> 本清单把报告每一节的文字结论映射到**已存在**的图与表，供后续报告（含第一标注论文）直接取用。
> 全部路径相对项目根；图均为 PDF + PNG 双格式。gene_set v1.0；人侧为 13 bin 主口径（`collapse_adult_bins=true`）。
> 生成依据：`results/v1.0/logs/{s01,s02,s03,s04,s07,s09}_report_*.tsv`（三道人工门均已签署，见 `reviews/`）。

## 0. 报告正文配图（3 张，A4 版式；`scripts/a7_report_figures.py`，2026-11-24）

位置：`results/v1.0/figures/report/report_fig{1,2,3}.{pdf,png,eps}`（PDF/EPS 矢量、PNG 300 dpi；EPS 字体以 Type42 嵌入）。
字体：首选微软雅黑；本机无该字体时回退 Noto Sans CJK SC（不强制；将 msyh.ttc 放入 `fonts/` 后重跑即得雅黑版）。
版面：按 A4 印刷实际尺寸 1:1 设计（页边距 25 mm → 正文宽 160 mm）；图 1 = 160×105 mm、图 2 = 160×118 mm、
图 3 = 160×122 mm，各占约半页；字号 4.2–9 pt（1:1 下即所得）。
**图内只含面板字母、简短标题、轴名/刻度、图例与数据标签；不含口径说明文字与假设检验结果**——
图题、图注、统计量（p/q/δ）、判读提示统一见 `docs/执行情况报告_图注.md`。
配色（dataviz 六项检验实测；**不考虑黑白印刷**，但仍须色盲安全）：机制类别 7 色改用 Okabe-Ito 色系
（项目原 `mech_colors` 未过 CVD 分离：#ff7f0e↔#2ca02c 红色盲 ΔE=0.7）；轨迹三色 #0072B2/#D55E00/#009E73、
人/鼠双系列 #1f77b4/#d62728 均六项全过；有序量用单色相顺序色（Blues / viridis），缺失用斜纹（`na_style`）。

| 图 | 面板 | 内容 | 数据来源（无需重跑管线） |
| --- | --- | --- | --- |
| 图1 图谱构成与证据分布 | a | 小鼠 10 类群 × 12 分期细胞数矩阵（斜纹＝低于 30 核门槛，不作估计） | `intermediate/mouse/cube_bin.tsv` |
| | b | 人源同结构矩阵（13 分期） | `intermediate/human/cube_bin.tsv` |
| | c | 临床异质性四联横条：机制类/遗传方式/发病年龄/病程特征 | `genes/ataxia_genes.tsv`、`genes/ataxia_gene_disease.tsv` |
| 图2 共性在时序而非细胞类型 | a | 机制类表达高峰分布：行内彩色点＝各基因高峰分期、同色细线＝范围、黑粗线＝中位（右侧标中位分期）；按机制类着色（Okabe-Ito 7 色） | `tables/cross/peak_rank_by_class.tsv`（mouse_only 轴） |
| | b | 每个基因"相对表达最高类群"的归属计数（人/鼠并列）；**标题给结论**"致病基因不集中于单一细胞类群"，标题下一行一句大白话给计数规则，横轴仅"基因数" | `tables/{human,mouse}/specificity_global.tsv` |
| | c | 随发育呈现明确变化模式（非表达平稳）的基因占比（人/鼠并列）；**标题给结论**"离子通道与突触类基因几乎都随发育变化"，占比口径放纵轴 | `tables/cross/module_assignment.tsv`（`shape`）、`enrichment_class_x_module.tsv`（test=class_x_shape，q 值入图注） |
| 图3 代表基因与疾病级一致性 | a–d | GRM1 与 ATXN3 的人/鼠轨迹（三色实线：全小脑/浦肯野/颗粒谱系） | `intermediate/{sp}/cube_bin.tsv` 的 `log_GRM1`、`log_ATXN3` |
| | e | 疾病级一致性：**全部 117 个疾病**的综合评分分布（按评分降序的秩–评分散点，viridis 色标；中位与 0 参考线；仅标注 AT/SCA44/SCA3/SCA10/JBTS9 五个方位参照，不做 top/bottom 截选） | `tables/cross/fidelity_by_disease.tsv` |

说明：原"DNA 修复类 10 基因类群归属"面板（免疫细胞 8/10）按 2026-11-24 评审意见**撤下**——该信息量单一，
由正文第 5 段与图注承载即可；图2 因此为 a/b/c 三面板，b、c 更宽更清晰。

与正文（`docs/执行情况报告_修改后版本.txt`）的对应：图1↔第2–3段；图2↔第5–6段（b 对应第5段细胞类型分布与
小胶质细胞观察，c 对应第5段末"离子通道与突触类更倾向呈现明确发育变化"，a 对应第6段时序层级）；
图3↔第7–8段（a–d 对应第8段 GRM1 与 ATXN3，e 对应第7段疾病级评分）。
管线既有图形（下表各节）仍作为附录/论文素材，不进入正文这 3 张。

## 0. 图形清单总览

| 层 | 目录 | 图 | 数量 |
| --- | --- | --- | --- |
| 物种层（鼠） | `results/v1.0/figures/mouse/` | a1 ALL 趋势热图、b1 分组热图×9 组、b2 轨迹小倍数×9 组、b3 组×基因 ρ 热图、b4 时序模式堆叠条、c1 特异性 dotplot、c2 期×组占比热图、d1 emergence 时间轴 | 108 茎 / 216 文件 |
| 物种层（人） | `results/v1.0/figures/human/` | 同上（前缀 `human_`） | 108 茎 / 216 文件 |
| 跨物种正式层 | `results/v1.0/figures/cross/` | fig3a 富集热图、fig3b 机制类×shape 计数、fig3c 三套峰时轴、fig3d 代表轨迹、x1(=Fig5a) ALL 状态热图、x2(=Fig5b) 分组状态热图、x3(=Fig5c) 人鼠并排轨迹、x4(=Fig5d) ρ 一致性散点、fig5e 保真度表热图 | 9 茎 / 18 文件 |
| 探索性补充 | `results/v1.0/analysis/a4_common_mechanism_{mouse,human}/`、`.../a6_cross_species/` | a4_1 模块热图、a4_2 模块 eigen 轨迹、a4_3 富集、a4_4 跨组一致性；a6_1 模块对应、a6_2 谱相关散点 | 6 茎 ×2 物种（a4）+ 2 茎（a6） |

分页规则（用于定位单个基因）：a1 每页 60 行、按 ALL 组 ρ 降序；b1/b2 每页 12 格（3 列×4 行）、组内按 ρ 降序；c1 行序 = 主组→τ 降序。

---

## 一、第一节「任务完成情况」

| 文字结论 | 图 | 表/日志 |
| --- | --- | --- |
| 空间维 = 8 主组 + UBC 描述组 + ALL 参考组；时间维 = 鼠 12 期 / 人 13 期 | 无需图（口径） | `config/config.yaml`（`cell_groups`、`group_tiers`、`species.*.time_labels`、`periods`） |
| 一键复现、三道人工门通过 | — | `scripts/run_all.sh`；`reviews/gate1_mouse.md`、`gate2_human.md`、`gate3_cross.md` |

## 二、第二节「基因集与临床异质性」

| 文字结论 | 图 | 表 |
| --- | --- | --- |
| 177 基因 / 188 疾病行 / 186 病种；core 108 / extended 54 / peripheral 15 | — | `genes/ataxia_genes.tsv`、`genes/ataxia_gene_disease.tsv` |
| 8 个异质性维度分布（机制 8、遗传 4、起病档 4、突变 7、病程 5、病理 24、临床分组 8、polyQ 7） | **建议新做一张"异质性构成"总览图**（目前无对应图） | 同上两表按列汇总 |
| 验收 5 项通过、2 项顺延 v1.1 | — | `genes/curation_log.md` §4 处置表；s00b 验收自检输出 |

## 三、第三节「啮齿类图谱主体结果」

| 文字结论 | 图（首选 → 备选） | 表 |
| --- | --- | --- |
| 立方 120 格 / 可行 92 格；不可行记缺失不插值 | 物种报告 §2 可行性矩阵表（**建议后续做成 S1 热图**） | `results/v1.0/logs/s01_feasibility_mouse.tsv`、`intermediate/mouse/cube_bin.tsv` |
| 136 个显著「基因×细胞组」趋势（q<0.05） | `figures/mouse/mouse_a1_ALL_p1.pdf`（ALL 组趋势热图，ρ 降序）→ `mouse_b3_p1/p2.pdf`（组×基因 ρ 热图，带显著性星号） | `tables/mouse/trend_by_group.tsv` |
| 三锚点方向全对（ITPR1↑ / MBP↑ / ATOH1↓） | ITPR1：`mouse_b2_Purkinje_p2.pdf` **第 2 格**；人侧 `human_b2_Purkinje_p1.pdf` **第 2 格**。**MBP / ATOH1 无图**（二者属 sanity marker，不在统计 panel，故不进 a1/b1/b2）→ 数值见 `logs/s02_report_mouse.tsv` 与物种报告 §5 | `logs/s02_report_{mouse,human}.tsv` |
| 7 类时序模式分布（flat 727 / endpoint_high 93 / monotonic_up 79 / monotonic_down 51 / window_high 35 / transient_peak 10 / 缺失 75） | `figures/mouse/mouse_b4.pdf`（按组堆叠条） | `tables/mouse/temporal_patterns.tsv` |
| 6 个强特异基因（τ=0.865–0.954） | `figures/mouse/mouse_c1.pdf`（dotplot，行序即 τ 降序：ITPR1、DAB1、CYP27A1、SLC1A3、RELN、PRKCG）→ `mouse_c2_p1/p2.pdf`（期×组占比热图） | `tables/mouse/specificity_global.tsv` |
| 首次出现期与峰期分布 | `figures/mouse/mouse_d1.pdf`（emergence 时间轴） | `tables/mouse/emergence.tsv`、`peak_by_group.tsv` |

## 四、第四节「共性规律的正式检验」（报告核心）

| 文字结论 | 图 | 表 |
| --- | --- | --- |
| 87/101 基因广谱表达；机制类×模块富集 0 格显著 | `figures/cross/fig3a_enrichment_heatmap.pdf`（人鼠并排 log2FE 热图） | `tables/cross/module_assignment.tsv`（`g_star=broad` 计数）、`enrichment_class_x_module.tsv` |
| 唯一跨物种重复的显著格：channel_synapse × flat 缺失（鼠 −1.11/q=0.021；人 −2.31/q=0.0042） | `figures/cross/fig3b_shape_mosaic.pdf`（机制类×shape 计数，人鼠并排） | `tables/cross/enrichment_class_x_module.tsv`（`test=class_x_shape`）、`fisher_class_x_shape.tsv` |
| **机制类峰时层级**：dna_repair 0.5 → 中段 6.0 → channel_synapse 14.5；KW p=1.3×10⁻³；4 对显著、delta 0.52–0.76 | `figures/cross/fig3c_peak_ridge.pdf`（**三套轴并排**，中间面板 = human_only，右面板 = mouse_only＝正文采信轴）→ `fig3d_representative_traj.pdf`（代表基因轨迹） | `tables/cross/peak_rank_by_class.tsv`（逐基因秩）、`peak_rank_pairwise.tsv`（两两 MW + Cliff's delta，按 `axis` 筛选） |

## 五、第五节「人源数据补充与物种保守性」

| 文字结论 | 图 | 表 |
| --- | --- | --- |
| 人侧图谱建成（107 基因全可测、185 显著趋势） | `figures/human/human_a1_ALL_p1/p2.pdf`、`human_b3_p1/p2.pdf`、`human_c1.pdf`、`human_d1.pdf` | `tables/human/*.tsv`（7 张）、`reports/human/species_report.md` |
| 26 个共享 dev_state 对齐（不做年龄换算） | `figures/cross/x1_state_heatmap_ALL.pdf`（= Fig5a，人鼠双联、物种内 z、独立色标）→ `x2_state_heatmap_groups.pdf`（= Fig5b，Purkinje/GranuleLineage/OligoLineage/Astro 四组） | `tables/cross/state_profiles.tsv`（10,100 行长表，含 `n_cells`/`feasible`/`consensus_rank`） |
| 方向同向率 70.0% | `figures/cross/x4_rho_concordance.pdf`（= Fig5d，含象限计数与同向率） | `tables/cross/cross_species_gene_group.tsv`（808 行） |
| 轨迹形状相关中位 +0.563 | `figures/cross/x3_paired_trajectories.pdf`（= Fig5c；当前选定 ITPR1、ATOH1 + traj_r 最高 3 基因 PNPLA6/PTF1A/SLC1A3，均 Purkinje） | 同上（`traj_r` 列） |
| 保真分：high 76 / medium 25 / low 6；按机制类 cytoskeleton 0.989 > channel_synapse 0.792 > … > proteostasis_RNA 0.292 | `figures/cross/fig5e_fidelity_table.pdf`（疾病级保真度表热图，前 15 + 后 10） | `tables/cross/conservation_gene.tsv`（107 行基因级）、`conservation_gene_by_group.tsv`（856 行）、`fidelity_by_disease.tsv`（186 行） |

## 六、第六节「GRM1 与 SCA3 专段」（与专项内其他课题呼应）

### GRM1（SCA44 / SCAR13）

| 文字结论 | 图（精确格位） |
| --- | --- |
| 鼠 6 组显著单调上升 | `mouse_a1_ALL_p1.pdf` **第 4 行**；`mouse_b2_Purkinje_p1.pdf` **第 5 格**；`mouse_b2_GranuleLineage_p1.pdf` **第 4 格**；`mouse_b2_DN_p1.pdf` **第 7 格**；`mouse_b2_Interneuron_p1.pdf` **第 10 格**；`mouse_b2_Progenitor_p1.pdf` **第 4 格** |
| 人 5 组显著上升 | `human_a1_ALL_p1.pdf` **第 13 行**；`human_b2_Purkinje_p1.pdf` **第 9 格**；`human_b2_GranuleLineage_p1.pdf` **第 8 格**；`human_b2_DN_p3.pdf` **第 9 格**；`human_b2_Interneuron_p1.pdf` **第 2 格** |
| 高峰在最晚档（统一秩 14.5） | `figures/cross/fig3c_peak_ridge.pdf`（channel_synapse 行最右端的点之一） |
| 人鼠轨迹形状相关 0.991（Purkinje）/ 0.986（GranuleLineage） | **当前 x3 图未包含 GRM1** → 若后续报告要展示，需在 `scripts/s04_figures_crossspecies.py` 的 x3 选基因列表中加入 `GRM1`（一行改动，重跑 s04 约 15 秒） |
| 综合保真分 0.714（high） | `figures/cross/fig5e_fidelity_table.pdf`（SCA44 行）；`tables/cross/conservation_gene.tsv` GRM1 行 |
| 稳健性（多口径一致上升） | `results/v1.0/sensitivity/robust_core_genes.tsv`（GRM1：人 3 组 / 鼠 4 组，方向 up）※正文未引用，备查 |

### ATXN3（SCA3 / Machado-Joseph，临床课题主要病种）

| 文字结论 | 图（精确格位） |
| --- | --- |
| 鼠：全小脑 ρ=−0.762、颗粒谱系 ρ=−0.699 | `mouse_a1_ALL_p2.pdf` **第 35 行**；`mouse_b2_GranuleLineage_p7.pdf` **第 5 格** |
| 人：全小脑 ρ=−0.703 | `human_a1_ALL_p2.pdf` **第 38 行** |
| 胚胎期最高（统一峰时秩 2.5，鼠峰 E12.5） | `figures/cross/fig3c_peak_ridge.pdf`（polyQ_repeat 行最左端的点）；`mouse_b2_GranuleLineage_p7.pdf` 第 5 格可见 E12.5 峰形 |
| 4 个可判组方向一致 | `figures/cross/x4_rho_concordance.pdf`（左下"−−"象限内的 ATXN3 点，按组着色） |
| 综合保真分 0.500（medium） | `figures/cross/fig5e_fidelity_table.pdf`（SCA3 行）；`tables/cross/conservation_gene.tsv` ATXN3 行 |
| 四个外部数据集均可测（为疾病投影预留） | — （`genes/ataxia_genes.tsv` 的 `in_human_sepp/in_mouse_sepp/in_aldinger/in_castonguay` 全为 true） |

## 七、第七节「局限与下一阶段」

| 文字结论 | 现状 |
| --- | --- |
| 疾病组学未接入（疾病投影层、共享靶点网络） | **无图无表**：v3 层 b00/b01、n00–n04 未实现；疾病数据已在 `/mnt/data/xyy/Ataxia_CB` 完成盘点（`outputs/disease_registry.tsv`，308 KB） |
| 外部验证待补（BrainSpan / GTEx / 模型鼠 DEG） | **无图**：`external/` 三文件缺席，s10 未实现（门 1/门 2 已就此记 WARN 并判非阻断） |
| 基因集 v1.1（来源核实 #2、双标注 #6、CL-024） | — 见 `genes/curation_log.md` §4「v1.1 待办」 |
| 第一标注论文 | 图骨架已具备（a1/b1–b4/c1/c2/d1 + fig3a–d + fig5a–e）；缺 s05 主表导出、s08 起病关联、s10 外部验证、s11 图表组装 |

---

## 附：后续若要出图，建议补做的 4 张（目前缺）

1. **基因集临床异质性总览**（第二节）：机制类 × 遗传方式 × 起病档 × 临床分组的组合条图/旭日图——最能一眼说明"异质性"三个字，专项办友好。
2. **可行性矩阵热图（S1）**（第三节）：组×bin 细胞数热图 + NA 斜纹，用于说明"缺失不插值"的严谨性。
3. **机制类峰时序示意图**（第四节）：把 fig3c 简化成"机制类 → 发育时间轴"的横向条带图（含 Cliff's delta 标注），比 ridge 点图更适合汇报。
4. **GRM1 / ATXN3 双基因专题图**（第六节）：两基因在人鼠、各细胞组的轨迹并排（含 dev_state 对齐版），直接服务与其他课题的呼应。
