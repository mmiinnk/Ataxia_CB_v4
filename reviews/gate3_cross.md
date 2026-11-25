# review 门 3 — cross/spec §9(s04+s07+s09 自动项; 人工项勾选后 commit)

自动项总状态: 通过(FAIL 0 / WARN 1 / 日志行 26)

> **签署记录**：reviewer **XYY**，2026-11-24，结论 **过门**，四条清单项全部同意。
> 唯一 WARN（s07「三套峰轴方向一致性 19/24 = 79.2%」）经判读为**非阻断**，处置口径（写入后续报告与图注）：
> 按 spec §6，涉及机制类峰时次序的正文表述**以 `mouse_only` 轴为准**，并明确标注——
> 「两端结论（`dna_repair` 最早、`channel_synapse`/`polyQ_repeat` 最晚）在三套轴下一致；
> 中间类（`cilia_congenital_dev` / `cytoskeleton_transport_glia`）次序不稳定，因 `stage_map` 中段为
> 插值锚点且 P4/P7 同秩、adult 统一记 14.5，粒度较粗」。
> 勾选与签名由 Claude Code **转录自对话**（判定为人工做出，脚本无跳过开关，硬约束 9）。
> 本层范围：s04 + s07 + s09（s05 主表导出 = 论文阶段；s08 起病关联待 CL-024 裁决；s10 待 `external/` 三文件）。
> 三道门现状：门 1（鼠）过门 · 门 2（人）整改后复检→过门 · 门 3（跨物种）过门，均 reviewer XYY 2026-11-24。
> 本仓库尚未 `git init`，标题所称"commit"一步以本文件签署记录留痕代替。

## (1) logs 无 assert 失败且 WARNING 已逐条看过
- 自动: ✓ 通过 — FAIL 0 条
- WARN 1 条(判读后勾选):
| source | check | observed | note |
|---|---|---|---|
| s07_report | 三套峰轴方向一致性 | 19/24 = 79.2% | spec §6: 三者方向一致才可写正文; 不一致时**以 mouse_only 为准**(单一物种、无映射混杂) |

- [x] 人工确认: WARNING 已逐条看过且非阻断 —— XYY 2026-11-24（唯一 WARN = 三套峰轴一致率 79.2%，判为非阻断；处置口径见顶部签署记录：正文以 mouse_only 为准 + 标注中间类次序不稳）

## (2) 跨物种层口径核查(锚点表的跨物种等价项)
- 自动项(逐条见 logs/s04_report.tsv、s07_report.tsv):
| source | check | expected | observed | level |
|---|---|---|---|---|
| s04_report | 跨物种共有统计 panel | core & in_panel 且两物种 cube 均命中 | 101 基因; 仅单物种命中 6: ['ATN1', 'CA8', 'KIAA0586', 'PDYN', 'SAMD9L', 'TCTN2'] | PASS |
| s04_report | ALL 状态剖面可行性 | 每状态 ≥30 核 | 人 26/26 状态可行, 鼠 26/26; 共识序长度 26 | PASS |
| s04_report | 趋势方向同向率(组×基因) | 两侧 /rho_bin/ ≥ 0.2 才可判(spec §8.1) | 可判 487/808 行, 同向 70.0% | PASS |
| s04_report | traj_r 可算率 | ≥3 个共享可行状态 | 299/808 行, 中位 +0.563 | PASS |
| s04_report | 产物 | tables/cross + figures/cross | state_profiles(10100 行), cross_species_gene_group(808 行); 图 x1–x4 ×2 格式 | PASS |
| s04_report | Fig5 对应关系 | spec §8.1: fig5a/b/c 用 s04 x1–x3, d 用 x4 | x1→fig5a, x2→fig5b, x3→fig5c, x4→fig5d; fig5e(保真度表热图)由 s09 产出 | PASS |
| s07_report | 统计 panel | core & in_panel | 107 基因 | PASS |
| s07_report | 模块指派 human | status=ok 者进检验 | 107/107 ok; 其余 {'ok': np.int64(107)} | PASS |
| s07_report | 模块指派 mouse | status=ok 者进检验 | 101/107 ok; 其余 {'ok': np.int64(101), 'missing_gene': np.int64(6)} | PASS |
| s07_report | CL-027 残差类排除 human | 排除 ['other', 'unknown'] | 3 基因: ['CCDC88C', 'MME', 'WWOX']; 入富集 104 | PASS |
| s07_report | NA 行排除 human | 模块/shape 为 NA 者不进检验(NA 政策) | class×module 入 104/104; class×shape 入 104/104(shape=NA 0 行, 多为 not_expressed 的掩膜伪标签) | PASS |
| s07_report | 富集显著格 human | q<0.05 且 obs≥3 | 0 格 | PASS |
| s07_report | CL-027 残差类排除 mouse | 排除 ['other', 'unknown'] | 3 基因: ['CCDC88C', 'MME', 'WWOX']; 入富集 98 | PASS |
| s07_report | NA 行排除 mouse | 模块/shape 为 NA 者不进检验(NA 政策) | class×module 入 98/98; class×shape 入 96/98(shape=NA 2 行, 多为 not_expressed 的掩膜伪标签) | PASS |
| s07_report | 富集显著格 mouse | q<0.05 且 obs≥3 | 1 格; 最强 class_x_shape channel_synapse×flat log2FE=-1.11 q=0.021 | PASS |
| s07_report | 峰时轴 mixed | KW + 两两 MW(BH) + Cliff's delta | n=101 基因 / 7 类, KW p=0.00016, 显著对(q<0.05) 7 | PASS |
| s07_report | 峰时轴 human_only | KW + 两两 MW(BH) + Cliff's delta | n=69 基因 / 7 类, KW p=0.0134, 显著对(q<0.05) 0 | PASS |
| s07_report | 峰时轴 mouse_only | KW + 两两 MW(BH) + Cliff's delta | n=98 基因 / 7 类, KW p=0.00131, 显著对(q<0.05) 4 | PASS |
| s07_report | 三套峰轴方向一致性 | 两两类别中位差符号 vs mouse_only 一致率 | 19/24 = 79.2% | WARN |
| s07_report | 产物 | tables/cross + figures/cross | module_assignment(214 行), enrichment(287), fisher(77), peak_rank(107), pairwise(63); 图 4 张 ×2 格式 | PASS |

- [x] 人工确认: dev_state 对齐口径(§4.4/§4.5A)、CL-027 残差类排除、NA 政策均无异议 —— XYY 2026-11-24
      判读注: 跨物种时间轴用 26 个共享 dev_state(按细胞分化进度对齐, **不做人鼠年龄换算**), 两物种 26/26 状态均 ≥30 核可行;
      残差类 other/unknown 已按 CL-027 排除于富集网格之外(人鼠各 3 基因: CCDC88C/MME/WWOX, 改单基因描述);
      不可行格一律 NA、不插值不填 0; traj_r 需 ≥3 个共享可行状态, 算不出即 NA(不当 0, 故 low 证据 ≠ 不保守)。

## (3) marker 校验
- 自动: 跨物种层不重复 marker 校验(marker 属物种内 check, 已在门 1/门 2 判读)。
- [x] 人工确认: 认可本层免检(如需, 回门 1/门 2 复核) —— XYY 2026-11-24（marker 属物种内 check，已在门 1（鼠）/门 2（人）逐条判读并签署）

## (4) fig3a 富集与 fig5 保真度统计(门 3 焦点)
- fig3a 显著富集格(q<0.05 且 obs≥3): **1** 格 — 最强: mouse channel_synapse×flat log2FE=-1.11 q=0.021
- 保真度: 基因级 107 行; 证据等级 {'high': np.int64(76), 'medium': np.int64(25), 'low': np.int64(6)}; 保真分中位 +0.685(区间 [−0.3, 1.0])
- 疾病级 186 行; 中位保真分 +0.714
- [x] 人工确认: fig3a(富集热图)/fig3c(三套峰轴)/fig5a–e 的统计与图面一致、结论表述与证据等级相符 —— XYY 2026-11-24
      判读注（结论强度分级，报告须照此表述）:
      ① fig3a 机制类×模块富集 **0 格显著**（人/鼠均然）→ 属信息性结果，只可写"未发现机制类偏好某细胞类型×时期模块"，
         **不得**写"无共性机制"；根因是 86/107（人）、88/101（鼠）基因 G*=broad（广谱表达）。
      ② 唯一跨物种重复的显著格 = `channel_synapse × shape=flat` 的**缺失**（人 log2FE −2.31 q=0.0042 / 鼠 −1.11 q=0.021）
         → 可写"通道/突触类基因几乎不呈平坦表达，均有真实发育动态"。
      ③ fig3c 峰时轴 KW：mixed p=1.6e-4（7 对显著）、mouse_only p=1.3e-3（4 对）、human_only p=0.013（0 对）
         → 机制类的发育**时序**差异可写（两端稳健、中间类待定，见顶部处置口径）。
      ④ fig5 保真分区间 [−0.3, 1.0]，**不是百分比**；`fidelity=1.0 且 traj_r=NA` 只表示三个二元指标一致，
         不表示轨迹形状已验证（27 个基因的人侧主组落在 Immune/UBC 薄数据组）；low 证据（6 基因）= 证据不足，非阴性结论。

---
reviewer: XYY  日期: 2026-11-24  结论: [x] 过门 / [ ] 整改后复检
(自动项生成 2026-11-24 04:31; gene_set v1.0; 产物 results/v1.0/{tables,figures,logs}/cross 与 /{human,mouse})
(签署方式: reviewer 在对话中裁决"门3过门, 四条同意, WARN 判为非阻断", 由 Claude Code 转录)
(三道门全部通过 → 按硬约束 9, s01–s09 正式层产物可用于阶段报告; 后续 s10(外部验证)/s11(图表组装)
 另需门 4 与论文阶段工作, 不在本次范围)
