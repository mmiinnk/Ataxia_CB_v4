# review 门 3 — cross/spec §9(s04+s07+s09 自动项; 人工项勾选后 commit)

自动项总状态: 通过(FAIL 0 / WARN 1 / 日志行 26)

> **签署记录**：reviewer **XYY** 于 **2026-11-26** 判读后裁决：**四条清单项全部同意，过门**。
> 本次为可靠性修复（commit c05c594）后重签：原 2026-11-24 签署件归档为
> `gate3_cross.2026-11-24_修复前.md`；重签材料见 `gate3_重签材料_2026-11-26.md`。
> 勾选与签名由 Claude Code **转录自该对话**（判定为人工做出：用户 2026-11-26 明确"全部确认，
> reviewer XYY"；脚本无跳过开关，硬约束 9）。
> **关键口径确认（材料 §3）**：
> · **B 跨物种峰时限定**：修复恢复 30 个人峰=AD 基因后，三套峰轴一致率 79.2%→50.0%，人侧 DNA 修复类
>   峰呈双峰（5 早期 + 5 成年）。原"两端结论三套轴一致"不再成立。**正文时序层级以 mouse_only 为权威轴**
>   （KW p=1.3×10⁻³、4 对显著不变），明确声明峰时排序跨物种不全自动保守；共表达模块层面仍保守（a6）。
> · **C AD 峰秩 interim**：合并 AD bin 赋质心秩 14.5（与鼠 adult 映射一致）接受为 interim 口径；
>   v1.1 须与 stage_map 按 13-bin 重记统一。
> · **A NA 政策**：6 个鼠缺失基因疾病由假 0 转 NA（不可测，非阴性），疾病级可评分 117→111，中位 0.714 不变。

## (1) logs 无 assert 失败且 WARNING 已逐条看过
- 自动: ✓ 通过 — FAIL 0 条
- WARN 1 条(判读后勾选):
| source | check | observed | note |
|---|---|---|---|
| s07_report | 三套峰轴方向一致性 | 12/24 = 50.0% | spec §6: 三者方向一致才可写正文; 不一致时**以 mouse_only 为准**(单一物种、无映射混杂) |

- [x] 人工确认: WARNING 已逐条看过且非阻断 —— XYY 2026-11-26（唯一 WARN = 三套峰轴一致率 50.0%，
      判为非阻断；处置：正文以 mouse_only 为准 + 声明峰时排序跨物种不全自动保守，见顶部口径 B）

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
| s07_report | 统计 panel | tiers_main & in_panel(load_panel) | 107 基因 | PASS |
| s07_report | 模块指派 human | status=ok 者进检验 | 107/107 ok; 其余 {'ok': np.int64(107)} | PASS |
| s07_report | 模块指派 mouse | status=ok 者进检验 | 101/107 ok; 其余 {'ok': np.int64(101), 'missing_gene': np.int64(6)} | PASS |
| s07_report | CL-027 残差类排除 human | 排除 ['other', 'unknown'] | 3 基因: ['CCDC88C', 'MME', 'WWOX']; 入富集 104 | PASS |
| s07_report | NA 行排除 human | 模块/shape 为 NA 者不进检验(NA 政策) | class×module 入 104/104; class×shape 入 104/104(shape=NA 0 行, 多为 not_expressed 的掩膜伪标签) | PASS |
| s07_report | 富集显著格 human | q<0.05 且 obs≥3 | 0 格 | PASS |
| s07_report | CL-027 残差类排除 mouse | 排除 ['other', 'unknown'] | 3 基因: ['CCDC88C', 'MME', 'WWOX']; 入富集 98 | PASS |
| s07_report | NA 行排除 mouse | 模块/shape 为 NA 者不进检验(NA 政策) | class×module 入 98/98; class×shape 入 96/98(shape=NA 2 行, 多为 not_expressed 的掩膜伪标签) | PASS |
| s07_report | 富集显著格 mouse | q<0.05 且 obs≥3 | 1 格; 最强 class_x_shape channel_synapse×flat log2FE=-1.11 q=0.021 | PASS |
| s07_report | 峰时轴 mixed | KW + 两两 MW(BH) + Cliff's delta | n=103 基因 / 7 类, KW p=0.0111, 显著对(q<0.05) 2 | PASS |
| s07_report | 峰时轴 human_only | KW + 两两 MW(BH) + Cliff's delta | n=98 基因 / 7 类, KW p=0.0266, 显著对(q<0.05) 2 | PASS |
| s07_report | 峰时轴 mouse_only | KW + 两两 MW(BH) + Cliff's delta | n=98 基因 / 7 类, KW p=0.00131, 显著对(q<0.05) 4 | PASS |
| s07_report | 三套峰轴方向一致性 | 两两类别中位差符号 vs mouse_only 一致率 | 12/24 = 50.0% | WARN |
| s07_report | 产物 | tables/cross + figures/cross | module_assignment(214 行), enrichment(287), fisher(77), peak_rank(107), pairwise(63); 图 4 张 ×2 格式 | PASS |

- [x] 人工确认: dev_state 对齐口径(§4.4/§4.5A)、CL-027 残差类排除、NA 政策均无异议 —— XYY 2026-11-26
      判读注: 跨物种时间轴用 26 个共享 dev_state（按细胞分化进度对齐，不做人鼠年龄换算）；残差类
      other/unknown 排除于富集网格（人鼠各 3 基因，改单基因描述）；不可行格一律 NA、不插值不填 0。
      **NA 政策修正接受**：6 个鼠缺失基因 max_group_match 记 NA（原误为 0），疾病级可评分 117→111。
      **AD 峰秩 14.5 interim 口径接受**（顶部口径 C），v1.1 与 stage_map 统一。

## (3) marker 校验
- 自动: 跨物种层不重复 marker 校验(marker 属物种内 check, 已在门 1/门 2 判读)。
- [x] 人工确认: 认可本层免检(如需, 回门 1/门 2 复核) —— XYY 2026-11-26（marker 属物种内 check，
      门 1/2 已签；本次 sanity.marker_table 补全 UBC 描述组 marker 校验，主组 verdict 不变，
      物种内 marker 日志待下次 s01 重跑更新）

## (4) fig3a 富集与 fig5 保真度统计(门 3 焦点)
- fig3a 显著富集格(q<0.05 且 obs≥3): **1** 格 — 最强: mouse channel_synapse×flat log2FE=-1.11 q=0.021
- 保真度: 基因级 107 行; 证据等级 {'high': 76, 'medium': 25, 'none': 6}; 保真分中位 +0.714(区间 [−0.3, 1.0])
- 疾病级 186 行（其中 **111 行可评分**，fidelity 非 NA）; 可评分中位保真分 +0.714
- [x] 人工确认: fig3a/fig3c/fig5a–e 的统计与图面一致、结论表述与证据等级相符 —— XYY 2026-11-26
      判读注（结论强度分级，报告须照此表述）:
      ① fig3a 机制类×模块富集 **0 格显著**（人/鼠均然）→ 信息性结果，只可写"未发现机制类偏好某细胞类型×时期模块"，
         **不得**写"无共性机制"；根因是多数基因 G*=broad（广谱表达）。
      ② 唯一跨物种重复的显著格 = `channel_synapse × shape=flat` 的**缺失**（人 log2FE −2.31 q=0.0042 / 鼠 −1.11 q=0.021）
         → 可写"通道/突触类基因几乎不呈平坦表达，均有真实发育动态"。
      ③ fig3c 峰时轴 KW：mouse_only p=1.3×10⁻³（4 对显著，**正文权威**）、mixed p=0.0111（2 对）、human_only p=0.0266（2 对）
         → 机制类发育时序差异可写，**以 mouse_only 为准**；峰时排序跨物种不全自动保守（顶部口径 B）。
      ④ fig5 保真分区间 [−0.3, 1.0]，**不是百分比**；fig3e 展示 111 个可评分疾病；low/none（6 基因）= 证据不足，非阴性结论。
      ⑤ 高证据（76 基因）口径机制类保真中位：cytoskeleton_transport 0.989 / channel_synapse 0.792 / proteostasis_RNA 0.292（不变）。

---
reviewer: XYY  日期: 2026-11-26  结论: [x] 过门 / [ ] 整改后复检
(自动项生成 2026-11-25 18:29，签署转录 2026-11-26; gene_set v1.0; 产物 results/v1.0/)
(签署方式: reviewer 在对话中裁决"全部确认，reviewer XYY"，由 Claude Code 转录; 修复后重签，原 2026-11-24 件归档)
(三道门现状: 门 1 鼠过门 · 门 2 人过门 · 门 3 跨物种过门〔2026-11-26 修复后重签〕; s10/s11 另需门 4 与论文阶段)
