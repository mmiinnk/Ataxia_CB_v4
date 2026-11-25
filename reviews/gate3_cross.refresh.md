# review 门 3 — cross/spec §9(s04+s07+s09 自动项; 人工项勾选后 commit)

自动项总状态: 通过(FAIL 0 / WARN 1 / 日志行 26)

## (1) logs 无 assert 失败且 WARNING 已逐条看过
- 自动: ✓ 通过 — FAIL 0 条
- WARN 1 条(判读后勾选):
| source | check | observed | note |
|---|---|---|---|
| s07_report | 三套峰轴方向一致性 | 12/24 = 50.0% | spec §6: 三者方向一致才可写正文; 不一致时**以 mouse_only 为准**(单一物种、无映射混杂) |

- [ ] 人工确认: WARNING 已逐条看过且非阻断

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

- [ ] 人工确认: dev_state 对齐口径(§4.4/§4.5A)、CL-027 残差类排除、NA 政策均无异议

## (3) marker 校验
- 自动: 跨物种层不重复 marker 校验(marker 属物种内 check, 已在门 1/门 2 判读)。
- [ ] 人工确认: 认可本层免检(如需, 回门 1/门 2 复核)

## (4) fig3a 富集与 fig5 保真度统计(门 3 焦点)
- fig3a 显著富集格(q<0.05 且 obs≥3): **1** 格 — 最强: mouse channel_synapse×flat log2FE=-1.11 q=0.021
- 保真度: 基因级 107 行; 证据等级 {'high': np.int64(76), 'medium': np.int64(25), 'none': np.int64(6)}; 保真分中位 +0.714(区间 [−0.3, 1.0])
- 疾病级 186 行; 中位保真分 +0.714
- [ ] 人工确认: fig3a(富集热图)/fig3c(三套峰轴)/fig5a–e 的统计与图面一致、结论表述与证据等级相符

---
reviewer: ____________  日期: ____________  结论: [ ] 过门 / [ ] 整改后复检
(自动生成 2026-11-25 18:29; gene_set v1.0; 产物 results/v1.0/)
