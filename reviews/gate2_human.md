# review 门 2 — human(spec §9 四条清单; s06 自动生成自动项, 人工项勾选后 commit)

自动项总状态: 通过(FAIL 0 / WARN 7 / 锚点 19 / marker 1 条)

> **签署记录（复核改判）**：reviewer **XYY**，2026-11-24。
> 首签（同日）结论为「整改后复检」，唯一整改项 = §11 未决 7（人成人 4 bin 42/44/46/52y 各 1 供体、
> 生物学同期，却被 Spearman 当作 4 个有序时间点）。
> **整改已完成并验证**：`collapse_adult_bins=true` 升为人侧**主口径**（config.yaml + spec §11-7 + s01
> `cells_default`/`bins_default` 判据），人侧 s01–s06 已按 **13 bin** 重跑；旧 16 bin 口径降为敏感性
> 对照变体 `adult_16bin`（`run_all.sh sensitivity adult_16bin`）。
> 整改效果实测：显著趋势行 **256→185（−28%）**，保留行**方向一致率 100%**（无一例翻转）；
> 旧口径可复现新口径 **94.6%** 的结论，但多出 **81 行**仅靠"成人排序"才显著 → 即被砍掉的可疑信号。
> 单供体可行 bin 由 52 降至 **29**。reviewer 复核后裁决：**过门**。
> 勾选与签名由 Claude Code **转录自对话**（判定为人工做出，脚本无跳过开关，硬约束 9）。
> 本次自动项 = 13 bin 新主口径重跑结果（FAIL 0 / **WARN 7**）；较首签多的 1 条为
> 「§4.2 可行性锚点表按 16 bin 记录、与 13 bin 不可比 → 跳过比对」，判读为非阻断，
> 并列为 **v1.1 待办**（按新主口径重记 §4.2 锚点表）。
> 首签版本留档 `gate2_human.2026-11-24_整改后复检.md`；本仓库尚未 `git init`，"commit" 一步以本文件留痕代替。

## (1) logs 无 assert 失败且 WARNING 已逐条看过
- 自动: ✓ 自动判定通过 — FAIL 0 条
- WARN 7 条(判读后勾选):
| 来源 | check | observed | note |
|---|---|---|---|
| s00_report | castonguay: normal 供体/PC 锚点 §4.3 | {'n_donors_normal': 93, 'n_pc_normal': 31089, 'n_donors_ge30pc': 89, 'median_pc_ge30': 257} → {'n_donors_normal': 93, 'n | median_pc_ge30: 实测 264.0 != 锚点 257 |
| s00_report | mouse: 鼠 dropout 现场重算 §3.12 | 锚点±容差内(WARN=0) → WARN=2, ABSENT=0 | 详见 logs/s00_dropout.tsv |
| s00_report | external: external 文件 brainspan | /home/victor/Projects/Ataxia_CB_v4/external/brainspan_cerebellum_rpkm.tsv → 缺席 | s10 对缺失文件只 WARN 并出占位 panel(§11 未决 6) |
| s00_report | external: external 文件 gtex | /home/victor/Projects/Ataxia_CB_v4/external/gtex_v8_cerebellum_tpm.tsv → 缺席 | s10 对缺失文件只 WARN 并出占位 panel(§11 未决 6) |
| s00_report | external: external 文件 model_mouse_deg | /home/victor/Projects/Ataxia_CB_v4/external/model_mouse_deg_signatures.tsv → 缺席 | s10 对缺失文件只 WARN 并出占位 panel(§11 未决 6) |
| s01_report_human | 组×bin 可行性矩阵 §4.2 | - → 人侧主口径已合并成人 bin(collapse_adult_bins=true; 2026-11-24 用户裁决, 门 2 整改项): spec §4.2 可行性锚点表按 16 bin 记录, 与合并后 13 bin 不可比 → 跳 |  |
| s01_report_human | marker 校验(spec §9 清单) | 期望组=argmax(或≥50%argmax) → FAIL=0(清单内) / 3(清单外), ABSENT=4 |  |

- [x] 人工确认: WARNING 已逐条看过且非阻断 —— XYY 2026-11-24 **复核**（7 条）：前 5 条与门 1 同源（castonguay 中位数 264≠257、鼠 dropout Ptprc/Csf1r 4.31%>4、external 三文件缺席）已判读非阻断；第 6 条 = §4.2 可行性锚点表按 16 bin 记录、与新主口径 13 bin 不可比 → 跳过比对（**整改的直接连带项**，非数据问题，v1.1 按 13 bin 重记锚点表）；第 7 条 = 人侧 marker 清单内 FAIL=0（详见勾选(3)）

## (2) §4 锚点全部匹配
- 自动: ⚠ 有 WARN, 待人工裁决(FAIL 0)(锚点行 19 条; FAIL→不过, WARN→人工裁决)
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

- [x] 人工确认: 锚点数值与 spec §4 逐项核对 —— XYY 2026-11-24 **复核**（19 条 FAIL=0）：human 形状 (180956, 27410) 精确、size_factor↔X 行和 r=1.00000、**§4.1 9 组细胞数在 13 bin 合并口径下仍精确匹配 PASS**（s01 已拆出 `cells_default` 判据——该锚点按组计数、与 bin 切法无关，过去随 `default_mode` 一并跳过属过宽，本轮修正）；mouse/aldinger/castonguay 同门 1

## (3) marker 校验通过
- 自动: ⚠ 有 WARN, 待人工裁决(FAIL 0)(marker 行 1 条; FAIL→不过, WARN→人工裁决)
| 来源 | check | observed | level |
|---|---|---|---|
| s01_report_human | marker 校验(spec §9 清单) | 期望组=argmax(或≥50%argmax) → FAIL=0(清单内) / 3(清单外), ABSENT=4 | WARN |

- [x] 人工确认: ITPR1/CA8→Purkinje, GABRA6/CBLN3→GranuleLineage, MBP/PLP1→OligoLineage, AQP4→Astro, SOX2/MKI67→Progenitor, PTPRC→Immune(鼠免疫豁免如上)
      —— XYY 2026-11-24 **复核**（与首签同判，13 bin 重跑后 marker 结果不变: 清单内 FAIL=0 / 清单外 3 / ABSENT 4）。
      清单外 3 条判读: PCP2(argmax=Astro, PC 组均值仅 0.024)、HOMER3(argmax=Interneuron, 0.004) 系**人侧 PC
      仅覆盖产前(CS18–NB)、胎龄 PC 未成熟**所致(成熟 PC 标志尚未上调); 对照鼠侧两者均 PASS(PCP2 在 Purkinje 0.668)
      → 属覆盖面缺陷, 非细胞注释错误。PAX6 argmax=GranuleLineage 与鼠侧一致(真实生物学: 菱形唇/颗粒祖细胞)。
      ABSENT 4 条 = TBR1(DN)、EOMES/LMX1A(UBC)、SOX2(Progenitor): 人侧 var 无该基因或按 species_valid 过滤(CL-018),
      故本条文字中 SOX2→Progenitor 在人侧不适用。
      **遗留影响(须写入报告)**: 人侧无法验证出生后 PC 特异性, 该段只能由鼠承担; v3 PC 网络层(b00/b01、n00–n04)
      若拟用人侧 PC 数据, 必须先接受"人 PC 仅产前、成熟标志不表达"这一事实。

## (4) 趋势锚点方向正确
- 自动: ✓ 自动判定通过(正式趋势锚点 3 条; §9(4) 期望: ITPR1×Purkinje ρ>0; MBP×OligoLineage ρ>0; ATOH1×GranuleLineage ρ<0; 方向不符记 WARN 亦判不过)
| check | observed | level |
|---|---|---|
| 趋势锚点 ITPR1×Purkinje(正式) | rho + | rho_bin=+1.000 q_bin=0 status=ok | PASS |
| 趋势锚点 MBP×OligoLineage(正式) | rho + | rho_bin=+0.929 q_bin=0.000863 status=ok | PASS |
| 趋势锚点 ATOH1×GranuleLineage(正式) | rho - | rho_bin=-0.564 q_bin=0.0448 status=ok | PASS |

- 锚点来源(config trend_anchors.source):
  - ITPR1×Purkinje(+): Gold 2003 Neuron(RORα 依赖 PC 成熟程序); Aldinger2021 胎窗实测 rho=+0.643(w9–21, H-PC, n=8 周)
  - MBP×OligoLineage(+): Zečević 1998 Dev Brain Res(人脑 MBP 随髓鞘化上升); Ozarkar 2025 Front Neuroanat(小脑 P10–P23 最快上升); 髓鞘化为出生后事件→Aldinger 胎窗不适用(rho=-0.143, 记 N/A_WINDOW)
  - ATOH1×GranuleLineage(-): Haldipur 2022 Annu Rev Neurosci(RL/EGL 祖细胞标记, 分化下调); Aldinger2021 胎窗实测 rho=-0.790 + 祖(RL+GCP)>已分化(GN)横向对照 75% 周阳性
- [x] 人工确认: 方向与量级合理 —— XYY 2026-11-24 **复核**（13 bin 新主口径）：ITPR1×Purkinje **+1.000**(q=0, n=8 产前 bin)、MBP×OligoLineage **+0.929**(q=8.63e-4, n=8)、ATOH1×GranuleLineage **−0.564**(q=0.045, n=13)，方向全对。
      判读注: 三值较 16 bin 口径略降(+0.964→+0.929、−0.601→−0.564)，系合并后时间点减少导致的功效下降，**方向与显著性均未变**；
      ITPR1 ρ=+1.000 为产前 8 bin 的完美单调(q=0 来自渐近近似，8 bin > exact_perm_max_bins=6)，不可解读为"全生命周期上升"；
      ATOH1 弱于鼠侧(−0.907)因人侧轴延伸至成人、ATOH1 长期低位。首签所附"跨成人趋势待敏感性确认"的条件**已解除**
      (成人段已合并为单一 AD bin；且 `adult_16bin`/`exclude_disease`/`min_cells_50` 三套对照方向一致率均 100%)。

---
reviewer: XYY  日期: 2026-11-24  结论: [x] 过门 / [ ] 整改后复检
(签署过程: 首签「整改后复检」→ 整改项 collapse_adult_bins 升为人侧主口径并重跑 s01–s06 → 同日复核改判「过门」;
 两次判定均由 reviewer 在对话中做出, 由 Claude Code 转录; 首签版本留档 gate2_human.2026-11-24_整改后复检.md)
(口径: gene_set v1.0; 人侧 13 bin 主口径(collapse_adult_bins=true); 产物 results/v1.0/;
 敏感性对照 results/v1.0/sensitivity/{adult_16bin, exclude_disease, adult_16bin+exclude_disease, min_cells_50, min_cells_100}_human)
(同批附带裁决: `bins.min_cells` **维持 30 为主口径** —— 人侧 30↔50 仅差 1 行(保留率 99.5%),
 100 档人侧保留率 53.0%/鼠侧 71.3% 会砍掉近半结论 → 50/100 仅作敏感性对照, 不得作主口径)
(下一步: 门 1、门 2 均已过 → 按硬约束 9 可进跨物种层(s04/s05/s07/s08/s09, 均未实现))
