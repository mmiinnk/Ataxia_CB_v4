# review 门 2 — human(spec §9 四条清单; s06 自动生成自动项, 人工项勾选后 commit)

自动项总状态: 通过(FAIL 0 / WARN 6 / 锚点 19 / marker 1 条)

> **签署记录**：reviewer **XYY** 于 **2026-11-24** 判读后裁决：**四条清单项同意通过，但 §11 未决条目 7
> （人成人 4 个 bin 42/44/46/52y 各 1 供体、生物学同期却被 Spearman 当作 4 个有序时间点）不同意按现状通过**，
> 要求把成人段合并为**单一时间点**后复检 → 故本门结论为 **[ ] 过门 / [x] 整改后复检**。
> 勾选与签名由 Claude Code **转录自该对话**（判定为人工做出，脚本无跳过开关，硬约束 9）。
> **整改项（唯一）**：`sensitivity.collapse_adult_bins` —— 人 42/44/46/52y → 单一 `AD` bin（4 供体）。
> 整改依据已量化：该变体已跑（`results/v1.0/sensitivity/collapse_adult_human/`），
> 合并后**显著趋势行 256 → 185（−71 行，−28%）**，即默认口径下约三成"显著发育趋势"来自把同期成人
> 当成有序时间点；故 reviewer 的科学判断得到数据支持。待裁决：是否把合并口径升为人侧**主分析口径**
> （需同步修订 spec §11-7"主分析保持 16 bin"一句，并重跑人侧 s01–s06 + 重签本门）。
> 自动项内容 = gene_set v1.0 重跑结果（`results/v1.0/`）。本仓库尚未 `git init`，"commit" 一步暂以本文件留痕代替。

## (1) logs 无 assert 失败且 WARNING 已逐条看过
- 自动: ✓ 自动判定通过 — FAIL 0 条
- WARN 6 条(判读后勾选):
| 来源 | check | observed | note |
|---|---|---|---|
| s00_report | castonguay: normal 供体/PC 锚点 §4.3 | {'n_donors_normal': 93, 'n_pc_normal': 31089, 'n_donors_ge30pc': 89, 'median_pc_ge30': 257} → {'n_donors_normal': 93, 'n | median_pc_ge30: 实测 264.0 != 锚点 257 |
| s00_report | mouse: 鼠 dropout 现场重算 §3.12 | 锚点±容差内(WARN=0) → WARN=2, ABSENT=0 | 详见 logs/s00_dropout.tsv |
| s00_report | external: external 文件 brainspan | /home/victor/Projects/Ataxia_CB_v4/external/brainspan_cerebellum_rpkm.tsv → 缺席 | s10 对缺失文件只 WARN 并出占位 panel(§11 未决 6) |
| s00_report | external: external 文件 gtex | /home/victor/Projects/Ataxia_CB_v4/external/gtex_v8_cerebellum_tpm.tsv → 缺席 | s10 对缺失文件只 WARN 并出占位 panel(§11 未决 6) |
| s00_report | external: external 文件 model_mouse_deg | /home/victor/Projects/Ataxia_CB_v4/external/model_mouse_deg_signatures.tsv → 缺席 | s10 对缺失文件只 WARN 并出占位 panel(§11 未决 6) |
| s01_report_human | marker 校验(spec §9 清单) | 期望组=argmax(或≥50%argmax) → FAIL=0(清单内) / 3(清单外), ABSENT=4 |  |

- [x] 人工确认: WARNING 已逐条看过且非阻断 —— XYY 2026-11-24(6 条: 前 5 条与门 1 同源已判读; 第 6 条人侧 marker 清单内 FAIL=0)

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
| s01_report_human | 趋势锚点 MBP×OligoLineage | rho + → rho=+0.964(n=11 bin, 粗算) | PASS |
| s01_report_human | 趋势锚点 ATOH1×GranuleLineage | rho - → rho=-0.601(n=16 bin, 粗算) | PASS |
| s02_report_human | 趋势锚点 ITPR1×Purkinje(正式) | rho + | rho_bin=+1.000 q_bin=0 status=ok | PASS |
| s02_report_human | 趋势锚点 MBP×OligoLineage(正式) | rho + | rho_bin=+0.964 q_bin=1.85e-06 status=ok | PASS |
| s02_report_human | 趋势锚点 ATOH1×GranuleLineage(正式) | rho - | rho_bin=-0.601 q_bin=0.0139 status=ok | PASS |

- [x] 人工确认: 锚点数值与 spec §4 逐项核对 —— XYY 2026-11-24(19 条 FAIL=0; human 形状 (180956,27410) 精确、9 组细胞数精确匹配、size_factor r=1.00000; 唯一 WARN 为 castonguay 中位数, 非阻断)

## (3) marker 校验通过
- 自动: ⚠ 有 WARN, 待人工裁决(FAIL 0)(marker 行 1 条; FAIL→不过, WARN→人工裁决)
| 来源 | check | observed | level |
|---|---|---|---|
| s01_report_human | marker 校验(spec §9 清单) | 期望组=argmax(或≥50%argmax) → FAIL=0(清单内) / 3(清单外), ABSENT=4 | WARN |

- [x] 人工确认: ITPR1/CA8→Purkinje, GABRA6/CBLN3→GranuleLineage, MBP/PLP1→OligoLineage, AQP4→Astro, SOX2/MKI67→Progenitor, PTPRC→Immune(鼠免疫豁免如上)
      —— XYY 2026-11-24。判读注(人侧): 清单内 FAIL=0。清单外 FAIL 3 条已判读为非阻断——
      PCP2(argmax=Astro, PC 组均值仅 0.024)、HOMER3(argmax=Interneuron, 0.004) 系**人侧 PC 数据仅覆盖产前
      (CS18–NB)、胎龄 PC 未成熟**所致(成熟 PC 标志尚未上调); 对照鼠侧两者均 PASS(PCP2 在 Purkinje 0.668),
      故属覆盖面缺陷而非细胞注释错误。PAX6 argmax=GranuleLineage 与鼠侧一致(真实生物学: 菱形唇/颗粒祖细胞)。
      ABSENT 4 条 = TBR1(DN)、EOMES/LMX1A(UBC)、SOX2(Progenitor) —— 人侧 var 无该基因或按 species_valid 过滤,
      故本勾选框文字中 SOX2→Progenitor、CA8 在人侧不适用(CA8 在 var 中且 specific→Purkinje τ=0.962, 但非本清单判定项)。
      **遗留影响(须写入报告)**: 人侧无法验证出生后 PC 特异性, 该段只能由鼠承担; 对 v3 PC 网络层(b00/b01、n00–n04)
      若拟用人侧 PC 数据, 必须先接受"人 PC 仅产前、成熟标志不表达"这一事实。

## (4) 趋势锚点方向正确
- 自动: ✓ 自动判定通过(正式趋势锚点 3 条; §9(4) 期望: ITPR1×Purkinje ρ>0; MBP×OligoLineage ρ>0; ATOH1×GranuleLineage ρ<0; 方向不符记 WARN 亦判不过)
| check | observed | level |
|---|---|---|
| 趋势锚点 ITPR1×Purkinje(正式) | rho + | rho_bin=+1.000 q_bin=0 status=ok | PASS |
| 趋势锚点 MBP×OligoLineage(正式) | rho + | rho_bin=+0.964 q_bin=1.85e-06 status=ok | PASS |
| 趋势锚点 ATOH1×GranuleLineage(正式) | rho - | rho_bin=-0.601 q_bin=0.0139 status=ok | PASS |

- 锚点来源(config trend_anchors.source):
  - ITPR1×Purkinje(+): Gold 2003 Neuron(RORα 依赖 PC 成熟程序); Aldinger2021 胎窗实测 rho=+0.643(w9–21, H-PC, n=8 周)
  - MBP×OligoLineage(+): Zečević 1998 Dev Brain Res(人脑 MBP 随髓鞘化上升); Ozarkar 2025 Front Neuroanat(小脑 P10–P23 最快上升); 髓鞘化为出生后事件→Aldinger 胎窗不适用(rho=-0.143, 记 N/A_WINDOW)
  - ATOH1×GranuleLineage(-): Haldipur 2022 Annu Rev Neurosci(RL/EGL 祖细胞标记, 分化下调); Aldinger2021 胎窗实测 rho=-0.790 + 祖(RL+GCP)>已分化(GN)横向对照 75% 周阳性
- [x] 人工确认: 方向与量级合理 —— XYY 2026-11-24(ITPR1×Purkinje +1.000 / MBP×OligoLineage +0.964 / ATOH1×GranuleLineage −0.601, 方向全对)
      判读注: ITPR1 ρ=+1.000 系**产前 8 个 bin 的完美单调**(q_bin=0 来自渐近近似, 8 bin > exact_perm_max_bins=6),
      不可解读为"全生命周期上升"; ATOH1 弱于鼠侧(−0.907) 因人侧轴延伸至成人、ATOH1 长期低位。
      **附条件**: 凡跨越成人 4 bin 的趋势, 其量级受 §11-7 影响 → 见本门整改项, 合并成人口径后须复核。

---
reviewer: XYY  日期: 2026-11-24  结论: [ ] 过门 / [x] 整改后复检
(签署方式: 对话中逐条判读后裁决, 由 Claude Code 转录; gene_set v1.0, 产物 results/v1.0/)
整改项(1 条, 详见顶部签署记录): §11 未决 7 —— 人成人 4 bin(42/44/46/52y, 各 1 供体)生物学同期,
不得当作 4 个有序时间点; 须合并为单一 AD 时间点(`collapse_adult_bins`)后复检。
四条清单项本身均已同意通过(见上勾选与判读注)。整改完成后由 reviewer 复核, 再把结论改为「过门」。
