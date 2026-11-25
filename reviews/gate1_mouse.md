# review 门 1 — mouse(spec §9 四条清单; s06 自动生成自动项, 人工项勾选后 commit)

自动项总状态: 通过(FAIL 0 / WARN 6 / 锚点 19 / marker 1 条)

> **签署记录**：reviewer **XYY** 于 **2026-11-24** 逐条判读 6 条 WARN 后裁决「**过门**，四项全部同意」
> （含：castonguay `median_pc_ge30` 264≠257 记为非阻断、v1.1 复核；PAX6 属清单外、EOMES/LMX1A 为
> UBC marker 不适用）。勾选与签名由 Claude Code **转录自该对话**——判定为人工做出，脚本无跳过开关（硬约束 9）。
> 自动项内容 = gene_set v1.0 重跑结果（`results/v1.0/`），自 `gate1_mouse.refresh.md` 同步；
> draft 期旧清单留档 `gate1_mouse.v0.9.1-draft.md`。本仓库尚未 `git init`，故清单标题所称"commit"
> 一步暂以本文件签署记录留痕代替（`git init` 见 curation_log §4 待办）。

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
| s01_report_mouse | marker 校验(spec §9 清单) | 期望组=argmax(或≥50%argmax) → FAIL=0(清单内) / 1(清单外), ABSENT=2 |  |

- [x] 人工确认: WARNING 已逐条看过且非阻断 —— XYY 2026-11-24(6 条逐条判读, 均不触及鼠侧统计路径)

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
| s01_report_mouse | 形状锚点 | (115282, 20409) → (115282, 20409) | PASS |
| s01_report_mouse | 9 组细胞数锚点 §4.1 | 精确匹配 → 匹配 | PASS |
| s01_report_mouse | 趋势锚点 ITPR1×Purkinje | rho + → rho=+0.915(n=10 bin, 粗算) | PASS |
| s01_report_mouse | 趋势锚点 MBP×OligoLineage | rho + → rho=+0.943(n=6 bin, 粗算) | PASS |
| s01_report_mouse | 趋势锚点 ATOH1×GranuleLineage | rho - → rho=-0.907(n=12 bin, 粗算) | PASS |
| s02_report_mouse | 趋势锚点 ITPR1×Purkinje(正式) | rho + | rho_bin=+0.915 q_bin=0.000204 status=ok | PASS |
| s02_report_mouse | 趋势锚点 MBP×OligoLineage(正式) | rho + | rho_bin=+0.943 q_bin=0.0167 status=ok | PASS |
| s02_report_mouse | 趋势锚点 ATOH1×GranuleLineage(正式) | rho - | rho_bin=-0.907 q_bin=4.63e-05 status=ok | PASS |

- [x] 人工确认: 锚点数值与 spec §4 逐项核对 —— XYY 2026-11-24(19 条 FAIL=0; 唯一 WARN 为 castonguay 中位数, 非阻断)

## (3) marker 校验通过
- 自动: ⚠ 有 WARN, 待人工裁决(FAIL 0)(marker 行 1 条; FAIL→不过, WARN→人工裁决)
| 来源 | check | observed | level |
|---|---|---|---|
| s01_report_mouse | marker 校验(spec §9 清单) | 期望组=argmax(或≥50%argmax) → FAIL=0(清单内) / 1(清单外), ABSENT=2 | WARN |

- 注: 鼠免疫 marker 按 dropout 豁免(spec §11.1: Neurod1/Ptprc/Aif1/Csf1r/Tmem119; 免疫锚点改用 Cx3cr1/P2ry12/Hexb)。
- [x] 人工确认: ITPR1/CA8→Purkinje, GABRA6/CBLN3→GranuleLineage, MBP/PLP1→OligoLineage, AQP4→Astro, SOX2/MKI67→Progenitor, PTPRC→Immune(鼠免疫豁免如上)
      —— XYY 2026-11-24。判读注: 清单内 FAIL=0; PAX6(清单外)argmax=GranuleLineage 属真实生物学(菱形唇/颗粒祖细胞), 非标注错误;
      EOMES/LMX1A 为 UBC(descriptive 层)marker, 不入 8 主组校验 → ABSENT 属设计行为; CA8/SOX2 在鼠侧按 species_valid(CL-018)不适用。

## (4) 趋势锚点方向正确
- 自动: ✓ 自动判定通过(正式趋势锚点 3 条; §9(4) 期望: ITPR1×Purkinje ρ>0; MBP×OligoLineage ρ>0; ATOH1×GranuleLineage ρ<0; 方向不符记 WARN 亦判不过)
| check | observed | level |
|---|---|---|
| 趋势锚点 ITPR1×Purkinje(正式) | rho + | rho_bin=+0.915 q_bin=0.000204 status=ok | PASS |
| 趋势锚点 MBP×OligoLineage(正式) | rho + | rho_bin=+0.943 q_bin=0.0167 status=ok | PASS |
| 趋势锚点 ATOH1×GranuleLineage(正式) | rho - | rho_bin=-0.907 q_bin=4.63e-05 status=ok | PASS |

- 锚点来源(config trend_anchors.source):
  - ITPR1×Purkinje(+): spec §9(4); Sepp 2024 鼠轴 PC 成熟程序; 草稿 v0.9.1 实测 rho=+0.915(s01 粗算)
  - MBP×OligoLineage(+): spec §9(4); 髓鞘化上升; 草稿实测 rho=+0.943
  - ATOH1×GranuleLineage(-): spec §9(4); 祖细胞标记分化下调; 草稿实测 rho=-0.907
- [x] 人工确认: 方向与量级合理 —— XYY 2026-11-24(ITPR1 +0.915 / MBP +0.943 / ATOH1 −0.907, 与 draft 逐位一致)

---
reviewer: XYY  日期: 2026-11-24  结论: [x] 过门 / [ ] 整改后复检
(签署方式: 对话中逐条判读后裁决, 由 Claude Code 转录; gene_set v1.0, 产物 results/v1.0/)
