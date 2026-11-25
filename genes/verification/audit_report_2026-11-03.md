# checklist_filled.tsv 首过件审计报告（2026-11-03）

> 对象：`/home/victor/Projects/Ataxia_CB/checklist_filled.tsv`（项目根目录，927 KB，mtime 2026-11-03 19:11）。
> 基准：`genes/verification/checklist.tsv`（空白模板）+ `README.md`（填写规范）+ `gene_curation_plan.md` §2/§3。
> 复现命令：
> `python -I scripts/apply_verification_results.py <root> checklist_filled.tsv --check-only`
> （逐行处置报告：`genes/verification/apply_report_2026-11-03.tsv`；早期探索性审计脚本 /tmp/audit_checklist/audit.py，结论已被前者取代）。

## 0. 总体判定

**结构完整，规范可修复，内容待人工复核——尚不能直接回填。**

- 填写者为 AI 辅助会话：190/190 行 `curator = claude-code(AI-assisted · 待人工复核)`，
  `verification_date` 2026-10-29→2026-11-03。按 README §7，**人工 20% 抽查（quote ↔ URL 对照）
  未做之前，本件不构成验收 #2 的合格来源轮**。
- 行仍在补：PUM1::SCA47 verdict 空（notes 自述等待 Pass3/PanelApp 兜底）；多行 notes 含
  Pass2/Pass3 过程痕迹。本件按“进行中快照”审计。

## 1. 结构完整性（全部通过）

| 项 | 结果 |
| --- | --- |
| 行数/行集合 | 190/190，与模板一致，无缺/多/重复 |
| 27 个预填列（*_cur 等） | **零篡改** |
| 未动过的行 | 0 |
| verdict | corrected 161 / confirmed 28 / 空 1（PUM1）/ unverifiable 0 |
| source_type | OMIM 117 / GeneReviews 72 / PanelApp 0 |
| 引文长度 | 全部 ≤200 字 |
| light 行 + SPECIAL9 的 ataxia_core_feature_v | 全部已显式填写 |

## 2. 规范偏差分级（G1 门处置）

**A. 机械可归一（apply 自动处理并记录，无需退回）——约 180 处：**
- `source_primary_v` 带尾注（如 “OMIM:614650 Clinical Synopsis”）117 处 → 归一为 `OMIM:614650`，
  节段信息移入表 note；
- `onset_basis_v` 填自由文本而非四值词表 ~40 处 → 含唯一关键词者映射（midpoint/range→range_midpoint，
  mean/average→reported_mean，median→reported_median）；**congenital/at-birth 措辞支持的 onset=0 行**
  词表本无对应项（规范缺口），归一为 `expert_estimate` 并逐行记录；
- `inheritance_v`：XLR→XL、“AR/AD”→AD+AR、“AR; AD for specific variants(…)”→取首段+限定语移 note；
- `course_v` 空格→下划线；
- `mech_class_v = gain_of_function`（词表外，属 mutation 机制语义）13 处 → **丢弃不回填**
  （README §4 规则 5：本轮不裁 mech_class），建议值记 note 转 CL-014 双标注轮。

**B. 退回标注者修正（8 行，apply 干跑 rejected）：**
- confirmed 误标 6 行（值有语义差异应标 corrected）：ITPR1::SCA29、THAP11::SCA51、GRM1::SCA44、
  PTF1A::PTF1A、MT-ATP6::NARP（primary_pathology/min-max 被改）、ATP7B::WD（onset 15→36.5 + acf 翻转）；
- basis 歧义 2 行（自由文本同时命中 mean+range，无法机械裁决）：SACS::ARSACS、SETX::AOA2。

**C. 人工裁决门（非格式问题）：**
- **G2/CL-024**：39 个 full 行「草稿有数值→核实留空」（core 主病 30 行）。s08 可用数 99 → 76
  （policy=na）或 → 106（policy=keep-flagged）。政策未裁决前正式回填 ABORT。
- **G3/CL-025**：COQ6::COQ6、OPHN1::OPHN1 acf true→false × 策展 tier=core × 单疾病行 →
  回填即触发 s00b 防呆断言。已行级 BLOCKED，二选一裁决见 CL-025。
- acf 翻转其余 12 行属 CL-023 自动重裁：GBE1/PGAP1/SURF1/B3GALNT2/CAPN1/CYP7B1（+ATP7B 修复后）
  → P1 降 peripheral；ITM2B/GOSR2/KCNJ10 → 升回 extended；MT-TL1/TREX1 规则标签 P2→P1（成员不变）。
  干跑预测构成：core 110 / extended 51 / peripheral 18（软窗口 [120,180] 内，union=161）。

## 3. 内容级红旗（20% 人工抽查优先名单；格式合法但语义存疑）

| 行 | 问题 | 证据（本件内） |
| --- | --- | --- |
| CA8::CAMRQ3 | onset 0→9.57，但引文是“mean age of **walking** was 9.57 years”——运动里程碑≠起病年龄 | quote 列自证 |
| COQ2::COQ2 | onset 2.0→34.5 = “birth to seventh decade”全距中点；COQ2 以婴儿型多系统病为主，中点 34.5 生物学误导 | quote+notes 自证（min=0 依 birth 措辞） |
| COQ6::COQ6 | onset=1.2y 取自**蛋白尿**中位（肾病表型），非 ataxia 起病；且 acf 翻转触发 G3 | quote 自证 |
| ATP7B::WD | onset 15→36.5 = 全病谱（肝/神经/精神）范围 3–70 中点；notes 自认“非单纯神经起病”；verdict 又误标 confirmed | notes 自证 |
| FAT2::SCA45 | 引文本身支持 “after age 40”，onset 却留空（“无 GR 出处”）——OMIM 是合格来源（plan §2），处理自相矛盾 | quote vs onset_v 空 |
| VAMP1::SPAX1 | evidence_grade=strong，违 special_note“CL-015 报道家系少，定级从严” | special_note vs grade |
| CACNA1G::SCA42 | 25→43.5 = “9–78 全距中点”，反转 CL-004 保留值；引文含“most had onset in mid-adulthood”，规则合规但需人工认可口径 | quote 自证 |
| 全表 evidence_grade | strong 162 / limited 28，**零 definitive/moderate**——分布可疑，是否逐条实查 ClinGen 待抽查（CL-012 不能据此销项） | 干跑统计 |
| 引文复用 | “Onset at birth”（OMIM CS 套语）跨 ATP2B3/SLC16A2/ARMC9/ITPR1 四基因；Joubert 婴儿期句跨 TMEM67/ARL13B——各有独立 URL/MIM，可接受但属弱引文 | quote 比对 |

## 4. 正向结论（本轮已产生的价值）

- 45 个 NA-onset full 行救回 29 个（含 core 主病 7 个）；**UBA5::SCAR24 = 6.5（5–8，Duan 2016，
  OMIM:617133）终结 CL-005 冲突**；ZFHX3/NPTX1/SAMD9L/TDP2 等 2021+ 新基因首次落值。
- 66 行 onset 数值获得带引文修正；CL-003 四基因（TMEM240/ITPR1/TTBK2/ATXN10）复核为
  16.5/38/29.5/30（与 plan §8 已核值有小偏移，抽查时优先对照）。
- notes 普遍记录了检索过程（888 章 GeneReviews 本地索引、eutils 旁证、Pass2/Pass3），
  “无专篇”结论可追溯——unverifiable=0 但检索留痕充分。

## 5. 到「直接运行得到结果」还差什么（按序）

1. 标注方补完 PUM1::SCA47（verdict + 来源五列或 unverifiable）与 8 个 rejected 行（§2B）；
2. 用户裁决 CL-024（onset 政策）与 CL-025（COQ6/OPHN1）；
3. 人工 20% 抽查（优先 §3 红旗 9 项 + evidence_grade 定级依据）；
4. `apply_verification_results.py`（正式模式）→ `s00b_derive_gene_table.py` 重跑 → CL 销项 → 版本升级。
