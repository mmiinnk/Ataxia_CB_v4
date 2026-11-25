# 来源核实轮·人工任务说明（决定二·方案 B）

> 生成：2026-10-29。配套文件：`checklist.tsv`（190 行 = 基因×疾病行全量，由
> `scripts/make_verification_checklist.py` 从 `genes/ataxia_gene_disease.tsv` + `ataxia_genes.tsv` 派生）。
> 对应验收标准：**plan §0 #2**（每基因至少一条 GeneReviews/OMIM/PanelApp 来源）；
> 销项目标：`curation_log.md` CL-003/004/005/008/015/021 与冻结清单第 1 条。

## 0. 这项工作在做什么、为什么

当前基因表 179 个基因中 **174 个的 `source_primary=unverified`**：值（起病年龄、遗传方式、机制类、
病程、病理）来自 xlsx 草稿与文档转引，没有落到合格文献出处。这些值会直接进入统计——
`onset_typical_y` 是第 4 层（起病年龄×表达窗口 Spearman）的自变量，`mech_class` 是 Fig3
机制类富集的分组变量——**一行无出处的错误值就是一个潜在的假生物学结论**；同时这张表将来是
论文的 Supplementary Table，审稿人抽查任何一行的成本极低。本任务就是逐行把出处补上、把值复核。

## 1. 你要做的三件事（每行）

1. **找来源**：按优先级（plan §2）GeneReviews 专篇 > OMIM Clinical Synopsis > PanelApp 条目；
   evidence 定级另查 ClinGen gene-disease validity。
2. **摘原文**：把支持该行关键值的原文句子**逐字摘录**到 `evidence_quote_v`（英文原文即可，≤200 字）。
   起病年龄必须摘到含 mean/median/range 的那句话；light 行摘基因–疾病关系确立句即可。
3. **对值 + 裁决**：文献值与预填 `*_cur` 列逐列比对，按 §3 规则填 `*_v` 列；
   `verdict` 填 confirmed / corrected / unverifiable 之一，签名（`curator`）+ 日期（`verification_date`）。

## 2. 文件使用规则（先读，防返工）

- **`*_cur` 与全部预填列是冻结草稿记录，不得修改**——核实结果一律写进右侧空白 `*_v` 列；
  这样"草稿值 vs 核实值"的差异本身就是审计记录。
- 填写方式任选：LibreOffice Calc / Excel / WPS 打开（UTF-8 TSV），或纯文本编辑器；
  **完成后另存为 `checklist_filled.tsv`（同目录），不要覆盖 `checklist.tsv` 空白模板**。
  多人分工可按批次拆分文件（如 `checklist_filled_batchA.tsv`），行结构保持一致。
- 可以分批交：每批交回后我方抽查 20% 行（quote ↔ 来源 URL 对照），通过后并入。
- **一次只处理一行**（plan §6.2 明令禁止批量猜测）；同基因多疾病行（CACNA1A 3 行、ITPR1 3 行、
  KCNC3/SPTBN2/STUB1/GRM1/FGF14/AFG3L2/PNPT1/TREX1 各 2 行）要按**该行疾病**分别找出处，
  不得用同一段引文覆盖多行（综述级 Overview 引文只能作 `source_secondary` 级补充，写在 notes）。

## 3. 各填写列的格式与规则

**填写要求速查（20 个空白列分三档——不是每列都要填）**：

| 档 | 列 | 何时填 |
| --- | --- | --- |
| **每行必填（8）** | `source_type`, `source_primary_v`, `source_url_v`, `accessed_date_v`, `evidence_quote_v`, `verdict`, `curator`, `verification_date` | 所有行（unverifiable 行只填后 4 列 + notes，前 5 列留空） |
| **full 行条件必填（5）** | `onset_min/max/typical/basis_v`：文献给数值就必填（与 `*_cur` 相同也照抄）；确无数值 → 4 列留空 + notes 写「NA 合理」。`evidence_grade_v`：full 行实为必填（现值几乎全是 unverified，CL-012 靠本轮定级） | full 行（175 行） |
| **仅指定行必填（1）** | `ataxia_core_feature_v`：① 全部 light 行（15）；② 9 个指定基因（CL-023 矛盾 7 基因 + ITM2B/SLC16A2）的 full 行——无论改不改都显式填 true/false | 见左 |
| **仅当与 `*_cur` 不同才填（5）** | `inheritance_v`, `course_v`, `mutation_class_v`, `mech_class_v`, `primary_pathology_v` | corrected 时；一致留空 = "已确认与文献一致" |
| **情境必填（1）** | `notes`：unverifiable（写明检索过什么）/ NA-onset 结论 / 月龄换算 / expert_estimate 理由 / 增行建议 | 触发即填 |

最少工作量：**light 行 = 9 列；full 且全部一致的 confirmed 行 ≈ 13 列**；corrected 行按需增加。


| 列 | 规则 |
| --- | --- |
| `source_type` | GeneReviews / OMIM / PanelApp（三选一，主来源；ClinGen 定级另记入 `evidence_grade_v`） |
| `source_primary_v` | 格式见 plan §2 记录规则：`GeneReviews:NBK1757`、`OMIM:164400`、`PanelApp:Hereditary ataxia v3.x (Green)`。GeneReviews 须记 NBK 号；能同时看到最后修订日期更好（记入 notes） |
| `source_url_v` | 完整 URL |
| `accessed_date_v` | ISO 日期（你实际打开页面的日期） |
| `evidence_quote_v` | 原文逐字摘录 ≤200 字；含起病数值的句子必摘；不得转述、不得翻译后冒充原文 |
| `onset_min_y_v` / `onset_max_y_v` / `onset_typical_y_v` | **只要文献给出数值就必须填**（与 `*_cur` 相同也要照抄——这是显式确认）；典型值取原文 mean/median/typical；只给范围 → 填范围中点并把 `onset_basis_v` 记 `range_midpoint`；先天/出生起病填 0（须原文明确 congenital/at birth）；确无数值 → 三列留空 + notes 注明「NA 合理」 |
| `onset_basis_v` | reported_median / reported_mean / range_midpoint（以上三选一优先）；expert_estimate 不建议，若必须用请在 notes 说明理由 |
| `inheritance_v` / `ataxia_core_feature_v` / `course_v` / `mutation_class_v` / `mech_class_v` / `primary_pathology_v` | **仅当文献值与 `*_cur` 不同才填**（受控词表见 plan §3/§4；`ataxia_core_feature_v` 例外——CL-023 矛盾 7 基因与 ITM2B/SLC16A2 无论是否改值都必须显式填 true/false，见 §5） |
| `evidence_grade_v` | ClinGen gene-disease validity：definitive / strong / moderate / limited / disputed；无 ClinGen 记录时可给等价判断并在 notes 说明依据 |
| `verdict` | confirmed（文献与草稿一致）/ corrected（至少一个 `*_v` 与 `*_cur` 不同）/ unverifiable（找不到合格来源——`*_v` 全留空，notes 写明检索过哪些来源） |
| `curator` / `verification_date` | 签名 + ISO 日期 |
| `notes` | 自由文本：分歧记录、检索过程、需要增行/删行的建议（**不要自行增删行**） |

## 4. 硬性规则（与项目宪法对齐，违反即返工）

1. **无来源不填值**（硬约束 10）：每个 `*_v` 数值都必须能被同行 `evidence_quote_v` 支持；
   凭记忆/印象填值是本轮要消灭的错误本身。
2. **NA 政策**（硬约束 5）：找不到数值就保持 NA——绝不填 0、绝不插值、绝不拿相近疾病的值顶替；
   0 只用于原文明确的先天/出生起病。
3. **unverifiable 是合法结果**：标出来比硬凑一个来源更有价值；这些行将留在 unverified，
   冻结前按 plan §1.2 处置（evidence=limited 或出局）。
4. `onset_class` **不要填**——它由 `onset_typical_y_v` 按 plan §4.3 阈值 (1,13,30) 程序派生。
5. 不改 `mech_class_v` 除非文献明确矛盾；CL-014 四个机制类争议（VAMP1/SLC16A2/IFIH1/mito 亚类）
   属双标注轮（决定三），本轮不裁决。

## 5. 批次、深度与工作量

按 `batch` 列顺序做（a→e），批内按行序：

| 批 | 内容 | 行数 | 深度 | 预估 |
| --- | --- | --- | --- | --- |
| **a** | 起病值有变更史（CL-003/004/005：TMEM240、ITPR1、TTBK2、ATXN10、CACNA1G、FAT2、SPTBN2、ANO10、UBA5）+ core NA-onset（CL-008 其余基因） | 22 | full | 3–5 h |
| **b** | 旧清单 21 基因（CL-015） | 21（17 full + 4 light） | 混合 | 2–4 h |
| **c** | 新近发现（CL-021：RAB3A、THAP11、PIK3R5、CLCN2、PNPT1） | 6 | full | ~1 h |
| **d** | 其余 effective core | 88 | full | 10–15 h |
| **e** | 其余 extended（full）+ peripheral（light） | 53（42 full + 11 light） | 混合 | 4–7 h |
| 合计 | | **190** | full 175 / light 15 | **约 20–30 h**（full 5–15 min/行，light 2–5 min/行；可多人分批并行） |

**light 行（15 行，全部 peripheral）只需三件事**：① 一条合格来源证明基因–疾病关系成立
（`source_type/primary_v/url_v/date_v/quote_v` 五列）；② `ataxia_core_feature_v` 显式复核；
③ `verdict` + 签名。onset/mech/course 等不要求（peripheral 不进任何统计，CL-023）。

**full 且 `onset_typical_y_cur=NA` 的 45 行**是本轮最高价值的行（`special_note` 已标）：
s08 连续统计目前因 NA 自动排除这些行，补上数值即直接提升统计功效；确无数值也要给
「NA 合理」的结论。

## 6. 特殊行的处理要点（`special_note` 列已逐行标注）

- **UBA5::SCAR24**（批次 a 第一优先）：plan 种子 0 岁 vs xlsx 6 岁冲突，现为 NA——
  以 GeneReviews SCAR24 专篇为裁决点，摘到数值即终结 CL-005。
- **CL-023 矛盾 7 基因**（PRPS1、MTPAP、GBE1、PLP1、SURF1、RRM2B、PGAP1）：`ataxia_core_feature_v`
  **必须显式填**（无论 true/false）。填 false 的行不需要你改任何表——重跑 s00b 后该基因自动降
  peripheral（成员是标准的推论，CL-023）；填 true 则维持 extended。MTPAP 行同时核病名缩写
  （SPAX4？）与 MIM 号（现为 NA）。
- **ITM2B、SLC16A2**：其 peripheral 降级依赖 `ataxia_core_feature=false` 这个**初判**（CL-013 待核），
  同样必须显式填 `ataxia_core_feature_v`。
- **PNPT1::SCA25**：基因归属本身待核（单家系连锁，CL-007）——若 ClinGen/文献判 disputed，
  `evidence_grade_v=disputed` 并按 plan §1.2 在 notes 给出处置建议。
- **RUBCN、SYT14、WWOX、VAMP1**：报道家系少，evidence 定级从严（CL-015）。
- **5 个已有「plan §8 已核」转引的行**（TMEM240、ITPR1、TTBK2、ATXN10 等）：数值已有据，
  只需回填 NBK/MIM 编号并复核原文与当前值一致。
- **ATXN8OS、MT-\*、TREX1**（light）：见 §5；TREX1 若你认为应增 RVCL 行，写 notes，不要自行加行。

## 7. 交回之后发生什么（我方流程，供了解）

`checklist_filled*.tsv` 交回 → 抽查 20%（quote ↔ 来源 URL 对照）→ 运行回填脚本（**已写，2026-11-03**）：

```bash
# 先干跑（只出报告不写盘；报告 genes/verification/apply_report_<date>.tsv）
python -I scripts/apply_verification_results.py <root> <checklist_filled.tsv> --check-only
# 正式回填（G2 政策必须显式裁决, CL-024）
python -I scripts/apply_verification_results.py <root> <checklist_filled.tsv> \
    --onset-unverified-policy na|keep-flagged
# 然后重跑 s00b 重裁 effective tier（CL-023）
python -I scripts/s00b_derive_gene_table.py <root>   # --strict 下验收 #2 须转 PASS
```

三道门（详见脚本 docstring）：**G1** 行级规范门——机械偏差自动归一并记录（source_primary 尾注、
basis 自由文本、XLR/复合格式、course 空格），归一不了或 verdict 误标的行 reject 退回；
**G2** onset「草稿有数值→不可证」政策门（CL-024）；**G3** s00b 预演门——acf 翻转若使策展 core
基因命中 peripheral 规则即行级 BLOCKED（现：COQ6/OPHN1，CL-025）。回填 confirmed/corrected 的
source_primary/url/date/quote 与 corrected 值，unverifiable 行保持 unverified；写盘前自动备份
`.bak_preapply_<date>`。之后 → curation_log 销项（CL-003/004/005/008/015/021）→ 版本升级。
全部批次完成且 #2 PASS 后，冻结清单只剩双标注（#6）、BioMart 映射与 P2c 检验。

## 8. 常见问题

- **没有 GeneReviews 专篇**：用 OMIM Clinical Synopsis（记 MIM 号）；再没有 → PanelApp 条目
  （记面板名+版本+颜色）；三者都无 → unverifiable。
- **专篇内起病描述不一致**（不同段落/家系）：以 Natural History/Clinical Description 段为主，
  `onset_min/max_y_v` 记录报道的全距，notes 说明差异来源。
- **文献只有月龄**（婴儿期疾病常见）：换算为年（如 8 个月 → 0.67），notes 记录原文月龄表述。
- **值与草稿一致**：onset 四列照抄（显式确认），其余 `*_v` 留空，verdict=confirmed。
- **发现新的疾病行/别名问题**：不自行改表，notes 记录，我方走 s00a/aliases 流程。
