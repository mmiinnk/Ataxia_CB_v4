# curation_log.md — 基因集策展日志

> 版本：**v1.0**（2026-11-23 冻结；权威来源 = `inputs/merged.xlsx` G3/YYB 人工复核层，用户裁决
> 「基因集就按照 merged.xlsx 中的信息进行冻结」。接线 `scripts/s00a_from_merged_xlsx.py` →
> 派生 `scripts/s00b_derive_gene_table.py`。规模：177 基因 / 188 疾病行；effective core 108 /
> extended 54 / peripheral 15（CL-023 派生层），core∪extended=162）。
> **v1.0 是「集合冻结」而非「验收全过」**：#2（来源逐条核实）/#6（双标注）经用户裁决顺延 **v1.1**。
> #3（mech_class 每类 ≥5）在冻结当日于 G3 权威值下 FAIL（残差类 `other` 仅 3：WWOX/CCDC88C/MME），
> **原样留痕、未以放宽阈值或改判数据绕过**；2026-11-24 经用户裁决认定 `other` 为兜底残差类，
> 按 plan §0-3/§4.1 原文对齐 s00b 检查（残差类豁免、真实机制类仍硬门槛）→ FAIL 0，见 **CL-027**。
> COQ8B/PUM1 两基因的行级剔除经用户确认（YYB 黄标 + OMIM 依据），见 **CL-028**。
> CL-024（onset 未核实政策）/CL-025（COQ6/OPHN1 acf 翻转）顺延 v1.1：merged.xlsx
> `core_grade_pending_h` 全 188 行 = 「YYB未覆盖·待G3终判」，故 acf/evidence/tier 一律沿用 `*_cur`，
> AI 核验层 `*_v` 的翻转与随之而来的 CL-023 自动重裁本轮**不采纳**。
> 历史版本：v0.9.1-draft（2026-10-29 CL-023 三层制，effective core∪extended=164；2026-11-03 增补
> 来源核实轮首过件审计 + `apply_verification_results.py`，新增 CL-024/CL-025）；
> v0.9-draft（2026-10-28 LLM 整理 + 待人工 review）。
> 本文件是 `ataxia_gene_disease.tsv` / `ataxia_genes.tsv` / `gene_aliases.tsv` / `sanity_markers.tsv` 的
> 裁决留痕（方案 §0 产出要求）。`config.yaml gene_set.version` 必须与本版本一致（硬约束 6，s00 版本门）。
> 规则依据：`gene_curation_plan.md`（下称 plan）§1–§6；工程口径：`implementation_spec.md`（下称 spec）。

---

## 1. 输入与产物（可复算）

**输入（v1.0 起）**
0. **`inputs/merged.xlsx` `main` 表（188 行 × 73 列；sha256 见 `inputs/sha256.tsv`）——v1.0 唯一权威来源。**
   四个层：策展层 0–10（`hgnc_symbol`/`tier_curated`/`is_primary_disease`…）、机器草稿层 `*_cur` 11–25、
   **AI 核验层 26–46**（`special_note`/`source_type`/`*_v`/`verdict`/`curator`/`notes`；
   curator 全为 `claude-code(AI-assisted · 待人工复核)`，人工 20% 抽查未做（audit §0）→ **不作权威值层**，
   仅用于 merged.xlsx 别无来源的溯源字段）、**G3/YYB 人工复核层 `*_h` 47–72（权威）**。
   取值一律**按列名**（非位置）由 `scripts/s00a_from_merged_xlsx.py` 完成，运行时打印「字段 ← 列: 行数」
   溯源表；列名缺失/重复即 AssertionError 停（防工作簿改版静默错值）。
1. ~~`共济失调基因分类表(1).xlsx` v2026-10-17（179 基因草稿；快照 `genes/source/gene_classification_table_v20261017.tsv`，
   解析脚本 `scripts/parse_xlsx_snapshot.py`）~~ —— **v1.0 起不再是权威来源**（其内容已经 `*_cur` 层并入
   merged.xlsx）；快照与 `s00a_build_gene_disease_table.py` 保留供 draft 复算/追溯，draft→v1.0 差异见 §5。
2. plan §8 种子清单（含 2026-10-14「已核」条目——本日志中标 `via plan §8 已核` 的数值均转引自该次核实，
   原始出处以 plan 记录为准，v1.0 冻结前仍需回填 NBK/MIM 编号）。
3. 用户旧清单《共济失调致病基因11.6.xlsx》——**仅经 xlsx sheet3 间接获知，非合格来源**（plan §0.2），
   涉及 21 基因见 CL-015。
4. 四个 h5ad 的 `var`（只读实测，2026-10-28）：`/mnt/data/xyy/Ataxia_CB/{human,mouse,aldinger2021_fetal_snRNAseq,castonguay2026_adult_eQTL}.h5ad`。

**产物与生成命令**（幂等，重跑覆盖）
```bash
P=/home/victor/.conda/envs/omics/bin/python          # 单核分析依赖
R=/home/victor/.conda/envs/pd_tremor_cl/bin/python   # 含 openpyxl(读 merged.xlsx)
$P -I scripts/extract_ensembl_mapping.py  $PWD   # → genes/ensembl_mapping_measured.tsv, genes/gene_aliases.tsv
$R -I scripts/s00a_from_merged_xlsx.py    $PWD   # v1.0 冻结接线: merged.xlsx → genes/ataxia_gene_disease.tsv
                                                 #   (首次运行把旧表备份为 *.bak_premerge_<date>, 不覆盖)
$P -I scripts/s00b_derive_gene_table.py   $PWD --onset-thresholds 1,13,30 --span-threshold 20
                                           # → genes/ataxia_genes.tsv + sanity_markers.tsv 回填映射 + 验收自检
$P -I scripts/make_verification_checklist.py $PWD  # → genes/verification/checklist.tsv(190 行人工核实工作表,
                                           #    批次 a–e; 说明书 genes/verification/README.md; 2026-10-29 决定二·方案B)
# draft 链路(留档复算, 非 v1.0): $P -I scripts/s00a_build_gene_disease_table.py $PWD --onset-thresholds 1,13,30
```
注：plan 说 s00b 由疾病表派生基因表——本仓库的 `scripts/s00b_derive_gene_table.py` 即该脚本的引导版。
v1.0 起疾病表的上游是 `s00a_from_merged_xlsx.py`（幂等，可随工作簿更新重跑；重跑会以 merged.xlsx 为准
覆盖 `ataxia_gene_disease.tsv`，故**人工手改应改 merged.xlsx，不要直接改派生 tsv**）。

**规模（v1.0）**：177 基因 / 188 疾病行 / 42 个 sanity marker。策展层 `tier_curated`（主病行计）
core 108 / extended 69；effective tier（CL-023 派生层，s00b 计算）core 108 / extended 54 /
**peripheral 15**（P1 9 / P2a 3 / P2b 1 / P2c 2）；effective core∪extended=**162**（软窗口 120–180 内）。
**vs draft（179 基因 / 190 行 / core 110 / union 164）**：少 **COQ8B、PUM1** 两基因——merged.xlsx 已把
二者列为排除记录（`docs/legacy/README_v3.md`：「188 行对应 177 个基因，另有两条排除记录保存在原工作簿」；
`PUM1::SCA47` 即 audit §0 中 verdict 空缺、等待 Pass3/PanelApp 兜底的那行），core 110→108 全部由此二基因
退出，**无任何基因发生 curated tier 变更**；`COA7::COA7`→`COA7::SCAN3`、`POLG::MIRAS`→`POLG::MIRAS / SANDO`
为同基因同行的缩写更名。字段级差异（186 个共同行）见 §5 v1.0 条。

---

## 2. 决策日志（CL-xxx；s00a/s00b 脚本内引用同编号）

- **CL-001 溯源链约定**（2026-10-28）：每疾病行 `source_secondary` 记录来源链
  （xlsx 快照 + plan §8 种子 / 旧清单 via xlsx / 检索补入 via xlsx）；`source_primary` 仅当有合格出处
  （GeneReviews/OMIM/PanelApp，含 via plan 已核）时填写，否则 `unverified:pending…`。本表 **174/179 基因
  source_primary=unverified** → 验收 #2 未过，为 v1.0 冻结的首要待办（§4 清单）。
- **CL-002 总数上限偏离 →已裁决（2026-10-29）**：core∪extended=179 > plan §0.1 上限 160。根因：xlsx 补入
  旧清单 21 基因 + 检索补入约 30 基因。**用户裁决**：上限 160 系设计期估计、非生物学真理，不作为必须
  满足的硬性条件；采纳三层结构，但 peripheral 成员**不照抄任何名单**（xlsx sheet3 的 11 基因清单与
  表内 `ataxia_core_feature` 实测注释互相矛盾——其中 7 个标注为 true），改按成文标准逐基因派生。
  执行记录见 **CL-023**；plan §0.1 已修订为软窗口 120–180（WARN 语义）。
- **CL-003 恢复 plan「已核」起病值（xlsx 曾无出处改动）**：
  TMEM240 10→**17**（ADCA I 综述均值 17.4，范围 6–30；onset_class 重算 childhood→juvenile）；
  ITPR1/SCA15-16 35→**40**（GeneReviews 均值 39.6，范围 20–66）；
  TTBK2 25→**30**（GeneReviews NBK1757 各家系 9–50，range_midpoint；重算 juvenile→adult）；
  ATXN10 35→**32**（均值 32.2，范围 18–45）。
  以上出处均为 plan §8 已核记录（2026-10-14）；xlsx 的改动无任何引注，按硬约束 10 恢复。
- **CL-004 保留 xlsx 有据修正**：CACNA1G 40→**25**（xlsx 记录文献范围 13–32）；FAT2 40→**60**
  （SCA45 晚发 50–73）；SPTBN2 代表行 onset_class=juvenile(25)；ANO10=30（plan §7 修正方向一致，
  min/max=20/50）。这些修正方向正确，但**原文引注仍未落表**（basis=expert_estimate），冻结前补引用。
- **CL-005 UBA5 冲突置 NA**：plan 种子 0 vs xlsx 6（引"用户表 2 例 5、8 岁"，非合格来源）→
  `onset_typical_y=NA`，类标签暂 childhood，待 GeneReviews SCAR24 专篇裁决。
- **CL-006 多疾病基因拆行**（基因×疾病层，plan §5.1）：CACNA1A（SCA6 代表/EA2/FHM1）、
  ITPR1（SCA15/16 代表/SCA29/Gillespie）、KCNC3（R420H 成人型代表 ~40/R423H 婴儿型 0，取代 xlsx
  文本值 variable）、SPTBN2（SCA5 代表/SCAR14）、STUB1（SCA48 代表/SCAR16）、GRM1（SCAR13 代表/SCA44）、
  FGF14（SCA27B 代表/SCA27A）、AFG3L2（SCA28 代表/SPAX5）。次要行 onset 多为 expert_estimate/NA，
  冻结前逐条核实；FHM1 行 `ataxia_core_feature=false`。
- **CL-007 PNPT1 降级并重做**：xlsx 行自相矛盾（代表疾病 SCA25 为原始家系成人起病，起病类却填
  childhood——后者是 PNPT1 儿童期 MDDS/DEE 表型）。拆为 SCA25（AD，evidence=limited，起病 NA）与
  MDDS（AR，ataxia_core_feature=false）两行；tier core→**extended**。依据 plan §1.2（单家系连锁证据
  不入 core）。SCA25 基因归属本身待 ClinGen/原文核实。
- **CL-008 文本起病值 → NA 政策**：xlsx 中 52 个非数值起病（adult/childhood/variable 等），其中
  core 11 个（CWF19L1、CYP27A1、MTTP、NPTX1、PCNA、PHYH、PUM1、SAMD9L、TDP2、UBA5、ZFHX3）。
  处理：`onset_typical_y=NA`（绝不插值/填 0，硬约束 5 的口径延伸）、`onset_class` 暂沿用 xlsx 标签
  （basis=expert_estimate）、连续统计（log1p Spearman，s08）自动排除 NA，类别统计（Wilcoxon）保留。
  此为对 plan §4.3「onset_class 程序化派生」的**记录在案的偏离**：NA 行无法派生，标签为过渡态。
- **CL-009 in_panel=false 名单**：ATXN8OS（plan §1.2 非编码 RNA，xlsx 备注同）；MT-TK、MT-TI
  （实测四数据集 var 均无此特征，2026-10-28）。三者 tier 不变，但默认统计（core & in_panel）排除。
- **CL-010 mtDNA 基因与新增词表项**：MT-ATP6/MT-ND5/MT-TK/MT-TL1/MT-TI 保留 extended
  （2026-10-29 更新：五者经 CL-023 规则 P2a/P2b/P2c 进入 **peripheral 派生层**，策展层 tier 不变）；
  实测仅 MT-ATP6/MT-ND5/MT-TL1 存在于 Sepp human var，Aldinger/Castonguay 全缺 → 只能做主数据集
  描述性分析，无外部验证。**新增词表项待回写 plan**：`mutation_class=mtDNA_variant`（§4.2）、
  `primary_pathology += cerebellar_atrophy/PC_dysfunction/bergmann_glia/motor_neuron`（§4.4；
  前三者本就在 plan §8 种子中使用，属形式化）、辅助列 `clinical_group_draft`（xlsx 临床大类，
  未注册进 plan §3 schema，暂以 `_draft` 后缀标记，图表使用前须注册进 config）。
- **CL-011 遗传方式规范化**：`AD/AR`、`AR/AD`→`AD+AR`（plan §3 词表；本次拆行后基因层已无复合值，
  仅疾病行历史备注留痕）；`XL(前突变)`→`XL`（FMR1，前突变机制由 mutation_class=repeat_noncoding 表达）；
  `Mito(母系)`→`mito`。
- **CL-012 evidence 初值**：仅按现有文档依据填 limited——TGM6（plan §1.2 争议）、TRPC3/CACNB4/
  CCDC88C（xlsx 备注"证据有限"）、PNPT1-SCA25（CL-007）。**其余 174 基因 evidence=unverified**，
  冻结前按 ClinGen/GeneReviews 逐条定级。
- **CL-013 ataxia_core_feature=false 名单**：SCN2A/SCN8A/KCNJ10/GOSR2/PRRT2（plan §1.2 + xlsx 临床大类
  "伴发癫痫(ataxia次要)"）；PEX1/PEX6（xlsx 备注"ataxia非核心"）；ITM2B/SLC16A2（本次初判：
  CAA-痴呆/AHDS 以其他表现为主，**待核**）；疾病行级：FHM1、PNPT1-MDDS。
  2026-10-29 注：基因级 9 个（前两组）现为 CL-023 规则 P1 的自动降级对象（peripheral）；ITM2B/SLC16A2
  的 false 标注若被决定二核实轮推翻，重跑 s00b 即自动升回 extended，无需新裁决。
- **CL-014 机制类待仲裁清单（未改值，仅挂 note）**：VAMP1（SNARE，与 GOSR2 归
  cytoskeleton_transport_glia 的先例不一致——两者应统一口径）；SLC16A2（MCT8 激素转运，归
  cytoskeleton_transport_glia 牵强，备选 mito_metabolic/other）；IFIH1（胞内 RNA 感受器，other vs
  proteostasis_RNA）；mito_metabolic 过载（core 22/110=20%，含线粒体+溶酶体+脂代谢+糖基化+
  过氧化物酶体，Fig3 解读时注意异质性，或在 v1.1 考虑亚类）。双标注轮（CL-016）一并仲裁。
- **CL-015 旧清单 21 基因**（ATXN8OS、RUBCN、SYT14、WWOX、VAMP1、COQ2、COQ6、COQ8B、PDSS1、PDSS2
  为 core；PEX1、PEX6、SURF1、RRM2B、PLP1、PRPS1、GBE1、ITM2B、PGAP1、MTPAP、SLC16A2 为 extended）：
  来源为"用户旧清单 11.6.xlsx"，非合格 source_primary → 逐条补 GeneReviews/OMIM/PanelApp 出处后才可
  留在 core。其中 RUBCN/SYT14/WWOX/VAMP1 报道家系数少，evidence 定级从严。
- **CL-016 双标注未执行**（验收 #6）：本次为单轮 LLM 整理。冻结前需第二轮独立标注
  （只标 mech_class 与 onset_typical_y，不看本轮结果），一致率 ≥85%，分歧仲裁记录于此。
- **CL-017 Ensembl 映射为实测占位**：`ensg_id/ensmusg_id` 提取自四数据集 `var`
  （ENSG 命中 212/214 query；仅 MT-TK/MT-TI 无人侧 ID；ENSMUSG 201/214）。
  core&in_panel（109 基因）命中率实测：human 1.000 / mouse 0.945（缺 Atn1/Pdyn/Samd9l/Ca8/Kiaa0586/Tctn2）/
  aldinger 0.954（缺 MRE11/TWNK/COQ8A/RUBCN/COQ8B）/ castonguay 1.000，全部过 spec `min_hit_rate` 门禁；
  缺失基因按 NA 政策处理（不插值，s01 全骨架 reindex 承接）。`ortholog_type` 填
  `symbol_match(measured)/none_found(measured)`，**非 BioMart one2one 判定**——联网后按 plan §6.5
  用 BioMart 重导并覆盖。实测重要发现：① Aldinger var 缺失 TWNK/COQ8A/COQ8B/RUBCN/CPLANE1/MRE11
  （ENSG 锚定确认非旧符号问题，是真缺失）；② TREX1 仅 Aldinger 有，GOSR2/MARS2 经 Aldinger/
  Castonguay 补到 ENSG；③ 鼠 MRE11 以旧名 **Mre11a** 存在（已入 gene_aliases.tsv）；
  ④ 鼠 var 缺 Atn1/Pdyn/Samd9l/Ca8/Kiaa0586/Tctn2/Atxn8os/Trex1/Prrt2（作者过滤，跨物种比较时
  这些基因鼠侧按 NA 政策处理，s01 起全骨架 reindex 自然承接）。
- **CL-018 sanity markers**：42 行按 plan §8 转录，组名用 config 正式词（Purkinje/GranuleLineage/DN/
  Interneuron/UBC/OligoLineage/Astro/Progenitor/Immune）。双标（is_also_panel=true）：CA8/ITPR1/FAT2/
  SLC1A3（plan 带*）+ **PLP1/CSF1R**（xlsx extended panel 成员，plan 未标，本次补）。实测修正：
  **CA8、SOX2 在 Sepp mouse var 缺失** → species_valid 改 human（鼠侧校验豁免，与 spec §3.12 免疫
  dropout 豁免同口径）；NEUROD1/PTPRC/CSF1R/AIF1/TMEM119 维持 plan 的 human-only（鼠 dropout 实测值
  见 spec §3.12）。
- **CL-019 polyQ 口径**：xlsx sheet3 称"不满足 §0.3 每类 core ≥8"系**误引**——plan §0.3 阈值为 ≥5
  （≥8 是 onset congenital 类阈值）。实测 core mech_class=polyQ_repeat 共 6（ATXN1/2/3/7、TBP、ATN1），
  **达标**，按 spec `min_genes_per_class=5` 保留、图中标 low n。Fig4e polyQ 集按 plan §4.2 用
  `mutation_class==repeat_polyQ`（contains 语义）：拆行后由"任一疾病行"规则派生
  `is_polyQ_set`，core=7（6+CACNA1A）+extended THAP11=8，矛盾消除。
- **CL-020 引导期阈值口径**：config.yaml 尚未建立，onset 阈值 (1,13,30) 与 span 阈值 20 经 CLI 传入
  （默认=plan §4.3/§5.2 值），敏感性网格 (1,10,25)/(2,16,35) 同 spec `onset_thresholds_grid`。
  config 建立后**必须以 config 为唯一来源**（硬约束 6），脚本默认值仅兜底。
- **CL-021 新近发现基因**：ZFHX3/SCA4（2023 确立，填补 plan 种子缺口，core）、NPTX1/SCA50（MIM
  620158 via xlsx 备注，core）、RAB3A/SCA52（xlsx 称 PanelApp 绿 MIM 621535，extended）、
  THAP11/SCA51（OMIM 未正式分配，extended）、PIK3R5/AOA3、CLCN2（extended）。MIM 号与 PanelApp
  状态均需在线复核后才可去 unverified。
- **CL-022 course 列整体为规则派生**（xlsx sheet3 §一自述：结构性畸形→nonprogressive、EA→episodic、
  单纯小脑型 SCA→slowly_progressive、其余→progressive），**无逐条文献来源** → 整列视为草稿，
  冻结前按 GeneReviews 病程描述核实（slowly/progressive 边界尤其）。
- **CL-023 三层制：peripheral 派生层**（2026-10-29，执行用户裁决，了结 CL-002）：
  **原则**：peripheral 成员**不照抄任何名单**（xlsx sheet3 的 11 基因清单、原草案的 21 基因清单均不
  采用），由下列成文标准逐基因派生；名单是标准的**推论**而非输入。标准实现于 s00b
  `peripheral_rule_of()`，按顺序判定、命中即记入 `peripheral_rule` 列（审计留痕）：
  - **P1 ataxia 非核心表现**：该基因**无任何**疾病行 `ataxia_core_feature=true`（ANY-row 语义，与
    CL-019 polyQ 规则同构）；
  - **P2a 主数据集不可测**：Sepp 人、鼠两主数据集 `var` 实测均无此特征；
  - **P2b 仅单一主数据集**：仅一个主数据集可测，且两个外部验证集（Aldinger/Castonguay）均缺 →
    跨物种与外部验证均不可行；
  - **P2c mtDNA 编码基因**（`MT-` 前缀）：核文库中线粒体转录本的测量效度存疑（ambient/胞质携带）。
    **本批唯一非机械规则**，附行动项（§4 第 8 条）：s01 实测 MT 基因检出率与信号真实性；若信号真实
    且结构化，MT-ATP6/MT-ND5 可升回 extended。
  **派生结果**（2026-10-29 实测复跑，15 基因）：
  - P1（9）：SCN2A、SCN8A、KCNJ10、GOSR2、PRRT2、PEX1、PEX6、ITM2B、SLC16A2；
  - P2a（3）：MT-TK、MT-TI（0/4 实测，CL-009 同源）、**TREX1**（仅 Aldinger 有、两主数据集均缺，
    实测 var 扫描发现——xlsx sheet3 与本日志此前各版名单均未覆盖它，规则化派生优于手工名单的直接
    证据，CL-017②）；
  - P2b（1）：MT-TL1（仅 Sepp human）；
  - P2c（2）：MT-ATP6、MT-ND5。
  effective 规模：core 110 / extended 54 / peripheral 15；**effective core∪extended = 164**，
  落于修订后软窗口 [120,180]（plan §0.1）。core 计数、类规模（#3/#4）、命中率（#5）均不受影响。
  **名单与标准矛盾的处理**：xlsx sheet3 的 11 个"综合征型"中 7 个（PRPS1、MTPAP、GBE1、PLP1、SURF1、
  RRM2B、PGAP1）在表内标注 `ataxia_core_feature=true`——按 P1 **不进 peripheral，保留 extended**。
  其中 MTPAP（SPAX4 痉挛性共济失调；缩写/MIM 待核实轮坐实，omim_phenotype 现为 NA）本身就是命名性
  ataxia，说明 sheet3 名单确与注释列打架，
  以标准（而非名单）为准正确。这 7 个基因 true 标注的真实性交决定二（来源核实轮）按 GeneReviews/OMIM
  裁决：若被推翻为 false，重跑 s00b 即经 P1 自动进入 peripheral，**无需新的人工裁决**。
  **语义与工程口径**：peripheral 是派生层 effective tier（s00b 写入 `ataxia_genes.tsv` 的 `tier` 列；
  策展层 `ataxia_gene_disease.tsv` 的 tier 保持 core/extended，两表 tier 语义自此不同）。peripheral
  永不进富集统计（s07 过滤 `tier==core & in_panel` 不变），也不进敏感性层
  （`sensitivity.include_extended` 仅覆盖 effective core/extended，config 落地时注明）；仅描述性展示。
  **防呆断言**：peripheral 规则命中策展 core 基因 → s00b 直接 AssertionError（防注释矛盾静默降级）。
  **plan 修订**（2026-10-29，已执行）：§0 验收 #1（软窗口 120–180 + peripheral 不计数）、§1.2
  （癫痫 5 基因降级注记）、§3 tier 行、§5.2（peripheral 派生规则行）、§9（敏感性层口径）。

- **CL-024（未决）：onset「草稿有数值→合格来源不可证」政策**（2026-11-03 登记，来源核实轮审计）：
  `checklist_filled.tsv` 首过件（AI 辅助，curator 自标待人工复核）中 **39 个 full 行**的草稿数值经核实
  留空（README 严格口径 → 回填即 NA），其中 **core 主病行 30 个**；另救回原 NA 行 29 个（core 主病 7 个：
  CYP27A1/NPTX1/PHYH/SAMD9L/TDP2/UBA5/ZFHX3，UBA5=6.5 终结 CL-005 冲突）。s08 净影响实测
  （`apply_verification_results.py --check-only`）：core 主病 onset 可用行 **99 → 76**（政策 `na`）
  或 **→ 106**（政策 `keep-flagged`：保留草稿值 + note 标 `onset_draft_unverified` + 敏感性分析排除，
  即此前提议的 A+ 方案）。脚本已把该裁决实现为 **G2 强制门**：不给 `--onset-unverified-policy`
  旗标即 ABORT，不默认。待用户裁决；选 keep-flagged 则 plan §3 须注册审计列与敏感性开关。
- **CL-025（未决）：COQ6/OPHN1 策展 core × ataxia_core_feature=false 矛盾**（2026-11-03，核实轮 G3 拦截）：
  两行被复核者翻转 acf true→false（COQ6::COQ6 引 OMIM:614650 CS，onset 亦改记肾病中位 1.2y；
  OPHN1::OPHN1 引 OMIM:300486），且均为单疾病行基因、策展 tier=core → 回填后将命中 P1，
  触发 s00b 防呆断言（CL-023 设计行为，非 bug）。apply 已行级 BLOCKED（不写 acf 改动）。
  人工裁决二选一：(a) 复核推翻翻转（若文献支持 ataxia 属核心表现）；(b) 接受 false → 新 CL 条目把
  两基因策展 tier core→extended（附出处与理由），重跑 s00b 自动降 peripheral。冻结前必须了结。
  **v1.0 处置（2026-11-23）**：用户裁决 CL-025 顺延 v1.1。依据 = merged.xlsx `core_grade_pending_h`
  **全 188 行**均为「YYB未覆盖·待G3终判」（G3 未对任何行的 core/acf 终判），故本轮 `ataxia_core_feature`
  一律取 `*_cur`，AI 层 `ataxia_core_feature_v` 的翻转不采纳 → COQ6/OPHN1 保持 acf=true ∩ tier=core
  （与 draft 连续），s00b CL-023 防呆断言不触发。v1.1 待 G3 终判后按 (a)/(b) 了结。
- **CL-026 v1.0 冻结：G3 人工层取代 CL-003 期起病值（5 行，硬约束 10 留痕）**（2026-11-23）：
  冻结以 merged.xlsx G3 `*_h` 为权威，其中 4 个基因的起病值与 draft 期 CL-003 恢复值不同，且 G3
  **自带枚举 basis 与引文**（`onset_basis_enum_h` + `evidence_mech_h`/`source_urls_h` 内 PMID/PMC），
  故非「无出处改动」（硬约束 10 不适用回退）：

  | 行 | draft（CL-003，依据 plan §8 已核） | v1.0（G3 `onset_typical_y_h`） | G3 basis | G3 min/max |
  | --- | --- | --- | --- | --- |
  | `ITPR1::SCA15/16` | 40（GR 均值 39.6，范围 20–66） | **35** | reported_mean | 18 / 66 |
  | `TMEM240::SCA21` | 17（ADCA I 综述均值 17.4，6–30） | **18** | reported_median | 1 / 61 |
  | `TTBK2::SCA11` | 30（NBK1757 各家系 9–50，range_midpoint） | **25** | reported_mean | 4 / 70 |
  | `ATXN10::SCA10` | 32（均值 32.2，范围 18–45） | **32.2** | reported_mean | 12 / 83 |
  | `AFG3L2::SCA28` | 20 | **26.5** | （见行内 G3 引文） | — |

  裁决理由：G3/YYB 人工复核（2026-11-19）晚于 CL-003（2026-10-29）且有出处，按用户「基因集按
  merged.xlsx 冻结」采纳 G3 值；CL-004 两项（`CACNA1G`=25、`FAT2`=60）因 G3 未给 typical 值、
  由 `*_cur` 承接而**完整保留**（无回退）。v1.1 双标注轮（#6 / CL-014）应把这 5 行列为优先对照项。
  附带发现：`UBA5::SCAR24` G3 有 min 5 / max 8 但 typical 空 → v1.0 仍 NA（audit §4 曾救回 6.5，
  因 `apply_verification_results.py` 未正式回填而未入表）；全表 `onset_basis` NA 84 行 → 与 CL-024 一并处置。
- **CL-027 `other` 为残差兜底类：豁免验收 #3 的 ≥5 门槛 + 明定下游使用规则**（2026-11-24，用户裁决）：
  v1.0 冻结集下 `mech_class=other` 的 core 基因仅 3 个（WWOX/CCDC88C/MME），触发 s00b「#3 每类 ≥5」FAIL。
  **裁决**：认定 `other` 为残差兜底类（不是一种机制），不适用 ≥5；硬门槛只作用于 7 个真实机制类。
  **依据**：(i) plan §0 验收 3 / §4.1 原文即规定「类内 core <5 → 并入 `other`（允许，须在本日志记录）」，
  `other` 本就是该规则的落脚点，s00b 的硬 FAIL 实现严于 plan；(ii) 工作簿 `分子机制` sheet 定义
  `other`=「无法归入上述类别」、`unknown`=「仅当正常功能确实无权威注释」，二者不可互换——3 个基因的功能
  均有权威注释（CCDC88C: UniProt Q9P219；WWOX: PMID 42128308；MME: `mech_class_alt=mito_metabolic`），
  故正确取值是 `other`，改判 `unknown` 属事实错误；(iii) 7 个真实类实测 7–24 个基因全部达标
  （mito_metabolic 24 / cilia_congenital_dev 19 / channel_synapse 19 / proteostasis_RNA 15 /
  cytoskeleton_transport_glia 11 / dna_repair 10 / polyQ_repeat 7），规则目的（每类有足够基因支撑共性结论）已满足。
  **不采纳的替代方案**：删 `other` 类 / 删这 3 个基因——前者需推翻 G3 带引文的人工判定或为其强开新类，
  后者把 SCAR12、SCA40、SCA43 三个真实致病基因踢出图谱，直接损害指标 1.1 的覆盖面。
  **落地**：`s00b_derive_gene_table.py` #3 改为只检查真实类 + 新增「残差类成员」信息项；plan §0-3 补写
  豁免条款与下游规则；**基因数据一字未改**（`ataxia_genes.tsv` md5 前后一致）。
  **下游使用规则（s07/s08 实现时必须遵守）**：残差类**不进入** s07 机制类×模块富集网格（"other 富集于某模块"
  无生物学解释；其成员改做单基因描述）；s08 分层置换中残差类**自成 stratum 并标 `low n`**（n=3 → 层内仅
  6 种排列，`atx/stats_enrich.py` 的 `_MIN_PERM_OBS=3` 为结构下限），且必须同时报未分层 `p_naive` 对照。
  **影响面**：s01–s06 不使用 `mech_class`（仅 s01 作注释列携带），已产出的 v1.0 鼠侧图谱零影响。
- **CL-028 COQ8B / PUM1 行级剔除确认**（2026-11-24，用户裁决「两个基因都按照黄标进剔除」）：
  merged.xlsx `excluded` 附表（2 行 × 74 列，比 main 多一列 `exclusion_reason`）记录 YYB 黄标剔除，
  两行 `merge_notes` 均含「行级剔除至附表(YYB黄标建议)」：

  | 行 | `review_flag_h` | `exclusion_reason`（原文照录） |
  | --- | --- | --- |
  | `PUM1::SCA47` | 争议基因，建议删除 | unknown争议基因，建议删除（YYB黄标原文）；另见机器notes OMIM-2025 VUS/Cuillerier链：SCA47重分类为VUS、家系归FGF14/SCA27B |
  | `COQ8B::COQ8B` | 主要表现为肾病，建议剔除 | unknown 主要表现为肾病，建议剔除（YYB黄标原文） |

  即：**PUM1** 的基因—疾病因果关系已不被承认（OMIM 2025 将 SCA47 重分类为 VUS，原家系归因 FGF14/SCA27B）；
  **COQ8B** 主症为肾病而非共济失调（与 acf「共济失调是否核心表现」判据同源，参 CL-025）。
  二者原均为 `tier_curated=core`，故 v1.0 的 core 110→108、鼠统计 panel 109→107 **全部**由此而来，
  无任何其他基因发生 tier 变更。两行仍完整保存在工作簿 `excluded` 表（含 G3 值层：PUM1 onset 30–50、
  mech `proteostasis_RNA`；COQ8B onset 0.03–32 / typical 12 `reported_median`、mech `mito_metabolic`），
  若将来 OMIM/PanelApp 恢复其因果关系，可凭该表按新 CL 条目恢复，无需重新策展。
  注：`s00a_from_merged_xlsx.py` 只读 `main` 表，`excluded` 表**不进**派生链（这是剔除生效的机制）。

---

## 3. 与 plan 的差异汇总（冻结时须修订 plan 或回改表）

| # | 差异 | 建议 |
| --- | --- | --- |
| 1 | ~~总数 179 > §0.1 上限 160~~ **已解决**（2026-10-29） | plan §0.1 已修订为软窗口 120–180 + 三层制；effective 164（CL-002→CL-023） |
| 2 | `mutation_class` 新值 `mtDNA_variant` | 回写 §4.2 |
| 3 | `primary_pathology` 扩充 4 值 | 回写 §4.4（种子里本就在用） |
| 4 | 辅助列 `clinical_group_draft` | 若图表要用 → 注册进 §3 schema + config 词表；否则删列 |
| 5 | onset NA 行的 onset_class 为标签非派生 | §4.3 增补 NA 政策一句（CL-008） |
| 6 | `ortholog_type` 为实测占位值 | §6.5 BioMart 步骤完成后覆盖（CL-017） |
| 7 | xlsx sheet3 "polyQ ≥8" 误引 | 无需改 plan；以 CL-019 为准 |

## 4. 冻结清单（v1.0 处置 → v1.1 待办）

**v1.0 处置（2026-11-23；用户裁决「基因集就按照 merged.xlsx 中的信息进行冻结」，
并接受本阶段不强制第一标注论文）**

| 项 | v1.0 状态 | 依据 / v1.1 待办 |
| --- | --- | --- |
| 权威值层 | ✅ 落地：G3 `*_h` > `*_cur`；AI `*_v` 不作权威值层 | `s00a_from_merged_xlsx.py` 按列名取值 + 运行时溯源表 + schema 断言 |
| 基因集范围：COQ8B / PUM1 剔除 | ✅ **用户确认剔除（2026-11-24，CL-028）** | merged.xlsx `excluded` 附表 2 行（`merge_notes`：「行级剔除至附表(YYB黄标建议)」）。core 110→108、鼠统计 panel 109→107 **全部**由此二基因退出，无其他 curated tier 变更 |
| #2 来源核实（172/177 行 `source_primary=unverified`） | ⏸ **顺延 v1.1**（WARN） | 用户裁决。AI 首过件的候选出处已保留可查（`source_primary_v`；`source_url` ← G3 `source_urls_h` 188/188、`evidence_quote` ← G3 `evidence_mech_h` 188/188）；v1.1 做 README §7 人工 20% 抽查后逐条坐实 NBK/MIM |
| #6 双标注一致率 | ⏸ **顺延 v1.1**（WARN，未执行） | 用户裁决；与 CL-014 机制类争议同轮 |
| #3 mech_class core 每类 ≥5 | ✅ **已了结（2026-11-24，CL-027）** → 全 PASS | 冻结当日 FAIL（残差类 `other` 仅 3：WWOX/CCDC88C/MME），**原样留痕未绕过**；用户裁决认定 `other` 为兜底残差类，按 plan §0-3/§4.1 原文（"类内 core <5 → 并入 `other`，允许，须留痕"）对齐 `s00b` 检查：残差类 `other`/`unknown` 豁免 ≥5，**7 个真实机制类仍为硬门槛**（实测 7–24，全达标）。**基因数据一字未改**（`ataxia_genes.tsv` md5 前后一致）→ `results/v1.0` 无需重跑。下游规则已写入 plan §0-3：s07 富集排除残差类、s08 残差类自成 stratum 标 `low n` |
| CL-024 onset 未核实政策 | ⏸ 顺延 v1.1 | `onset_basis` 本轮只取 `onset_basis_enum_h`(62 行)/`onset_basis_cur`(42 行)，84 行 NA；自由文本引文列（`onset_basis_h` 54 / `onset_basis_v` 35）**未**入表 |
| CL-025 COQ6/OPHN1 acf 翻转 | ⏸ 顺延 v1.1 | `core_grade_pending_h` 全 188 行「YYB未覆盖·待G3终判」→ acf 一律取 `*_cur`（COQ6/OPHN1 保持 true∩core，与 draft 连续），s00b CL-023 core 防呆断言不触发 |
| CL-023 中 acf 驱动的自动重裁（audit §2C 其余 12 行） | ⏸ 顺延 v1.1 | 同上：不采纳 AI 层 acf 翻转 → peripheral 名单与 draft 完全一致（P1 9 基因不变） |
| evidence 定级（CL-012） | ⏸ 顺延 v1.1 | `evidence_grade_v` 分布可疑（audit §3：strong 162/limited 28，零 definitive/moderate）→ 本轮 `evidence` 取 `*_cur` |
| BioMart 同源（CL-017） | ⏸ 顺延 v1.1 | `ortholog_type` 仍为实测占位值 |
| course 核实（CL-022） | ⚠ 部分改善：G3 `course_h` 覆盖 188/188 行，词表比 draft 更细（`slowly_progressive`/`static_nonprogressive` 等） | v1.1 回写 plan §4 词表 |
| P2c MT-ATP6/MT-ND5 实测检验 | ⏸ 仍 peripheral（派生层未变） | v1.1 在 s01 阶段给结论 |
| plan §3 差异表 + **v1.0 新发现词表分歧** | ⏸ 顺延 v1.1（**下游有影响，优先**） | `mutation_class`：G3 `mixed`/`missense` vs draft `LoF`/`missense_GoF_DN`（115 行变更）；`inheritance`：G3 `XLR`/`XLD`/`AR/AD` vs draft 归一 `XL`/`AR`（14 行）；`primary_pathology`：G3 `cerebellar_atrophy` vs draft `PC_loss`（86 行）。详见下方「词表权威」注 |
| `git init`（CLAUDE.md 下一步） | ⏸ 未执行 | 工程项 |

**词表权威已转移（v1.0 起）**：merged.xlsx 自带 6 张定义 sheet——`course `、`inheritance `、`onset `、
`突变类型`、`分子机制`、`病理类型`（每个取值都有中文名、定义与判定要点）→ 这 6 张 sheet 为词表权威，
plan §4 须按之修订；s08 / 疾病投影层按 `mutation_class`/`inheritance`/`primary_pathology`/`course` 分组前
**必须先统一词表**（否则同一含义会因两套取值而分裂成两组）。
已发现的 plan↔G3 实质分歧：plan §4.1 明写「**CACNA1A 归 `channel_synapse`**」，G3 人工层判 `polyQ_repeat`
（与 CL-019 的 Fig4e polyQ 集一致，s00b 实测 `is_polyQ_set` core 7 含 CACNA1A）→ **以 G3 为准**，v1.1 修订 plan 该句。

**v1.1 待办（原 v1.0 冻结清单，按优先级）**

1. **来源核实（验收 #2，最大工作量）**：174 个 unverified 基因逐条补 GeneReviews NBK/OMIM MIM/
   PanelApp 版本+颜色，摘录 `evidence_quote`（≤200 字）。优先 core 110；其中最先：
   (a) CL-003/004/005 涉及的 8 个起病值补原文引注；(b) 旧清单 21 基因（CL-015）；
   (c) 检索补入基因（CL-021 + sheet3 §四"本表新增条目"全部）；(d) core 11 个 NA-onset（CL-008）。
   peripheral 15 基因（CL-023）只做**轻核实**（足以证明收录合理，无需精确起病）；其中 7 个
   sheet3 矛盾基因（CL-023）的 `ataxia_core_feature=true` 标注优先复核——结果直接决定其 effective tier。
   **2026-10-29：人工核实清单已生成（方案 B）**——`genes/verification/checklist.tsv`（190 行×47 列，
   批次 a 22 / b 21 / c 6 / d 88 / e 53；full 175 / light 15；full 且 onset=NA 共 45 行）+
   `README.md`（任务说明/格式规则/硬规则/工作量约 20–30h）。
   **2026-11-03：首过件已回收并审计**——`checklist_filled.tsv`（190 行全填，curator 全部
   =claude-code(AI-assisted · 待人工复核)，即 AI 辅助初填、人工复核未做）。审计与干跑详情见
   `genes/verification/audit_report_2026-11-03.md`；`scripts/apply_verification_results.py`
   已写（三道门 G1/G2/G3，原子写盘+备份+逐行报告），干跑：179 staged / 8 rejected /
   2 BLOCKED（COQ6/OPHN1→CL-025）/ 1 incomplete（PUM1）；G2 政策待裁决（CL-024）。
   人工 20% 抽查（README §7）通过后方可正式回填并重跑 s00b。
   注意本机 WebFetch 对 nature/pmc/biorxiv 被阻——GeneReviews(Bookshelf)/OMIM/PanelApp 的可达性
   需先探测；不可达则该步为人工/在线步骤（plan 定位本就如此）。
2. **双标注（验收 #6）**：第二轮独立标 mech_class+onset，比对、仲裁、留痕（CL-016）；
   一并裁决 CL-014 四个机制类争议。
3. **evidence 定级**（CL-012）：ClinGen gene-disease validity 逐条，disputed/limited 者按 §1.2 处置。
4. **BioMart 映射**（CL-017）：ENSG↔ENSMUSG 同源类型，覆盖 ortholog_type；核对鼠侧缺失 9 基因的
   同源存在性（是真无同源还是 var 过滤）。
5. **course 列核实**（CL-022）。
6. **plan 修订**（§3 表剩余六项；第 1 项已于 2026-10-29 完成，CL-023）+ 本日志补记 +
   版本号升 **v1.0** + `config.yaml gene_set.version` 对齐。
7. 前置工程项：`git init`（CLAUDE.md 当前状态所列下一步，尚未执行）。
8. **P2c 实测检验**（CL-023 行动项，s01 阶段）：核文库中 MT-ATP6/MT-ND5 的检出率与信号真实性
   （ambient/胞质携带评估）；若真实且结构化 → 升回 extended 并重跑 s00b。非冻结阻塞项，
   但须在 v1.0 冻结前给出结论（影响 peripheral 名单终态）。

## 5. 变更历史

- **v1.0 补充（2026-11-24，两项用户裁决落地）**：
  **(1) #3 了结（CL-027）**——用户裁决认定 `other` 为兜底残差类。核对 plan 原文发现 §0 验收 3 与 §4.1
  本就规定「类内 core <5 → 该类在统计中并入 `other`（允许，须留痕）」，即 `other` 正是该规则的落脚点，
  **s00b 把它实现成硬 FAIL 严于 plan**。据此对齐：#3 的 ≥5 硬门槛只作用于 7 个真实机制类（实测 7–24，
  全达标），残差类 `other`/`unknown` 豁免并降为信息项（列出成员 + 下游用法）；plan §0-3 同步补写豁免条款与
  下游规则（s07 富集排除残差类；s08 分层置换中残差类自成 stratum 且标 `low n`）。**未删 `other` 类、
  未改判 3 个基因、未删任何基因**——删类要么推翻 G3 带引文的判定、要么把 SCAR12/SCA40/SCA43 三个真实
  致病基因踢出图谱；且 `unknown`（"正常功能确实无权威注释"）对这 3 个基因不成立。影响面核查：
  s02/s03/s06 **完全不使用** `mech_class`（s01 仅作注释列携带）→ 已产出的 v1.0 鼠侧图谱零影响，
  受影响的只有尚未实现的 s07/s08，规则已前置写入 plan 与 CL-027。
  **(2) COQ8B / PUM1 剔除确认（CL-028）**——用户裁决两基因均按 YYB 黄标剔除（依据原文见 CL-028）。
  **验证**：`s00b` 重跑 FAIL 0 / 退出码 0，`ataxia_genes.tsv` md5 与裁决前逐字节一致
  （`53466a74ee226a150ad731f00428ddb9`）——本轮只改检查口径与文档，**基因数据未动**，
  故 `results/v1.0`（gate 31 PASS / 鼠侧 216 图 / 物种报告 / 门 1 清单）全部有效，无需重跑。
- **v1.0（2026-11-23，集合冻结）**：权威来源切换为 `inputs/merged.xlsx`（用户裁决）。新增接线
  `scripts/s00a_from_merged_xlsx.py`（openpyxl 环境；**按列名取值**，启动断言列名唯一且必需列齐全，
  运行时打印逐字段溯源）。取值层级：G3 `*_h`(47–72) > `*_cur`(11–25)；AI `*_v`(26–46) 仅用于
  merged.xlsx 别无来源的溯源字段且 G3 优先——实测溯源：`inheritance`/`mutation_class`/`mech_class`/
  `course`/`primary_pathology`/`disease_name`/`disease_abbrev` ← `*_h` 188/188，`onset_min_y` ← `*_h` 136
  （`*_cur` 1、NA 51）、`onset_max_y` ← `*_h` 117（NA 71）、`onset_typical_y` ← `*_h` 97（`*_cur` 63、NA 28）、
  `mech_class_alt` ← `*_h` 1（`*_cur` 4、NA 183）、`source_url` ← `source_urls_h` 188、
  `evidence_quote` ← `evidence_mech_h` 188、`curator` ← `g3_curator` 188、`accessed/curation_date` ← `g3_date` 188、
  `note` ← `notes` 188（并追加 `review_flag_h`/`onset_derivation_note`/`g3_verdict`/`quote_check_status_h`/
  `evidence_path_h` 溯源片段）；`onset_basis` ← `onset_basis_enum_h` 62 / `onset_basis_cur` 42 / NA 84；
  `onset_class`/`ataxia_core_feature`/`evidence`/`source_primary`/`source_secondary` ← `*_cur`（理由见 §4 处置表）。
  `clinical_group_draft` 不在工作簿 → 自冻结前旧表备份携带（carry 源固定用 `*.bak_premerge_<date>`，
  避免把本轮派生值当旧值回灌）。
  **事故留痕（硬约束 10）**：首版接线凭记忆硬编码列位置，误把 AI 核验层(26–46)当作「旧 `*_h` 层」，
  造成字段错位：`source_primary` ← `course_v`（写出 "slowly_progressive"/"progressive"）、`note` ←
  `mutation_class_v`（"missense"）、`onset_class` ← `onset_typical_y_v`（数值混入，验收 #4 出现
  `'12.5','36','6.1','6.5','38.5'` 等伪类别）、`onset_basis` 二级回退 ← `onset_max_y_v`；同时误把
  `ataxia_core_feature_v`(37) 当人工层而触发 CL-025 断言。该轮产物已作废并重写（`ataxia_gene_disease.tsv`
  被覆盖；冻结前旧表备份 `*.bak_premerge_2026-11-23` 未受影响，`ataxia_genes.tsv` 因断言先于写盘而未被污染）。
  修正 = 改按列名定位 + 启动断言 + 溯源打印。教训：**工作簿列序/列名不得凭记忆硬编码，必须现场读取并断言**。
  规模与差异（vs draft，186 个共同疾病行）：177 基因 / 188 行；curated core 108（−COQ8B、−PUM1）/ extended 69；
  effective core 108 / extended 54 / peripheral 15，union 162。G3 人工层带来的字段变更行数：
  `onset_min_y` 134、`mutation_class` 115、`onset_max_y` 113、`onset_typical_y` 94（大量 NA→有值，
  起病覆盖显著改善）、`course` 90、`primary_pathology` 86、`onset_basis` 59、`mech_class` 25、
  `inheritance` 14、`disease_name` 3、`mech_class_alt` 1；**无基因发生 curated tier 变更**。
  `mech_class` 25 行变更中 8 个 core 基因由 `other` 判入实类 → 残留 `other` core=3 → 当日验收 **#3 FAIL**
  （原样留痕、未绕过；2026-11-24 经用户裁决按 plan §0-3/§4.1 对齐为残差类豁免 → PASS，见 **CL-027**）。
  验收自检（`s00b`，strict=False）：#1a core=108 PASS；#1b union=162∈[120,180] PASS；#3 **FAIL**（other=3）；
  #4 core onset_class congenital 31 / childhood 26 / juvenile 17 / adult 34 PASS；#5a human 1.000 /
  #5b mouse 0.944 / #5c aldinger 0.963、castonguay 1.000 PASS；#7 PASS；#2 WARN（172/177 行 unverified）；
  #6 WARN（未执行）；NA-onset core 8（CWF19L1/CYP27A1/MTTP/PCNA/SAMD9L/TDP2/UBA5/ZFHX3）；
  is_polyQ_set core 7（ATN1/ATXN1/ATXN2/ATXN3/ATXN7/CACNA1A/TBP）。**当日 s00b 退出码 1（仅因 #3）**；
  `run_all.sh gate/run_mouse/run_human` 不调用 s00b，故不阻断流水线。#3 已于 2026-11-24 了结（CL-027）
  → s00b 退出码 0 / FAIL 0，且 `ataxia_genes.tsv` md5 前后一致（`53466a74ee226a150ad731f00428ddb9`），
  **基因数据未变，results/v1.0 全部产物无需重跑**。
  起病值层「G3 人工层取代 CL-003 期已核值」的 5 行逐条留痕见 **CL-026**（含 CL-004 两项完整保留的说明）。
  `config.yaml gene_set.version` 同步升 `1.0`（`draft=False` → 产物目录 `results/v1.0`，不再打 DRAFT 标；
  s00 版本门要求与本文件头部版本一致）。
  **冻结后全链路重跑（同日）**：`run_all.sh gate` → PASS 31 / WARN 5 / FAIL 0（131 s）；
  `run_all.sh run_mouse` → s01–s06 完成（统计 panel **107** 基因、cube 命中 101、缺失 6；216 个图形文件；
  物种报告 134 行；门 1 自动项 FAIL 0 / WARN 6 → `reviews/gate1_mouse.refresh.md`，人工勾选后方可过门）。
  三个趋势锚点与 draft **完全一致**（ITPR1×Purkinje ρ=+0.915 q=2.04e-4、MBP×OligoLineage ρ=+0.943
  q=0.0167、ATOH1×GranuleLineage ρ=−0.907 q=4.63e-5），确认换表未扰动统计层；panel 109→107
  即 COQ8B/PUM1 退出所致。
- **v0.9.1-draft 增补（2026-11-03）**：来源核实轮首过件审计 + 回填脚本落地。结构审计：190/190 行、
  27 预填列零篡改、verdict 161 corrected / 28 confirmed / 1 空（PUM1::SCA47）/ 0 unverifiable；
  source_type OMIM 117 / GeneReviews 72 / PanelApp 0；evidence_grade strong 162 / limited 28
  （零 definitive——定级是否实查 ClinGen 列入抽查重点）。`scripts/apply_verification_results.py`
  写入：G1 行级规范门（词表/格式/verdict 一致性；机械归一 source_primary 尾注 117 处、basis 自由
  文本、XLR→XL、AD/AR 复合、course 空格；词表外 mech_class_v 13 处丢弃转 CL-014）、G2 onset 政策门
  （CL-024，39 行强制显式裁决）、G3 s00b 预演门（CL-025，COQ6/OPHN1 行级 BLOCKED）。干跑处置：
  179 staged / 8 rejected（6 confirmed 误标 + SACS/SETX basis 歧义）/ 2 BLOCKED / 1 incomplete；
  回填后 effective 构成预测 core 110 / extended 51 / peripheral 18（终态取决于 rejected 修复与
  CL-024/025 裁决）。产物：`genes/verification/{audit_report_2026-11-03.md, apply_report_2026-11-03.tsv}`。
- **v0.9.1-draft（2026-10-29）**：CL-023 三层制执行（用户裁决 CL-002）——peripheral 派生层 15 基因
  （P1 9 / P2a 3 含 TREX1 / P2b 1 / P2c 2），effective core∪extended=164；plan §0/§1.2/§3/§5.2/§9
  同步修订；s00b 新增 `peripheral_rule_of()`、`peripheral_rule` 列、`--union-window` 与 core 防呆断言。
  验收自检：#1a/**#1b（软窗口 120–180）**/#3/#4/#5/#7 全 PASS，#2/#6 仍为冻结门禁（WARN）。
  **可靠性复核（同日）**：全链路重跑（含四 h5ad var 重读）五产物字节级一致；独立重算脚本
  （/tmp/verify_tiers/verify_tiers.py，11 项检查，不 import s00b）全 PASS；P2a/P2b/P2c 六基因存在性
  经全新 h5py 探针直证（TREX1=ENSG00000213689 仅 Aldinger，human/mouse var 仅有旁系同源 MTREX/TREX2
  非别名；鼠 var 的 mt- 基因无 mt-Tl1；MT-TK/MT-TI 四库全无）。决定二·方案 B 人工核实清单生成
  （`genes/verification/checklist.tsv` + `README.md`，`scripts/make_verification_checklist.py`）。
  更正：本日早先 CL-023 草案与对话中曾误记 MTPAP 为"SPAX6"——表内缩写为 SPAX4（xlsx 转引，
  未核），已回改；MIM 号入核实清单批次 b 坐实（硬约束 10 留痕）。
- **v0.9-draft（2026-10-28）**：由 xlsx v2026-10-17 + plan §8 + 四 h5ad var 实测整理生成
  （CL-001–CL-022）；验收自检：#1a/#3/#4/#5/#7 PASS，#1b WARN（CL-002），#2/#6 未过（冻结门禁）。
