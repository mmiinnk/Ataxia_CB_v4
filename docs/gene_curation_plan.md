# Ataxia 基因集与临床注释策展方案（实施版）

> 读者：负责执行策展的人员或 LLM 代理。本文件定义**做什么、按什么规则做、产出什么、怎样算合格**。
> 与 `analysis_protocol.md`（导师汇报版）和 `implementation_spec.md`（流水线工程规范）配套；本文件产出的
> 表是流水线 s00 的输入，字段名在三份文件中保持一致。
> 日期：2026-10-13。

---

## 0. 目标与验收标准（先读）

**目标**：产出一张可追溯、可复算的 ataxia 致病基因主表，使第 1–6 层分析能按"机制类"、"发病年龄"、"病理靶细胞"三个临床维度对基因分组并做统计检验。

**产出文件**（全部放 `genes/`，UTF-8、制表符分隔）：

| 文件 | 粒度 | 作用 |
| --- | --- | --- |
| `ataxia_gene_disease.tsv` | 一行 = 一个 **基因 × 疾病** 对 | 原始策展层，保留所有证据与来源，不做任何合并 |
| `ataxia_genes.tsv` | 一行 = 一个 **基因** | 由上表按 §5 规则**程序化派生**（脚本 `scripts/s00b_derive_gene_table.py`），是流水线直接读取的表 |
| `gene_aliases.tsv` | 一行 = 一个别名映射 | 符号别名 → HGNC 现行符号 / ENSG / 鼠同源 |
| `curation_log.md` | 自由文本 | 每条裁决、每个争议基因的决定与理由，附日期 |
| `sanity_markers.tsv` | 一行 = 一个标志基因 | 非致病基因，只用于流水线 marker 校验，不进统计 |

**验收标准**（s00 门禁逐条硬检查，不满足即停）：

1. `ataxia_genes.tsv` 中 `tier == core` 的基因 ≥ 100；effective `core ∪ extended` 落于**软窗口 120–180**（原 120–160 为设计期估计、非硬性阈值；2026-10-29 修订，CL-002/CL-023——超出窗口记 WARN 并在 `curation_log.md` 说明，不阻塞）。`peripheral` 派生层（§5.2，CL-023）不计入该总数，且不进任何富集/敏感性统计。（§8 种子清单 core 约 99 个，需从 PanelApp 绿色基因补足；补入基因同样逐条核实。）
2. 每个基因至少一条来源（`source_primary` 非空，且为 GeneReviews / OMIM / PanelApp 三者之一）。
3. `mech_class` 受控词表内**每个真实机制类**在 core 中 ≥ 5 个基因，否则该类在统计中并入 `other`（允许，但必须在 `curation_log.md` 记录）。
   **残差类 `other` / `unknown` 不适用 ≥5 门槛**（2026-11-24 明确，CL-027）：按本条与 §4.1 的定义，`other`
   正是"小类并入"与"无法归入现有类"的落脚点，它不是一种机制，对它要求 ≥5 无科学含义。残差类成员的义务改为：
   (i) 逐个有依据（`evidence_mech_h` / `mech_class_alt` / `note` 引 UniProt、PMID 等）；
   (ii) **下游使用规则**——s07 机制类×模块富集**排除残差类**（其成员改做单基因描述，不产生"残差类富集"结论）；
   s08 分层置换中残差类**自成 stratum 并标 `low n`**（n=3 时层内仅 6 种排列，必须同时看未分层 `p_naive` 对照）。
   `polyQ_repeat` 预计只有 6 个（ATXN1/2/3/7、TBP、ATN1），属核心假设类，保留但在所有统计图中标 `low n`。
4. `onset_class` 四类在 core 中各 ≥ 10 个基因（congenital 类允许 ≥ 8）。
5. 每个基因有 `ensg_id`，且在 `/mnt/data/xyy/Ataxia_CB/human.h5ad` 的 `var.index` 中命中率 ≥ 95%；鼠 `ensmusg_id` 命中率 ≥ 90%（允许少量无同源）。
6. 两名独立标注者（或两次独立 LLM 运行 + 人工仲裁）对 `mech_class` 与 `onset_class` 的一致率 ≥ 85%，分歧条目全部在 `curation_log.md` 留痕。
7. §7 中列出的"已知错误/排除"基因不得出现在 core 中。

---

## 1. 纳入与排除标准

### 1.1 纳入（任一满足即进入 `ataxia_gene_disease.tsv`）

- **A 类（core 候选）**：GeneReviews《Hereditary Ataxia Overview》（含 SCA、AR ataxia、EA、X-linked 小节）列出的基因；或 OMIM 表型系列 *Spinocerebellar ataxia*（PS164400）、*Spinocerebellar ataxia, autosomal recessive*（PS213200）、*Episodic ataxia*（PS160120）、*Cerebellar ataxia, mental retardation, and dysequilibrium syndrome*（CAMRQ, PS224050）、*Joubert syndrome*（PS213300）中有确定基因–表型关系的基因；或 Genomics England PanelApp "Hereditary ataxia"（成人+儿童面板）**绿色**基因。
- **B 类（extended 候选）**：ataxia 为**突出但非核心**表现的综合征：脑白质病性共济失调（POLR3A/B、TUBB4A、AARS2、ABCD1、CSF1R、EIF2B1–5、DARS2、GFAP、LMNB1）、I 型干扰素病（TREX1、RNASEH2A/B/C、SAMHD1、ADAR、IFIH1）、可治疗代谢病（CYP27A1、PHYH、PEX7、GBA2、SLC2A1、NPC1）、PanelApp 橙色基因。

### 1.2 排除（记录于 `curation_log.md` 的"排除清单"）

- 基因–表型关系被 ClinGen 或后续文献判定为 *disputed/refuted*（例：ZNF592–SCAR5；TGM6–SCA35 争议较大 → 标 `evidence=limited`，不入 core）。
- 单家系、无功能验证、无后续独立报道（`evidence=limited`，可留 extended，不入 core）。
- 非编码 RNA 基因（ATXN8OS、SCA8）：保留在表中，`in_panel=false`（h5ad 可能无此特征），不参与统计。
- 以共济失调为**次要**表现的广谱神经发育/癫痫基因（SCN2A、SCN8A、KCNJ10、GOSR2、PRRT2）：放 extended，`ataxia_core_feature=false`。（2026-10-29：此类基因经 §5.2 派生规则 P1 降为 `peripheral` effective tier，仅描述性展示，CL-023。）
- 被误标进示意脚本的非 ataxia 基因（§7）。

---

## 2. 数据来源与优先级

| 优先级 | 来源 | 用途 | 访问方式 |
| --- | --- | --- | --- |
| 1 | **GeneReviews**（NCBI Bookshelf）：Hereditary Ataxia Overview + 各病专篇 | 疾病名、遗传模式、**发病年龄范围与典型值**、病理描述 | 网页逐篇阅读；记录 NBK 编号与最后修订日期 |
| 2 | **OMIM** | 基因–表型确认、MIM 号、表型系列 | 网页/API；记录 phenotype MIM 与 gene MIM |
| 3 | **Genomics England PanelApp** "Hereditary ataxia" / "Cerebellar hypoplasia" / "Joubert" 面板 | 证据等级（绿/橙/红）、遗传模式 | 网页或 API 导出 TSV |
| 4 | **ClinGen** gene–disease validity | 争议基因裁决 | 网页 |
| 5 | 综述：Synofzik & Németh 2018 *Handb Clin Neurol*；Klockgether, Mariotti & Paulson 2019 *Nat Rev Dis Primers*（SCA）；Beaudin et al. 2019 *Cerebellum*（AR ataxia 分类）；Coutelier et al. 2017 *Brain*（基因面板实证） | 机制类分组依据、AR ataxia 系统分类 | PDF |
| 6 | **HGNC**（符号/别名/ENSG）、**Ensembl BioMart**（人–鼠同源 ENSMUSG + MGI 符号） | 别名与跨物种映射 | 导出表 |
| 7 | **MGI / IMPC** 小鼠模型表型 | 第 6 层模型鼠签名候选（可选列 `mouse_model_available`） | 网页 |

**记录规则**：每条 `source_primary` 写成 `GeneReviews:NBK1138` / `OMIM:164400` / `PanelApp:Hereditary ataxia v3.x (Green)` 的形式；`source_url` 填完整 URL；`accessed_date` 填 ISO 日期。LLM 执行时**不得凭记忆填写发病年龄与机制**，每条必须由上述来源页面的原文支持，并把支持原文摘录到 `evidence_quote`（≤ 200 字）。

---

## 3. 表结构：`ataxia_gene_disease.tsv`（基因 × 疾病层）

| 列 | 类型 / 受控词表 | 说明 |
| --- | --- | --- |
| `hgnc_symbol` | 字符串 | HGNC 现行符号（大写） |
| `ensg_id` | ENSG… | 人 Ensembl gene ID（与 h5ad `var.index` 匹配的主键） |
| `mgi_symbol` | 字符串 | 鼠同源符号（首字母大写；无同源留空） |
| `ensmusg_id` | ENSMUSG… | 鼠 Ensembl ID |
| `ortholog_type` | one2one / one2many / none | BioMart 同源关系 |
| `disease_name` | 字符串 | 规范病名（SCA1、ARSACS、AOA2、Joubert syndrome 9 …） |
| `disease_abbrev` | 字符串 | 短名（图中用） |
| `omim_phenotype` | 数字 | 表型 MIM 号 |
| `inheritance` | AD / AR / XL / XLD / mito / AD+AR | 遗传模式 |
| `mutation_class` | repeat_polyQ / repeat_noncoding / LoF / missense_GoF_DN / mixed / unknown | **突变机制**（与 `mech_class` 分开记录，见 §4.2） |
| `mech_class` | 见 §4.1 八类 | 蛋白功能/致病通路类别 |
| `mech_class_alt` | 同上或空 | 次选类别（双归属基因，敏感性分析用） |
| `onset_min_y`, `onset_max_y` | 数值（年） | 文献报道范围（先天/出生时记 0；胎儿期记 0） |
| `onset_typical_y` | 数值（年） | 典型/中位发病年龄（GeneReviews 原文"mean/median/typical"；若只给范围则取范围中点，并在 `onset_basis` 标 `range_midpoint`） |
| `onset_basis` | reported_median / reported_mean / range_midpoint / expert_estimate | 典型值的来源类型 |
| `onset_class` | congenital / childhood / juvenile / adult | 由 `onset_typical_y` 按 §4.3 阈值**程序化**派生，不手填 |
| `course` | nonprogressive / slowly_progressive / progressive / episodic | 病程 |
| `primary_pathology` | 多值，分号分隔，词表见 §4.4 | 主要病理靶细胞/结构 |
| `ataxia_core_feature` | true / false | 共济失调是否核心表现 |
| `evidence` | definitive / strong / moderate / limited / disputed | ClinGen 等级或等价判断 |
| `tier` | core / extended / excluded（策展层；派生层 effective tier 另有 `peripheral`，见 §5.2 / CL-023） | §1 规则 |
| `is_primary_disease` | true / false | 该基因的多个疾病中选哪一个代表基因（§5.1） |
| `source_primary`, `source_secondary`, `source_url`, `accessed_date`, `evidence_quote` | 字符串 | 溯源 |
| `curator`, `curation_date`, `note` | 字符串 | 留痕 |

---

## 4. 受控词表与判定规则

### 4.1 `mech_class`（八类；统计时类内 core 基因 < 5 的并入 `other`）

| 代码 | 含义 | 典型基因 | 判定要点 |
| --- | --- | --- | --- |
| `polyQ_repeat` | 编码区 CAG/polyQ 扩展，蛋白功能为转录/RNA 调控或未归类 | ATXN1、ATXN2、ATXN3、ATXN7、TBP、ATN1 | 只放"蛋白功能本身不属于其他类"的 polyQ 基因；**CACNA1A 归 `channel_synapse`**（`mutation_class=repeat_polyQ` 保留 polyQ 身份），敏感性分析可切换 |
| `channel_synapse` | 离子通道、Ca²⁺ 信号、谷氨酸受体、突触支架、神经元兴奋性 | ITPR1、CACNA1A、CACNA1G、KCNC3、KCND3、KCNA1、CACNB4、GRID2、GRM1、PRKCG、TRPC3、FGF14、SLC1A3、CA8、CASK、ATP2B3、ANO10 | 以蛋白功能为准；SLC1A3（Bergmann 胶质谷氨酸转运体）标 `mech_class_alt=cytoskeleton_transport_glia`；CASK（突触支架，但 MICPCH 为先天小脑发育不良）标 `mech_class_alt=cilia_congenital_dev`；CA8（结合 ITPR1 调控 Ca²⁺）与 ANO10（TMEM16 家族，功能有争议）标 `mech_class_alt=other` |
| `dna_repair` | DNA 单/双链修复、复制、拓扑异构酶相关 | ATM、APTX、SETX、PNKP、XRCC1、TDP1、TDP2、MRE11、PCNA、RFC1 | RFC1（CANVAS，内含子重复）`mutation_class=repeat_noncoding` |
| `mito_metabolic` | 线粒体、脂代谢、维生素、溶酶体、糖基化、过氧化物酶体 | FXN、POLG、TWNK、SPG7、AFG3L2、COQ8A、PMPCA、TTPA、MTTP、ELOVL4、ELOVL5、TPP1、PNPLA6、ABCB7、PMM2、CYP27A1、GBA2、PHYH | |
| `cilia_congenital_dev` | 纤毛病、神经元迁移/reelin 通路、发育转录因子、结构性小脑畸形 | AHI1、CEP290、TMEM67、CC2D2A、INPP5E、ARL13B、RPGRIP1L、NPHP1、KIAA0586、TCTN2、CSPP1、OFD1、VLDLR、RELN、DAB1、WDR81、PTF1A、TTBK2 | DAB1（SCA37 成人起病但 reelin 通路）保留在此类——这是"发育通路基因成人起病"的有意义特例 |
| `cytoskeleton_transport_glia` | 细胞骨架、膜运输、自噬、细胞器动力学、胶质功能 | SYNE1、SACS、SPTBN2、KIF1C、VPS13D、SNX14、SCYL1、FAT2、OPHN1 | SPTBN2 为 PC 骨架蛋白（`mech_class_alt=channel_synapse`） |
| `proteostasis_RNA` | 分子伴侣/泛素系统、翻译、RNA 加工、非编码重复 RNA 毒性、转录因子（非发育） | STUB1、SIL1、RNF216、PUM1、EEF2、CWF19L1、NOP56、BEAN1、ATXN10、PPP2R2B、CAMTA1、FMR1、UBA5 | 非编码重复（NOP56、BEAN1、ATXN10、PPP2R2B、FMR1 前突变）归此类并标 `mutation_class=repeat_noncoding` |
| `other` | 无法归入上述类或功能不明 | TMEM240、PDYN、MME、PLD3、SAMD9L、CCDC88C、ATCAY、TGM6 | |

**双归属裁决顺序**：(1) 蛋白的分子功能 → (2) 公认的疾病机制综述分组（Klockgether 2019 / Beaudin 2019）→ (3) 若仍两可，主类取与小脑病理最直接的一类，次类进 `mech_class_alt`。所有双归属基因在 `curation_log.md` 列出。

### 4.2 `mutation_class`

`repeat_polyQ`（编码 CAG）、`repeat_noncoding`（内含子/UTR 重复：FXN、RFC1、NOP56、BEAN1、ATXN10、PPP2R2B、FGF14-SCA27B、DAB1、FMR1）、`LoF`（无义/移码/缺失/剪接，含 AR 双等位）、`missense_GoF_DN`（显性错义、功能获得或显性负性）、`mixed`、`unknown`。第 4 层 polyQ 特例（Fig 4e）取 `mutation_class == repeat_polyQ` 的所有基因（含 CACNA1A）。

### 4.3 `onset_class` 阈值（程序化派生；阈值写入 config，可做敏感性分析）

| 类 | `onset_typical_y` | 序数 |
| --- | --- | --- |
| `congenital` | < 1（出生时/婴儿期起病，含非进行性先天性共济失调与小脑畸形） | 0 |
| `childhood` | 1 ≤ y < 13 | 1 |
| `juvenile` | 13 ≤ y < 30 | 2 |
| `adult` | ≥ 30 | 3 |

同时保留连续变量 `onset_typical_y`（统计时用 `log1p(onset_typical_y)` 做 Spearman），类别只用于分组展示与 Wilcoxon。敏感性：阈值改为 (1, 10, 25) 与 (2, 16, 35) 各跑一次。

### 4.4 `primary_pathology` 词表

`PC_loss`（Purkinje 变性）、`DN_loss`（齿状核/深部核团）、`granule_loss`、`olivopontocerebellar`、`spinocerebellar_tract`、`DRG_sensory`（背根神经节/感觉神经元）、`white_matter`（髓鞘/脑白质）、`cerebellar_hypoplasia`（发育不良/体积小）、`vermis_hypoplasia_molar_tooth`（Joubert 型）、`cerebellar_agenesis`、`retina`、`peripheral_nerve`、`unknown`。来源为 GeneReviews/OMIM 的病理或影像描述。此列用于第 1 层验证："预测主组 = Purkinje"的基因是否富集 `PC_loss`。

---

## 5. 基因层派生表 `ataxia_genes.tsv`（脚本生成，不手填）

### 5.1 代表疾病选择（`is_primary_disease`）

一个基因有多条疾病行时，代表疾病按以下顺序选：(1) `ataxia_core_feature=true` 且 `evidence ∈ {definitive,strong}`；(2) 患病率/文献报道最多者（GeneReviews 专篇存在者优先）；(3) 仍并列时取发病年龄**较晚**者并在 `note` 记录（理由：避免把罕见婴儿型当代表，使成人起病 SCA 不被低估）。例：CACNA1A → SCA6（代表）、EA2 与 FHM1 作次要行；ITPR1 → SCA15/16（代表），SCA29/Gillespie 次要；STUB1 → SCA48（AD）代表，SCAR16 次要；GRM1 → SCAR13（AR，证据更强）代表；SPTBN2 → SCA5 代表；FGF14 → SCA27B（GAA 重复，患病率高）代表；AFG3L2 → SCA28 代表。

### 5.2 派生列

| 列 | 规则 |
| --- | --- |
| `disease_primary`, `inheritance`, `mutation_class`, `mech_class`, `mech_class_alt`, `onset_typical_y`, `onset_class`, `primary_pathology`, `tier`（策展层初值） | 取代表疾病行；`tier` 的 effective 值另按下行 `peripheral_rule` 规则覆写（CL-023） |
| `onset_min_all_y` / `onset_max_all_y` | 跨该基因所有 core 疾病行的最小/最大 |
| `onset_span_flag` | `onset_max_all_y - onset_min_all_y > 20` 时为 true（多表型基因，Fig 4 中用空心点） |
| `n_diseases` | 该基因疾病行数 |
| `in_panel` | 默认 true；非编码 RNA、`tier=excluded` 为 false |
| `is_polyQ_set` | `mutation_class == repeat_polyQ`（任一疾病行，contains 语义，CL-019） |
| `tier`（effective）+ `peripheral_rule` | **peripheral 派生层（CL-023，2026-10-29）**：按顺序判定——**P1** 无任何疾病行 `ataxia_core_feature=true`（ataxia 非核心表现）；**P2a** 人、鼠两主数据集实测 `var` 均缺；**P2b** 仅一个主数据集可测且两外部验证集均缺；**P2c** mtDNA 编码基因（`MT-` 前缀，附 s01 实测升回条款）。命中 → `tier=peripheral`（仅描述性展示，不进富集与敏感性统计），规则记入 `peripheral_rule`；策展层 tier 不变；命中 core 基因 → s00b 断言失败（防注释矛盾）。成员是标准的推论而非名单：注释或实测存在性变化时重跑 s00b 自动重裁 |

---

## 6. 执行步骤（给策展代理的逐步指令）

1. **建骨架**：从 GeneReviews Hereditary Ataxia Overview 的全部表格、OMIM 五个表型系列、PanelApp 两个面板的绿/橙基因，合并去重得到候选符号列表（预期 180–220 个）。每个符号经 HGNC 规范化（别名 → 现行符号），写入 `gene_aliases.tsv`。
2. **逐条填写**：对每个候选基因，打开 GeneReviews 专篇（无专篇则 OMIM 临床概要），填 §3 全部列；`evidence_quote` 必须是原文摘录。一次只处理一个基因，禁止批量猜测。
3. **机制分类**：按 §4.1 规则填 `mech_class`；对不确定者查 UniProt 功能描述，并在 `note` 写明依据。
4. **双标注**：第二名标注者（或第二次独立运行、不看第一轮结果）只标 `mech_class` 与 `onset_typical_y`；程序比较两轮，分歧条目人工仲裁并记录。
5. **映射**：用 Ensembl BioMart 导出（人 ENSG ↔ 符号 ↔ 鼠 ENSMUSG/MGI 符号 + 同源类型），合并到表；再跑 `scripts/s00_verify_inputs.py` 的基因命中检查，对未命中者查 h5ad 的 `feature_name` 是否使用旧符号（已知：Castonguay 中 KIAA0226 → RUBCNL；SOGA1 在 Castonguay 中无此符号——需按 ENSG 匹配并记录）。
6. **派生与验收**：运行 `scripts/s00b_derive_gene_table.py` 生成 `ataxia_genes.tsv`，运行 §0 验收检查；不达标回到第 2 步补基因或调整 tier。
7. **冻结**：验收通过后在 `curation_log.md` 写"冻结版本 v1.0 + 日期 + 基因数"；后续任何更改升版本号，流水线产物目录带版本号。

---

## 7. 已知需修正的问题（来自示意脚本 `make_sample_figures.py` 的占位面板，不得照抄）

| 示意脚本中的条目 | 问题 | 处理 |
| --- | --- | --- |
| `DNAJC5` 标为 SCA48 | 错误。SCA48 = **STUB1**；DNAJC5 为成人型神经元蜡样脂褐质沉积症（Kufs 病），非 ataxia 核心基因 | 排除 DNAJC5；加入 STUB1 |
| `HPCA` | HPCA 为 DYT2 肌张力障碍基因，非 ataxia | 排除 |
| `ZNF592`（SCAR5） | 基因–表型关系已被质疑 | `evidence=disputed`，排除 |
| `SLC1A6`、`PCP2`、`CALB1`、`GABRA6` | 为细胞标志基因，非 Mendelian ataxia 基因 | 移至 `sanity_markers.tsv` |
| `ANO10` 标"青少年" | SCAR10 典型起病 20–50 岁 | 以 GeneReviews/OMIM 原文为准重填 |
| `TBP` 主组标为少突 | 无依据的示意值 | 由数据决定，不预设 |
| `ATXN3` 标"中间神经元"、`FXN`/`STUB1` 标"深部核团" | 同上 | 同上 |

---

## 8. 种子清单（草稿；每条必须经 §2 来源核实后才可进入正式表）

> 下表只用于加速第 6 节第 1–2 步，`onset_typical_y` 为策展者初步印象值，**不是**来源值；核实时以 GeneReviews 原文覆盖。标注「(已核:...)」的条目已于 2026-10-14 对照 GeneReviews/OMIM/综述核实，可直接填入 `evidence_quote`；其余仍需逐条核实。列：基因 | 代表疾病 | 遗传 | mutation_class | mech_class | 初步典型起病(年) | 主要病理。

**显性 SCA / EA**
ATXN1|SCA1|AD|repeat_polyQ|polyQ_repeat|35|PC_loss;olivopontocerebellar ·
ATXN2|SCA2|AD|repeat_polyQ|polyQ_repeat|35|PC_loss;olivopontocerebellar ·
ATXN3|SCA3|AD|repeat_polyQ|polyQ_repeat|40|DN_loss;spinocerebellar_tract ·
ATXN7|SCA7|AD|repeat_polyQ|polyQ_repeat|30|PC_loss;retina ·
TBP|SCA17|AD|repeat_polyQ|polyQ_repeat|35|PC_loss ·
ATN1|DRPLA|AD|repeat_polyQ|polyQ_repeat|30|DN_loss ·
CACNA1A|SCA6 (代表) / EA2 / FHM1|AD|repeat_polyQ (SCA6) / LoF (EA2)|channel_synapse|45 / 10|PC_loss ·
ITPR1|SCA15/16 (代表) / SCA29 / Gillespie|AD / AD / AR|LoF / missense_GoF_DN|channel_synapse|40 (已核:GeneReviews 均值 39.6,范围 20–66) / 0 / 0|PC_loss ·
KCNC3|SCA13|AD|missense_GoF_DN|channel_synapse|两型分行:R420H 型成人起病(代表,~40) / R423H 型婴儿起病伴智力障碍(次要,0)|PC_loss ·
KCND3|SCA19/22|AD|missense_GoF_DN|channel_synapse|30|PC_loss ·
CACNA1G|SCA42|AD|missense_GoF_DN|channel_synapse|40|cerebellar atrophy ·
PRKCG|SCA14|AD|missense_GoF_DN|channel_synapse|30|PC_loss ·
TRPC3|SCA41|AD|missense_GoF_DN|channel_synapse|40|limited evidence ·
FGF14|SCA27B (代表) / SCA27A|AD|repeat_noncoding / LoF|channel_synapse|60 / 15|PC ·
SPTBN2|SCA5 (代表) / SCAR14|AD / AR|missense_GoF_DN / LoF|cytoskeleton_transport_glia|25 / 0|PC_loss ·
PPP2R2B|SCA12|AD|repeat_noncoding|proteostasis_RNA|35|PC_loss ·
TTBK2|SCA11|AD|LoF (frameshift)|cilia_congenital_dev|30 (已核:GeneReviews NBK1757 各家系 9–50 岁,取范围中点,onset_basis=range_midpoint)|PC_loss ·
TMEM240|SCA21|AD|missense_GoF_DN|other|17 (已核:ADCA I 综述均值 17.4,范围 6–30)|cerebellar atrophy ·
PDYN|SCA23|AD|missense_GoF_DN|other|50|PC_loss ·
EEF2|SCA26|AD|missense_GoF_DN|proteostasis_RNA|40|cerebellar atrophy ·
AFG3L2|SCA28 (代表) / SPAX5|AD / AR|missense_GoF_DN / LoF|mito_metabolic|20 (已核:均值 19.5) / 0|PC_loss ·
BEAN1|SCA31|AD|repeat_noncoding|proteostasis_RNA|58|PC_loss ·
NOP56|SCA36|AD|repeat_noncoding|proteostasis_RNA|50|PC_loss;motor neuron ·
ATXN10|SCA10|AD|repeat_noncoding|proteostasis_RNA|32 (已核:均值 32.2,范围 18–45)|PC_loss ·
DAB1|SCA37|AD|repeat_noncoding|cilia_congenital_dev|48|PC_loss ·
ELOVL4|SCA34|AD|missense_GoF_DN|mito_metabolic|30|cerebellar atrophy ·
ELOVL5|SCA38|AD|missense_GoF_DN|mito_metabolic|40|cerebellar atrophy ·
MME|SCA43|AD|missense_GoF_DN|other|50|peripheral_nerve ·
GRM1|SCAR13 (代表) / SCA44|AR / AD|LoF / missense_GoF_DN|channel_synapse|0 / adult|PC ·
FAT2|SCA45|AD|missense_GoF_DN|cytoskeleton_transport_glia|40|cerebellar atrophy ·
PLD3|SCA46|AD|missense_GoF_DN|other|50|peripheral_nerve ·
PUM1|SCA47 (PRCA)|AD|missense_GoF_DN|proteostasis_RNA|adult|PC_loss ·
STUB1|SCA48 (代表) / SCAR16|AD / AR|missense_GoF_DN / LoF|proteostasis_RNA|42 / 20|PC_loss ·
SAMD9L|SCA49 / ataxia-pancytopenia|AD|missense_GoF_DN|other|childhood|cerebellar atrophy ·
CCDC88C|SCA40|AD|missense_GoF_DN|other|40|limited ·
KCNA1|EA1|AD|missense_GoF_DN|channel_synapse|7|PC dysfunction ·
CACNB4|EA5|AD|missense_GoF_DN|channel_synapse|20|limited ·
SLC1A3|EA6|AD|missense_GoF_DN|channel_synapse (alt: glia)|5|Bergmann glia

**隐性 / X 连锁**
FXN|FRDA|AR|repeat_noncoding|mito_metabolic|12|DRG_sensory;DN_loss;spinocerebellar_tract ·
ATM|AT|AR|LoF|dna_repair|2|PC_loss ·
APTX|AOA1|AR|LoF|dna_repair|5|PC_loss;peripheral_nerve ·
SETX|AOA2|AR|LoF|dna_repair|15|PC_loss;peripheral_nerve ·
PNKP|AOA4|AR|LoF|dna_repair|4|PC_loss ·
XRCC1|AOA5|AR|LoF|dna_repair|30|cerebellar atrophy ·
TDP1|SCAN1|AR|missense_GoF_DN|dna_repair|15|cerebellar atrophy;peripheral_nerve ·
TDP2|SCAR23|AR|LoF|dna_repair|childhood|cerebellar atrophy ·
MRE11|ATLD|AR|LoF|dna_repair|3|PC_loss ·
PCNA|ATLD2|AR|missense_GoF_DN|dna_repair|childhood|cerebellar atrophy ·
RFC1|CANVAS|AR|repeat_noncoding|dna_repair|55|PC_loss;DRG_sensory ·
SACS|ARSACS|AR|LoF|cytoskeleton_transport_glia|2|PC_loss;white_matter ·
SYNE1|SCAR8/ARCA1|AR|LoF|cytoskeleton_transport_glia|30|cerebellar atrophy ·
ANO10|SCAR10|AR|LoF|channel_synapse|30|cerebellar atrophy ·
SPG7|SPG7|AR|LoF|mito_metabolic|40|cerebellar atrophy;spinocerebellar_tract ·
POLG|MIRAS/SANDO|AR|missense_GoF_DN|mito_metabolic|25|cerebellar atrophy;peripheral_nerve ·
TWNK|IOSCA|AR|missense_GoF_DN|mito_metabolic|1|cerebellar atrophy;peripheral_nerve ·
COQ8A|ARCA2|AR|LoF|mito_metabolic|5|cerebellar atrophy ·
PMPCA|SCAR2|AR|missense_GoF_DN|mito_metabolic|0|cerebellar_hypoplasia ·
TTPA|AVED|AR|LoF|mito_metabolic|13|spinocerebellar_tract;DRG_sensory ·
MTTP|abetalipoproteinemia|AR|LoF|mito_metabolic|childhood|spinocerebellar_tract ·
PNPLA6|Boucher-Neuhäuser|AR|LoF|mito_metabolic|15|cerebellar atrophy ·
TPP1|SCAR7|AR|LoF (hypomorphic)|mito_metabolic|10|cerebellar atrophy ·
SIL1|Marinesco-Sjögren|AR|LoF|proteostasis_RNA|1|PC_loss ·
RNF216|Gordon Holmes|AR|LoF|proteostasis_RNA|20|cerebellar atrophy ·
CWF19L1|SCAR17|AR|LoF|proteostasis_RNA|childhood|cerebellar_hypoplasia ·
UBA5|SCAR24|AR|LoF|proteostasis_RNA|0|cerebellar_hypoplasia ·
VPS13D|SCAR4|AR|missense_GoF_DN|cytoskeleton_transport_glia|15|cerebellar atrophy ·
SNX14|SCAR20|AR|LoF|cytoskeleton_transport_glia|0|cerebellar atrophy ·
SCYL1|SCAR21|AR|LoF|cytoskeleton_transport_glia|2|cerebellar atrophy ·
KIF1C|SPAX2|AR|LoF|cytoskeleton_transport_glia|10|white_matter ·
GRID2|SCAR18|AR|LoF|channel_synapse|0|PC;cerebellar atrophy ·
CA8|CAMRQ3|AR|LoF|channel_synapse|0|cerebellar_hypoplasia ·
VLDLR|CAMRQ1|AR|LoF|cilia_congenital_dev|0|cerebellar_hypoplasia ·
WDR81|CAMRQ2|AR|missense_GoF_DN|cilia_congenital_dev|0|cerebellar_hypoplasia ·
ATCAY|Cayman ataxia|AR|LoF|other|0|cerebellar_hypoplasia ·
RELN|LCH (Norman-Roberts)|AR|LoF|cilia_congenital_dev|0|cerebellar_hypoplasia ·
PTF1A|cerebellar agenesis|AR|LoF|cilia_congenital_dev|0|cerebellar_agenesis ·
PMM2|PMM2-CDG|AR|missense_GoF_DN|mito_metabolic|0|cerebellar_hypoplasia ·
CAMTA1|CANPMR|AD|LoF|proteostasis_RNA|0|cerebellar atrophy (nonprogressive) ·
CASK|MICPCH|XL|LoF|channel_synapse|0|cerebellar_hypoplasia ·
OPHN1|XL ID with cerebellar hypoplasia|XL|LoF|cytoskeleton_transport_glia|0|cerebellar_hypoplasia ·
ATP2B3|SCAX1|XL|missense_GoF_DN|channel_synapse|0|cerebellar_hypoplasia ·
ABCB7|XLSA/A|XL|missense_GoF_DN|mito_metabolic|0|cerebellar atrophy ·
SLC9A6|Christianson|XL|LoF|other|0|cerebellar atrophy ·
FMR1|FXTAS|XL (premutation)|repeat_noncoding|proteostasis_RNA|60|PC_loss;white_matter ·
Joubert（全部 congenital, cilia_congenital_dev, LoF, vermis_hypoplasia_molar_tooth）：AHI1、CEP290、TMEM67、CC2D2A、INPP5E、ARL13B、RPGRIP1L、NPHP1、KIAA0586、TCTN2、CSPP1、OFD1(XL)、TMEM216、CPLANE1、MKS1

**extended 候选**：POLR3A、POLR3B、TUBB4A、AARS2、ABCD1、CSF1R、EIF2B5、DARS2、GFAP、LMNB1（白质病）；TREX1、RNASEH2B、SAMHD1、ADAR、IFIH1（AGS）；CYP27A1、PHYH、GBA2、NPC1、SLC2A1（代谢/可治疗）；SCN2A、SCN8A、KCNJ10、GOSR2、PRRT2（癫痫/发作性，ataxia 非核心）；TGM6（争议）。

**sanity_markers.tsv**（非致病，仅校验）：Purkinje = CALB1、PCP2、CA8*、ITPR1*、FOXP2、SKOR2、ESRRB、HOMER3；颗粒 = GABRA6、CBLN3、RBFOX3、NEUROD1（鼠 dropout，仅人用）、FAT2*、GABRD；深部核团 = SLC17A6、LHX9、TBR1；中间神经元 = PVALB、GAD1、SORCS3；UBC = EOMES、LMX1A；少突 = MBP、PLP1、OLIG2、PDGFRA；星形/Bergmann = AQP4、GDF10、SLC1A3*、HOPX；前体 = SOX2、ATOH1、PAX6、NES、MKI67；免疫 = PTPRC、CSF1R、AIF1、TMEM119（**仅人用**；鼠全库 dropout，2026-10-14 在 mouse.h5ad 实测）+ CX3CR1、P2RY12、HEXB（两物种均可用，鼠免疫组检出 42–52%）。表中加列 `species_valid`（human / mouse / both）。（带 * 者同时是 panel 基因，`in_panel=true, is_sanity=true` 双标。）

---

## 9. 与流水线的接口

- s00 读取 `genes/ataxia_genes.tsv`（列：`hgnc_symbol, ensg_id, mgi_symbol, ensmusg_id, tier, peripheral_rule, in_panel, mech_class, mech_class_alt, mutation_class, onset_typical_y, onset_class, primary_pathology, disease_primary, inheritance, is_polyQ_set, onset_span_flag`；2026-10-29 增 `peripheral_rule`，CL-023）；按 `ensg_id` 匹配人 h5ad 与两个外部人类数据集，按 `ensmusg_id` 匹配鼠（失败时回退符号匹配 + `gene_aliases.tsv`）。
- `config.yaml` 的 `gene_set.version` 必须等于 `curation_log.md` 的冻结版本号；不一致 s00 报错。
- 统计默认只用 `tier == core & in_panel`；`extended` 作敏感性重跑（`sensitivity.include_extended: true`）。`peripheral`（派生层，CL-023）在默认与敏感性统计中**均排除**，仅描述性展示；基因表 `tier` 列为 effective 值，规则见 `peripheral_rule` 列。
