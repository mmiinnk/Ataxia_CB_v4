# 小脑发育时间窗与共济失调致病基因 · 研究方案（导师汇报版）

## 拟题（2026-09-22）

- 推荐：*A cross-species single-nucleus atlas of cerebellar development links spatiotemporal expression of ataxia genes to cell-type selectivity and age of onset*
- 短版：*When and where ataxia genes act: developmental time windows across the human and mouse cerebellum*
- 备选（模型价值向）：*Cell-type-resolved developmental windows of cerebellar ataxia genes explain clinical onset age and mouse-model fidelity*

> 数据：Sepp et al., *Nature* 2024 人/鼠小脑发育单核测序图谱（本组已有同数据的成熟分析流水线 Depression_CB，本项目**只参考其已验证的实现思路与数据陷阱，代码独立重写**（见 `implementation_spec.md` §3）；细胞分组、基因集与统计口径均按 ataxia 疾病生物学独立设计，不要求与 Depression_CB 结果可比）。
> 方法学参照：Lee et al., *Nature Genetics* 2026（把"时间"作为主变量、按细胞类型分层建模、用独立队列互验、把外部疾病信息投影到"细胞×时间"分辨率的思路；该文为衰老时间轴，本项目为发育时间轴，互验队列换成人-鼠跨物种）。

---

## 一、一句话研究问题

在已知共济失调（ataxia）致病基因清单上，逐一测量每个基因"**在小脑哪类细胞中表达、在哪个发育时间段高表达**"，检验这张"时空表达地图"能否解释共济失调的两个核心临床现象：**为什么病变选择性地发生在特定细胞（如 Purkinje 细胞）**、**为什么发病年龄从先天跨越到中老年**。

## 二、为什么这组数据合适

- **时间轴完整**：人 16 个发育阶段（胚胎 CS18 → 52 岁），鼠 12 个时间点（E10.5 → 成年）；
- **细胞注释齐全**：Purkinje、颗粒细胞、大胶质、小胶质等可直接提取；
- **人鼠配对**：同一研究、同一注释体系 → 天然回答"小鼠模型对每种共济失调的参考价值有多大"；
- **本组已有先例**：Depression_CB 项目（抑郁基因面板）已用同一数据跑通并人工复核全流程，本项目直接复用，风险低。

### 人类数据的补充数据集（评估结论，2026-09-22）

- **可用：Aldinger et al., Nat Neurosci 2021**（CELLxGENE collection 1b014f39-…，人胎儿小脑 SPLiT-seq snRNA-seq + LCM 空间转录组，含菱唇/脑室区等亚区）→ 定位：**胎儿期独立队列复核**（不同建库平台）+ **空间/亚区维度**（菱唇=颗粒/UBC 起源、脑室区=Purkinje 起源，与先天性共济失调/畸形的发育起源相关）。局限：仅胎儿期，**不能**补产后/成人空洞；SPLiT-seq 灵敏度较低，只做复核与定性，不与主数据合并定量。
- **不用：CNGB CNP0002781**（元数据已核实，43 样本 = 8 例正常胎儿小脑 scRNA-seq（8–17 PCW）+ 35 例儿童髓母细胞瘤组织（scRNA/snRNA/scATAC/bulk，1.5–12.3 岁）；对应 Luo et al., Nature 2022，该文已于 2025-08 被 Nature 撤稿）→ 三重理由排除：① 撤稿数据不作新文章支柱；② 正常小脑仅 8 例胎儿样本，时间区间（8–17 PCW）与主数据 Sepp 人胎期（9–20 wpc）完全重叠，无增量；③ 其余为肿瘤组织，与 ataxia 主题无关。
- **产后/成人空洞的真正补强**（第 6 层验证扩展，链接供人工核查，按优先级）：
  1. Castonguay et al., Nat Genet 2026（PMID 41857318）：成人小脑皮质 snRNA-seq eQTL 图谱，>100 万细胞/109 例 —— https://www.nature.com/articles/s41588-026-02544-8 ；CELLxGENE https://cellxgene.cziscience.com/collections/b3ab2f02-ba90-4abf-87c2-f9d5bd850abc ；Zenodo https://zenodo.org/records/20404733 ；Broad SCP3177 https://singlecell.broadinstitute.org/single_cell/study/SCP3177
     - **平台注记（2026-09-22 核实）**：该数据为 **SPLiT-seq（组合索引）**而非 10x——无微液滴的核尺寸限制，因此成人 Purkinje 可被捕获（原文报告 PC 占比在样本间稳定）。但单核灵敏度低于 10x、dropout 更重 → 定位为**成人 PC 表达模式的秩基/方向性验证**，不与主数据定量合并；与 Aldinger 胎儿 SPLiT-seq 同平台，跨数据集复核在方法学上更干净。
  2. ~~Purkinje 富集 snRNA-seq（SCA7 成人小脑，bioRxiv 2023.03.19.533345）~~：2026-09-22 人工核查，**未提供数据下载，排除**。
  3. Allen WHB-10Xv3 / Human Brain Cell Atlas v1.0（成人全脑含小脑 CB- 解剖块；仅 3 供体；注意 10x 对大核系统欠捕获）：https://alleninstitute.github.io/abc_atlas_access/descriptions/WHB-10Xv3.html ；https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443
  4. bulk 时间序列：BrainSpan https://www.brainspan.org/ ；GTEx https://gtexportal.org/home/
  5. 鼠成年/老化窗：Kozareva 2021 成年鼠小脑 https://www.nature.com/articles/s41586-021-03220-z ；Tabula Muris Senis https://tabula-muris-senis.ds.czbiohub.org/
  注意：Hao et al. Science 2024（Stereo-seq 小脑皮质）物种为猴/狨/鼠、**无人**，不作人类补强。

### 项目数据目录（2026-09-22 已统一到位：/mnt/data/xyy/Ataxia_CB/）

| 文件 | 内容 | 来源 | 入库核验结果 |
| --- | --- | --- | --- |
| human.h5ad | Sepp 2024 人小脑：180,956 核 × 27,410 基因，CS18→52 岁 16 阶段 27 供体 | 拷贝自 /mnt/data/xyy/Nature2024/（原库只读） | md5 与原库一致（f08ae756…） |
| mouse.h5ad | Sepp 2024 鼠小脑：115,282 核 × 20,409 基因，E10.5→9 周 12 时间点 | 同上 | 大小一致；分析起点即此文件 |
| castonguay2026_adult_eQTL.h5ad | **成人**（63–100 岁）人小脑皮质 SPLiT-seq snRNA-seq：1,004,112 核 × 109 供体（正常 750,036 核 + 特发性震颤 254,076 核） | CELLxGENE 731e83e8-…（6.2 GB） | **Purkinje 36,911 核（3.7%）；只用 normal 供体：93 供体 / 31,089 PC 核 / 89 供体 ≥30 PC（中位 257）→ 可做供体级 pseudobulk；16 例特发性震颤供体（5,822 PC）排除**；`X` 为 log 归一化浮点，计数在 `raw/X`；PC 标签经 marker 抽查证实（CA8 2.9 / ITPR1 7.7 / CALB1 3.2，显著高于 Granule/Bergmann/MLI 的 0.1–1.2）；中位 ~1,700 counts/核 → 低灵敏度，**只做秩基/方向验证**；ataxia 面板 43/46 基因可检出（KIAA0226、SOGA1 需别名映射）；另有 Bergmann 97,153、MLI 75,727、Golgi 22,424、UBC 15,561 等成人组 |
| aldinger2021_fetal_snRNAseq.h5ad | 胎儿（9–21 wpc，13 供体）人小脑 SPLiT-seq snRNA-seq：69,174 核，21 个精细标签（H-PC 25,724、H-GCP、H-RL、H-BG、H-MLI 等） | CELLxGENE 99f7a3d3-…（565 MB） | 胎儿期独立复核可用；与 Castonguay 同平台（SPLiT-seq），跨数据集比较方法学干净 |

## 三、分析思路：六层递进

每层按【问题 → 做法（通俗）→ 产出图 → 图怎么读（生物学意义）→ 预期结果】说明。

### 第 1 层　致病基因在哪些细胞中表达（细胞特异性）

- **做法**：细胞按 ataxia 病理分为**九个主组**（Purkinje / 颗粒谱系 / 深部核团 / 皮质中间神经元 / 少突谱系 / 星形胶质(含 Bergmann) / 神经前体 / 小胶质免疫 / UBC，另加 ALL 整体参考；定义与依据见"细胞组设计"）。每个基因算各组平均表达（同一供体、同一阶段、同一类细胞的所有细胞核平均成一点，误差棒=供体间差异，避免"几万细胞当几万样本"的假统计）。再算两个指标：该基因表达量在各组中的**占比**，以及**特异性指数 tau**（0–1，越大=越集中于单一细胞类型）。
- **产出图**：Fig 1——特异性点图（点大小=表达细胞比例、颜色=表达量、附 tau 列）+ 人/鼠分组热图。
- **生物学意义**：给出每个致病基因的"作用位置"。Purkinje 特异的基因，其突变直接破坏 Purkinje 功能 → 对应" Purkinje 变性"型共济失调；广泛表达的基因（如 DNA 修复类）→ 突变可能经多细胞/多系统起病。此层直接对应临床的**病变细胞选择性**。
- **预期**：离子通道/突触类（ITPR1、CACNA1A、CACNB4、SPTBN2、GRID2）集中于 Purkinje 或颗粒细胞；DNA 修复类（ATM、SETX、APTX）与线粒体类（FXN、POLG）广谱；部分隐性基因（SYNE1、SACS）偏胶质。预扫描与 Depression_CB 实测（ITPR1 的 tau≈0.9、Purkinje 特异）已支持。

### 第 2 层　致病基因何时表达（发育时间窗）

- **做法**：每个基因×细胞组，沿发育时间轴画平均表达曲线，做三种描述：①**单调趋势**（随发育升/降，秩相关+统计检验）；②**时间窗型**（仅某段高、前后都低；要求峰两侧都有实测点，避免把"后面没测到"误判成窗口）；③**特异性建立时间**（最早出现"集中于某细胞"的时间点）。
- **产出图**：Fig 2——带供体散点与误差棒的轨迹曲线、时间窗基因清单、特异性建立时间轴。
- **生物学意义**：时间窗=基因参与的发育事件（增殖、迁移、突触形成、髓鞘化、成熟维持）。只在胚胎增殖期高的基因，突变靠"破坏小脑搭建"致病 → 先天/早童起病、常伴结构异常；出生后上升并维持的基因负责成熟细胞功能 → 突变表现为晚发、进行性变性。
- **预期**：方法阳性对照——髓鞘基因 MBP/CNP 应呈出生后（鼠 P14）窗口（Depression_CB 已验证该检测有效）；ataxia 基因中离子通道类出生后上升并维持（鼠 Itpr1 由胚胎近 0 升至成年高表达），DNA 修复/线粒体类全期平直，polyQ 类（ATXN1 等）胚胎即表达、围产期达峰。

### 第 3 层　不同机制类别的致病基因，窗口是否系统性不同（分组富集）

- **做法**：基因集按**机制类**（polyQ、离子通道/突触、DNA 修复、线粒体、胞骨架/胶质等）分组；统计每类基因是否**显著超额**集中于某细胞组的某发育期窗口（比例比较+置换检验，即"比随机抽同样多的基因更集中"的程度）。
- **产出图**：Fig 3——机制类 ×（细胞组 × 发育期）富集热图，星号=统计显著。
- **生物学意义**：把结论从单基因提升到**疾病机制类**：若离子通道类整体落在"出生后 Purkinje 成熟窗"、DNA 修复类整体"胚胎期广谱"，则说明**共济失调的致病路径由其基因的正常发育表达程序预先决定**——这是本文的核心论点。
- **预期**：离子通道/突触类富集于 Purkinje/颗粒的出生后成熟窗；DNA 修复/线粒体类富集于早期广谱或"不特异"；SYNE1、SACS 等富集于星形/少突谱系出生后窗。

### 第 4 层　窗口时序能否解释临床发病年龄（基因水平关联）

- **做法**：给每个基因加临床注释（先天/儿童/青少年/成人起病，来源 GeneReviews/OMIM）；检验"表达峰/窗口出现的时间"在不同起病年龄的基因之间是否系统性不同（分组比较+置换检验）；并画"发育时间轴 × 疾病组"泳道示意。
- **产出图**：Fig 4——左：散点/箱线（x=窗口时序，y=起病年龄类别）；右：疾病组时间窗泳道图。
- **生物学意义**：直接回答"同是小脑病变，为何发病从婴儿期到五十岁"：若成人起病基因的窗口显著晚于先天起病基因，说明**发病年龄很大程度上取决于基因何时开始在该细胞中担当重要角色**——即突变效应"暴露"的时间。此层把图谱从描述提升为对临床表型的解释。
- **预期**：正关联显著（置换 P<0.01）；先天类=胚胎窗/广谱，成人起病类（SCA6、EA2、SCA15 等）=出生后成熟窗。polyQ 类（成人起病但胚胎即表达）作为特例讨论：其 **Purkinje 特异性高表达部分**仍落在成熟期，支持"成年起病=突变毒性需要成熟细胞环境"的经典假说。

### 第 5 层　人鼠是否保守（小鼠模型的参考价值）

- **做法**：用人鼠共享的 33 种"发育状态"（细胞分化状态注释，而非日历年龄）作共同坐标；逐基因比较人鼠轨迹（趋势方向一致性+轨迹相关性）得**保守性评分**；按疾病汇总成"小鼠模型保真度表"。
- **产出图**：Fig 5——人鼠状态对齐并排热图、趋势方向一致性散点、保真度表。
- **生物学意义**：回答实验同事最关心的问题：**哪些共济失调的小鼠模型，其基因表达时空背景与人高度一致**（结论可外推到人）；哪些基因人鼠分歧大（小鼠模型需谨慎、或提示人特有调控）。
- **预期**：成熟窗基因（ITPR1、GRID2、CACNA1A）高保守（Depression_CB 在此数据上实测双物种显著趋势 65 对中 57 对同方向）；少数基因（如 KCNA1 在人核数据低检出、个别胶质基因模式偏移）列入分歧清单。

### 第 6 层　外部独立数据验证

- **做法**：① 与 **BrainSpan** 人脑发育测序（独立样本）比对时间曲线相关性；② 在 **Allen/HDBR 原位杂交图谱**抽查关键基因的时空模式；③ 把已发表的共济失调模型鼠（SCA1、EA2、FRDA 等）转录组差异基因"**回投**"到我们的正常时间窗：看疾病扰动是否集中于我们定义的窗口。
- **产出图**：Fig 6——BrainSpan 相关散点、ISH 示例、模型签名富集条图。
- **生物学意义**：①②证明窗口不是技术假象；③把"正常发育程序"与"疾病状态"连接：若模型鼠紊乱的基因正是成熟窗基因，则支持"**共济失调是发育程序维持/失调病**"，并提示干预的目标时间窗。
- **预期**：BrainSpan 相关 r>0.6；ISH 与预测的细胞类型/时段吻合；模型签名富集于对应窗口。

> 可选扩展层：上游开关——哪些转录因子按时打开这些窗口（共表达/调控网络分析，借鉴 Lee 2026 的调控网络思路）。结果好则作 Fig 7，否则放入展望。

### 细胞组设计（按 ataxia 病理定义，2026-09-22 修订）

本项目为独立的初步分析，**不沿用 Depression_CB 的四组方案、不要求可比**。分组对应文章的三个问题：① 病变累及哪些细胞（小脑回路：颗粒细胞 → Purkinje → 深部核团）；② 各机制类致病基因作用于哪类细胞（通道/突触 → Purkinje/颗粒；髓鞘/脑白质病 → 少突；干扰素病/炎症 → 小胶质）；③ 基因何时起作用（前体增殖迁移窗 → 先天结构性共济失调；成熟维持窗 → 晚发变性）。

**主组（9 组 + ALL 整体参考）**，组归属由 `author_cell_type` + `subtype` 唯一映射：

| 组 | 标签来源 | ataxia 依据 | 可行性（2026-09-22 实测 crosstab，`results_ct_stage/`） |
| --- | --- | --- | --- |
| Purkinje | Purkinje（人亚型 EB/LB；鼠 Foxp1/Cdh9/Etv1/Rorb） | 绝大多数 SCA/EA/FRDA 的原发变性靶细胞；ITPR1、GRID2、SPTBN2、CA8、CACNA1A 直接作用位点 | 人仅产前（CS22=25,993，产后 0–53）；鼠 E12.5→成年全序列 |
| 颗粒谱系 | GC、GC/UBC + 菱唇前体（progenitor_RL(_early)）+ **glutamatergic_uncertain**（已核实为成熟颗粒细胞，见下）；亚型 GC_diff_1(早)/GC_diff_2(晚) 可再拆 | GCP 增殖窗=先天小脑发育不良类（CAMTA1、CASK、PTF1A）作用期；CACNA1A 亦在颗粒表达 | 两物种最大群（人 69k/鼠 31k），全窗可行，统计锚；GU 并入后 44 岁 bin 补至 496 核 → 可行 |
| 深部核团神经元 | glut_DN + GABA_DN（亚型 ventral/posterior 可再拆） | SCA1/2/3/7 齿状核变性（核内包涵体、神经元丢失） | 人仅产前（CS22≈3.8k）；鼠至 P0–P7 后 <30 → NA |
| 皮质中间神经元 | interneuron（GL=高尔基 / ML1、ML2=篮·星状 / PL / early） | GABA 能回路障碍；篮状细胞病理见于多个 SCA | 人产后薄但多数 bin 可行（27–930）；鼠全序列（38–833） |
| 少突谱系 | oligo（OPC/pre_OPC/COP；人分 OPC_early/late） | 脑白质病性共济失调（POLR3A/B、TUBB4A、AARS2、SPG7、ABCD1、CSF1R、RNASEH2B）；FRDA 白质受累；ET 的少突易感性（Castonguay 2026） | 人产后丰富（398–5,264）、产前薄；鼠 P17.5 起 90–158、胚胎 NA |
| 星形胶质（含 Bergmann） | astroglia 中非前体亚型：astro_Bergmann、astro_parenh、astroblast、glioblast_PWM | Bergmann 胶质谷氨酸清除与反应性胶质增生（SCA1/SCA7）；SYNE1 | 人可行（含 Bergmann 1,314）；鼠成熟星形极薄（Bergmann 181/parenh 79）→ 鼠该组由早期亚型主导，成熟期自动 NA |
| 神经前体 | VZ_neuroblast、NTZ_neuroblast、NTZ_mixed + astroglia 内 progenitor_*（VZ/bipotent/gliogenic/isthmic/Nckap5_neg/RP） | 先天/结构性共济失调（Joubert、Dandy-Walker、小脑发育不良）的靶细胞与时间窗 | 两物种胚胎期高丰度（人 CS18–22 达 576–2,848；鼠 E10.5–P0） |
| 小胶质/免疫 | immune（microglia、nonparenh_macrophage、T-cell） | AGS 等 I 型干扰素病伴共济失调（TREX1、RNASEH2B、IFIH1）；神经炎症 | 人可行（45–675，个别婴儿 bin NA）；鼠 6–56/bin 全 NA + 免疫基因文库 dropout → 仅描述 |
| UBC | UBC（亚型 Trpc3/Hcrtr2） | 前庭小脑（绒球小结叶）；部分先天 ataxia 相关 | 人仅胎期+NB（峰 1,920）；鼠产后 <30 → NA；**定位为描述/敏感性分析组** |
| ALL | 全部细胞核 | 参考分母 | — |

**不入组（非小脑或过薄群体，仅保留在 ALL 中）**：isth_N、isthmic_neuroblast、parabrachial、noradrenergic、GABA_MB、MBO（鼠）、motorneuron、ependymal、neural_crest_progenitor、mural/endoth、meningeal、erythroid。

**glutamatergic_uncertain 抽查结论（2026-09-22 完成）**：这 1,107 个核（44 岁供体 #215=496、52 岁供体 #138=611，均来自齿状核组织块 CrDn）**是成熟颗粒细胞，不是漏标的 Purkinje**——① 全转录组 pseudobulk 相关：vs 成人 GC r=0.96–0.98，vs 胎期 PC 仅 0.07–0.18；② GC 特征基因齐备：CBLN3（82%/75% 细胞表达）、GABRA6、RBFOX3、GRID2（98%，均值 13–15）、FAT2、GABRD 均达成人 GC 水平；③ PC 标志全阴：CA8 均值 0.2（胎期 PC 7.4–12.4）、FOXP2 0.1（胎期 PC 29–49）、SKOR2/ESRRB/HOMER3/PCP2≈0；④ GABA 能标志阴性（GAD1/2、SLC32A1≈0），排除中间神经元。→ **归入颗粒谱系组**（配置留"剔除 GU"敏感性开关）。成人 PC 空洞无法由此挽救，改由 Castonguay 数据集承担成人 PC 验证（见数据目录）。

**统计口径**：tau/type_fraction、富集检验在 **8 主组**上计算（UBC 为描述组、不入分母，ALL 为参考）；可行组 <4 的 bin/期 tau、fraction 记 NA，且不同 bin 的 tau/fraction 只在可行组集合相同时可比（输出 `groups_used`）；bin 内 <30 细胞自动 NA、不插值；亚型水平（人 37 / 鼠 34）仅作精细映射补充热图，不做主统计。

### 图面 panel 构成（每 Figure 多 panel，一图一科学信息；逐 panel 详述）

**样例图（2026-09-22）**：`figures_mock/fig{1..6}_sample.png`，由 `scripts/make_sample_figures.py` 生成；
全部为示意数据，仅演示各 Figure 的 panel 组构与画法规则（NA 斜纹、donor 散点、rho 星号、泳道图等），
正式数值待第 1–6 层流水线与外部数据接入后替换。
**示意脚本中与正式口径不一致之处（不沿用）**：9 wpc 被放在胚胎期（正式为胎儿期）；UBC 作为主组（正式为描述组）；
机制类 6 类（正式 8 类）；部分基因与疾病对应错误（见 `gene_curation_plan.md` §7）。

**配套实施文件（2026-10-13）**：`gene_curation_plan.md`（基因集与临床注释策展：纳入/排除标准、受控词表、
表结构、验收标准、种子清单、示意脚本中需修正的错误）与 `implementation_spec.md`（流水线工程规范：完整
config、对 Depression_CB 引擎的逐文件改动、9 组实测锚点与可行性矩阵、s07–s11 新脚本的算法与产物、外部数据
集的 `raw/X` 读取与 normal-only 过滤、执行顺序与 review 门）。后续实现以这两份文件为准。

**读图指南（2026-09-23）**：`figure_reading_guide.md`——按图注尺度给出每张 Fig 的主要问题与逐 panel
"画什么/怎么算/怎么读/caveat"，含通用约定（伪bulk、NA 政策、rho/tau/fraction/emergence、temporal_pattern、
dev_state 对齐、秩基验证）与 S1–S6 简注；下文表格为 panel 组构的工程口径，两者互补。

#### Fig 1｜致病基因在哪里：图谱总览与细胞类型特异性

| panel | 内容与画法 | 怎么读 |
| --- | --- | --- |
| a | 研究设计示意（新绘）：人/鼠双时间轴（人 CS18→52 岁 16 阶段 27 供体；鼠 E10.5→9 周 12 时间点）、小脑回路简图（颗粒→Purkinje→深部核团）、9 细胞组图例、ataxia 机制类色带 | 一页看懂数据规模与设计逻辑 |
| b | 细胞组成堆叠条形图（新）：x=阶段（上人/下鼠两排），y=各组核占比 | GC 主导、PC 集中胎期、少突产后上升；直观暴露取样空洞 |
| c | 人 9 组特异性 dotplot（c1 式）：行=ataxia 基因（按主组→tau 排序），列=9 组；点大小=表达核比例、色=均值 CP10K；边列=tau | 每个致病基因的"作用位置" |
| d | 鼠 9 组特异性 dotplot（c1 式，同排布） | 与 c 对照，人鼠作用位置不一致的基因在此暴露 |
| e | 发育期×组占比热图（c2 式）：时间轴分 4 期（胚胎 CS18–CS22 / 胎儿 9–20 wpc / 婴幼儿 NB–3 岁 / 成人 42–52 岁），值=期内均值占比；缺组期 NA 斜纹；列头标注该期可行组数 n（人婴幼儿/成人期仅 5 组，占比与胎期 8 组不可直接比较） | 特异性何时建立（主组随发育期的切换） |
| f | 特异性建立时间轴（d1 式）：x=阶段，行=基因按 emergence 排序，色=主组 | 各基因获得细胞类型身份的时间先后 |

#### Fig 2｜致病基因何时起作用：发育时间窗

| panel | 内容与画法 | 怎么读 |
| --- | --- | --- |
| a | 人 ALL 组趋势热图（a1 式）：行=基因按 Spearman rho 降序，列=16 bin，色=CP10K z；不可行 bin 白底斜纹 | 全小脑尺度的升/降基因总览 |
| b | 鼠 ALL 组趋势热图（a1 式，12 bin） | 与 a 对照；出生后成熟波最清楚 |
| c | 分组热图选页（b1 式：Purkinje + 颗粒谱系） | 组内时间模式，不被组成变化稀释 |
| d | 轨迹小倍数（b2 式）：每基因一格，x=bin，y=CP10K；灰点=供体，实线=均值±sem，缺口断线不桥接，单供体 bin ▲ | 窗口的真实形状与证据强度 |
| e | 时间窗分类汇总（b4 式+新）：temporal_pattern（transient_peak / window_high / monotonic↑↓ / endpoint_high / flat）按组的堆叠条 + 代表性窗口基因轨迹 | 多少基因是"窗口型"、集中在哪组 |
| f | 组×基因 rho 热图（b3 式）：色=Spearman rho，星=显著（可行 bin≤6 用精确置换 p） | 一图看全趋势方向，供第 3/4 层统计取数 |

#### Fig 3｜机制类别的窗口是否系统性不同

| panel | 内容与画法 | 怎么读 |
| --- | --- | --- |
| a | 富集热图（新）：行=机制类（polyQ、离子通道/突触、DNA 修复、线粒体、脂代谢、胞骨架/胶质、纤毛/Joubert 类），列="组×期"模块（来自第 2 层窗口判定），色=富集倍数，星=置换 P | 哪一"类"致病基因整体落在哪类细胞的哪个发育期——本文核心统计图 |
| b | 窗类型×机制类列联（新，马赛克/热图）+ Fisher/置换检验 | 例如"通道类显著超额出现在出生后 window_high" |
| c | 峰时/开关点 ridgeplot（新）：x=发育时间，y=密度，每机制类一色（取每个基因的峰 bin 或 emergence bin） | 类别间峰时的系统性错位：先天类在前、成人起病类在后 |
| d | 代表基因轨迹小倍数（b2 式）：每类 1–2 个（如 ITPR1、CACNA1A、ATM、FXN、SYNE1、CAMTA1） | a–c 统计结论的原始证据 |

#### Fig 4｜时间窗能否解释临床发病年龄

| panel | 内容与画法 | 怎么读 |
| --- | --- | --- |
| a | 加权散点（新）：x=基因在其主组的窗口/峰时，y=发病年龄类别（先天/儿童/青少年/成人，序数），点大小=证据强度（可行 bin 数），趋势检验+置换 P | 窗口越晚 → 发病越晚 |
| b | 发病类别×峰时箱线（新）+ 两两 Wilcoxon + 校正 | 组间差异的显著性与效应量 |
| c | 疾病组×模块归属马赛克（新） | 先天类→前体/胚胎模块，成人类→Purkinje 成熟模块 |
| d | 疾病组泳道示意（新）：x=发育时间轴，行=疾病（SCA1/2/3/6/7、EA2、FRDA、AT、AOA1/2、Joubert、小脑发育不良类…），条带=致病基因窗口、色=细胞组 | 全文核心论点的主示意图 |
| e | polyQ 特例（新）：ATXN1/2/3/7、TBP 双轨迹——全组织表达（胚胎已高）vs Purkinje 特异占比（成熟期才升） | 解释"胚胎期就表达、为何成人才发病"：毒性需要成熟 PC 环境 |

#### Fig 5｜人鼠保守性与小鼠模型参考价值

| panel | 内容与画法 | 怎么读 |
| --- | --- | --- |
| a | dev_state 对齐 ALL 热图（x1 式）：33 个共享分化状态为列，人/鼠双联 | 用"分化状态"统一坐标代替日历年龄的跨物种轨迹比较 |
| b | 状态对齐分组热图（x2 式，Purkinje、颗粒谱系等 2–4 组） | 组内保守性 |
| c | 并排轨迹（x3 式）：选定基因人/鼠并排 | 单基因保守性目检 |
| d | rho 一致性散点（x4 式）：x=人 rho、y=鼠 rho、色=组；象限计数+总体同向率 | 趋势方向的一致性（引擎实测参考值：65 对中 57 对同向） |
| e | 保真度表热图（新）：行=疾病（或基因），列=保守性指标（rho 同向、轨迹 r、tau 主组一致、峰状态±1 一致），综合评分 | 每种 ataxia 的小鼠模型可信度评级——实验同事最关心的一张表 |

#### Fig 6｜外部验证

| panel | 内容与画法 | 怎么读 |
| --- | --- | --- |
| a | BrainSpan bulk 轨迹相关（新）：我们的组 pseudobulk 轨迹 vs BrainSpan 人小脑独立样本轨迹，逐基因/模块散点+r | 时间窗不是本数据集的技术假象 |
| b | ISH 示例（新）：Allen HDBR（胎）/成人脑原位杂交，2–4 个关键基因，与预测的细胞类型/时段对照；可选 Aldinger LCM 亚区（菱唇/脑室区）小 panel | 解剖学层面的独立佐证 |
| c | 模型鼠签名回投（新）：已发表 SCA1/SCA7/EA2/FRDA 模型差异基因在我们"组×期"模块中的富集条图 | 疾病扰动集中于对应正常窗口 → 支持"发育程序维持/失调病" |
| d | **Castonguay 成人 PC 秩基验证（新；数据已入库）**：109 供体成人 Purkinje（36,911 核）的 ataxia 基因平均表达秩 vs 我们鼠成年 PC / 人胎期 PC 的秩，Spearman rho + 方向一致率；重点核查成熟窗基因（ITPR1、CA8、GRID2 等）在成人是否维持高表达 | 成人 PC 表达状态的直接证据——最大空洞的补强 |
| e | 图形摘要（新）：时间轴×细胞组×疾病总览 | 一图总结全文 |

**补充图（详细清单）**：S1 QC 与可行性矩阵（组×bin 细胞数热图，NA 标注）；S2 全量组热图与轨迹（b1/b2 所有组）；S3 敏感性分析（剔除疾病供体、min-cell 阈值梯度、GU 剔除与否、GC 谱系拆/并、九组 vs 粗分组口径对比）；S4 glutamatergic_uncertain marker 抽查图（vs 成人 GC/胎期 PC 对照）；S5 Aldinger 胎儿期独立复核（rho 散点+方向一致率+峰阶段吻合+LCM 亚区定性）；S6 基因主表与保真度表全量。

## 四、计算主干：独立实现（参考 Depression_CB 的思路，不拷贝代码）

| 内容 | 实现方式 |
| --- | --- |
| 输入门禁、pseudobulk（供体×阶段×组）、趋势/特异性/时间窗/特异性建立统计、单物种图、人鼠 dev_state 对齐图、质控锚点、"样本量不足即标 NA 不插值"规则、绘图规范 | **按 `implementation_spec.md` §3 的模块设计重写**（`scripts/atx/` 包 + s00–s06）。借鉴 Depression_CB 已验证的口径（size_factor 归一化、稀疏指示矩阵聚合、全骨架 reindex、精确置换 p、缺组期不重归一、内部峰判定）与环境陷阱（pandas 3、scipy、backed 切片），但组分三层（8 主组 / UBC 描述组 / ALL）、基因键为 Ensembl ID、鼠时间列为 `author_stage`、emergence 定义为"主组占比 ≥ 0.6"，与旧代码的隐性假设不兼容，故不改造旧脚本 |
| 基因集（`genes/ataxia_genes.tsv`，带机制类/发病年龄/病理注释；别名与 Ensembl 映射） | **按 `gene_curation_plan.md` 策展**；数据输入只来自 /mnt/data/xyy/Ataxia_CB/（4 个 h5ad），不再读取 /mnt/data/xyy/Nature2024/ |
| 机制类×细胞×期富集（第 3 层）、发病年龄关联（第 4 层）、保守性评分与保真度表（第 5 层）、外部验证（第 6 层：Castonguay 成人 PC 秩基验证、Aldinger 胎儿复核、BrainSpan/GTEx、模型鼠签名） | **新开发**（s07–s11；算法见 `implementation_spec.md` §6–8；数据集间细胞类别与时间对齐规则见其 §4.5） |

注：抑郁面板与 ataxia 基因集有重叠（ITPR1、ATXN1、PPP2R2B=SCA12、SYNE1=ARCA1、GRID2、CAMTA1 等），这些基因在 Depression_CB 中的已有结果仅用作**方法学 sanity-check 锚点**（新实现在同一数据上应复现其方向性结论，如 ITPR1 的 Purkinje 特异性与鼠 Itpr1 产后上升），不构成结果可比性约束。

## 五、分析目标可行性核查（2026-09-22，基于实测数据）

| 目标（层） | 实测数据支撑 | 缺口 | 处理 |
| --- | --- | --- | --- |
| 1 细胞特异性 | 人胎期 9 组全部可行；人产后颗粒/中间神经元/少突/星形/免疫可行（GU 并入后 44 岁颗粒 bin 也可行）；鼠除免疫/成熟星形外全程可行 | 人成人 PC/DN/UBC 空白；鼠成年 PC 仅 21 核（9 周） | 成人 PC 用 Castonguay（normal 31,089 核、89 供体 ≥30）做秩基验证；DN 成人在鼠也仅个位数 → 结论止于 P0–P7，不外推 |
| 2 时间窗 | 颗粒谱系/中间神经元/少突/前体双物种窗口完整；PC 人胎期+鼠 E12.5→P7 密集 | 人 PC 窗止于胎期+NB；鼠 PC 成人点薄（P14=30、9wk=21，单 bin） | 人胎期窗与鼠产后窗拼接叙述；成人维持态由 Castonguay + BrainSpan/GTEx bulk 佐证 |
| 3 机制类富集 | 建立在第 1/2 层结果之上，9 组均有可行窗口覆盖；髓鞘类所需的少突组人产后 398–5,264/bin、鼠 P17.5 起 90–158/bin | 无 | 无需补数据 |
| 4 发病年龄关联 | 基因临床注释（GeneReviews/OMIM）+ 各基因窗口判定均产出 | 成人起病基因的"成熟维持窗"在人 PC 不可直接测 | 三角互证：鼠 P0→成年 + Castonguay 成人 PC 秩基 + BrainSpan/GTEx bulk |
| 5 跨物种保守 | 33 个共享 dev_state；引擎实测双物种显著趋势 65 对中 57 对同向 | 鼠免疫（6–56/bin 全 NA）与成熟星形（Bergmann 181）薄 | 这两组鼠侧仅描述；保真度表中相应基因标注证据等级 |
| 6 外部验证 | **验证数据集已全部入库**：Castonguay（成人 PC）、Aldinger（胎儿 SPLiT-seq，同平台）；BrainSpan/GTEx 为公开小体量下载；模型鼠 DEG 表可从文献/GEO 整理 | SCA7 Purkinje 富集库无数据（已排除）；**无公开的 ataxia 患者小脑 snRNA-seq** | 患者组织层面用"模型鼠签名回投"替代；作为局限性主动声明 |

**结论：六层目标全部可达成，无需新增测序或申请受限数据。** 剩余待办均为小体量下载/整理：① BrainSpan 人小脑发育 bulk 矩阵；② GTEx 成人小脑表达；③ 各 ataxia 模型鼠已发表 DEG 表；④ ataxia 基因集策展（~120–150 基因，注意 Castonguay 中 KIAA0226→RUBCNL、SOGA1 等别名映射）。

## 六、局限性（汇报时主动说明）

1. **人 Purkinje、深部核团与 UBC 出生后基本测不到**（该图谱公认的技术空洞；人产后实际可行的组仅颗粒谱系/中间神经元/少突/星形/免疫）→ 相应结论止于胚胎期；产后/成年窗以鼠 + Castonguay 成人 SPLiT-seq（秩基验证）+ BrainSpan/GTEx 佐证，不宣称过头。**鼠成年 Purkinje 亦薄**（9 周仅 21 核、P14 30 核）→ PC"成年维持"结论主要由 Castonguay 与 bulk 数据支撑。
2. 人部分阶段单供体（成人 42/44/46/52 岁各 1 例）、鼠早胚为合并文库 → 这些点无重复误差棒，图中明确标记。
2b. Castonguay 供体全部为 63–100 岁（老龄成人）且为 SPLiT-seq（中位 ~1,700 counts/核）→ 只能回答"成人（老龄）PC 是否维持该表达状态"，不能覆盖产后→青壮年的连续轨迹；中段轨迹靠 BrainSpan/GTEx bulk 衔接，且只做秩基/方向性结论。
3. 核测序低估胞质/线粒体基因（如 FXN）→ 线粒体类结论以外部 bulk 数据佐证。
4. 鼠免疫组每 bin 仅 6–56 个细胞（多数 bin 不可行）且免疫基因存在文库级 dropout；鼠成熟星形/Bergmann 极薄（79/181 个）→ 鼠的免疫组与星形成熟窗结论仅作描述。
5. 全部关联为"表达-表型"水平、非因果验证 → 定位为"图谱+假说生成"，验证靠第 6 层与后续实验。

## 七、时间与产出

- 第 0 周（已完成，2026-09-22）：数据目录统一入库（/mnt/data/xyy/Ataxia_CB/，4 个 h5ad）并核验；细胞组按 ataxia 病理重定义；GU 细胞身份查清；六层目标可行性核查通过。
- 第 1 周：基因集策展与临床注释（含别名映射）；流水线改参重跑（第 1/2/5 层图全部产出）；
- 第 2 周：第 3/4 层统计与图；
- 第 3 周：第 5 层保真度表 + 第 6 层验证（BrainSpan/GTEx 下载、Castonguay 成人 PC 秩基验证、Aldinger 胎儿复核、模型鼠 DEG 表整理）；
- 第 4 周：组图与初稿（目标 Nat Commun / Brain / Nat Neurosci 档）。
- 交付：6 张主图+补充图（S1–S6）；基因主表（每基因×人/鼠×全部统计量×临床注释）；小鼠模型保真度表。
