# Ataxia_CB 流水线实施规范（给实现者 LLM / 工程师）

> 目的：读完本文件即可在不询问的情况下把 `analysis_protocol.md` 的六层分析实现并跑通。
> 本文件只写**工程口径**（输入、字段、公式、阈值、产物、校验锚点、执行顺序）；科学理由见
> `analysis_protocol.md`，读图规则见 `figure_reading_guide.md`，基因表见 `gene_curation_plan.md`。
> 所有数字锚点均于 2026-10-13 在真实数据上实测。

---

## 0. 硬性约束（违反即视为实现错误）

1. **输入只读且唯一**：流水线只读取 `/mnt/data/xyy/Ataxia_CB/` 下的 4 个 h5ad，不得修改、不得在该目录写任何文件。**不再读取 `/mnt/data/xyy/Nature2024/`**：其 `human.h5ad`/`mouse.h5ad` 与本项目目录中的副本 md5 一致（f08ae756… / 354a7f77…），`mouse_tp*.h5ad` 是 Depression_CB 的再加工版本、本项目不用，`41586_2023_6884_MOESM3_ESM.xlsx` 为论文补充表仅供人工查阅（环境无 openpyxl）。任何脚本中出现 `Nature2024` 路径视为错误。
1b. **不拷贝 Depression_CB 代码**：本项目按 §3 的模块设计**独立重写**；Depression_CB 只作为"已验证可行的实现思路"参考（哪些陷阱、哪些口径有效），不得 `cp` 其脚本或整段粘贴函数。原因：组数（4→9+分层）、基因表结构（纯符号列表→带临床注释的表）、时间列（鼠 `timepoint`→`author_stage`）、外部数据集、第 3–6 层统计均不同，改造旧代码的隐性假设（4 组分母、固定锚点、固定组名）比重写更易出错。
2. **Python**：`/home/victor/.conda/envs/omics/bin/python`（3.14；scanpy 1.12.4、anndata 0.13.3、pandas 3.0.6、matplotlib 3.11、statsmodels 0.15）。**无 pyarrow → 禁 parquet**，中间产物只用 h5ad + TSV。无 openpyxl/pypdf。
3. **内存策略**：这些 h5ad 的 backed 模式列切片是 O(nnz)，会卡死；人/鼠/Aldinger 全载入内存（峰值 < 20 GB，机器 62 GB）。Castonguay（6.2 GB 文件、100 万核）**不要全载**：用 h5py 按 `raw/X` 的 `indptr` 做行切片（见 §8.1）。
4. **NA 政策**：组×bin < `min_cells`（30）→ 表中 NA、图中白底 45° 斜纹；**绝不插值、绝不填 0**；单供体或 pooled bin 的 sem 一律 NA。
5. **单一事实来源**：所有参数在 `config/config.yaml`；脚本不得硬编码阈值、组名、颜色。
6. **幂等**：每个脚本可重复运行，覆盖同名产物；产物目录带基因集版本号 `results/v{gene_set.version}/`。
7. **pandas 3 陷阱**（Depression_CB 在同一环境实测过，重写时同样适用）：`groupby([单键])` 返回 tuple key 需拆包；`observed=True` 不产生空组合 → 必须对全骨架 `reindex`；`sum(axis=1)` 默认 skipna → 缺组时显式置 NA；`idxmax` 遇全 NA 行抛错 → 用 `_idxmax_na_safe`；bool 列含 NA 用 nullable `"boolean"`；写入用 `.loc`。
8. **scipy `spearmanr` 在 rho=±1 时 p=0**：可行 bin ≤ 6 时用精确置换 p（`cb_stats._exact_perm_p`）。

---

## 1. 目录与代码来源

```
/home/victor/Projects/Ataxia_CB/
├── config/config.yaml              # §2
├── genes/                          # gene_curation_plan.md 的产物
│   ├── ataxia_genes.tsv  ataxia_gene_disease.tsv  gene_aliases.tsv  sanity_markers.tsv  curation_log.md
├── external/                       # §8 下载的小体量外部表(BrainSpan/GTEx/模型鼠 DEG)
├── scripts/
│   ├── atx/                         # 可复用模块(包),按 §3 独立实现
│   │   ├── __init__.py  config.py  genes.py  io.py  groups.py  pseudobulk.py
│   │   ├── stats_trend.py  stats_specificity.py  stats_pattern.py  stats_enrich.py
│   │   ├── align.py  external.py  sanity.py  plots.py
│   ├── s00_verify_inputs.py  s00b_derive_gene_table.py
│   ├── s01_build_cubes.py  s02_stats.py  s03_figures_species.py  s04_figures_crossspecies.py
│   ├── s05_export_tables.py  s06_report.py
│   ├── s07_category_enrichment.py  s08_onset_association.py  s09_conservation_fidelity.py
│   ├── s10_external_validation.py  s11_assemble_figures.py
│   ├── make_sample_figures.py       # 既有示意图脚本,保留不动
│   └── run_all.sh
├── tests/                           # pytest:合成小数据上的单元测试(§3.13)
└── results/v1.0/{intermediate,tables,figures,reports,logs}/{human,mouse,cross,external}
```

**第一步**：`git init`；创建 `scripts/atx/` 包并按 §3 自顶向下实现；先用 `tests/` 的合成数据（200 核 × 20 基因 × 3 组 × 5 bin）把 `pseudobulk`/`stats_*` 的口径跑对，再接真实数据。

---

## 2. `config/config.yaml`（完整字段；值为实测或决策值）

```yaml
meta:
  project: Ataxia_CB
  seed: 20261013
  paper_doi: "10.1038/s41586-023-06884-x"

paths:
  project_root: /home/victor/Projects/Ataxia_CB
  genes_table: /home/victor/Projects/Ataxia_CB/genes/ataxia_genes.tsv
  sanity_table: /home/victor/Projects/Ataxia_CB/genes/sanity_markers.tsv
  aliases_table: /home/victor/Projects/Ataxia_CB/genes/gene_aliases.tsv
  results_root: /home/victor/Projects/Ataxia_CB/results    # 实际写入 results/v{gene_set.version}
  human_h5ad: /mnt/data/xyy/Ataxia_CB/human.h5ad
  mouse_h5ad: /mnt/data/xyy/Ataxia_CB/mouse.h5ad            # 以原始 mouse.h5ad 为准(md5 354a7f77…)
  aldinger_h5ad: /mnt/data/xyy/Ataxia_CB/aldinger2021_fetal_snRNAseq.h5ad
  castonguay_h5ad: /mnt/data/xyy/Ataxia_CB/castonguay2026_adult_eQTL.h5ad
  external_dir: /home/victor/Projects/Ataxia_CB/external
  python: /home/victor/.conda/envs/omics/bin/python

gene_set:
  version: "1.0"                 # 必须等于 genes/curation_log.md 冻结版本
  tiers_main: [core]             # peripheral(genes/curation_log.md CL-023 派生层)不进任何统计层; effective tier=extended 仅经 sensitivity.include_extended_tier 纳入
  match_key: ensg_id             # 人及外部人类数据按 ENSG;鼠按 ensmusg_id;失败回退 symbol+aliases
  min_hit_rate: {human: 0.95, mouse: 0.90, aldinger: 0.90, castonguay: 0.90}

normalization: {method: size_factor, target_sum: 1.0e4}   # Sepp: X 原始计数; obs.size_factor = 作者基因过滤前的总 UMI(与 X 行和 r=0.99997,行和低≈1%),统一用 size_factor

species:
  human:
    h5ad_key: human_h5ad
    n_cells_expected: 180956
    n_vars_expected: 27410
    time_column: development_stage          # 存储 categories 乱序,顺序以 time_order 为准
    time_order: ["Carnegie stage 18","Carnegie stage 19","Carnegie stage 22",
      "9th week post-fertilization stage","11th week post-fertilization stage",
      "17th week post-fertilization stage","20th week post-fertilization stage",
      "newborn stage (0-28 days)","6-month-old stage","7-month-old stage","9-month-old stage",
      "3-year-old stage","42-year-old stage","44-year-old stage","46-year-old stage","52-year-old stage"]
    time_labels: [CS18,CS19,CS22,9w,11w,17w,20w,NB,6mo,7mo,9mo,3y,42y,44y,46y,52y]
    # 人等价坐标(用于 Fig3c/Fig4 的峰时轴;单位 = 粗略"受精后周"的对数刻度序号,只用于秩)
    periods:
      # 胚胎期止于 CS23(≈8 wpc),9 wpc 已属胎儿期(示意脚本把 9w 放胚胎期是错误的,不沿用)
      embryonic: ["Carnegie stage 18","Carnegie stage 19","Carnegie stage 22"]
      fetal:     ["9th week post-fertilization stage","11th week post-fertilization stage","17th week post-fertilization stage","20th week post-fertilization stage"]
      infant:    ["newborn stage (0-28 days)","6-month-old stage","7-month-old stage","9-month-old stage","3-year-old stage"]
      adult:     ["42-year-old stage","44-year-old stage","46-year-old stage","52-year-old stage"]
    donor_column: donor_id
    disease_column: disease       # normal / acute myocardial infarction(46y) / heart failure(44y) / congenital heart disease(NB 供体)
    exclude_disease_donors: false
    single_donor_bins: [CS19, 20w, 6mo, 7mo, 9mo, 42y, 44y, 46y, 52y]   # 实测;s00 校验
  mouse:
    h5ad_key: mouse_h5ad
    n_cells_expected: 115282
    n_vars_expected: 20409
    time_column: author_stage     # 注意:原始 mouse.h5ad 无 'timepoint' 列,用 author_stage
    time_order: [E10.5,E11.5,E12.5,E13.5,E14.5,E15.5,E17.5,P0,P4,P7,P14,"adult (9 weeks)"]
    time_labels: [E10.5,E11.5,E12.5,E13.5,E14.5,E15.5,E17.5,P0,P4,P7,P14,9wk]
    periods:
      embryonic: [E10.5,E11.5,E12.5,E13.5,E14.5,E15.5]
      perinatal: [E17.5,P0]
      postnatal: [P4,P7,P14]
      adult:     ["adult (9 weeks)"]
    donor_column: donor_id
    pooled_timepoints: [E10.5, E11.5]   # donor_id=='pooled' → sem NA、图 ▲
    exclude_disease_donors: false

# 鼠 bin → 人等价坐标(Fig3c/Fig4a 峰时统一轴;只做秩比较,不做定量)。人 bin 序号尺度:CS18=1,CS19=2,CS22=3,
# 9w=4,11w=5,17w=6,20w=7,NB=8,6mo=9,7mo=10,9mo=11,3y=12,42y=13 … 52y=16;非整数 = 落在两 bin 之间。
# 来源(2026-10-14 核对,Sepp 2024 Nature 正文/Extended Data Fig.2 与 bioRxiv 2021.12.20.473443 方法:动态时间规整,以鼠为焦点物种):
#   sourced: CS18–19(7 wpc)↔E11.5;CS22(8 wpc)↔E13.5;11 wpc↔E15.5;EGL/PL 出现 人 18 wpc↔鼠 E17.5;新生儿↔1 周龄鼠(P4 与 P7 互为重复);
#   interpolated(两 sourced 锚点之间按序插值,未在原文逐字确认):E10.5(早于 CS18)、E12.5、E14.5↔9 wpc、P0↔20 wpc、P14↔婴儿(6–9 mo)、P63↔成人;
#   人 3y(toddler)无鼠对应(介于 P14 与 P63 之间)。
stage_map_mouse_to_human_rank:
  E10.5: 0.5            # interpolated
  E11.5: 1.5            # sourced: CS18–19
  E12.5: 2.5            # interpolated
  E13.5: 3              # sourced: CS22
  E14.5: 4              # interpolated: 9 wpc
  E15.5: 5              # sourced: 11 wpc
  E17.5: 6              # sourced(组织学): 17–18 wpc
  P0: 7                 # interpolated: 20 wpc
  P4: 8                 # sourced: newborn
  P7: 8                 # sourced: newborn(与 P4 同级)
  P14: 10               # interpolated: infant 6–9 mo
  "adult (9 weeks)": 14.5   # interpolated: adult 42–52y
stage_map_source: {sourced: [E11.5, E13.5, E15.5, E17.5, P4, P7], interpolated: [E10.5, E12.5, E14.5, P0, P14, "adult (9 weeks)"]}
peak_source_rule:
  human_censored_groups: [Purkinje, DN, UBC]   # 人侧这些组产后全空 → 其峰时右删失,必须用鼠峰映射(peak_source=mouse)
  # 其余组(GranuleLineage/Interneuron/OligoLineage/Astro/Immune/Progenitor)人侧覆盖全轴 → 用人峰(peak_source=human);
  # Progenitor 两物种都在产后消失,峰早是生物学而非删失,不视为 censored

cell_groups:   # 两物种同一表达式;pandas eval(engine=python);列已 NA 安全化为 *_f;互斥(实测 overlap=0)
  Purkinje:       {expr: "author_cell_type == 'Purkinje'"}
  GranuleLineage: {expr: "author_cell_type in ['GC','GC/UBC','glutamatergic_uncertain'] | subtype in ['progenitor_RL','progenitor_RL_early']"}
  DN:             {expr: "author_cell_type in ['glut_DN','GABA_DN']"}
  Interneuron:    {expr: "author_cell_type == 'interneuron'"}
  UBC:            {expr: "author_cell_type == 'UBC'"}
  OligoLineage:   {expr: "author_cell_type == 'oligo'"}
  Astro:          {expr: "author_cell_type == 'astroglia' & subtype in ['astro_Bergmann','astro_parenh','astroblast','glioblast_PWM']"}
  Progenitor:     {expr: "author_cell_type in ['VZ_neuroblast','NTZ_neuroblast','NTZ_mixed'] | (author_cell_type == 'astroglia' & subtype.str.startswith('progenitor_') & ~subtype.isin(['progenitor_RL','progenitor_RL_early'])) | (author_cell_type == 'astroglia' & subtype == '__NA__' & dev_state == 'progenitor')"}
  Immune:         {expr: "author_cell_type == 'immune'"}
  ALL:            {expr: null}
group_tiers:          # 统计口径
  main: [Purkinje, GranuleLineage, DN, Interneuron, OligoLineage, Astro, Progenitor, Immune]   # 进 tau/fraction 分母、富集检验
  descriptive: [UBC]  # 只出图与表,不进分母与检验(人产后空、鼠仅 5 bin)
  reference: [ALL]

sensitivity:
  drop_glutamatergic_uncertain: false     # true → GranuleLineage 表达式去掉 glutamatergic_uncertain
  split_granule_progenitor: false          # true → 把 progenitor_RL(_early) 拆为 GranuleProgenitor 组
  include_ubc_in_main: false
  min_cells_grid: [30, 50, 100]
  exclude_disease_donors: false
  include_extended_tier: false              # true → 敏感性层纳入 effective tier=extended; peripheral(genes/curation_log.md CL-023 派生层)恒排除, 仅描述性展示
  onset_thresholds_grid: [[1,13,30],[1,10,25],[2,16,35]]
  cacna1a_as_polyQ: false                  # true → CACNA1A 的 mech_class 改 polyQ_repeat
  collapse_adult_bins: false               # true → 人 42/44/46/52y 合并为单一 adult bin(4 供体)
  peak_axis: mixed                         # mixed | human_only | mouse_only(§6 三套结果)

bins:
  min_cells: 30
  min_cells_total_for_trend: 100
  min_bins_for_trend: 5
  per_group_override:
    UBC:    {min_bins_for_trend: 4}
    Immune: {min_bins_for_trend: 4}
  na_policy: blank_hatch

stats:
  trend: {method: spearman, levels: [bin, donor], fdr: benjamini-hochberg, fdr_family: species_group, exact_perm_max_bins: 6}
  specificity:
    scale: mean_cp10k_linear
    pseudocount: 1.0e-3
    aggregation: bin_equal_weight
    background: main_groups          # 分母 = group_tiers.main(8 组);UBC/ALL 不入
    tau_threshold: 0.6
    fraction_threshold: 0.5
    min_expr: 0.1
  tau_by_bin: {min_groups: 4, tau_threshold: 0.6, fraction_threshold: 0.6}   # emergence = 首个主组 fraction ≥ 0.6 的 bin
  temporal_pattern: {pulse_min_peak_vs_flank: 0.5, window_min_ratio: 2.0, mono_rho: 0.4}
  enrichment:        # s07
    n_perm: 10000
    min_genes_per_class: 5            # polyQ_repeat 预计仅 6 个基因(ATXN1/2/3/7、TBP、ATN1),是核心假设类,阈值取 5 保留之;其结果须标低功效
    modules: group_x_period            # 模块 = 主组 × period(人 4 期/鼠 4 期)
    module_assign_rule: peak_period_in_max_group   # 基因→模块:主组(tau argmax)× 该组内峰 bin 所在 period
  onset:             # s08
    continuous: log1p_onset_typical_y
    categorical: onset_class
    n_perm: 10000
    stratify_by: mech_class            # 分层置换:只在同机制类内打乱 onset
  conservation:      # s09
    rho_sign_weight: 0.3
    traj_r_weight: 0.3
    max_group_weight: 0.2
    peak_state_weight: 0.2
    peak_state_tolerance: 1

state_alignment:
  vocab: dev_state                     # 实测共享 33 类(见 §4.4)
  min_cells: 30
  order_rule: consensus_mean_time_rank
  group_states:
    Purkinje:       [Purkinje_diff, Purkinje_defined, Purkinje_maturing]
    GranuleLineage: [GCP, GCP/UBCP, GC_diff_1, GC_diff_2, GC_defined]
    DN:             [glut_DN_defined, GABA_DN_defined]
    Interneuron:    [interneuron_diff, interneuron_defined]
    UBC:            [UBC_diff, UBC_defined]
    OligoLineage:   [oligo_progenitor, oligodendrocyte]
    Astro:          [glioblast, astrocyte]
    Progenitor:     [progenitor, VZ_neuroblast_1, VZ_neuroblast_2, VZ_neuroblast_3, NTZ_neuroblast_1, NTZ_neuroblast_2, NTZ_neuroblast_3]
    Immune:         [immune]

external:
  castonguay:
    use_disease: [normal]              # 排除 essential tremor(16 供体/254,076 核)
    cell_type_column: author_cell_type
    purkinje_label: Purkinje
    other_groups: {Granule: GranuleLineage, Bergmann: Astro, Astrocytes: Astro, MLI_1: Interneuron, MLI_2: Interneuron, Golgi: Interneuron, PLI: Interneuron, Oligodendrocytes: OligoLineage, OPC: OligoLineage, Microglia: Immune, UBC: UBC}
    counts_layer: raw/X                # X 为 log-normalized 浮点;raw/X 为整数计数(已验证)
    min_pc_per_donor: 30
    subsample_pc_per_donor: 300        # 供体级 pseudobulk 前每供体最多抽 300 PC(平衡供体权重;seed 固定)
  aldinger:
    cell_type_column: fig_cell_type
    # 以作者 fig_cell_type 为准(CELLxGENE cell_type 本体列把 eCN/UBC、iCN 都标成 interneuron,不可用);
    # 依据 Aldinger 2021 Fig.1 簇定义:iCN=inhibitory cerebellar nuclei(→DN),eCN/UBC=兴奋性核团与 UBC 混合簇(不入主组),
    # PIP=PAX2+ interneuron progenitors(→Interneuron,谱系前体)。实现前对照原文复核一次。
    group_map: {H-PC: Purkinje, H-RL: GranuleLineage, H-GCP: GranuleLineage, H-GN: GranuleLineage, H-iCN: DN, H-eCN/UBC: mixed_DN_UBC, H-PIP: Interneuron, H-MLI: Interneuron, H-BG: Astro, H-Ast: Astro, H-Glia: Astro, H-OPC: OligoLineage, H-Committed OPC: OligoLineage, H-Microglia: Immune}
    descriptive_only: [mixed_DN_UBC]
    all_only: [H-Brainstem, H-Choroid, H-Meninges, H-Endothelial, H-Pericytes, H-Ast/Ependymal, H-BS Choroid/Ependymal]
    sepp_bin_pairs: {"9-10w": "9th week post-fertilization stage", "11-12w": "11th week post-fertilization stage", "17-18w": "17th week post-fertilization stage", "20-21w": "20th week post-fertilization stage"}   # 14w 无对应
    counts_layer: raw/X
    time_column: development_stage     # 9 个 wpc 阶段,13 供体;合并为 bins: [9-10w, 11-12w, 14w, 17-18w, 20-21w]
    time_bins: {"9-10w": ["9th week post-fertilization stage","10th week post-fertilization stage"], "11-12w": ["11th week post-fertilization stage","12th week post-fertilization stage"], "14w": ["14th week post-fertilization stage"], "17-18w": ["17th week post-fertilization stage","18th week post-fertilization stage"], "20-21w": ["20th week post-fertilization stage","21st week post-fertilization stage"]}
  brainspan: {file: external/brainspan_cerebellum_rpkm.tsv, structure: CBC, min_samples_per_stage: 2}
  gtex:      {file: external/gtex_v8_cerebellum_tpm.tsv, tissues: ["Brain - Cerebellum", "Brain - Cerebellar Hemisphere"]}
  model_mouse_deg: external/model_mouse_deg_signatures.tsv   # 列: model, gene_mouse, direction, log2fc, padj, source_pmid

figures:
  formats: [pdf, png]
  dpi: 300
  sequential_cmap: Blues
  diverging_cmap: RdBu_r
  tau_cmap: YlOrBr
  na_style: {facecolor: white, hatch: "///", edgecolor: "#c3c2b7"}
  group_colors: {Purkinje: "#d62728", GranuleLineage: "#1f77b4", DN: "#9467bd", Interneuron: "#2ca02c", UBC: "#8c564b", OligoLineage: "#ff7f0e", Astro: "#e377c2", Progenitor: "#7f7f7f", Immune: "#17becf", ALL: "#898781"}
  mech_colors: {polyQ_repeat: "#d62728", channel_synapse: "#1f77b4", dna_repair: "#2ca02c", mito_metabolic: "#ff7f0e", cilia_congenital_dev: "#9467bd", cytoskeleton_transport_glia: "#8c564b", proteostasis_RNA: "#e377c2", other: "#7f7f7f"}
  heatmap_page_size: 60
  row_order_rule: rho_desc
```

---

## 3. 模块设计与实现思路（独立实现；Depression_CB 仅作思路参考）

每个模块给出：职责 → 借鉴的思路（Depression_CB 已验证有效的口径/陷阱）→ 本项目的差异化设计 → 关键函数签名。实现者按签名落地即可；不得从 Depression_CB 复制代码。

### 3.1 `atx/config.py`
- 职责：加载 `config.yaml`，派生每物种的 `bin_order`（1..n）、`period_of`、短标签、组层级；解析 `results/v{version}/` 目录。
- 借鉴：所有参数集中于 YAML；时间顺序显式声明（h5ad 中 categories 乱序）。
- 差异：组分三层（`main`/`descriptive`/`reference`），所有统计函数通过 `main_groups(cfg)` 取分母组，`plot_groups(cfg)` 取出图组；`gene_set.version` 校验。
- 签名：`load(path) -> Cfg`；`Cfg.species(sp) -> SpeciesCfg`（属性 `time_order, bin_order, period_of, labels, donor_col, time_col, pooled`）；`Cfg.main_groups / plot_groups / results(*parts)`。

### 3.2 `atx/genes.py`
- 职责：读 `ataxia_genes.tsv` + `sanity_markers.tsv` + `gene_aliases.tsv`；建立 panel（core & in_panel ∪ sanity）；对任意数据集 `var` 做映射。
- 借鉴：映射结果固化为表（`gene_mapping.tsv`）并写明 `match_method`；panel 基因不得落在重复符号上。
- 差异：主键是 **ENSG/ENSMUSG**（四个数据集 `var.index` 都是 Ensembl ID），符号只是回退；返回 `GeneMap`（DataFrame：`gene_key, ensg, ensmusg, row_<dataset>, status_<dataset>, match_method_<dataset>, in_panel, is_sanity, mech_class, onset_*`）。
- 签名：`load_panel(cfg) -> pd.DataFrame`；`map_to_var(panel, var, id_col='index', symbol_col='feature_name', species) -> pd.DataFrame`；`hit_rate(gm, dataset) -> float`。

### 3.3 `atx/io.py`
- 职责：三种读取模式。(a) `read_full(path)`：anndata 全载（Sepp 人/鼠、Aldinger）；(b) `read_obs_var(path)`：只用 h5py + `anndata.io.read_elem` 读 obs/var，不触矩阵；(c) `read_rows_csr(path, rows, layer='raw/X', chunk=5000)`：按行号从 CSR 的 `indptr/indices/data` 切片（Castonguay 专用）。
- 借鉴：backed 列切片 O(nnz) 会卡死 → 永远先切行/全载再切列。
- 差异：新增行切片读取器，并对 `raw/X` 做整数性断言（`assert np.allclose(data, np.round(data))`）。
- 签名：如上；`read_rows_csr` 返回 `(scipy.sparse.csr_matrix, var_index)`，内部按 `rows` 排序后分块 `h5py` 读取，再按原顺序还原。

### 3.4 `atx/groups.py`
- 职责：obs 派生列 + 细胞组指派。
- 借鉴：分组表达式用 `pandas.eval(engine='python')`；分类列先做 NA 安全化（`<NA>` → 哨兵字符串）；敏感性开关只打标 `use_cell`，由 s01 同步过滤 X 与 obs（避免行错位）。
- 差异：① 要安全化的列扩展为 `author_cell_type, subtype, dev_state, cell_type`；② `subtype.str.startswith(...)` 之类表达式**先预计算布尔辅助列**（`sub_is_progenitor`, `sub_is_rl`, `sub_na`）再 eval，避免依赖 `eval` 的 `.str` 支持；③ 断言组互斥（除 ALL）与"ungrouped 只进 ALL"，并把每组计数与 §4.1 锚点比对；④ 鼠时间列 `author_stage`，人 `development_stage`；⑤ 敏感性开关 `drop_glutamatergic_uncertain` / `split_granule_progenitor` 在此处改写表达式。
- 签名：`derive_obs(obs, spcfg, cfg) -> obs`；`assign(obs, cfg) -> pd.DataFrame[cell_pos, cell_group]`；`check_anchors(obs_long, species, cfg)`。

### 3.5 `atx/pseudobulk.py`
- 职责：细胞 × 基因 → (donor × bin × group) → (group × bin) 两级聚合。
- 借鉴：① 归一化 = `X / obs.size_factor × 1e4`（Sepp 的 `size_factor` ≡ 总计数，s00 复验）；② 用稀疏指示矩阵 `S @ X` 做一次性聚合，层 `mean_log`、`mean_cp10k`、`pct_pos`；③ 全骨架 reindex 保证不可行 bin 显式存在；④ 单供体/pooled bin 的 sem = NA（硬断言）。
- 差异：① 骨架为 `plot_groups × bins`；② `feasible` 由 `min_cells` 参数化，便于 `min_cells_grid` 敏感性复用同一 donor cube（**只重算 bin 表，不重算 donor cube**）；③ 对外部数据集（Aldinger/Castonguay）复用同一函数，输入改为 `size_factor = 行和`。
- 签名：`normalize(X, size_factor, target) -> (X_log, X_cp)`；`donor_cube(X_log, X_cp, obs, obs_long, var) -> AnnData`；`bin_table(cube, spcfg, cfg, min_cells) -> pd.DataFrame`；`spot_check(X_log, obs, obs_long, bin_table, n=3)`（随机抽 3 个 (gene, group, bin) 用 numpy 直接重算比对）。

### 3.6 `atx/stats_trend.py`
- 职责：Spearman 趋势（bin 级主、donor 级辅）、BH-FDR。
- 借鉴：可行 bin ≤ 6 用精确置换 p（scipy 在 rho=±1 时 p=0）；FDR 族 = 物种 × 组；时间轴统一用 `bin_order`。
- 差异：`min_bins_for_trend` 按组可覆盖（UBC/Immune 4）；输出 `low_power` 标记；置换 p 的排列表缓存按 n 建一次。
- 签名：`spearman_exact(t, x, exact_max_n=6) -> (rho, p)`；`trend_table(bin_table, cube, genes, group, spcfg, cfg) -> pd.DataFrame`。

### 3.7 `atx/stats_specificity.py`
- 职责：type_fraction、tau、max_group、specific_call（global 与 by-period）；tau×bin 与 emergence。
- 借鉴：线性 CP10K、bin 等权；缺组期不隐式重归一；`idxmax` NA 安全。
- 差异：① 分母 = 8 主组（UBC/ALL 不入）；② 可行组数 < `min_groups`(4) 的期/bin → fraction/tau 一律 NA，并输出 `groups_used` 字符串列与 `n_groups_used`（让读者知道分母由哪些组构成——人产后期必然缺 PC/DN/Progenitor）；**tau 与 fraction 的数值只在 `groups_used` 相同的 bin/期之间可比**（tau 的分母 n−1 随组数变；fraction 的分母随组集变），因此 Fig 1e 每期列头必须标注 n_groups，emergence 只在 `n_groups_used ≥ 4` 的 bin 上判定，且 s02 输出 `fraction_8group`（固定 8 组分母、缺组按 NA 不参与、仅当 8 组全可行才非 NA）作严格可比版本供补充表；③ **emergence 定义改为"首个主组 fraction ≥ 0.6 的 bin"**（主组 = global `max_group`），与读图指南一致；④ 同时给出 `fraction_in_feasible`（在可行组内的占比）与 `max_group_period`（分期主组），供 Fig 1e。
- 签名：`aggregate(bin_table, genes, cfg, period=None) -> long`；`scores(agg, genes, cfg) -> pd.DataFrame`；`tau_by_bin(bin_table, genes, cfg)`；`emergence(tau_bin, max_group_map, cfg)`。

### 3.8 `atx/stats_pattern.py`
- 职责：temporal_pattern 标签与 `peak_by_group`。
- 借鉴：内部峰判定（峰两侧都有可行 bin）、`peak_vs_flank`、`period_ratio`、标签优先级 transient > window > monotonic > endpoint_high > flat。
- 差异：① 轨迹存在中间缺口（人 Interneuron 6mo/44y、Astro 11w）时，"两侧有可行 bin"按**可行 bin 序列**判定而非相邻 bin；② 新增 `half_max_window`（CP10K ≥ 峰值 50% 的连续可行 bin 区间），供 Fig 4d 泳道宽度；③ 输出 `peak_by_group.tsv`（gene × group：peak_bin, peak_period, peak_cp, n_feasible_bins, half_max_lo, half_max_hi）。
- 签名：`patterns(bin_table, genes, group, cfg) -> pd.DataFrame`；`label(df, cfg) -> pd.Series`；`peak_table(bin_table, genes, cfg)`。

### 3.9 `atx/stats_enrich.py`（第 3/4 层，全新）
- 职责：模块指派、机制类 × 模块置换富集、分层置换相关、Kruskal/Mann-Whitney + Cliff's delta、Fisher 列联。
- 思路：富集 = 观测/期望的 log2 倍数 + 标签置换（打乱机制类，固定模块）；发病年龄关联的置换在 `mech_class` 内分层打乱 y，消除"机制类与发病年龄共线"的混杂；所有置换用 `np.random.default_rng(cfg.seed)`。
- 签名：`assign_modules(spec_global, peak_by_group, cfg) -> pd.DataFrame`；`perm_enrichment(labels_a, labels_b, n_perm) -> pd.DataFrame[a, b, obs, exp, log2fe, p, q]`；`stratified_perm_spearman(x, y, strata, n_perm) -> (rho, p_strat, p_naive)`；`group_compare(values, groups) -> (kw_p, pairwise_df)`；`cliffs_delta(a, b)`.

### 3.10 `atx/align.py`（第 5 层）
- 职责：跨物种 dev_state 剖面、共识状态排序、组内状态子集、保真度评分。
- 借鉴：物种内 z-score、独立 colorbar、只比形状；共识序 = 各物种内状态的平均 `bin_order` 秩再取两物种秩均值。
- 差异：状态剖面也按 `min_cells` 判可行；保真度 4 指标 + 权重 + `n_metrics` 证据等级（§8.1）。
- 签名：`state_profiles(X_log, obs, genes, states, min_cells, group_mask=None)`；`consensus_order(profiles_h, profiles_m)`；`fidelity(trend_h, trend_m, spec_h, spec_m, prof_h, prof_m, cfg)`.

### 3.11 `atx/external.py`（第 6 层）
- 职责：Castonguay 行切片 + 供体级 pseudobulk；Aldinger 全载 + 分组分 bin pseudobulk；BrainSpan/GTEx 表读取；模型鼠签名富集。
- 设计：两套外部单核数据都走 `pseudobulk.py` 的同一函数（`size_factor = raw 行和`），只是 obs 派生不同（§4.5 的映射表）。
- 签名：`castonguay_pc_pseudobulk(cfg, gm) -> pd.DataFrame[gene, mean_cp10k, n_donors, n_donors_detected]`；`aldinger_cube(cfg, gm) -> (cube, bin_table)`；`rank_compare(a, b) -> dict(rho, p, concordance, outliers)`.

### 3.12 `atx/sanity.py` 与 `atx/plots.py`
- sanity：锚点字典 = §4.1/§4.2 实测值；marker 期望组从 `sanity_markers.tsv` 读；鼠 dropout 基因在 s00 现场计算（全库 pct+ < 1% 即豁免）并写入 `logs/s00_dropout.tsv`——**2026-10-14 已在原始 `mouse.h5ad` 实测**：Neurod1 0.66% / Ptprc 0.07% / Aif1 0.06% / Csf1r 0.08% / Tmem119 0.05%（免疫组内 pct+ ≤ 4%）→ 豁免；Cx3cr1 / P2ry12 / Hexb 在免疫组 pct+ 52% / 42% / 50% → 作鼠免疫锚点；人侧 PTPRC/CSF1R/CX3CR1/P2RY12 在免疫组 59–68% 正常。`size_factor`：与当前 X 行和相关 0.99997 但行和约低 1%（作者在基因过滤前计算总 UMI）→ 归一化统一用作者 `size_factor`，外部数据集无此列时用 `raw/X` 行和；趋势方向锚点 WARNING 级。
- plots：热图（imshow + NA 斜纹 + 分页 + 全局色标）、轨迹小倍数（donor 灰点、sem 带只在连续可行段内、单供体 ▲、疾病供体 △）、dotplot（大小=pct、色=log 均值、tau 边列）、emergence 时间轴、ridge、mosaic、swimlane、rho 散点、保真度表。所有颜色/斜纹样式从 config 读。

### 3.13 测试（`tests/`）
- `test_pseudobulk.py`：合成数据上验证 bin 均值 = 手算；单供体 sem 为 NA；不可行 bin 为 NA。
- `test_specificity.py`：3 组可行时 tau/fraction 为 NA（min_groups=4）；8 组时 tau 公式 = 手算。
- `test_pattern.py`：构造内部峰/端点峰/中间缺口三条轨迹，断言标签。
- `test_enrich.py`：固定种子下 `perm_enrichment` 对无关标签 p 接近均匀；`stratified_perm_spearman` 在"y 完全由 strata 决定"时 `p_strat` 不显著而 `p_naive` 显著。
- `test_io.py`：`read_rows_csr` 与全载后切行结果一致（用 Aldinger 前 2,000 行）。

## 4. 数据事实与锚点（s00/s01 硬校验值）

### 4.1 9 组细胞数（互斥；ungrouped 只进 ALL）

| 组 | 人 | 鼠 |
| --- | ---: | ---: |
| Purkinje | 45,341 | 18,271 |
| GranuleLineage | 72,112 | 35,247 |
| DN | 6,352 | 6,656 |
| Interneuron | 6,805 | 5,286 |
| UBC | 3,847 | 1,208 |
| OligoLineage | 11,064 | 881 |
| Astro | 3,150 | 1,843 |
| Progenitor | 25,543 | 36,419 |
| Immune | 1,770 | 371 |
| ungrouped（仅 ALL） | 4,972 | 9,100 |
| 合计 | 180,956 | 115,282 |

`sensitivity.drop_glutamatergic_uncertain=true` 时人 GranuleLineage = 71,005。鼠 `astroglia` 且 `subtype` 为 NA 的 937 核全部 `dev_state=='progenitor'`，已被 Progenitor 表达式第三项收纳。

### 4.2 组 × bin 细胞数与可行性（≥30）

**人**（列 = CS18…52y）：

```
Purkinje        80  656 25993 9300 6136 1887 1231   53    0    0    0    4    1    0    0    0   → 可行 8 (CS18–NB)
GranuleLineage 354  206   875  762 2155 3144 2198 13187 2969 4606 4070 13432 14102 496 6214 3342 → 16
DN              58  463  3807 1163  525  336    0    0    0    0    0    0    0    0    0    0   → 6 (CS18–17w)
Interneuron      2    6   282  812 1220 1594  933  930   17   78   62   87  214   27  399  142   → 12 (6mo/44y NA, CS18/19 NA)
UBC              0    0    51   95 1279 1920  444   58    0    0    0    0    0    0    0    0   → 6
OligoLineage     0    0     2    0   73  289  124  208   15   37   71  398  466 3209  908 5264   → 11
Astro            0    0     2    0    8  453   73  585  134  132  199  161  583  167  392  261   → 11
Progenitor    4832 4436  6851 1301 1658 3251 2582  619    5    3    3    2    0    0    0    0   → 8
Immune           6    9   166   45   47  110   46  142    3   11    8   46  115  675   49  292   → 11
```

**鼠**（列 = E10.5…9wk）：

```
Purkinje         21   38 2420 7664 2272 2756 1517 1111  231  190   30   21  → 10 (E12.5–P14;P14 恰 30;9wk NA)
GranuleLineage 1477  876  688 1461  657 1474 2727 4107 4692 4176 4781 8131  → 12
DN               16  101 1597 2557  825  741  437  284   37   45    9    7  → 9 (E11.5–P7)
Interneuron      38   60  211  833  311  552  698  759  560  498  400  366  → 12
UBC               0    0    0    0   27  102  311  536  115   86   17   14  → 5
OligoLineage      6    4    0    1    6    6   90  158  155  155  149  151  → 6 (E17.5–9wk)
Astro             0    0    0   13   70   76  245  458  403  307  125  146  → 8
Progenitor     7595 10634 4760 5656 975 1574 1487 1732 1013  896   67   30  → 12
Immune           56   41   14   27    6   37   16   41   41   52   21   19  → 6 (低功效,per_group_override)
```

供体数：人 CS19/20w/6mo/7mo/9mo/42y/44y/46y/52y 各 1 供体（sem NA）；鼠 E10.5/E11.5 为 pooled（sem NA）；其余 2–4 供体。

### 4.3 外部数据集事实

| 项 | Aldinger 2021 | Castonguay 2026 |
| --- | --- | --- |
| 形状 | 69,174 × 26,326 | 1,004,112 × 47,789 |
| `X` | **log-normalized 浮点** | **log-normalized 浮点** |
| `raw/X` | 整数计数，同形状，`raw/var` 索引 = `var` 索引 | 同左 |
| 基因键 | `var.index` = ENSG；`var.feature_name` 符号 | 同左；KIAA0226 已更名 RUBCNL；SOGA1 符号缺失（按 ENSG 匹配） |
| 细胞标签列 | `fig_cell_type`（21 类，H-PC 25,724） | `author_cell_type`（14 类，Purkinje 36,911） |
| 供体/年龄 | 13 供体，9–21 wpc；PC/供体 4–12,182 | 109 供体 63–100 岁；**normal 93 供体 / 31,089 PC / 89 供体 ≥30 PC（中位 257）**；ET 16 供体排除 |
| 每核中位基因数 | ~1,146 | ~719（低灵敏度 → 只做秩） |

### 4.4 共享 `dev_state`（33 类，跨物种对齐轴）

GABA_DN_defined, GCP, GCP/UBCP, GC_defined, GC_diff_1, GC_diff_2, NTZ_neuroblast_1/2/3, Purkinje_defined, Purkinje_diff, Purkinje_maturing, UBC_defined, UBC_diff, VZ_neuroblast_1/2/3, astrocyte, erythroid, glioblast, glut_DN_defined, immune, interneuron_defined, interneuron_diff, isth_N_defined, isth_N_diff, meningeal, mural/endoth, noradrenergic, oligo_progenitor, oligodendrocyte, parabrachial, progenitor。人独有 GABA_MB/glut_DN_maturing/isthmic_neuroblast，鼠独有 GC/UBC_diff/MBO/ependymal/motorneuron/neural_crest_progenitor → 不进对齐。

---

### 4.5 数据集之间细胞类别与发育时间点的对齐规则

**总原则**：① 任何两个数据集之间**不合并表达矩阵**，只在"组 × 时间 bin 的 pseudobulk 统计量"层面比较；② 细胞类别对齐 = 把每个数据集作者自己的标签**映射到本项目 9 组**（映射表固定在 config，不做重新聚类/标签迁移）；③ 时间对齐分三档——**同一尺度直接对齐**（Aldinger ↔ Sepp 人，都是 wpc）、**用共享分化状态对齐**（Sepp 人 ↔ 鼠，`dev_state`）、**不对齐只做单点秩比较**（Castonguay 成人/老年）。

#### A. Sepp 人 ↔ Sepp 鼠（主数据，第 5 层）

| 维度 | 规则 |
| --- | --- |
| 细胞类别 | 两物种用**同一套** `author_cell_type` + `subtype` 词表（Sepp 作者联合注释），9 组表达式对两物种一字不差（§2 `cell_groups`，实测互斥、无歧义）。人独有标签 `glutamatergic_uncertain`（已核实为成熟 GC → 颗粒谱系）；鼠独有 `MBO`、`ependymal`、`motorneuron`、`neural_crest_progenitor` 只进 ALL。 |
| 时间：主分析 | **不对齐**。两物种各自用自己的 bin 序号做 Spearman、峰判定、temporal_pattern；人 16 bin、鼠 12 bin 完全独立；period 名称也各自定义（人 embryonic/fetal/infant/adult，鼠 embryonic/perinatal/postnatal/adult）。 |
| 时间：跨物种比较（Fig 5a–c） | 用 **`dev_state`**（33 个共享分化状态，§4.4）作共同横轴：对每个状态取该物种中该状态所有核的均值（ALL）或组内核的均值（`group_states`），物种内逐基因 z-score；状态顺序 = 共识序（每物种内按状态的平均 `bin_order` 排秩，两物种秩平均）。这是"按细胞分化进度而非日历年龄对齐"的做法，不依赖任何人鼠年龄换算。 |
| 时间：峰时合并到一根轴（Fig 3c、Fig 4a/d） | 只在需要把两物种基因放到同一横轴时才用 `stage_map_mouse_to_human_rank`（§2）：鼠 bin → 人 bin 序号的近似秩。**已于 2026-10-14 对照 Sepp 2024 核定**：作者用动态时间规整（以鼠为焦点物种）给出 CS18–19↔E11.5、CS22↔E13.5、11 wpc↔E15.5、新生儿↔P4/P7，组织学上 EGL/PL 在人 18 wpc 与鼠 E17.5 同期出现；其余鼠 bin 在锚点间按序插值并在 config 标 `interpolated`。图注须写明"锚点来自 Sepp 2024，其余插值"。规则（`peak_source_rule`）：主组 ∈ {Purkinje, DN, UBC} 的基因人侧峰时右删失 → 用鼠峰映射（`peak_source=mouse`，图中虚线边框）；其余组用人峰。同时报告仅人/仅鼠两套作为对照；统计只做秩（Spearman/Kruskal），不做定量。 |
| 结果口径 | 人产后 PC/DN/UBC/Progenitor 空 → 这些组的"成熟期"只能来自鼠；Fig 5 的保真度评分即是"鼠能否代人"的校准，凡保真度 low 的基因在 Fig 4 中标注。 |

#### B. Aldinger 2021（胎儿 SPLiT-seq）↔ Sepp 人（S5 复核）

| 维度 | 规则 |
| --- | --- |
| 细胞类别 | `fig_cell_type`（21 类）→ 9 组映射（config `external.aldinger.group_map`，**本版已修正**）：H-PC→Purkinje；H-RL、H-GCP、H-GN→GranuleLineage；**H-iCN（inhibitory cerebellar nuclei）→DN**；**H-eCN/UBC（excitatory cerebellar nuclei 与 UBC 混合簇）→ 不入主组（`mixed_DN_UBC`，仅描述）**；H-PIP（PAX2⁺ 中间神经元前体）→Interneuron（注 `lineage_progenitor`）；H-MLI→Interneuron；H-BG、H-Ast、H-Glia→Astro；H-OPC、H-Committed OPC→OligoLineage；H-Microglia→Immune；H-Brainstem、H-Choroid、H-Meninges、H-Endothelial、H-Pericytes、H-Ast/Ependymal、H-BS Choroid/Ependymal→ 只进 ALL。Aldinger 没有与 Sepp `Progenitor`（VZ/NTZ neuroblast）直接对应的簇 → 该组在 S5 中 NA。CELLxGENE 的 `cell_type` 本体列把 H-eCN/UBC 标成 "interneuron"、H-iCN 标成 "inhibitory interneuron"，**不可用**，必须用 `fig_cell_type`。**2026-10-14 用 `raw/X` 实测 marker（均值 CP10K）**：H-PIP PAX2 12.9 / GAD1 7.8 / GAD2 6.7（PAX2⁺ 中间神经元前体 ✓）；H-eCN/UBC EOMES 3.9 / LMX1A 11.6（UBC 与兴奋性核团混合 ✓，不入主组）；H-PC FOXP2 35.6 / SKOR2 2.0 ✓；H-GN RBFOX3 16.8 ✓；H-RL LMX1A 5.9 / MKI67 5.5 ✓；**H-iCN 的 GAD1/GAD2 仅 0.8–0.9、SLC17A6 1.1，抑制性身份在本数据中不明确** → Aldinger 的 DN 组降为描述级，S5 的定量复核只用 Purkinje / GranuleLineage / Interneuron / Astro / OligoLineage / Immune 六组。 |
| 时间 | 同一尺度（受精后周）。Aldinger 9/10/11/12/14/17/18/20/21 wpc 合并为 5 个 bin（§2 `time_bins`），与 Sepp 人的 9w/11w/17w/20w 按下表配对：`9-10w↔9w`、`11-12w↔11w`、`17-18w↔17w`、`20-21w↔20w`；Aldinger 的 `14w` bin 无 Sepp 对应，只用于 Aldinger 内部轨迹的峰判定。Sepp 的 CS18–CS22 与产后 bin 在 Aldinger 中无覆盖 → 不比较。**实测组×bin 核数（9-10w / 11-12w / 14w / 17-18w / 20-21w）**：Purkinje 4,867 / 15,751 / 1,346 / 3,389 / 371；GranuleLineage 153 / 3,040 / 1,295 / 4,061 / 855；Interneuron 427 / 2,270 / 377 / 2,115 / 762；Astro 425 / 2,414 / 1,272 / 3,773 / 1,019；OligoLineage 28 / 227 / 337 / 327 / 67；Immune 57 / 174 / 159 / 401 / 76；DN(描述) 749 / 1,837 / 888 / 371 / 296。供体数每 bin 2 / 4 / 2 / 3 / 2，但 **Purkinje 在 20-21w 仅 1 供体、17-18w 2 供体**（sem 按规则 NA）；OligoLineage 9-10w = 28 核不可行。每核中位基因数 ~900–1,700。 |
| 比较方式 | 同组同 bin 的供体级 pseudobulk CP10K：跨基因 Spearman（平台不同、灵敏度不同，只比秩）；每基因 4-bin 趋势方向一致率；峰 bin 吻合率（±1 bin）。 |

#### C. Castonguay 2026（成人/老年 SPLiT-seq）↔ Sepp（Fig 6d）

| 维度 | 规则 |
| --- | --- |
| 细胞类别 | `author_cell_type`（14 类）→ 9 组（config `external.castonguay.other_groups`）：Purkinje→Purkinje；Granule→GranuleLineage；Bergmann、Astrocytes→Astro；MLI_1、MLI_2、Golgi、PLI→Interneuron；Oligodendrocytes、OPC→OligoLineage；Microglia→Immune；UBC→UBC；Endocytes（内皮）、Pericytes→ 只进 ALL。成人无 DN（取样为小脑皮质）、无 Progenitor。 |
| 时间 | **不对齐**：109 供体全部 63–100 岁，视为单一"老年成人"点；不与 Sepp 的任何 bin 配对，也不进入人时间轴。仅用 normal 供体（93 人；ET 供体排除）。 |
| 比较方式 | 跨基因秩比较：Castonguay 成人 PC（供体级 pseudobulk 均值）vs Sepp 人 PC CS22（人 PC 最大 bin）、vs Sepp 鼠 PC P14（最后可行 bin）与 9wk（n=21，描述性）。回答的问题只有一个：胎儿/鼠成年建立的 PC 基因间排序在老年人 PC 中是否保持。其他成人组（Granule/Bergmann/MLI）同法作补充图。 |

#### D. 外部 bulk（BrainSpan / GTEx，Fig 6a）

- 无细胞类别：只与 Sepp **ALL 组**比较。
- 时间：BrainSpan 样本年龄（pcw / 月 / 岁）按人 4 期归并（embryonic ≤ 8 pcw（BrainSpan 最早样本 8 pcw，仅 1–2 例，只画不检验）、fetal 9–37 pcw、infant 0–3 岁、adult ≥ 18 岁；中间年龄 4–17 岁作第 5 期 `child_adolescent`，Sepp 无对应，只画不检验）；GTEx 全为成人，归 adult 期。比较 = 逐基因跨期 Spearman + 跨基因相关。

#### E. 模型鼠签名（Fig 6c）

- 签名来自模型鼠小脑（多为 bulk 或 PC 富集）的 DEG，用 `ensmusg_id` 映射到本项目鼠 panel；模型鼠取样年龄（如 SCA1 154Q 5 周/12 周）记录在签名表 `model_age_weeks` 列，用于解释（成年期），不参与对齐。

#### F. 对齐决策的留痕

所有映射表都在 `config.yaml`；s00 把每个外部数据集"作者标签 → 本项目组 → 细胞数"的交叉表写到 `logs/s00_label_mapping_<dataset>.tsv`，S1 补充图展示。未映射标签必须在表中显式列为 `unmapped`，不能静默丢弃。

## 5. 第 1–2 层（s01–s03，复用引擎）产物定义

所有表以 core 基因全集为行（缺失基因 `status=missing`），列名前缀约定：`log_<gene>`（bin 内供体均值的 log1p CP10K）、`cp_<gene>`（线性 CP10K）、`pct_<gene>`、`logsem_<gene>`。

| 产物 | 内容 |
| --- | --- |
| `intermediate/{sp}/cube_donor.h5ad` | donor × bin × group pseudobulk；`X`=mean log1p CP10K，layers `mean_cp10k`、`pct_pos` |
| `intermediate/{sp}/cube_bin.tsv` | 全骨架 10 组 × bin；`feasible`、`single_donor`、`n_donors`、`n_cells_total` |
| `intermediate/{sp}/cells_panel.h5ad` | 细胞级 log 面板矩阵（s04/s09 状态剖面用） |
| `tables/{sp}/trend_by_group.tsv` | gene × group：`rho_bin, p_bin, q_bin, rho_donor, p_donor, n_bins, low_power` |
| `tables/{sp}/temporal_patterns.tsv` | gene × group：`peak_bin, peak_bin_order, peak_period, peak_internal, peak_vs_flank, period_ratio, temporal_pattern` |
| `tables/{sp}/specificity_global.tsv` | gene：8 组 `cp10k_/pct_/fraction_`、`tau, max_group, max_fraction, specific_call, n_groups_used` |
| `tables/{sp}/specificity_by_period.tsv` | 同上 × period（可行组 < 4 的期 → tau/fraction NA） |
| `tables/{sp}/tau_by_bin.tsv`、`emergence.tsv` | 每 bin tau 与各组 fraction；首个主组 fraction ≥ 0.6 的 bin |
| `tables/{sp}/peak_by_group.tsv` | gene × group：峰 bin、峰 period、峰 CP10K、可行 bin 数（s07/s08 输入） |
| 图 | 单 panel 文件命名 `{species}_{view}_{group}.pdf`（view ∈ a1 ALL 热图 / b1 组热图 / b2 轨迹 / b3 rho 热图 / b4 窗口基因 / c1 dotplot / c2 分期占比 / d1 emergence）；Fig 组装在 s11 |

**temporal_pattern 标签算法**（已在引擎实现）：对组内可行 bin 的 `log_` 轨迹取 argmax；`peak_internal` = 峰两侧均有可行 bin；`peak_vs_flank` = 峰 − max(左均, 右均)；`period_ratio` = 峰所在 period 的 CP10K bin 均值 / 全轴 bin 均值。标签优先级：`transient_peak`（peak_vs_flank ≥ 0.5 且内部峰）> `window_high`（period_ratio ≥ 2 且内部峰）> `monotonic_up/down`（q < 0.05 且 |rho| ≥ 0.4）> `endpoint_high`（period_ratio ≥ 2 但端点峰）> `flat`。

**低表达伪标签掩膜（P7 裁决①，用户确认 2026-11-05）**：`specific_call = not_expressed` 的基因，其 `window_high/transient_peak` 标签一律在 s02 掩膜为 `NA`（`stats_pattern.mask_lowexpr_labels`；真实鼠草稿 61 行窗口/峰标签中 15 行属此类，Purkinje 组仅有的 2 个窗口标签 BEAN1/TTBK2 均系接近零表达上的假峰）。依据＝NA 政策同构（不可靠证据不给标签）；图上按 pattern_colors.NA 灰色呈现，s07/s08 峰时输入自动排除。掩膜行可识别：`status=='ok'` 且 `temporal_pattern=='NA'`。`endpoint_high` 同类残余不在本次裁决，s06 报告计数留档。

### 5A. `s06_report.py` 物种报告与门清单（P7 提取补注；spec 原无 s06 专节，R1 落点）

输入：`logs/*_report*.tsv`（统一 schema `check/observed/level/note`）、s02 六表、`cube_bin.tsv`、
`cube_donor.h5ad`（仅取供体数/疾病供体数）、`figures/<sp>/` 目录清单。未跑到的层（s07 及以后）
对应小节整节跳过，不产占位。无时间戳（确定性、可 diff）。

产物一 `results/v{ver}/reports/<sp>/species_report.md`，小节固定为：
① 头部（基因集版本/draft 标、基因全集与 cube 命中数、供体数、疾病供体数）；
② 证据结构（组×bin 可行性矩阵、单供体/pooled 清单、NA 清单计数）；
③ 第 1 层趋势（每组 q<0.05 升/降计数表）；
④ 第 2 层（pattern 分布、specific_call 分布与 specific/not_expressed 基因清单、
低表达伪标签掩膜计数（P7 裁决①）与 endpoint_high 残余留档、emergence 计数、峰期×组分布）；
⑤ 锚点与门自动核查（§4 锚点与趋势锚点的 PASS/WARN/FAIL 逐条带值）；
⑥ 本物种相关 §11 未决条目编号（鼠 3/8，人 2/7/9；条文以 §11 为准，代码只列编号防漂移）；
⑦ 图清单（每视图文件茎与件数）；⑧ review 焦点（b2/c1 首页路径）；
⑨ 口径注（NA 政策、b2 的 log1p 尺度注记、draft footer、低表达掩膜待裁决状态）。

产物二 `reviews/gate{N}_<sp>.md`（N：mouse=1、human=2）＝§9 四条清单：
(1)(2)(4) 自动填充（FAIL 计数判级 + WARN/锚点逐条列出），(3) 取自 `s01_report_<sp>` 的 marker 行
（鼠免疫 dropout 豁免注记自动附，§11.1）；每条留人工确认复选框与签名栏。
条目 (4)：锚点定义读 `config species.<sp>.trend_anchors`（2026-11-05 裁决人鼠均设，§9(4)）；
期望方向按 config 动态生成，config 未设 → "—"，设了但日志无行（s02 未跑）→ ✗。
**防覆盖规则**：目标文件已存在时改写 `gate{N}_<sp>.refresh.md`，人工勾选记录不被重跑破坏。

### 5B. `s05_export_tables.py` 基因主表列清单（P7 提取草案；P13 按确认结果实现）

行粒度 = 基因（analysis_protocol L239"每基因×人/鼠×全部统计量×临床注释"）。四文件：
1. `tables/cross/gene_master.tsv`（每基因一行）：身份与临床列 = `genes/ataxia_genes.tsv` 全列；
   每物种 global 统计（`<sp>_tau/_max_group/_max_fraction/_specific_call/_n_groups_used/_emergence_bin`）；
   主组统计（`<sp>_mg_rho_bin/_q_bin/_temporal_pattern/_peak_bin/_peak_period`，主组 = 该物种 max_group）；
   跨物种（s04：`rho_human/rho_mouse/rho_sign_concordant/traj_r`；s09：四保守性指标+综合分+`n_metrics`+证据等级；
   s10：`castonguay_rho/castonguay_dir/aldinger_rho/aldinger_peak_match/brainspan_r`，缺席外部文件 → NA+status）；
   s07 模块（`g_star/module/模块富集 q`）与 s08 onset 统计列。
2. `tables/cross/gene_by_group_long.tsv`（gene×species×group 长表）：rho_bin/p_bin/q_bin/rho_donor/
   temporal_pattern/peak_bin/peak_period/peak_log/peak_cp/status —— S6 全量与 Fig2f 的完整来源。
3. `tables/cross/gene_disease_clinical.tsv`（gene×disease 长表）：`ataxia_gene_disease.tsv` 直连
   （保真度表行来源，§8.1 fidelity_by_disease 的输入）。
4. `tables/cross/gene_master_codebook.tsv`：列名 → 来源脚本/层 → 口径注 对照（硬约束 10：每列有来源）。
NA 政策同硬约束 5；cube 缺失基因 → 各统计列 NA + status=missing。**裁决记录（P7 报告）**：
② 已定（用户确认 2026-11-05）＝宽表仅载主组统计 + `gene_by_group_long.tsv` 全量伴随（本草案方案）。
① 已定（用户确认 2026-11-05）＝not_expressed 基因的 window/peak 伪标签在 s02 掩膜为 NA（§5 掩膜规则；
s05 主表/长表转录掩膜后的值，掩膜行以 `status=='ok'` 且 `temporal_pattern=='NA'` 识别）。
onset 的 CL-024 开关（`onset_include_unverified_draft`）s05 只转录不裁决。

---

## 6. 第 3 层 `s07_category_enrichment.py`

**输入**：两物种 `specificity_global.tsv`、`peak_by_group.tsv`、`temporal_patterns.tsv`；`ataxia_genes.tsv`。

**基因 → 模块指派**（每物种独立）：
1. 主组 `G* = max_group`（global tau argmax，8 主组内）。若 `max_fraction < 0.35` 则 `G* = broad`（广谱）。
2. 期 `P* = peak_period`（`peak_by_group.tsv` 中 G* 行的峰 period）。若 G* 在该物种可行 bin < `min_bins_for_trend` → 模块 = `NA`，该基因不进检验。
3. 模块 = `G* × P*`；另记 `shape = temporal_pattern(G*)`。

**检验**（每物种，机制类 c、模块 m）：观测 = c 中落在 m 的基因数 `O_cm`；期望 = `n_c × n_m / N`；`log2FE = log2((O+0.5)/(E+0.5))`；置换：固定每基因模块，随机打乱机制类标签 10,000 次，双侧 p = 比例(|log2FE_perm| ≥ |log2FE_obs|)；BH 校正族 = 物种内全部 (c, m) 格。`min_genes_per_class`（5）不足的类并入 `other`；基因数 < 8 的类在图中标 `low n`。

**残差类不进入检验网格（CL-027，2026-11-24 用户裁决；plan §0-3/§4.1）**：`other` / `unknown` 是「小类并入」与「无法归入现有类」的**兜底桶，不是一种机制**，「残差类富集于某模块」没有生物学解释。因此机制类 × 模块（Fig3a）与机制类 × shape（Fig3b）的**全部检验格必须排除** `mech_class ∈ {other, unknown}` 的基因；这些基因改做单基因描述（Fig3d 代表轨迹可保留）。落地三处：配置项 `stats.enrichment.exclude_mech_classes`；`atx.stats_enrich.perm_enrichment()` 内置**硬保护**（`labels_a` 含残差类即 `ValueError`，须显式传 `allow_residual=True` 才能绕过且在报告注明）；s07 驱动脚本须在构造标签前过滤并打印被排除的基因清单。v1.0 冻结集下受影响的是 3 个 core 基因：WWOX(SCAR12)、CCDC88C(SCA40)、MME(SCA43)。

**列联（Fig 3b）**：机制类 × shape 六类计数；Fisher 精确（2×k 退化为逐列 2×2 + BH）与同样的置换。

**峰时 ridge（Fig 3c）与统一峰时轴**：每基因 `peak_rank_human_equiv` 按 `peak_source_rule` 确定——G* ∉ `human_censored_groups` 时取人数据 G* 的峰 bin 序号（`peak_source=human`）；G* ∈ {Purkinje, DN, UBC} 时人侧轨迹产后全空、峰时右删失，改取**鼠** G* 的峰 bin 经 `stage_map_mouse_to_human_rank` 映射（`peak_source=mouse`）；若鼠 G* 也不可行（如 DN 的成年段）记 NA。此规则是确定性的，不依赖人侧 temporal_pattern 判定。Kruskal–Wallis + 两两 Mann–Whitney（BH）。**必须同时报告三套结果**：混合轴（主）、仅 `peak_source=human` 的基因、仅鼠峰（所有基因都用鼠峰映射）——三者结论方向一致才可写入正文；不一致时以"仅鼠峰"为准（单一物种、无映射混杂）并在讨论中说明。

**产物**：`tables/cross/module_assignment.tsv`（gene, species, max_group, max_fraction, peak_period, module, shape, peak_rank_human_equiv, peak_source）、`tables/cross/enrichment_class_x_module.tsv`、`enrichment_class_x_shape.tsv`、`peak_rank_by_class.tsv`；图 `fig3a_enrichment_heatmap`、`fig3b_shape_mosaic`、`fig3c_peak_ridge`、`fig3d_representative_traj`。

---

## 7. 第 4 层 `s08_onset_association.py`

**输入**：`module_assignment.tsv`、`ataxia_genes.tsv`（`onset_typical_y, onset_class, mech_class, is_polyQ_set, onset_span_flag`）。

**主检验**：x = `peak_rank_human_equiv`（按 §6 的 `peak_source_rule`；同样报告混合轴 / 仅人 / 仅鼠三套），y = `log1p(onset_typical_y)`；Spearman rho；**分层置换**：在 `mech_class` 内打乱 y（10,000 次），报告 `p_strat`；同时报告未分层 p 作对照。**残差类处理（CL-027）**：与 §6 富集**不同**，分层置换**不排除**残差类（排除会改变样本集、损失功效），而是让 `other`/`unknown` **自成一个 stratum 并在图表标 `low n`**——v1.0 下 `other` 仅 3 个 core 基因（WWOX/CCDC88C/MME），层内只有 3!=6 种排列（`atx/stats_enrich.py` 的 `_MIN_PERM_OBS=3` 为结构下限），该层的分层 p 无解释力；结论必须以未分层 `p_naive` 与「剔除残差层」的敏感性结果为准。排除 `onset_span_flag=true` 的基因做敏感性（主图用空心点标出）。

**类别检验**：`onset_class` 四组 `peak_rank` 的 Kruskal–Wallis；两两 Mann–Whitney + BH；效应量 Cliff's delta。

**位置维度（Fig 4c）**：`onset_class × max_group` 列联，Fisher/置换。

**泳道图（Fig 4d）**：行 = `disease_primary`（按 `peak_rank` 排序）；条带 = `[peak_rank − w_left, peak_rank + w_right]`，其中 w 取该基因在 G* 中 CP10K ≥ 峰值 50% 的连续 bin 范围（**由数据决定，不是示意宽度**）；色 = G*；右注 = 疾病名 + onset_class；`peak_source=mouse` 的行用虚线边框。

**polyQ 特例（Fig 4e）**：对 `is_polyQ_set` 基因，画鼠 ALL `log_` 轨迹（黑）与 Purkinje `fraction_Purkinje`（来自 `tau_by_bin.tsv`，红）；人侧同图但产后 PC 为斜纹。

**敏感性**：`onset_thresholds_grid` 三组阈值、`include_extended_tier`、`cacna1a_as_polyQ`、`drop_glutamatergic_uncertain` 各重跑一次主检验，汇总为 `tables/cross/onset_sensitivity.tsv`。

**产物**：`tables/cross/onset_association.tsv`、`onset_pairwise.tsv`、`onset_sensitivity.tsv`、`swimlane_windows.tsv`；图 fig4a–e。

---

## 8. 第 5–6 层

### 8.1 `s09_conservation_fidelity.py`（第 5 层）

输入：两物种 `trend_by_group.tsv`、`specificity_global.tsv`、`peak_by_group.tsv`、s04 的 `state_profiles.tsv`（ALL 与各组 × 共享 dev_state 的 mean log）。

每基因 × 组（两物种均可行）计算：
- `rho_sign_concordant`：`sign(rho_h) == sign(rho_m)`（两者均 |rho| ≥ 0.2；否则 NA）；
- `traj_r`：人鼠在共享 dev_state 子集（`group_states[g]`）上 z-score 轨迹的 Pearson r（≥ 3 个共享可行状态才算）；
- `max_group_match`：global `max_group` 一致；
- `peak_state_match`：各自峰 dev_state 在共识序中相差 ≤ `peak_state_tolerance`。

基因保真分 = 加权均值（权重见 config，NA 项按剩余权重重归一并记录 `n_metrics`）。疾病保真分 = 该疾病基因的均值。证据等级：`n_metrics ≥ 3` 为 high，2 为 medium，≤ 1 为 low。

产物：`tables/cross/conservation_gene.tsv`、`fidelity_by_disease.tsv`；图 fig5a–e（a/b/c 用 s04 x1–x3，d 用 x4，e 新）。

### 8.2 `s10_external_validation.py`（第 6 层）

**(d) Castonguay 成人 PC 秩基验证**（优先实现，数据已在库）：
1. h5py 打开，读 `obs`（`anndata.io.read_elem`），取 `disease=='normal' & author_cell_type=='Purkinje'`；按供体分组，每供体随机抽 ≤ 300 核（seed）；供体 PC < 30 的剔除。
2. 用 `raw/X` 的 `indptr/indices/data` 做**行切片**（按选中行号排序后分块读取，每块 ≤ 5,000 行），得 CSR；每核 CP10K = counts / 行和 × 1e4；按 ENSG 取 panel 列。
3. 供体级 pseudobulk = 供体内均值；再跨供体均值 → `cast_pc_mean_cp10k`，并记 `n_donors_detected`（供体均值 > 0 的供体数）。
4. 与 Sepp 人 PC CS22 bin 的 `cp_` 均值、鼠 Purkinje 9wk 直接均值（n=21，描述性）、鼠 P14 bin 均值分别做跨基因 Spearman；报告 rho、p、同向率（以各自中位数为界的象限一致率）；偏离对角线 |residual| 最大的 6 基因标注。
5. 同法对 Bergmann/Granule/MLI 组做"成人细胞组 × 基因"dotplot（补充图）。

**(S5) Aldinger 胎儿复核**：`raw/X` 计数 → CP10K；按 `group_map` 分组、`time_bins` 分 bin、供体级 pseudobulk；与 Sepp 人同组同期（9w/11w/17w/20w）的 `cp_` 做：跨基因 Spearman（每组每 bin）、每基因组内趋势方向一致率、峰 bin 吻合率（±1 bin）。

**(a) BrainSpan/GTEx**：需人工下载 `external/brainspan_cerebellum_rpkm.tsv`（BrainSpan RNA-seq Gencode v10，结构 CBC，列：gene_symbol, ensg, age, rpkm）与 GTEx v8 小脑 TPM 的 panel 基因子集。按年龄合并为与人 periods 对应的 4 期；与 Sepp ALL 组的 `cp_` 期均值做逐基因跨期 Spearman 与跨基因相关。文件缺失时脚本跳过并在报告标"待下载"。

**(c) 模型鼠签名回投**：`external/model_mouse_deg_signatures.tsv` 由人工整理（候选：SCA1 Atxn1-154Q PC 转录组 Ingram 2016 Neuron；SCA7 266Q；EA2 tottering；FRDA KIKO/YG8；SCA3 Q84）。签名表每行一个 DEG（`model, gene_mouse, ensmusg_id, direction, log2fc, padj, model_age_weeks, tissue, source_pmid`），取 `padj < 0.05`。分两档实现：**档 1（首轮必做）**：签名 ∩ ataxia panel 的基因在本项目"组×期"模块中的分布，与 panel 全体的模块分布比较（超几何 + 标签置换）——回答"模型鼠扰动的 ataxia 基因是否集中在某模块"；因交集通常只有几个到十几个基因，只作描述与方向性证据。**档 2（可选，证据更强）**：对鼠**全转录组**（20,409 基因）跑一次 pseudobulk（复用 `pseudobulk.py`；115,282 × 20,409 dense float32 ≈ 9.4 GB，可承受），对每个基因做同样的模块指派，得到"全基因 × 模块"背景；再做签名基因的模块富集（超几何 / 秩基 GSEA）。档 2 的结论才能写成"SCA1 模型的转录扰动集中于 Purkinje×出生后模块"。

**(b) ISH**：不自动化；s11 预留 panel 位置并读取 `external/ish/*.png`（人工放置），缺失则绘占位框。

产物：`tables/external/castonguay_pc_rank.tsv`、`aldinger_replication.tsv`、`brainspan_corr.tsv`、`model_signature_enrichment.tsv`；图 fig6a–d。

### 8.3 `s11_assemble_figures.py`

按 `analysis_protocol.md` 的 panel 表把 s03/s04/s07–s10 的单 panel PDF 组装为 Fig 1–6 与 S1–S6（matplotlib `GridSpec` + 读入 PNG，或直接重绘）。每个 panel 左上角字母、每图底部一行"数据口径"注（bin 阈值、基因集版本、NA 规则）。

---

## 9. 执行顺序与人工 review 门

```
./scripts/run_all.sh gate        # s00 + s00b:基因表版本/命中率/锚点/外部文件;失败即停
./scripts/run_all.sh run_mouse   # s01→s02→s03→s06(鼠)        ── 人工 review 门 1:鼠 b2 轨迹与 c1 dotplot
./scripts/run_all.sh run_human   # 同上(人)                     ── 人工 review 门 2
./scripts/run_all.sh run_cross   # s04→s05→s07→s08→s09           ── 人工 review 门 3:fig3a/fig4a 统计
./scripts/run_all.sh run_external# s10                            ── 门 4
./scripts/run_all.sh assemble    # s11
./scripts/run_all.sh sensitivity # 对 sensitivity 各开关循环跑 s01–s02、s07–s08,产物到 results/v1.0/sensitivity/<switch>/
```

每个 review 门的检查清单：(1) `logs/` 无 assert 失败且 WARNING 已逐条看过；(2) 锚点表（§4）全部匹配；(3) marker 校验通过（ITPR1/CA8 → Purkinje，GABRA6/CBLN3 → GranuleLineage，MBP/PLP1 → OligoLineage，AQP4 → Astro，SOX2/MKI67 → Progenitor，PTPRC → Immune；鼠免疫 marker 按 dropout 豁免）；(4) 趋势锚点方向正确。**趋势锚点定义 = `config species.<sp>.trend_anchors`（单一事实来源，硬约束 6；2026-11-05 用户裁决人侧增设并自 s01 常量迁入）**：鼠/人同三条方向 — Itpr1(ITPR1)×Purkinje ρ>0、Mbp(MBP)×OligoLineage ρ>0、Atoh1(ATOH1)×GranuleLineage ρ<0。人侧来源（硬约束 10）＝文献（Gold 2003 Neuron：Itpr1 属 RORα 依赖 PC 成熟程序；Zečević 1998 Dev Brain Res + Ozarkar 2025 Front Neuroanat：MBP 随髓鞘化上升；Haldipur 2022 Annu Rev Neurosci：ATOH1 为 RL/EGL 祖标记、分化下调）＋独立人胎数据 Aldinger 2021 实测（`scripts/dev_check_human_anchors.py` → `logs/dev_human_anchor_check.tsv`，2026-11-05：ITPR1 ρ=+0.643、ATOH1 ρ=−0.790、祖>GN 横向对照 75% 周阳性；MBP 系出生后事件、胎窗不适用 → N/A_WINDOW，方向由文献支撑、门 2 在人 Sepp 数据上正式核验）。人侧可行性（乙）在门 2 用真实人类 cube 核对。

---

## 10. 预期运行资源

| 步骤 | 时间 | 峰值内存 |
| --- | --- | --- |
| s01 人 / 鼠 | ~6 / ~3 min | ~16 / ~8 GB |
| s02–s06 | 各 < 2 min | < 4 GB |
| s07–s09（10,000 置换） | < 10 min | < 4 GB |
| s10 Castonguay（行切片 ~27k 核） | ~5 min | < 6 GB |
| s10 Aldinger | ~2 min | ~4 GB |
| 敏感性全套 | ~1.5 h | — |

---

## 11. 已知未决事项（实现时须显式处理，不能默默假设）

1. ~~鼠 4 基因 dropout 需在 `mouse.h5ad` 重验~~ → 已验（§3.12）：Neurod1/Ptprc/Aif1/Csf1r/Tmem119 豁免，鼠免疫锚点改用 Cx3cr1/P2ry12/Hexb；s00 仍现场重算并与这些数值比对。
2. 人 `interneuron` 在 6mo（17）与 44y（27）不可行；人 `Astro` 在 11w（8）不可行——这些是"中间缺口"，`peak_internal` 判定与轨迹断线要正确处理不连续可行段（引擎已支持，但需在 S1 可行性矩阵中显示）。
3. 鼠 Purkinje P14 恰好 30 核：`min_cells_grid` 的 50/100 档会把它变 NA，使鼠 PC 轨迹止于 P7；敏感性报告要单独指出。**2026-11-24 用户裁决：`bins.min_cells` 维持 30 为主口径**——实测人侧 30↔50 仅差 1 行显著趋势（保留率 99.5%，判级"稳健"），而 100 档人侧保留率 53.0%、鼠侧 71.3%，会砍掉近半结论，故 **50/100 仅作敏感性对照，不得作主口径**；全部 7 个变体的方向一致率均为 100%（无一例升降翻转）。变体产物见 `results/v1.0/sensitivity/`，比较见同目录 `COMPARISON.md`（由 `scripts/a5_sensitivity_compare.py` 生成）。
4. `stage_map_mouse_to_human_rank` 的 6 个锚点来自 Sepp 2024（§4.5 A），其余 6 个为插值；任何用它得到的 p 值须在图注标"鼠峰时按 Sepp 2024 阶段对应映射（部分插值）"，并以"仅鼠峰"结果作稳健性对照。
5. Castonguay 基因符号以 ENSG 匹配；若 panel 某基因的 ENSG 在 Castonguay `var.index` 缺失（如版本差异），记 `status=missing_in_castonguay`，不能回退到符号匹配后静默使用。
6. `external/` 三个人工文件（BrainSpan、GTEx、模型鼠 DEG）在首轮实现时可缺席；s10 对缺失文件只 WARN 并产出占位 panel，不得报错中断。
7. 人成人 4 个 bin（42/44/46/52 岁，各 1 供体）在 `bin_order` 上是 4 个有序点，但生物学上同属成人、无发育先后；Spearman 趋势会把它们当 4 个时间点。敏感性开关 `collapse_adult_bins: true`（把 4 bin 合并为 1 个 adult bin、4 供体）必须实现并在 S3 报告；~~主分析保持 16 bin 以与可行性矩阵一致~~ → **2026-11-24 用户裁决改判：`collapse_adult_bins: true` 为人侧主口径（13 bin），16 bin 降为敏感性对照（变体名 `adult_16bin`）**。依据：四个成人 bin 生物学同期、各仅 1 供体（44y 唯一供体为心衰、46y 唯一供体为急性心梗），当作有序时间点会制造假趋势——实测合并后人侧显著趋势行 256→185（−28%），而保留行**方向一致率 100%**（无一例翻转），说明被去掉的是错误时间排序产生的信号而非真实生物学。连带影响两条：(a) 本节与 §4.2 的**可行性锚点表按 16 bin 记录**，与新主口径不可比 → s01 在合并口径下跳过该比对并写明原因（`bins_default=False`），锚点表待 v1.1 按 13 bin 重记；(b) §4.1「9 组细胞数锚点」按**组**计数、与 bin 切法无关，已改为在合并口径下**照常比对**（s01 新增 `cells_default` 判据，过去随 `default_mode` 一并跳过属过宽）。
8. 鼠 `author_stage` 的 `adult (9 weeks)` 带空格与括号，所有以组名/bin 名构造文件名的地方须做安全化（映射到 `time_labels` 的 `9wk`）。
9. 人 `development_stage` 的 `newborn stage (0-28 days)` 合并了 `author_stage` 的 `newborn`（11,141 核，2 供体）与 `newborn (5 days)`（5,046 核，先心病供体）；bin 级按 `development_stage` 合并为 NB（3 供体），疾病供体在 donor cube 中保留 `disease_flag`。
