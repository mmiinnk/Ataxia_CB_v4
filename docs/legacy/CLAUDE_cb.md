# CLAUDE.md · Ataxia_CB

人/鼠小脑发育 snRNA-seq（Sepp 2024）上的共济失调致病基因时空表达研究。本文件只列硬约束与文档分工；细节一律以下列文档为准，不要凭记忆填写数字或映射。

## 文档分工（先读）

| 文件 | 用途 |
| --- | --- |
| `analysis_protocol.md` | 科学方案（导师汇报版）：研究问题、六层分析、图面 panel 构成、局限性 |
| `implementation_spec.md` | **流水线工程规范**：完整 `config.yaml`、`scripts/atx/` 模块设计、实测锚点与可行性矩阵、数据集对齐规则（§4.5）、s07–s11 算法、review 门、已知未决事项（§11） |
| `gene_curation_plan.md` | 基因集与临床注释策展：纳入/排除、受控词表、表结构、验收标准、种子清单（标"已核"者才可直接用） |
| `figure_reading_guide.md` | 逐 panel 读图规则与统计约定 |

实现任务从 `implementation_spec.md` 开始；策展任务从 `gene_curation_plan.md` §6 开始。

## 硬约束

1. **输入只读且唯一**：只读 `/mnt/data/xyy/Ataxia_CB/{human,mouse,aldinger2021_fetal_snRNAseq,castonguay2026_adult_eQTL}.h5ad`。不读 `/mnt/data/xyy/Nature2024/`，不在数据目录写文件。
2. **不拷贝 Depression_CB 代码**（`/home/victor/Projects/Depression_CB`）。可读其 README/CLAUDE.md 了解已验证的口径与 pandas 3 / scipy 陷阱，但本项目按 `implementation_spec.md` §3 独立实现。
3. **Python**：`/home/victor/.conda/envs/omics/bin/python`（3.14）。无 pyarrow → 禁 parquet（中间产物 h5ad + TSV）；无 openpyxl / pypdf。
4. **数据读取**：backed 列切片会卡死——Sepp/Aldinger 全载；Castonguay（100 万核）用 h5py 对 `raw/X` 按行切片。Aldinger/Castonguay 的 `X` 是 log 归一化值，**计数在 `raw/X`**。Sepp 归一化用作者 `obs.size_factor`。
5. **NA 政策**：组×bin < 30 核 → 表 NA、图斜纹；绝不插值、绝不填 0；单供体/pooled bin 的 sem = NA。
6. **参数单一事实来源**：`config/config.yaml`；脚本不得硬编码阈值、组名、颜色。基因表版本 `gene_set.version` 必须与 `genes/curation_log.md` 冻结版本一致。
7. **统计口径**：8 主组（UBC 为描述组、ALL 为参考）；人 9 wpc 属胎儿期；跨物种只用共享 `dev_state` 对齐，日历年龄不对齐；峰时统一轴按 `peak_source_rule`，并同时报告混合/仅人/仅鼠三套结果。
8. **示意图不是口径**：`scripts/make_sample_figures.py` 与 `figures_mock/` 为合成数据，其中的分期、基因–疾病对应与机制类数目有已知错误（见 `analysis_protocol.md` 样例图说明与 `gene_curation_plan.md` §7），不得照抄。
9. **人工 review 门**：鼠跑完 → 人跑完 → 跨物种/第 3–4 层 → 外部验证，每一门经人工确认后再继续（`implementation_spec.md` §9）。
10. 任何写入文档的数字、映射、临床注释必须有来源：实测（给出脚本/命令）或文献（给出 GeneReviews/OMIM/论文出处）；无法核实的标 `unverified`，不得默默假设。

## 当前状态（2026-10-14）

设计阶段：只有文档、示意图与 `results_ct_stage/` 交叉表；无流水线代码、无基因表、非 git 仓库。下一步：`git init` → 基因集策展 → `scripts/atx/` + s00–s06。
