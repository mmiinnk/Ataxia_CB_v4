# Purkinje 细胞跨病因共享靶点研究

本目录含研究方案、四份参考 H5AD 的独立副本、经人工复核的基因–疾病工作簿、数据登记表、首轮脚本和已经算出的结果。可以将整个目录移动到另一台机器。脚本默认从自身位置定位输入，疾病矩阵另行取得。

[完整研究方案](ataxia_network_study_proposal.md) · [疾病与扩展队列登记表](config/cohorts.tsv) · [188 条疾病记录及当前资源状态](outputs/disease_registry.tsv)

## 从这里开始

1. 查看 [四库供体覆盖汇总](outputs/reference_summary.tsv) 和 [逐供体、阶段的 PC 核数](outputs/donor_stage_coverage.tsv)。判断哪些阶段适合建网、哪些只能作描述。Sepp 人数据有 14 位供体至少在一个阶段达到 30 个 PC 核，Aldinger 有 9 位。Castonguay 有 89 位正常供体达到该标准。Sepp 鼠早期合并样本保留在覆盖表，个体推断单列。
2. 使用目录内原始矩阵，从 Sepp 人 PC 构建供体内、分阶段和亚型的 metacell，开始方案的正常模块发现。现成的 [人 PC 聚合样本表](outputs/sepp_human_samples.tsv)、[全基因索引](outputs/sepp_human_features.tsv) 和 `outputs/sepp_human_counts.npz` 支持立即核对表达、计数和样本覆盖。44 个聚合行覆盖全部 45,341 个 PC 核与 27,410 个基因。建网仍须回到核级数据，亚型、metacell 和重复分化批次共享供体时要一起重采样。
3. 并行核对人 SCA1、A-T 的样本和疾病数据可得性。给每个供体或动物填写 [样本对应表](templates/sample_manifest.tsv)，确认组织、技术重复、计数层、PC 数及与其他研究的脑库重叠，再做队列内疾病效应比较。同步盘点本地 SCA3 独立家系、随访和患者来源细胞，并评估人体小脑数据的获取条件。新资源登记到 [队列表](config/cohorts.tsv)。预留的 SCA7 先仅核对技术元数据。
4. 完成正常模块复核与两个发现队列的分析后，复制 [预测冻结表](templates/prediction_lock.tsv) 到带日期的分析目录，填写成员和权重文件校验值、方向、主终点、效应下限、组织、发育阶段及排除规则。保存后再运行 SCA7 预留检验。比较同一基因在两个不同机制背景中的功能救援是实验阶段的主要任务。

## 复现已运行的步骤

在本目录内执行以下命令。本机有两套已验证的环境，第一套包含读取 Excel 所需的 openpyxl，第二套包含单核数据分析依赖。

```bash
/home/victor/.conda/envs/pd_tremor_cl/bin/python scripts/import_inputs.py --verify
/home/victor/.conda/envs/pd_tremor_cl/bin/python scripts/00_registry.py
/home/victor/.conda/envs/omics/bin/python scripts/01_reference.py --pseudobulk sepp_human
```

换机器时可用 Python 3.14 的虚拟环境安装启动脚本使用的包，然后运行同样的脚本。下列依赖用于启动步骤。后续 hdWGCNA、pySCENIC 和对应物种的基序库应在各自兼容环境中独立记录版本。

```bash
python -m venv .venv
.venv/bin/python -m pip install -r requirements-start.txt
.venv/bin/python scripts/import_inputs.py --verify
.venv/bin/python scripts/00_registry.py
.venv/bin/python scripts/01_reference.py --pseudobulk sepp_human
```

`01_reference.py` 无参数时仅更新四库覆盖表。可把 `--pseudobulk` 换成 `sepp_mouse`、`aldinger` 或 `castonguay`，或指定 `--out` 写至新的分析目录。重算只覆盖相应的自动生成文件。个体层人工复核意见应另建文件，以 `row_id_h` 或样本编号关联。若用 Git 分发，H5AD 已列入本目录忽略规则，搬迁时应一并传输 `inputs/reference/`，并用 `--verify` 校验。如果手上只有方案目录的代码，`scripts/import_inputs.py --reference-dir /path/to/four_h5ad --workbook /path/to/merged.xlsx` 可重新导入源文件。

## 输入、计数与适用范围

| 文件或结果 | 用法 |
|---|---|
| [SHA-256 清单](inputs/sha256.tsv) 与 [参考输入配置](config/reference.tsv) | 前者标识四个 H5AD 和工作簿的实际字节版本。后者指定原作者 PC 标签、供体和时间列、计数层与参考集定义。Sepp 主文件保留临床状态标签，按方案另做剔除合并症供体的敏感性分析。 |
| [G3 工作簿](inputs/merged.xlsx) 与 [逐病种登记表](outputs/disease_registry.tsv) | `main` 表的 188 行对应 177 个基因，另有两条排除记录保存在原工作簿。导出保留人工层 `*_h` 和待裁定信息。每行匹配已登记资源或标记“待逐型核查”。这是数据盘点，不能代替个体病例表。 |
| [鼠面板前期报告](inputs/mouse_pilot.md) | 验证细胞定位和发育趋势的历史草稿。已知阳性基因须在新网络中标明。 |
| `donor_stage_coverage.tsv` | 每个有细胞的供体、阶段和疾病状态各一行。真实无 PC 的样本显示 0，低覆盖和 pooled 样本分别标识。未注释细胞计入总核数，汇总表单列其数量。 |
| `*_counts.npz`、`*_samples.tsv`、`*_features.tsv` | SciPy CSR 原始计数和。矩阵行依样本表的零基 `row` 排列，列依特征表的零基 `column` 排列。Sepp 额外按 `dev_state`、`subtype` 分层，并保存每行作者 size factor 之和。稀少的行保留供质控，低覆盖时分开解释。 |
| [样本对应表](templates/sample_manifest.tsv) 与 [预测冻结表](templates/prediction_lock.tsv) | 两者仅含表头，分别用于登记新队列和在留出检验之前固定预测。空白代表待落实的资料。 |

## 取得疾病数据与外部支持

| 用途 | 文件与存放位置 |
|---|---|
| 人 SCA1 发现、鼠病程 | 从 GEO 的 Supplementary files 取人和鼠 H5AD 压缩文件及 Table S1。保存于 `inputs/disease/GSE246183/`。先合并同一个人的技术重复，审查各时间点单只动物的键，再计算 PC 评分。[GSE246183](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246183) |
| 人 A-T 发现 | 从 Broad 的小脑项目获取处理后表达、细胞注释与供体对应。保存于 `inputs/disease/SCP1300/`。单独核对入组 6 位病例、7 位对照中有多少人能支持目标 PC 亚型估计。[SCP1300](https://singlecell.broadinstitute.org/single_cell/study/SCP1300/single-cell-atlas-of-the-human-cerebellum-in-ataxia-telangiectasia) |
| 小鼠 SCA7 留出检验 | 取 RDS 归档、MULTI-seq 条码和 GSM 元数据，存放于 `inputs/disease/GSE269430/`。以“时间点、条码、动物”的映射还原个体，先验证每组有效动物数。候选冻结后查看表达效应。[GSE269430](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269430) |
| SCA3 人体扩展与本地临床对接 | 本地 SCA3 患者资源相对丰富，优先盘点可随访的独立家系。用锁定的模块评估可取得的人体小脑组学，所得方向和读数用于设计患者来源 PC 验证。EGA 尸检小脑 RNA-seq 用于组织层比较，公开的分选 PC 数据用于核对基因可测性。具备匹配的 PC 对照后可估计 PC 内疾病效应。[完整研究方案](ataxia_network_study_proposal.md#病种选择与可用数据) · [EGA 人小脑数据](https://www.ega-archive.org/datasets/EGAD00001009317) · [人分选 PC 原始研究](https://doi.org/10.1002/mds.70465) |
| 其他亚型 | 队列表已登记多病种人脑资料、SCA14 和 ARSACS 组织证据，以及 SCA27B 患者来源细胞。各自标明公开入口、申请条件与证据层次。新增病例矩阵保存原始来源、获取日期、校验值、注释和样本对应表。[队列登记表](config/cohorts.tsv) |

## 后续交付顺序

| 环节 | 输出 | 完成条件 |
|---|---|---|
| 正常网络 | 模块成员与权重、独立胎儿复核、老年 PC 检出和删供体稳定性 | 冻结至多三个通过复核的模块及其评分，记录哪些阶段无法估计 |
| 调控与疾病发现 | TF–靶基因逐边证据表、人 SCA1 与 A-T 的供体散点及效应区间、SCA1 鼠早期轨迹 | 选出一个方向可检验的主节点和至多一个备选，归档冻结表 |
| 独立亚型 | 预留病种效应、节点与模块各自的结论、亚型构成和损伤标志对照 | 按预设标准归为支持、方向相反或精度不足 |
| 转化验证 | 两个不同病因的匹配模型、单一功能主终点、同向干预与健康加干预对照 | 各遗传背景分别检验靶点作用及功能改善，再决定是否进入动物疗效、耐受性和临床开发 |

本包已执行输入核验、G3 导出、四库供体覆盖和人 PC 全转录组聚合。后续网络建模、病例投影、动物及患者研究的执行条件和判据见完整方案。
