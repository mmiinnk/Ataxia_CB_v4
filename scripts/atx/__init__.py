"""atx — Ataxia_CB 流水线可复用模块包(implementation_spec.md §3)。

独立实现(硬约束 2:不拷贝 Depression_CB 代码;仅借鉴其已验证口径与陷阱,见各模块 docstring)。
参数唯一事实来源 = config/config.yaml(硬约束 6);本包任何模块不得硬编码阈值/组名/颜色。

模块一览(按依赖序):
  config      Cfg/SpeciesCfg 加载与派生(bin_order/period_of/组层级/版本门)
  sanity      §4 实测锚点字典(组计数/可行性矩阵/dropout/形状) + marker 校验
  genes       panel 装载与 var 映射(ENSG/ENSMUSG 主键, symbol+alias 回退)
  io          三种读取模式(全载 / obs+var / CSR 行切片) + 行和流式计算
  groups      obs 派生列 + 细胞组指派(互斥断言) + 锚点核对
  pseudobulk  归一化 + donor cube + bin 表(全骨架, NA 政策) + spot_check
  stats_trend  Spearman 趋势(bin 级主/donor 级辅, ≤6 bin 精确置换 p)+ BH-FDR
  stats_specificity  type_fraction/tau/specific_call(global+by-period)/tau×bin/emergence
  stats_pattern  temporal_pattern 五类标签+NA / peak_by_group / half_max_window(可跨缺口)
  stats_enrich  模块指派(G*×P*) / 置换富集 / 分层置换 Spearman / KW+MW+Cliff's δ / Fisher
  align / external / plots (后续阶段实现)
"""
__version__ = "0.1.0"
