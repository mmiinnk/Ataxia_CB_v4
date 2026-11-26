# a8 — M1 相对富集/剂量阈值假说检验（鼠）

gene_set v1.0；共同分母 7 组：Astro, DN, GranuleLineage, Immune, Interneuron, Progenitor, Purkinje。
PC 占比在两期共有组上重算（去除 postnatal 多 OligoLineage 的分母混杂）；ΔPC=postnatal−embryonic。
另报 postnatal 8 组 fraction（期内同分母，跨基因可比）。PC_loss 来自临床注释；onset 剔除 NA（CL-024）。

## 检验结果

| 检验 | 指标 | n | 统计量 | p | 效应 |
|---|---|---|---|---|---|
| PC_loss vs 非 | delta_pc_c7 | 17/90 | MW | 0.759 | δ=-0.051 |
| PC_loss vs 非 | pcfrac_post_8g | 17/90 | MW | 0.780 | δ=0.045 |
| onset 分层 Spearman | delta_pc_c7 | 92 | ρ=0.186 | p_strat=0.023 | p_naive=0.076 |
| onset 分层 Spearman | pcfrac_post_8g | 94 | ρ=0.187 | p_strat=0.025 | p_naive=0.072 |

## 代表基因（共同分母 PC 占比）

| 基因 | 机制类 | PC_loss | embryonic | postnatal | ΔPC | postnatal 8组 |
|---|---|---|---|---|---|---|
| ATXN1 | polyQ_repeat | 是 | 0.146 | 0.099 | -0.047 | 0.090 |
| ATXN3 | polyQ_repeat | 否 | 0.125 | 0.123 | -0.002 | 0.104 |
| ATXN7 | polyQ_repeat | 否 | 0.111 | 0.144 | 0.034 | 0.131 |
| CACNA1A | polyQ_repeat | 是 | 0.158 | 0.316 | 0.158 | 0.300 |
| ITPR1 | channel_synapse | 否 | 0.047 | 0.911 | 0.864 | 0.898 |
| TBP | polyQ_repeat | 是 | 0.149 | 0.087 | -0.063 | 0.079 |

## 判读边界

- 关联非因果；属 M1 假说检验，须疾病线验证（核心假说文档 §4）。
- 仅鼠（人 PC 产前-only）；PC 占比为相对量，受其他组表达影响。
- onset 关联为初步（CL-024 未核实，NA 剔除）；残差类 other 自成 stratum（low n）。
- ΔPC 用共同 7 组分母校正；fraction_8group 严格版 embryonic 全 NA 不可用。
