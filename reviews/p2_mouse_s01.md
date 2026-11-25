# P2 核对单 — 鼠 s01 真跑(v0.9.1-draft,2026-11-04)

运行:`run_all.sh run_mouse` 的 s00 门 + `s01_build_cubes.py . --species mouse`。
产物:`results/v0.9.1-draft/intermediate/mouse/`(不入 git,可再生);报告 `logs/s01_report_mouse.tsv`。

## 结果总览

- s00 门:PASS=31 / WARN=5 / FAIL=0(与 P1 基线一致,WARN 处置见 commit e3599f2)
- s01 鼠:**PASS=15 / WARN=1 / FAIL=0**,退出码 0;耗时 6.8s,峰值 RSS 1.93GB(全载 330MB h5ad 后早切列)
- 本次真跑暴露并修复 1 个 bug:`s01` 中 `sp=args.species` 遮蔽 `scipy.sparse as sp`(改别名 `sparse`)

## 锚点核对(全部精确)

| 核对项 | 期望 | 实测 | 判定 |
| --- | --- | --- | --- |
| 形状 | 115282×20409 | 同 | ✓ |
| X 整数计数 | 全整数(前 1M nnz) | 全整数 | ✓ |
| 9 组细胞数 §4.1 | 11 值精确 | 全对(含 RL 裁定) | ✓ |
| 组×bin 可行性矩阵 §4.2 | 9×12 全对 | 全对 | ✓ |
| bin 表骨架 | 10 组×12 bin=120 行 | 120 行,可行 92 | ✓ |
| spot_check 两级聚合重算 | 3/3 | 3/3 PASS(seed=20261013) | ✓ |
| 趋势锚点 Itpr1×PKJ | ρ>0 | +0.915(10 bin 粗算) | ✓ |
| 趋势锚点 Mbp×Oligo | ρ>0 | +0.943(6 bin) | ✓ |
| 趋势锚点 Atoh1×GL | ρ<0 | −0.907(12 bin) | ✓ |
| Purkinje×P14 恰 30 核(§11 未决 3) | feasible=True | 30 核/2 供体,可行 | ✓(敏感性 50/100 档将变 NA,P18 须指出) |
| `adult (9 weeks)`→`9wk` 安全化(§11 未决 8) | 无空格标签 | bin 标签含 9wk,无原始串 | ✓ |
| pooled bin(E10.5/E11.5)sem=NA | 硬断言 | bin_table 内断言通过 | ✓ |

## 基因映射

- 统计 panel 命中 103/109=0.9450(≥0.90 过门);缺失 6:ATN1/CA8/KIAA0586/PDYN/SAMD9L/TCTN2(Sepp 鼠 var 无,实测)
- sanity-only 缺失 1:SOX2(CL-018 已豁免,species_valid=human)
- cube 基因数 = 147 − 7 = **140**

## WARN=1 判读(marker 校验,清单外)

- **PAX6**:期望 Progenitor,实测 argmax=GranuleLineage(3.59 vs 1.26,比值 0.35<0.5→FAIL)。
  判读:非流水线错误——PAX6 高表达于 RL 来源颗粒前体,而 RL 细胞经 2026-11-04 实测裁定归
  GranuleLineage(groups.py docstring,§4.1 锚点复现)。属分组口径的已知生物学后果,留门 1 人工确认。
- 边缘 PASS 留档:HOMER3(argmax=Immune 0.55 vs PKJ 0.29,过 50% 线)、SORCS3(argmax=Oligo 2.47 vs IN 2.40)。
- ABSENT=2:EOMES/LMX1A(UBC marker;UBC 为描述组不在 8 主组分母,无法校验,预期行为)。
- spec §9 清单内 marker(ITPR1/CA8→PKJ 等)FAIL=0;鼠免疫按 dropout 豁免口径,PTPRC 不作 ERROR。

## P2c MT 检出留档(CL-023 冻结行动项)

鼠 14 个 mt- 基因,pct+ 中位 1.21%(最高 mt-Rnr2 16.9%,最低 mt-Nd4l 0.07%);
mt-Tl1 不在 Sepp 鼠 var(作者过滤),MT-TL1 在草稿表中 in_panel=false。
线粒体基因整体低检出为实测事实,供冻结时复核 MT 基因处置。

## 产物清单(draft,不入 git)

| 文件 | 规模 | 内容 |
| --- | --- | --- |
| cube_donor.h5ad | 1.5MB | 226 (donor×bin×group) × 140 基因;X=mean log1p CP10K;layers mean_cp10k/pct_pos |
| cube_bin.tsv | 430KB | 120 行全骨架 × 569 列(log_/cp_/pct_/logsem_<gene> + 元数据) |
| cells_panel.h5ad | 20MB | 115282 × 140 细胞级 log 面板(s04/s09 用) |
