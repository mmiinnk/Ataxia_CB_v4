#!/usr/bin/env bash
# run_all.sh — 流水线总入口(implementation_spec §9;人工 review 门在目标之间,门内步骤失败即停)。
#
# 用法: bash scripts/run_all.sh <target>
#   gate         s00 输入门(四数据集形状/整数性/命中率/锚点)——任何后续目标的前置
#   run_mouse    鼠:s01→s02→s03→s06(物种报告+门清单)→ 人工 review 门 1
#   run_human    人:s01→s02→s03→s06 → 人工 review 门 2;必须在门 1 通过后跑
#   run_cross    跨物种层: s04(状态剖面/图 x1–x4)→ s07(机制类富集/fig3a–d)→ s09(保真度/fig5e)→ 人工 review 门 3
#                (未接入: s05 主表导出=论文阶段; s08 起病关联=待 CL-024 裁决; s10=external/ 文件缺席)
#   run_external 外部验证(s10 Aldinger/Castonguay/Bulk)→ 门 4     [未实现]
#   assemble     s11 图表组装                                      [未实现]
#   sensitivity  敏感性套件: 逐变体重跑 s01–s02, 产物隔离到 results/v{ver}/sensitivity/<变体>_<物种>/
#                (方案A; config.yaml 不动, 开关经 --sens 进程内覆盖) 用法: sensitivity [变体名]
#   --- v4 合并管线新增目标 ---
#   ingest       v4_d00 疾病数据盘点/统一清单/样本登记预填(只读 /mnt)
#   registry     00_registry:merged.xlsx(*_h 人工层,权威)→ disease_registry.tsv [openpyxl 环境]
#   reference    01_reference:四库供体覆盖(+可选 --pseudobulk) [--data-root 直读 /mnt]
#   network       n00–n02: 正常 PC metacell → 共表达模块+LODO/bootstrap 稳定性 → Aldinger/Castonguay
#                 独立复核 → module_freeze.tsv 冻结(阶段1; n01 含重采样重建, 约 15–25 min)
#   bridge|project|lock   b00/b01 成熟度轴与预测锁定阶段        [未实现(阶段 2)]
#
# 约定:python -I(隔离模式,脚本自行插入 sys.path);产物版本目录 results/v{gene_set.version}/;
# FAIL 退出码非 0 即停,不带 --force 类旁路(硬约束 9:门只能人工确认,脚本不设跳过开关)。
# 权威注释来源:inputs/merged.xlsx(G3 *_h 层);genes/ 表为派生物,冻结 v1.0 前全部 DRAFT。
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${ATX_PYTHON:-/home/victor/.conda/envs/omics/bin/python}"
REG_PY="${REG_PYTHON:-/home/victor/.conda/envs/pd_tremor_cl/bin/python}"
DATA_ROOT="${ATX_DATA_ROOT:-/mnt/data/xyy/Ataxia_CB}"
TARGET="${1:-}"

not_impl() { echo "目标 '$1' 依赖的脚本尚未实现(见 implementation_spec §3);当前已接入: gate/run_mouse/run_human 的 s00–s03+s06。" >&2; exit 2; }

case "$TARGET" in
  gate)
    "$PY" -I "$ROOT/scripts/s00_verify_inputs.py" "$ROOT"
    ;;
  run_mouse)
    "$PY" -I "$ROOT/scripts/s00_verify_inputs.py" "$ROOT"
    "$PY" -I "$ROOT/scripts/s01_build_cubes.py" "$ROOT" --species mouse
    "$PY" -I "$ROOT/scripts/s02_stats.py" "$ROOT" --species mouse
    "$PY" -I "$ROOT/scripts/s03_figures_species.py" "$ROOT" --species mouse
    "$PY" -I "$ROOT/scripts/s06_report.py" "$ROOT" --species mouse
    echo "== 鼠 s01–s06 完成。下一步: 人工 review 门 1(reviews/gate1_mouse.md 勾选后 commit)。"
    ;;
  run_human)
    "$PY" -I "$ROOT/scripts/s00_verify_inputs.py" "$ROOT"
    "$PY" -I "$ROOT/scripts/s01_build_cubes.py" "$ROOT" --species human
    "$PY" -I "$ROOT/scripts/s02_stats.py" "$ROOT" --species human
    "$PY" -I "$ROOT/scripts/s03_figures_species.py" "$ROOT" --species human
    "$PY" -I "$ROOT/scripts/s06_report.py" "$ROOT" --species human
    echo "== 人 s01–s06 完成。人工 review 门 2(reviews/gate2_human.md)之后再进跨物种层(硬约束 9)。"
    ;;
  ingest)
    "$PY" -I "$ROOT/scripts/v4_d00_ingest_disease.py" --inventory --data-root "$DATA_ROOT"
    ;;
  registry)
    "$REG_PY" -I "$ROOT/scripts/00_registry.py"
    ;;
  reference)
    "$PY" -I "$ROOT/scripts/01_reference.py" --data-root "$DATA_ROOT"
    ;;
  reference_pb)
    "$PY" -I "$ROOT/scripts/01_reference.py" --data-root "$DATA_ROOT" --pseudobulk "${2:-sepp_human}"
    ;;
  sensitivity)
    # 方案A(spec §9 敏感性套件): 每个变体产物隔离到 results/v{ver}/sensitivity/<变体>_<物种>/,
    # **不覆盖**正式产物; 开关经 --sens 在进程内覆盖, config/config.yaml 一字不动。
    # 每变体跑 s01→s02(不出图, 与 spec "循环跑 s01–s02" 一致; 需要图时手工加跑 s03 --out-subdir 同值)。
    # 用法: run_all.sh sensitivity [变体名]  —— 不给变体名则跑全部。
    ONLY="${2:-}"
    # 格式: 变体名|物种|仅 s01 的参数|s01 与 s02 共用的参数
    # (--min-cells 只属 s01: 它决定 cube_bin.tsv 的 feasible 标记, s02/s03 从 cube 读可行性,
    #  s02 另用 bins.min_cells_total_for_trend, 与 --min-cells 无关)
    VARIANTS=(
      "adult_16bin|human||--sens collapse_adult_bins=false"
      "exclude_disease|human||--sens exclude_disease_donors=true"
      "adult_16bin+exclude_disease|human||--sens collapse_adult_bins=false --sens exclude_disease_donors=true"
      "min_cells_50|human|--min-cells 50|"
      "min_cells_100|human|--min-cells 100|"
      "min_cells_50|mouse|--min-cells 50|"
      "min_cells_100|mouse|--min-cells 100|"
    )
    n_run=0
    for v in "${VARIANTS[@]}"; do
      IFS='|' read -r name sp s01only shared <<< "$v"
      if [ -n "$ONLY" ] && [ "$ONLY" != "$name" ] && [ "$ONLY" != "${name}_${sp}" ]; then continue; fi
      sub="sensitivity/${name}_${sp}"
      echo "== 敏感性变体 ${name} (${sp}) → results/v1.0/${sub}/"
      mkdir -p "$ROOT/results/v1.0/$sub"
      { echo "variant=$name"; echo "species=$sp"; echo "s01_args=$s01only"; echo "shared_args=$shared";
        echo "baseline=results/v1.0/tables/$sp/ (默认口径)"; echo "run_at=$(date -Is)"; } \
        > "$ROOT/results/v1.0/$sub/VARIANT.txt"
      # shellcheck disable=SC2086
      "$PY" -I "$ROOT/scripts/s01_build_cubes.py" "$ROOT" --species "$sp" --out-subdir "$sub" $s01only $shared
      # shellcheck disable=SC2086
      "$PY" -I "$ROOT/scripts/s02_stats.py"       "$ROOT" --species "$sp" --out-subdir "$sub" $shared
      n_run=$((n_run + 1))
    done
    if [ "$n_run" -eq 0 ]; then echo "无匹配变体: '$ONLY'(可选: collapse_adult / exclude_disease / collapse_adult+exclude_disease / min_cells_50 / min_cells_100)" >&2; exit 2; fi
    echo "== 敏感性套件完成(${n_run} 个变体); 产物在 results/v1.0/sensitivity/<变体>_<物种>/{logs,intermediate,tables}/"
    echo "== 对照基准 = results/v1.0/tables/<物种>/(默认口径); 比较时以 trend_by_group / temporal_patterns / specificity_global 为主。"
    ;;
  run_cross)
    # spec §9: run_cross = s04→s05→s07→s08→s09 → 人工 review 门 3。
    # 本仓库当前接入 s04→s07→s09(用户裁决 2026-11-24 的范围):
    #   · s05_export_tables(基因主表导出)为纯汇总层, 推迟到论文阶段;
    #   · s08_onset_association 被 CL-024(onset 未核实政策)未裁决阻塞, 且 v1.0 有 8 个 core 基因 onset=NA、
    #     onset_basis 仅 104/188 行有值 → 待裁决后再进(不在本次关键路径)。
    "$PY" -I "$ROOT/scripts/s04_figures_crossspecies.py" "$ROOT"
    "$PY" -I "$ROOT/scripts/s07_category_enrichment.py" "$ROOT"
    "$PY" -I "$ROOT/scripts/s09_conservation_fidelity.py" "$ROOT"
    echo "== 跨物种层 s04+s07+s09 完成。下一步: 人工 review 门 3(reviews/gate3_cross.md 勾选后 commit)。"
    echo "   未接入: s05(主表导出, 论文阶段)、s08(起病关联, 待 CL-024 裁决)、s10(外部验证, external/ 文件缺席)。"
    ;;
  network)
    # 阶段1: 正常 PC 网络(n00 metacell → n01 模块发现+LODO/bootstrap 稳定性 → n02 独立复核+冻结)。
    # n01 含供体重采样重建(约 15–25 min), 建议后台跑; 冻结门与口径见 config network.* 与 n02 文档串。
    "$PY" -I "$ROOT/scripts/n00_build_metacells.py" "$ROOT"
    "$PY" -I "$ROOT/scripts/n01_discover_modules.py" "$ROOT"
    "$PY" -I "$ROOT/scripts/n02_replication.py" "$ROOT"
    echo "== 正常 PC 网络 n00–n02 完成: 产物 results/v1.0/network/(module_freeze.tsv 为冻结契约)。"
    echo "   bridge(b00 成熟度轴)/project/lock 未实现(阶段 2); 人工复核见 logs/n0*_report.tsv。"
    ;;
  run_external|assemble|bridge|project|lock)
    not_impl "$TARGET"
    ;;
  *)
    echo "用法: bash scripts/run_all.sh {gate|run_mouse|run_human|run_cross|run_external|assemble|sensitivity|ingest|registry|reference|reference_pb|bridge|network|project|lock}" >&2
    exit 2
    ;;
esac
