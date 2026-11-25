#!/usr/bin/env python3
"""make_verification_checklist.py — 生成决定二·方案B的人工来源核实清单(验收 #2, CL-015/021 销项前置)。

一行 = 一个基因×疾病行(与 ataxia_gene_disease.tsv 同粒度), 预填全部草稿值(*_cur, 冻结记录, 核实者不得改),
右侧为空白 *_v 填写列 + verdict/签名列。批次(a→e)与轻/全深度按 curation_log §4 第 1 条与 CL-023:
  a: 起病值有变更史或 NA 的基因(CL-003/004/005/008) — 直接进 s08 连续统计, 最先;
  b: 旧清单 21 基因(CL-015) — core/extended 资格悬在未核来源上;
  c: 新近发现基因(CL-021, 未被 a 覆盖者) — MIM/PanelApp 状态需坐实;
  d: 其余 effective core;
  e: 其余(extended=full, peripheral=light, CL-023)。
输出: genes/verification/checklist.tsv (幂等重跑覆盖; 填写后的回收文件另存 checklist_filled.tsv)
用法: python -I scripts/make_verification_checklist.py <project_root>
"""
import csv
import os
import sys
from collections import Counter

BATCH_A = {  # CL-003/004/005 起病值变更史 + CL-008 core NA-onset
    'TMEM240', 'ITPR1', 'TTBK2', 'ATXN10',            # CL-003 恢复已核值
    'CACNA1G', 'FAT2', 'SPTBN2', 'ANO10',             # CL-004 有据修正无引注
    'UBA5',                                            # CL-005 冲突置 NA
    'CWF19L1', 'CYP27A1', 'MTTP', 'NPTX1', 'PCNA',    # CL-008 core NA-onset
    'PHYH', 'PUM1', 'SAMD9L', 'TDP2', 'ZFHX3',
}
BATCH_B = {  # CL-015 旧清单 21 基因
    'ATXN8OS', 'RUBCN', 'SYT14', 'WWOX', 'VAMP1', 'COQ2', 'COQ6', 'COQ8B', 'PDSS1', 'PDSS2',
    'PEX1', 'PEX6', 'SURF1', 'RRM2B', 'PLP1', 'PRPS1', 'GBE1', 'ITM2B', 'PGAP1', 'MTPAP', 'SLC16A2',
}
BATCH_C = {'RAB3A', 'THAP11', 'PIK3R5', 'CLCN2', 'PNPT1'}  # CL-021(ZFHX3/NPTX1 已入 a)

CONFLICT7 = ['PRPS1', 'MTPAP', 'GBE1', 'PLP1', 'SURF1', 'RRM2B', 'PGAP1']  # CL-023 sheet3 名单与注释矛盾

SPECIAL_GENE = {
    'TMEM240': 'CL-003: xlsx 曾无出处改为10 → 已恢复 plan §8 已核值 17(综述均值17.4, 范围6–30); 请回填原始出处编号并复核原文',
    'ITPR1': 'CL-003: SCA15/16 行已恢复 40(GeneReviews 均值39.6, 范围20–66); 请回填 NBK 编号并复核; SCA29/Gillespie 行为 expert_estimate, 一并核',
    'TTBK2': 'CL-003: 已恢复 30(NBK1757 各家系9–50, range_midpoint); 请复核',
    'ATXN10': 'CL-003: 已恢复 32(均值32.2, 范围18–45); 请回填出处',
    'CACNA1G': 'CL-004: 25(xlsx 引文献范围13–32)方向正确但原文引注未落表',
    'FAT2': 'CL-004: 60(SCA45 晚发50–73)原文引注未落表',
    'SPTBN2': 'CL-004: SCA5 代表行(25)/SCAR14 行(0)分别复核',
    'ANO10': 'CL-004: 30(min/max=20/50)原文引注未落表',
    'UBA5': 'CL-005: plan 种子 0 vs xlsx 6 冲突 → typical=NA; 以 GeneReviews SCAR24 专篇为裁决点, 本行批次 a 最优先',
    'PNPT1': 'CL-007: SCA25 基因归属本身待核(单家系连锁证据); MDDS 行 core_feature=false; tier 已降 extended',
    'ITM2B': 'CL-013/CL-023: core_feature=false 系初判待核 — 该标注决定 P1 降级是否成立; 若文献推翻, 填 ataxia_core_feature_v=true 即可, 重跑 s00b 自动升回 extended',
    'SLC16A2': 'CL-013/CL-023: 同 ITM2B, P1 降级依赖 false 标注; AHDS 的 ataxia 地位以 GeneReviews 为准',
    'ATXN8OS': 'CL-009: 非编码 RNA, in_panel=false; light 核实(收录证明即可)',
    'RUBCN': 'CL-015: 报道家系少, evidence 定级从严; 注意 Castonguay 中 KIAA0226→RUBCNL 口径差异(plan §6.5), 本行只核文献',
    'SYT14': 'CL-015: 报道家系少, evidence 定级从严',
    'WWOX': 'CL-015: 报道家系少, evidence 定级从严',
    'VAMP1': 'CL-015: 报道家系少, evidence 定级从严; CL-014 机制类争议(channel_synapse vs cytoskeleton_transport_glia)留待双标注轮, 本轮不改 mech_class 除非文献明确',
    'ZFHX3': 'CL-008+CL-021: SCA4 2023 确立, onset=NA 需找数值; MIM 号需在线坐实',
    'NPTX1': 'CL-008+CL-021: SCA50(MIM 620158 via xlsx 备注), onset=NA 需找数值; MIM 需坐实',
    'RAB3A': 'CL-021: SCA52, xlsx 称 PanelApp 绿(MIM 621535)需在线复核',
    'THAP11': 'CL-021: SCA51, OMIM 未正式分配 — 核实基因-疾病关系到何种证据强度',
    'PIK3R5': 'CL-021: AOA3, MIM/证据等级需坐实',
    'CLCN2': 'CL-021:  leukoencephalopathy 伴 ataxia, 证据等级需坐实',
    'TREX1': 'CL-023-P2a: 实测仅 Aldinger var 有(两主数据集均缺)→ peripheral; light 核实(AGS1 行, 表内唯一行; 等位病 RVCL 未入表, 如认为应增行记入 notes, 勿自行加行)',
    'MT-ATP6': 'CL-023-P2c: 核文库 mtDNA 测量效度存疑(s01 行动项); light 核实(NARP/MILS)',
    'MT-ND5': 'CL-023-P2c: 同 MT-ATP6; light 核实(MELAS)',
    'MT-TL1': 'CL-023-P2b: 仅 Sepp human 可测; light 核实(MELAS/MIDD)',
    'MT-TK': 'CL-009/CL-023-P2a: 四数据集均无, in_panel=false; light 核实(MERRF)',
    'MT-TI': 'CL-009/CL-023-P2a: 同 MT-TK; light 核实',
}
SPECIAL_GENE.update({g: 'CL-023: sheet3「综合征型」名单与注释矛盾基因 — 重点复核 ataxia_core_feature; '
                        'true 维持 extended, 推翻为 false 则重跑 s00b 自动降 peripheral, 两者都无需手改表'
                     for g in CONFLICT7})
SPECIAL_GENE['MTPAP'] += '; 缩写 SPAX4 与 MIM 号均未坐实(omim_phenotype=NA), 本轮一并核病名/MIM'

PREFILL = ['batch', 'check_depth', 'row_id', 'hgnc_symbol', 'disease_abbrev', 'disease_name',
           'omim_phenotype', 'tier_curated', 'effective_tier', 'peripheral_rule', 'is_primary_disease',
           'inheritance_cur', 'mutation_class_cur', 'mech_class_cur', 'mech_class_alt_cur',
           'onset_min_y_cur', 'onset_max_y_cur', 'onset_typical_y_cur', 'onset_basis_cur',
           'onset_class_cur', 'course_cur', 'primary_pathology_cur', 'ataxia_core_feature_cur',
           'evidence_cur', 'source_primary_cur', 'source_secondary_cur', 'special_note']
FILLIN = ['source_type', 'source_primary_v', 'source_url_v', 'accessed_date_v', 'evidence_quote_v',
          'onset_min_y_v', 'onset_max_y_v', 'onset_typical_y_v', 'onset_basis_v', 'inheritance_v',
          'ataxia_core_feature_v', 'evidence_grade_v', 'course_v', 'mutation_class_v', 'mech_class_v',
          'primary_pathology_v', 'verdict', 'curator', 'verification_date', 'notes']


def main():
    root = sys.argv[1]
    gd = list(csv.DictReader(open(f'{root}/genes/ataxia_gene_disease.tsv', encoding='utf-8'), delimiter='\t'))
    gt = {r['hgnc_symbol']: r for r in csv.DictReader(open(f'{root}/genes/ataxia_genes.tsv', encoding='utf-8'), delimiter='\t')}

    out = []
    for r in gd:
        g = r['hgnc_symbol']
        eff = gt[g]
        batch = ('a' if g in BATCH_A else 'b' if g in BATCH_B else 'c' if g in BATCH_C
                 else 'd' if eff['tier'] == 'core' else 'e')
        depth = 'light' if eff['tier'] == 'peripheral' else 'full'
        notes = []
        if g in SPECIAL_GENE:
            notes.append(SPECIAL_GENE[g])
        if g in BATCH_B and g not in SPECIAL_GENE:
            notes.append('CL-015: 来源为用户旧清单(非合格 source_primary), 补齐出处后方可留在当前 tier')
        if r['onset_typical_y'] == 'NA' and depth == 'full':
            notes.append('onset=NA(CL-008): 在文献中找数值范围/典型值; 确无数值 → verdict 照填, notes 注明「NA 合理」; 绝不填 0/插值')
        if not r['source_primary'].startswith('unverified'):
            notes.append('已有 plan §8 已核转引: 仅需回填 NBK/MIM 编号并复核原文与当前值一致')
        row = {
            'batch': batch, 'check_depth': depth, 'row_id': f"{g}::{r['disease_abbrev']}",
            'hgnc_symbol': g, 'disease_abbrev': r['disease_abbrev'], 'disease_name': r['disease_name'],
            'omim_phenotype': r['omim_phenotype'], 'tier_curated': r['tier'],
            'effective_tier': eff['tier'], 'peripheral_rule': eff['peripheral_rule'],
            'is_primary_disease': r['is_primary_disease'],
            'inheritance_cur': r['inheritance'], 'mutation_class_cur': r['mutation_class'],
            'mech_class_cur': r['mech_class'], 'mech_class_alt_cur': r['mech_class_alt'],
            'onset_min_y_cur': r['onset_min_y'], 'onset_max_y_cur': r['onset_max_y'],
            'onset_typical_y_cur': r['onset_typical_y'], 'onset_basis_cur': r['onset_basis'],
            'onset_class_cur': r['onset_class'], 'course_cur': r['course'],
            'primary_pathology_cur': r['primary_pathology'],
            'ataxia_core_feature_cur': r['ataxia_core_feature'], 'evidence_cur': r['evidence'],
            'source_primary_cur': r['source_primary'], 'source_secondary_cur': r['source_secondary'],
            'special_note': ' | '.join(notes),
        }
        row.update({k: '' for k in FILLIN})
        out.append(row)

    out.sort(key=lambda x: (x['batch'], x['hgnc_symbol'], x['is_primary_disease'] != 'true', x['disease_abbrev']))
    os.makedirs(f'{root}/genes/verification', exist_ok=True)
    path = f'{root}/genes/verification/checklist.tsv'
    with open(path, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=PREFILL + FILLIN, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(out)

    bc = Counter(r['batch'] for r in out)
    dc = Counter(r['check_depth'] for r in out)
    na = sum(1 for r in out if r['onset_typical_y_cur'] == 'NA' and r['check_depth'] == 'full')
    genes = {b: len({r['hgnc_symbol'] for r in out if r['batch'] == b}) for b in 'abcde'}
    print(f'[checklist] 写入 {path}: {len(out)} 行(190 疾病行全覆盖), {len(gt)} 基因')
    print(f'[checklist] 批次行数: {dict(sorted(bc.items()))}; 批次基因数: {genes}')
    print(f'[checklist] 深度: {dict(dc)}; full 且 onset=NA 的行: {na}')


if __name__ == '__main__':
    main()
