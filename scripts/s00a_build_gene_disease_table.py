#!/usr/bin/env python3
"""s00a_build_gene_disease_table.py — 构建 genes/ataxia_gene_disease.tsv（基因×疾病层）。

输入（只读）:
  1. genes/source/gene_classification_table_v20261017.tsv  — xlsx sheet2 快照（179 基因, 基因粒度）
  2. genes/ensembl_mapping_measured.tsv                    — 实测 ENSG/ENSMUSG 映射（extract_ensembl_mapping.py 产出; 可缺省）
输出:
  genes/ataxia_gene_disease.tsv                            — 方案 §3 全列, 多疾病基因按 §5.1 拆行

裁决全部带 CL-xxx 编号, 见 genes/curation_log.md（2026-10-28, v0.9-draft）。
规则来源: gene_curation_plan.md §1–§5; 起病阈值默认 §4.3 (1,13,30), 经 CLI 传入（config.yaml 建立前为引导期口径, CL-020）。
本脚本幂等: 重跑覆盖输出。用法:
  python -I scripts/s00a_build_gene_disease_table.py <project_root> [--onset-thresholds 1,13,30]
"""
import csv
import sys
import os

CURATOR = 'LLM-draft(claude,victor-review-pending)'
CURATION_DATE = '2026-10-28'
XLSX_REF = '共济失调基因分类表(1).xlsx v2026-10-17 sheet2'
UNVERIFIED = 'unverified:pending GeneReviews/OMIM/PanelApp'

# ---------------------------------------------------------------- 裁决数据 ----
# 多疾病基因拆行（方案 §5.1 全部点名基因 + PNPT1/KCNC3, 见 CL-006/CL-007）。
# 未给出的字段继承 xlsx 基因行。typ=None → NA（CL-008）。prim=True → is_primary_disease。
SPLITS = {
    'CACNA1A': [
        dict(d='Spinocerebellar ataxia type 6', a='SCA6', inh='AD', mut='repeat_polyQ',
             typ=45, basis='expert_estimate', course='slowly_progressive', prim=True,
             ref='plan §5.1 代表疾病; §8 值45(种子未核)'),
        dict(d='Episodic ataxia type 2', a='EA2', inh='AD', mut='LoF',
             typ=10, basis='expert_estimate', course='episodic', prim=False,
             ref='plan §5.1 次要行; §8 值10'),
        dict(d='Familial hemiplegic migraine type 1', a='FHM1', inh='AD', mut='missense_GoF_DN',
             typ=None, basis='NA', course='episodic', prim=False, core=False,
             ref='plan §5.1 次要行; 起病/病理待核'),
    ],
    'ITPR1': [
        dict(d='Spinocerebellar ataxia 15/16', a='SCA15/16', inh='AD', mut='LoF',
             typ=40, omin=20, omax=66, basis='reported_mean', course='slowly_progressive', prim=True,
             sp='GeneReviews (via plan §8 已核 2026-10-14)',
             quote='GeneReviews 均值 39.6, 范围 20–66（plan §8 已核记录）',
             note='CL-003: 恢复已核值40（xlsx 曾改35, 无出处）', ref='plan §5.1/§8'),
        dict(d='SCA29 (先天性非进行性共济失调)', a='SCA29', inh='AD', mut='missense_GoF_DN',
             typ=0, basis='expert_estimate', course='nonprogressive', prim=False,
             note='course 待核', ref='plan §8 (SCA29 起病0)'),
        dict(d='Gillespie syndrome', a='Gillespie', inh='AR', mut='missense_GoF_DN',
             typ=0, basis='expert_estimate', course='nonprogressive', prim=False,
             note='plan §8 仅记 AR/起病0; 突变机制(显性负性错义为主)待核', ref='plan §8'),
    ],
    'KCNC3': [
        dict(d='SCA13 (R420H 成人起病型)', a='SCA13', inh='AD', mut='missense_GoF_DN',
             typ=40, basis='expert_estimate', course='progressive', prim=True,
             ref='CL-006: plan §8 两型分行(代表 ~40); 取代 xlsx 文本值 variable'),
        dict(d='SCA13 (R423H 婴儿起病型, 伴智力障碍)', a='SCA13-infantile', inh='AD', mut='missense_GoF_DN',
             typ=0, basis='expert_estimate', course='progressive', prim=False,
             ref='plan §8 两型分行(次要 0)'),
    ],
    'SPTBN2': [
        dict(d='Spinocerebellar ataxia type 5', a='SCA5', inh='AD', mut='missense_GoF_DN',
             typ=25, basis='expert_estimate', course='slowly_progressive', prim=True,
             ref='plan §5.1/§8'),
        dict(d='SCAR14 (婴儿型)', a='SCAR14', inh='AR', mut='LoF',
             typ=0, basis='expert_estimate', course='progressive', prim=False,
             ref='plan §8'),
    ],
    'STUB1': [
        dict(d='Spinocerebellar ataxia type 48', a='SCA48', inh='AD', mut='missense_GoF_DN',
             typ=42, basis='expert_estimate', course='progressive', prim=True,
             note='plan §7: SCA48=STUB1, 示意脚本 DNAJC5 误标已排除', ref='plan §5.1/§8'),
        dict(d='SCAR16', a='SCAR16', inh='AR', mut='LoF',
             typ=20, basis='expert_estimate', course='progressive', prim=False,
             ref='plan §8'),
    ],
    'GRM1': [
        dict(d='SCAR13', a='SCAR13', inh='AR', mut='LoF',
             typ=0, basis='expert_estimate', course='progressive', prim=True,
             ref='plan §5.1 (AR 证据更强为代表)'),
        dict(d='SCA44', a='SCA44', inh='AD', mut='missense_GoF_DN',
             typ=None, ocls='adult', basis='expert_estimate', course='progressive', prim=False,
             ref='plan §8 (adult)'),
    ],
    'FGF14': [
        dict(d='SCA27B (内含子 GAA 重复)', a='SCA27B', inh='AD', mut='repeat_noncoding',
             typ=60, basis='expert_estimate', course='progressive', prim=True,
             note='4-氨基吡啶反应良好(xlsx 备注)', ref='plan §5.1 (患病率高为代表)'),
        dict(d='SCA27A', a='SCA27A', inh='AD', mut='LoF',
             typ=15, basis='expert_estimate', course='progressive', prim=False,
             ref='plan §8'),
    ],
    'AFG3L2': [
        dict(d='Spinocerebellar ataxia type 28', a='SCA28', inh='AD', mut='missense_GoF_DN',
             typ=20, basis='reported_mean', course='progressive', prim=True,
             sp='GeneReviews/文献 (via plan §8 已核 2026-10-14)',
             quote='均值 19.5（plan §8 已核记录）', ref='plan §5.1/§8'),
        dict(d='SPAX5 (线粒体型痉挛性共济失调5型)', a='SPAX5', inh='AR', mut='LoF',
             typ=0, basis='expert_estimate', course='progressive', prim=False,
             ref='plan §8'),
    ],
    'PNPT1': [
        dict(d='Spinocerebellar ataxia type 25', a='SCA25', inh='AD', mut='LoF',
             typ=None, ocls='NA', basis='NA', course='progressive', prim=True,
             ev='limited', tier='extended',
             note='CL-007: xlsx 起病类 childhood 与 SCA25 原始家系成人起病记载矛盾, 起病置 NA 待核; '
                  '单家系连锁证据 → evidence=limited, tier 降 extended; OMIM 608703 为基因 MIM(xlsx 备注), 表型 MIM 待补',
             ref='CL-007'),
        dict(d='PNPT1 相关线粒体 DNA 耗竭综合征 (MDDS/DEE)', a='MDDS', inh='AR', mut='LoF',
             typ=None, ocls='childhood', basis='expert_estimate', course='progressive', prim=False,
             core=False, ev='unverified', tier='extended',
             note='CL-007: ataxia 非核心表现', ref='CL-007'),
    ],
}

# 基因级覆盖（起病恢复/修正/置 NA）, 见 CL-003/004/005。
OVERRIDES = {
    'TMEM240': dict(typ=17, omin=6, omax=30, basis='reported_mean',
                    sp='GeneReviews/综述 (via plan §8 已核 2026-10-14)',
                    quote='ADCA I 综述均值 17.4, 范围 6–30（plan §8 已核记录）',
                    note='CL-003: 恢复已核值17（xlsx 曾改10, 无出处）; onset_class 由程序重算 childhood→juvenile'),
    'TTBK2': dict(typ=30, omin=9, omax=50, basis='range_midpoint',
                  sp='GeneReviews:NBK1757 (via plan §8 已核 2026-10-14)',
                  quote='GeneReviews NBK1757 各家系 9–50 岁, 取范围中点（plan §8 已核记录）',
                  note='CL-003: 恢复已核值30（xlsx 曾改25, 无出处）; onset_class 重算 juvenile→adult'),
    'ATXN10': dict(typ=32, omin=18, omax=45, basis='reported_mean',
                   sp='GeneReviews/文献 (via plan §8 已核 2026-10-14)',
                   quote='均值 32.2, 范围 18–45（plan §8 已核记录）',
                   note='CL-003: 恢复已核值32（xlsx 曾改35, 无出处）'),
    'CACNA1G': dict(typ=25, omin=13, omax=32, basis='expert_estimate',
                    note='CL-004: 保留 xlsx 有据修正(40→25, 文献起病13–32岁); 原文引注待补'),
    'FAT2': dict(typ=60, omin=50, omax=73, basis='expert_estimate',
                 note='CL-004: 保留 xlsx 有据修正(40→60, SCA45 晚发50–73岁); 原文引注待补'),
    'ANO10': dict(typ=30, omin=20, omax=50, basis='expert_estimate',
                  note='CL-004: plan §7 修正方向一致(青少年→典型20–50岁); 取30, 待原文核实'),
    'UBA5': dict(typ=None, ocls='childhood', basis='expert_estimate',
                 note='CL-005: plan 种子0 vs xlsx 6(用户表2例5、8岁)冲突, 均非合格来源 → onset=NA, 类标签暂 childhood'),
}

# 备注型注释（不改值, 只挂 note）, 见 CL-009/010/014/021。
NOTES = {
    'NPTX1': 'CL-021: 新近报道(xlsx: MIM 620158); mech_class 为推断(xlsx 备注); 待 OMIM/PanelApp 复核',
    'RAB3A': 'CL-021: xlsx 称 PanelApp 绿(MIM 621535); mech_class 为推断; 待复核',
    'THAP11': 'CL-019/021: OMIM 未正式分配(xlsx 备注); extended; is_polyQ_set=true(派生层)',
    'ZFHX3': 'CL-021: SCA4 GGC 重复, 2023 长读长测序确立(xlsx 备注); 原文引注待补; 填补 plan 种子缺口',
    'VAMP1': 'CL-014: SNARE 蛋白; 与 GOSR2(归 cytoskeleton_transport_glia)先例一致性待仲裁',
    'SLC16A2': 'CL-013/014: AHDS 以 ID/运动障碍为主, ataxia_core_feature=false(待核); MCT8 膜转运, mech 归属待仲裁',
    'IFIH1': 'CL-014: 胞内 RNA 感受器, mech=other; 备选 proteostasis_RNA 待仲裁',
    'ATXN8OS': 'CL-009: plan §1.2 非编码 RNA → in_panel=false(s00b 强制), 不参与统计; SCA8 证据等级待核',
    'ITM2B': 'CL-013: 脑淀粉样血管病/痴呆为主, ataxia_core_feature=false(待核)',
    'TGM6': 'CL-012: plan §1.2 SCA35 争议较大 → evidence=limited, 不入 core',
    'TRPC3': 'CL-012: xlsx 备注"证据有限" → evidence=limited',
    'CACNB4': 'CL-012: xlsx 备注"证据有限" → evidence=limited',
    'CCDC88C': 'CL-012: xlsx 备注"证据有限" → evidence=limited',
    'RUBCN': 'CL-021: 旧名 KIAA0226(见 gene_aliases.tsv); SCAR15 家系数少, 证据等级待核(ClinGen)',
    'SYT14': 'CL-021: SCAR11 报道家系少, 证据等级待核',
    'WWOX': 'CL-021: SCAR12 证据等级待核',
    'MT-ATP6': 'CL-010: mtDNA 母系; mutation_class=mtDNA_variant 为新增词表项(待回写 plan §4.2)',
    'MT-ND5': 'CL-010: mtDNA 母系; 同 MT-ATP6',
    'MT-TK': 'CL-009/010: mtDNA; 实测四数据集 var 均无此特征 → in_panel=false',
    'MT-TL1': 'CL-010: mtDNA 母系; 同 MT-ATP6',
    'MT-TI': 'CL-009/010: mtDNA; 实测四数据集 var 均无此特征 → in_panel=false',
}
for _mt in ('MT-ATP6', 'MT-ND5', 'MT-TL1'):
    NOTES[_mt] += '; 实测仅 Sepp human var 存在, Aldinger/Castonguay 全缺 → 仅描述性分析'

TIER_CHANGES = {'PNPT1': ('extended', 'CL-007')}
EVIDENCE_MAP = {'TGM6': 'limited', 'TRPC3': 'limited', 'CACNB4': 'limited', 'CCDC88C': 'limited'}
CORE_FALSE = {'SCN2A', 'SCN8A', 'KCNJ10', 'GOSR2', 'PRRT2',   # plan §1.2 / xlsx 临床大类"伴发癫痫(ataxia次要)"
              'PEX1', 'PEX6',                                   # xlsx 备注"ataxia非核心"
              'ITM2B', 'SLC16A2'}                               # CL-013 初判(待核)
OMIM_PHENOTYPE = {'NPTX1': ('620158', 'via xlsx 备注, 待复核'), 'RAB3A': ('621535', 'via xlsx 备注, 待复核')}
INH_NORM = {'AD/AR': 'AD+AR', 'AR/AD': 'AD+AR', 'XL(前突变)': 'XL', 'Mito(母系)': 'mito'}

# sheet3 §四"本次补入 21 个基因"（来源=用户旧清单, 非合格 source_primary, CL-001/015）
OLD_LIST_21 = {'ATXN8OS', 'RUBCN', 'SYT14', 'WWOX', 'VAMP1', 'COQ2', 'COQ6', 'COQ8B', 'PDSS1', 'PDSS2',
               'PEX1', 'PEX6', 'SURF1', 'RRM2B', 'PLP1', 'PRPS1', 'GBE1', 'ITM2B', 'PGAP1', 'MTPAP', 'SLC16A2'}

COLUMNS = ['hgnc_symbol', 'ensg_id', 'mgi_symbol', 'ensmusg_id', 'ortholog_type',
           'disease_name', 'disease_abbrev', 'omim_phenotype', 'inheritance', 'mutation_class',
           'mech_class', 'mech_class_alt', 'onset_min_y', 'onset_max_y', 'onset_typical_y',
           'onset_basis', 'onset_class', 'course', 'primary_pathology', 'ataxia_core_feature',
           'evidence', 'tier', 'is_primary_disease', 'source_primary', 'source_secondary',
           'source_url', 'accessed_date', 'evidence_quote', 'curator', 'curation_date', 'note']


def onset_class_of(y, thr):
    c, j, a = thr
    if y < c:
        return 'congenital'
    if y < j:
        return 'childhood'
    if y < a:
        return 'juvenile'
    return 'adult'


def main():
    root = sys.argv[1]
    thr = (1, 13, 30)
    if '--onset-thresholds' in sys.argv:
        thr = tuple(float(x) for x in sys.argv[sys.argv.index('--onset-thresholds') + 1].split(','))
    snap = os.path.join(root, 'genes/source/gene_classification_table_v20261017.tsv')
    mapf = os.path.join(root, 'genes/ensembl_mapping_measured.tsv')
    outf = os.path.join(root, 'genes/ataxia_gene_disease.tsv')

    genes = list(csv.DictReader(open(snap, encoding='utf-8'), delimiter='\t'))
    mapping = {}
    if os.path.exists(mapf):
        for m in csv.DictReader(open(mapf, encoding='utf-8'), delimiter='\t'):
            mapping[m['query_symbol']] = m

    out_rows = []
    n_split_genes = 0
    for g in genes:
        sym = g['基因符号'].strip()
        notes = [g['备注'].strip()] if g['备注'].strip() else []
        tier = g['层级(tier)'].strip()
        if sym in TIER_CHANGES:
            tier, ref = TIER_CHANGES[sym]
            notes.append(f'tier 调整: {g["层级(tier)"]}→{tier} ({ref})')
        ev = EVIDENCE_MAP.get(sym, 'unverified')
        core_feat = 'false' if sym in CORE_FALSE else 'true'
        inh = INH_NORM.get(g['遗传方式'].strip(), g['遗传方式'].strip())
        if g['遗传方式'].strip() == 'XL(前突变)':
            notes.append('遗传方式规范化 XL(前突变)→XL; 前突变机制见 mutation_class=repeat_noncoding (CL-011)')
        if g['遗传方式'].strip() in ('AD/AR', 'AR/AD'):
            notes.append(f"遗传方式规范化 {g['遗传方式']}→AD+AR (CL-011)")

        # 溯源链（CL-001/015）
        sec = [XLSX_REF]
        if sym in OLD_LIST_21:
            sec.append('用户旧清单11.6.xlsx(via xlsx sheet3 §四, 非合格来源, 待逐条核实)')
        elif '检索补入' in g['备注'] or '补漏' in g['备注']:
            sec.append('xlsx 检索补入(sheet3 §四, 待核)')
        else:
            sec.append('gene_curation_plan.md §8 种子清单')
        source_secondary = ' + '.join(sec)

        m = mapping.get(sym, {})
        base = dict(
            hgnc_symbol=sym,
            ensg_id=m.get('ensg_id', 'NA') or 'NA',
            mgi_symbol=m.get('mgi_symbol', 'NA') or 'NA',
            ensmusg_id=m.get('ensmusg_id', 'NA') or 'NA',
            ortholog_type=m.get('ortholog_type', 'NA') or 'NA',
            mech_class=g['分子机制类'].strip(),
            mech_class_alt=(g['次选机制类'].strip() or 'NA'),
            course=g['病程(course)'].strip(),
            primary_pathology=g['主要病理'].strip(),
            ataxia_core_feature=core_feat,
            evidence=ev,
            tier=tier,
            source_primary=UNVERIFIED,
            source_secondary=source_secondary,
            source_url='NA',
            accessed_date='2026-10-17',
            evidence_quote='NA',
            curator=CURATOR,
            curation_date=CURATION_DATE,
            clinical_group_draft=g['临床大类'].strip(),   # 辅助列, 未注册进 plan 词表 (CL-010)
        )
        if sym in OMIM_PHENOTYPE:
            base['omim_phenotype'] = OMIM_PHENOTYPE[sym][0]
            notes.append(f"omim_phenotype={OMIM_PHENOTYPE[sym][0]} ({OMIM_PHENOTYPE[sym][1]})")
        else:
            base['omim_phenotype'] = 'NA'
        if sym in NOTES:
            notes.append(NOTES[sym])

        ov = OVERRIDES.get(sym, {})

        def make_row(spec=None):
            r = dict(base)
            spec = spec or {}
            r['disease_name'] = (spec.get('d') or g['代表疾病'].replace(' (代表)', '').replace('(代表)', '')).strip()
            r['disease_abbrev'] = (spec.get('a') or g['疾病缩写']).strip()
            r['inheritance'] = spec.get('inh', inh)
            r['mutation_class'] = spec.get('mut', g['突变机制'].strip())
            r['course'] = spec.get('course', base['course'])
            r['ataxia_core_feature'] = {True: 'true', False: 'false'}.get(spec.get('core'), base['ataxia_core_feature']) if 'core' in spec else base['ataxia_core_feature']
            r['evidence'] = spec.get('ev', base['evidence'])
            r['tier'] = spec.get('tier', base['tier'])
            r['is_primary_disease'] = 'true' if (spec.get('prim', True) if spec else True) else 'false'
            r['source_primary'] = spec.get('sp', ov.get('sp', UNVERIFIED))
            r['evidence_quote'] = spec.get('quote', ov.get('quote', 'NA'))
            if r['source_primary'].startswith('GeneReviews'):
                r['accessed_date'] = '2026-10-14'   # plan 已核日期
            # 起病
            typ = spec.get('typ', ov.get('typ', g['典型起病(岁)'].strip())) if (spec or ov) else g['典型起病(岁)'].strip()
            if typ is None or typ == 'NA':
                r['onset_typical_y'] = 'NA'
            else:
                try:
                    r['onset_typical_y'] = str(float(typ)) if not isinstance(typ, (int, float)) else str(typ)
                except ValueError:
                    r['onset_typical_y'] = 'NA'   # 文本值 → NA (CL-008)
            r['onset_min_y'] = str(spec.get('omin', ov.get('omin', 'NA')))
            r['onset_max_y'] = str(spec.get('omax', ov.get('omax', 'NA')))
            r['onset_basis'] = spec.get('basis', ov.get('basis', 'NA')) if r['onset_typical_y'] != 'NA' else \
                spec.get('basis', ov.get('basis', 'expert_estimate'))
            # onset_class: 数值 → 程序派生(§4.3); NA → 沿用标签并标记(CL-008)
            if r['onset_typical_y'] != 'NA':
                derived = onset_class_of(float(r['onset_typical_y']), thr)
                xl_label = g['发病年龄类'].strip()
                if 'ocls' in spec:
                    derived = spec['ocls'] if spec['ocls'] != 'NA' else derived
                if derived != xl_label and 'ocls' not in spec:
                    notes_row = f"onset_class 重算 {xl_label}→{derived}(阈值{thr}, §4.3)"
                    r['_extra_note'] = notes_row
                r['onset_class'] = derived
            else:
                r['onset_class'] = spec.get('ocls', ov.get('ocls', g['发病年龄类'].strip()))
                if r['onset_class'] not in ('NA', ''):
                    r['_extra_note'] = f"CL-008: 文本起病({g['典型起病(岁)'].strip()})→onset_typical_y=NA, 类标签暂沿用, 待核"
                else:
                    r['_extra_note'] = 'onset_class=NA: 记载矛盾或无来源, 待核'
            # note 汇总
            row_notes = list(notes)
            if spec and spec.get('note'):
                row_notes.append(spec['note'])
            if not spec and ov.get('note'):
                row_notes.append(ov['note'])
            if r.get('_extra_note'):
                row_notes.append(r.pop('_extra_note'))
            if spec and spec.get('ref'):
                row_notes.append(f"[{spec['ref']}]")
            # 去掉与基因级重复的 note（拆行基因每条疾病行保留自身裁决）
            r['note'] = ' | '.join(dict.fromkeys(row_notes)) if spec is not None or True else ''
            return r

        if sym in SPLITS:
            n_split_genes += 1
            for sp_ in SPLITS[sym]:
                out_rows.append(make_row(sp_))
        else:
            out_rows.append(make_row(None))

    # 去重 note（基因级 note 在拆行时会在每行重复, 保留——单行审计自足）
    with open(outf, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS + ['clinical_group_draft'], delimiter='\t',
                           extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        for r in out_rows:
            r.pop('_extra_note', None)
            w.writerow(r)

    # ---------------------------------------------------------------- 摘要 ----
    syms = sorted({r['hgnc_symbol'] for r in out_rows})
    tiers = {}
    for r in out_rows:
        if r['is_primary_disease'] == 'true':
            tiers[r['hgnc_symbol']] = r['tier']
    n_core = sum(1 for t in tiers.values() if t == 'core')
    na_onset = [r['hgnc_symbol'] for r in out_rows if r['is_primary_disease'] == 'true' and r['onset_typical_y'] == 'NA']
    na_core = sorted({s for s in na_onset if tiers.get(s) == 'core'})
    print(f'[s00a] 写入 {outf}')
    print(f'[s00a] 疾病行数={len(out_rows)}  基因数={len(syms)}  拆行基因={n_split_genes}')
    print(f'[s00a] tier(基因级): core={n_core}  extended={len(tiers)-n_core}')
    print(f'[s00a] 代表行 onset=NA: {len(set(na_onset))} 基因 (core {len(na_core)}): {na_core}')
    # §7 硬断言: 已知错误基因不得出现在表中/core
    forbidden_core = {'DNAJC5', 'HPCA', 'ZNF592', 'SLC1A6', 'PCP2', 'CALB1', 'GABRA6'}
    bad = forbidden_core & set(tiers)
    assert not bad, f'§7 违规基因出现在表中: {bad}'
    print('[s00a] §7 校验: DNAJC5/HPCA/ZNF592/标志基因均未出现 PASS')


if __name__ == '__main__':
    main()
