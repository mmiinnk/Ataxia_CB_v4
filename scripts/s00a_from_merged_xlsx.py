#!/usr/bin/env python3
"""s00a_from_merged_xlsx.py — v1.0 冻结接线: 自 inputs/merged.xlsx main 派生基因×疾病表。

用户裁决(2026-11-23): 基因集按 merged.xlsx 冻结。

merged.xlsx main 的真实分层(73 列, 列名唯一; 本脚本一律按列名定位, 不硬编码位置):
  0-10   行标识/策展层: batch, check_depth, row_id, hgnc_symbol, disease_abbrev,
         disease_name, omim_phenotype, tier_curated, effective_tier, peripheral_rule,
         is_primary_disease
  11-25  *_cur 机器/草稿层(含 source_secondary_cur)
  26-46  AI 核验层(*_v + special_note/source_type/verdict/curator/verification_date/notes)
         —— curator='claude-code(AI-assisted · 待人工复核)', 人工 20% 抽查未做
         (audit_report_2026-11-03 §0), 故本次冻结**不**作为权威值层
  47-72  G3/YYB 人工复核层(*_h + g3_verdict/g3_curator/g3_date/merge_notes)——权威

取值策略(逐字段, 前者非 NA 即取):
  G3 人工 *_h > *_cur 机器层;  AI *_v 仅用于 merged.xlsx 别无来源的溯源字段
  (source_primary/source_url/accessed_date/evidence_quote/curator/date), 且 G3 优先。
  onset_basis 只取枚举列 onset_basis_enum_h > onset_basis_cur;
  **绝不**取 onset_basis_h(54, 自由文本引文)或 onset_basis_v(35, AI 自由文本)。
  onset_class/ataxia_core_feature/evidence 无 G3 列 → 取 *_cur:
    · core_grade_pending_h 全 188 行 = 'YYB未覆盖·待G3终判' → G3 未对 core/acf 终判,
      故 AI 层 ataxia_core_feature_v 的翻转(CL-025: COQ6/OPHN1 等)与随之而来的
      CL-023 自动重裁**全部顺延 v1.1**(用户已裁决 CL-024/CL-025 顺延);
    · evidence_grade_v 分布可疑(audit §3: strong162/limited28, 零 definitive/moderate)→ 顺延。
  clinical_group_draft 不在 merged.xlsx → 自冻结前疾病表(.bak_premerge 优先)携带。
  tier=tier_curated(人工策展); effective_tier/peripheral_rule 不用工作簿值, 交 s00b CL-023 重裁。

用法: <openpyxl python> -I scripts/s00a_from_merged_xlsx.py <root>
"""
import csv, os, sys, datetime
from collections import Counter
import openpyxl

COLS = ['hgnc_symbol', 'ensg_id', 'mgi_symbol', 'ensmusg_id', 'ortholog_type',
        'disease_name', 'disease_abbrev', 'omim_phenotype', 'inheritance', 'mutation_class',
        'mech_class', 'mech_class_alt', 'onset_min_y', 'onset_max_y', 'onset_typical_y',
        'onset_basis', 'onset_class', 'course', 'primary_pathology', 'ataxia_core_feature',
        'evidence', 'tier', 'is_primary_disease', 'source_primary', 'source_secondary',
        'source_url', 'accessed_date', 'evidence_quote', 'curator', 'curation_date', 'note']

# 字段 → 来源列名(按优先级)。列名缺失即 AssertionError(防工作簿改版静默错值)。
SRC = {
    'disease_name':         ['disease_name_h', 'disease_name'],
    'disease_abbrev':       ['disease_abbrev_h', 'disease_abbrev'],
    'omim_phenotype':       ['omim_phenotype'],
    'inheritance':          ['inheritance_h', 'inheritance_cur'],
    'mutation_class':       ['mutation_class_h', 'mutation_class_cur'],
    'mech_class':           ['mech_class_h', 'mech_class_cur'],
    'mech_class_alt':       ['mech_class_alt_h', 'mech_class_alt_cur'],
    'onset_min_y':          ['onset_min_y_h', 'onset_min_y_cur'],
    'onset_max_y':          ['onset_max_y_h', 'onset_max_y_cur'],
    'onset_typical_y':      ['onset_typical_y_h', 'onset_typical_y_cur'],
    'onset_basis':          ['onset_basis_enum_h', 'onset_basis_cur'],
    'onset_class':          ['onset_class_cur'],
    'course':               ['course_h', 'course_cur'],
    'primary_pathology':    ['primary_pathology_h', 'primary_pathology_cur'],
    'ataxia_core_feature':  ['ataxia_core_feature_cur'],
    'evidence':             ['evidence_cur'],
    'tier':                 ['tier_curated'],
    'is_primary_disease':   ['is_primary_disease'],
    'source_primary':       ['source_primary_cur', 'source_primary_v'],
    'source_secondary':     ['source_secondary_cur'],
    'source_url':           ['source_urls_h', 'source_url_v'],
    'accessed_date':        ['g3_date', 'accessed_date_v'],
    'evidence_quote':       ['evidence_mech_h', 'evidence_quote_v'],
    'curator':              ['g3_curator', 'curator'],
    'curation_date':        ['g3_date', 'verification_date'],
    'note':                 ['notes', 'special_note'],
}
# 溯源用附加列(写入 note, 不参与统计)
TRACE = ['review_flag_h', 'onset_derivation_note', 'g3_verdict', 'quote_check_status_h']
INTERIM = ('CL-023/CL-025 interim: G3 core/acf 未终判(core_grade_pending_h 全表待判) → '
           'acf/evidence/tier 沿用 *_cur, AI *_v 翻转与重裁顺延 v1.1')


def norm(v):
    if v is None:
        return 'NA'
    if isinstance(v, bool):
        return 'true' if v else 'false'
    s = str(v).strip()
    return s or 'NA'


def main():
    root = sys.argv[1]
    wb = openpyxl.load_workbook(os.path.join(root, 'inputs/merged.xlsx'),
                                read_only=True, data_only=True)
    ws = wb['main']
    it = ws.iter_rows(values_only=True)
    hdr = [norm(h) for h in next(it)]
    assert len(set(hdr)) == len(hdr), '工作簿出现重复列名, 按名定位失效, 请人工核对'
    ix = {h: i for i, h in enumerate(hdr)}
    need = ({c for cs in SRC.values() for c in cs} | set(TRACE)
            | {'hgnc_symbol', 'disease_abbrev_h', 'core_grade_pending_h',
               'ataxia_core_feature_v', 'evidence_grade_v', 'evidence_path_h', 'merge_notes'})
    missing = sorted(need - set(ix))
    assert not missing, f'merged.xlsx 缺列(schema 变动?): {missing}'

    rows = [r for r in it if not all(v is None for v in r)]
    mp = {m['query_symbol']: m for m in csv.DictReader(
        open(os.path.join(root, 'genes/ensembl_mapping_measured.tsv'),
             encoding='utf-8'), delimiter='\t')}
    gdf = os.path.join(root, 'genes/ataxia_gene_disease.tsv')
    today = datetime.date.today().isoformat()
    bak = gdf + '.bak_premerge_' + today
    # 携带源优先用冻结前备份(纯净旧表), 避免把本次派生值当旧值回灌
    carry_src = bak if os.path.exists(bak) else gdf
    carry = {(r['hgnc_symbol'], r['disease_abbrev']): r for r in csv.DictReader(
        open(carry_src, encoding='utf-8'), delimiter='\t')}

    def pick(r, names):
        for n in names:
            v = norm(r[ix[n]])
            if v != 'NA':
                return v, n
        return 'NA', None

    out, prov = [], Counter()
    pending_core = 0
    for r in rows:
        sym = norm(r[ix['hgnc_symbol']])
        m = mp.get(sym, {})
        old = carry.get((sym, norm(r[ix['disease_abbrev_h']])), {})
        rec = {'hgnc_symbol': sym}
        for f in ('ensg_id', 'mgi_symbol', 'ensmusg_id', 'ortholog_type'):
            rec[f] = m.get(f, 'NA')
        for f, names in SRC.items():
            v, used = pick(r, names)
            rec[f] = v
            prov[(f, used or 'NA')] += 1
        rec['is_primary_disease'] = rec['is_primary_disease'].lower()
        for f in ('clinical_group_draft', 'source_secondary', 'evidence_quote', 'source_url'):
            if rec.get(f, 'NA') == 'NA':
                rec[f] = old.get(f, 'NA')
        # note: 基线 + G3 溯源列 + 冻结戳
        bits = [] if rec['note'] == 'NA' else [rec['note']]
        for t in TRACE:
            v = norm(r[ix[t]])
            if v != 'NA':
                bits.append(f'{t}={v}')
        path_q = norm(r[ix['evidence_path_h']])
        if path_q != 'NA':
            bits.append('evidence_path_h=' + path_q[:120])
        if norm(r[ix['core_grade_pending_h']]) != 'NA':
            pending_core += 1
            bits.append(INTERIM)
        bits.append('v1.0 freeze: merged.xlsx G3 *_h authoritative ' + today)
        rec['note'] = ' | '.join(bits)
        rec['accessed_date'] = rec['accessed_date'] if rec['accessed_date'] != 'NA' \
            else rec['curation_date']
        out.append(rec)

    if os.path.exists(gdf) and not os.path.exists(bak):
        os.replace(gdf, bak)
    with open(gdf, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS + ['clinical_group_draft'],
                           delimiter='\t', extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        for rec in out:
            w.writerow(rec)

    print('[s00a_xlsx] rows', len(out), 'genes', len({r['hgnc_symbol'] for r in out}),
          'tier_curated', dict(Counter(r['tier'] for r in out)))
    print('[s00a_xlsx] core/acf 待 G3 终判行(core_grade_pending_h 非空):', pending_core, '/', len(out))
    print('[s00a_xlsx] 取值溯源(字段←列: 行数):')
    for (f, used), n in sorted(prov.items()):
        print(f'    {f:20s} ← {used:24s} {n}')
    print('[s00a_xlsx] wrote', gdf, '| carry from', os.path.basename(carry_src),
          '| backup', bak if os.path.exists(bak) else '(已存在, 未覆盖)')


if __name__ == '__main__':
    main()
