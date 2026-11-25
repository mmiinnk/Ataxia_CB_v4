#!/usr/bin/env python3
"""apply_verification_results.py — 把来源核实轮(决定二·方案B)结果回填 genes/ataxia_gene_disease.tsv。

输入: checklist_filled.tsv(核实者交回件) + genes/verification/checklist.tsv(空白模板, 防篡改基准)
      + genes/ataxia_gene_disease.tsv(回填目标) + genes/ataxia_genes.tsv(s00b 预演用实测旗标/策展 tier)。
输出: 原子重写 ataxia_gene_disease.tsv(先备份 .bak_preapply_<date>), 并写逐行处置报告
      genes/verification/apply_report_<date>.tsv。回填后需重跑 s00b 重裁 effective tier(CL-023)。

用法:
  python -I scripts/apply_verification_results.py <project_root> <checklist_filled.tsv>
         [--check-only] [--batch a,b] [--onset-unverified-policy na|keep-flagged]
         [--onset-thresholds 1,13,30] [--report PATH]

设计: 可反复运行(分批交回场景)——未完成/不合格的行只进报告, 不触碰表格; 合格的行才回填。
三道门(任一触发即中止或行级拦截, 绝不静默写入):
  G1 行级规范门: 词表/格式/verdict 一致性/来源五件套/引文≤200字。可机械归一化的偏差
     (source_primary 尾注、basis 自由文本含唯一词表关键词、XLR→XL、'AD,AR'→AD+AR、
     course 空格→下划线)默认归一并记录(--no-normalize 关闭); 归一不了的行 reject。
     mech_class_v 词表外值(如 gain_of_function)一律**丢弃不回填**(README §4 规则5: 本轮不裁
     mech_class, CL-014 双标注轮另议), 建议值记入该行 note。
  G2 onset 政策门: cur 有数值但核实后无数值(留空)的行 =「草稿值在合格来源不可证」。
     存在此类行且非 --check-only 时必须显式给 --onset-unverified-policy, 否则 ABORT:
       na           → 置 NA(README 严格口径; 草稿值仅存于 checklist 审计记录)
       keep-flagged → 保留草稿值 + note 标记 onset_draft_unverified(敏感性分析开关, 待 plan §3 注册)
  G3 s00b 预演门: 用回填后的 ataxia_core_feature 预演 CL-023 peripheral 规则; 若规则命中
     策展 tier=core 的基因(如 COQ6/OPHN1 翻转为 false), s00b 将 AssertionError —— 本脚本
     提前行级拦截该基因的全部 acf 翻转行(BLOCKED), 交人工裁决(推翻翻转, 或走 CL 条目降策展
     tier)。其余基因的 P1 降级/升回 extended 属 CL-023 自动重裁, 放行并打印预测构成。
verdict 语义: confirmed=全列一致(onset 照抄或留空均视为确认; 任何语义差异→G1 reject 要求改标
  corrected); corrected=应用 *_v; unverifiable=表不动(保持 unverified), 仅计数。
"""
import csv
import datetime
import os
import re
import sys
from collections import Counter, defaultdict

PREFILL = ['batch', 'check_depth', 'row_id', 'hgnc_symbol', 'disease_abbrev', 'disease_name',
           'omim_phenotype', 'tier_curated', 'effective_tier', 'peripheral_rule', 'is_primary_disease',
           'inheritance_cur', 'mutation_class_cur', 'mech_class_cur', 'mech_class_alt_cur',
           'onset_min_y_cur', 'onset_max_y_cur', 'onset_typical_y_cur', 'onset_basis_cur',
           'onset_class_cur', 'course_cur', 'primary_pathology_cur', 'ataxia_core_feature_cur',
           'evidence_cur', 'source_primary_cur', 'source_secondary_cur', 'special_note']
SPECIAL9 = {'PRPS1', 'MTPAP', 'GBE1', 'PLP1', 'SURF1', 'RRM2B', 'PGAP1', 'ITM2B', 'SLC16A2'}
VOC = {
    'source_type': {'GeneReviews', 'OMIM', 'PanelApp'},
    'verdict': {'confirmed', 'corrected', 'unverifiable'},
    'evidence_grade': {'definitive', 'strong', 'moderate', 'limited', 'disputed'},
    'inheritance': {'AD', 'AR', 'XL', 'XLD', 'mito', 'AD+AR'},
    'basis': {'reported_median', 'reported_mean', 'range_midpoint', 'expert_estimate'},
    'course': {'nonprogressive', 'slowly_progressive', 'progressive', 'episodic'},
    'mech_class': {'polyQ_repeat', 'channel_synapse', 'mito_metabolic', 'dna_repair',
                   'proteostasis_RNA', 'cilia_congenital_dev', 'cytoskeleton_transport_glia', 'other'},
    'bool': {'true', 'false'},
}
ISO = re.compile(r'^\d{4}-\d{2}-\d{2}$')
NUM = re.compile(r'^\d+(\.\d+)?$')
SRC_RE = [(re.compile(r'GeneReviews:\s*(NBK\d+)'), 'GeneReviews:{}'),
          (re.compile(r'OMIM:\s*(\d+)'), 'OMIM:{}'),
          (re.compile(r'^(PanelApp:.+)$'), '{}')]
BASIS_KEYS = [('reported_median', 'reported_median'), ('median', 'reported_median'),
              ('reported_mean', 'reported_mean'), ('average', 'reported_mean'), ('mean', 'reported_mean'),
              ('range_midpoint', 'range_midpoint'), ('midpoint', 'range_midpoint'), ('range', 'range_midpoint'),
              ('expert_estimate', 'expert_estimate'), ('expert', 'expert_estimate')]
CONGENITAL_RE = re.compile(r'congenital|at birth|from birth|since birth|\bbirth\b|neonatal|perinatal|prenatal|in utero', re.I)
# *_v 列 → 疾病表列(仅当值与 cur 不同才写)
PAIRS = [('inheritance_v', 'inheritance'), ('course_v', 'course'),
         ('mutation_class_v', 'mutation_class'), ('primary_pathology_v', 'primary_pathology'),
         ('ataxia_core_feature_v', 'ataxia_core_feature')]


def onset_class_of(y, thr):
    c, j, a = thr
    return 'congenital' if y < c else 'childhood' if y < j else 'juvenile' if y < a else 'adult'


def peripheral_rule_of(sym, acf_values, m):
    """与 s00b 同构的 CL-023 规则; acf_values = 该基因各行回填后的 ataxia_core_feature 可迭代。"""
    if not any(v == 'true' for v in acf_values):
        return 'P1'
    if m.get('in_human_sepp') == 'false' and m.get('in_mouse_sepp') == 'false':
        return 'P2a'
    if (m.get('in_human_sepp') == 'true') + (m.get('in_mouse_sepp') == 'true') <= 1 \
            and m.get('in_aldinger') == 'false' and m.get('in_castonguay') == 'false':
        return 'P2b'
    if sym.startswith('MT-'):
        return 'P2c'
    return ''


def num_eq(a, b):
    a, b = a.strip(), b.strip()
    if NUM.match(a) and NUM.match(b):
        return abs(float(a) - float(b)) < 1e-9
    return a == b


def norm_source(v):
    v = v.strip()
    for rx, fmt in SRC_RE:
        m = rx.search(v)
        if m:
            canon = fmt.format(m.group(1))
            tail = (v[:m.start()] + v[m.end():]).strip(' ,;:()')
            return canon, tail
    return None, v


def norm_basis(v, typ_v=''):
    v = v.strip()
    if v in VOC['basis']:
        return v, ''
    hits = {canon for key, canon in BASIS_KEYS if key in v.lower()}
    if len(hits) == 1:
        return hits.pop(), f'basis 自由文本归一: "{v[:60]}"'
    # 规范缺口兜底: congenital/at-birth 措辞支持的 onset=0 行, 四值词表无天然对应
    # (README §3: 0 须原文明确 congenital/at birth) → 按 plan 定义归 expert_estimate, 记录在案
    if typ_v and NUM.match(typ_v) and float(typ_v) == 0 and CONGENITAL_RE.search(v):
        return 'expert_estimate', f'congenital 措辞 0 值行 basis→expert_estimate: "{v[:60]}"'
    return None, f'basis 无法归一(命中 {len(hits)} 个词表关键词): "{v[:80]}"'


def norm_inh(v):
    v = v.strip()
    if v in VOC['inheritance']:
        return v, ''
    if v in ('XLR', 'X-linked recessive'):
        return 'XL', 'XLR→XL'
    toks = [t.strip() for t in re.split(r'[,;/]', v)]
    if set(toks[:2]) == {'AD', 'AR'}:
        return 'AD+AR', f'"{v[:40]}"→AD+AR'
    if toks[0] in VOC['inheritance']:
        return toks[0], f'取首段 {toks[0]}, 限定语移 note: "{v[len(toks[0]):][:60]}"'
    return None, f'inheritance 无法归一: "{v[:60]}"'


def norm_course(v):
    w = v.strip().replace(' ', '_').lower()
    return (w, '') if w in VOC['course'] else (None, f'course 非法: "{v[:40]}"')


class Report:
    FIELDS = ['row_id', 'batch', 'status', 'reasons', 'normalized', 'applied_fields']

    def __init__(self):
        self.rows = []

    def add(self, rid, batch, status, reasons=(), norm=(), applied=()):
        self.rows.append({'row_id': rid, 'batch': batch, 'status': status,
                          'reasons': '; '.join(reasons), 'normalized': '; '.join(norm),
                          'applied_fields': '; '.join(applied)})

    def drop_staged(self, rid):
        self.rows = [x for x in self.rows if not (x['row_id'] == rid and x['status'] == 'staged')]

    def write(self, path):
        with open(path, 'w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=self.FIELDS, delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(self.rows)


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print(__doc__)
        sys.exit(2)
    root, filled_path = args[0], args[1]
    check_only = '--check-only' in args
    normalize = '--no-normalize' not in args
    batches = set(args[args.index('--batch') + 1].split(',')) if '--batch' in args else None
    policy = None
    if '--onset-unverified-policy' in args:
        policy = args[args.index('--onset-unverified-policy') + 1]
        assert policy in ('na', 'keep-flagged'), f'非法政策: {policy}'
    thr = (1.0, 13.0, 30.0)
    if '--onset-thresholds' in args:
        thr = tuple(float(x) for x in args[args.index('--onset-thresholds') + 1].split(','))
    today = datetime.date.today().isoformat()

    tpl = {r['row_id']: r for r in csv.DictReader(
        open(f'{root}/genes/verification/checklist.tsv', encoding='utf-8'), delimiter='\t')}
    fil = list(csv.DictReader(open(filled_path, encoding='utf-8'), delimiter='\t'))
    gd_path = f'{root}/genes/ataxia_gene_disease.tsv'
    gd = list(csv.DictReader(open(gd_path, encoding='utf-8'), delimiter='\t'))
    gt = {r['hgnc_symbol']: r for r in csv.DictReader(
        open(f'{root}/genes/ataxia_genes.tsv', encoding='utf-8'), delimiter='\t')}
    gd_key = {(r['hgnc_symbol'], r['disease_abbrev']): r for r in gd}
    rows_of_gene = defaultdict(list)
    for r in gd:
        rows_of_gene[r['hgnc_symbol']].append(r)
    rep = Report()

    # ---- 0. 结构完整性(硬前提) ----
    fids = [r['row_id'] for r in fil]
    assert set(fids) == set(tpl), '行集合与模板不一致——拒绝处理'
    assert len(fids) == len(set(fids)), 'row_id 重复——拒绝处理'
    tampered = [(r['row_id'], c) for r in fil for c in PREFILL if r[c] != tpl[r['row_id']][c]]
    assert not tampered, f'预填列被修改 {len(tampered)} 处(前5): {tampered[:5]}——拒绝处理'
    usable_before = sum(1 for r in gd if r['is_primary_disease'] == 'true'
                        and r['tier'] == 'core' and NUM.match(r['onset_typical_y']))

    # ---- G1 行级规范门 ----
    plans, incomplete, rejected = {}, [], []
    for r in fil:
        rid, g, depth = r['row_id'], r['hgnc_symbol'], r['check_depth']
        if batches and r['batch'] not in batches:
            rep.add(rid, r['batch'], 'skipped', ['--batch 过滤'])
            continue
        v = r['verdict'].strip()
        if not v:
            incomplete.append(rid)
            rep.add(rid, r['batch'], 'incomplete', ['verdict 空(该行核实未完成)'])
            continue
        reasons, norm, new = [], [], {}
        note_extra = []
        if v not in VOC['verdict']:
            reasons.append(f'verdict 非法: {v}')
        if not r['curator'].strip():
            reasons.append('curator 空')
        if not ISO.match(r['verification_date'].strip()):
            reasons.append(f"verification_date 非 ISO: '{r['verification_date']}'")
        if v == 'unverifiable':
            for c in ('source_type', 'source_primary_v', 'source_url_v', 'accessed_date_v', 'evidence_quote_v'):
                if r[c].strip():
                    reasons.append(f'unverifiable 但 {c} 非空')
            if not r['notes'].strip():
                reasons.append('unverifiable 无 notes(检索过程)')
        else:
            st = r['source_type'].strip()
            if st not in VOC['source_type']:
                reasons.append(f'source_type 非法/空: {st!r}')
            sp, tail = norm_source(r['source_primary_v'])
            if sp is None:
                reasons.append(f'source_primary_v 无法归一: {r["source_primary_v"][:50]!r}')
            else:
                new['source_primary'] = sp
                if tail:
                    note_extra.append(f'来源节段: {tail[:60]}')
                    if normalize:
                        norm.append(f'source_primary 尾注剥离: "{tail[:40]}"')
            if not r['source_url_v'].strip().startswith('http'):
                reasons.append(f'source_url_v 非 URL: {r["source_url_v"][:40]!r}')
            else:
                new['source_url'] = r['source_url_v'].strip()
            if ISO.match(r['accessed_date_v'].strip()):
                new['accessed_date'] = r['accessed_date_v'].strip()
            else:
                reasons.append(f'accessed_date_v 非 ISO: {r["accessed_date_v"]!r}')
            q = r['evidence_quote_v'].strip()
            if not q:
                reasons.append('evidence_quote_v 空')
            elif len(q) > 200:
                reasons.append(f'quote 超长({len(q)}>200字)')
            else:
                new['evidence_quote'] = q
            if depth == 'full':
                eg = r['evidence_grade_v'].strip()
                if eg not in VOC['evidence_grade']:
                    reasons.append(f'evidence_grade_v 非法/空: {eg!r}')
                else:
                    new['evidence'] = eg
        if (depth == 'light' or g in SPECIAL9) and v != 'unverifiable':
            if r['ataxia_core_feature_v'].strip() not in VOC['bool']:
                reasons.append('强制行 ataxia_core_feature_v 缺失/非法')

        # ---- onset 与词表列 ----
        typ_v = r['onset_typical_y_v'].strip()
        cur_num = bool(NUM.match(r['onset_typical_y_cur'].strip()))
        onset_to_na = cur_num and not typ_v and v == 'corrected'
        if v != 'unverifiable':
            if typ_v:
                if not NUM.match(typ_v):
                    reasons.append(f'onset_typical_y_v 非数值: {typ_v!r}')
                b, msg = norm_basis(r['onset_basis_v'], typ_v)
                if b is None:
                    reasons.append(msg)
                else:
                    new['onset_basis'] = b
                    if msg and normalize:
                        norm.append(msg)
                for vc, cc in (('onset_min_y_v', 'onset_min_y'), ('onset_max_y_v', 'onset_max_y'),
                               ('onset_typical_y_v', 'onset_typical_y')):
                    if r[vc].strip():
                        if NUM.match(r[vc].strip()):
                            new[cc] = r[vc].strip()
                        else:
                            reasons.append(f'{vc} 非数值: {r[vc][:20]!r}')
            for vc, cc in PAIRS:
                val = r[vc].strip()
                if not val:
                    continue
                val_cur = r[cc + '_cur'].strip()
                if cc == 'course':
                    w, msg = norm_course(val)
                    if w is None:
                        reasons.append(msg)
                        continue
                    val, val_cur = w, val_cur.replace(' ', '_').lower()
                elif cc == 'inheritance':
                    w, msg = norm_inh(val)
                    if w is None:
                        reasons.append(msg)
                        continue
                    if msg:
                        norm.append(msg)
                        note_extra.append(msg)
                    val = w
                elif cc in ('mutation_class', 'primary_pathology') and val != val_cur:
                    norm.append(f'{cc} 新值(词表外请人工抽查): "{val[:30]}"')
                if val != val_cur:
                    new[cc] = val
            acf_v = r['ataxia_core_feature_v'].strip()
            if acf_v in VOC['bool'] and acf_v != r['ataxia_core_feature_cur'].strip():
                new['ataxia_core_feature'] = acf_v
            mv = r['mech_class_v'].strip()
            if mv:
                if mv in VOC['mech_class']:
                    if mv != r['mech_class_cur'].strip():
                        new['mech_class'] = mv
                        norm.append(f'mech_class 本轮被改 {r["mech_class_cur"]}→{mv}(核实者声明文献矛盾)')
                else:
                    msg = f'mech_class_v 词表外值丢弃(本轮不裁 mech_class, 转 CL-014): "{mv}"'
                    norm.append(msg)
                    note_extra.append(msg)

        # ---- verdict 一致性(数值比较) ----
        onset_changed = False
        if typ_v and NUM.match(typ_v):
            onset_changed = (not cur_num) or not num_eq(typ_v, r['onset_typical_y_cur'])
        elif onset_to_na:
            onset_changed = True
        mm_changed = any(r[vc].strip() and not num_eq(r[vc].strip(), r[cc]) for vc, cc in
                         (('onset_min_y_v', 'onset_min_y_cur'), ('onset_max_y_v', 'onset_max_y_cur')))
        cat_changed = [cc for vc, cc in PAIRS if cc in new] + (['mech_class'] if 'mech_class' in new else [])
        # basis NA→有值 = 出处补全(同 source 回填), 不算值修正; 实值被改才算
        if 'onset_basis' in new and r['onset_basis_cur'].strip() not in ('', 'NA') \
                and new['onset_basis'] != r['onset_basis_cur'].strip():
            cat_changed.append('onset_basis')
        # source 回填(unverified→实源)本身构成 corrected 的合法差异; confirmed 行允许同时回填 source
        source_changed = 'source_primary' in new and new['source_primary'] != r['source_primary_cur'].strip()
        if v == 'confirmed' and (onset_changed or mm_changed or cat_changed):
            reasons.append(f'confirmed 但存在语义差异(onset={onset_changed or mm_changed}, {cat_changed})——应改标 corrected')
        if v == 'corrected' and not (onset_changed or mm_changed or cat_changed or source_changed):
            reasons.append('corrected 但无任何 *_v 与 cur 不同——应改标 confirmed')

        if reasons:
            rejected.append(rid)
            rep.add(rid, r['batch'], 'rejected', reasons, norm)
        else:
            plans[rid] = {'row': r, 'new': new, 'norm': norm, 'note_extra': note_extra,
                          'onset_to_na': onset_to_na, 'verdict': v}
            rep.add(rid, r['batch'], 'staged', [], norm)

    # ---- G2 onset 政策门 ----
    na_rows = sorted(rid for rid, p in plans.items() if p['onset_to_na'])
    core_na = [rid for rid in na_rows if tpl[rid]['tier_curated'] == 'core'
               and tpl[rid]['is_primary_disease'] == 'true']
    resc_core = [rid for rid, p in plans.items()
                 if p['row']['onset_typical_y_cur'].strip() == 'NA'
                 and NUM.match(p['new'].get('onset_typical_y', ''))
                 and tpl[rid]['tier_curated'] == 'core' and tpl[rid]['is_primary_disease'] == 'true']
    usable_after = usable_before - len(core_na) + len(resc_core)
    if na_rows and policy is None and not check_only:
        print(f'[G2] ABORT: {len(na_rows)} 行「草稿有数值→核实无出处」(其中 core 主病行 {len(core_na)} 个; '
              f'另救回 core 主病 NA 行 {len(resc_core)} 个)。s08 样本量: {usable_before} → {usable_after}。')
        print('       这是未决政策裁决(curation_log 未决事项), 必须显式指定其一:')
        print('       --onset-unverified-policy na            (置 NA, README 严格口径)')
        print('       --onset-unverified-policy keep-flagged  (保留草稿值+unverified 标记, 敏感性分析用)')
        print('       core 主病行清单:', ', '.join(core_na))
        rep.write(f'{root}/genes/verification/apply_report_{today}.tsv')
        sys.exit(1)
    if na_rows and policy == 'keep-flagged':
        for rid in na_rows:
            for k in ('onset_min_y', 'onset_max_y', 'onset_typical_y', 'onset_basis'):
                plans[rid]['new'].pop(k, None)
            plans[rid]['note_extra'].append(
                'onset_draft_unverified: 合格来源未证数值, 保留草稿值(敏感性分析需排除)')
            plans[rid]['onset_to_na'] = False

    # ---- G3 s00b 预演门(迭代到不动点) ----
    def simulate(blocked):
        hits = {}
        for sym, grows in rows_of_gene.items():
            acf_vals = []
            for r in grows:
                rid = f"{sym}::{r['disease_abbrev']}"
                val = r['ataxia_core_feature']
                if rid in plans and rid not in blocked and 'ataxia_core_feature' in plans[rid]['new']:
                    val = plans[rid]['new']['ataxia_core_feature']
                acf_vals.append(val)
            prim = [r for r in grows if r['is_primary_disease'] == 'true'][0]
            rule = peripheral_rule_of(sym, acf_vals, gt.get(sym, {}))
            if rule and prim['tier'] == 'core':
                hits[sym] = rule
        return hits

    blocked = set()
    hits = {}
    for _ in range(5):
        hits = simulate(blocked)
        new_block = {rid for sym in hits for rid, p in plans.items()
                     if tpl[rid]['hgnc_symbol'] == sym and 'ataxia_core_feature' in p['new']}
        if new_block <= blocked:
            break
        blocked |= new_block
    for rid in sorted(blocked):
        p = plans.pop(rid)
        sym = tpl[rid]['hgnc_symbol']
        rep.drop_staged(rid)
        rep.add(rid, p['row']['batch'], 'BLOCKED',
                [f'G3: acf 翻转将使策展 core 基因 {sym} 命中 peripheral 规则({hits.get(sym, "?")}) '
                 f'→ s00b 断言失败; 人工裁决: 推翻翻转, 或走 CL 条目降 {sym} 策展 tier'], p['norm'])

    # ---- 预测回填后的 effective tier 构成(信息项) ----
    comp = Counter()
    for sym, grows in rows_of_gene.items():
        acf_vals = []
        for r in grows:
            rid = f"{sym}::{r['disease_abbrev']}"
            val = r['ataxia_core_feature']
            if rid in plans and 'ataxia_core_feature' in plans[rid]['new']:
                val = plans[rid]['new']['ataxia_core_feature']
            acf_vals.append(val)
        prim = [r for r in grows if r['is_primary_disease'] == 'true'][0]
        rule = peripheral_rule_of(sym, acf_vals, gt.get(sym, {}))
        comp['peripheral' if rule else prim['tier']] += 1

    # ---- 回填 ----
    n_applied = Counter()
    if not check_only:
        for rid, p in sorted(plans.items()):
            r, new = p['row'], p['new']
            rep.drop_staged(rid)
            if p['verdict'] == 'unverifiable':
                rep.add(rid, r['batch'], 'unverifiable-untouched', ['保持 unverified(plan §1.2 处置)'], p['norm'])
                continue
            tr = gd_key[(tpl[rid]['hgnc_symbol'], tpl[rid]['disease_abbrev'])]
            applied = []
            for col in ('source_primary', 'source_url', 'accessed_date', 'evidence_quote', 'evidence',
                        'inheritance', 'course', 'mutation_class', 'primary_pathology', 'mech_class',
                        'ataxia_core_feature', 'onset_basis'):
                if col in new:
                    tr[col] = new[col]
                    applied.append(col)
            if p['onset_to_na'] and policy == 'na':
                for col in ('onset_min_y', 'onset_max_y', 'onset_typical_y', 'onset_class'):
                    tr[col] = 'NA'
                applied.append('onset→NA(G2 政策)')
            else:
                for col in ('onset_min_y', 'onset_max_y', 'onset_typical_y'):
                    if col in new:
                        tr[col] = new[col]
                        applied.append(col)
                if new.get('onset_typical_y') and NUM.match(new['onset_typical_y']):
                    oc = onset_class_of(float(new['onset_typical_y']), thr)
                    if oc != tr['onset_class']:
                        applied.append(f'onset_class {tr["onset_class"]}→{oc}')
                    tr['onset_class'] = oc
            tr['curator'] = r['curator'].strip()
            tr['curation_date'] = r['verification_date'].strip()
            tag = f'来源核实轮{r["verification_date"].strip()}: verdict={p["verdict"]}'
            if p['note_extra']:
                tag += '; ' + '; '.join(p['note_extra'])
            tr['note'] = (tr['note'] + ' | ' if tr['note'] and tr['note'] != 'NA' else '') + tag
            rep.add(rid, r['batch'], 'applied', [], p['norm'], applied)
            n_applied[p['verdict']] += 1
        if n_applied:
            bak = f'{gd_path}.bak_preapply_{today}'
            if not os.path.exists(bak):
                with open(bak, 'w', encoding='utf-8') as fb, open(gd_path, encoding='utf-8') as fo:
                    fb.write(fo.read())
            tmp = gd_path + '.tmp'
            with open(tmp, 'w', encoding='utf-8', newline='') as f:
                w = csv.DictWriter(f, fieldnames=list(gd[0].keys()), delimiter='\t', lineterminator='\n')
                w.writeheader()
                w.writerows(gd)
            os.replace(tmp, gd_path)

    # ---- 汇总 ----
    rpath = args[args.index('--report') + 1] if '--report' in args \
        else f'{root}/genes/verification/apply_report_{today}.tsv'
    rep.write(rpath)
    st = Counter(x['status'] for x in rep.rows)
    print(f'[apply] 模式: {"check-only(未写盘)" if check_only else "回填"}; 报告: {rpath}')
    print(f'[apply] 行处置: {dict(sorted(st.items()))}')
    if incomplete:
        print(f'[apply] 未完成行({len(incomplete)}): {", ".join(incomplete)}')
    if rejected:
        print(f'[apply] G1 拒绝行({len(rejected)}):')
        for rid in rejected:
            rs = [x for x in rep.rows if x['row_id'] == rid and x['status'] == 'rejected'][0]['reasons']
            print(f'    {rid}: {rs[:150]}')
    if na_rows or resc_core:
        print(f'[G2] onset 数值→NA 候选 {len(na_rows)} 行(core 主病 {len(core_na)}); 救回 core 主病 NA 行 {len(resc_core)}; '
              f's08 可用数 {usable_before}→{usable_after}; 政策: '
              f'{policy or "未指定(check-only 仅统计, 正式回填必须先裁决)"}')
    if blocked:
        print(f'[G3] BLOCKED {len(blocked)} 行(策展 core × acf 翻转 → s00b 断言): {", ".join(sorted(blocked))}')
    print(f'[G3] 回填后 effective tier 构成预测(以重跑 s00b 为准): {dict(sorted(comp.items()))}')
    print(f'[apply] s08 core 主病行 onset 可用数: 回填前 {usable_before}')
    if not check_only and n_applied:
        print(f'[apply] 已写盘 {dict(n_applied)}; 备份: {gd_path}.bak_preapply_{today}')
        print(f'[apply] 下一步: python -I scripts/s00b_derive_gene_table.py {root} '
              '(重裁 effective tier; --strict 下验收 #2 应转 PASS)')


if __name__ == '__main__':
    main()
