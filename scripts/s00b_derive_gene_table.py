#!/usr/bin/env python3
"""s00b_derive_gene_table.py — 由 ataxia_gene_disease.tsv 程序化派生 ataxia_genes.tsv。

方案 gene_curation_plan.md §5 派生规则:
  - 代表疾病行(is_primary_disease=true)提供基因级注释列;
  - onset_min_all_y/onset_max_all_y = 该基因全部 ataxia_core_feature=true 疾病行的最小/最大;
  - onset_span_flag = (max-min) > 20 (阈值经 CLI, config 建立前为引导期口径 CL-020);
  - is_polyQ_set = 任一疾病行 mutation_class == repeat_polyQ (contains 规则, 修复 sheet3 指出的
    CACNA1A 复合值与 plan §4.2 Fig4e 的矛盾, CL-019);
  - in_panel: 默认 true; ATXN8OS(plan §1.2 非编码RNA)与 MT-TK/MT-TI(实测四数据集 var 缺失)为 false (CL-009);
  - ensg/ensmusg/存在性旗标 join 自 genes/ensembl_mapping_measured.tsv (实测, CL-017)。
  - peripheral 派生层(CL-023 三层制, 2026-10-29): effective tier 按成文标准逐基因派生, 规则见
    peripheral_rule_of(); 策展层(疾病表)tier 不变; 验收 #1b 改为软窗口(默认 120-180, WARN 语义,
    plan §0.1 2026-10-29 修订; 上限为设计期估计, 非硬性阈值)。
同时对 sanity_markers.tsv 回填 ensg/ensmusg 列, 并跑方案 §0 验收自检(草稿模式=WARN, --strict=v1.0 门禁)。
用法: python -I scripts/s00b_derive_gene_table.py <project_root> [--onset-thresholds 1,13,30]
      [--span-threshold 20] [--union-window 120,180] [--strict]
"""
import csv
import os
import sys

IN_PANEL_FALSE = {
    'ATXN8OS': 'plan §1.2 非编码RNA (CL-009)',
    'MT-TK': '实测四数据集 var 均无 (CL-009)',
    'MT-TI': '实测四数据集 var 均无 (CL-009)',
}
FORBIDDEN_CORE = {'DNAJC5', 'HPCA', 'ZNF592', 'SLC1A6', 'PCP2', 'CALB1', 'GABRA6'}  # plan §7

GENE_COLS = ['hgnc_symbol', 'ensg_id', 'mgi_symbol', 'ensmusg_id', 'ortholog_type', 'tier',
             'peripheral_rule', 'in_panel',
             'mech_class', 'mech_class_alt', 'mutation_class', 'onset_typical_y', 'onset_basis',
             'onset_class', 'onset_min_all_y', 'onset_max_all_y', 'onset_span_flag', 'n_diseases',
             'disease_primary', 'inheritance', 'course', 'primary_pathology', 'ataxia_core_feature',
             'evidence', 'is_polyQ_set', 'clinical_group_draft', 'in_human_sepp', 'in_mouse_sepp',
             'in_aldinger', 'in_castonguay', 'source_primary', 'note']


def onset_class_of(y, thr):
    c, j, a = thr
    return 'congenital' if y < c else 'childhood' if y < j else 'juvenile' if y < a else 'adult'


def peripheral_rule_of(sym, gr, m):
    """CL-023 三层制: peripheral 派生规则(2026-10-29 用户裁决执行)。

    peripheral 是**派生层 effective tier**, 不照抄任何名单(xlsx sheet3 / 旧清单均不采用),
    由下列成文标准逐基因判定, 按顺序命中即记; 规则输入(ataxia_core_feature 注释或实测存在性)
    变化时重跑本脚本即自动重裁成员, 无需新的人工裁决。策展层(疾病表)tier 保持 core/extended。
      P1  ataxia 非核心表现: 无任何疾病行 ataxia_core_feature=true (ANY-row 语义, 与 CL-019 同构);
      P2a 主数据集不可测: Sepp human 与 mouse var 实测均无此特征;
      P2b 仅单一主数据集可测, 且两个外部验证集(Aldinger/Castonguay)均缺 → 无法跨物种/外部验证;
      P2c mtDNA 编码基因(MT- 前缀): 核文库中线粒体转录本测量效度存疑(ambient/胞质携带)。
          唯一非机械规则, 附行动项: s01 实测检出率/信号真实性, 若真实且结构化可升回 extended。
    """
    if not any(r['ataxia_core_feature'] == 'true' for r in gr):
        return 'P1'
    if m.get('in_human_sepp') == 'false' and m.get('in_mouse_sepp') == 'false':
        return 'P2a'
    if (m.get('in_human_sepp') == 'true') + (m.get('in_mouse_sepp') == 'true') <= 1 \
            and m.get('in_aldinger') == 'false' and m.get('in_castonguay') == 'false':
        return 'P2b'
    if sym.startswith('MT-'):
        return 'P2c'
    return ''


def main():
    root = sys.argv[1]
    strict = '--strict' in sys.argv
    thr = (1, 13, 30)
    span_thr = 20.0
    union_win = (120, 180)   # plan §0.1 软窗口(2026-10-29 修订, CL-002/CL-023); config 建立后以 config 为准
    if '--onset-thresholds' in sys.argv:
        thr = tuple(float(x) for x in sys.argv[sys.argv.index('--onset-thresholds') + 1].split(','))
    if '--span-threshold' in sys.argv:
        span_thr = float(sys.argv[sys.argv.index('--span-threshold') + 1])
    if '--union-window' in sys.argv:
        union_win = tuple(int(x) for x in sys.argv[sys.argv.index('--union-window') + 1].split(','))

    gd = os.path.join(root, 'genes/ataxia_gene_disease.tsv')
    mapf = os.path.join(root, 'genes/ensembl_mapping_measured.tsv')
    outf = os.path.join(root, 'genes/ataxia_genes.tsv')
    rows = list(csv.DictReader(open(gd, encoding='utf-8'), delimiter='\t'))
    mapping = {m['query_symbol']: m for m in csv.DictReader(open(mapf, encoding='utf-8'), delimiter='\t')}

    genes = {}
    for r in rows:
        genes.setdefault(r['hgnc_symbol'], []).append(r)

    out = []
    for sym in sorted(genes):
        gr = genes[sym]
        prim = [r for r in gr if r['is_primary_disease'] == 'true']
        assert len(prim) == 1, f'{sym}: is_primary_disease 行数 {len(prim)} != 1'
        p = prim[0]
        m = mapping.get(sym, {})
        # onset_class 数值行程序重算(§4.3), NA 行沿用标签
        ocls = p['onset_class']
        if p['onset_typical_y'] != 'NA':
            ocls = onset_class_of(float(p['onset_typical_y']), thr)
        # span: 跨 ataxia_core_feature=true 的疾病行
        vals = [float(r[k]) for r in gr if r['ataxia_core_feature'] == 'true'
                for k in ('onset_typical_y',) if r[k] != 'NA']
        mn = str(min(vals)) if vals else 'NA'
        mx = str(max(vals)) if vals else 'NA'
        span = 'true' if (vals and (max(vals) - min(vals)) > span_thr) else 'false'
        polyq = 'true' if any(r['mutation_class'] == 'repeat_polyQ' for r in gr) else 'false'
        in_panel = 'false' if sym in IN_PANEL_FALSE else 'true'
        # CL-023: peripheral 派生(effective tier); 命中 core 即注释矛盾, 直接报错而非静默降级
        prule = peripheral_rule_of(sym, gr, m)
        if prule and p['tier'] == 'core':
            raise AssertionError(f'{sym}: peripheral 规则 {prule} 命中 core 基因——'
                                 f'ataxia_core_feature/存在性旗标与 tier 矛盾, 请复核注释')
        eff_tier = 'peripheral' if prule else p['tier']
        note = p['note']
        if sym in IN_PANEL_FALSE:
            note = (note + ' | ' if note and note != 'NA' else '') + f'in_panel=false: {IN_PANEL_FALSE[sym]}'
        out.append(dict(
            hgnc_symbol=sym,
            ensg_id=m.get('ensg_id', 'NA'), mgi_symbol=m.get('mgi_symbol', 'NA'),
            ensmusg_id=m.get('ensmusg_id', 'NA'), ortholog_type=m.get('ortholog_type', 'NA'),
            tier=eff_tier, peripheral_rule=prule, in_panel=in_panel,
            mech_class=p['mech_class'], mech_class_alt=p['mech_class_alt'],
            mutation_class=p['mutation_class'],
            onset_typical_y=p['onset_typical_y'], onset_basis=p['onset_basis'], onset_class=ocls,
            onset_min_all_y=mn, onset_max_all_y=mx, onset_span_flag=span, n_diseases=str(len(gr)),
            disease_primary=p['disease_abbrev'], inheritance=p['inheritance'], course=p['course'],
            primary_pathology=p['primary_pathology'], ataxia_core_feature=p['ataxia_core_feature'],
            evidence=p['evidence'], is_polyQ_set=polyq,
            clinical_group_draft=p['clinical_group_draft'],
            in_human_sepp=m.get('in_human_sepp', 'NA'), in_mouse_sepp=m.get('in_mouse_sepp', 'NA'),
            in_aldinger=m.get('in_aldinger', 'NA'), in_castonguay=m.get('in_castonguay', 'NA'),
            source_primary=p['source_primary'], note=note))

    with open(outf, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=GENE_COLS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(out)
    print(f'[s00b] 写入 {outf} ({len(out)} 基因, {len(rows)} 疾病行)')

    # ---- sanity markers 回填映射列 ----
    mkf = os.path.join(root, 'genes/sanity_markers.tsv')
    mk = list(csv.DictReader(open(mkf, encoding='utf-8'), delimiter='\t'))
    for r in mk:
        m = mapping.get(r['gene_symbol'], {})
        r['ensg_id'] = m.get('ensg_id', 'NA')
        r['ensmusg_id'] = m.get('ensmusg_id', 'NA')
        r['mgi_symbol'] = m.get('mgi_symbol', 'NA')
    with open(mkf, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(mk[0].keys()), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(mk)
    print(f'[s00b] sanity_markers.tsv 回填映射列 ({len(mk)} 行)')

    # ---------------------------------------------------------------- 验收自检 ----
    results = []

    def check(name, ok, detail, gate=True):
        results.append((name, 'PASS' if ok else ('FAIL' if gate else 'WARN'), detail))

    core = [r for r in out if r['tier'] == 'core']
    core_panel = [r for r in core if r['in_panel'] == 'true']
    ext = [r for r in out if r['tier'] == 'extended']
    check('#1a core ≥ 100', len(core) >= 100, f'core={len(core)}')
    eff = [r for r in out if r['tier'] in ('core', 'extended')]
    lo, hi = union_win
    check(f'#1b effective core∪extended ∈ [{lo},{hi}] (软窗口, plan §0.1 2026-10-29 修订, CL-002/CL-023)',
          lo <= len(eff) <= hi,
          f'effective 总数={len(eff)}; peripheral={len(out) - len(eff)}(派生层, 不计数, CL-023)', gate=False)
    per_by_rule = {}
    for r in out:
        if r['tier'] == 'peripheral':
            per_by_rule.setdefault(r['peripheral_rule'], []).append(r['hgnc_symbol'])
    check('peripheral 派生层构成(信息项)', True,
          '; '.join(f'{k}({len(v)}): {sorted(v)}' for k, v in sorted(per_by_rule.items())) or '空', gate=False)
    from collections import Counter
    mech_c = Counter(r['mech_class'] for r in core)
    # #3 与 plan §0 验收 3 / §4.1 对齐(用户裁决 2026-11-24, CL-027):
    # plan 原文 = "每类在 core 中 ≥5, 否则该类在统计中并入 other(允许, 须留痕)" —— other 正是该规则的
    # 落脚点(残差兜底类), 故 ≥5 硬门槛只作用于**真实机制类**; 残差类 other/unknown 豁免, 改为信息项留痕,
    # 并写明下游用法(s07 富集排除残差类; s08 分层置换中残差类自成 stratum 且标 low n)。
    RESIDUAL_MECH = ('other', 'unknown')
    bad_mech = {k: v for k, v in mech_c.items() if v < 5 and k not in RESIDUAL_MECH}
    check('#3 mech_class core 每个真实机制类 ≥5(残差类 other/unknown 豁免; plan §0-3/§4.1, CL-027)',
          not bad_mech, f'{dict(mech_c)}; 不足的真实类={bad_mech or "无"}')
    resid = {k: sorted(r['hgnc_symbol'] for r in core if r['mech_class'] == k)
             for k in RESIDUAL_MECH if mech_c.get(k)}
    check('残差类成员(信息项; 须逐个有依据; s07 富集排除, s08 独立 stratum 标 low n)', True,
          '; '.join(f'{k}({len(v)}): {v}' for k, v in resid.items()) or '空', gate=False)
    oc = Counter(r['onset_class'] for r in core)
    ok_oc = all(oc.get(k, 0) >= (8 if k == 'congenital' else 10) for k in ('congenital', 'childhood', 'juvenile', 'adult'))
    check('#4 onset_class core 各类 ≥10(congenital ≥8)', ok_oc, f'{dict(oc)}')
    # #5 命中率: 人按 ensg (统计口径=core&in_panel), 鼠按 ensmusg
    def rate(rs, key):
        return sum(1 for r in rs if r[key] == 'true') / len(rs) if rs else float('nan')
    hs = rate(core_panel, 'in_human_sepp')
    ms = rate(core_panel, 'in_mouse_sepp')
    al = rate(core_panel, 'in_aldinger')
    ca = rate(core_panel, 'in_castonguay')
    check('#5a human 命中率 ≥0.95 (core&in_panel, ENSG 锚定实测)', hs >= 0.95, f'{hs:.3f}')
    check('#5b mouse 命中率 ≥0.90', ms >= 0.90, f'{ms:.3f}')
    check('#5c aldinger ≥0.90 / castonguay ≥0.90', al >= 0.90 and ca >= 0.90, f'aldinger={al:.3f} castonguay={ca:.3f}')
    # #7 §7 违规基因
    check('#7 §7 已知错误基因不入 core', not (FORBIDDEN_CORE & {r['hgnc_symbol'] for r in core}), '硬断言')
    # #2 来源(草稿期 WARN, strict 时 FAIL)
    unv = sum(1 for r in out if r['source_primary'].startswith('unverified'))
    check('#2 每基因合格 source_primary', unv == 0, f'{unv}/{len(out)} 行为 unverified (CL-015, 冻结前须逐条核实)',
          gate=strict)
    # #6 双标注(草稿期未做)
    check('#6 双标注一致率 ≥85%', False, '未执行 (CL-016, 冻结门禁)', gate=strict)
    # NA onset 统计
    na_core = [r['hgnc_symbol'] for r in core if r['onset_typical_y'] == 'NA']
    check('NA-onset core 数量(信息项)', True, f'{len(na_core)}: {na_core}', gate=False)
    polyq_core = [r['hgnc_symbol'] for r in core if r['is_polyQ_set'] == 'true']
    print(f'[s00b] is_polyQ_set core={len(polyq_core)}: {polyq_core} (+extended THAP11 → Fig4e 集含 CACNA1A, CL-019)')

    print(f'\n[s00b] 验收自检 (strict={strict}, 阈值 onset={thr} span>{span_thr}):')
    n_fail = 0
    for name, status, detail in results:
        print(f'  [{status}] {name}: {detail}')
        n_fail += status == 'FAIL'
    if n_fail:
        sys.exit(1)
    print('[s00b] 非 strict 门禁通过(FAIL 0)。仍为 WARN 的冻结门禁项 #2 来源核实 / #6 双标注按用户裁决'
          '顺延 v1.1; #3 残差类豁免见 CL-027 —— 处置总表 genes/curation_log.md §4')


if __name__ == '__main__':
    main()
