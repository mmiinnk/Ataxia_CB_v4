#!/usr/bin/env python3
"""extract_ensembl_mapping.py — 从四个 h5ad 的 var 实测提取 Ensembl ID 映射与跨数据集别名。

方法（全部实测, 满足硬约束10; 只读 var, 不触矩阵）:
  1. 读每个数据集 var/_index(Ensembl ID, spec §3.2: 四数据集 var.index 都是 Ensembl ID)与 var/feature_name(符号)。
  2. 人侧: query 符号先在 Sepp human feature_name 精确匹配; 未命中者用 Aldinger/Castonguay 补(记录来源数据集)。
  3. 别名发现: 以 ENSG 为锚, 若同一 ENSG 在另一数据集的 feature_name 与现行符号不同 → 别名(实测)。
  4. 鼠侧: 符号小写化匹配 mouse.h5ad; ENSMUSG 取 var/_index。BioMart ortholog_type 为待办(离线不可得),
     本表 ortholog_type=symbol_match(measured)/none_found(measured), 仅为占位(CL-017)。
输入: genes/source/gene_classification_table_v20261017.tsv + genes/sanity_markers.tsv 的符号并集
输出: genes/ensembl_mapping_measured.tsv, genes/gene_aliases.tsv
用法: python -I scripts/extract_ensembl_mapping.py <project_root>
"""
import csv
import os
import sys

import h5py

DATA = {  # key → (path, species)
    'human_sepp': ('/mnt/data/xyy/Ataxia_CB/human.h5ad', 'human'),
    'mouse_sepp': ('/mnt/data/xyy/Ataxia_CB/mouse.h5ad', 'mouse'),
    'aldinger': ('/mnt/data/xyy/Ataxia_CB/aldinger2021_fetal_snRNAseq.h5ad', 'human'),
    'castonguay': ('/mnt/data/xyy/Ataxia_CB/castonguay2026_adult_eQTL.h5ad', 'human'),
}
TODAY = '2026-10-28'


def _decode(arr):
    return [x.decode('utf-8') if isinstance(x, bytes) else x for x in arr]


def read_var(path):
    """Return (ids, symbols): parallel lists from var/_index and var/feature_name."""
    with h5py.File(path, 'r') as f:
        v = f['var']
        idx_key = v.attrs.get('_index', '_index')
        if isinstance(idx_key, bytes):
            idx_key = idx_key.decode()
        ids = _decode(v[idx_key][:])
        d = v['feature_name']
        if isinstance(d, h5py.Group):          # categorical
            cats = _decode(d['categories'][:])
            syms = [cats[i] for i in d['codes'][:]]
        else:
            syms = _decode(d[:])
    return ids, syms


def main():
    root = sys.argv[1]
    genes_dir = os.path.join(root, 'genes')

    # ---- query 符号集合: 疾病表快照 + sanity markers ----
    snap = os.path.join(genes_dir, 'source/gene_classification_table_v20261017.tsv')
    queries = [r['基因符号'].strip() for r in csv.DictReader(open(snap, encoding='utf-8'), delimiter='\t')]
    markers = os.path.join(genes_dir, 'sanity_markers.tsv')
    if os.path.exists(markers):
        queries += [r['gene_symbol'].strip() for r in csv.DictReader(open(markers, encoding='utf-8'), delimiter='\t')]
    queries = sorted(set(queries))

    # ---- 读四个数据集 var ----
    vars_ = {}
    for key, (path, sp) in DATA.items():
        ids, syms = read_var(path)
        assert all(str(i).startswith('ENSG') or str(i).startswith('ENSMUSG') for i in ids[:50]), \
            f'{key}: var index 不是 Ensembl ID'
        by_sym, by_id = {}, {}
        for i, s in zip(ids, syms):
            by_sym.setdefault(s if sp == 'human' else s.lower(), i)
            by_id.setdefault(i, s)
        vars_[key] = dict(ids=set(ids), by_sym=by_sym, by_id=by_id, species=sp)
        print(f'[map] {key}: {len(ids)} vars')

    alias_rows = []   # (alias, hgnc_symbol, ensg, mgi, ensmusg, source, note)
    out_rows = []
    for q in queries:
        row = dict(query_symbol=q, ensg_id='NA', ensg_source='NA', mgi_symbol='NA', ensmusg_id='NA',
                   ortholog_type='NA', symbol_human_sepp='NA', symbol_aldinger='NA',
                   symbol_castonguay='NA', symbol_mouse_sepp='NA', in_human_sepp='false',
                   in_mouse_sepp='false', in_aldinger='false', in_castonguay='false',
                   match_method='NA', note='')
        # 人侧 ENSG: Sepp → Aldinger → Castonguay 依次尝试(符号精确匹配)
        ensg, src = None, None
        for key in ('human_sepp', 'aldinger', 'castonguay'):
            v = vars_[key]
            if q in v['by_sym']:
                ensg, src = v['by_sym'][q], key
                break
        if ensg:
            row['ensg_id'], row['ensg_source'] = ensg, src
            row['match_method'] = f'symbol_exact({src})'
            # ENSG 锚定: 各人类数据集里这个 ENSG 用的符号
            for key, col in (('human_sepp', 'symbol_human_sepp'), ('aldinger', 'symbol_aldinger'),
                             ('castonguay', 'symbol_castonguay')):
                v = vars_[key]
                if ensg in v['by_id']:
                    s = v['by_id'][ensg]
                    row[col] = s
                    if s != q:
                        alias_rows.append((s, q, ensg, 'NA', 'NA',
                                           f'measured:{key} var ENSG-anchored {TODAY}',
                                           '同一 ENSG 在该数据集使用不同符号'))
                row['in_' + ('human_sepp' if key == 'human_sepp' else key)] = \
                    'true' if ensg in v['by_id'] else 'false'
        else:
            row['match_method'] = 'not_found_human(all datasets)'
        # 鼠侧: 符号小写匹配(占位 ortholog_type, 待 BioMart, CL-017)
        vm = vars_['mouse_sepp']
        if q.lower() in vm['by_sym']:
            row['ensmusg_id'] = vm['by_sym'][q.lower()]
            row['mgi_symbol'] = vm['by_id'][row['ensmusg_id']]
            row['symbol_mouse_sepp'] = row['mgi_symbol']
            row['in_mouse_sepp'] = 'true'
            row['ortholog_type'] = 'symbol_match(measured)'
        else:
            row['ortholog_type'] = 'none_found(measured)' if ensg else 'NA'
        out_rows.append(row)

    # 鼠侧别名特例(实测探针 2026-10-28): MRE11 → Mre11a
    if 'mre11a' in vars_['mouse_sepp']['by_sym']:
        ensmusg = vars_['mouse_sepp']['by_sym']['mre11a']
        mgi = vars_['mouse_sepp']['by_id'][ensmusg]
        alias_rows.append(('Mre11a', 'MRE11', 'NA', mgi, ensmusg,
                           f'measured:mouse_sepp var symbol probe {TODAY}',
                           '人 MRE11 的鼠同源在 Sepp mouse var 中以旧名 Mre11a 存在; 修正 in_mouse_sepp'))
        for r in out_rows:
            if r['query_symbol'] == 'MRE11':
                r['ensmusg_id'], r['mgi_symbol'], r['symbol_mouse_sepp'] = ensmusg, mgi, mgi
                r['in_mouse_sepp'] = 'true'
                r['ortholog_type'] = 'symbol_match(measured,alias Mre11a)'
                r['note'] = '鼠侧以 Mre11a 命中(实测)'

    # 手工种子别名(xlsx 备注/plan §6.5, 未经数据集证实者标 pending)
    manual = [('C10orf2', 'TWNK', 'xlsx 备注"旧名C10orf2"'), ('PEO1', 'TWNK', '常用旧名(unverified)'),
              ('ADCK3', 'COQ8A', 'xlsx 备注"旧名ADCK3"'), ('ADCK4', 'COQ8B', 'xlsx 备注"旧名ADCK4"'),
              ('C5orf42', 'CPLANE1', 'xlsx 备注"旧名C5orf42"'), ('KIAA0226', 'RUBCN', 'xlsx 备注"旧名KIAA0226"; plan §6.5 警示 Castonguay 中 KIAA0226→RUBCNL 的口径差异, s00 时按 ENSG 核对')]
    measured_aliases = {(a, g) for a, g, *_ in alias_rows}
    for a, g, why in manual:
        if (a, g) not in measured_aliases:
            alias_rows.append((a, g, 'NA', 'NA', 'NA', f'documented:{why} (pending dataset confirmation)', ''))

    # ---- 写映射表 ----
    mapf = os.path.join(genes_dir, 'ensembl_mapping_measured.tsv')
    cols = ['query_symbol', 'ensg_id', 'ensg_source', 'mgi_symbol', 'ensmusg_id', 'ortholog_type',
            'symbol_human_sepp', 'symbol_aldinger', 'symbol_castonguay', 'symbol_mouse_sepp',
            'in_human_sepp', 'in_mouse_sepp', 'in_aldinger', 'in_castonguay', 'match_method', 'note']
    with open(mapf, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(out_rows)

    # ---- 写别名表 ----
    alf = os.path.join(genes_dir, 'gene_aliases.tsv')
    with open(alf, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['alias', 'hgnc_symbol', 'ensg_id', 'mgi_symbol', 'ensmusg_id', 'source', 'note'])
        seen = set()
        for row_ in alias_rows:
            if (row_[0], row_[1]) in seen:
                continue
            seen.add((row_[0], row_[1]))
            w.writerow(row_)

    # ---- 摘要 ----
    n_ensg = sum(1 for r in out_rows if r['ensg_id'] != 'NA')
    n_ensmusg = sum(1 for r in out_rows if r['ensmusg_id'] != 'NA')
    miss_human = [r['query_symbol'] for r in out_rows if r['ensg_id'] == 'NA']
    print(f'[map] 写入 {mapf} ({len(out_rows)} 行) 与 {alf} ({len(alias_rows)} 别名)')
    print(f'[map] ENSG 命中: {n_ensg}/{len(out_rows)}; ENSMUSG 命中: {n_ensmusg}/{len(out_rows)}')
    print(f'[map] 人侧完全未命中: {miss_human}')
    print(f'[map] 实测发现的跨数据集别名: {sorted(measured_aliases)}')


if __name__ == '__main__':
    main()
