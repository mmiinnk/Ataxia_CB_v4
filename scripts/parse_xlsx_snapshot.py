#!/usr/bin/env python3
"""Parse xlsx without openpyxl: zipfile + ElementTree. Read-only."""
import sys, zipfile, re
import xml.etree.ElementTree as ET

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
      'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}

def col_to_idx(ref):
    m = re.match(r'([A-Z]+)(\d+)', ref)
    col, row = m.group(1), int(m.group(2))
    ci = 0
    for ch in col:
        ci = ci * 26 + (ord(ch) - 64)
    return row, ci - 1

def main(path, out_dir):
    z = zipfile.ZipFile(path)
    # sheet names + order
    wb = ET.fromstring(z.read('xl/workbook.xml'))
    rels = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    rel_map = {rel.get('Id'): rel.get('Target') for rel in rels}
    sheets = []
    for sh in wb.find('m:sheets', NS):
        rid = sh.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
        target = rel_map[rid].lstrip('/')
        if not target.startswith('xl/'):
            target = 'xl/' + target
        if target not in z.namelist():
            target = rel_map[rid].lstrip('/')
        sheets.append((sh.get('name'), target))

    # shared strings (may not exist)
    shared = []
    if 'xl/sharedStrings.xml' in z.namelist():
        sst = ET.fromstring(z.read('xl/sharedStrings.xml'))
        for si in sst:
            shared.append(''.join(t.text or '' for t in si.iter('{%s}t' % NS['m'])))

    for name, target in sheets:
        print(f'===== SHEET: {name} ({target}) =====')
        root = ET.fromstring(z.read(target))
        rows = {}
        maxc = 0
        for c in root.iter('{%s}c' % NS['m']):
            ref = c.get('r')
            if not ref:
                continue
            row, ci = col_to_idx(ref)
            t = c.get('t')
            v = None
            if t == 'inlineStr':
                is_el = c.find('m:is', NS)
                if is_el is not None:
                    v = ''.join(tt.text or '' for tt in is_el.iter('{%s}t' % NS['m']))
            else:
                v_el = c.find('m:v', NS)
                if v_el is not None:
                    v = v_el.text
                    if t == 's' and v is not None:
                        v = shared[int(v)]
            if v is not None:
                rows.setdefault(row, {})[ci] = v
                maxc = max(maxc, ci)
        safe = re.sub(r'[^\w\-]', '_', name)
        out = f'{out_dir}/sheet__{safe}.tsv'
        with open(out, 'w', encoding='utf-8') as f:
            for r in sorted(rows):
                cells = rows[r]
                line = '\t'.join((cells.get(i, '') or '').replace('\n', ' ').replace('\t', ' ') for i in range(maxc + 1))
                f.write(f'{line}\n')
        print(f'  -> {len(rows)} rows, {maxc+1} cols, written to {out}')

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
