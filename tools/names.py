#!/usr/bin/env python3
"""names.py — give every sequence its OEIS name.

    python3 names.py <raw_dir>        e.g.  python3 names.py raw     (before compact.py)

Rewrites <raw_dir>/catalog.csv in place: 'name' becomes the OEIS name (%N line, without
its final period, from oeis.json), and the generator's own short name moves to a new
'alias' column when it differs, so the page still finds it by search. A 'keywords' column
keeps the OEIS keywords the gallery filters on: core and nice (%K line).
"""
import csv, json, pathlib, sys

HERE = pathlib.Path(__file__).resolve().parent
path = pathlib.Path(sys.argv[1]) / 'catalog.csv'
oeis = json.loads((HERE / 'oeis.json').read_text())
rows = list(csv.DictReader(open(path, newline='')))
fields = list(rows[0].keys())
if 'alias' not in fields: fields.insert(fields.index('name') + 1, 'alias')
if 'keywords' not in fields: fields.insert(fields.index('family') + 1, 'keywords')
KEEP = ('core', 'nice')
for r in rows:
    short = r.get('alias') or r['name']
    full = oeis[r['anumber']]['name']
    if full.endswith('.'): full = full[:-1]
    r['name'] = full
    r['alias'] = short if short.lower() != full.lower() else ''
    kw = oeis[r['anumber']].get('keywords', [])
    r['keywords'] = ' '.join(k for k in KEEP if k in kw)
with open(path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields, lineterminator='\n')
    w.writeheader(); w.writerows(rows)
print(f'{len(rows)} names from the OEIS')
