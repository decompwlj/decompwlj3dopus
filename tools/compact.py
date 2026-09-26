#!/usr/bin/env python3
"""compact.py — turn the generator's a,d,k,L chunks into the compact chunks the page loads.

    python3 compact.py <raw_dir> <site_data_dir>      e.g.  python3 compact.py raw ../data

Each compact chunk is '#a0=<first term>', 'd,k', then one 'd,k' row per term.  The page
rebuilds a as a running sum of d and L = (a - d)/k, and checks both.  Every chunk is
read back here and compared with its source row by row before the script finishes.
"""
import csv, pathlib, shutil, sys
src, dst = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
(dst / 'seq').mkdir(parents=True, exist_ok=True)
shutil.copy(src / 'catalog.csv', dst / 'catalog.csv')
total = 0
for rec in csv.DictReader(open(src / 'catalog.csv')):
    nxt = None
    for c in range(int(rec['chunks'])):
        rows = [tuple(map(int, l.split(','))) for l in open(src / 'seq' / rec['id'] / f'chunk-{c:03d}.csv').read().split('\n')[1:] if l]
        p = dst / 'seq' / rec['id'] / f'chunk-{c:03d}.csv'
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text('\n'.join([f'#a0={rows[0][0]}', 'd,k'] + [f'{d},{k}' for a, d, k, L in rows]) + '\n')
        total += p.stat().st_size
        lines = p.read_text().split('\n'); a = int(lines[0][4:])
        assert nxt is None or a == nxt, (rec['id'], c)
        for (a0, d0, k0, L0), line in zip(rows, lines[2:]):
            d, k = map(int, line.split(',')); L = 0 if k == 0 else (a - d) // k
            assert (a, d, k, L) == (a0, d0, k0, L0) and (k == 0 or (a - d) % k == 0), (rec['id'], a)
            a += d
        nxt = a
print(f'{total:,} bytes written; every row read back exactly')
