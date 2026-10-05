#!/usr/bin/env python3
"""verify_all.py — run verify_all.c on every sequence of data/, in parallel, and check the
first terms against the OEIS records (oeis.json).

    cc -O2 -o tools/verify_all tools/verify_all.c -lz
    python3 tools/verify_all.py [--jobs 4] [--only id,id,...] [--bin tools/verify_all]

Every row of every sequence is checked (see verify_all.c): the chunk format, the identity
a = k·L + d, k > d, k the least divisor of a − d above d, decomposition exactly when a > 2d,
and the catalogue's counts. Exit status 1 on any failure.
"""
import csv, json, pathlib, subprocess, sys, time, concurrent.futures as cf

ROOT = pathlib.Path(__file__).resolve().parent.parent
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
jobs, binp = int(arg('--jobs', 4)), arg('--bin', str(ROOT / 'tools' / 'verify_all'))
rows = list(csv.DictReader(open(ROOT / 'data' / 'catalog.csv', newline='')))
if '--only' in sys.argv:
    want = set(arg('--only', '').split(',')); rows = [r for r in rows if r['id'] in want or r['anumber'] in want]
oeis = json.load(open(ROOT / 'tools' / 'oeis.json'))
F = 'chunks chunk_rows terms decomposable level weight ties level_one forced amin amax kmax Lmax dmin dmax'.split()

def run(part):
    text = ''.join(f"{r['id']} " + ' '.join(r[f] for f in F) + '\n' for r in part)
    p = subprocess.run([binp], input=text, capture_output=True, text=True, cwd=ROOT)
    return p.stdout.splitlines()

t0 = time.time()
parts = [rows[i::jobs * 8] for i in range(jobs * 8)]
by = {r['id']: r for r in rows}
ok = fails = total = weights = 0
with cf.ThreadPoolExecutor(jobs) as ex:
    for out in ex.map(run, parts):
        for line in out:
            f = line.split(' ', 4)
            if f[1] != 'ok': fails += 1; print(line); continue
            r = by[f[0]]; total += int(f[2]); weights += int(f[3])
            o = oeis[r['anumber']]; want = o['terms'][int(r['n0']) - o['offset']:][:40]
            head = [int(x) for x in f[4].split(',')][:len(want)]
            if len(want) < 12 or head != want: fails += 1; print(f"{f[0]} FAIL first terms {head[:8]}… differ from the OEIS {want[:8]}…"); continue
            ok += 1
missing = len(rows) - ok - fails
if missing > 0: fails += missing; print(f'{missing} sequences gave no answer')
print(f"{ok} of {len(rows)} sequences verified: {total:,} rows, {weights:,} weights re-derived as the least divisor, "
      f"first terms equal to the OEIS; {fails} failures ({time.time() - t0:.0f} s)")
sys.exit(1 if fails else 0)
