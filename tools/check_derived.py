#!/usr/bin/env python3
"""check_derived.py — what the site shows of each sequence, checked against its data.

    python3 tools/check_derived.py [--shard i/n] [--jobs N] [--only id,id,...]

Each sequence of data/catalog.csv (or of the shard: every n-th from the i-th) is decoded once,
by compact.py's own decoder, and then:
  1  its CSV (the file the "Download CSV" buttons build) hashes to the SHA-256 published in
     data/sha256.txt;
  2  its gallery plate thumbs/<id>.webp, a lossless image, holds exactly the pixels thumbs.py
     draws from these rows (a few edge pixels may move with the platform's log; at most 0.05 %).
The pages, their OEIS links and the fingerprint on their buttons are checked by check_site.py.
Exit status 1 on any failure. Run in shards by .github/workflows/checks.yml.
"""
import csv, hashlib, math, pathlib, sys, time, concurrent.futures as cf
import numpy as np
from PIL import Image

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'tools'))
from compact import decode

arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d

# the plate, exactly as tools/thumbs.py draws it
S = 480; M = 0.03; W = (0x3b, 0x7f, 0xe6); LV = (0xe0, 0x62, 0x1a); A0 = 0.55
def plate(a, k, L):
    A = a.max(); m = k > 0
    X = np.log(k[m]) / math.log(A); Y = np.log(L[m]) / math.log(A); lev = k[m] > L[m]
    px = np.clip(((X + M) / (1 + 2 * M) * S).astype(int), 0, S - 2)
    py = np.clip(((1 - (Y + M) / (1 + 2 * M)) * S).astype(int), 0, S - 2)
    layers = []
    for sel in (~lev, lev):
        cnt = np.zeros((S, S), np.int32)
        for dx in (0, 1):
            for dy in (0, 1):
                np.add.at(cnt, (py[sel] + dy, px[sel] + dx), 1)
        layers.append(1 - (1 - A0) ** cnt)
    aw, al = layers
    alpha = al + aw * (1 - al)
    rgb = np.zeros((S, S, 3))
    for ch in range(3):
        rgb[..., ch] = np.where(alpha > 0, (LV[ch] * al + W[ch] * aw * (1 - al)) / np.maximum(alpha, 1e-9), 0)
    return np.dstack([rgb, alpha * 255]).round().astype(np.uint8)

def check(job):
    rec, sha = job
    i, A = rec['id'], rec['anumber']
    fails = []
    rows = []
    for c in range(int(rec['chunks'])):
        rows += decode((ROOT / 'data' / 'seq' / i / f'chunk-{c:03d}.bin.gz').read_bytes())
    n0 = int(rec['n0'])
    if len(rows) != int(rec['terms']): fails.append(f'{len(rows)} rows, catalogue says {rec["terms"]}')
    # 1 the CSV, as the buttons build it
    h = hashlib.sha256(b'n;a;weight;level;jump\n')
    h.update(''.join(f'{n0 + j};{a};{k if k else ""};{L if k else ""};{d}\n' for j, (a, d, k, L) in enumerate(rows)).encode())
    if h.hexdigest() != sha: fails.append(f'CSV SHA-256 {h.hexdigest()[:16]}… is not the published {str(sha)[:16]}…')
    # 2 the plate
    arr = np.array(rows, dtype=np.float64)
    want = plate(arr[:, 0], arr[:, 2], arr[:, 3])
    got = np.asarray(Image.open(ROOT / 'thumbs' / f'{i}.webp').convert('RGBA'))
    if got.shape != want.shape: fails.append(f'plate is {got.shape}, expected {want.shape}')
    else:
        diff = int(np.count_nonzero(np.any(got != want, axis=2)))
        if diff > S * S * 0.0005: fails.append(f'plate differs from the data in {diff} pixels')
    return i, A, len(rows), fails

if __name__ == '__main__':
    rows = list(csv.DictReader(open(ROOT / 'data' / 'catalog.csv', newline='')))
    si, sn = (int(x) for x in arg('--shard', '0/1').split('/'))
    rows = rows[si::sn]
    if '--only' in sys.argv:
        want = set(arg('--only', '').split(',')); rows = [r for r in rows if r['id'] in want or r['anumber'] in want]
    sha = {l.split()[1][len('decompwlj_'):-len('.csv')]: l.split()[0]
           for l in (ROOT / 'data' / 'sha256.txt').read_text().splitlines() if l.strip()}
    t0 = time.time(); bad = 0; terms = 0
    with cf.ProcessPoolExecutor(int(arg('--jobs', 4))) as ex:
        for i, A, n, fails in ex.map(check, [(r, sha.get(r['anumber'])) for r in rows], chunksize=4):
            terms += n
            for f in fails: bad += 1; print(f'  FAIL {A} ({i}): {f}')
    print(f'{len(rows)} sequences, {terms:,} rows: every CSV fingerprint and plate checked against the data; '
          f'{bad} failures ({time.time() - t0:.0f} s)')
    sys.exit(1 if bad else 0)
