#!/usr/bin/env python3
"""fingerprint.py — the SHA-256 of every sequence's CSV, the site's integrity reference.

    python3 tools/fingerprint.py              write data/sha256.txt for every sequence
    python3 tools/fingerprint.py --check      recompute and compare with data/sha256.txt
    python3 tools/fingerprint.py --check --every 25     … for one sequence in 25

data/sha256.txt is in the format of sha256sum, one line per sequence:

    <sha256>  decompwlj_A000040.csv

where the hash is that of the exact file the site's "Download CSV" buttons produce:
"n;a;weight;level;jump\\n", then one line "n;a;k;L;d\\n" per term (weight and level empty
when the term does not decompose), decoded here by compact.py's own decoder. The page checks
its CSV against this hash before saving it, the viewer checks the numbers it shows, and a
downloaded file can be checked with:  sha256sum -c sha256.txt
"""
import csv, hashlib, pathlib, sys, concurrent.futures as cf

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from compact import decode

def csv_sha(rec):
    h = hashlib.sha256(b'n;a;weight;level;jump\n')
    n = int(rec['n0']); out = []
    for c in range(int(rec['chunks'])):
        rows = decode((ROOT / 'data' / 'seq' / rec['id'] / f'chunk-{c:03d}.bin.gz').read_bytes())
        for a, d, k, L in rows:
            out.append(f'{n};{a};{k if k else ""};{L if k else ""};{d}\n'); n += 1
        h.update(''.join(out).encode()); out = []
    assert n - int(rec['n0']) == int(rec['terms']), rec['id']
    return rec['anumber'], h.hexdigest()

if __name__ == '__main__':
    rows = list(csv.DictReader(open(ROOT / 'data' / 'catalog.csv', newline='')))
    path = ROOT / 'data' / 'sha256.txt'
    if '--every' in sys.argv: rows = rows[::int(sys.argv[sys.argv.index('--every') + 1])]
    with cf.ProcessPoolExecutor() as ex:
        got = dict(ex.map(csv_sha, rows, chunksize=8))
    if '--check' in sys.argv:
        ref = {l.split()[1][len('decompwlj_'):-len('.csv')]: l.split()[0] for l in path.read_text().splitlines() if l.strip()}
        bad = [a for a in got if ref.get(a) != got[a]]
        for a in bad: print(f'  FAIL {a}: data give {got[a][:16]}…, sha256.txt has {(ref.get(a) or "nothing")[:16]}…')
        print(f'{len(got) - len(bad)} of {len(got)} sequences match data/sha256.txt')
        sys.exit(1 if bad else 0)
    path.write_text(''.join(f'{got[a]}  decompwlj_{a}.csv\n' for a in sorted(got)))
    print(f'{len(got)} fingerprints written to {path.relative_to(ROOT)}')
