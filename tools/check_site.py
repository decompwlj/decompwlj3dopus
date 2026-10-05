#!/usr/bin/env python3
"""check_site.py — is the published site complete and in step with its catalogue?

    python3 tools/check_site.py                       the file checks below
    python3 tools/check_site.py --regen GEN id,id     also rebuild these sequences with the compiled
                                                      generator GEN and compare, byte for byte, with
                                                      the chunks in data/

Checks: every sequence of data/catalog.csv has its chunks, its gallery plate (thumbs/), its share
image (share/), its page (seq/<A-number>/) and its OEIS record (tools/oeis.json); nothing is left
over from a sequence that is gone; the counts written in index.html, the sitemap and the page list
agree with the catalogue.  The data itself is checked by audit.py.  Run by the GitHub workflow
.github/workflows/checks.yml on every pull request.
"""
import csv, json, pathlib, re, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
errors = []
def bad(msg):
    errors.append(msg)
    if len(errors) <= 40: print('  FAIL', msg)

rows = list(csv.DictReader(open(ROOT / 'data' / 'catalog.csv', newline='')))
N = len(rows)
ids = [r['id'] for r in rows]; ans = [r['anumber'] for r in rows]
for what, v in (('id', ids), ('A-number', ans)):
    dup = sorted({x for x in v if v.count(x) > 1}) if len(set(v)) != len(v) else []
    if dup: bad(f'duplicate {what}s in catalog.csv: {dup[:10]}')

oeis = json.load(open(ROOT / 'tools' / 'oeis.json'))
for r in rows:
    i, A = r['id'], r['anumber']
    d = ROOT / 'data' / 'seq' / i
    want = {f'chunk-{c:03d}.bin.gz' for c in range(int(r['chunks']))}
    have = {p.name for p in d.iterdir()} if d.is_dir() else set()
    if want - have: bad(f'{i}: missing {sorted(want - have)}')
    if have - want: bad(f'{i}: unexpected files {sorted(have - want)}')
    for p in (ROOT / 'thumbs' / f'{i}.webp', ROOT / 'share' / f'{i}.jpg'):
        if not p.is_file() or p.stat().st_size == 0: bad(f'{i}: missing {p.relative_to(ROOT)}')
    page = ROOT / 'seq' / A / 'index.html'
    if not page.is_file(): bad(f'{i}: missing seq/{A}/index.html')
    else:
        t = page.read_text()
        if f'share/{i}.jpg' not in t or f'thumbs/{i}.webp' not in t: bad(f'seq/{A}/: does not point to its images')
    if A not in oeis: bad(f'{i}: {A} missing from tools/oeis.json')

# leftovers of sequences that are gone
known, knownA = set(ids), set(ans)
for p in (ROOT / 'data' / 'seq').iterdir():
    if p.name not in known: bad(f'data/seq/{p.name}: not in the catalogue')
for sub, ext in (('thumbs', '.webp'), ('share', '.jpg')):
    for p in (ROOT / sub).iterdir():
        if p.suffix != ext or p.stem not in known: bad(f'{sub}/{p.name}: not in the catalogue')
for p in (ROOT / 'seq').iterdir():
    if p.is_dir() and p.name not in knownA: bad(f'seq/{p.name}/: not in the catalogue')

# counts written in the pages
html = (ROOT / 'index.html').read_text()
for pat in (r'<span id="homeTotal">(\d+)</span>', r'Search (\d+) sequences', r'for (\d+) integer sequences</title>'):
    m = re.search(pat, html)
    if not m: bad(f'index.html: {pat} not found')
    elif int(m.group(1)) != N: bad(f'index.html says {m.group(1)} sequences, the catalogue has {N}')
fams = {r['family'] for r in rows}
sm = (ROOT / 'sitemap.xml').read_text().count('<url>')
if sm != N + len(fams) + 3: bad(f'sitemap.xml has {sm} URLs, expected {N + len(fams) + 3}')
lst = (ROOT / 'seq' / 'index.html').read_text().count('<tr><td>')
if lst != N: bad(f'seq/index.html lists {lst} sequences, the catalogue has {N}')
# the published fingerprints: one per sequence, and on its page
shaf = ROOT / 'data' / 'sha256.txt'
if not shaf.is_file(): bad('missing data/sha256.txt (run tools/fingerprint.py)')
else:
    sha = {l.split()[1][10:-4]: l.split()[0] for l in shaf.read_text().splitlines() if l.strip()}
    for A in ans:
        if not re.fullmatch(r'[0-9a-f]{64}', sha.get(A, '')): bad(f'{A}: no SHA-256 in data/sha256.txt')
        elif f'data-sha256="{sha[A]}"' not in (ROOT / 'seq' / A / 'index.html').read_text(): bad(f'seq/{A}/: does not carry its SHA-256')
    for A in set(sha) - set(ans): bad(f'data/sha256.txt: {A} is not in the catalogue')
for p in ('learn/index.html', 'learn/learn.js', 'seq/csv.js', '404.html', 'robots.txt', 'og.png',
          'favicon.svg', 'favicon.ico', 'apple-touch-icon.png', 'vendor/three.module.min.js', 'vendor/addons/controls/OrbitControls.js', 'vendor/mp4-muxer.min.mjs'):
    if not (ROOT / p).is_file(): bad(f'missing {p}')

# rebuild some sequences with the C generator and compare byte for byte
if '--regen' in sys.argv:
    sys.path.insert(0, str(ROOT / 'tools'))
    from compact import encode
    gen, pick = sys.argv[sys.argv.index('--regen') + 1], sys.argv[sys.argv.index('--regen') + 2].split(',')
    by = {r['id']: r for r in rows}
    with tempfile.TemporaryDirectory() as tmp:
        for i in pick:
            subprocess.run([gen, tmp, i], check=True, stdout=subprocess.DEVNULL)
            for c in range(int(by[i]['chunks'])):
                lines = open(f'{tmp}/seq/{i}/chunk-{c:03d}.csv').read().split('\n')[1:]
                new = encode([tuple(map(int, l.split(','))) for l in lines if l])
                if new != (ROOT / 'data' / 'seq' / i / f'chunk-{c:03d}.bin.gz').read_bytes():
                    bad(f'{i}: chunk {c} differs from a fresh run of the generator')
            print(f'  regenerated {i}: identical' if not any(e.startswith(i + ':') for e in errors) else f'  regenerated {i}')

if errors:
    sys.exit(f'{len(errors)} problem(s)')
print(f'site complete: {N} sequences, {len(fams)} families, every page, plate, share image and OEIS record in place')
