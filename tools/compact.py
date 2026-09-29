#!/usr/bin/env python3
"""compact.py — turn the generator's a,d,k,L chunks into the compact chunks the page loads.

    python3 compact.py <raw_dir> <site_data_dir>      e.g.  python3 compact.py raw ../data

Each chunk is written as chunk-NNN.bin.gz: a small binary file, gzip-compressed (the page
inflates it; the files carry no timestamp, so a rebuild is byte-identical).

    'dwj2'                      4 bytes, the format tag
    m                           1 byte, how the jumps are coded (below)
    n, a0                       rows in the chunk, first term
    j[0] … j[n-1]               the jumps d, coded by m
    s[0] … s[n-1]               the smaller factor of a - d = k·L, and which one it is:
                                  0         the term does not decompose (k = 0)
                                  2k        k ≤ L  (weight class, ties included)
                                  2L + 1    L < k  (level class)

every number an unsigned LEB128 varint (7 bits a byte, low first). a is the running sum of d
from a0; the factor not stored is (a - d) divided by the one stored. Storing the smaller factor,
at most √(a - d), roughly halves the digits of the large weights and levels.

The jumps are coded three ways and the chunk keeps whichever compresses smallest:
    m = 0   d itself
    m = 1   d - d_prev                  (zigzag: v >= 0 -> 2v, v < 0 -> -2v - 1)
    m = 2   d - 2 d_prev + d_prevprev   (zigzag), with d_prev = d_prevprev = 0 before the first
Differences suit sequences whose jumps grow regularly: for a polynomial of degree 2 the second
differences are constant. A mode is used only if every coded value stays below 2^53, so the page
decodes it exactly. The older 'dwj1' chunks (m = 0, no mode byte) are still read.

Every chunk is read back here and compared with its source row by row before the script finishes.
"""
import csv, gzip, io, pathlib, shutil, sys

MAGIC1, MAGIC2 = b'dwj1', b'dwj2'
LIMIT = 1 << 53

def _uv(x, out):
    while x >= 0x80:
        out.append((x & 0x7f) | 0x80); x >>= 7
    out.append(x)

def _zz(v): return 2 * v if v >= 0 else -2 * v - 1
def _unzz(z): return z // 2 if z % 2 == 0 else -(z + 1) // 2

def _gzip(b):
    buf = io.BytesIO()
    with gzip.GzipFile(filename='', mode='wb', compresslevel=9, fileobj=buf, mtime=0) as g:
        g.write(bytes(b))
    return buf.getvalue()

def _jumps(ds, m):
    """the coded jumps for mode m, or None if a value would reach 2^53"""
    out, p, pp = [], 0, 0
    for d in ds:
        v = d if m == 0 else _zz(d - p) if m == 1 else _zz(d - 2 * p + pp)
        if v >= LIMIT: return None
        out.append(v); pp, p = p, d
    return out

def encode(rows):
    """rows: [(a, d, k, L), ...] consecutive terms -> gzip bytes (the smallest of the modes)."""
    ds = [d for a, d, k, L in rows]
    tail = bytearray()
    for a, d, k, L in rows:
        _uv(0 if k == 0 else (2 * k if k <= L else 2 * L + 1), tail)
    best = None
    for m in (0, 1, 2):
        js = _jumps(ds, m)
        if js is None: continue
        b = bytearray(MAGIC2); b.append(m)
        _uv(len(rows), b); _uv(rows[0][0], b)
        for v in js: _uv(v, b)
        g = _gzip(b + tail)
        if best is None or len(g) < len(best): best = g
    return best

def decode(data):
    """gzip bytes -> [(a, d, k, L), ...], checking every division."""
    b = gzip.decompress(data)
    if b[:4] == MAGIC2: m, pos = b[4], 5
    elif b[:4] == MAGIC1: m, pos = 0, 4
    else: raise ValueError('not a dwj1/dwj2 chunk')
    def rd():
        nonlocal pos
        v = sh = 0
        while True:
            c = b[pos]; pos += 1
            v |= (c & 0x7f) << sh; sh += 7
            if c < 0x80: return v
    n = rd(); a = rd()
    ds, p, pp = [], 0, 0
    for _ in range(n):
        v = rd()
        d = v if m == 0 else p + _unzz(v) if m == 1 else 2 * p - pp + _unzz(v)
        ds.append(d); pp, p = p, d
    rows = []
    for d in ds:
        s = rd()
        if s == 0: k = L = 0
        elif s % 2 == 0: k = s // 2; L, r = divmod(a - d, k); assert r == 0 and k <= L
        else: L = s // 2; k, r = divmod(a - d, L); assert r == 0 and L < k
        rows.append((a, d, k, L))
        a += d
    assert pos == len(b), 'trailing bytes'
    return rows

if __name__ == '__main__':
    src, dst = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    (dst / 'seq').mkdir(parents=True, exist_ok=True)
    shutil.copy(src / 'catalog.csv', dst / 'catalog.csv')
    total = 0
    for rec in csv.DictReader(open(src / 'catalog.csv')):
        nxt = None
        for c in range(int(rec['chunks'])):
            rows = [tuple(map(int, l.split(','))) for l in open(src / 'seq' / rec['id'] / f'chunk-{c:03d}.csv').read().split('\n')[1:] if l]
            p = dst / 'seq' / rec['id'] / f'chunk-{c:03d}.bin.gz'
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(encode(rows))
            total += p.stat().st_size
            back = decode(p.read_bytes())
            assert back == rows, (rec['id'], c)
            assert nxt is None or back[0][0] == nxt, (rec['id'], c)
            nxt = back[-1][0] + back[-1][1]
    print(f'{total:,} bytes written; every row read back exactly')
