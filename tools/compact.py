#!/usr/bin/env python3
"""compact.py — turn the generator's a,d,k,L chunks into the compact chunks the page loads.

    python3 compact.py <raw_dir> <site_data_dir>      e.g.  python3 compact.py raw ../data

Each chunk is written as chunk-NNN.bin.gz: a small binary file, gzip-compressed (the page
inflates it; the files carry no timestamp, so a rebuild is byte-identical).

    'dwj1'                      4 bytes, the format tag
    n, a0                       rows in the chunk, first term
    d[0] … d[n-1]               the jumps
    s[0] … s[n-1]               the smaller factor of a - d = k·L, and which one it is:
                                  0         the term does not decompose (k = 0)
                                  2k        k ≤ L  (weight class, ties included)
                                  2L + 1    L < k  (level class)

every number an unsigned LEB128 varint (7 bits a byte, low first). a is the running sum of d
from a0; the factor not stored is (a - d) divided by the one stored. Storing the smaller factor,
at most √(a - d), roughly halves the digits of the large weights and levels.

Every chunk is read back here and compared with its source row by row before the script finishes.
"""
import csv, gzip, io, pathlib, shutil, sys

MAGIC = b'dwj1'

def _uv(x, out):
    while x >= 0x80:
        out.append((x & 0x7f) | 0x80); x >>= 7
    out.append(x)

def encode(rows):
    """rows: [(a, d, k, L), ...] consecutive terms -> gzip bytes."""
    b = bytearray(MAGIC)
    _uv(len(rows), b); _uv(rows[0][0], b)
    for a, d, k, L in rows: _uv(d, b)
    for a, d, k, L in rows:
        _uv(0 if k == 0 else (2 * k if k <= L else 2 * L + 1), b)
    buf = io.BytesIO()
    with gzip.GzipFile(filename='', mode='wb', compresslevel=9, fileobj=buf, mtime=0) as g:
        g.write(bytes(b))
    return buf.getvalue()

def decode(data):
    """gzip bytes -> [(a, d, k, L), ...], checking every division."""
    b = gzip.decompress(data)
    assert b[:4] == MAGIC, 'not a dwj1 chunk'
    pos = 4
    def rd():
        nonlocal pos
        v = sh = 0
        while True:
            c = b[pos]; pos += 1
            v |= (c & 0x7f) << sh; sh += 7
            if c < 0x80: return v
    n = rd(); a = rd()
    ds = [rd() for _ in range(n)]
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
