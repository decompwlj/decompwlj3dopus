#!/usr/bin/env python3
"""thumbs.py — the gallery plates, one 480 px transparent PNG per sequence.

    python3 thumbs.py <raw_dir> <thumbs_dir>          e.g.  python3 thumbs.py raw ../thumbs

The weight-level plane on a square with equal aspect: X = log k / log A, Y = log L / log A,
A the largest term, so every sequence shares the diagonal k = L and the edge kL = A.
Points only; the page draws the frame.  Needs numpy and Pillow.
"""
import csv, math, pathlib, sys
import numpy as np
from PIL import Image
src, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
S = 480; M = 0.03; W = (0x3b, 0x7f, 0xe6); LV = (0xe0, 0x62, 0x1a); A0 = 0.55
for rec in csv.DictReader(open(src / 'catalog.csv')):
    rows = np.concatenate([np.loadtxt(src / 'seq' / rec['id'] / f'chunk-{c:03d}.csv', delimiter=',', skiprows=1, dtype=np.float64, ndmin=2)
                           for c in range(int(rec['chunks']))])
    a, d, k, L = rows.T
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
    Image.fromarray(np.dstack([rgb, alpha * 255]).round().astype(np.uint8), 'RGBA').save(out / f"{rec['id']}.png", optimize=True)
print('plates written to', out)
