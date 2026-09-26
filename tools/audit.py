#!/usr/bin/env python3
"""audit.py — independent verification of the generated CSVs.

Nothing here shares code with the C generator.  The weight is recomputed the
slow, obvious way — scan k = d+1, d+2, ... until one divides l — so the two
programmes can only agree by both being right.

Checks, per sequence:
  1  row count and chunk boundaries match the catalogue
  2  the terms are strictly increasing across the whole concatenation
  3  d(n) = a(n+1) - a(n) on every consecutive pair inside the file
  4  decomposability: k = 0 exactly when a <= 2d
  5  a = k*L + d on every decomposable row
  6  k really is the LEAST divisor of l exceeding d  (sampled, exhaustive)
  7  the catalogue's counts (decomposable / level / weight / ties / L=1,
     and the min-max columns) are reproduced from the rows
  8  the first 12 rows equal the OEIS DATA line read from index n0 (checks n0 too)
  9  every value is below 2^53, so the browser's doubles hold it exactly

Usage:  python3 tools/audit.py [data_dir] [--sample N]
"""

import sys, os, csv, random, math

# DATA lines as printed on oeis.org (fetched 25 Sep 2026), each with its OEIS offset.
# The rows must equal DATA[n0 - offset :], which checks n0 as well as the terms.
OEIS = {
 'A000027': (1, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]),
 'A000037': (1, [2,3,5,6,7,8,10,11,12,13,14,15,17,18,19]),
 'A000040': (1, [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47]),
 'A000201': (1, [1,3,4,6,8,9,11,12,14,16,17,19,21,22,24]),
 'A000217': (0, [0,1,3,6,10,15,21,28,36,45,55,66,78,91,105]),
 'A000290': (0, [0,1,4,9,16,25,36,49,64,81,100,121,144,169,196]),
 'A000292': (0, [0,1,4,10,20,35,56,84,120,165,220,286,364,455,560]),
 'A000326': (0, [0,1,5,12,22,35,51,70,92,117,145,176,210,247,287]),
 'A000330': (0, [0,1,5,14,30,55,91,140,204,285,385,506,650,819,1015]),
 'A000384': (0, [0,1,6,15,28,45,66,91,120,153,190,231,276,325,378]),
 'A000578': (0, [0,1,8,27,64,125,216,343,512,729,1000,1331,1728,2197,2744]),
 'A000959': (1, [1,3,7,9,13,15,21,25,31,33,37,43,49,51,63]),
 'A000960': (1, [1,3,7,13,19,27,39,49,63,79,91,109,133,147,181]),
 'A000961': (1, [1,2,3,4,5,7,8,9,11,13,16,17,19,23,25]),
 'A001248': (1, [4,9,25,49,121,169,289,361,529,841,961,1369,1681,1849,2209]),
 'A001358': (1, [4,6,9,10,14,15,21,22,25,26,33,34,35,38,39]),
 'A001359': (1, [3,5,11,17,29,41,59,71,101,107,137,149,179,191,197]),
 'A001481': (1, [0,1,2,4,5,8,9,10,13,16,17,18,20,25,26]),
 'A001651': (1, [1,2,4,5,7,8,10,11,13,14,16,17,19,20,22]),
 'A001694': (1, [1,4,8,9,16,25,27,32,36,49,64,72,81,100,108]),
 'A001855': (1, [0,1,3,5,8,11,14,17,21,25,29,33,37,41,45]),
 'A002088': (0, [0,1,2,4,6,10,12,18,22,28,32,42,46,58,64]),
 'A002113': (1, [0,1,2,3,4,5,6,7,8,9,11,22,33,44,55]),
 'A002145': (1, [3,7,11,19,23,31,43,47,59,67,71,79,83,103,107]),
 'A002378': (0, [0,2,6,12,20,30,42,56,72,90,110,132,156,182,210]),
 'A002620': (0, [0,0,1,2,4,6,9,12,16,20,25,30,36,42,49]),
 'A002808': (1, [4,6,8,9,10,12,14,15,16,18,20,21,22,24,25]),
 'A002858': (1, [1,2,3,4,6,8,11,13,16,18,26,28,36,38,47]),
 'A003052': (1, [1,3,5,7,9,20,31,42,53,64,75,86,97,108,110]),
 'A003714': (0, [0,1,2,4,5,8,9,10,16,17,18,20,21,32,33]),
 'A004202': (1, [2,5,6,10,11,12,17,18,19,20,26,27,28,29,30]),
 'A004207': (0, [1,1,2,4,8,16,23,28,38,49,62,70,77,91,101]),
 'A005101': (1, [12,18,20,24,30,36,40,42,48,54,56,60,66,70,72]),
 'A005117': (1, [1,2,3,5,6,7,10,11,13,14,15,17,19,21,22]),
 'A005153': (1, [1,2,4,6,8,12,16,18,20,24,28,30,32,36,40]),
 'A005214': (1, [1,3,4,6,9,10,15,16,21,25,28,36,45,49,55]),
 'A005349': (1, [1,2,3,4,5,6,7,8,9,10,12,18,20,21,24]),
 'A005408': (0, [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29]),
 'A005843': (0, [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28]),
 'A006218': (0, [0,1,3,5,8,10,14,16,20,23,27,29,35,37,41]),
 'A006446': (1, [1,2,3,4,6,8,9,12,15,16,20,24,25,30,35]),
 'A006512': (1, [5,7,13,19,31,43,61,73,103,109,139,151,181,193,199]),
 'A006881': (1, [6,10,14,15,21,22,26,33,34,35,38,39,46,51,55]),
 'A007088': (0, [0,1,10,11,100,101,110,111,1000,1001,1010,1011,1100,1101,1110]),
 'A007504': (0, [0,2,5,10,17,28,41,58,77,100,129,160,197,238,281]),
 'A008864': (1, [3,4,6,8,12,14,18,20,24,30,32,38,42,44,48]),
 'A024619': (1, [6,10,12,14,15,18,20,21,22,24,26,28,30,33,34]),
 'A024916': (1, [1,4,8,15,21,33,41,56,69,87,99,127,141,165,189]),
 'A026424': (1, [2,3,5,7,8,11,12,13,17,18,19,20,23,27,28]),
 'A052382': (1, [1,2,3,4,5,6,7,8,9,11,12,13,14,15,16]),
}

def least_divisor_above(l, d):
    """The definition, done the slow way: the least k > d dividing l."""
    k = d + 1
    while k * k <= l:
        if l % k == 0:
            return k
        k += 1
    # No divisor in (d, sqrt(l)], so every divisor above d exceeds sqrt(l) and its
    # cofactor e = l/k satisfies e < l/d, i.e. e <= (l-1)//d.  The least such k
    # comes from the LARGEST admissible e.
    e = min(int(math.isqrt(l)), (l - 1) // d)
    while e >= 1:
        if l % e == 0:
            return l // e
        e -= 1
    return l

def read_rows(base, rec):
    rows = []
    for c in range(int(rec['chunks'])):
        path = os.path.join(base, 'seq', rec['id'], f'chunk-{c:03d}.csv')
        with open(path) as f:
            head = f.readline().strip()
            assert head == 'a,d,k,L', f'{path}: header is {head!r}'
            n = 0
            for line in f:
                a, d, k, L = (int(x) for x in line.split(','))
                rows.append((a, d, k, L))
                n += 1
        expect = min(int(rec['chunk_rows']), int(rec['terms']) - c * int(rec['chunk_rows']))
        assert n == expect, f'{path}: {n} rows, expected {expect}'
    return rows

def main():
    base = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith('-') else 'data'
    sample = 400
    if '--sample' in sys.argv:
        sample = int(sys.argv[sys.argv.index('--sample') + 1])
    rng = random.Random(20260826)

    with open(os.path.join(base, 'catalog.csv'), newline='') as f:
        cat = list(csv.DictReader(f))

    total_rows = total_sampled = 0
    for rec in cat:
        rows = read_rows(base, rec)
        assert len(rows) == int(rec['terms']), f"{rec['id']}: row count"

        off, data = OEIS[rec['anumber']]
        want = data[int(rec['n0']) - off:][:12]
        head = [r[0] for r in rows[:len(want)]]
        assert len(want) >= 12 and head == want, f"{rec['id']}: rows {head} vs OEIS {want} (n0 {rec['n0']}, offset {off})"

        dec = lev = wgt = tie = one = forced = 0
        amin = kmax = Lmax = dmax = None
        amax = 0; dmin = None
        for i, (a, d, k, L) in enumerate(rows):
            assert a + d < 2 ** 53, f"{rec['id']}: {a} + {d} is not exact in a double"
            if i + 1 < len(rows):
                assert rows[i + 1][0] > a, f"{rec['id']}: not increasing at row {i}"
                assert rows[i + 1][0] - a == d, f"{rec['id']}: gap mismatch at row {i}"
            if a > 2 * d:
                assert k != 0, f"{rec['id']}: a={a} d={d} should decompose"
                assert k * L + d == a, f"{rec['id']}: identity fails at a={a}"
                assert (a - d) == k * L, f"{rec['id']}: l != k*L at a={a}"
                assert k > d, f"{rec['id']}: k <= d at a={a}"
                dec += 1
                if k > L: lev += 1
                else:     wgt += 1
                if k == L: tie += 1
                if L == 1: one += 1
                # forced: the window (d, sqrt l] is empty, so the term is level
                # whatever l factors into.  Exactly alpha = log d / log l >= 1/2.
                if a - d <= d * d:
                    assert k > L, f"{rec['id']}: l <= d^2 at a={a} but not level"
                    forced += 1
                kmax = k if kmax is None else max(kmax, k)
                Lmax = L if Lmax is None else max(Lmax, L)
            else:
                assert k == 0 and L == 0, f"{rec['id']}: a={a} d={d} should not decompose"
            amin = a if amin is None else min(amin, a)
            amax = max(amax, a)
            dmin = d if dmin is None else min(dmin, d)
            dmax = d if dmax is None else max(dmax, d)

        for name, got in (('decomposable', dec), ('level', lev), ('weight', wgt), ('forced', forced),
                          ('ties', tie), ('level_one', one), ('amin', amin), ('amax', amax),
                          ('kmax', kmax), ('Lmax', Lmax), ('dmin', dmin), ('dmax', dmax)):
            if name not in rec: continue        # older catalogues have no 'forced'
            assert int(rec[name]) == got, f"{rec['id']}: catalogue {name} {rec[name]} != {got}"

        # re-derive the weight from scratch on a random sample
        idx = [i for i in range(len(rows)) if rows[i][2] != 0]
        pick = rng.sample(idx, min(sample, len(idx)))
        for i in pick:
            a, d, k, L = rows[i]
            kk = least_divisor_above(a - d, d)
            assert kk == k, f"{rec['id']}: a={a} d={d} weight {k} != {kk}"
            assert (a - d) // kk == L, f"{rec['id']}: a={a} level mismatch"

        total_rows += len(rows)
        total_sampled += len(pick)
        print(f"  {rec['id']:<12} {len(rows):>7} rows  ok   "
              f"({dec} decomposable, {lev} level, {tie} ties, {len(pick)} re-derived)")

    print(f"\nall clear — {total_rows:,} rows checked, "
          f"{total_sampled:,} weights re-derived by exhaustive search")

if __name__ == '__main__':
    main()
