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
  8  the first rows (up to 40) equal the OEIS terms read from index n0 (checks n0 too)
  9  every value is below 2^53, so the browser's doubles hold it exactly

Usage:  python3 tools/audit.py [data_dir] [--sample N]
"""

import sys, os, csv, random, math, json

# The OEIS records (name, offset, first terms), from github.com/oeis/oeisdata by
# fetch_oeis.py.  The rows must equal terms[n0 - offset :], which checks n0 as well.
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'oeis.json')) as _f:
    OEIS = {k: (v['offset'], v['terms']) for k, v in json.load(_f).items()}

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
        want = data[int(rec['n0']) - off:][:40]
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
