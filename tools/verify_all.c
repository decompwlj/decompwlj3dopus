/* verify_all.c — check EVERY row of the published data, independently of the generator.
 *
 *     cc -O2 -o verify_all verify_all.c -lz
 *     python3 verify_all.py                      (runs it on all of data/, in parallel)
 *
 * Shares no code with decompwlj_gen.c: its own reader for the compact chunks (dwj1/dwj2 in
 * gzip), its own factorisation (trial division, Miller-Rabin with the 12 prime bases that are
 * exact below 3.3e24, Pollard-Brent), its own divisor walk. For every row of every sequence:
 *
 *   1  the chunk is well formed: tag, mode, row count, no trailing bytes, values < 2^53
 *   2  a(n) increases, chunk c starts where chunk c-1 ends, d(n) = a(n+1) - a(n)
 *   3  a term decomposes exactly when a > 2d; then k·L = a - d exactly and k > d
 *   4  k is the LEAST divisor of l = a - d exceeding d (all divisors of l are walked)
 *   5  the catalogue's counts and ranges are reproduced
 *
 * Input, one line per sequence on stdin:
 *   id chunks chunk_rows terms decomposable level weight ties level_one forced amin amax kmax Lmax dmin dmax
 * Output: one line per sequence: "<id> ok <rows> <weights re-derived> <first 40 terms>", or
 * "<id> FAIL <the first failure>"; exit status 1 if any failed.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>

typedef unsigned long long u64;
typedef unsigned __int128 u128;

static u64 mulm(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powm(u64 a, u64 e, u64 m) { u64 r = 1; a %= m; while (e) { if (e & 1) r = mulm(r, a, m); a = mulm(a, a, m); e >>= 1; } return r; }

static int prime(u64 n)
{
    static const u64 B[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < 2) return 0;
    for (int i = 0; i < 12; i++) if (n % B[i] == 0) return n == B[i];
    u64 d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        u64 x = powm(B[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int ok = 0;
        for (int r = 1; r < s && !ok; r++) { x = mulm(x, x, n); if (x == n - 1) ok = 1; }
        if (!ok) return 0;
    }
    return 1;
}

static u64 gcd(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

/* Pollard-Brent: a non-trivial factor of the odd composite n */
static u64 brent(u64 n)
{
    for (u64 c = 1;; c++) {
        u64 y = 2, x, g = 1, q = 1, ys = 2, m = 128, r = 1;
        do {
            x = y;
            for (u64 i = 0; i < r; i++) y = (mulm(y, y, n) + c) % n;
            u64 k = 0;
            do {
                ys = y;
                for (u64 i = 0; i < m && i < r - k; i++) { y = (mulm(y, y, n) + c) % n; q = mulm(q, x > y ? x - y : y - x, n); }
                g = gcd(q, n); k += m;
            } while (k < r && g == 1);
            r *= 2;
        } while (g == 1);
        if (g == n) { do { ys = (mulm(ys, ys, n) + c) % n; g = gcd(x > ys ? x - ys : ys - x, n); } while (g == 1); }
        if (g != n) return g;
    }
}

static u64 P[64]; static int E[64], NP;
static void addp(u64 p) { for (int i = 0; i < NP; i++) if (P[i] == p) { E[i]++; return; } P[NP] = p; E[NP++] = 1; }
static void split(u64 n)
{
    if (n == 1) return;
    if (prime(n)) { addp(n); return; }
    u64 f = brent(n); split(f); split(n / f);
}
static void factor(u64 n)
{
    NP = 0;
    static const int small[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
    for (int i = 0; i < 25; i++) while (n % small[i] == 0) { addp(small[i]); n /= small[i]; }
    split(n);
}

/* least divisor of l exceeding d, over every divisor of l */
static u64 best, lim;
static void walk(int i, u64 v)
{
    if (v >= best) return;                                     /* divisors only grow from here */
    if (i == NP) { if (v > lim) best = v; return; }
    u64 p = 1;
    for (int e = 0; e <= E[i]; e++) { walk(i + 1, v * p); if (e < E[i]) p *= P[i]; }   /* v·p divides l < 2^53 */
}
static u64 least_above(u64 l, u64 d) { factor(l); best = l + 1; lim = d; walk(0, 1); return best; }

/* ── reading a chunk ── */
static unsigned char *buf; static size_t cap;
static long readgz(const char *path)
{
    gzFile f = gzopen(path, "rb"); if (!f) return -1;
    size_t n = 0;
    for (;;) {
        if (n + 65536 > cap) { cap = cap ? cap * 2 : 1 << 22; buf = realloc(buf, cap); }
        int r = gzread(f, buf + n, 65536);
        if (r < 0) { gzclose(f); return -2; }
        if (r == 0) break;
        n += (size_t)r;
    }
    if (gzclose(f) != Z_OK) return -2;                         /* a CRC or length mismatch shows here */
    return (long)n;
}
static size_t pos, len;
static int bad;
static u64 rd(void)
{
    u64 v = 0; int sh = 0; unsigned char c;
    do { if (pos >= len || sh > 56) { bad = 1; return 0; } c = buf[pos++]; v |= (u64)(c & 127) << sh; sh += 7; } while (c & 128);
    return v;
}

#define FAIL(...) do { printf("%s FAIL ", id); printf(__VA_ARGS__); printf("\n"); return 1; } while (0)

static int verify(char *line)
{
    char id[64]; long chunks, crow, terms;
    u64 cat[12];
    if (sscanf(line, "%63s %ld %ld %ld %llu %llu %llu %llu %llu %llu %llu %llu %llu %llu %llu %llu", id, &chunks, &crow, &terms,
               &cat[0], &cat[1], &cat[2], &cat[3], &cat[4], &cat[5], &cat[6], &cat[7], &cat[8], &cat[9], &cat[10], &cat[11]) != 16) {
        printf("bad input line\n"); return 1;
    }
    const u64 LIMIT = 1ULL << 53;
    u64 dec = 0, lev = 0, wgt = 0, tie = 0, one = 0, forced = 0, amin = ~0ULL, amax = 0, kmax = 0, Lmax = 0, dmin = ~0ULL, dmax = 0;
    u64 expect_a = 0; int have = 0; long rows = 0;
    char head[1024]; int hl = 0; head[0] = 0;                 /* the first 40 terms, for the OEIS check */
    u64 *D = NULL; size_t dcap = 0;
    for (long c = 0; c < chunks; c++) {
        char path[512]; snprintf(path, sizeof path, "data/seq/%s/chunk-%03ld.bin.gz", id, c);
        long n = readgz(path);
        if (n < 0) FAIL("%s: cannot read (%s)", path, n == -1 ? "missing" : "gzip error or bad CRC");
        len = (size_t)n; pos = 0; bad = 0;
        if (len < 5 || memcmp(buf, "dwj", 3) || (buf[3] != '1' && buf[3] != '2')) FAIL("%s: not a dwj chunk", path);
        int m = buf[3] == '2' ? buf[4] : 0; pos = buf[3] == '2' ? 5 : 4;
        if (m > 2) FAIL("%s: unknown mode %d", path, m);
        u64 nr = rd(), a = rd();
        long want = crow < terms - c * crow ? crow : terms - c * crow;
        if (bad || (long)nr != want) FAIL("%s: %llu rows, expected %ld", path, nr, want);
        if (have && a != expect_a) FAIL("%s: starts at %llu, previous chunk ends at %llu", path, a, expect_a);
        if (dcap < nr) { dcap = nr; D = realloc(D, dcap * sizeof *D); }
        long long p = 0, pp = 0;
        for (u64 r = 0; r < nr; r++) {
            u64 v = rd();
            long long dz = v % 2 ? -(long long)((v + 1) / 2) : (long long)(v / 2);    /* zigzag */
            long long d = m == 0 ? (long long)v : m == 1 ? p + dz : 2 * p - pp + dz;
            if (bad || d <= 0 || (u64)d >= LIMIT) FAIL("%s: jump %lld out of range at row %llu", path, d, r);
            D[r] = (u64)d; pp = p; p = d;
        }
        for (u64 r = 0; r < nr; r++) {
            u64 s = rd(), d = D[r], k = 0, L = 0;
            if (bad) FAIL("%s: truncated", path);
            if (a >= LIMIT || a + d >= LIMIT) FAIL("%s: a = %llu passes 2^53", path, a);
            if (a > 2 * d) {
                if (!s) FAIL("%s: a = %llu, d = %llu should decompose", path, a, d);
                u64 l = a - d, h = s / 2;
                if (!h || l % h) FAIL("%s: a = %llu: stored factor %llu does not divide %llu", path, a, h, l);
                k = s % 2 ? l / h : h; L = s % 2 ? h : l / h;
                if ((u128)k * L != l || k <= d) FAIL("%s: a = %llu: k = %llu, L = %llu, d = %llu", path, a, k, L, d);
                if (s % 2 ? !(L < k) : !(k <= L)) FAIL("%s: a = %llu: the stored factor is not the smaller one", path, a);
                u64 kk = least_above(l, d);
                if (kk != k) FAIL("%s: a = %llu, d = %llu: k = %llu but the least divisor of %llu above d is %llu", path, a, d, k, l, kk);
                dec++; if (k > L) lev++; else wgt++;
                if (k == L) tie++;
                if (L == 1) one++;
                if ((u128)l <= (u128)d * d) {                  /* d² passes 2^64 when d > 2^32 */ if (!(k > L)) FAIL("%s: a = %llu: l <= d^2 but not level", path, a); forced++; }
                if (k > kmax) kmax = k;
                if (L > Lmax) Lmax = L;
            } else if (s) FAIL("%s: a = %llu, d = %llu should not decompose", path, a, d);
            if (rows < 40 && hl < 960) hl += snprintf(head + hl, sizeof head - hl, rows ? ",%llu" : "%llu", a);
            if (a < amin) amin = a;
            if (a > amax) amax = a;
            if (d < dmin) dmin = d;
            if (d > dmax) dmax = d;
            a += d; rows++;
        }
        if (pos != len) FAIL("%s: %zu trailing bytes", path, len - pos);
        expect_a = a; have = 1;
    }
    free(D);
    if (rows != terms) FAIL("%ld rows, catalogue says %ld", rows, terms);
    u64 got[12] = {dec, lev, wgt, tie, one, forced, amin, amax, kmax, Lmax, dmin, dmax};
    static const char *nm[12] = {"decomposable", "level", "weight", "ties", "level_one", "forced", "amin", "amax", "kmax", "Lmax", "dmin", "dmax"};
    for (int i = 0; i < 12; i++) if (got[i] != cat[i]) FAIL("catalogue %s = %llu, the rows give %llu", nm[i], cat[i], got[i]);
    printf("%s ok %ld %llu %s\n", id, rows, dec, head);
    return 0;
}

int main(void)
{
    char line[4096]; int fails = 0;
    while (fgets(line, sizeof line, stdin)) { fails += verify(line); fflush(stdout); }
    return fails ? 1 : 0;
}
