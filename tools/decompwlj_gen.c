/* decompwlj_gen.c — data generator for the decompwlj 3-D atlas.
 *
 *   a(n) = k(n) * L(n) + d(n)          weight x level + jump
 *
 *   d(n) = a(n+1) - a(n)
 *   l(n) = a(n) - d(n)   if a(n) - d(n) > d(n),  else 0      [ i.e. a(n) > 2 d(n) ]
 *   k(n) = least divisor of l(n) strictly exceeding d(n)      [ 0 if l(n) = 0 ]
 *   L(n) = l(n) / k(n)                                        [ 0 if k(n) = 0 ]
 *
 * Reference kernel (PARI/GP, project canon, tools/decompwlj.gp):
 *   decomp(a,b)={my(d=b-a,l); if(a<=2*d,return([0,0,d])); l=a-d;
 *                fordiv(l,k, if(k>d, return([k,l/k,d])))}
 *
 * This C engine is the factor-based variant (decomp_fact): l is factored, the
 * divisor lattice is walked, and the least divisor exceeding d is kept — O(tau(l)).
 *
 * Output, per sequence <id>:
 *   data/seq/<id>/chunk-000.csv ...   header "a,d,k,L", CHUNK_ROWS rows each
 *   data/catalog.csv                  one row per sequence
 *
 * IMPORTANT — every term is written, decomposable or not.  A non-decomposable
 * term is written with k = L = 0.  Row i of the concatenated chunks is therefore
 * always sequence index n = n0 + i.  (See finding 6 of the plates record: the
 * census engines that *skip* non-decomposable terms silently shift every index.)
 *
 * History
 *   26 Aug 2026  ten sequences (naturals ... digitadd), 2e5 terms
 *   24 Sep 2026  twenty (odd ... wythoff added)
 *   25 Sep 2026  fifty: the whole 50-sequence census, every sequence at 1e5 terms
 *                (A007088 at 65,535: a(65536) = 10^16 passes 2^53).  Row 0 is
 *                a(OEIS offset), except A002620 and A004207, which repeat their
 *                first value and so start one index later; n0 is always the
 *                true OEIS index of row 0, checked against the DATA lines.
 *   26 Sep 2026  one hundred: fifty more (odious ... binpal), for the GitHub Pages
 *                site.  The 24 of them also in wlj-atlas agree with it row for row.
 *   26 Sep 2026  two hundred: a hundred more (nnp3 ... n4p3); the 12 of them in
 *                wlj-atlas agree with it row for row.
 *   26 Sep 2026  four hundred: two hundred more (gon9 ... r5_024), every one
 *                checked against its OEIS terms and offset (tools/oeis.json).
 *   26 Sep 2026  eight hundred: four hundred more, found by family in the OEIS
 *                data (prime classes, linear and quadratic prime conditions,
 *                binary forms, residue classes, polynomials, Beatty sequences)
 *                and checked the same way.
 *   26 Sep 2026  one thousand: two hundred more (primes k^2 + c, n^2 + n + c,
 *                k ones in binary, divisor-function pairs, sums of consecutive
 *                primes, power sieves, Beatty constants, and from the same pools).
 *   26 Sep 2026  twelve hundred: one hundred base-dependent sequences (OEIS keyword
 *                "base") and one hundred classic ones (keywords "core", "nice").
 *   26 Sep 2026  sixteen hundred: four hundred more, mostly polynomial prime
 *                conditions read from the OEIS names (n with f(n) prime, primes of
 *                the form f(k), primes p with f(p) prime) and semiprime and
 *                squarefree variants, plus more from the verified pools.
 *
 * Build:  cc -O2 -o decompwlj_gen decompwlj_gen.c -lm
 * Run:    ./decompwlj_gen <outdir>            all sequences
 *         ./decompwlj_gen <outdir> primes     one sequence
 *         ./decompwlj_gen <outdir> primes 50000
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint8_t  u8;
typedef unsigned __int128 u128;

#define CHUNK_ROWS 50000

/* ------------------------------------------------------------------ */
/* smallest-prime-factor sieve + small prime list                      */
/* ------------------------------------------------------------------ */

#define SPF_MAX (1u << 23)          /* 8 388 608 */
static u32 *spf = NULL;
static u32 *smallp = NULL;          /* primes below 65536 */
static int  nsmallp = 0;

static void sieve_init(void)
{
    spf = calloc(SPF_MAX + 1, sizeof(u32));
    if (!spf) { fprintf(stderr, "spf alloc failed\n"); exit(1); }
    for (u32 i = 2; i <= SPF_MAX; i++) {
        if (!spf[i]) for (u64 j = i; j <= SPF_MAX; j += i) if (!spf[j]) spf[j] = i;
    }
    smallp = malloc(7000 * sizeof(u32));
    for (u32 i = 2; i < 65536; i++) if (spf[i] == i) smallp[nsmallp++] = i;
}

/* ------------------------------------------------------------------ */
/* Miller-Rabin (deterministic for 64-bit) and Brent's rho             */
/* ------------------------------------------------------------------ */

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }

static u64 powmod(u64 a, u64 e, u64 m)
{
    u64 r = 1; a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}

static const u64 mr_bases[12] = {2,3,5,7,11,13,17,19,23,29,31,37};

static int is_prime(u64 n)
{
    if (n < 2) return 0;
    for (int i = 0; i < 12; i++) { u64 p = mr_bases[i]; if (n % p == 0) return n == p; }
    u64 d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (int bi = 0; bi < 12; bi++) {
        u64 a = mr_bases[bi];
        u64 x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int ok = 0;
        for (int i = 1; i < s; i++) { x = mulmod(x, x, n); if (x == n - 1) { ok = 1; break; } }
        if (!ok) return 0;
    }
    return 1;
}

static u64 gcd_u64(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

static u64 rho(u64 n)
{
    if (!(n & 1)) return 2;
    for (u64 c = 1; ; c++) {
        u64 x = 2, y = 2, d = 1;
        do {
            x = (mulmod(x, x, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            d = gcd_u64(x > y ? x - y : y - x, n);
        } while (d == 1);
        if (d != n) return d;
    }
}

/* ------------------------------------------------------------------ */
/* factorisation                                                       */
/* ------------------------------------------------------------------ */

static u64 g_pf[64]; static int g_pe[64]; static int g_npf;

static void push_factor(u64 p)
{
    for (int i = 0; i < g_npf; i++) if (g_pf[i] == p) { g_pe[i]++; return; }
    g_pf[g_npf] = p; g_pe[g_npf] = 1; g_npf++;
}

static void factor_rec(u64 n)
{
    if (n == 1) return;
    if (n <= SPF_MAX) { while (n > 1) { u64 p = spf[n]; push_factor(p); n /= p; } return; }
    if (is_prime(n)) { push_factor(n); return; }
    u64 d = rho(n);
    factor_rec(d);
    factor_rec(n / d);
}

static void factor(u64 n)
{
    g_npf = 0;
    if (n <= SPF_MAX) { while (n > 1) { u64 p = spf[n]; push_factor(p); n /= p; } return; }
    for (int i = 0; i < nsmallp; i++) {
        u32 p = smallp[i];
        if ((u64)p * p > n) break;
        while (n % p == 0) { push_factor(p); n /= p; }
    }
    factor_rec(n);
}

/* least divisor of l exceeding d, walking the divisor lattice */
static u64 g_d, g_best;

static void walk(int i, u64 cur)
{
    if (i == g_npf) { if (cur > g_d && cur < g_best) g_best = cur; return; }
    u64 v = cur;
    for (int e = 0; e <= g_pe[i]; e++) { walk(i + 1, v); if (e < g_pe[i]) v *= g_pf[i]; }
}

/* decompose a with successor b; writes k, L, d.  returns 1 if decomposable. */
static int decomp(u64 a, u64 b, u64 *ok, u64 *oL, u64 *od)
{
    u64 d = b - a;
    *od = d;
    if (a <= 2 * d) { *ok = 0; *oL = 0; return 0; }
    u64 l = a - d;
    factor(l);
    g_d = d; g_best = ~0ULL;
    walk(0, 1);
    *ok = g_best; *oL = l / g_best;
    return 1;
}

/* ------------------------------------------------------------------ */
/* sequence generators — each fills t[0..cnt-1] with strictly           */
/* increasing terms; one extra term is always generated so that the     */
/* last requested term has a gap.                                       */
/* ------------------------------------------------------------------ */

static void gen_naturals(u64 *t, long cnt)     /* A000027 */
{ for (long i = 0; i < cnt; i++) t[i] = (u64)(i + 1); }

static void gen_primes(u64 *t, long cnt)       /* A000040 */
{
    /* upper bound for the cnt-th prime */
    double n = (double)cnt;
    u64 lim = (n < 6) ? 20 : (u64)(n * (log(n) + log(log(n))) + 100);
    u8 *c = calloc(lim + 1, 1);
    long k = 0;
    for (u64 i = 2; i <= lim && k < cnt; i++) {
        if (!c[i]) { t[k++] = i; for (u64 j = i * i; j <= lim; j += i) c[j] = 1; }
    }
    free(c);
    if (k < cnt) { fprintf(stderr, "prime bound too small\n"); exit(1); }
}

static void gen_triangular(u64 *t, long cnt)   /* A000217, from a(0) = 0 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (n + 1) / 2; } }

static void gen_palindromes(u64 *t, long cnt)  /* A002113, from a(1) = 0 */
{
    long k = 0;
    t[k++] = 0;
    for (int m = 1; k < cnt; m++) {
        int h = (m + 1) / 2;
        u64 lo = 1; for (int i = 1; i < h; i++) lo *= 10;
        u64 hi = lo * 10;
        for (u64 half = lo; half < hi && k < cnt; half++) {
            char buf[24]; int n = snprintf(buf, sizeof buf, "%llu", (unsigned long long)half);
            char pal[48]; int q = 0;
            for (int i = 0; i < n; i++) pal[q++] = buf[i];
            for (int i = (m & 1) ? n - 2 : n - 1; i >= 0; i--) pal[q++] = buf[i];
            pal[q] = 0;
            t[k++] = strtoull(pal, NULL, 10);
        }
    }
}

static void gen_ulam(u64 *t, long cnt)         /* A002858 */
{
    /* a(n) ~ 13.5 n for large n; take a generous limit and grow if needed */
    u64 lim = (u64)(14.5 * cnt + 1000);
    for (;;) {
        size_t W = lim / 64 + 2;
        u64 *ulam = calloc(W, 8), *one = calloc(W, 8), *many = calloc(W, 8);
        if (!ulam || !one || !many) { fprintf(stderr, "ulam alloc failed\n"); exit(1); }
        long k = 0;
        u64 frontier = 0;
        int overflow = 0;
        for (;;) {
            u64 u;
            if (k == 0) u = 1;
            else if (k == 1) u = 2;
            else {
                /* smallest index > frontier carrying exactly one representation */
                u64 i = frontier + 1;
                while (i <= lim && !((one[i >> 6] >> (i & 63)) & 1)) i++;
                if (i > lim) { overflow = 1; break; }
                u = i;
            }
            /* fold in the sums u + x for every earlier Ulam number x */
            u64 w = u >> 6, b = u & 63;
            for (size_t i = w; i < W; i++) {
                u64 nh;
                if (b == 0) nh = ulam[i - w];
                else nh = (ulam[i - w] << b) | (i - w >= 1 ? (ulam[i - w - 1] >> (64 - b)) : 0);
                if (!nh) continue;
                u64 o = one[i], m = many[i];
                many[i] = m | (o & nh);
                one[i]  = (o & ~nh) | (nh & ~o & ~m);
            }
            ulam[u >> 6] |= 1ULL << (u & 63);
            frontier = u;
            t[k++] = u;
            if (k == cnt) break;
        }
        free(ulam); free(one); free(many);
        if (!overflow) return;
        lim *= 2;
    }
}

static int digitsum(u64 n) { int s = 0; while (n) { s += n % 10; n /= 10; } return s; }

static void gen_harshad(u64 *t, long cnt)      /* A005349 */
{
    long k = 0;
    for (u64 n = 1; k < cnt; n++) if (n % (u64)digitsum(n) == 0) t[k++] = n;
}

static void gen_skiptake(u64 *t, long cnt)     /* A004202 = { m^2 + j : 1 <= j <= m } */
{
    long k = 0;
    for (u64 m = 1; k < cnt; m++)
        for (u64 j = 1; j <= m && k < cnt; j++) t[k++] = m * m + j;
}

static void gen_flavius(u64 *t, long cnt)      /* A000960 */
{
    /* The n-th survivor is found by inverting the sieve.  After step k the
     * surviving position p maps to p - floor(p/(k+1)); inverting, a term at
     * position q after step k sat at position q + floor((q-1)/k) before it.
     * Steps k >= n never move position n, so n-1 inversions suffice.        */
    for (long i = 0; i < cnt; i++) {
        u64 p = (u64)(i + 1);
        for (u64 k = (u64)i; k >= 1; k--) p += (p - 1) / k;
        t[i] = p;
    }
}

static void gen_sorting(u64 *t, long cnt)      /* A001855, from index 1 */
{
    u64 a = 0; long k = 0;
    t[k++] = 0;
    for (u64 n = 2; k < cnt; n++) {
        int c = 0; u64 v = 1; while (v < n) { v <<= 1; c++; }   /* ceil(log2 n) */
        a += c; t[k++] = a;
    }
}

static void gen_digitadd(u64 *t, long cnt)     /* A004207, strictly increasing part */
{
    u64 a = 1; long k = 0;
    while (k < cnt) { t[k++] = a; a += (u64)digitsum(a); }
}

/* ---------------------------- the ten added 24 Sep 2026 ------------- */

static void gen_odd(u64 *t, long cnt)          /* A005408, from a(0) = 1 */
{ for (long i = 0; i < cnt; i++) t[i] = 2 * (u64)i + 1; }

static void gen_non3(u64 *t, long cnt)         /* A001651 */
{ long k = 0; for (u64 n = 1; k < cnt; n++) if (n % 3) t[k++] = n; }

static void gen_twinlo(u64 *t, long cnt)       /* A001359 */
{
    for (u64 lim = 1ULL << 26; ; lim *= 2) {
        u8 *c = calloc(lim + 3, 1);
        if (!c) { fprintf(stderr, "twin alloc failed\n"); exit(1); }
        for (u64 i = 2; i * i <= lim + 2; i++)
            if (!c[i]) for (u64 j = i * i; j <= lim + 2; j += i) c[j] = 1;
        long k = 0;
        for (u64 p = 3; p <= lim && k < cnt; p += 2) if (!c[p] && !c[p + 2]) t[k++] = p;
        free(c);
        if (k == cnt) return;
    }
}

static void gen_primesq(u64 *t, long cnt)      /* A001248 */
{ gen_primes(t, cnt); for (long i = 0; i < cnt; i++) t[i] *= t[i]; }

static void gen_semiprimes(u64 *t, long cnt)   /* A001358: Omega(n) = 2 */
{
    long k = 0;
    for (u64 n = 4; k < cnt; n++) {
        if (n > SPF_MAX) { fprintf(stderr, "semiprimes past spf table\n"); exit(1); }
        int om = 0; u64 m = n;
        while (m > 1 && om < 3) { m /= spf[m]; om++; }
        if (om == 2 && m == 1) t[k++] = n;
    }
}

static int squarefree_small(u64 b)             /* b <= SPF_MAX */
{
    while (b > 1) { u64 p = spf[b]; b /= p; if (b % p == 0) return 0; }
    return 1;
}

static int cmp_u64(const void *x, const void *y)
{ u64 a = *(const u64 *)x, b = *(const u64 *)y; return (a > b) - (a < b); }

static void gen_powerful(u64 *t, long cnt)     /* A001694 = { a^2 b^3 : b squarefree }, unique form */
{
    for (u64 lim = 10000000000ULL; ; lim *= 2) {
        size_t cap = (size_t)(3.0 * sqrt((double)lim)) + 1000, m = 0;
        u64 *v = malloc(cap * sizeof(u64));
        if (!v) { fprintf(stderr, "powerful alloc failed\n"); exit(1); }
        for (u64 b = 1; b * b * b <= lim; b++) {
            if (!squarefree_small(b)) continue;
            u64 b3 = b * b * b;
            for (u64 a = 1; a * a * b3 <= lim; a++) {
                if (m == cap) { fprintf(stderr, "powerful cap\n"); exit(1); }
                v[m++] = a * a * b3;
            }
        }
        if ((long)m >= cnt) {
            qsort(v, m, sizeof(u64), cmp_u64);
            for (size_t i = 1; i < m; i++)
                if (v[i] == v[i - 1]) { fprintf(stderr, "powerful: duplicate %llu\n", (unsigned long long)v[i]); exit(1); }
            memcpy(t, v, (size_t)cnt * sizeof(u64));
            free(v);
            return;
        }
        free(v);
    }
}

static void gen_twosq(u64 *t, long cnt)        /* A001481, from a(1) = 0 */
{
    /* n = x^2 + y^2  iff  every prime = 3 (mod 4) divides n to an even power */
    long k = 0;
    t[k++] = 0;
    for (u64 n = 1; k < cnt; n++) {
        if (n > SPF_MAX) { fprintf(stderr, "twosq past spf table\n"); exit(1); }
        u64 m = n; int ok = 1;
        while (m > 1) {
            u64 p = spf[m]; int e = 0;
            while (m % p == 0) { m /= p; e++; }
            if ((p & 3) == 3 && (e & 1)) { ok = 0; break; }
        }
        if (ok) t[k++] = n;
    }
}

static void gen_lucky(u64 *t, long cnt)        /* A000959 */
{
    /* Inverting the sieve, as for A000960.  Start from the odd numbers.  The
     * step with m = t[j] (j >= 1: 3, 7, 9, ...) deletes every m-th survivor, so
     * a survivor at position q after it sat at q + floor((q-1)/(m-1)) before
     * it.  Steps with m > n never move position n, and every step value m <= n
     * is an earlier term, so the terms can be produced in order.             */
    long jmax = 0;
    for (long i = 0; i < cnt; i++) {
        u64 n = (u64)(i + 1);
        while (jmax + 1 < i && t[jmax + 1] <= n) jmax++;
        u64 q = n;
        for (long j = jmax; j >= 1; j--) q += (q - 1) / (t[j] - 1);
        t[i] = 2 * q - 1;
    }
}

static void gen_totsum(u64 *t, long cnt)       /* A002088, from a(0) = 0 */
{
    t[0] = 0;
    for (long i = 1; i < cnt; i++) {
        u64 n = (u64)i, m = n, ph = n;
        if (n > SPF_MAX) { fprintf(stderr, "totsum past spf table\n"); exit(1); }
        while (m > 1) { u64 p = spf[m]; ph = ph / p * (p - 1); while (m % p == 0) m /= p; }
        t[i] = t[i - 1] + ph;
    }
}

static u64 isqrt_u64(u64 x)
{
    u64 r = (u64)sqrtl((long double)x);
    while (r * r > x) r--;
    while ((r + 1) * (r + 1) <= x) r++;
    return r;
}

static void gen_wythoff(u64 *t, long cnt)      /* A000201: floor(n phi) = floor((n + isqrt(5 n^2)) / 2) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = (n + isqrt_u64(5 * n * n)) / 2; } }

/* ---------------------------- the thirty added 25 Sep 2026 ----------- */
/* The rest of the 50-sequence census (decompwlj_atlas_50_record_2026_08_28). */

static u8 *composite_flags(u64 lim)            /* c[i] = 1 iff i is composite (i >= 2) */
{
    u8 *c = calloc(lim + 3, 1);
    if (!c) { fprintf(stderr, "sieve alloc failed\n"); exit(1); }
    for (u64 i = 2; i * i <= lim + 2; i++)
        if (!c[i]) for (u64 j = i * i; j <= lim + 2; j += i) c[j] = 1;
    return c;
}

static void need_spf(u64 n, const char *who)
{ if (n > SPF_MAX) { fprintf(stderr, "%s: past the spf table\n", who); exit(1); } }

static int omega_big(u64 n)                    /* Omega(n), n <= SPF_MAX */
{ int c = 0; while (n > 1) { n /= spf[n]; c++; } return c; }

static int is_prime_power(u64 n)               /* n >= 2 */
{ u64 p = spf[n]; while (n % p == 0) n /= p; return n == 1; }

static u64 sigma_small(u64 n)                  /* sum of divisors, n <= SPF_MAX */
{
    u64 s = 1;
    while (n > 1) { u64 p = spf[n], pk = 1, t = 1; while (n % p == 0) { n /= p; pk *= p; t += pk; } s *= t; }
    return s;
}

static u64 tau_small(u64 n)                    /* number of divisors */
{
    u64 s = 1;
    while (n > 1) { u64 p = spf[n]; int e = 0; while (n % p == 0) { n /= p; e++; } s *= (u64)(e + 1); }
    return s;
}

static void gen_even(u64 *t, long cnt)         /* A005843, from a(0) = 0 */
{ for (long i = 0; i < cnt; i++) t[i] = 2 * (u64)i; }

static void gen_composites(u64 *t, long cnt)   /* A002808 */
{ long k = 0; for (u64 n = 4; k < cnt; n++) { need_spf(n, "composites"); if (spf[n] != n) t[k++] = n; } }

static void gen_nonsquares(u64 *t, long cnt)   /* A000037 */
{
    long k = 0; u64 r = 1;
    for (u64 n = 2; k < cnt; n++) { while ((r + 1) * (r + 1) <= n) r++; if (r * r != n) t[k++] = n; }
}

static void gen_primes3mod4(u64 *t, long cnt)  /* A002145 */
{
    for (u64 lim = 1ULL << 22; ; lim *= 2) {
        u8 *c = composite_flags(lim); long k = 0;
        for (u64 p = 3; p <= lim && k < cnt; p += 4) if (!c[p]) t[k++] = p;
        free(c); if (k == cnt) return;
    }
}

static void gen_twinhi(u64 *t, long cnt)       /* A006512 */
{
    for (u64 lim = 1ULL << 25; ; lim *= 2) {
        u8 *c = composite_flags(lim); long k = 0;
        for (u64 p = 3; p + 2 <= lim && k < cnt; p += 2) if (!c[p] && !c[p + 2]) t[k++] = p + 2;
        free(c); if (k == cnt) return;
    }
}

static void gen_primepowers(u64 *t, long cnt)  /* A000961: 1 and p^e, e >= 1 */
{
    long k = 0; t[k++] = 1;
    for (u64 n = 2; k < cnt; n++) { need_spf(n, "primepowers"); if (is_prime_power(n)) t[k++] = n; }
}

static void gen_notpp(u64 *t, long cnt)        /* A024619: n >= 2 with at least two distinct prime factors */
{ long k = 0; for (u64 n = 2; k < cnt; n++) { need_spf(n, "notpp"); if (!is_prime_power(n)) t[k++] = n; } }

static void gen_pplus1(u64 *t, long cnt)       /* A008864 */
{ gen_primes(t, cnt); for (long i = 0; i < cnt; i++) t[i] += 1; }

static void gen_squarefree(u64 *t, long cnt)   /* A005117 */
{ long k = 0; for (u64 n = 1; k < cnt; n++) { need_spf(n, "squarefree"); if (squarefree_small(n)) t[k++] = n; } }

static void gen_sfsemiprimes(u64 *t, long cnt) /* A006881: pq, p < q */
{
    long k = 0;
    for (u64 n = 6; k < cnt; n++) { need_spf(n, "sfsemiprimes"); if (omega_big(n) == 2 && squarefree_small(n)) t[k++] = n; }
}

static void gen_oddomega(u64 *t, long cnt)     /* A026424: Omega(n) odd */
{ long k = 0; for (u64 n = 2; k < cnt; n++) { need_spf(n, "oddomega"); if (omega_big(n) & 1) t[k++] = n; } }

static void gen_abundant(u64 *t, long cnt)     /* A005101: sigma(n) > 2n */
{ long k = 0; for (u64 n = 1; k < cnt; n++) { need_spf(n, "abundant"); if (sigma_small(n) > 2 * n) t[k++] = n; } }

static void gen_practical(u64 *t, long cnt)    /* A005153, Stewart-Sierpinski criterion */
{
    /* n = p1^e1 ... pk^ek, p1 < ... < pk, is practical iff p1 = 2 (or n = 1) and
       p_i <= 1 + sigma(p1^e1 ... p_{i-1}^e_{i-1}) for every i >= 2.            */
    long k = 0; t[k++] = 1;
    for (u64 n = 2; k < cnt; n += 2) {
        need_spf(n, "practical");
        u64 m = n, s = 1; int ok = 1;
        while (m > 1) {
            u64 p = spf[m], pk = 1, sp = 1;
            if (p > 1 + s) { ok = 0; break; }
            while (m % p == 0) { m /= p; pk *= p; sp += pk; }
            s *= sp;
        }
        if (ok) t[k++] = n;
    }
}

static void gen_squares(u64 *t, long cnt)      /* A000290, from a(0) = 0 */
{ for (long i = 0; i < cnt; i++) t[i] = (u64)i * (u64)i; }

static void gen_pentagonal(u64 *t, long cnt)   /* A000326, from a(0) = 0 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (3 * n - 1) / 2; } }

static void gen_hexagonal(u64 *t, long cnt)    /* A000384, from a(0) = 0 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (2 * n - 1); } }

static void gen_oblong(u64 *t, long cnt)       /* A002378, from a(0) = 0 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (n + 1); } }

static void gen_quartersq(u64 *t, long cnt)    /* A002620 from a(1) = 0 (a(0) = a(1) = 0) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = n * n / 4; } }

static void gen_cubes(u64 *t, long cnt)        /* A000578, from a(0) = 0 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * n * n; } }

static void gen_tetrahedral(u64 *t, long cnt)  /* A000292, from a(0) = 0 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (n + 1) * (n + 2) / 6; } }

static void gen_pyramidal(u64 *t, long cnt)    /* A000330, from a(0) = 0 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (n + 1) * (2 * n + 1) / 6; } }

static void gen_trisq(u64 *t, long cnt)        /* A005214: triangular numbers and squares, merged, no 0 */
{
    u64 j = 1, m = 1; long k = 0;
    while (k < cnt) {
        u64 tri = j * (j + 1) / 2, sq = m * m;
        if (tri < sq) { t[k++] = tri; j++; }
        else if (sq < tri) { t[k++] = sq; m++; }
        else { t[k++] = sq; j++; m++; }
    }
}

static void gen_floorsqrt(u64 *t, long cnt)    /* A006446: floor(sqrt k) | k  =  { m^2, m^2+m, m^2+2m } */
{
    long k = 0;
    for (u64 m = 1; k < cnt; m++)
        for (u64 j = 0; j < 3 && k < cnt; j++) t[k++] = m * m + j * m;
}

static void gen_zeroless(u64 *t, long cnt)     /* A052382 */
{
    long k = 0;
    for (u64 n = 1; k < cnt; n++) { u64 v = n; int ok = 1; while (v) { if (v % 10 == 0) { ok = 0; break; } v /= 10; } if (ok) t[k++] = n; }
}

static void gen_self(u64 *t, long cnt)         /* A003052: not of the form m + digitsum(m) */
{
    for (u64 lim = 1ULL << 20; ; lim *= 2) {
        u8 *hit = calloc(lim + 1, 1); long k = 0;
        if (!hit) { fprintf(stderr, "self alloc failed\n"); exit(1); }
        for (u64 m = 1; m <= lim; m++) { u64 g = m + (u64)digitsum(m); if (g <= lim) hit[g] = 1; }
        for (u64 n = 1; n <= lim && k < cnt; n++) if (!hit[n]) t[k++] = n;
        free(hit); if (k == cnt) return;
    }
}

static void gen_fibbinary(u64 *t, long cnt)    /* A003714, from a(0) = 0 */
{ long k = 0; for (u64 n = 0; k < cnt; n++) if (!(n & (n >> 1))) t[k++] = n; }

static void gen_binary(u64 *t, long cnt)       /* A007088, from a(0) = 0: n in base 2, read in base 10 */
{
    for (long i = 0; i < cnt; i++) {
        u64 n = (u64)i, v = 0, p = 1;
        while (n) { v += (n & 1) * p; p *= 10; n >>= 1; }
        t[i] = v;
    }
}

static void gen_primesum(u64 *t, long cnt)     /* A007504, from a(0) = 0 */
{
    gen_primes(t + 1, cnt - 1);
    t[0] = 0;
    for (long i = 1; i < cnt; i++) t[i] += t[i - 1];
}

static void gen_divsum(u64 *t, long cnt)       /* A006218, from a(0) = 0: sum of d(k), k <= n */
{ t[0] = 0; for (long i = 1; i < cnt; i++) { need_spf((u64)i, "divsum"); t[i] = t[i - 1] + tau_small((u64)i); } }

static void gen_sigmasum(u64 *t, long cnt)     /* A024916, from a(1) = 1: sum of sigma(k), k <= n */
{ u64 s = 0; for (long i = 0; i < cnt; i++) { need_spf((u64)(i + 1), "sigmasum"); s += sigma_small((u64)(i + 1)); t[i] = s; } }

/* ---------------------------- the fifty added 26 Sep 2026 ----------- */
/* The second fifty, for the 100-sequence GitHub Pages site.            */

static int popcount_u64(u64 n) { int c = 0; while (n) { n &= n - 1; c++; } return c; }

static void gen_odious(u64 *t, long cnt)       /* A000069, from a(1) = 1: odd binary weight */
{ long k = 0; for (u64 n = 1; k < cnt; n++) if (popcount_u64(n) & 1) t[k++] = n; }

static void gen_evil(u64 *t, long cnt)         /* A001969, from a(1) = 0: even binary weight */
{ long k = 0; for (u64 n = 0; k < cnt; n++) if (!(popcount_u64(n) & 1)) t[k++] = n; }

static void gen_heptagonal(u64 *t, long cnt)   /* A000566, from a(0) = 0 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (5 * n - 3) / 2; } }

static void gen_octagonal(u64 *t, long cnt)    /* A000567, from a(0) = 0 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (3 * n - 2); } }

static void gen_censquare(u64 *t, long cnt)    /* A001844, from a(0) = 1 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = 2 * n * (n + 1) + 1; } }

static void gen_cenhex(u64 *t, long cnt)       /* A003215, from a(0) = 1 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = 3 * n * (n + 1) + 1; } }

static void gen_centri(u64 *t, long cnt)       /* A005448, from a(1) = 1 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = 3 * n * (n - 1) / 2 + 1; } }

static void gen_lazycaterer(u64 *t, long cnt)  /* A000124, from a(0) = 1 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (n + 1) / 2 + 1; } }

static void gen_sqplus1(u64 *t, long cnt)      /* A002522, from a(0) = 1 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * n + 1; } }

static void gen_genpent(u64 *t, long cnt)      /* A001318, from a(0) = 0: m(3m-1)/2, m = 0, 1, -1, 2, -2, ... */
{
    for (long i = 0; i < cnt; i++) {
        u64 j = (u64)(i + 1) / 2;
        t[i] = (i & 1) ? j * (3 * j - 1) / 2 : j * (3 * j + 1) / 2;
    }
}

static void gen_upperwythoff(u64 *t, long cnt) /* A001950: floor(n phi^2) = floor(n phi) + n */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = (n + isqrt_u64(5 * n * n)) / 2 + n; } }

static void gen_beattysqrt2(u64 *t, long cnt)  /* A001951, from a(0) = 0: floor(n sqrt 2) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = isqrt_u64(2 * n * n); } }

static void gen_beattysqrt3(u64 *t, long cnt)  /* A022838: floor(n sqrt 3) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = isqrt_u64(3 * n * n); } }

/* primes p (in order) with pred(p, c), c the composite flags up to lim; the limit
   doubles until cnt terms are found.  pred may look at c[] up to 2 lim + 1.      */
static void gen_prime_filter(u64 *t, long cnt, u64 lim, int (*pred)(u64, const u8 *))
{
    for (; ; lim *= 2) {
        u8 *c = composite_flags(2 * lim + 8); long k = 0;
        for (u64 p = 2; p <= lim && k < cnt; p++) if (!c[p] && pred(p, c)) t[k++] = p;
        free(c); if (k == cnt) return;
    }
}

static int pr(u64 n, const u8 *c) { return n >= 2 && !c[n]; }
static int p_1mod4(u64 p, const u8 *c)  { (void)c; return p % 4 == 1; }
static int p_1mod6(u64 p, const u8 *c)  { (void)c; return p % 6 == 1; }
static int p_5mod6(u64 p, const u8 *c)  { (void)c; return p % 6 == 5; }
static int p_sophie(u64 p, const u8 *c) { return pr(2 * p + 1, c); }
static int p_safe(u64 p, const u8 *c)   { return p > 2 && pr((p - 1) / 2, c); }
static int p_isolated(u64 p, const u8 *c) { return !(p >= 2 && pr(p - 2, c)) && !pr(p + 2, c); }
static int p_twin(u64 p, const u8 *c)   { return (p >= 2 && pr(p - 2, c)) || pr(p + 2, c); }
static int p_cousin(u64 p, const u8 *c) { return pr(p + 4, c); }
static int p_sexy(u64 p, const u8 *c)   { return pr(p + 6, c); }

static void gen_primes1mod4(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 22, p_1mod4); }  /* A002144 */
static void gen_primes1mod6(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 22, p_1mod6); }  /* A002476 */
static void gen_primes5mod6(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 22, p_5mod6); }  /* A007528 */
static void gen_sophie(u64 *t, long cnt)      { gen_prime_filter(t, cnt, 1ULL << 24, p_sophie); } /* A005384 */
static void gen_safe(u64 *t, long cnt)        { gen_prime_filter(t, cnt, 1ULL << 25, p_safe); }   /* A005385 */
static void gen_isolated(u64 *t, long cnt)    { gen_prime_filter(t, cnt, 1ULL << 21, p_isolated); } /* A007510 */
static void gen_twins(u64 *t, long cnt)       { gen_prime_filter(t, cnt, 1ULL << 23, p_twin); }   /* A001097 */
static void gen_cousin(u64 *t, long cnt)      { gen_prime_filter(t, cnt, 1ULL << 24, p_cousin); } /* A023200 */
static void gen_sexy(u64 *t, long cnt)        { gen_prime_filter(t, cnt, 1ULL << 23, p_sexy); }   /* A023201 */

static void gen_primeidx(u64 *t, long cnt)     /* A006450: prime(prime(n)) */
{
    u64 *p = malloc(sizeof(u64) * 1500000);
    gen_primes(p, 1500000);
    for (long i = 0; i < cnt; i++) {
        if (p[i] > 1500000) { fprintf(stderr, "primeidx: prime table too short\n"); exit(1); }
        t[i] = p[p[i] - 1];
    }
    free(p);
}

static u64 reverse10(u64 n) { u64 r = 0; while (n) { r = r * 10 + n % 10; n /= 10; } return r; }

static void gen_emirps(u64 *t, long cnt)       /* A006567: primes whose reversal is a different prime */
{
    /* the reversal has no more digits than p, so a sieve to a power of ten suffices */
    for (u64 lim = 10000000ULL; ; lim *= 10) {
        u8 *c = composite_flags(lim); long k = 0;
        for (u64 p = 2; p < lim && k < cnt; p++) {
            if (c[p]) continue;
            u64 r = reverse10(p);
            if (r != p && r >= 2 && !c[r]) t[k++] = p;
        }
        free(c); if (k == cnt) return;
    }
}

static void gen_primesminus1(u64 *t, long cnt) /* A006093 */
{ gen_primes(t, cnt); for (long i = 0; i < cnt; i++) t[i] -= 1; }

static void gen_almost3(u64 *t, long cnt)      /* A014612: Omega(n) = 3 */
{ long k = 0; for (u64 n = 8; k < cnt; n++) { need_spf(n, "almost3"); if (omega_big(n) == 3) t[k++] = n; } }

static void gen_almost4(u64 *t, long cnt)      /* A014613: Omega(n) = 4 */
{ long k = 0; for (u64 n = 16; k < cnt; n++) { need_spf(n, "almost4"); if (omega_big(n) == 4) t[k++] = n; } }

static void gen_sphenic(u64 *t, long cnt)      /* A007304: pqr, p < q < r */
{ long k = 0; for (u64 n = 30; k < cnt; n++) { need_spf(n, "sphenic"); if (omega_big(n) == 3 && squarefree_small(n)) t[k++] = n; } }

static void gen_nonsqfree(u64 *t, long cnt)    /* A013929 */
{ long k = 0; for (u64 n = 4; k < cnt; n++) { need_spf(n, "nonsqfree"); if (!squarefree_small(n)) t[k++] = n; } }

static void gen_deficient(u64 *t, long cnt)    /* A005100: sigma(n) < 2n */
{ long k = 0; for (u64 n = 1; k < cnt; n++) { need_spf(n, "deficient"); if (sigma_small(n) < 2 * n) t[k++] = n; } }

static int happy(u64 n)
{
    while (n != 1 && n != 4) { u64 s = 0; while (n) { u64 q = n % 10; s += q * q; n /= 10; } n = s; }
    return n == 1;
}

static void gen_happy(u64 *t, long cnt)        /* A007770 */
{ long k = 0; for (u64 n = 1; k < cnt; n++) if (happy(n)) t[k++] = n; }

static void gen_haszero(u64 *t, long cnt)      /* A011540, from a(1) = 0: a digit 0 */
{
    long k = 0; t[k++] = 0;
    for (u64 n = 1; k < cnt; n++) { u64 v = n; while (v) { if (v % 10 == 0) { t[k++] = n; break; } v /= 10; } }
}

static void gen_odddigits(u64 *t, long cnt)    /* A014261: every digit odd */
{
    /* i + 1 in bijective base 5, digits 1..5 written as 1, 3, 5, 7, 9: this keeps the order */
    for (long i = 0; i < cnt; i++) {
        u64 m = (u64)(i + 1), v = 0, p = 1;
        while (m) { u64 r = (m - 1) % 5; v += (2 * r + 1) * p; p *= 10; m = (m - 1) / 5; }
        t[i] = v;
    }
}

/* n > 0 is kept iff every prime p with bad(p) divides n to an even power */
static void gen_evenpower(u64 *t, long cnt, int from0, int (*bad)(u64))
{
    long k = 0;
    if (from0) t[k++] = 0;
    for (u64 n = 1; k < cnt; n++) {
        need_spf(n, "evenpower");
        u64 m = n; int ok = 1;
        while (m > 1) {
            u64 p = spf[m]; int e = 0;
            while (m % p == 0) { m /= p; e++; }
            if (bad(p) && (e & 1)) { ok = 0; break; }
        }
        if (ok) t[k++] = n;
    }
}
static int bad_loesch(u64 p) { return p % 3 == 2; }
static int bad_x2p2y2(u64 p) { return p % 8 == 5 || p % 8 == 7; }

static void gen_loeschian(u64 *t, long cnt) { gen_evenpower(t, cnt, 1, bad_loesch); }   /* A003136 = x^2 + xy + y^2 */
static void gen_x2p2y2(u64 *t, long cnt)    { gen_evenpower(t, cnt, 1, bad_x2p2y2); }   /* A002479 = x^2 + 2 y^2 */

static void gen_twopossq(u64 *t, long cnt)     /* A000404: x^2 + y^2, x, y >= 1 */
{
    for (u64 lim = 1ULL << 18; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 x = 1; 2 * x * x <= lim; x++)
            for (u64 y = x; x * x + y * y <= lim; y++) h[x * x + y * y] = 1;
        for (u64 n = 1; n <= lim && k < cnt; n++) if (h[n]) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}

static int is_4a8b7(u64 n)                     /* n = 4^a (8b + 7) */
{ if (!n) return 0; while (n % 4 == 0) n /= 4; return n % 8 == 7; }

static void gen_threesq(u64 *t, long cnt)      /* A000378, from a(1) = 0: sums of three squares */
{ long k = 0; for (u64 n = 0; k < cnt; n++) if (!is_4a8b7(n)) t[k++] = n; }

static void gen_notthreesq(u64 *t, long cnt)   /* A004215: 4^a (8b + 7) */
{ long k = 0; for (u64 n = 1; k < cnt; n++) if (is_4a8b7(n)) t[k++] = n; }

static void gen_coprime30(u64 *t, long cnt)    /* A007775: gcd(n, 30) = 1 */
{ long k = 0; for (u64 n = 1; k < cnt; n++) if (n % 2 && n % 3 && n % 5) t[k++] = n; }

static void gen_coprime6(u64 *t, long cnt)     /* A007310: n = 1 or 5 (mod 6) */
{ long k = 0; for (u64 n = 1; k < cnt; n++) if (n % 2 && n % 3) t[k++] = n; }

static void gen_mult3(u64 *t, long cnt)        /* A008585, from a(0) = 0 */
{ for (long i = 0; i < cnt; i++) t[i] = 3 * (u64)i; }

static void gen_oddnonprimes(u64 *t, long cnt) /* A014076 */
{ long k = 0; for (u64 n = 1; k < cnt; n += 2) { need_spf(n, "oddnonprimes"); if (n == 1 || spf[n] != n) t[k++] = n; } }

static void gen_cantor(u64 *t, long cnt)       /* A005836, from a(1) = 0: no digit 2 in base 3 */
{
    for (long i = 0; i < cnt; i++) {
        u64 n = (u64)i, v = 0, p = 1;
        while (n) { v += (n & 1) * p; p *= 3; n >>= 1; }
        t[i] = v;
    }
}

static void gen_moser(u64 *t, long cnt)        /* A000695, from a(0) = 0: sums of distinct powers of 4 */
{
    for (long i = 0; i < cnt; i++) {
        u64 n = (u64)i, v = 0, p = 1;
        while (n) { v += (n & 1) * p; p *= 4; n >>= 1; }
        t[i] = v;
    }
}

static void gen_evenzeros(u64 *t, long cnt)    /* A003159: even number of trailing 0 bits */
{
    long k = 0;
    for (u64 n = 1; k < cnt; n++) { int z = 0; u64 v = n; while (!(v & 1)) { v >>= 1; z++; } if (!(z & 1)) t[k++] = n; }
}

static void gen_binpal(u64 *t, long cnt)       /* A006995, from a(1) = 0: binary palindromes */
{
    long k = 0; t[k++] = 0;
    for (u64 n = 1; k < cnt; n += 2) {         /* palindromes > 0 are odd */
        u64 r = 0, v = n; while (v) { r = (r << 1) | (v & 1); v >>= 1; }
        if (r == n) t[k++] = n;
    }
}

static void gen_ternary(u64 *t, long cnt)      /* A007089, from a(0) = 0: n in base 3, read in base 10 */
{
    for (long i = 0; i < cnt; i++) {
        u64 n = (u64)i, v = 0, p = 1;
        while (n) { v += (n % 3) * p; p *= 10; n /= 3; }
        t[i] = v;
    }
}

static void gen_pentpyr(u64 *t, long cnt)      /* A002411, from a(0) = 0: n^2 (n+1)/2 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * n * (n + 1) / 2; } }

static void gen_octahedral(u64 *t, long cnt)   /* A005900, from a(0) = 0: n (2 n^2 + 1)/3 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (2 * n * n + 1) / 3; } }

static int powerful_small(u64 n)
{
    while (n > 1) { u64 p = spf[n]; int e = 0; while (n % p == 0) { n /= p; e++; } if (e == 1) return 0; }
    return 1;
}

static void gen_weak(u64 *t, long cnt)         /* A052485: not powerful */
{ long k = 0; for (u64 n = 2; k < cnt; n++) { need_spf(n, "weak"); if (!powerful_small(n)) t[k++] = n; } }

/* ---------------------------- the hundred added 26 Sep 2026 ---------- */
/* The third and fourth fifty, for 200 sequences.                       */

/* polynomial and figurate */
static void gen_nnp3(u64 *t, long cnt)         /* A000096, from a(0) = 0: n(n+3)/2 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (n + 3) / 2; } }
static void gen_nnp2(u64 *t, long cnt)         /* A005563, from a(0) = 0: n(n+2) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (n + 2); } }
static void gen_twicesq(u64 *t, long cnt)      /* A001105, from a(0) = 0: 2n^2 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = 2 * n * n; } }
static void gen_oddsq(u64 *t, long cnt)        /* A016754, from a(0) = 1: (2n+1)^2 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = (2 * n + 1) * (2 * n + 1); } }
static void gen_genoct(u64 *t, long cnt)       /* A001082, from a(1) = 0: m(3m-2), m = 0, 1, -1, 2, -2, ... */
{
    for (long i = 0; i < cnt; i++) {
        u64 j = (u64)(i + 1) / 2;
        t[i] = (i & 1) ? j * (3 * j - 2) : j * (3 * j + 2);
    }
}
static void gen_cenpent(u64 *t, long cnt)      /* A005891, from a(0) = 1: (5n^2 + 5n + 2)/2 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = (5 * n * n + 5 * n + 2) / 2; } }
static void gen_cenocta(u64 *t, long cnt)      /* A001845, from a(0) = 1: (2n+1)(2n^2+2n+3)/3 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = (2 * n + 1) * (2 * n * n + 2 * n + 3) / 3; } }
static void gen_cencube(u64 *t, long cnt)      /* A005898, from a(0) = 1: n^3 + (n+1)^3 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * n * n + (n + 1) * (n + 1) * (n + 1); } }
static void gen_stella(u64 *t, long cnt)       /* A007588, from a(0) = 0: n(2n^2 - 1) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n ? n * (2 * n * n - 1) : 0; } }
static void gen_oddsqsum(u64 *t, long cnt)     /* A000447, from a(0) = 0: n(2n-1)(2n+1)/3 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n ? n * (2 * n - 1) * (2 * n + 1) / 3 : 0; } }
static void gen_hexpyr(u64 *t, long cnt)       /* A002412, from a(0) = 0: n(n+1)(4n-1)/6 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n ? n * (n + 1) * (4 * n - 1) / 6 : 0; } }
static void gen_cake(u64 *t, long cnt)         /* A000125, from a(0) = 1: (n^3 + 5n + 6)/6 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = (n * n * n + 5 * n + 6) / 6; } }
static void gen_magic(u64 *t, long cnt)        /* A006003, from a(0) = 0: n(n^2 + 1)/2 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (n * n + 1) / 2; } }
static void gen_sechex(u64 *t, long cnt)       /* A014105, from a(0) = 0: n(2n+1) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (2 * n + 1); } }
static void gen_secpent(u64 *t, long cnt)      /* A005449, from a(0) = 0: n(3n+1)/2 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = n * (3 * n + 1) / 2; } }
static void gen_trimatch(u64 *t, long cnt)     /* A045943, from a(0) = 0: 3n(n+1)/2 */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = 3 * n * (n + 1) / 2; } }
static void gen_thirdsq(u64 *t, long cnt)      /* A000212 from a(1) = 0 (a(0) = a(1) = 0): floor(n^2/3) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = n * n / 3; } }
static void gen_halfsq(u64 *t, long cnt)       /* A007590 from a(1) = 0 (a(0) = a(1) = 0): floor(n^2/2) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = n * n / 2; } }
static void gen_pow32(u64 *t, long cnt)        /* A000093, from a(0) = 0: floor(n^(3/2)) = isqrt(n^3) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i; t[i] = isqrt_u64(n * n * n); } }

/* primes */
static int p_1mod8(u64 p, const u8 *c) { (void)c; return p % 8 == 1; }
static int p_3mod8(u64 p, const u8 *c) { (void)c; return p % 8 == 3; }
static int p_5mod8(u64 p, const u8 *c) { (void)c; return p % 8 == 5; }
static int p_7mod8(u64 p, const u8 *c) { (void)c; return p % 8 == 7; }
static int p_end1(u64 p, const u8 *c)  { (void)c; return p % 10 == 1; }
static int p_end3(u64 p, const u8 *c)  { (void)c; return p % 10 == 3; }
static int p_end7(u64 p, const u8 *c)  { (void)c; return p % 10 == 7; }
static int p_end9(u64 p, const u8 *c)  { (void)c; return p % 10 == 9; }
static int p_2pm1(u64 p, const u8 *c)  { return pr(2 * p - 1, c); }
static int p_halfp1(u64 p, const u8 *c) { return p > 2 && pr((p + 1) / 2, c); }
static int p_plus8(u64 p, const u8 *c) { return pr(p + 8, c); }
static int p_trip26(u64 p, const u8 *c) { return pr(p + 2, c) && pr(p + 6, c); }
static int p_trip46(u64 p, const u8 *c) { return pr(p + 4, c) && pr(p + 6, c); }
static int p_gap6(u64 p, const u8 *c)  { return pr(p + 6, c) && !pr(p + 2, c) && !pr(p + 4, c); }

static void gen_p1mod8(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 23, p_1mod8); }   /* A007519 */
static void gen_p3mod8(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 23, p_3mod8); }   /* A007520 */
static void gen_p5mod8(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 23, p_5mod8); }   /* A007521 */
static void gen_p7mod8(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 23, p_7mod8); }   /* A007522 */
static void gen_pend1(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 23, p_end1); }    /* A030430 */
static void gen_pend3(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 23, p_end3); }    /* A030431 */
static void gen_pend7(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 23, p_end7); }    /* A030432 */
static void gen_pend9(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 23, p_end9); }    /* A030433 */
static void gen_p2pm1(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 24, p_2pm1); }    /* A005382 */
static void gen_phalf(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 24, p_halfp1); }  /* A005383 */
static void gen_pplus8(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 24, p_plus8); }   /* A023202 */
static void gen_trip26(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 27, p_trip26); }  /* A022004 */
static void gen_trip46(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 27, p_trip46); }  /* A022005 */
static void gen_gap6(u64 *t, long cnt)   { gen_prime_filter(t, cnt, 1ULL << 24, p_gap6); }    /* A031924 */

/* primes p(n) with 1 < n and a condition on p(n-1), p(n), p(n+1) */
static void gen_prime_nbr(u64 *t, long cnt, int (*keep)(u64, u64, u64))
{
    for (u64 lim = 1ULL << 24; ; lim *= 2) {
        u8 *c = composite_flags(lim); long k = 0;
        u64 a = 2, b = 3;
        for (u64 q = 5; q <= lim && k < cnt; q += 2) {
            if (c[q]) continue;
            if (keep(a, b, q)) t[k++] = b;
            a = b; b = q;
        }
        free(c); if (k == cnt) return;
    }
}
static int k_balanced(u64 a, u64 b, u64 c) { return 2 * b == a + c; }
static int k_strong(u64 a, u64 b, u64 c)   { return 2 * b > a + c; }
static int k_weak(u64 a, u64 b, u64 c)     { return 2 * b < a + c; }
static void gen_balanced(u64 *t, long cnt) { gen_prime_nbr(t, cnt, k_balanced); }   /* A006562 */
static void gen_strong(u64 *t, long cnt)   { gen_prime_nbr(t, cnt, k_strong); }     /* A051634 */
static void gen_weakp(u64 *t, long cnt)    { gen_prime_nbr(t, cnt, k_weak); }       /* A051635 */

static void gen_pprod(u64 *t, long cnt)        /* A006094: p(n) p(n+1) */
{ u64 *p = malloc((cnt + 1) * sizeof(u64)); gen_primes(p, cnt + 1); for (long i = 0; i < cnt; i++) t[i] = p[i] * p[i + 1]; free(p); }
static void gen_npn(u64 *t, long cnt)          /* A033286: n p(n) */
{ gen_primes(t, cnt); for (long i = 0; i < cnt; i++) t[i] *= (u64)(i + 1); }
static void gen_pnplusn(u64 *t, long cnt)      /* A014688: p(n) + n */
{ gen_primes(t, cnt); for (long i = 0; i < cnt; i++) t[i] += (u64)(i + 1); }
static void gen_psum2(u64 *t, long cnt)        /* A001043: p(n) + p(n+1) */
{ u64 *p = malloc((cnt + 1) * sizeof(u64)); gen_primes(p, cnt + 1); for (long i = 0; i < cnt; i++) t[i] = p[i] + p[i + 1]; free(p); }
static void gen_twop(u64 *t, long cnt)         /* A100484: 2p */
{ gen_primes(t, cnt); for (long i = 0; i < cnt; i++) t[i] *= 2; }
static void gen_threep(u64 *t, long cnt)       /* A001748: 3p */
{ gen_primes(t, cnt); for (long i = 0; i < cnt; i++) t[i] *= 3; }

/* multiplicative */
static int omega_distinct(u64 n) { int c = 0; while (n > 1) { u64 p = spf[n]; while (n % p == 0) n /= p; c++; } return c; }
static int cubefree_small(u64 n)
{ while (n > 1) { u64 p = spf[n]; int e = 0; while (n % p == 0) { n /= p; e++; } if (e >= 3) return 0; } return 1; }
static int pq2(u64 n)                          /* n = p q^2, p != q */
{
    u64 p = spf[n]; int e1 = 0; while (n % p == 0) { n /= p; e1++; }
    if (n == 1) return 0;
    u64 q = spf[n]; int e2 = 0; while (n % q == 0) { n /= q; e2++; }
    return n == 1 && e1 + e2 == 3 && e1 != e2;
}

#define SPF_FILTER(name, from, cond) \
static void name(u64 *t, long cnt) \
{ long k = 0; for (u64 n = (from); k < cnt; n++) { need_spf(n, #name); if (cond) t[k++] = n; } }

SPF_FILTER(gen_pq2,        12, pq2(n))                                               /* A054753 */
SPF_FILTER(gen_fourdist,   210, omega_distinct(n) == 4 && squarefree_small(n))       /* A046386 */
SPF_FILTER(gen_almost5,    32, omega_big(n) == 5)                                    /* A014614 */
SPF_FILTER(gen_cubefree,   1, cubefree_small(n))                                     /* A004709 */
SPF_FILTER(gen_noncubefree, 8, !cubefree_small(n))                                   /* A046099 */
SPF_FILTER(gen_omega2,     6, omega_distinct(n) == 2)                                /* A007774 */
SPF_FILTER(gen_omega3,     30, omega_distinct(n) == 3)                               /* A033992 */
SPF_FILTER(gen_omega4,     210, omega_distinct(n) == 4)                              /* A033993 */
SPF_FILTER(gen_mu1,        1, squarefree_small(n) && !(omega_big(n) & 1))            /* A030229 */
SPF_FILTER(gen_mum1,       2, squarefree_small(n) && (omega_big(n) & 1))             /* A030059 */
SPF_FILTER(gen_evenomega,  1, !(omega_big(n) & 1))                                   /* A028260 */
SPF_FILTER(gen_arith,      1, sigma_small(n) % tau_small(n) == 0)                    /* A003601 */

/* values of Euler's phi.  phi(m) > m / (e^gamma lnln m + 3 / lnln m) for m >= 3
   (Rosser-Schoenfeld), so every value n <= N comes from some m <= M once that
   bound at M exceeds N; the bound increases with m.                          */
static u8 *totient_flags(u64 N)
{
    u64 M = 8 * N;
    double ll = log(log((double)M));
    if ((double)M / (1.7810724179901979 * ll + 3.0 / ll) <= (double)N) { fprintf(stderr, "totient bound\n"); exit(1); }
    need_spf(M, "totients");
    u8 *f = calloc(N + 1, 1);
    for (u64 m = 1; m <= M; m++) {
        u64 x = m, ph = m;
        while (x > 1) { u64 p = spf[x]; ph = ph / p * (p - 1); while (x % p == 0) x /= p; }
        if (ph <= N) f[ph] = 1;
    }
    return f;
}
static void gen_totients(u64 *t, long cnt)     /* A002202 */
{
    for (u64 N = 1ULL << 18; ; N *= 2) {
        u8 *f = totient_flags(N); long k = 0;
        for (u64 n = 1; n <= N && k < cnt; n++) if (f[n]) t[k++] = n;
        free(f); if (k == cnt) return;
    }
}
static void gen_nontotients(u64 *t, long cnt)  /* A007617 */
{
    for (u64 N = 1ULL << 18; ; N *= 2) {
        u8 *f = totient_flags(N); long k = 0;
        for (u64 n = 1; n <= N && k < cnt; n++) if (!f[n]) t[k++] = n;
        free(f); if (k == cnt) return;
    }
}

/* digits */
static void gen_base_read(u64 *t, long cnt, u64 b)   /* n in base b, read in decimal, from a(0) = 0 */
{
    for (long i = 0; i < cnt; i++) {
        u64 n = (u64)i, v = 0, p = 1;
        while (n) { v += (n % b) * p; p *= 10; n /= b; }
        t[i] = v;
    }
}
static void gen_base4(u64 *t, long cnt) { gen_base_read(t, cnt, 4); }   /* A007090 */
static void gen_base5(u64 *t, long cnt) { gen_base_read(t, cnt, 5); }   /* A007091 */
static void gen_base6(u64 *t, long cnt) { gen_base_read(t, cnt, 6); }   /* A007092 */

static int nozero_base(u64 n, u64 b) { while (n) { if (n % b == 0) return 0; n /= b; } return 1; }
static int has_digit(u64 n, int d) { do { if ((int)(n % 10) == d) return 1; n /= 10; } while (n); return 0; }
static int all_even_digits(u64 n) { do { if ((n % 10) & 1) return 0; n /= 10; } while (n); return 1; }
static int ndigits_bin(u64 n) { int c = 0; while (n) { c++; n >>= 1; } return c; }
static int is_palin10(u64 n) { return n == reverse10(n); }

#define N_FILTER(name, from, cond) \
static void name(u64 *t, long cnt) { long k = 0; for (u64 n = (from); k < cnt; n++) if (cond) t[k++] = n; }

N_FILTER(gen_nozero4,  1, nozero_base(n, 4))                                  /* A023705 */
N_FILTER(gen_nozero3,  1, nozero_base(n, 3))                                  /* A032924 */
N_FILTER(gen_has1,     1, has_digit(n, 1))                                    /* A011531 */
N_FILTER(gen_has9,     1, has_digit(n, 9))                                    /* A011539 */
N_FILTER(gen_no1,      0, !has_digit(n, 1))                                   /* A052383 */
N_FILTER(gen_no2,      0, !has_digit(n, 2))                                   /* A052404 */
N_FILTER(gen_evendig,  0, all_even_digits(n))                                 /* A014263 */
N_FILTER(gen_dsprime,  1, is_prime((u64)digitsum(n)))                         /* A028834 */
N_FILTER(gen_dseven,   0, !(digitsum(n) & 1))                                 /* A054683 */
N_FILTER(gen_binbal,   1, 2 * popcount_u64(n) == ndigits_bin(n))              /* A031443 */
N_FILTER(gen_oddval,   1, (__builtin_ctzll(n) & 1))                           /* A036554 */
N_FILTER(gen_unhappy,  1, !happy(n))                                          /* A031177 */
N_FILTER(gen_nonniven, 1, n % (u64)digitsum(n) != 0)                          /* A065877 */
N_FILTER(gen_nonpal,   1, !is_palin10(n))                                     /* A029742 */

/* numbers with monotone digits, in order: by length, then lexicographically */
static long mono_k, mono_cnt; static u64 *mono_t;
static void mono_rec(int left, int prev, u64 v, int up)
{
    if (mono_k >= mono_cnt) return;
    if (!left) { mono_t[mono_k++] = v; return; }
    if (up) for (int d = prev; d <= 9; d++) mono_rec(left - 1, d, v * 10 + d, up);
    else    for (int d = 0; d <= prev; d++) mono_rec(left - 1, d, v * 10 + d, up);
}
static void gen_mono(u64 *t, long cnt, int up)
{
    mono_t = t; mono_cnt = cnt; mono_k = 0;
    t[mono_k++] = 0;
    for (int len = 1; mono_k < cnt; len++)
        for (int d = 1; d <= 9; d++) mono_rec(len - 1, d, (u64)d, up);
}
static void gen_nondec(u64 *t, long cnt) { gen_mono(t, cnt, 1); }   /* A009994, from a(1) = 0 */
static void gen_noninc(u64 *t, long cnt) { gen_mono(t, cnt, 0); }   /* A009996, from a(1) = 0 */

/* Beatty */
static void gen_beatty_1s2(u64 *t, long cnt)   /* A003151: floor(n(1 + sqrt 2)) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = n + isqrt_u64(2 * n * n); } }
static void gen_beatty_1h2(u64 *t, long cnt)   /* A003152: floor(n(1 + 1/sqrt 2)) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = n + isqrt_u64(n * n / 2); } }
static void gen_beatty_2s2(u64 *t, long cnt)   /* A001952: floor(n(2 + sqrt 2)) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = 2 * n + isqrt_u64(2 * n * n); } }
static void gen_beatty_s5(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = isqrt_u64(5 * n * n); } } /* A022839 */
static void gen_beatty_s6(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = isqrt_u64(6 * n * n); } } /* A022840 */
static void gen_beatty_s7(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = isqrt_u64(7 * n * n); } } /* A022841 */
static void gen_beatty_s8(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = isqrt_u64(8 * n * n); } } /* A022842 */

/* squares, powers, complements */
static void gen_twodistsq(u64 *t, long cnt)    /* A004431: x^2 + y^2, 0 < x < y */
{
    for (u64 lim = 1ULL << 18; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 x = 1; 2 * x * x < lim; x++)
            for (u64 y = x + 1; x * x + y * y <= lim; y++) h[x * x + y * y] = 1;
        for (u64 n = 1; n <= lim && k < cnt; n++) if (h[n]) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static void gen_threeposq(u64 *t, long cnt)    /* A000408: x^2 + y^2 + z^2, x, y, z >= 1 */
{
    for (u64 lim = 1ULL << 18; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 x = 1; 3 * x * x <= lim; x++)
            for (u64 y = x; x * x + 2 * y * y <= lim; y++)
                for (u64 z = y; x * x + y * y + z * z <= lim; z++) h[x * x + y * y + z * z] = 1;
        for (u64 n = 1; n <= lim && k < cnt; n++) if (h[n]) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static int twosq_small(u64 n)
{
    while (n > 1) { u64 p = spf[n]; int e = 0; while (n % p == 0) { n /= p; e++; } if ((p & 3) == 3 && (e & 1)) return 0; }
    return 1;
}
SPF_FILTER(gen_nottwosq, 1, !twosq_small(n))                                        /* A022544 */

static int is_triangular(u64 n) { u64 r = isqrt_u64(8 * n + 1); return r * r == 8 * n + 1; }
static int is_cube(u64 n) { u64 r = (u64)cbrtl((long double)n); while (r * r * r > n) r--; while ((r + 1) * (r + 1) * (r + 1) <= n) r++; return r * r * r == n; }
N_FILTER(gen_nontri,  1, !is_triangular(n))                                   /* A014132 */
N_FILTER(gen_noncube, 2, !is_cube(n))                                         /* A007412 */

static int is_square(u64 x) { u64 r = isqrt_u64(x); return r * r == x; }
static int is_fib(u64 n) { return is_square(5 * n * n + 4) || (n && is_square(5 * n * n - 4)); }
N_FILTER(gen_nonfib,  4, !is_fib(n))                                          /* A001690 */

static u8 *perfect_power_flags(u64 lim)        /* f[n] = 1 iff n = m^e, m >= 2, e >= 2 */
{
    u8 *f = calloc(lim + 1, 1);
    for (u64 m = 2; m * m <= lim; m++) for (u64 v = m * m; v <= lim; v *= m) { f[v] = 1; if (v > lim / m) break; }
    return f;
}
static void gen_perfpow(u64 *t, long cnt)      /* A001597: 1 and m^e, e >= 2 */
{
    /* about sqrt(x) of them below x: collect the powers instead of flagging */
    for (u64 lim = 1ULL << 34; ; lim *= 4) {
        size_t cap = 2 * (size_t)sqrtl((long double)lim) + 1000, m = 0;
        u64 *v = malloc(cap * sizeof(u64));
        v[m++] = 1;
        for (u64 b = 2; b * b <= lim; b++) for (u64 x = b * b; ; x *= b) { v[m++] = x; if (x > lim / b) break; }
        qsort(v, m, sizeof(u64), cmp_u64);
        long k = 0;
        for (size_t i = 0; i < m && k < cnt; i++) if (i == 0 || v[i] != v[i - 1]) t[k++] = v[i];
        free(v); if (k == cnt) return;
    }
}
static void gen_nonperfpow(u64 *t, long cnt)   /* A007916 */
{
    for (u64 lim = 1ULL << 18; ; lim *= 2) {
        u8 *f = perfect_power_flags(lim); long k = 0;
        for (u64 n = 2; n <= lim && k < cnt; n++) if (!f[n]) t[k++] = n;
        free(f); if (k == cnt) return;
    }
}

static void gen_figfig(u64 *t, long cnt)       /* A005228: a(n+1) = a(n) + (n-th number not in the sequence) */
{
    u64 lim = 64; while (lim < (u64)cnt * (u64)cnt) lim *= 2;
    u8 *in = calloc(lim / 8 + 2, 1);
    #define IN(x) ((in[(x) >> 3] >> ((x) & 7)) & 1)
    t[0] = 1; in[0] |= 2;
    u64 c = 1;                                  /* last complement term used */
    for (long i = 1; i < cnt; i++) {
        do c++; while (IN(c));
        t[i] = t[i - 1] + c;
        if (t[i] >= lim) { fprintf(stderr, "figfig limit\n"); exit(1); }
        in[t[i] >> 3] |= (u8)(1 << (t[i] & 7));
    }
    #undef IN
    free(in);
}

/* arithmetic progressions */
N_FILTER(gen_non5,     1, n % 5 != 0)                                         /* A047201 */
static void gen_ap(u64 *t, long cnt, u64 a, u64 b) { for (long i = 0; i < cnt; i++) t[i] = a * (u64)i + b; }
static void gen_mult4(u64 *t, long cnt) { gen_ap(t, cnt, 4, 0); }   /* A008586 */
static void gen_mult5(u64 *t, long cnt) { gen_ap(t, cnt, 5, 0); }   /* A008587 */
static void gen_mult6(u64 *t, long cnt) { gen_ap(t, cnt, 6, 0); }   /* A008588 */
static void gen_mult7(u64 *t, long cnt) { gen_ap(t, cnt, 7, 0); }   /* A008589 */
static void gen_3np1(u64 *t, long cnt)  { gen_ap(t, cnt, 3, 1); }   /* A016777 */
static void gen_3np2(u64 *t, long cnt)  { gen_ap(t, cnt, 3, 2); }   /* A016789 */
static void gen_4np1(u64 *t, long cnt)  { gen_ap(t, cnt, 4, 1); }   /* A016813 */
static void gen_4np3(u64 *t, long cnt)  { gen_ap(t, cnt, 4, 3); }   /* A004767 */

/* ---------------------------- the two hundred added 26 Sep 2026 ------ */
/* Sequences 201-400.  Every one is checked against its OEIS terms.      */

#define POLY(name, n_first, expr) \
static void name(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { u64 n = (u64)i + (n_first); t[i] = (expr); } }

POLY(gen_gon9,   0, n * (7 * n - 5) / 2)                 /* A001106 */
POLY(gen_gon10,  0, n * (4 * n - 3))                     /* A001107 */
POLY(gen_gon11,  0, n * (9 * n - 7) / 2)                 /* A051682 */
POLY(gen_gon12,  0, n * (5 * n - 4))                     /* A051624 */
POLY(gen_gon13,  0, n * (11 * n - 9) / 2)                /* A051865 */
POLY(gen_gon14,  0, n * (6 * n - 5))                     /* A051866 */
POLY(gen_gon15,  0, n * (13 * n - 11) / 2)               /* A051867 */
POLY(gen_gon16,  0, n * (7 * n - 6))                     /* A051868 */
POLY(gen_cen9,   1, (9 * n * n - 9 * n + 2) / 2)         /* A060544 */
POLY(gen_cen10,  1, 5 * n * n - 5 * n + 1)               /* A062786 */
POLY(gen_cen11,  1, (11 * n * n - 11 * n + 2) / 2)       /* A069125 */
POLY(gen_cen12,  1, 6 * n * (n - 1) + 1)                 /* A003154 */
POLY(gen_cen7,   1, (7 * n * n - 7 * n + 2) / 2)         /* A069099 */
POLY(gen_p4n13,  0, (4 * n + 1) * (4 * n + 3))           /* A001539 */
POLY(gen_2n2np1, 0, 2 * n * (2 * n + 1))                 /* A002943 */
POLY(gen_n4nm1,  0, n * (4 * n - 1))                     /* A033991 */
POLY(gen_n4np1,  0, n * (4 * n + 1))                     /* A007742 */
POLY(gen_2n2nm1, 0, 2 * n * (2 * n - 1))                 /* A002939 */
POLY(gen_cpoly,  1, n * n - n + 1)                       /* A002061, from a(1) = 1 (a(0) = a(1) = 1) */
POLY(gen_2n2p1,  0, 2 * n * n + 1)                       /* A058331 */
POLY(gen_2n2m1,  1, 2 * n * n - 1)                       /* A056220, from a(1) = 1 (a(0) = -1) */
POLY(gen_n2np2,  0, n * n + n + 2)                       /* A014206 */
POLY(gen_nnp1sq, 0, n + (n + 1) * (n + 1))               /* A028387 */
POLY(gen_n2m2,   2, n * n - 2)                           /* A008865, from a(2) = 2 (a(1) = -1) */
POLY(gen_2pyr,   0, n * (n + 1) * (2 * n + 1) / 3)       /* A006331 */
POLY(gen_pyr7,   0, n * (n + 1) * (5 * n - 2) / 6 + 0 * (n == 0))  /* A002413 */
POLY(gen_pyr8,   1, n * (n + 1) * (2 * n - 1) / 2)       /* A002414 */
POLY(gen_pyr9,   0, n ? n * (n + 1) * (7 * n - 4) / 6 : 0)  /* A007584 */
POLY(gen_pyr10,  0, n ? n * (n + 1) * (8 * n - 5) / 6 : 0)  /* A007585 */
POLY(gen_centet, 0, (2 * n + 1) * (n * n + n + 3) / 3)   /* A005894 */
POLY(gen_dodeca, 0, n ? n * (3 * n - 1) * (3 * n - 2) / 2 : 0)  /* A006566 */
POLY(gen_icosa,  1, n * (5 * n * n - 5 * n + 2) / 2)     /* A006564 */
POLY(gen_rhdod,  1, n * n * n * n - (n - 1) * (n - 1) * (n - 1) * (n - 1))  /* A005917 */
POLY(gen_cubep1, 0, n * n * n + 1)                       /* A001093, from a(0) = 1 (a(-1) = 0) */
POLY(gen_cubem1, 1, n * n * n - 1)                       /* A068601 */
POLY(gen_cubepn, 0, n * n * n + n)                       /* A034262 */
POLY(gen_hexprism, 0, (n + 1) * (3 * n * n + 3 * n + 1)) /* A005915 */
POLY(gen_cenicosa, 0, (2 * n + 1) * (5 * n * n + 5 * n + 3) / 3)  /* A005902 */
POLY(gen_trunctet, 0, (n + 1) * (23 * n * n + 19 * n + 6) / 6)    /* A005906 */
POLY(gen_secpentodd, 0, n ? (2 * n - 1) * (3 * n - 1) : 1)        /* A033568 */

/* Beatty and nearest-integer sequences */
static const long double E_L  = 2.718281828459045235360287471352662498L;
static const long double PI_L = 3.141592653589793238462643383279502884L;
POLY(gen_beattye,   0, (u64)floorl((long double)n * E_L))                    /* A022843 */
POLY(gen_beattypi,  0, (u64)floorl((long double)n * PI_L))                   /* A022844 */
POLY(gen_beattyee1, 1, (u64)floorl((long double)n * E_L / (E_L - 1.0L)))     /* A054385 */
POLY(gen_beatty2s3, 1, 2 * n + isqrt_u64(3 * n * n))                         /* A003512 */
POLY(gen_halfs2,    0, isqrt_u64((2 * n + 1) * (2 * n + 1) / 2))             /* A001953 */
POLY(gen_halfs22,   0, 2 * n + 1 + isqrt_u64((2 * n + 1) * (2 * n + 1) / 2)) /* A001954 */
POLY(gen_beattys55, 1, (5 * n + isqrt_u64(5 * n * n)) / 2)                   /* A003231 */
POLY(gen_zeckeven,  1, (n + isqrt_u64(5 * n * n)) / 2 - 1)                   /* A022342 */
POLY(gen_wythAA,    1, (n + isqrt_u64(5 * n * n)) / 2 + n - 1)               /* A003622 */
POLY(gen_near2,     0, (1 + isqrt_u64(8 * n * n)) / 2)                       /* A022846 */
POLY(gen_near3,     0, (1 + isqrt_u64(12 * n * n)) / 2)                      /* A022847 */
POLY(gen_near5,     0, (1 + isqrt_u64(20 * n * n)) / 2)                      /* A022848 */
POLY(gen_beattys5m1, 1, isqrt_u64(5 * n * n) - n)                            /* A001961 */

/* primes */
static void gen_prime_filter_m(u64 *t, long cnt, u64 lim, int (*pred)(u64, const u8 *), u64 mult)
{
    for (; ; lim *= 2) {
        u8 *c = composite_flags(mult * lim + 8); long k = 0;
        for (u64 p = 2; p <= lim && k < cnt; p++) if (!c[p] && pred(p, c)) t[k++] = p;
        free(c); if (k == cnt) return;
    }
}
static int p_1mod12(u64 p, const u8 *c)  { (void)c; return p % 12 == 1; }
static int p_5mod12(u64 p, const u8 *c)  { (void)c; return p % 12 == 5; }
static int p_7mod12(u64 p, const u8 *c)  { (void)c; return p % 12 == 7; }
static int p_11mod12(u64 p, const u8 *c) { (void)c; return p % 12 == 11; }
static int p_14mod5(u64 p, const u8 *c)  { (void)c; return p % 5 == 1 || p % 5 == 4; }
static int p_23mod5(u64 p, const u8 *c)  { (void)c; return p % 5 == 2 || p % 5 == 3; }
static int p_127mod8(u64 p, const u8 *c) { (void)c; return p % 8 == 1 || p % 8 == 2 || p % 8 == 7; }
static int p_pm3mod8(u64 p, const u8 *c) { (void)c; return p % 8 == 3 || p % 8 == 5; }
static int p_pm1mod8(u64 p, const u8 *c) { (void)c; return p % 8 == 1 || p % 8 == 7; }
static int p_13mod8(u64 p, const u8 *c)  { (void)c; return p % 8 == 1 || p % 8 == 3; }
static int p_first1(u64 p, const u8 *c)  { (void)c; while (p >= 10) p /= 10; return p == 1; }
static int p_2pp3(u64 p, const u8 *c)    { return pr(2 * p + 3, c); }
static int p_3pp2(u64 p, const u8 *c)    { return pr(3 * p + 2, c); }
static int p_4pp1(u64 p, const u8 *c)    { return pr(4 * p + 1, c); }
static int p_chen(u64 p, const u8 *c)    { if (pr(p + 2, c)) return 1; need_spf(p + 2, "chen"); return omega_big(p + 2) == 2; }
static int p_triple(u64 p, const u8 *c)  { return pr(p + 6, c) && (pr(p + 2, c) || pr(p + 4, c)); }

static int primroot(u64 g, u64 p)              /* g is a primitive root mod p (p odd prime, p <= SPF_MAX) */
{
    if (g % p == 0) return 0;
    u64 m = p - 1;
    while (m > 1) { u64 q = spf[m]; while (m % q == 0) m /= q; if (powmod(g, (p - 1) / q, p) == 1) return 0; }
    return 1;
}
static int p_root2(u64 p, const u8 *c)  { (void)c; return p > 2 && primroot(2, p); }
static int p_root10(u64 p, const u8 *c) { (void)c; return p > 5 && primroot(10, p); }

static void gen_p1mod12(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 23, p_1mod12); }   /* A068228 */
static void gen_p5mod12(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 23, p_5mod12); }   /* A040117 */
static void gen_p7mod12(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 23, p_7mod12); }   /* A068229 */
static void gen_p11mod12(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 23, p_11mod12); }  /* A068231 */
static void gen_p14mod5(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 22, p_14mod5); }   /* A045468 */
static void gen_p23mod5(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 22, p_23mod5); }   /* A003631 */
static void gen_p127mod8(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 22, p_127mod8); }  /* A038873 */
static void gen_ppm3mod8(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 22, p_pm3mod8); }  /* A003629 */
static void gen_ppm1mod8(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 22, p_pm1mod8); }  /* A001132 */
static void gen_p13mod8(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 22, p_13mod8); }   /* A033200 */
static void gen_pfirst1(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 22, p_first1); }   /* A045707 */
static void gen_p2pp3(u64 *t, long cnt)    { gen_prime_filter(t, cnt, 1ULL << 23, p_2pp3); }     /* A023204 */
static void gen_p3pp2(u64 *t, long cnt)    { gen_prime_filter_m(t, cnt, 1ULL << 23, p_3pp2, 3); } /* A023208 */
static void gen_p4pp1(u64 *t, long cnt)    { gen_prime_filter_m(t, cnt, 1ULL << 23, p_4pp1, 4); } /* A023212 */
static void gen_chen(u64 *t, long cnt)     { gen_prime_filter(t, cnt, 1ULL << 21, p_chen); }     /* A109611 */
static void gen_triples(u64 *t, long cnt)  { gen_prime_filter(t, cnt, 1ULL << 27, p_triple); }   /* A007529 */
static void gen_root2(u64 *t, long cnt)    { gen_prime_filter(t, cnt, 1ULL << 21, p_root2); }    /* A001122 */
static void gen_root10(u64 *t, long cnt)   { gen_prime_filter(t, cnt, 1ULL << 21, p_root10); }   /* A001913 */

static int k_up6(u64 a, u64 b, u64 c)   { (void)c; return b - a == 6; }
static int k_lo8(u64 a, u64 b, u64 c)   { (void)a; return c - b == 8; }
static int k_lo10(u64 a, u64 b, u64 c)  { (void)a; return c - b == 10; }
static int k_lo12(u64 a, u64 b, u64 c)  { (void)a; return c - b == 12; }
static void gen_up6(u64 *t, long cnt)  { gen_prime_nbr(t, cnt, k_up6); }    /* A031925 */
static void gen_lo8(u64 *t, long cnt)  { gen_prime_nbr(t, cnt, k_lo8); }    /* A031926 */
static void gen_lo10(u64 *t, long cnt) { gen_prime_nbr(t, cnt, k_lo10); }   /* A031928 */
static void gen_lo12(u64 *t, long cnt) { gen_prime_nbr(t, cnt, k_lo12); }   /* A031930 */

static void gen_ppp(u64 *t, long cnt)          /* A038580: prime(prime(prime(n))) */
{
    long np = 22000000;
    u64 *p = malloc((size_t)np * sizeof(u64));
    gen_primes(p, np);
    for (long i = 0; i < cnt; i++) {
        u64 j = p[p[i] - 1] - 1;
        if ((long)j >= np) { fprintf(stderr, "ppp: table too short\n"); exit(1); }
        t[i] = p[j];
    }
    free(p);
}
static void gen_pnonprimeidx(u64 *t, long cnt) /* A007821: prime(k), k = 1 or composite */
{
    long np = 3 * cnt + 1000;
    u64 *p = malloc((size_t)np * sizeof(u64));
    gen_primes(p, np);
    long k = 0;
    for (long i = 0; i < np && k < cnt; i++) { u64 idx = (u64)i + 1; need_spf(idx, "pnonprimeidx"); if (idx == 1 || spf[idx] != idx) t[k++] = p[i]; }
    free(p);
    if (k < cnt) { fprintf(stderr, "pnonprimeidx short\n"); exit(1); }
}

#define PRIME_MAP(name, extra, expr) \
static void name(u64 *t, long cnt) \
{ u64 *p = malloc((size_t)(cnt + (extra)) * sizeof(u64)); gen_primes(p, cnt + (extra)); \
  for (long i = 0; i < cnt; i++) t[i] = (expr); \
  free(p); }

PRIME_MAP(gen_fourp,   0, 4 * p[i])                          /* A001749 */
PRIME_MAP(gen_ppm1,    0, p[i] * (p[i] - 1))                 /* A036689 */
PRIME_MAP(gen_ppp1,    0, p[i] * (p[i] + 1))                 /* A036690 */
PRIME_MAP(gen_psqm2,   0, p[i] * p[i] - 2)                   /* A049001 */
PRIME_MAP(gen_psum3,   2, p[i] + p[i + 1] + p[i + 2])        /* A034961 */
PRIME_MAP(gen_psum4,   3, p[i] + p[i + 1] + p[i + 2] + p[i + 3])  /* A034963 */
PRIME_MAP(gen_pavg,    2, (p[i + 1] + p[i + 2]) / 2)         /* A024675 */
PRIME_MAP(gen_pminus2, 0, p[i] - 2)                          /* A040976 */
PRIME_MAP(gen_phalfm,  1, (p[i + 1] - 1) / 2)                /* A005097 */
PRIME_MAP(gen_phalfp,  1, (p[i + 1] + 1) / 2)                /* A006254 */

static void gen_twin6m(u64 *t, long cnt)       /* A002822: 6m - 1, 6m + 1 both prime */
{
    for (u64 lim = 1ULL << 24; ; lim *= 2) {
        u8 *c = composite_flags(lim + 2); long k = 0;
        for (u64 m = 1; 6 * m + 1 <= lim && k < cnt; m++) if (!c[6 * m - 1] && !c[6 * m + 1]) t[k++] = m;
        free(c); if (k == cnt) return;
    }
}
static void gen_sum2p(u64 *t, long cnt)        /* A014091: sums of two primes (Goldbach holds far beyond this range) */
{
    for (u64 lim = 1ULL << 18; ; lim *= 2) {
        u8 *c = composite_flags(lim); long k = 0;
        for (u64 n = 4; n <= lim && k < cnt; n++) if (!(n & 1) || !c[n - 2]) t[k++] = n;
        free(c); if (k == cnt) return;
    }
}
static void gen_notsum2p(u64 *t, long cnt)     /* A014092 */
{
    for (u64 lim = 1ULL << 18; ; lim *= 2) {
        u8 *c = composite_flags(lim); long k = 0;
        for (u64 n = 1; n <= lim && k < cnt; n++) {
            int s = (n >= 4 && !(n & 1)) || (n >= 4 && (n & 1) && !c[n - 2]);
            if (!s) t[k++] = n;
        }
        free(c); if (k == cnt) return;
    }
}

/* positive-definite binary forms and small sums, by enumeration */
static void gen_form_pr(u64 *t, long cnt, u64 A, u64 B, u64 C, int primes_only, int from0)
{   /* values A x^2 + B x y + C y^2 with x, y >= 0 */
    for (u64 lim = 1ULL << 20; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 y = 0; C * y * y <= lim; y++)
            for (u64 x = 0; A * x * x + B * x * y + C * y * y <= lim; x++) h[A * x * x + B * x * y + C * y * y] = 1;
        u8 *c = primes_only ? composite_flags(lim) : NULL;
        for (u64 n = from0 ? 0 : 1; n <= lim && k < cnt; n++)
            if (h[n] && (!primes_only || (n >= 2 && !c[n]))) t[k++] = n;
        free(h); free(c); if (k == cnt) return;
    }
}
static void gen_pform712(u64 *t, long cnt) { gen_form_pr(t, cnt, 1, 1, 2, 1, 0); }   /* A106856 */
static void gen_pform15(u64 *t, long cnt)  { gen_form_pr(t, cnt, 1, 0, 5, 1, 0); }   /* A033205 */
static void gen_form16(u64 *t, long cnt)   { gen_form_pr(t, cnt, 1, 0, 6, 0, 1); }   /* A002481 */
static void gen_form23(u64 *t, long cnt)   { gen_form_pr(t, cnt, 2, 0, 3, 0, 1); }   /* A002480 */
static void gen_form15(u64 *t, long cnt)   { gen_form_pr(t, cnt, 1, 0, 5, 0, 1); }   /* A020669 */

static void gen_form112(u64 *t, long cnt)      /* A000401: x^2 + y^2 + 2 z^2 */
{
    for (u64 lim = 1ULL << 18; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 z = 0; 2 * z * z <= lim; z++)
            for (u64 y = 0; 2 * z * z + y * y <= lim; y++)
                for (u64 x = 0; x <= y && 2 * z * z + y * y + x * x <= lim; x++) h[2 * z * z + y * y + x * x] = 1;
        for (u64 n = 0; n <= lim && k < cnt; n++) if (h[n]) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static void gen_sq2sq(u64 *t, long cnt)        /* A028982: squares and twice squares */
{
    long k = 0; u64 a = 1, b = 1;
    while (k < cnt) { u64 s = a * a, d = 2 * b * b; if (s < d) { t[k++] = s; a++; } else { t[k++] = d; b++; } }
}
static int is_sq2sq(u64 n) { if (is_square(n)) return 1; return !(n & 1) && is_square(n / 2); }
N_FILTER(gen_sigmaeven, 1, !is_sq2sq(n))                                      /* A028983 */
N_FILTER(gen_exactly3sq, 1, !is_4a8b7(n) && !is_square(n) && !twosq_small(n)) /* A000419 */
static void gen_twocubes(u64 *t, long cnt)     /* A003325: x^3 + y^3, x, y >= 1 */
{
    for (u64 lim = 1ULL << 22; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 x = 1; 2 * x * x * x <= lim; x++)
            for (u64 y = x; x * x * x + y * y * y <= lim; y++) h[x * x * x + y * y * y] = 1;
        for (u64 n = 1; n <= lim && k < cnt; n++) if (h[n]) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static void gen_twotri_any(u64 *t, long cnt, int want)
{
    for (u64 lim = 1ULL << 18; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 i = 0; i * (i + 1) <= lim; i++)
            for (u64 j = i; i * (i + 1) / 2 + j * (j + 1) / 2 <= lim; j++) h[i * (i + 1) / 2 + j * (j + 1) / 2] = 1;
        for (u64 n = 0; n <= lim && k < cnt; n++) if (h[n] == want) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static void gen_twotri(u64 *t, long cnt)    { gen_twotri_any(t, cnt, 1); }   /* A020756 */
static void gen_nottwotri(u64 *t, long cnt) { gen_twotri_any(t, cnt, 0); }   /* A020757 */

/* multiplicative */
static int bitexp_even(u64 n)                  /* A000379: the 1 bits of all exponents, even in total */
{ int b = 0; while (n > 1) { u64 p = spf[n]; int e = 0; while (n % p == 0) { n /= p; e++; } b += popcount_u64((u64)e); } return !(b & 1); }
static int biqfree(u64 n) { while (n > 1) { u64 p = spf[n]; int e = 0; while (n % p == 0) { n /= p; e++; } if (e >= 4) return 0; } return 1; }
static u64 phi_small(u64 n) { u64 m = n, ph = n; while (m > 1) { u64 p = spf[m]; ph = ph / p * (p - 1); while (m % p == 0) m /= p; } return ph; }
static int blum(u64 n)
{
    if (omega_big(n) != 2 || !squarefree_small(n)) return 0;
    u64 p = spf[n], q = n / p; return p % 4 == 3 && q % 4 == 3;
}
static int ndigits10(u64 n) { int c = 0; do { c++; n /= 10; } while (n); return c; }
static int factor_digits(u64 n)                /* digits of the prime factorization, exponents > 1 counted */
{ int c = 0; while (n > 1) { u64 p = spf[n]; int e = 0; while (n % p == 0) { n /= p; e++; } c += ndigits10(p) + (e > 1 ? ndigits10((u64)e) : 0); } return c; }
static int smith(u64 n)
{ if (spf[n] == n) return 0; int s = 0; u64 m = n; while (m > 1) { u64 p = spf[m]; m /= p; s += digitsum(p); } return s == digitsum(n); }
static int hoax(u64 n)
{ if (spf[n] == n) return 0; int s = 0; u64 m = n; while (m > 1) { u64 p = spf[m]; s += digitsum(p); while (m % p == 0) m /= p; } return s == digitsum(n); }
static int hati(u64 n) { int i = 0; while (n % 2 == 0) { n /= 2; i++; } while (n % 3 == 0) { n /= 3; i++; } return !(i & 1); }

SPF_FILTER(gen_almost6,  64, omega_big(n) == 6)                              /* A046306 */
SPF_FILTER(gen_almost7,  128, omega_big(n) == 7)                             /* A046308 */
SPF_FILTER(gen_almost8,  256, omega_big(n) == 8)                             /* A046310 */
SPF_FILTER(gen_oddsemi,  9, (n & 1) && omega_big(n) == 2)                    /* A046315 */
SPF_FILTER(gen_oddpq,    15, (n & 1) && omega_big(n) == 2 && squarefree_small(n))  /* A046388 */
SPF_FILTER(gen_atleast3, 30, omega_distinct(n) >= 3)                         /* A000977 */
SPF_FILTER(gen_bitexp,   1, bitexp_even(n))                                  /* A000379 */
SPF_FILTER(gen_biqfree,  1, biqfree(n))                                      /* A046100 */
SPF_FILTER(gen_biqful,   16, !biqfree(n))                                    /* A046101 */
SPF_FILTER(gen_cyclic,   1, gcd_u64(n, phi_small(n)) == 1)                   /* A003277 */
SPF_FILTER(gen_rough11,  1, n == 1 || spf[n] >= 11)                          /* A008364 */
SPF_FILTER(gen_rough13,  1, n == 1 || spf[n] >= 13)                          /* A008365 */
SPF_FILTER(gen_rough17,  1, n == 1 || spf[n] >= 17)                          /* A008366 */
SPF_FILTER(gen_blum,     21, blum(n))                                        /* A016105 */
SPF_FILTER(gen_tau4,     6, tau_small(n) == 4)                               /* A030513 */
SPF_FILTER(gen_tau6,     12, tau_small(n) == 6)                              /* A030515 */
SPF_FILTER(gen_tau8,     24, tau_small(n) == 8)                              /* A030626 */
static void gen_tau10(u64 *t, long cnt)        /* A030628: 1, p q^4 (p != q) and p^9, i.e. 10 divisors */
{
    for (u64 lim = 1ULL << 26; ; lim *= 2) {
        u8 *c = composite_flags(lim / 16 + 2);
        size_t cap = 4000000, m = 0; u64 *v = malloc(cap * sizeof(u64));
        v[m++] = 1;
        for (u64 q = 2; q * q * q * q * 2 <= lim; q++) {
            if (c[q]) continue;
            u64 q4 = q * q * q * q;
            for (u64 pp = 2; pp * q4 <= lim; pp++) if (!c[pp] && pp != q) { if (m == cap) { fprintf(stderr, "tau10 cap\n"); exit(1); } v[m++] = pp * q4; }
            if (q4 * q4 * q <= lim) v[m++] = q4 * q4 * q;
        }
        free(c);
        if ((long)m >= cnt) { qsort(v, m, sizeof(u64), cmp_u64); memcpy(t, v, (size_t)cnt * sizeof(u64)); free(v); return; }
        free(v);
    }
}
SPF_FILTER(gen_tau12,    60, tau_small(n) == 12)                             /* A030630 */
SPF_FILTER(gen_equidig,  1, n == 1 || factor_digits(n) == ndigits10(n))      /* A046758 */
SPF_FILTER(gen_wasteful, 4, n > 1 && factor_digits(n) > ndigits10(n))        /* A046760 */
SPF_FILTER(gen_smith,    4, smith(n))                                        /* A006753 */
SPF_FILTER(gen_hoax,     22, hoax(n))                                        /* A019506 */
SPF_FILTER(gen_hati,     1, hati(n))                                         /* A036668 */
SPF_FILTER(gen_sqfcomp,  6, squarefree_small(n) && spf[n] != n)              /* A120944 */

static void gen_sqfsq(u64 *t, long cnt)        /* A062503: squares of squarefree numbers */
{ long k = 0; for (u64 m = 1; k < cnt; m++) { need_spf(m, "sqfsq"); if (squarefree_small(m)) t[k++] = m * m; } }

static void gen_evennontot(u64 *t, long cnt)   /* A005277 */
{
    for (u64 N = 1ULL << 18; ; N *= 2) {
        u8 *f = totient_flags(N); long k = 0;
        for (u64 n = 2; n <= N && k < cnt; n += 2) if (!f[n]) t[k++] = n;
        free(f); if (k == cnt) return;
    }
}

static void gen_brilliant(u64 *t, long cnt)    /* A078972: pq, p <= q primes with the same number of digits */
{
    u64 *p = malloc(1300 * sizeof(u64)); gen_primes(p, 1229);   /* the primes below 10^4 */
    size_t cap = 700000, m = 0; u64 *v = malloc(cap * sizeof(u64));
    for (int i = 0; i < 1229; i++)
        for (int j = i; j < 1229 && ndigits10(p[j]) == ndigits10(p[i]); j++) v[m++] = p[i] * p[j];
    qsort(v, m, sizeof(u64), cmp_u64);
    if ((long)m < cnt) { fprintf(stderr, "brilliant short\n"); exit(1); }
    memcpy(t, v, (size_t)cnt * sizeof(u64));
    free(p); free(v);
}

/* smooth numbers: the smallest cnt of them, sorted (the bound doubles until there are enough) */
static size_t smooth_m, smooth_cap; static u64 *smooth_v;
static void smooth_rec(const u64 *ps, int np, int i, u64 v, u64 lim)
{
    if (i == np) {
        if (smooth_m == smooth_cap) { smooth_cap *= 2; smooth_v = realloc(smooth_v, smooth_cap * sizeof(u64)); }
        smooth_v[smooth_m++] = v; return;
    }
    for (u64 x = v; ; x *= ps[i]) { smooth_rec(ps, np, i + 1, x, lim); if (x > lim / ps[i]) break; }
}
static long smooth_all(u64 *t, long cnt, const u64 *ps, int np)
{
    for (u64 lim = 1ULL << 20; ; lim *= 2) {
        smooth_cap = 1 << 16; smooth_m = 0; smooth_v = malloc(smooth_cap * sizeof(u64));
        smooth_rec(ps, np, 0, 1, lim);
        if ((long)smooth_m >= cnt || lim >= (1ULL << 53)) {
            qsort(smooth_v, smooth_m, sizeof(u64), cmp_u64);
            long k = (long)smooth_m < cnt ? (long)smooth_m : cnt;
            memcpy(t, smooth_v, (size_t)k * sizeof(u64));
            free(smooth_v);
            return k;
        }
        free(smooth_v);
    }
}
static const u64 PS17[7] = {2, 3, 5, 7, 11, 13, 17};
static void gen_smooth11(u64 *t, long cnt) { if (smooth_all(t, cnt, PS17, 5) < cnt) { fprintf(stderr, "smooth11 short\n"); exit(1); } }  /* A051038 */
static void gen_smooth13(u64 *t, long cnt) { if (smooth_all(t, cnt, PS17, 6) < cnt) { fprintf(stderr, "smooth13 short\n"); exit(1); } }  /* A080197 */
static void gen_smooth17(u64 *t, long cnt) { if (smooth_all(t, cnt, PS17, 7) < cnt) { fprintf(stderr, "smooth17 short\n"); exit(1); } }  /* A080681 */

/* digits */
N_FILTER(gen_has2, 1, has_digit(n, 2))                                        /* A011532 */
N_FILTER(gen_has3, 1, has_digit(n, 3))                                        /* A011533 */
N_FILTER(gen_has4, 1, has_digit(n, 4))                                        /* A011534 */
N_FILTER(gen_has5, 1, has_digit(n, 5))                                        /* A011535 */
N_FILTER(gen_has6, 1, has_digit(n, 6))                                        /* A011536 */
N_FILTER(gen_has7, 1, has_digit(n, 7))                                        /* A011537 */
N_FILTER(gen_has8, 1, has_digit(n, 8))                                        /* A011538 */
N_FILTER(gen_no3,  0, !has_digit(n, 3))                                       /* A052405 */
N_FILTER(gen_no4,  0, !has_digit(n, 4))                                       /* A052406 */
N_FILTER(gen_no5,  0, !has_digit(n, 5))                                       /* A052413 */
N_FILTER(gen_no6,  0, !has_digit(n, 6))                                       /* A052414 */
N_FILTER(gen_no7,  0, !has_digit(n, 7))                                       /* A052419 */
N_FILTER(gen_no8,  0, !has_digit(n, 8))                                       /* A052421 */
static void gen_base7(u64 *t, long cnt) { gen_base_read(t, cnt, 7); }   /* A007093 */
static void gen_base8(u64 *t, long cnt) { gen_base_read(t, cnt, 8); }   /* A007094 */
static void gen_base9(u64 *t, long cnt) { gen_base_read(t, cnt, 9); }   /* A007095 */
static int is_fibval(u64 n) { u64 a = 0, b = 1; while (a < n) { u64 c = a + b; a = b; b = c; } return a == n; }
N_FILTER(gen_dsodd,  1, digitsum(n) & 1)                                      /* A054684 */
N_FILTER(gen_dsfib,  0, is_fibval((u64)digitsum(n)))                          /* A028840 */
static int has_zero_base(u64 n, u64 b) { if (!n) return 1; while (n) { if (n % b == 0) return 1; n /= b; } return 0; }
N_FILTER(gen_has0b3, 0, has_zero_base(n, 3))                                  /* A081605 */
static int no_prime_digit(u64 n) { do { int d = (int)(n % 10); if (d == 2 || d == 3 || d == 5 || d == 7) return 0; n /= 10; } while (n); return 1; }
N_FILTER(gen_noprimedig, 0, no_prime_digit(n))                                /* A084984 */
static int alt_parity(u64 n) { int last = -1; do { int d = (int)(n % 10) & 1; if (d == last) return 0; last = d; n /= 10; } while (n); return 1; }
N_FILTER(gen_altpar, 0, alt_parity(n))                                        /* A030141 */
static void gen_no3b4(u64 *t, long cnt)        /* A023717, from a(0) = 0: n in base 3, read in base 4 */
{
    for (long i = 0; i < cnt; i++) { u64 n = (u64)i, v = 0, p = 1; while (n) { v += (n % 3) * p; p *= 4; n /= 3; } t[i] = v; }
}

/* palindromes in base b, in order: by length, then by the first half */
static void gen_palb(u64 *t, long cnt, u64 b)
{
    long k = 0; t[k++] = 0;
    for (int len = 1; k < cnt; len++) {
        int h = (len + 1) / 2;
        u64 lo = 1; for (int i = 1; i < h; i++) lo *= b;
        for (u64 half = lo; half < lo * b && k < cnt; half++) {
            u64 v = half, r = (len & 1) ? half / b : half;
            while (r) { v = v * b + r % b; r /= b; }
            t[k++] = v;
        }
    }
}
static void gen_pal3(u64 *t, long cnt) { gen_palb(t, cnt, 3); }   /* A014190 */
static void gen_pal4(u64 *t, long cnt) { gen_palb(t, cnt, 4); }   /* A014192 */
static void gen_pal5(u64 *t, long cnt) { gen_palb(t, cnt, 5); }   /* A029952 */
static void gen_pal6(u64 *t, long cnt) { gen_palb(t, cnt, 6); }   /* A029953 */

/* binary patterns */
static int bitlen(u64 n) { int c = 0; while (n) { c++; n >>= 1; } return c; }
static int no00(u64 n) { if (!n) return 1; u64 m = (bitlen(n) >= 64) ? ~0ULL : ((1ULL << bitlen(n)) - 1); u64 z = ~n & m; return !(z & (z >> 1)); }
N_FILTER(gen_has11,  1, (n & (n >> 1)) != 0)                                  /* A004780 */
N_FILTER(gen_no00,   0, no00(n))                                              /* A003754 */
N_FILTER(gen_rsneg,  1, popcount_u64(n & (n >> 1)) & 1)                       /* A022155 */
N_FILTER(gen_rspos,  0, !(popcount_u64(n & (n >> 1)) & 1))                    /* A203463 */
static u64 oddpart(u64 n) { while (!(n & 1)) n >>= 1; return n; }
N_FILTER(gen_pf1,    1, oddpart(n) % 4 == 1)                                  /* A091072 */
N_FILTER(gen_pf3,    1, oddpart(n) % 4 == 3)                                  /* A091067 */
static int zeck_terms(u64 n)
{
    static u64 F[92]; static int nf = 0;
    if (!nf) { F[0] = 1; F[1] = 2; nf = 2; while (F[nf - 1] < (1ULL << 62)) { F[nf] = F[nf - 1] + F[nf - 2]; nf++; } }
    int c = 0; for (int i = nf - 1; i >= 0 && n; i--) if (F[i] <= n) { n -= F[i]; c++; }
    return c;
}
N_FILTER(gen_zeckodd, 1, zeck_terms(n) & 1)                                   /* A020899 */
static void gen_binself(u64 *t, long cnt)      /* A010061: not m + (binary weight of m) */
{
    for (u64 lim = 1ULL << 20; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 m = 0; m <= lim; m++) { u64 g = m + (u64)popcount_u64(m); if (g <= lim) h[g] = 1; }
        for (u64 n = 1; n <= lim && k < cnt; n++) if (!h[n]) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}

/* sieves and self-reference */
static void gen_ludic(u64 *t, long cnt)        /* A003309: 1, then the ludic sieve on 2, 3, 4, ... */
{
    /* Each step takes the first survivor x and removes the survivors at ranks 1, 1 + x,
       1 + 2x, ...  Ranks count from the front, so sieving the finite list 2..N gives
       exact terms: the cut-off only shortens the tail.  A Fenwick tree finds ranks.   */
    for (long N = 1L << 21; ; N *= 2) {
        int *tree = calloc((size_t)N + 1, sizeof(int));  /* position p holds the number p + 1 */
        for (long i = 1; i <= N; i++) { tree[i]++; long j = i + (i & -i); if (j <= N) tree[j] += tree[i]; }
        long LOG = 1; while ((LOG << 1) <= N) LOG <<= 1;
        long alive = N, k = 0;
        t[k++] = 1;
        while (alive > 0 && k < cnt) {
            u64 x = 0;
            long nrm = 0;
            for (long q = -1; q < nrm; q++) {
                long rank = q < 0 ? 1 : (nrm - q) * (long)x;          /* the head, then j x for j = nrm .. 1 */
                long pos = 0, rem = rank;
                for (long s2 = LOG; s2; s2 >>= 1) if (pos + s2 <= N && tree[pos + s2] < rem) { pos += s2; rem -= tree[pos]; }
                if (q < 0) {                                   /* the head: read it, remove it below */
                    x = (u64)pos + 2;
                    t[k++] = x;
                    nrm = (alive - 1) / (long)x;               /* ranks 1 + x, 1 + 2x, ... after the head */
                    for (long i = pos + 1; i <= N; i += i & -i) tree[i]--;
                    alive--;
                    /* with the head gone, rank 1 + j x becomes j x */
                    continue;
                }
                for (long i = pos + 1; i <= N; i += i & -i) tree[i]--;
                alive--;
            }
        }
        free(tree);
        if (k == cnt) return;
    }
}
static void gen_evenlucky(u64 *t, long cnt)    /* A045954: the lucky sieve on the even numbers */
{
    /* as for A000959: start from 2, 4, 6, ...; the step with m = t[j] (j >= 1: 4, 6, 10, ...)
       deletes every m-th survivor, so a survivor at q sat at q + floor((q-1)/(m-1)) before it */
    long jmax = 0;
    for (long i = 0; i < cnt; i++) {
        u64 n = (u64)(i + 1);
        while (jmax + 1 < i && t[jmax + 1] <= n) jmax++;
        u64 q = n;
        for (long j = jmax; j >= 1; j--) q += (q - 1) / (t[j] - 1);
        t[i] = 2 * q;
    }
}
static void gen_figcomp(u64 *t, long cnt)      /* A030124: the complement of A005228 */
{
    u64 *f = malloc((size_t)(cnt + 10) * sizeof(u64));
    gen_figfig(f, cnt + 10);
    long k = 0, j = 0;
    for (u64 n = 1; k < cnt; n++) { while (f[j] < n) j++; if (f[j] != n) t[k++] = n; }
    free(f);
}
static void gen_ulam_ab(u64 *t, long cnt, u64 u1, u64 u2)
{
    u64 lim = (u64)(16.0 * cnt + 1000);
    for (;;) {
        u8 *rep = calloc(lim + 1, 1);                /* 0, 1 or 2 (= more) representations */
        long k = 0; int over = 0;
        t[k++] = u1; t[k++] = u2; if (u1 + u2 <= lim) rep[u1 + u2] = 1;
        u64 last = u2;
        while (k < cnt) {
            u64 u = last + 1;
            while (u <= lim && rep[u] != 1) u++;
            if (u > lim) { over = 1; break; }
            for (long j = 0; j < k; j++) { u64 s = t[j] + u; if (s > lim) break; if (rep[s] < 2) rep[s]++; }
            t[k++] = u; last = u;
        }
        free(rep);
        if (!over) return;
        lim *= 2;
    }
}
static void gen_ulam13(u64 *t, long cnt) { gen_ulam_ab(t, cnt, 1, 3); }   /* A002859 */
static void gen_klarner(u64 *t, long cnt)      /* A002977: 1, and 2m + 1, 3m + 1 for every m in it */
{
    for (u64 lim = 1ULL << 20; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        h[1] = 1;
        for (u64 m = 1; m <= lim; m++) if (h[m]) { if (2 * m + 1 <= lim) h[2 * m + 1] = 1; if (3 * m + 1 <= lim) h[3 * m + 1] = 1; }
        for (u64 n = 1; n <= lim && k < cnt; n++) if (h[n]) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static void gen_addomega(u64 *t, long cnt)     /* A094222: 1, 2, then a(n+1) = a(n) + omega(a(n)) */
{ t[0] = 1; t[1] = 2; for (long i = 2; i < cnt; i++) { need_spf(t[i - 1], "addomega"); t[i] = t[i - 1] + (u64)omega_distinct(t[i - 1]); } }
static void gen_barriers(u64 *t, long cnt)     /* A005236: m + omega(m) <= n for every m < n */
{
    for (u64 lim = 1ULL << 24; ; lim *= 2) {
        u8 *om = calloc(lim + 1, 1);            /* omega by sieve */
        for (u64 p = 2; p <= lim; p++) if (!om[p]) for (u64 j = p; j <= lim; j += p) om[j]++;
        long k = 0; u64 reach = 0;              /* max of m + omega(m) over m < n */
        for (u64 n = 1; n <= lim && k < cnt; n++) {
            if (n >= 2 && reach <= n) t[k++] = n;
            if (n + om[n] > reach) reach = n + om[n];
        }
        free(om); if (k == cnt) return;
    }
}

/* partial sums */
static void gen_bitsum(u64 *t, long cnt)       /* A000788, from a(0) = 0 */
{ u64 s = 0; for (long i = 0; i < cnt; i++) { s += (u64)popcount_u64((u64)i); t[i] = s; } }
static void gen_digsum(u64 *t, long cnt)       /* A037123, from a(0) = 0 */
{ u64 s = 0; for (long i = 0; i < cnt; i++) { s += (u64)digitsum((u64)i); t[i] = s; } }
static void gen_omegasum(u64 *t, long cnt)     /* A013939, from a(1) = 0 */
{ u64 s = 0; for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); need_spf(n, "omegasum"); s += (u64)(n > 1 ? omega_distinct(n) : 0); t[i] = s; } }
static void gen_bigomegasum(u64 *t, long cnt)  /* A022559 from a(1) = 0 (a(0) = a(1) = 0) */
{ u64 s = 0; for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); need_spf(n, "bigomegasum"); s += (u64)(n > 1 ? omega_big(n) : 0); t[i] = s; } }
static void gen_halfsum(u64 *t, long cnt)      /* A005187, from a(0) = 0: a(n) = a(floor(n/2)) + n */
{ t[0] = 0; for (long i = 1; i < cnt; i++) t[i] = t[i / 2] + (u64)i; }
static void gen_pascalodd(u64 *t, long cnt)    /* A006046, from a(0) = 0: sum of 2^(binary weight of k), k < n */
{ u64 s = 0; for (long i = 0; i < cnt; i++) { t[i] = s; s += 1ULL << popcount_u64((u64)i); } }
static void gen_unitsum(u64 *t, long cnt)      /* A064608: sum of 2^omega(k), k <= n */
{ u64 s = 0; for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); need_spf(n, "unitsum"); s += 1ULL << (n > 1 ? omega_distinct(n) : 0); t[i] = s; } }

/* residue classes: n >= from with (mask >> (n % m)) & 1 */
static void gen_res(u64 *t, long cnt, u64 m, u64 mask, u64 from)
{ long k = 0; for (u64 n = from; k < cnt; n++) if ((mask >> (n % m)) & 1) t[k++] = n; }
#define RES(name, m, mask, from) static void name(u64 *t, long cnt) { gen_res(t, cnt, m, mask, from); }
RES(gen_r5_0234,  5, 0x1D, 0)     /* A047203 */
RES(gen_r5_14,    5, 0x12, 1)     /* A047209 */
RES(gen_r5_013,   5, 0x0B, 0)     /* A047220 */
RES(gen_r6_0234,  6, 0x1D, 0)     /* A047229 */
RES(gen_r6_02,    6, 0x05, 0)     /* A047238 */
RES(gen_r6_0123,  6, 0x0F, 0)     /* A047246 */
RES(gen_r6_1235,  6, 0x2E, 1)     /* A047255 */
RES(gen_r6_245,   6, 0x34, 2)     /* A047261 */
RES(gen_r6_015,   6, 0x23, 0)     /* A047266 */
RES(gen_r6_0135,  6, 0x2B, 0)     /* A047273 */
RES(gen_r10_1379, 10, 0x28A, 1)   /* A045572 */
RES(gen_r4_03,    4, 0x09, 0)     /* A014601 */
RES(gen_r4_12,    4, 0x06, 1)     /* A042963 */
RES(gen_r5_24,    5, 0x14, 2)     /* A047211 */
RES(gen_r3_02,    3, 0x05, 0)     /* A007494 */
RES(gen_r3_01,    3, 0x03, 0)     /* A032766 */
RES(gen_r5_024,   5, 0x15, 0)     /* A047212 */
N_FILTER(gen_cop21, 1, n % 3 && n % 7)                                        /* A160545 */

SPF_FILTER(gen_sqfpair, 1, squarefree_small(n) && squarefree_small(n + 1))   /* A007674 */

/* ---------------------------- the four hundred added 26 Sep 2026 ----- */
/* Sequences 401-800, found by family in the OEIS data and checked against
   their OEIS terms: generic generators with fixed parameters, one wrapper each. */

typedef __int128 i128;
static int isprime_s(i128 v) { return v >= 2 && is_prime((u64)v); }

static void gen_pres(u64 *t, long cnt, u64 M, u64 mask)        /* primes p with p mod M in a set */
{
    for (u64 lim = 1ULL << 22; ; lim *= 2) {
        u8 *c = composite_flags(lim); long k = 0;
        for (u64 p = 2; p <= lim && k < cnt; p++) if (!c[p] && ((mask >> (p % M)) & 1)) t[k++] = p;
        free(c); if (k == cnt) return;
    }
}
static void gen_plin(u64 *t, long cnt, long long A, long long B) /* primes p with A p + B prime */
{
    for (u64 lim = 1ULL << 22; ; lim *= 2) {
        u8 *c = composite_flags(lim); long k = 0;
        for (u64 p = 2; p <= lim && k < cnt; p++) if (!c[p] && isprime_s((i128)A * p + B)) t[k++] = p;
        free(c); if (k == cnt) return;
    }
}
static void gen_nlin(u64 *t, long cnt, long long A, long long B) /* n >= 0 with A n + B prime */
{ long k = 0; for (u64 n = 0; k < cnt; n++) if (isprime_s((i128)A * n + B)) t[k++] = n; }
static void gen_nsq(u64 *t, long cnt, long long C)              /* n >= 0 with n^2 + C prime */
{ long k = 0; for (u64 n = 0; k < cnt; n++) if (isprime_s((i128)n * n + C)) t[k++] = n; }
static void gen_bform(u64 *t, long cnt, u64 A, u64 B, u64 C, int primes)  /* A x^2 + B x y + C y^2, x, y >= 0 */
{
    for (u64 lim = 1ULL << 20; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 y = 0; C * y * y <= lim; y++)
            for (u64 x = 0; A * x * x + B * x * y + C * y * y <= lim; x++) h[A * x * x + B * x * y + C * y * y] = 1;
        for (u64 n = 0; n <= lim && k < cnt; n++) if (h[n] && (!primes || is_prime(n))) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static void gen_beattyl(u64 *t, long cnt, long double al, u64 n_first)  /* floor(n alpha) */
{ for (long i = 0; i < cnt; i++) t[i] = (u64)floorl((long double)(n_first + (u64)i) * al); }

/* the ones written out */
static int nodig_b(u64 n, u64 b, u64 d) { if (!n) return d != 0; while (n) { if (n % b == d) return 0; n /= b; } return 1; }
static int count_dig(u64 n, u64 b, u64 d) { int c = 0; do { if (n % b == d) c++; n /= b; } while (n); return c; }
static int div_by_a_digit(u64 n) { u64 m = n; while (m) { u64 d = m % 10; if (d && n % d == 0) return 1; m /= 10; } return 0; }
static int tern_sum(u64 n) { int s = 0; while (n) { s += (int)(n % 3); n /= 3; } return s; }
static int two_sq_primes(u64 n) { int c = 0; while (n > 1) { u64 p = spf[n]; int e = 0; while (n % p == 0) { n /= p; e++; } if (e >= 2) c++; } return c >= 2; }

N_FILTER(gen_a027697, 2, is_prime(n) && (popcount_u64(n) & 1))                        /* odious primes */
N_FILTER(gen_a027699, 2, is_prime(n) && !(popcount_u64(n) & 1))                       /* evil primes */
SPF_FILTER(gen_a007675, 1, squarefree_small(n) && squarefree_small(n + 1) && squarefree_small(n + 2))
SPF_FILTER(gen_a039955, 1, n % 4 == 1 && squarefree_small(n))
N_FILTER(gen_a001633, 0, ndigits10(n) & 1)
N_FILTER(gen_a001637, 10, !(ndigits10(n) & 1))
N_FILTER(gen_a034709, 1, n % 10 && n % (n % 10) == 0)
N_FILTER(gen_a038770, 1, div_by_a_digit(n))
N_FILTER(gen_a064150, 1, n % (u64)tern_sum(n) == 0)
SPF_FILTER(gen_a036785, 36, two_sq_primes(n))
N_FILTER(gen_a023709, 0, nodig_b(n, 4, 1))
N_FILTER(gen_a023713, 0, nodig_b(n, 4, 2))
N_FILTER(gen_a023721, 1, nodig_b(n, 5, 0))
N_FILTER(gen_a023725, 0, nodig_b(n, 5, 1))
N_FILTER(gen_a023729, 0, nodig_b(n, 5, 2))
N_FILTER(gen_a023733, 0, nodig_b(n, 5, 3))
N_FILTER(gen_a043493, 1, count_dig(n, 10, 1) == 1)
N_FILTER(gen_a023692, 1, count_dig(n, 3, 1) == 1)
SPF_FILTER(gen_a020893, 1, squarefree_small(n) && twosq_small(n))
SPF_FILTER(gen_a030634, 120, tau_small(n) == 16)
SPF_FILTER(gen_a030638, 240, tau_small(n) == 20)
static const u64 PS23[9] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
static void gen_a080682(u64 *t, long cnt) { if (smooth_all(t, cnt, PS23, 8) < cnt) { fprintf(stderr, "smooth19 short\n"); exit(1); } }
static void gen_a080683(u64 *t, long cnt) { if (smooth_all(t, cnt, PS23, 9) < cnt) { fprintf(stderr, "smooth23 short\n"); exit(1); } }

/* one wrapper per parametric sequence */
static void gen_a024913(u64 *t, long cnt) { gen_nlin(t, cnt, 10, -7); }
static void gen_a037030(u64 *t, long cnt) { gen_nlin(t, cnt, 666, 1); }
static void gen_a073085(u64 *t, long cnt) { gen_nlin(t, cnt, 210, 1); }
static void gen_a075745(u64 *t, long cnt) { gen_nlin(t, cnt, 210, 13); }
static void gen_a075746(u64 *t, long cnt) { gen_nlin(t, cnt, 210, -13); }
static void gen_a075747(u64 *t, long cnt) { gen_nlin(t, cnt, 210, 17); }
static void gen_a075748(u64 *t, long cnt) { gen_nlin(t, cnt, 210, -17); }
static void gen_a076354(u64 *t, long cnt) { gen_nlin(t, cnt, 210, -1); }
static void gen_a076355(u64 *t, long cnt) { gen_nlin(t, cnt, 210, 11); }
static void gen_a076356(u64 *t, long cnt) { gen_nlin(t, cnt, 210, -11); }
static void gen_a088958(u64 *t, long cnt) { gen_nlin(t, cnt, 60, 1); }
static void gen_a090614(u64 *t, long cnt) { gen_nlin(t, cnt, 14, 3); }
static void gen_a092022(u64 *t, long cnt) { gen_nlin(t, cnt, 16, 3); }
static void gen_a101084(u64 *t, long cnt) { gen_nlin(t, cnt, 97, 101); }
static void gen_a101503(u64 *t, long cnt) { gen_nlin(t, cnt, 11, 101); }
static void gen_a101557(u64 *t, long cnt) { gen_nlin(t, cnt, 101, 1009); }
static void gen_a102148(u64 *t, long cnt) { gen_nlin(t, cnt, 101, 11); }
static void gen_a102338(u64 *t, long cnt) { gen_nlin(t, cnt, 10, 3); }
static void gen_a102342(u64 *t, long cnt) { gen_nlin(t, cnt, 10, 7); }
static void gen_a102656(u64 *t, long cnt) { gen_nlin(t, cnt, 11, 1); }
static void gen_a014752(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 27, 1); }
static void gen_a033202(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 93, 1); }
static void gen_a033204(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 94, 1); }
static void gen_a033206(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 95, 1); }
static void gen_a033208(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 97, 1); }
static void gen_a033209(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 11, 1); }
static void gen_a033210(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 13, 1); }
static void gen_a033211(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 14, 1); }
static void gen_a033213(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 17, 1); }
static void gen_a033214(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 19, 1); }
static void gen_a033215(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 21, 1); }
static void gen_a033216(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 22, 1); }
static void gen_a051645(u64 *t, long cnt) { gen_plin(t, cnt, 30, 1); }
static void gen_a105961(u64 *t, long cnt) { gen_plin(t, cnt, 20, 3); }
static void gen_a112391(u64 *t, long cnt) { gen_plin(t, cnt, 23, 2); }
static void gen_a003625(u64 *t, long cnt) { gen_pres(t, cnt, 7, 104ULL); }
static void gen_a045372(u64 *t, long cnt) { gen_pres(t, cnt, 5, 6ULL); }
static void gen_a045429(u64 *t, long cnt) { gen_pres(t, cnt, 5, 10ULL); }
static void gen_a045378(u64 *t, long cnt) { gen_pres(t, cnt, 5, 20ULL); }
static void gen_a045435(u64 *t, long cnt) { gen_pres(t, cnt, 5, 24ULL); }
static void gen_a045321(u64 *t, long cnt) { gen_pres(t, cnt, 5, 14ULL); }
static void gen_a045371(u64 *t, long cnt) { gen_pres(t, cnt, 5, 22ULL); }
static void gen_a045428(u64 *t, long cnt) { gen_pres(t, cnt, 5, 26ULL); }
static void gen_a045327(u64 *t, long cnt) { gen_pres(t, cnt, 5, 28ULL); }
static void gen_a045392(u64 *t, long cnt) { gen_pres(t, cnt, 7, 4ULL); }
static void gen_a045437(u64 *t, long cnt) { gen_pres(t, cnt, 7, 8ULL); }
static void gen_a045471(u64 *t, long cnt) { gen_pres(t, cnt, 7, 16ULL); }
static void gen_a045458(u64 *t, long cnt) { gen_pres(t, cnt, 7, 32ULL); }
static void gen_a045473(u64 *t, long cnt) { gen_pres(t, cnt, 7, 64ULL); }
static void gen_a045465(u64 *t, long cnt) { gen_pres(t, cnt, 7, 3ULL); }
static void gen_a045391(u64 *t, long cnt) { gen_pres(t, cnt, 7, 6ULL); }
static void gen_a045436(u64 *t, long cnt) { gen_pres(t, cnt, 7, 10ULL); }
static void gen_a045343(u64 *t, long cnt) { gen_pres(t, cnt, 7, 12ULL); }
static void gen_a045469(u64 *t, long cnt) { gen_pres(t, cnt, 7, 18ULL); }
static void gen_a045387(u64 *t, long cnt) { gen_pres(t, cnt, 7, 20ULL); }
static void gen_a045432(u64 *t, long cnt) { gen_pres(t, cnt, 7, 24ULL); }
static void gen_a045456(u64 *t, long cnt) { gen_pres(t, cnt, 7, 34ULL); }
static void gen_a045368(u64 *t, long cnt) { gen_pres(t, cnt, 7, 36ULL); }
static void gen_a045416(u64 *t, long cnt) { gen_pres(t, cnt, 7, 40ULL); }
static void gen_a045452(u64 *t, long cnt) { gen_pres(t, cnt, 7, 48ULL); }
static void gen_a045472(u64 *t, long cnt) { gen_pres(t, cnt, 7, 66ULL); }
static void gen_a045389(u64 *t, long cnt) { gen_pres(t, cnt, 7, 68ULL); }
static void gen_a045434(u64 *t, long cnt) { gen_pres(t, cnt, 7, 72ULL); }
static void gen_a045467(u64 *t, long cnt) { gen_pres(t, cnt, 7, 80ULL); }
static void gen_a045455(u64 *t, long cnt) { gen_pres(t, cnt, 7, 96ULL); }
static void gen_a045342(u64 *t, long cnt) { gen_pres(t, cnt, 7, 14ULL); }
static void gen_a045386(u64 *t, long cnt) { gen_pres(t, cnt, 7, 22ULL); }
static void gen_a023203(u64 *t, long cnt) { gen_plin(t, cnt, 1, 10); }
static void gen_a046133(u64 *t, long cnt) { gen_plin(t, cnt, 1, 12); }
static void gen_a049488(u64 *t, long cnt) { gen_plin(t, cnt, 1, 16); }
static void gen_a049481(u64 *t, long cnt) { gen_plin(t, cnt, 1, 30); }
static void gen_a049489(u64 *t, long cnt) { gen_plin(t, cnt, 1, 32); }
static void gen_a062284(u64 *t, long cnt) { gen_plin(t, cnt, 1, 50); }
static void gen_a049482(u64 *t, long cnt) { gen_plin(t, cnt, 1, 210); }
static void gen_a063909(u64 *t, long cnt) { gen_plin(t, cnt, 2, -5); }
static void gen_a063910(u64 *t, long cnt) { gen_plin(t, cnt, 2, -7); }
static void gen_a063911(u64 *t, long cnt) { gen_plin(t, cnt, 2, -9); }
static void gen_a063912(u64 *t, long cnt) { gen_plin(t, cnt, 2, -11); }
static void gen_a063913(u64 *t, long cnt) { gen_plin(t, cnt, 2, -13); }
static void gen_a023209(u64 *t, long cnt) { gen_plin(t, cnt, 3, 4); }
static void gen_a023210(u64 *t, long cnt) { gen_plin(t, cnt, 3, 8); }
static void gen_a023211(u64 *t, long cnt) { gen_plin(t, cnt, 3, 10); }
static void gen_a062737(u64 *t, long cnt) { gen_plin(t, cnt, 4, -1); }
static void gen_a023213(u64 *t, long cnt) { gen_plin(t, cnt, 4, 3); }
static void gen_a023214(u64 *t, long cnt) { gen_plin(t, cnt, 4, 5); }
static void gen_a023215(u64 *t, long cnt) { gen_plin(t, cnt, 4, 7); }
static void gen_a023216(u64 *t, long cnt) { gen_plin(t, cnt, 4, 9); }
static void gen_a023217(u64 *t, long cnt) { gen_plin(t, cnt, 5, 2); }
static void gen_a023218(u64 *t, long cnt) { gen_plin(t, cnt, 5, 4); }
static void gen_a023220(u64 *t, long cnt) { gen_plin(t, cnt, 5, 8); }
static void gen_a007693(u64 *t, long cnt) { gen_plin(t, cnt, 6, 1); }
static void gen_a023221(u64 *t, long cnt) { gen_plin(t, cnt, 6, 5); }
static void gen_a023222(u64 *t, long cnt) { gen_plin(t, cnt, 6, 7); }
static void gen_a023223(u64 *t, long cnt) { gen_plin(t, cnt, 7, 2); }
static void gen_a023224(u64 *t, long cnt) { gen_plin(t, cnt, 7, 4); }
static void gen_a023225(u64 *t, long cnt) { gen_plin(t, cnt, 7, 6); }
static void gen_a023226(u64 *t, long cnt) { gen_plin(t, cnt, 7, 8); }
static void gen_a023227(u64 *t, long cnt) { gen_plin(t, cnt, 7, 10); }
static void gen_a023229(u64 *t, long cnt) { gen_plin(t, cnt, 8, 3); }
static void gen_a023231(u64 *t, long cnt) { gen_plin(t, cnt, 8, 7); }
static void gen_a023232(u64 *t, long cnt) { gen_plin(t, cnt, 8, 9); }
static void gen_a023233(u64 *t, long cnt) { gen_plin(t, cnt, 9, 2); }
static void gen_a023234(u64 *t, long cnt) { gen_plin(t, cnt, 9, 4); }
static void gen_a023235(u64 *t, long cnt) { gen_plin(t, cnt, 9, 8); }
static void gen_a023236(u64 *t, long cnt) { gen_plin(t, cnt, 9, 10); }
static void gen_a023237(u64 *t, long cnt) { gen_plin(t, cnt, 10, 1); }
static void gen_a023238(u64 *t, long cnt) { gen_plin(t, cnt, 10, 3); }
static void gen_a023239(u64 *t, long cnt) { gen_plin(t, cnt, 10, 7); }
static void gen_a023240(u64 *t, long cnt) { gen_plin(t, cnt, 10, 9); }
static void gen_a089443(u64 *t, long cnt) { gen_plin(t, cnt, 12, 13); }
static void gen_a113169(u64 *t, long cnt) { gen_plin(t, cnt, 13, 2); }
static void gen_a113115(u64 *t, long cnt) { gen_plin(t, cnt, 17, 2); }
static void gen_a067076(u64 *t, long cnt) { gen_nlin(t, cnt, 2, 3); }
static void gen_a098090(u64 *t, long cnt) { gen_nlin(t, cnt, 2, -3); }
static void gen_a089253(u64 *t, long cnt) { gen_nlin(t, cnt, 2, -5); }
static void gen_a089192(u64 *t, long cnt) { gen_nlin(t, cnt, 2, -7); }
static void gen_a102733(u64 *t, long cnt) { gen_nlin(t, cnt, 2, 101); }
static void gen_a024892(u64 *t, long cnt) { gen_nlin(t, cnt, 3, 1); }
static void gen_a087370(u64 *t, long cnt) { gen_nlin(t, cnt, 3, -1); }
static void gen_a024893(u64 *t, long cnt) { gen_nlin(t, cnt, 3, 2); }
static void gen_a034936(u64 *t, long cnt) { gen_nlin(t, cnt, 3, 4); }
static void gen_a089953(u64 *t, long cnt) { gen_nlin(t, cnt, 3, 7); }
static void gen_a005098(u64 *t, long cnt) { gen_nlin(t, cnt, 4, 1); }
static void gen_a095278(u64 *t, long cnt) { gen_nlin(t, cnt, 4, 3); }
static void gen_a111215(u64 *t, long cnt) { gen_nlin(t, cnt, 4, 5); }
static void gen_a111199(u64 *t, long cnt) { gen_nlin(t, cnt, 4, 9); }
static void gen_a024894(u64 *t, long cnt) { gen_nlin(t, cnt, 5, 1); }
static void gen_a024896(u64 *t, long cnt) { gen_nlin(t, cnt, 5, -2); }
static void gen_a111223(u64 *t, long cnt) { gen_nlin(t, cnt, 5, 2); }
static void gen_a024895(u64 *t, long cnt) { gen_nlin(t, cnt, 5, -3); }
static void gen_a087505(u64 *t, long cnt) { gen_nlin(t, cnt, 5, 3); }
static void gen_a024897(u64 *t, long cnt) { gen_nlin(t, cnt, 5, 4); }
static void gen_a081759(u64 *t, long cnt) { gen_nlin(t, cnt, 5, 6); }
static void gen_a107304(u64 *t, long cnt) { gen_nlin(t, cnt, 5, -7); }
static void gen_a111224(u64 *t, long cnt) { gen_nlin(t, cnt, 5, 7); }
static void gen_a111225(u64 *t, long cnt) { gen_nlin(t, cnt, 5, 8); }
static void gen_a111226(u64 *t, long cnt) { gen_nlin(t, cnt, 5, 12); }
static void gen_a111230(u64 *t, long cnt) { gen_nlin(t, cnt, 5, 14); }
static void gen_a024899(u64 *t, long cnt) { gen_nlin(t, cnt, 6, 1); }
static void gen_a059325(u64 *t, long cnt) { gen_nlin(t, cnt, 6, 5); }
static void gen_a024905(u64 *t, long cnt) { gen_nlin(t, cnt, 7, 1); }
static void gen_a024901(u64 *t, long cnt) { gen_nlin(t, cnt, 7, -2); }
static void gen_a105772(u64 *t, long cnt) { gen_nlin(t, cnt, 7, 2); }
static void gen_a089033(u64 *t, long cnt) { gen_nlin(t, cnt, 7, 3); }
static void gen_a024902(u64 *t, long cnt) { gen_nlin(t, cnt, 7, 4); }
static void gen_a024903(u64 *t, long cnt) { gen_nlin(t, cnt, 7, -4); }
static void gen_a024904(u64 *t, long cnt) { gen_nlin(t, cnt, 7, -5); }
static void gen_a111367(u64 *t, long cnt) { gen_nlin(t, cnt, 7, 5); }
static void gen_a024900(u64 *t, long cnt) { gen_nlin(t, cnt, 7, 6); }
static void gen_a111249(u64 *t, long cnt) { gen_nlin(t, cnt, 7, 8); }
static void gen_a111250(u64 *t, long cnt) { gen_nlin(t, cnt, 7, 10); }
static void gen_a033868(u64 *t, long cnt) { gen_nlin(t, cnt, 7, -11); }
static void gen_a089079(u64 *t, long cnt) { gen_nlin(t, cnt, 7, -23); }
static void gen_a108601(u64 *t, long cnt) { gen_nlin(t, cnt, 7, -911); }
static void gen_a108935(u64 *t, long cnt) { gen_nlin(t, cnt, 7, 911); }
static void gen_a005122(u64 *t, long cnt) { gen_nlin(t, cnt, 8, -1); }
static void gen_a005123(u64 *t, long cnt) { gen_nlin(t, cnt, 8, 1); }
static void gen_a005124(u64 *t, long cnt) { gen_nlin(t, cnt, 8, 3); }
static void gen_a005125(u64 *t, long cnt) { gen_nlin(t, cnt, 8, -3); }
static void gen_a105133(u64 *t, long cnt) { gen_nlin(t, cnt, 8, 5); }
static void gen_a024906(u64 *t, long cnt) { gen_nlin(t, cnt, 9, 1); }
static void gen_a024910(u64 *t, long cnt) { gen_nlin(t, cnt, 9, -2); }
static void gen_a024909(u64 *t, long cnt) { gen_nlin(t, cnt, 9, -4); }
static void gen_a024908(u64 *t, long cnt) { gen_nlin(t, cnt, 9, -5); }
static void gen_a024907(u64 *t, long cnt) { gen_nlin(t, cnt, 9, -7); }
static void gen_a024912(u64 *t, long cnt) { gen_nlin(t, cnt, 10, 1); }
static void gen_a105042(u64 *t, long cnt) { gen_nlin(t, cnt, 10, -1); }
static void gen_a024914(u64 *t, long cnt) { gen_nlin(t, cnt, 10, -3); }
static void gen_a005574(u64 *t, long cnt) { gen_nsq(t, cnt, 1); }
static void gen_a028870(u64 *t, long cnt) { gen_nsq(t, cnt, -2); }
static void gen_a067201(u64 *t, long cnt) { gen_nsq(t, cnt, 2); }
static void gen_a028873(u64 *t, long cnt) { gen_nsq(t, cnt, -3); }
static void gen_a049422(u64 *t, long cnt) { gen_nsq(t, cnt, 3); }
static void gen_a007591(u64 *t, long cnt) { gen_nsq(t, cnt, 4); }
static void gen_a028876(u64 *t, long cnt) { gen_nsq(t, cnt, -5); }
static void gen_a078402(u64 *t, long cnt) { gen_nsq(t, cnt, 5); }
static void gen_a028879(u64 *t, long cnt) { gen_nsq(t, cnt, -6); }
static void gen_a114269(u64 *t, long cnt) { gen_nsq(t, cnt, 6); }
static void gen_a028882(u64 *t, long cnt) { gen_nsq(t, cnt, -7); }
static void gen_a114270(u64 *t, long cnt) { gen_nsq(t, cnt, 7); }
static void gen_a028885(u64 *t, long cnt) { gen_nsq(t, cnt, -8); }
static void gen_a114271(u64 *t, long cnt) { gen_nsq(t, cnt, 8); }
static void gen_a114272(u64 *t, long cnt) { gen_nsq(t, cnt, 9); }
static void gen_a114273(u64 *t, long cnt) { gen_nsq(t, cnt, 10); }
static void gen_a114274(u64 *t, long cnt) { gen_nsq(t, cnt, 11); }
static void gen_a114275(u64 *t, long cnt) { gen_nsq(t, cnt, 12); }
static void gen_a113536(u64 *t, long cnt) { gen_nsq(t, cnt, 13); }
static void gen_a121250(u64 *t, long cnt) { gen_nsq(t, cnt, 14); }
static void gen_a121982(u64 *t, long cnt) { gen_nsq(t, cnt, 15); }
static void gen_a122062(u64 *t, long cnt) { gen_nsq(t, cnt, 16); }
static void gen_a020668(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 4, 0); }
static void gen_a020674(u64 *t, long cnt) { gen_bform(t, cnt, 2, 0, 5, 0); }
static void gen_a020677(u64 *t, long cnt) { gen_bform(t, cnt, 3, 0, 4, 0); }
static void gen_a020670(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 7, 0); }
static void gen_a020678(u64 *t, long cnt) { gen_bform(t, cnt, 3, 0, 5, 0); }
static void gen_a020671(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 8, 0); }
static void gen_a020675(u64 *t, long cnt) { gen_bform(t, cnt, 2, 0, 7, 0); }
static void gen_a020682(u64 *t, long cnt) { gen_bform(t, cnt, 4, 0, 5, 0); }
static void gen_a020672(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 9, 0); }
static void gen_a020679(u64 *t, long cnt) { gen_bform(t, cnt, 3, 0, 7, 0); }
static void gen_a020673(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 10, 0); }
static void gen_a020676(u64 *t, long cnt) { gen_bform(t, cnt, 2, 0, 9, 0); }
static void gen_a020680(u64 *t, long cnt) { gen_bform(t, cnt, 3, 0, 8, 0); }
static void gen_a020683(u64 *t, long cnt) { gen_bform(t, cnt, 4, 0, 7, 0); }
static void gen_a020685(u64 *t, long cnt) { gen_bform(t, cnt, 5, 0, 6, 0); }
static void gen_a020686(u64 *t, long cnt) { gen_bform(t, cnt, 5, 0, 7, 0); }
static void gen_a020681(u64 *t, long cnt) { gen_bform(t, cnt, 3, 0, 10, 0); }
static void gen_a020684(u64 *t, long cnt) { gen_bform(t, cnt, 4, 0, 9, 0); }
static void gen_a020687(u64 *t, long cnt) { gen_bform(t, cnt, 5, 0, 8, 0); }
static void gen_a020689(u64 *t, long cnt) { gen_bform(t, cnt, 6, 0, 7, 0); }
static void gen_a020688(u64 *t, long cnt) { gen_bform(t, cnt, 5, 0, 9, 0); }
static void gen_a020690(u64 *t, long cnt) { gen_bform(t, cnt, 7, 0, 8, 0); }
static void gen_a020691(u64 *t, long cnt) { gen_bform(t, cnt, 7, 0, 9, 0); }
static void gen_a020692(u64 *t, long cnt) { gen_bform(t, cnt, 7, 0, 10, 0); }
static void gen_a020693(u64 *t, long cnt) { gen_bform(t, cnt, 8, 0, 9, 0); }
static void gen_a020694(u64 *t, long cnt) { gen_bform(t, cnt, 9, 0, 10, 0); }
static void gen_a035121(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 82, 0); }
static void gen_a084865(u64 *t, long cnt) { gen_bform(t, cnt, 2, 0, 3, 1); }
static void gen_a106857(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 3, 1); }
static void gen_a106861(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 4, 1); }
static void gen_a106866(u64 *t, long cnt) { gen_bform(t, cnt, 2, 1, 3, 1); }
static void gen_a033199(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 6, 1); }
static void gen_a106862(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 5, 1); }
static void gen_a106871(u64 *t, long cnt) { gen_bform(t, cnt, 2, 1, 4, 1); }
static void gen_a106877(u64 *t, long cnt) { gen_bform(t, cnt, 3, 1, 3, 1); }
static void gen_a106889(u64 *t, long cnt) { gen_bform(t, cnt, 2, 0, 5, 1); }
static void gen_a106869(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 6, 1); }
static void gen_a106875(u64 *t, long cnt) { gen_bform(t, cnt, 3, 2, 3, 1); }
static void gen_a106885(u64 *t, long cnt) { gen_bform(t, cnt, 2, 1, 5, 1); }
static void gen_a106894(u64 *t, long cnt) { gen_bform(t, cnt, 3, 1, 4, 1); }
static void gen_a106870(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 7, 1); }
static void gen_a106882(u64 *t, long cnt) { gen_bform(t, cnt, 2, 2, 5, 1); }
static void gen_a106892(u64 *t, long cnt) { gen_bform(t, cnt, 3, 2, 4, 1); }
static void gen_a106897(u64 *t, long cnt) { gen_bform(t, cnt, 2, 1, 6, 1); }
static void gen_a106917(u64 *t, long cnt) { gen_bform(t, cnt, 2, 0, 7, 1); }
static void gen_a106918(u64 *t, long cnt) { gen_bform(t, cnt, 3, 1, 5, 1); }
static void gen_a106923(u64 *t, long cnt) { gen_bform(t, cnt, 4, 1, 4, 1); }
static void gen_a106963(u64 *t, long cnt) { gen_bform(t, cnt, 4, 0, 5, 1); }
static void gen_a102271(u64 *t, long cnt) { gen_bform(t, cnt, 3, 0, 7, 1); }
static void gen_a106874(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 8, 1); }
static void gen_a106883(u64 *t, long cnt) { gen_bform(t, cnt, 3, 3, 4, 1); }
static void gen_a106910(u64 *t, long cnt) { gen_bform(t, cnt, 2, 1, 7, 1); }
static void gen_a106914(u64 *t, long cnt) { gen_bform(t, cnt, 3, 2, 5, 1); }
static void gen_a106942(u64 *t, long cnt) { gen_bform(t, cnt, 3, 1, 6, 1); }
static void gen_a106956(u64 *t, long cnt) { gen_bform(t, cnt, 4, 1, 5, 1); }
static void gen_a033201(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 10, 1); }
static void gen_a047215(u64 *t, long cnt) { gen_res(t, cnt, 5, 5U, 0); }
static void gen_a047216(u64 *t, long cnt) { gen_res(t, cnt, 5, 6U, 1); }
static void gen_a047217(u64 *t, long cnt) { gen_res(t, cnt, 5, 7U, 0); }
static void gen_a047225(u64 *t, long cnt) { gen_res(t, cnt, 6, 3U, 0); }
static void gen_a047240(u64 *t, long cnt) { gen_res(t, cnt, 6, 7U, 0); }
static void gen_a047241(u64 *t, long cnt) { gen_res(t, cnt, 6, 10U, 1); }
static void gen_a047274(u64 *t, long cnt) { gen_res(t, cnt, 7, 3U, 0); }
static void gen_a047352(u64 *t, long cnt) { gen_res(t, cnt, 7, 5U, 0); }
static void gen_a047353(u64 *t, long cnt) { gen_res(t, cnt, 7, 6U, 1); }
static void gen_a008590(u64 *t, long cnt) { gen_res(t, cnt, 8, 1U, 0); }
static void gen_a047393(u64 *t, long cnt) { gen_res(t, cnt, 8, 3U, 0); }
static void gen_a047467(u64 *t, long cnt) { gen_res(t, cnt, 8, 5U, 0); }
static void gen_a090570(u64 *t, long cnt) { gen_res(t, cnt, 9, 3U, 0); }
static void gen_a087444(u64 *t, long cnt) { gen_res(t, cnt, 9, 18U, 1); }
static void gen_a054966(u64 *t, long cnt) { gen_res(t, cnt, 9, 259U, 0); }
static void gen_a090773(u64 *t, long cnt) { gen_res(t, cnt, 10, 80U, 4); }
static void gen_a078309(u64 *t, long cnt) { gen_res(t, cnt, 10, 146U, 1); }
static void gen_a090772(u64 *t, long cnt) { gen_res(t, cnt, 10, 260U, 2); }
static void gen_a008593(u64 *t, long cnt) { gen_res(t, cnt, 11, 1U, 0); }
static void gen_a008594(u64 *t, long cnt) { gen_res(t, cnt, 12, 1U, 0); }
static void gen_a083031(u64 *t, long cnt) { gen_res(t, cnt, 12, 137U, 0); }
static void gen_a083030(u64 *t, long cnt) { gen_res(t, cnt, 12, 145U, 0); }
static void gen_a008595(u64 *t, long cnt) { gen_res(t, cnt, 13, 1U, 0); }
static void gen_a092476(u64 *t, long cnt) { gen_res(t, cnt, 13, 522U, 1); }
static void gen_a008596(u64 *t, long cnt) { gen_res(t, cnt, 14, 1U, 0); }
static void gen_a113806(u64 *t, long cnt) { gen_res(t, cnt, 14, 320U, 6); }
static void gen_a113805(u64 *t, long cnt) { gen_res(t, cnt, 14, 544U, 5); }
static void gen_a008597(u64 *t, long cnt) { gen_res(t, cnt, 15, 1U, 0); }
static void gen_a087446(u64 *t, long cnt) { gen_res(t, cnt, 15, 66U, 1); }
static void gen_a008598(u64 *t, long cnt) { gen_res(t, cnt, 16, 1U, 0); }
static void gen_a106839(u64 *t, long cnt) { gen_res(t, cnt, 16, 2048U, 11); }
static void gen_a008599(u64 *t, long cnt) { gen_res(t, cnt, 17, 1U, 0); }
static void gen_a008600(u64 *t, long cnt) { gen_res(t, cnt, 18, 1U, 0); }
static void gen_a008601(u64 *t, long cnt) { gen_res(t, cnt, 19, 1U, 0); }
static void gen_a008602(u64 *t, long cnt) { gen_res(t, cnt, 20, 1U, 0); }
static void gen_a008603(u64 *t, long cnt) { gen_res(t, cnt, 21, 1U, 0); }
static void gen_a008604(u64 *t, long cnt) { gen_res(t, cnt, 22, 1U, 0); }
static void gen_a008605(u64 *t, long cnt) { gen_res(t, cnt, 23, 1U, 0); }
static void gen_a008606(u64 *t, long cnt) { gen_res(t, cnt, 24, 1U, 0); }
static void gen_a033430(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)4 * (N*N*N))); } }  /* 4*n^3 */
static void gen_a033431(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)2 * (N*N*N))); } }  /* 2*n^3 */
static void gen_a084377(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N*N) + (i128)7)); } }  /* n^3+7 */
static void gen_a084378(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N*N) + (i128)3)); } }  /* n^3+3 */
static void gen_a084380(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N*N) + (i128)2)); } }  /* n^3+2 */
static void gen_a084381(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N*N) + (i128)5)); } }  /* n^3+5 */
static void gen_a084382(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N*N) + (i128)6)); } }  /* n^3+6 */
static void gen_a117642(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)3 * (N*N*N))); } }  /* 3*n^3 */
static void gen_a084379(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N*N) + (i128)17)); } }  /* n^3+17 */
static void gen_a033562(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((i128)2 * (N*N*N)) + (i128)1)); } }  /* 2*n^3+1 */
static void gen_a100214(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((i128)4 * (N*N*N)) + (i128)4)); } }  /* 4*n^3+4 */
static void gen_a118465(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((i128)8 * (N*N*N)) + N)); } }  /* 8*n^3+n */
static void gen_a003777(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((((N*N*N) + (N*N)) - (i128)1)); } }  /* n^3+n^2-1 */
static void gen_a005491(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N*N*N) + ((i128)3 * N)) + (i128)1)); } }  /* n^3+3*n+1 */
static void gen_a011379(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N) * (N + (i128)1))); } }  /* n^2*(n+1) */
static void gen_a027444(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N*N*N) + (N*N)) + N)); } }  /* n^3+n^2+n */
static void gen_a098547(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N*N*N) + (N*N)) + (i128)1)); } }  /* n^3+n^2+1 */
static void gen_a105374(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((i128)4 * (N*N*N)) + ((i128)4 * N))); } }  /* 4*n^3+4*n */
static void gen_a114364(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((N * ((N + (i128)1)*(N + (i128)1)))); } }  /* n*(n+1)^2 */
static void gen_a119536(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((i128)3 * (N*N*N)) + ((i128)3 * N))); } }  /* 3*n^3+3*n */
static void gen_a122562(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)(((N*N*N) + ((i128)114 * N))); } }  /* n^3+114*n */
static void gen_a006002(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * ((N + (i128)1)*(N + (i128)1))) / (i128)2)); } }  /* n*(n+1)^2/2 */
static void gen_a006527(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N*N*N) + ((i128)2 * N)) / (i128)3)); } }  /* (n^3+2*n)/3 */
static void gen_a015237(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((((i128)2 * N) - (i128)1) * (N*N))); } }  /* (2*n-1)*n^2 */
static void gen_a053698(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((((N*N*N) + (N*N)) + N) + (i128)1)); } }  /* n^3+n^2+n+1 */
static void gen_a084367(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N * ((((i128)2 * N) + (i128)1)*(((i128)2 * N) + (i128)1)))); } }  /* n*(2*n+1)^2 */
static void gen_a089207(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((((i128)4 * (N*N*N)) + ((i128)2 * (N*N)))); } }  /* 4*n^3+2*n^2 */
static void gen_a099721(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N) * (((i128)2 * N) + (i128)1))); } }  /* n^2*(2*n+1) */
static void gen_a100109(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((((N*N*N) - ((i128)2 * (N*N))) + (i128)2)); } }  /* n^3-2*n^2+2 */
static void gen_a100705(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N*N) + ((N + (i128)1)*(N + (i128)1)))); } }  /* n^3+(n+1)^2 */
static void gen_a028347(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 2; t[i] = (u64)(((N*N) - (i128)4)); } }  /* n^2-4 */
static void gen_a028872(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 2; t[i] = (u64)(((N*N) - (i128)3)); } }  /* n^2-3 */
static void gen_a028881(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 3; t[i] = (u64)(((N*N) - (i128)7)); } }  /* n^2-7 */
static void gen_a033428(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)3 * (N*N))); } }  /* 3*n^2 */
static void gen_a033429(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)5 * (N*N))); } }  /* 5*n^2 */
static void gen_a033581(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)6 * (N*N))); } }  /* 6*n^2 */
static void gen_a033582(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)7 * (N*N))); } }  /* 7*n^2 */
static void gen_a059100(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N) + (i128)2)); } }  /* n^2+2 */
static void gen_a087475(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N) + (i128)4)); } }  /* n^2+4 */
static void gen_a114949(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N) + (i128)6)); } }  /* n^2+6 */
static void gen_a117619(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N) + (i128)7)); } }  /* n^2+7 */
static void gen_a117950(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N) + (i128)3)); } }  /* n^2+3 */
static void gen_a117951(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N) + (i128)5)); } }  /* n^2+5 */
static void gen_a033583(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)10 * (N*N))); } }  /* 10*n^2 */
static void gen_a033584(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)11 * (N*N))); } }  /* 11*n^2 */
static void gen_a064761(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)15 * (N*N))); } }  /* 15*n^2 */
static void gen_a064762(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)21 * (N*N))); } }  /* 21*n^2 */
static void gen_a064763(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)28 * (N*N))); } }  /* 28*n^2 */
static void gen_a114948(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N) + (i128)10)); } }  /* n^2+10 */
static void gen_a114962(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N) + (i128)14)); } }  /* n^2+14 */
static void gen_a114963(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N) + (i128)22)); } }  /* n^2+22 */
static void gen_a114964(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N) + (i128)30)); } }  /* n^2+30 */
static void gen_a114965(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N*N) + (i128)34)); } }  /* n^2+34 */
static void gen_a016766(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((i128)3 * N)*((i128)3 * N))); } }  /* (3*n)^2 */
static void gen_a016802(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((i128)4 * N)*((i128)4 * N))); } }  /* (4*n)^2 */
static void gen_a016850(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((i128)5 * N)*((i128)5 * N))); } }  /* (5*n)^2 */
static void gen_a016910(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((i128)6 * N)*((i128)6 * N))); } }  /* (6*n)^2 */
static void gen_a016982(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((i128)7 * N)*((i128)7 * N))); } }  /* (7*n)^2 */
static void gen_a017066(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((i128)8 * N)*((i128)8 * N))); } }  /* (8*n)^2 */
static void gen_a017162(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((i128)9 * N)*((i128)9 * N))); } }  /* (9*n)^2 */
static void gen_a027688(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N*N) + N) + (i128)3)); } }  /* n^2+n+3 */
static void gen_a027689(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N*N) + N) + (i128)4)); } }  /* n^2+n+4 */
static void gen_a004919(u64 *t, long cnt) { gen_beattyl(t, cnt, 6.85410196624968454461376050309691435L, 0); }
static void gen_a004920(u64 *t, long cnt) { gen_beattyl(t, cnt, 11.0901699437494742410229341718281906L, 0); }
static void gen_a004921(u64 *t, long cnt) { gen_beattyl(t, cnt, 17.9442719099991587856366946749251049L, 0); }
static void gen_a004922(u64 *t, long cnt) { gen_beattyl(t, cnt, 29.0344418537486330266596288467532955L, 0); }
static void gen_a004976(u64 *t, long cnt) { gen_beattyl(t, cnt, 4.23606797749978969640917366873127624L, 0); }
static void gen_a037085(u64 *t, long cnt) { gen_beattyl(t, cnt, 9.86960440108935861883449099987615114L, 0); }
static void gen_a037086(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.77245385090551602729816748334114518L, 0); }
static void gen_a037087(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.44466786100976613365833910859643022L, 0); }
static void gen_a038130(u64 *t, long cnt) { gen_beattyl(t, cnt, 6.28318530717958647692528676655900577L, 0); }
static void gen_a038152(u64 *t, long cnt) { gen_beattyl(t, cnt, 23.1406926327792690057290863679485474L, 1); }
static void gen_a038153(u64 *t, long cnt) { gen_beattyl(t, cnt, 22.459157718361045473427152204543735L, 1); }
static void gen_a054386(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.46694220692425985998339481323366757L, 1); }
static void gen_a054965(u64 *t, long cnt) { gen_beattyl(t, cnt, 2.09590327428938460429656752202140125L, 1); }
static void gen_a059531(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.31830988618379067153776752674502872L, 1); }
static void gen_a059532(u64 *t, long cnt) { gen_beattyl(t, cnt, 4.14159265358979323846264338327950288L, 1); }
static void gen_a059535(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.64493406684822643647241516664602519L, 1); }
static void gen_a059536(u64 *t, long cnt) { gen_beattyl(t, cnt, 2.5505460967304304402864869634760961L, 1); }
static void gen_a059537(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.20205690315959428539973816151144999L, 1); }
static void gen_a059538(u64 *t, long cnt) { gen_beattyl(t, cnt, 5.9491008936732628209342103204674846L, 1); }
static void gen_a059539(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.44224957030740830177251154964324087L, 1); }
static void gen_a059540(u64 *t, long cnt) { gen_beattyl(t, cnt, 3.26116669667965686230104438436683267L, 1); }
static void gen_a059541(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.69314718055994530941723212145817657L, 1); }
static void gen_a059542(u64 *t, long cnt) { gen_beattyl(t, cnt, 2.44269504088896340735992468100189214L, 1); }
static void gen_a059543(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.0986122886681096913952452369225257L, 1); }
static void gen_a059544(u64 *t, long cnt) { gen_beattyl(t, cnt, 11.1407239757471607801981033707929266L, 1); }
static void gen_a059545(u64 *t, long cnt) { gen_beattyl(t, cnt, 2.30258509299404568401799145468436421L, 1); }
static void gen_a059546(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.76770416411065987316179013970470072L, 1); }
static void gen_a059547(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.910239226626837393614240165736107L, 1); }
static void gen_a059548(u64 *t, long cnt) { gen_beattyl(t, cnt, 2.0986122886681096913952452369225257L, 1); }
static void gen_a059549(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.43429448190325182765112891891660508L, 1); }

/* ---------------------------- the two hundred added 26 Sep 2026 (801-1000) */

static void gen_pk2(u64 *t, long cnt, long long c)            /* primes k^2 + c, in order of k */
{ long k = 0; for (u64 n = 0; k < cnt; n++) { i128 v = (i128)n * n + c; if (isprime_s(v)) t[k++] = (u64)v; } }
static void gen_nn2(u64 *t, long cnt, long long c)            /* n >= 0 with n^2 + n + c prime */
{ long k = 0; for (u64 n = 0; k < cnt; n++) if (isprime_s((i128)n * n + n + c)) t[k++] = n; }
static void gen_kbits(u64 *t, long cnt, int bits)             /* exactly `bits` ones in binary (Gosper's hack) */
{
    u64 v = (1ULL << bits) - 1;
    for (long i = 0; i < cnt; i++) { t[i] = v; u64 c = v & -v, r = v + c; v = (((r ^ v) >> 2) / c) | r; }
}
static void gen_psq_c(u64 *t, long cnt, long long c)          /* primes p with p^2 + c prime */
{
    for (u64 lim = 1ULL << 22; ; lim *= 2) {
        u8 *cf = composite_flags(lim); long k = 0;
        for (u64 p = 2; p <= lim && k < cnt; p++) if (!cf[p] && isprime_s((i128)p * p + c)) t[k++] = p;
        free(cf); if (k == cnt) return;
    }
}
static void gen_pshift_semi(u64 *t, long cnt, int s)          /* primes p with p + s a semiprime */
{
    for (u64 lim = 1ULL << 22; ; lim *= 2) {
        u8 *cf = composite_flags(lim); long k = 0;
        for (u64 p = 2; p <= lim && k < cnt; p++) {
            if (cf[p] || (long long)p + s < 4) continue;
            u64 q = (u64)((long long)p + s); need_spf(q, "pshift_semi");
            if (omega_big(q) == 2) t[k++] = p;
        }
        free(cf); if (k == cnt) return;
    }
}
static void gen_consec_psum(u64 *t, long cnt, int j)          /* sums of j consecutive primes */
{
    u64 *p = malloc((size_t)(cnt + j) * sizeof(u64)); gen_primes(p, cnt + j);
    u64 s = 0; for (int i = 0; i < j; i++) s += p[i];
    for (long i = 0; i < cnt; i++) { t[i] = s; s += p[i + j] - p[i]; }
    free(p);
}
static void gen_powsieve(u64 *t, long cnt, u64 b)             /* delete every b-th, then every b^2-th, ... */
{
    /* inverting: a survivor at q after the step with m sat at q + floor((q-1)/(m-1)) before it;
       steps with m > the position never move it */
    for (long i = 0; i < cnt; i++) {
        u64 q = (u64)(i + 1), ms[64]; int nm = 0;
        for (u64 m = b; ; m *= b) { ms[nm++] = m; if (m > 4 * q + 64 || m > (1ULL << 60) / b) break; }
        for (int s = nm - 1; s >= 0; s--) q += (q - 1) / (ms[s] - 1);
        t[i] = q;
    }
}
static void gen_a031879(u64 *t, long cnt)                     /* nonprime lucky numbers */
{
    for (long L = cnt + cnt / 4; ; L += L / 2) {
        u64 *l = malloc((size_t)L * sizeof(u64)); gen_lucky(l, L);
        long k = 0; for (long i = 0; i < L && k < cnt; i++) if (!is_prime(l[i])) t[k++] = l[i];
        free(l); if (k == cnt) return;
    }
}
SPF_FILTER(gen_a005237, 1, tau_small(n) == tau_small(n + 1))
SPF_FILTER(gen_a006049, 1, omega_distinct(n) == omega_distinct(n + 1))
SPF_FILTER(gen_a045920, 1, omega_big(n) == omega_big(n + 1))
SPF_FILTER(gen_a063464, 1, omega_distinct(n) == omega_distinct(n + 2))
SPF_FILTER(gen_a069977, 1, squarefree_small(n) && squarefree_small(n + 2))
SPF_FILTER(gen_a121495, 2, squarefree_small(n) && squarefree_small(n + 1) && spf[n] != n && spf[n + 1] != n + 1)
SPF_FILTER(gen_a124940, 8, omega_big(n) == 3 && omega_big(n + 3) == 3)
SPF_FILTER(gen_a124941, 16, omega_big(n) == 4 && omega_big(n + 4) == 4)
SPF_FILTER(gen_a070552, 4, omega_big(n) == 2 && omega_big(n + 1) == 2)

static void gen_a002496(u64 *t, long cnt) { gen_pk2(t, cnt, 1); }
static void gen_a005473(u64 *t, long cnt) { gen_pk2(t, cnt, 4); }
static void gen_a028871(u64 *t, long cnt) { gen_pk2(t, cnt, -2); }
static void gen_a028874(u64 *t, long cnt) { gen_pk2(t, cnt, -3); }
static void gen_a028877(u64 *t, long cnt) { gen_pk2(t, cnt, -5); }
static void gen_a028880(u64 *t, long cnt) { gen_pk2(t, cnt, -6); }
static void gen_a028883(u64 *t, long cnt) { gen_pk2(t, cnt, -7); }
static void gen_a028886(u64 *t, long cnt) { gen_pk2(t, cnt, -8); }
static void gen_a049423(u64 *t, long cnt) { gen_pk2(t, cnt, 3); }
static void gen_a056899(u64 *t, long cnt) { gen_pk2(t, cnt, 2); }
static void gen_a056905(u64 *t, long cnt) { gen_pk2(t, cnt, 5); }
static void gen_a056909(u64 *t, long cnt) { gen_pk2(t, cnt, 6); }
static void gen_a079138(u64 *t, long cnt) { gen_pk2(t, cnt, 7); }
static void gen_a091272(u64 *t, long cnt) { gen_pk2(t, cnt, -11); }
static void gen_a027752(u64 *t, long cnt) { gen_nn2(t, cnt, 3); }
static void gen_a027754(u64 *t, long cnt) { gen_nn2(t, cnt, 5); }
static void gen_a027756(u64 *t, long cnt) { gen_nn2(t, cnt, 7); }
static void gen_a027757(u64 *t, long cnt) { gen_nn2(t, cnt, 9); }
static void gen_a028823(u64 *t, long cnt) { gen_nn2(t, cnt, 17); }
static void gen_a045546(u64 *t, long cnt) { gen_nn2(t, cnt, -1); }
static void gen_a048097(u64 *t, long cnt) { gen_nn2(t, cnt, 11); }
static void gen_a014312(u64 *t, long cnt) { gen_kbits(t, cnt, 4); }
static void gen_a014313(u64 *t, long cnt) { gen_kbits(t, cnt, 5); }
static void gen_a023688(u64 *t, long cnt) { gen_kbits(t, cnt, 6); }
static void gen_a023689(u64 *t, long cnt) { gen_kbits(t, cnt, 7); }
static void gen_a023690(u64 *t, long cnt) { gen_kbits(t, cnt, 8); }
static void gen_a023691(u64 *t, long cnt) { gen_kbits(t, cnt, 9); }
static void gen_a062324(u64 *t, long cnt) { gen_psq_c(t, cnt, 4); }
static void gen_a062326(u64 *t, long cnt) { gen_psq_c(t, cnt, -2); }
static void gen_a063637(u64 *t, long cnt) { gen_pshift_semi(t, cnt, 2); }
static void gen_a063638(u64 *t, long cnt) { gen_pshift_semi(t, cnt, -2); }
static void gen_a127333(u64 *t, long cnt) { gen_consec_psum(t, cnt, 6); }
static void gen_a127334(u64 *t, long cnt) { gen_consec_psum(t, cnt, 7); }
static void gen_a127336(u64 *t, long cnt) { gen_consec_psum(t, cnt, 9); }
static void gen_a127337(u64 *t, long cnt) { gen_consec_psum(t, cnt, 10); }
static void gen_a127338(u64 *t, long cnt) { gen_consec_psum(t, cnt, 11); }
static void gen_a127339(u64 *t, long cnt) { gen_consec_psum(t, cnt, 12); }
static void gen_a007950(u64 *t, long cnt) { gen_powsieve(t, cnt, 2); }
static void gen_a007951(u64 *t, long cnt) { gen_powsieve(t, cnt, 3); }
static void gen_a059550(u64 *t, long cnt) { gen_beattyl(t, cnt, 3.30258509299404568401799145468436421L, 1); }
static void gen_a059551(u64 *t, long cnt) { gen_beattyl(t, cnt, 2.67893853470774778891161190097964167L, 1); }
static void gen_a059552(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.59561441906750304645822162264394191L, 1); }
static void gen_a059553(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.35411793942640048300521854600597076L, 1); }
static void gen_a059554(u64 *t, long cnt) { gen_beattyl(t, cnt, 3.8239179342899090713720768726254423L, 1); }
static void gen_a059556(u64 *t, long cnt) { gen_beattyl(t, cnt, 2.73245471460063347358302531586082968L, 1); }
static void gen_a059557(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.33317792380771867431837613635524423L, 1); }
static void gen_a059559(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.54953931298164482233766176880290779L, 1); }
static void gen_a059560(u64 *t, long cnt) { gen_beattyl(t, cnt, 2.81970602717080047632906605775390654L, 1); }
static void gen_a059562(u64 *t, long cnt) { gen_beattyl(t, cnt, 7.90942298566142654468760580886779175L, 1); }
static void gen_a059563(u64 *t, long cnt) { gen_beattyl(t, cnt, 3.08616126963048755695581124151412337L, 1); }
static void gen_a059564(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.47934932670719437753878078263261546L, 1); }
static void gen_a059566(u64 *t, long cnt) { gen_beattyl(t, cnt, 2.28029101651436042828674681232510902L, 1); }
static void gen_a059567(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.36651292058166432701243915823266947L, 1); }
static void gen_a059568(u64 *t, long cnt) { gen_beattyl(t, cnt, 3.7284167729011498257169259985034365L, 1); }
static void gen_a066343(u64 *t, long cnt) { gen_beattyl(t, cnt, 3.32192809488736234787031942948939018L, 1); }
static void gen_a066344(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.43067655807339305067010656876396563L, 1); }
static void gen_a098005(u64 *t, long cnt) { gen_beattyl(t, cnt, 3.5496467783038448822263926847976956L, 1); }
static void gen_a108598(u64 *t, long cnt) { gen_beattyl(t, cnt, 1.80901699437494742410229341718281906L, 1); }
static void gen_a121283(u64 *t, long cnt) { gen_beattyl(t, cnt, 8.5397342226735670654635508695465745L, 0); }
static void gen_a102700(u64 *t, long cnt) { gen_nlin(t, cnt, 10, 9); }
static void gen_a102703(u64 *t, long cnt) { gen_nlin(t, cnt, 100, 99); }
static void gen_a102711(u64 *t, long cnt) { gen_nlin(t, cnt, 11, 7); }
static void gen_a102721(u64 *t, long cnt) { gen_nlin(t, cnt, 11, 13); }
static void gen_a102731(u64 *t, long cnt) { gen_nlin(t, cnt, 11, 23); }
static void gen_a102768(u64 *t, long cnt) { gen_nlin(t, cnt, 23, 11); }
static void gen_a103118(u64 *t, long cnt) { gen_nlin(t, cnt, 100, 57); }
static void gen_a103871(u64 *t, long cnt) { gen_nlin(t, cnt, 100, 69); }
static void gen_a105043(u64 *t, long cnt) { gen_nlin(t, cnt, 100, -1); }
static void gen_a105044(u64 *t, long cnt) { gen_nlin(t, cnt, 1000, -1); }
static void gen_a105134(u64 *t, long cnt) { gen_nlin(t, cnt, 16, 9); }
static void gen_a105135(u64 *t, long cnt) { gen_nlin(t, cnt, 32, 17); }
static void gen_a105136(u64 *t, long cnt) { gen_nlin(t, cnt, 64, 33); }
static void gen_a105137(u64 *t, long cnt) { gen_nlin(t, cnt, 128, 65); }
static void gen_a105138(u64 *t, long cnt) { gen_nlin(t, cnt, 256, 129); }
static void gen_a105139(u64 *t, long cnt) { gen_nlin(t, cnt, 512, 257); }
static void gen_a105583(u64 *t, long cnt) { gen_nlin(t, cnt, 101, 997); }
static void gen_a105679(u64 *t, long cnt) { gen_nlin(t, cnt, 997, 101); }
static void gen_a105773(u64 *t, long cnt) { gen_nlin(t, cnt, 11, 97); }
static void gen_a105775(u64 *t, long cnt) { gen_nlin(t, cnt, 97, 11); }
static void gen_a106690(u64 *t, long cnt) { gen_nlin(t, cnt, 11, -97); }
static void gen_a106692(u64 *t, long cnt) { gen_nlin(t, cnt, 97, -11); }
static void gen_a106695(u64 *t, long cnt) { gen_nlin(t, cnt, 101, -997); }
static void gen_a106697(u64 *t, long cnt) { gen_nlin(t, cnt, 997, -101); }
static void gen_a107305(u64 *t, long cnt) { gen_nlin(t, cnt, 11, -13); }
static void gen_a107366(u64 *t, long cnt) { gen_nlin(t, cnt, 101, 103); }
static void gen_a107369(u64 *t, long cnt) { gen_nlin(t, cnt, 103, 101); }
static void gen_a107371(u64 *t, long cnt) { gen_nlin(t, cnt, 101, -103); }
static void gen_a107372(u64 *t, long cnt) { gen_nlin(t, cnt, 103, -101); }
static void gen_a107400(u64 *t, long cnt) { gen_nlin(t, cnt, 107, 109); }
static void gen_a107405(u64 *t, long cnt) { gen_nlin(t, cnt, 109, 107); }
static void gen_a107406(u64 *t, long cnt) { gen_nlin(t, cnt, 107, -109); }
static void gen_a107407(u64 *t, long cnt) { gen_nlin(t, cnt, 109, -107); }
static void gen_a107960(u64 *t, long cnt) { gen_nlin(t, cnt, 11, -1); }
static void gen_a107992(u64 *t, long cnt) { gen_nlin(t, cnt, 11, -3); }
static void gen_a107994(u64 *t, long cnt) { gen_nlin(t, cnt, 11, -2); }
static void gen_a108027(u64 *t, long cnt) { gen_nlin(t, cnt, 137, 139); }
static void gen_a108028(u64 *t, long cnt) { gen_nlin(t, cnt, 139, 137); }
static void gen_a108029(u64 *t, long cnt) { gen_nlin(t, cnt, 149, 151); }
static void gen_a108030(u64 *t, long cnt) { gen_nlin(t, cnt, 151, 149); }
static void gen_a033217(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 23, 1); }
static void gen_a033218(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 26, 1); }
static void gen_a033219(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 29, 1); }
static void gen_a033220(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 30, 1); }
static void gen_a033221(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 31, 1); }
static void gen_a033222(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 33, 1); }
static void gen_a033223(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 34, 1); }
static void gen_a033224(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 35, 1); }
static void gen_a033225(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 37, 1); }
static void gen_a033226(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 38, 1); }
static void gen_a033227(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 39, 1); }
static void gen_a033228(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 41, 1); }
static void gen_a033229(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 42, 1); }
static void gen_a033230(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 43, 1); }
static void gen_a033231(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 46, 1); }
static void gen_a033232(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 47, 1); }
static void gen_a033233(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 51, 1); }
static void gen_a033234(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 53, 1); }
static void gen_a033235(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 55, 1); }
static void gen_a033236(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 57, 1); }
static void gen_a033237(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 58, 1); }
static void gen_a033238(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 59, 1); }
static void gen_a033239(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 61, 1); }
static void gen_a033240(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 62, 1); }
static void gen_a033241(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 65, 1); }
static void gen_a033242(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 66, 1); }
static void gen_a033243(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 67, 1); }
static void gen_a033244(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 69, 1); }
static void gen_a033245(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 70, 1); }
static void gen_a033246(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 71, 1); }
static void gen_a033247(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 73, 1); }
static void gen_a033248(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 74, 1); }
static void gen_a033249(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 77, 1); }
static void gen_a033250(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 78, 1); }
static void gen_a033251(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 79, 1); }
static void gen_a033252(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 82, 1); }
static void gen_a127989(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)(((((i128)2 * (N*N*N)) - ((i128)2 * N)) + (i128)9)); } }  /* 2*n^3-2*n+9 */
static void gen_a004126(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (((i128)7 * (N*N)) - (i128)1)) / (i128)6)); } }  /* n*(7*n^2-1)/6 */
static void gen_a004188(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (((i128)3 * (N*N)) - (i128)1)) / (i128)2)); } }  /* n*(3*n^2-1)/2 */
static void gen_a004466(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (((i128)5 * (N*N)) - (i128)2)) / (i128)3)); } }  /* n*(5*n^2-2)/3 */
static void gen_a006597(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N*N) * (((i128)5 * N) - (i128)3)) / (i128)2)); } }  /* n^2*(5*n-3)/2 */
static void gen_a061804(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((i128)2 * N) * (((i128)2 * (N*N)) + (i128)1))); } }  /* 2*n*(2*n^2+1) */
static void gen_a063521(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (((i128)7 * (N*N)) - (i128)4)) / (i128)3)); } }  /* n*(7*n^2-4)/3 */
static void gen_a063522(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (((i128)5 * (N*N)) - (i128)3)) / (i128)2)); } }  /* n*(5*n^2-3)/2 */
static void gen_a063523(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (((i128)8 * (N*N)) - (i128)5)) / (i128)3)); } }  /* n*(8*n^2-5)/3 */
static void gen_a067389(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((((i128)3 * (N*N*N)) + ((i128)2 * (N*N))) + N)); } }  /* 3*n^3+2*n^2+n */
static void gen_a117560(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 2; t[i] = (u64)((((N * ((N*N) - (i128)1)) / (i128)2) - (i128)1)); } }  /* n*(n^2-1)/2-1 */
static void gen_a004467(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (((i128)11 * (N*N)) - (i128)5)) / (i128)6)); } }  /* n*(11*n^2-5)/6 */
static void gen_a060163(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + -2; t[i] = (u64)(((((N*N*N) + ((i128)5 * N)) + (i128)18) / (i128)6)); } }  /* (n^3+5*n+18)/6 */
static void gen_a062025(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (((i128)13 * (N*N)) - (i128)7)) / (i128)6)); } }  /* n*(13*n^2-7)/6 */
static void gen_a005586(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N * (N + (i128)4)) * (N + (i128)5)) / (i128)6)); } }  /* n*(n+4)*(n+5)/6 */
static void gen_a006503(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N * (N + (i128)1)) * (N + (i128)8)) / (i128)6)); } }  /* n*(n+1)*(n+8)/6 */
static void gen_a027480(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N * (N + (i128)1)) * (N + (i128)2)) / (i128)2)); } }  /* n*(n+1)*(n+2)/2 */
static void gen_a027903(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (N + (i128)1)) * (((i128)3 * N) + (i128)1))); } }  /* n*(n+1)*(3*n+1) */
static void gen_a055112(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (N + (i128)1)) * (((i128)2 * N) + (i128)1))); } }  /* n*(n+1)*(2*n+1) */
static void gen_a059722(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N * ((((i128)2 * (N*N)) - ((i128)2 * N)) + (i128)1))); } }  /* n*(2*n^2-2*n+1) */
static void gen_a071233(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((((i128)2 * (N - (i128)1)) * ((N*N) + (i128)1))); } }  /* 2*(n-1)*(n^2+1) */
static void gen_a077414(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((((N * (N - (i128)1)) * (N + (i128)2)) / (i128)2)); } }  /* n*(n-1)*(n+2)/2 */
static void gen_a077415(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 2; t[i] = (u64)((((N * (N + (i128)2)) * (N - (i128)2)) / (i128)3)); } }  /* n*(n+2)*(n-2)/3 */
static void gen_a084990(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (((N*N) + ((i128)3 * N)) - (i128)1)) / (i128)3)); } }  /* n*(n^2+3*n-1)/3 */
static void gen_a085786(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)(((N * ((((i128)2 * (N*N)) + N) + (i128)1)) / (i128)2)); } }  /* n*(2*n^2+n+1)/2 */
static void gen_a090197(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((((N*N*N) + ((i128)6 * (N*N))) + ((i128)6 * N)) + (i128)1)); } }  /* n^3+6*n^2+6*n+1 */
static void gen_a094421(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((N * ((((i128)6 * (N*N)) + ((i128)6 * N)) + (i128)1))); } }  /* n*(6*n^2+6*n+1) */
static void gen_a110451(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N * ((((i128)4 * (N*N)) + ((i128)2 * N)) + (i128)1))); } }  /* n*(4*n^2+2*n+1) */
static void gen_a111396(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N * (N + (i128)7)) * (N + (i128)8)) / (i128)6)); } }  /* n*(n+7)*(n+8)/6 */
static void gen_a125200(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)(((N * ((((i128)4 * (N*N)) + N) - (i128)1)) / (i128)2)); } }  /* n*(4*n^2+n-1)/2 */
static void gen_a008607(u64 *t, long cnt) { gen_res(t, cnt, 25, 1ULL, 0); }
static void gen_a008854(u64 *t, long cnt) { gen_res(t, cnt, 5, 19ULL, 0); }
static void gen_a029739(u64 *t, long cnt) { gen_res(t, cnt, 6, 26ULL, 1); }
static void gen_a032769(u64 *t, long cnt) { gen_res(t, cnt, 5, 23ULL, 0); }
static void gen_a032775(u64 *t, long cnt) { gen_res(t, cnt, 7, 111ULL, 0); }
static void gen_a032793(u64 *t, long cnt) { gen_res(t, cnt, 5, 22ULL, 1); }
static void gen_a032796(u64 *t, long cnt) { gen_res(t, cnt, 7, 110ULL, 1); }
static void gen_a044102(u64 *t, long cnt) { gen_res(t, cnt, 36, 1ULL, 0); }
static void gen_a047202(u64 *t, long cnt) { gen_res(t, cnt, 5, 28ULL, 2); }
static void gen_a047204(u64 *t, long cnt) { gen_res(t, cnt, 5, 24ULL, 3); }
static void gen_a003628(u64 *t, long cnt) { gen_pres(t, cnt, 8, 160ULL); }
static void gen_a033212(u64 *t, long cnt) { gen_pres(t, cnt, 30, 524290ULL); }
static void gen_a039949(u64 *t, long cnt) { gen_pres(t, cnt, 30, 131072ULL); }
static void gen_a042987(u64 *t, long cnt) { gen_pres(t, cnt, 8, 172ULL); }
static void gen_a042988(u64 *t, long cnt) { gen_pres(t, cnt, 7, 63ULL); }
static void gen_a042989(u64 *t, long cnt) { gen_pres(t, cnt, 7, 61ULL); }
static void gen_a042990(u64 *t, long cnt) { gen_pres(t, cnt, 7, 111ULL); }
static void gen_a042992(u64 *t, long cnt) { gen_pres(t, cnt, 7, 109ULL); }
static void gen_a042994(u64 *t, long cnt) { gen_pres(t, cnt, 7, 47ULL); }
static void gen_a042995(u64 *t, long cnt) { gen_pres(t, cnt, 7, 45ULL); }
static void gen_a042997(u64 *t, long cnt) { gen_pres(t, cnt, 7, 124ULL); }
static void gen_a042998(u64 *t, long cnt) { gen_pres(t, cnt, 8, 46ULL); }
static void gen_a045320(u64 *t, long cnt) { gen_pres(t, cnt, 7, 95ULL); }
static void gen_a045322(u64 *t, long cnt) { gen_pres(t, cnt, 7, 93ULL); }
static void gen_a045323(u64 *t, long cnt) { gen_pres(t, cnt, 8, 142ULL); }

/* ---------------------------- the two hundred added 26 Sep 2026 (1001-1200) */
/* One hundred base-dependent sequences (OEIS keyword "base") and one hundred classic
   ones (keywords "core" and "nice"), each checked against its OEIS terms.            */

static u64 read_dec_base(u64 v, u64 b) { u64 r = 0, p = 1; while (v) { r += (v % b) * p; p *= 10; v /= b; } return r; }
static int dsum_b(u64 n, u64 b) { int s = 0; while (n) { s += (int)(n % b); n /= b; } return s; }
static u64 dprod10(u64 n) { u64 p = 1; do { p *= n % 10; n /= 10; } while (n); return p; }
static int all_digits_in(u64 n, unsigned mask) { do { if (!((mask >> (n % 10)) & 1)) return 0; n /= 10; } while (n); return 1; }
static int any_digit_in(u64 n, unsigned mask) { do { if ((mask >> (n % 10)) & 1) return 1; n /= 10; } while (n); return 0; }
static u64 pow10u(int k) { u64 p = 1; while (k--) p *= 10; return p; }

static void gen_sqb(u64 *t, long cnt, u64 b, u64 n_first)          /* n^2 written in base b */
{ for (long i = 0; i < cnt; i++) { u64 n = n_first + (u64)i; t[i] = read_dec_base(n * n, b); } }
static void gen_primeb(u64 *t, long cnt, u64 b)                    /* primes written in base b */
{ gen_primes(t, cnt); for (long i = 0; i < cnt; i++) t[i] = read_dec_base(t[i], b); }
static void gen_onedig(u64 *t, long cnt, u64 b, u64 d)             /* a single digit d in base b */
{ long k = 0; for (u64 n = 0; k < cnt; n++) if (count_dig(n, b, d) == 1) t[k++] = n; }
static void gen_digset(u64 *t, long cnt, unsigned mask, int mode, u64 from)  /* 0: all digits in mask, 1: some digit in mask, 2: no digit in mask */
{
    long k = 0;
    for (u64 n = from; k < cnt; n++) {
        int ok = mode == 0 ? all_digits_in(n, mask) : mode == 1 ? any_digit_in(n, mask) : !any_digit_in(n, mask);
        if (ok) t[k++] = n;
    }
}
static void gen_digadd_b(u64 *t, long cnt, u64 a0, u64 b)          /* a(n+1) = a(n) + digit sum of a(n) in base b */
{ u64 a = a0; for (long i = 0; i < cnt; i++) { t[i] = a; a += (u64)dsum_b(a, b); } }
static void gen_selfb(u64 *t, long cnt, u64 b)                     /* not of the form m + (digit sum of m in base b) */
{
    for (u64 lim = 1ULL << 20; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 m = 0; m <= lim; m++) { u64 g = m + (u64)dsum_b(m, b); if (g <= lim) h[g] = 1; }
        for (u64 n = 1; n <= lim && k < cnt; n++) if (!h[n]) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static void gen_concat_k(u64 *t, long cnt, u64 kmul, u64 add, u64 n_first)   /* n followed by kmul*n + add */
{ for (long i = 0; i < cnt; i++) { u64 n = n_first + (u64)i, m = kmul * n + add; t[i] = n * pow10u(ndigits10(m)) + m; } }
static int ra_steps(u64 n, int cap)                                /* reverse-and-add steps to a palindrome (at least one) */
{
    for (int s = 1; s <= cap; s++) { u64 r = reverse10(n); if (n > (1ULL << 62) - r) return -1; n += r; if (is_palin10(n)) return s; }
    return -1;
}
static void gen_rasteps(u64 *t, long cnt, int want)
{ long k = 0; for (u64 n = 1; k < cnt; n++) if (ra_steps(n, want) == want) t[k++] = n; }
static void gen_strobo(u64 *t, long cnt)                           /* A000787: the same upside down (0 1 8, 6 <-> 9) */
{
    static const int dg[5] = {0, 1, 6, 8, 9}, fl[5] = {0, 1, 9, 8, 6};
    long k = 0;
    for (int len = 1; k < cnt; len++) {
        int h = len / 2, mid = len & 1;
        long nh = 1; for (int i = 0; i < h; i++) nh *= 5;
        /* enumerate in increasing order: first half digits (leading nonzero unless len 1), then middle */
        for (long x = 0; x < nh && k < cnt; x++) {
            int d[16]; long y = x; for (int i = h - 1; i >= 0; i--) { d[i] = (int)(y % 5); y /= 5; }
            if (len > 1 && d[0] == 0) continue;
            for (int m = 0; m < (mid ? 3 : 1) && k < cnt; m++) {
                static const int md[3] = {0, 1, 8};
                if (len == 1) { t[k++] = (u64)md[m]; continue; }
                u64 v = 0;
                for (int i = 0; i < h; i++) v = v * 10 + (u64)dg[d[i]];
                if (mid) v = v * 10 + (u64)md[m];
                for (int i = h - 1; i >= 0; i--) v = v * 10 + (u64)fl[d[i]];
                t[k++] = v;
            }
        }
    }
}
static void gen_trailz_gaps(u64 *t, long cnt)                      /* A000966: n! never ends in this many 0's */
{
    long k = 0; u64 z = 0;
    for (u64 m = 5; k < cnt; m += 5) {
        u64 x = m, add = 0; while (x % 5 == 0) { x /= 5; add++; }
        for (u64 j = 1; j < add && k < cnt; j++) t[k++] = z + j;
        z += add;
    }
}
static void gen_dblbit(u64 *t, long cnt)                           /* A001196: every bit doubled */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i, v = 0, p = 1; while (n) { v += (n & 1) * 3 * p; p *= 4; n >>= 1; } t[i] = v; } }
static void gen_bin0pos(u64 *t, long cnt)                          /* A003607: positions of 0 in 1 10 11 100 ... */
{
    long k = 0; u64 pos = 1;                                       /* "0" comes first, at position 0 */
    t[k++] = 0;
    for (u64 n = 1; k < cnt; n++) { int L = bitlen(n); for (int i = L - 1; i >= 0 && k < cnt; i--) { if (!((n >> i) & 1)) t[k++] = pos; pos++; } }
}
static void gen_selfharshad(u64 *t, long cnt)                      /* A003219 */
{
    for (u64 lim = 1ULL << 22; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 m = 1; m <= lim; m++) { u64 g = m + (u64)digitsum(m); if (g <= lim) h[g] = 1; }
        for (u64 n = 1; n <= lim && k < cnt; n++) if (!h[n] && n % (u64)digitsum(n) == 0) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static int is_semiprime_small(u64 n) { return n >= 4 && omega_big(n) == 2; }
static int div_each_digit(u64 n) { u64 m = n; while (m) { u64 d = m % 10; if (d && n % d) return 0; m /= 10; } return 1; }
N_FILTER(gen_a001101, 1, n % (u64)digitsum(n) == 0 && is_prime(n / (u64)digitsum(n)))
N_FILTER(gen_a002796, 1, div_each_digit(n))
N_FILTER(gen_a085370, 1, n % 3 && n % (u64)digitsum(n) == 0)
N_FILTER(gen_a085371, 1, n % 3 == 0 && n % (u64)digitsum(n) != 0)
N_FILTER(gen_a070938, 1, n % (u64)digitsum(n) == 0 && n % pow10u(ndigits10((u64)digitsum(n))) == (u64)digitsum(n))
N_FILTER(gen_a038367, 1, dprod10(n) % (u64)digitsum(n) == 0)
N_FILTER(gen_a038368, 1, dprod10(n) != (u64)digitsum(n) && n % (dprod10(n) > (u64)digitsum(n) ? dprod10(n) - (u64)digitsum(n) : (u64)digitsum(n) - dprod10(n)) == 0)
N_FILTER(gen_a028838, 1, !(digitsum(n) & (digitsum(n) - 1)))
N_FILTER(gen_a059094, 1, digitsum(n) == 1 || digitsum(n) == 8 || digitsum(n) == 27 || digitsum(n) == 64 || digitsum(n) == 125)
N_FILTER(gen_a085802, 1, is_semiprime_small((u64)digitsum(n)))
N_FILTER(gen_a062713, 1, is_prime((u64)digitsum(n)) && n % (u64)digitsum(n) == 0)
N_FILTER(gen_a062996, 1, (u64)digitsum(n) >= dprod10(n))
N_FILTER(gen_a062997, 1, (u64)digitsum(n) > dprod10(n))
N_FILTER(gen_a062998, 1, (u64)digitsum(n) <= dprod10(n))
N_FILTER(gen_a006364, 0, !(popcount_u64(n >> 1) & 1))

/* classic sequences */
static int odd_exps_only(u64 n) { while (n > 1) { u64 p = spf[n]; int e = 0; while (n % p == 0) { n /= p; e++; } if (!(e & 1)) return 0; } return 1; }
static int all_1mod4(u64 n) { if (n < 2) return 0; while (n > 1) { u64 p = spf[n]; if (p % 4 != 1) return 0; n /= p; } return 1; }
static int refactorable(u64 n) { return n % tau_small(n) == 0; }
static int cop_sigma(u64 n) { return gcd_u64(n, sigma_small(n)) == 1; }
static int aliquot_prime(u64 n) { return n > 1 && is_prime(sigma_small(n) - n); }
static int odd_distinct(u64 n) { return omega_distinct(n) & 1; }
SPF_FILTER(gen_a000028, 1, !bitexp_even(n))
SPF_FILTER(gen_a000415, 2, twosq_small(n) && !is_square(n))
SPF_FILTER(gen_a000430, 2, spf[n] == n || (is_square(n) && spf[isqrt_u64(n)] == isqrt_u64(n)))
SPF_FILTER(gen_a002035, 2, odd_exps_only(n))
SPF_FILTER(gen_a005238, 1, tau_small(n) == tau_small(n + 1) && tau_small(n) == tau_small(n + 2))
SPF_FILTER(gen_a008846, 5, all_1mod4(n))
SPF_FILTER(gen_a014567, 1, cop_sigma(n))
SPF_FILTER(gen_a030230, 2, odd_distinct(n))
SPF_FILTER(gen_a030231, 1, !odd_distinct(n))
SPF_FILTER(gen_a033950, 1, refactorable(n))
SPF_FILTER(gen_a037020, 2, aliquot_prime(n))
SPF_FILTER(gen_a038509, 25, spf[n] != n && (n % 6 == 1 || n % 6 == 5))
SPF_FILTER(gen_a039956, 2, n % 2 == 0 && squarefree_small(n))
SPF_FILTER(gen_a039957, 3, n % 4 == 3 && squarefree_small(n))
SPF_FILTER(gen_a025583, 9, (n & 1) && spf[n] != n && spf[n - 2] != n - 2)
static int norm_phi_form(u64 n) {                                  /* x^2 + xy - y^2: primes 2, 3 mod 5 to even powers */
    u64 m = n; while (m > 1) { u64 p = spf[m]; int e = 0; while (m % p == 0) { m /= p; e++; } if ((p % 5 == 2 || p % 5 == 3) && (e & 1)) return 0; } return 1; }
SPF_FILTER(gen_a031363, 1, norm_phi_form(n))
N_FILTER(gen_a007921, 1, (n & 1) && !is_prime(n + 2))
static void gen_polyprime(u64 *t, long cnt, long long A, long long B, long long C, int out_value, u64 n_first)  /* A n^2 + B n + C prime */
{ long k = 0; for (u64 n = n_first; k < cnt; n++) { i128 v = (i128)A * n * n + (i128)B * n + C; if (isprime_s(v)) t[k++] = out_value ? (u64)v : n; } }
static void gen_a002731(u64 *t, long cnt) { long k = 0; for (u64 n = 1; k < cnt; n += 2) if (isprime_s(((i128)n * n + 1) / 2)) t[k++] = n; }
static void gen_a002815(u64 *t, long cnt)                          /* n + sum_{k <= n} pi(k) */
{ u64 s = 0, pi = 0; for (long i = 0; i < cnt; i++) { u64 n = (u64)i; if (is_prime(n)) pi++; s += pi; t[i] = n + s; } }
static void gen_a002821(u64 *t, long cnt)                          /* nearest integer to n^(3/2) */
{ for (long i = 0; i < cnt; i++) { u64 n = (u64)i, x = n * n * n, r = isqrt_u64(x); t[i] = (x > r * r + r) ? r + 1 : r; } }
static void gen_a002984(u64 *t, long cnt) { u64 a = 1; for (long i = 0; i < cnt; i++) { t[i] = a; a += isqrt_u64(a); } }
static void gen_sum3cubes(u64 *t, long cnt, int mode)              /* x^3 + y^3 + z^3, x <= y <= z >= 1; mode 0 any, 1 exactly one way, 2 more than one way */
{
    for (u64 lim = 1ULL << 20; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 x = 1; 3 * x * x * x <= lim; x++)
            for (u64 y = x; x * x * x + 2 * y * y * y <= lim; y++)
                for (u64 z = y; x * x * x + y * y * y + z * z * z <= lim; z++) { u64 v = x * x * x + y * y * y + z * z * z; if (h[v] < 2) h[v]++; }
        for (u64 n = 1; n <= lim && k < cnt; n++) if (mode == 0 ? h[n] > 0 : mode == 1 ? h[n] == 1 : h[n] == 2) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static void gen_a003072(u64 *t, long cnt) { gen_sum3cubes(t, cnt, 0); }
static void gen_a025395(u64 *t, long cnt) { gen_sum3cubes(t, cnt, 1); }
static void gen_a008917(u64 *t, long cnt) { gen_sum3cubes(t, cnt, 2); }
static u64 wyA(u64 n) { return (n + isqrt_u64(5 * n * n)) / 2; }
static void gen_a003623(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1); t[i] = wyA(wyA(n) + n); } }   /* floor(floor(n phi^2) phi) */
static void gen_a004201(u64 *t, long cnt)                          /* accept m, reject m, m = 1, 2, 3, ... */
{ long k = 0; u64 n = 1; for (u64 m = 1; k < cnt; m++) { for (u64 j = 0; j < m && k < cnt; j++) t[k++] = n++; n += m; } }
static void gen_distsq(u64 *t, long cnt, int terms, int zero_ok)   /* sums of `terms` distinct squares (3 or 4); zero allowed or not */
{
    for (u64 lim = 1ULL << 18; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0; u64 s0 = zero_ok ? 0 : 1;
        for (u64 a = s0; a * a <= lim; a++) for (u64 b = a + 1; a * a + b * b <= lim; b++) for (u64 c = b + 1; a * a + b * b + c * c <= lim; c++) {
            if (terms == 3) { h[a * a + b * b + c * c] = 1; continue; }
            for (u64 d = c + 1; a * a + b * b + c * c + d * d <= lim; d++) h[a * a + b * b + c * c + d * d] = 1;
        }
        for (u64 n = 0; n <= lim && k < cnt; n++) if (h[n]) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static void gen_a001974(u64 *t, long cnt) { gen_distsq(t, cnt, 3, 1); }
static void gen_a004432(u64 *t, long cnt) { gen_distsq(t, cnt, 3, 0); }
static void gen_a004433(u64 *t, long cnt) { gen_distsq(t, cnt, 4, 0); }
static void gen_a001983(u64 *t, long cnt)                          /* x^2 + y^2, 0 <= x < y */
{
    for (u64 lim = 1ULL << 18; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 x = 0; 2 * x * x < lim; x++) for (u64 y = x + 1; x * x + y * y <= lim; y++) h[x * x + y * y] = 1;
        for (u64 n = 0; n <= lim && k < cnt; n++) if (h[n]) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static void gen_a004999(u64 *t, long cnt)                          /* x^3 + y^3, x, y >= 0 */
{
    for (u64 lim = 1ULL << 22; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 x = 0; 2 * x * x * x <= lim; x++) for (u64 y = x; x * x * x + y * y * y <= lim; y++) h[x * x * x + y * y * y] = 1;
        for (u64 n = 0; n <= lim && k < cnt; n++) if (h[n]) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static void gen_a045980(u64 *t, long cnt)                          /* x^3 + y^3 or x^3 - y^3 */
{
    for (u64 lim = 1ULL << 22; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 x = 0; 2 * x * x * x <= lim; x++) for (u64 y = x; x * x * x + y * y * y <= lim; y++) h[x * x * x + y * y * y] = 1;
        for (u64 y = 0; 3 * y * y + 3 * y + 1 <= lim; y++) for (u64 x = y + 1; x * x * x - y * y * y <= lim; x++) h[x * x * x - y * y * y] = 1;
        for (u64 n = 0; n <= lim && k < cnt; n++) if (h[n]) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static void gen_a005658(u64 *t, long cnt)                          /* 1, and 2n, 3n + 2, 6n + 3 */
{
    for (u64 lim = 1ULL << 20; ; lim *= 2) {
        u8 *h = calloc(lim + 1, 1); long k = 0; h[1] = 1;
        for (u64 m = 1; m <= lim; m++) if (h[m]) { if (2 * m <= lim) h[2 * m] = 1; if (3 * m + 2 <= lim) h[3 * m + 2] = 1; if (6 * m + 3 <= lim) h[6 * m + 3] = 1; }
        for (u64 n = 1; n <= lim && k < cnt; n++) if (h[n]) t[k++] = n;
        free(h); if (k == cnt) return;
    }
}
static int p_mid6(u64 p, const u8 *c) { return p > 6 && pr(p - 6, c) && pr(p + 6, c); }
static void gen_a006489(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 22, p_mid6); }
static u64 rnd_phi(u64 n) { return (2 * n + isqrt_u64(5 * 4 * n * n) + 2) / 4; }   /* nearest integer to n phi */
static void gen_a007064(u64 *t, long cnt)                          /* not the nearest integer to n phi */
{ long k = 0; u64 j = 1, r = rnd_phi(1); for (u64 n = 1; k < cnt; n++) { while (r < n) r = rnd_phi(++j); if (r != n) t[k++] = n; } }
static void gen_a007066(u64 *t, long cnt)                          /* 1 + ceil((n-1) phi^2) */
{ for (long i = 0; i < cnt; i++) { u64 m = (u64)i; t[i] = m ? 1 + (wyA(m) + m) + 1 : 1; } }
static void gen_a007491(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { u64 n = (u64)(i + 1), v = n * n + 1; while (!is_prime(v)) v++; t[i] = v; } }   /* from a(1) = 2 */
static void gen_a013916(u64 *t, long cnt)                          /* k with p(1) + ... + p(k) prime */
{
    for (long L = 40L * cnt; ; L *= 2) {
        u64 *p = malloc((size_t)L * sizeof(u64)); gen_primes(p, L);
        long k = 0; u64 s = 0;
        for (long i = 0; i < L && k < cnt; i++) { s += p[i]; if (is_prime(s)) t[k++] = (u64)(i + 1); }
        free(p); if (k == cnt) return;
    }
}
static void gen_twinavg(u64 *t, long cnt) { u64 *q = malloc((size_t)cnt * sizeof(u64)); gen_twinlo(q, cnt); for (long i = 0; i < cnt; i++) t[i] = q[i] + 1; free(q); }
static int divides_2k1(u64 m)                                      /* m | 2^k + 1 for some k >= 0 */
{
    if (m == 1 || m == 2) return 1;
    if (!(m & 1)) return 0;
    need_spf(m, "divides_2k1");
    u64 o = phi_small(m), x = o;                                   /* multiplicative order of 2 mod m */
    while (x > 1) { u64 q = spf[x]; while (x % q == 0) x /= q; while (o % q == 0 && powmod(2, o / q, m) == 1) o /= q; }
    return !(o & 1) && powmod(2, o / 2, m) == m - 1;
}
N_FILTER(gen_a014657, 1, divides_2k1(n))
N_FILTER(gen_a014661, 1, n == 2 || !divides_2k1(n))      /* k > 0 here, so 2 is in */
static void gen_a022797(u64 *t, long cnt)                          /* p(n) + n-th nonprime */
{
    u64 *p = malloc((size_t)cnt * sizeof(u64)); gen_primes(p, cnt);
    long k = 0; for (u64 m = 1; k < cnt; m++) { need_spf(m, "a022797"); if (m == 1 || spf[m] != m) { t[k] = p[k] + m; k++; } }
    free(p);
}
static void gen_a025475(u64 *t, long cnt)                          /* 1 and p^m, m >= 2 */
{
    for (u64 lim = 1ULL << 36; ; lim *= 4) {
        u64 r = isqrt_u64(lim); u8 *c = composite_flags(r + 1);
        size_t cap = 2 * (size_t)r + 1000, m = 0; u64 *v = malloc(cap * sizeof(u64)); v[m++] = 1;
        for (u64 p = 2; p <= r; p++) if (!c[p]) for (u64 x = p * p; ; x *= p) { v[m++] = x; if (x > lim / p) break; }
        free(c); qsort(v, m, sizeof(u64), cmp_u64);
        if ((long)m >= cnt) { memcpy(t, v, (size_t)cnt * sizeof(u64)); free(v); return; }
        free(v);
    }
}
static void gen_a026430(u64 *t, long cnt)                          /* partial sums of 1 + Thue-Morse */
{ u64 s = 0; for (long i = 0; i < cnt; i++) { t[i] = s; s += 1 + (u64)(popcount_u64((u64)i) & 1); } }
static void gen_a027862(u64 *t, long cnt) { gen_polyprime(t, cnt, 2, 2, 1, 1, 0); }
static void gen_a034707(u64 *t, long cnt)                          /* sums of consecutive primes */
{
    for (u64 lim = 1ULL << 20; ; lim *= 2) {
        u8 *c = composite_flags(lim); u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 p = 2; p <= lim; p++) if (!c[p]) { u64 s = 0; for (u64 q = p; q <= lim; q++) if (!c[q]) { s += q; if (s > lim) break; h[s] = 1; } }
        for (u64 n = 1; n <= lim && k < cnt; n++) if (h[n]) t[k++] = n;
        free(c); free(h); if (k == cnt) return;
    }
}
static void gen_a035106(u64 *t, long cnt)                          /* 1, and k(k+1), k(k+2) */
{ long k = 0; t[k++] = 1; for (u64 m = 1; k < cnt; m++) { t[k++] = m * (m + 1); if (k < cnt) t[k++] = m * (m + 2); } }
static void gen_a038550(u64 *t, long cnt)                          /* an odd prime times a power of 2 */
{
    for (u64 lim = 1ULL << 22; ; lim *= 2) {
        u8 *c = composite_flags(lim); long k = 0;
        for (u64 n = 3; n <= lim && k < cnt; n++) { u64 o = n; while (!(o & 1)) o >>= 1; if (o > 2 && !c[o]) t[k++] = n; }
        free(c); if (k == cnt) return;
    }
}
static void gen_ppq(u64 *t, long cnt, int e2)                      /* p^2 + q^e2 over primes p, q */
{
    for (u64 lim = 1ULL << 22; ; lim *= 2) {
        u64 r = isqrt_u64(lim) + 2; u8 *c = composite_flags(r); u8 *h = calloc(lim + 1, 1); long k = 0;
        for (u64 p = 2; p * p <= lim; p++) if (!c[p]) for (u64 q = 2; ; q++) { if (c[q]) continue; u64 qe = e2 == 2 ? q * q : q * q * q; if (p * p + qe > lim) break; h[p * p + qe] = 1; }
        for (u64 n = 1; n <= lim && k < cnt; n++) if (h[n]) t[k++] = n;
        free(c); free(h); if (k == cnt) return;
    }
}
static void gen_a045636(u64 *t, long cnt) { gen_ppq(t, cnt, 2); }
static void gen_a045699(u64 *t, long cnt) { gen_ppq(t, cnt, 3); }
static void gen_a001463(u64 *t, long cnt)                          /* partial sums of Golomb's sequence */
{
    u32 *g = malloc((size_t)(cnt + 2) * sizeof(u32)); g[1] = 1;
    for (long n = 2; n <= cnt; n++) g[n] = 1 + g[n - g[g[n - 1]]];
    u64 s = 0; for (long i = 0; i < cnt; i++) { s += g[i + 1]; t[i] = s; }
    free(g);
}
static void gen_a001768(u64 *t, long cnt)                          /* sum_{k=1..n} ceil(log2(3k/4)) */
{ u64 s = 0; for (long i = 0; i < cnt; i++) { u64 k = (u64)(i + 1), c = 0, v = 1; while (4 * v < 3 * k) { v <<= 1; c++; } s += c; t[i] = s; } }
static void gen_a005598(u64 *t, long cnt)                          /* 1 + sum_{i<=n} (n-i+1) phi(i): partial sums of partial sums */
{ u64 s1 = 0, s2 = 0; for (long i = 0; i < cnt; i++) { u64 n = (u64)i; if (n) { need_spf(n, "a005598"); s1 += phi_small(n); s2 += s1; } t[i] = 1 + s2; } }
static void gen_a003278(u64 *t, long cnt) { gen_cantor(t, cnt); for (long i = 0; i < cnt; i++) t[i] += 1; }
static void gen_a006285(u64 *t, long cnt)                          /* odd numbers not p + 2^k */
{
    long k = 0;
    for (u64 n = 1; k < cnt; n += 2) { int ok = 1; for (u64 q = 1; q < n; q <<= 1) if (is_prime(n - q)) { ok = 0; break; } if (ok) t[k++] = n; }
}
static int x_pow_k_is_2(u64 p, u64 kk)                             /* x^kk = 2 solvable mod p */
{ if (p == 2) return 1; u64 g = gcd_u64(kk, p - 1); return powmod(2, (p - 1) / g, p) == 1; }
static int p_quart2(u64 p, const u8 *c) { (void)c; return x_pow_k_is_2(p, 4); }
static int p_oct2(u64 p, const u8 *c)   { (void)c; return x_pow_k_is_2(p, 8); }
static void gen_a040098(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 22, p_quart2); }
static void gen_a045315(u64 *t, long cnt) { gen_prime_filter(t, cnt, 1ULL << 22, p_oct2); }
static void gen_a002407(u64 *t, long cnt) { gen_polyprime(t, cnt, 3, 3, 1, 1, 1); }   /* cuban primes 3k^2 + 3k + 1 */

static void gen_a001840(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N % 3 == 0) ? (((i128)3 * N + (i128)1 * N * N) / (i128)6) : (N % 3 == 1) ? (((i128)2 + (i128)3 * N + (i128)1 * N * N) / (i128)6) : (((i128)2 + (i128)3 * N + (i128)1 * N * N) / (i128)6)); } }
static void gen_a001859(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N % 2 == 0) ? (((i128)4 * N + (i128)3 * N * N) / (i128)4) : (((i128)1 + (i128)4 * N + (i128)3 * N * N) / (i128)4)); } }
static void gen_a002492(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)2 * N + (i128)6 * N * N + (i128)4 * N * N * N) / (i128)3); } }
static void gen_a002623(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N % 2 == 0) ? (((i128)24 + (i128)34 * N + (i128)15 * N * N + (i128)2 * N * N * N) / (i128)24) : (((i128)21 + (i128)34 * N + (i128)15 * N * N + (i128)2 * N * N * N) / (i128)24)); } }
static void gen_a002717(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N % 2 == 0) ? (((i128)2 * N + (i128)5 * N * N + (i128)2 * N * N * N) / (i128)8) : (((i128)-1 + (i128)2 * N + (i128)5 * N * N + (i128)2 * N * N * N) / (i128)8)); } }
static void gen_a003600(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { if (i < 1) { t[i] = 1; continue; } i128 N = (i128)i + 0; t[i] = (u64)(((i128)8 * N + (i128)3 * N * N + (i128)1 * N * N * N) / (i128)6); } }
static void gen_a004006(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)5 * N + (i128)1 * N * N * N) / (i128)6); } }
static void gen_a005286(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)6 + (i128)20 * N + (i128)9 * N * N + (i128)1 * N * N * N) / (i128)6); } }
static void gen_a005744(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N % 2 == 0) ? (((i128)10 * N + (i128)15 * N * N + (i128)2 * N * N * N) / (i128)24) : (((i128)-3 + (i128)10 * N + (i128)15 * N * N + (i128)2 * N * N * N) / (i128)24)); } }
static void gen_a005893(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { if (i < 1) { t[i] = 1; continue; } i128 N = (i128)i + 0; t[i] = (u64)(((i128)2 + (i128)2 * N * N) / (i128)1); } }
static void gen_a005897(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { if (i < 1) { t[i] = 1; continue; } i128 N = (i128)i + 0; t[i] = (u64)(((i128)2 + (i128)6 * N * N) / (i128)1); } }
static void gen_a005899(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { if (i < 1) { t[i] = 1; continue; } i128 N = (i128)i + 0; t[i] = (u64)(((i128)2 + (i128)4 * N * N) / (i128)1); } }
static void gen_a005901(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { if (i < 1) { t[i] = 1; continue; } i128 N = (i128)i + 0; t[i] = (u64)(((i128)2 + (i128)10 * N * N) / (i128)1); } }
static void gen_a005914(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { if (i < 1) { t[i] = 1; continue; } i128 N = (i128)i + 0; t[i] = (u64)(((i128)2 + (i128)12 * N * N) / (i128)1); } }
static void gen_a005920(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)2 + (i128)6 * N + (i128)7 * N * N + (i128)3 * N * N * N) / (i128)2); } }
static void gen_a005993(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N % 2 == 0) ? (((i128)12 + (i128)14 * N + (i128)6 * N * N + (i128)1 * N * N * N) / (i128)12) : (((i128)6 + (i128)11 * N + (i128)6 * N * N + (i128)1 * N * N * N) / (i128)12)); } }
static void gen_a006004(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)(((i128)-2 + (i128)5 * N + (i128)-2 * N * N + (i128)1 * N * N * N) / (i128)2); } }
static void gen_a006918(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N % 2 == 0) ? (((i128)8 * N + (i128)6 * N * N + (i128)1 * N * N * N) / (i128)24) : (((i128)6 + (i128)11 * N + (i128)6 * N * N + (i128)1 * N * N * N) / (i128)24)); } }
static void gen_a007202(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N % 2 == 0) ? (((i128)4 + (i128)14 * N + (i128)21 * N * N + (i128)14 * N * N * N) / (i128)4) : (((i128)3 + (i128)14 * N + (i128)21 * N * N + (i128)14 * N * N * N) / (i128)4)); } }
static void gen_a007586(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)-2 * N + (i128)1 * N * N + (i128)3 * N * N * N) / (i128)2); } }
static void gen_a007587(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)-7 * N + (i128)3 * N * N + (i128)10 * N * N * N) / (i128)6); } }
static void gen_a007904(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N % 2 == 0) ? (((i128)12 + (i128)26 * N + (i128)15 * N * N + (i128)10 * N * N * N) / (i128)12) : (((i128)9 + (i128)26 * N + (i128)15 * N * N + (i128)10 * N * N * N) / (i128)12)); } }
static void gen_a008412(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { if (i < 1) { t[i] = 1; continue; } i128 N = (i128)i + 0; t[i] = (u64)(((i128)16 * N + (i128)8 * N * N * N) / (i128)3); } }
static void gen_a008577(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N % 3 == 0) ? (((i128)3 + (i128)4 * N + (i128)4 * N * N) / (i128)3) : (N % 3 == 1) ? (((i128)4 + (i128)4 * N + (i128)4 * N * N) / (i128)3) : (((i128)3 + (i128)4 * N + (i128)4 * N * N) / (i128)3)); } }
static void gen_a008580(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { if (i < 1) { t[i] = 1; continue; } i128 N = (i128)i + 0; t[i] = (u64)((N % 2 == 0) ? (((i128)-4 + (i128)10 * N + (i128)9 * N * N) / (i128)4) : (((i128)3 + (i128)8 * N + (i128)9 * N * N) / (i128)4)); } }
static void gen_a008804(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N % 4 == 0) ? (((i128)48 + (i128)44 * N + (i128)12 * N * N + (i128)1 * N * N * N) / (i128)48) : (N % 4 == 1) ? (((i128)42 + (i128)41 * N + (i128)12 * N * N + (i128)1 * N * N * N) / (i128)48) : (N % 4 == 2) ? (((i128)48 + (i128)44 * N + (i128)12 * N * N + (i128)1 * N * N * N) / (i128)48) : (((i128)30 + (i128)41 * N + (i128)12 * N * N + (i128)1 * N * N * N) / (i128)48)); } }
static void gen_a008810(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N % 3 == 0) ? (((i128)1 * N * N) / (i128)3) : (N % 3 == 1) ? (((i128)2 + (i128)1 * N * N) / (i128)3) : (((i128)2 + (i128)1 * N * N) / (i128)3)); } }
static void gen_a019298(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((N % 2 == 0) ? (((i128)2 * N + (i128)3 * N * N + (i128)2 * N * N * N) / (i128)8) : (((i128)1 + (i128)2 * N + (i128)3 * N * N + (i128)2 * N * N * N) / (i128)8)); } }
static void gen_a024206(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((N % 2 == 0) ? (((i128)-4 + (i128)2 * N + (i128)1 * N * N) / (i128)4) : (((i128)-3 + (i128)2 * N + (i128)1 * N * N) / (i128)4)); } }
static void gen_a028552(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)3 * N + (i128)1 * N * N) / (i128)1); } }
static void gen_a033586(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)4 * N + (i128)8 * N * N) / (i128)1); } }
static void gen_a035005(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)(((i128)2 * N + (i128)-12 * N * N + (i128)10 * N * N * N) / (i128)3); } }
static void gen_a035006(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)(((i128)-2 * N * N + (i128)2 * N * N * N) / (i128)1); } }
static void gen_a035008(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)8 * N + (i128)8 * N * N) / (i128)1); } }
static void gen_a045944(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)2 * N + (i128)3 * N * N) / (i128)1); } }
static void gen_a045946(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)6 * N + (i128)18 * N * N) / (i128)1); } }
static void gen_a001740(u64 *t, long cnt) { gen_sqb(t, cnt, 5, 0); }
static void gen_a001741(u64 *t, long cnt) { gen_sqb(t, cnt, 6, 0); }
static void gen_a002440(u64 *t, long cnt) { gen_sqb(t, cnt, 7, 1); }
static void gen_a002441(u64 *t, long cnt) { gen_sqb(t, cnt, 8, 1); }
static void gen_a002442(u64 *t, long cnt) { gen_sqb(t, cnt, 9, 1); }
static void gen_a001363(u64 *t, long cnt) { gen_primeb(t, cnt, 3); }
static void gen_a004678(u64 *t, long cnt) { gen_primeb(t, cnt, 4); }
static void gen_a004679(u64 *t, long cnt) { gen_primeb(t, cnt, 5); }
static void gen_a004680(u64 *t, long cnt) { gen_primeb(t, cnt, 6); }
static void gen_a004681(u64 *t, long cnt) { gen_primeb(t, cnt, 7); }
static void gen_a004682(u64 *t, long cnt) { gen_primeb(t, cnt, 8); }
static void gen_a004683(u64 *t, long cnt) { gen_primeb(t, cnt, 9); }
static void gen_a029954(u64 *t, long cnt) { gen_palb(t, cnt, 7); }
static void gen_a029803(u64 *t, long cnt) { gen_palb(t, cnt, 8); }
static void gen_a029955(u64 *t, long cnt) { gen_palb(t, cnt, 9); }
static void gen_a029956(u64 *t, long cnt) { gen_palb(t, cnt, 11); }
static void gen_a029957(u64 *t, long cnt) { gen_palb(t, cnt, 12); }
static void gen_a029958(u64 *t, long cnt) { gen_palb(t, cnt, 13); }
static void gen_a029959(u64 *t, long cnt) { gen_palb(t, cnt, 14); }
static void gen_a029960(u64 *t, long cnt) { gen_palb(t, cnt, 15); }
static void gen_a029730(u64 *t, long cnt) { gen_palb(t, cnt, 16); }
static void gen_a023699(u64 *t, long cnt) { gen_onedig(t, cnt, 3, 2); }
static void gen_a023706(u64 *t, long cnt) { gen_onedig(t, cnt, 4, 0); }
static void gen_a023710(u64 *t, long cnt) { gen_onedig(t, cnt, 4, 1); }
static void gen_a023714(u64 *t, long cnt) { gen_onedig(t, cnt, 4, 2); }
static void gen_a023718(u64 *t, long cnt) { gen_onedig(t, cnt, 4, 3); }
static void gen_a023722(u64 *t, long cnt) { gen_onedig(t, cnt, 5, 0); }
static void gen_a023726(u64 *t, long cnt) { gen_onedig(t, cnt, 5, 1); }
static void gen_a023730(u64 *t, long cnt) { gen_onedig(t, cnt, 5, 2); }
static void gen_a023734(u64 *t, long cnt) { gen_onedig(t, cnt, 5, 3); }
static void gen_a023738(u64 *t, long cnt) { gen_onedig(t, cnt, 5, 4); }
static void gen_a001729(u64 *t, long cnt) { gen_digset(t, cnt, 0xbe, 0, 0); }
static void gen_a001742(u64 *t, long cnt) { gen_digset(t, cnt, 0xae, 0, 0); }
static void gen_a001743(u64 *t, long cnt) { gen_digset(t, cnt, 0x341, 0, 0); }
static void gen_a001744(u64 *t, long cnt) { gen_digset(t, cnt, 0x351, 0, 0); }
static void gen_a001745(u64 *t, long cnt) { gen_digset(t, cnt, 0x351, 1, 0); }
static void gen_a001746(u64 *t, long cnt) { gen_digset(t, cnt, 0x341, 1, 0); }
static void gen_a046030(u64 *t, long cnt) { gen_digset(t, cnt, 0x213, 0, 0); }
static void gen_a046034(u64 *t, long cnt) { gen_digset(t, cnt, 0xac, 0, 0); }
static void gen_a028846(u64 *t, long cnt) { gen_digset(t, cnt, 0x116, 0, 0); }
static void gen_a006507(u64 *t, long cnt) { gen_digadd_b(t, cnt, 7, 10); }
static void gen_a007618(u64 *t, long cnt) { gen_digadd_b(t, cnt, 5, 10); }
static void gen_a010062(u64 *t, long cnt) { gen_digadd_b(t, cnt, 1, 2); }
static void gen_a010063(u64 *t, long cnt) { gen_digadd_b(t, cnt, 1, 3); }
static void gen_a010065(u64 *t, long cnt) { gen_digadd_b(t, cnt, 1, 4); }
static void gen_a010066(u64 *t, long cnt) { gen_digadd_b(t, cnt, 1, 5); }
static void gen_a010068(u64 *t, long cnt) { gen_digadd_b(t, cnt, 1, 6); }
static void gen_a010069(u64 *t, long cnt) { gen_digadd_b(t, cnt, 1, 7); }
static void gen_a010071(u64 *t, long cnt) { gen_digadd_b(t, cnt, 1, 8); }
static void gen_a010072(u64 *t, long cnt) { gen_digadd_b(t, cnt, 1, 9); }
static void gen_a010064(u64 *t, long cnt) { gen_selfb(t, cnt, 4); }
static void gen_a010067(u64 *t, long cnt) { gen_selfb(t, cnt, 6); }
static void gen_a010070(u64 *t, long cnt) { gen_selfb(t, cnt, 8); }
static void gen_a001704(u64 *t, long cnt) { gen_concat_k(t, cnt, 1, 1, 1); }
static void gen_a019550(u64 *t, long cnt) { gen_concat_k(t, cnt, 2, 0, 1); }
static void gen_a019553(u64 *t, long cnt) { gen_concat_k(t, cnt, 5, 0, 1); }
static void gen_a009440(u64 *t, long cnt) { gen_concat_k(t, cnt, 6, 0, 1); }
static void gen_a009441(u64 *t, long cnt) { gen_concat_k(t, cnt, 7, 0, 1); }
static void gen_a009470(u64 *t, long cnt) { gen_concat_k(t, cnt, 8, 0, 1); }
static void gen_a009474(u64 *t, long cnt) { gen_concat_k(t, cnt, 9, 0, 1); }
static void gen_a015976(u64 *t, long cnt) { gen_rasteps(t, cnt, 1); }
static void gen_a015977(u64 *t, long cnt) { gen_rasteps(t, cnt, 2); }
static void gen_a015979(u64 *t, long cnt) { gen_rasteps(t, cnt, 3); }
static void gen_a015980(u64 *t, long cnt) { gen_rasteps(t, cnt, 4); }
static void gen_a015982(u64 *t, long cnt) { gen_rasteps(t, cnt, 5); }
static void gen_a015984(u64 *t, long cnt) { gen_rasteps(t, cnt, 6); }
static void gen_a001912(u64 *t, long cnt) { gen_polyprime(t, cnt, 4, 0, 1, 0, 0); }
static void gen_a002837(u64 *t, long cnt) { gen_polyprime(t, cnt, 1, -1, 41, 0, 0); }
static void gen_a002081(u64 *t, long cnt) { gen_res(t, cnt, 20, 65812ULL, 2); }
static void gen_a016825(u64 *t, long cnt) { gen_res(t, cnt, 4, 4ULL, 2); }
static void gen_a042965(u64 *t, long cnt) { gen_res(t, cnt, 4, 11ULL, 0); }
static void gen_a003626(u64 *t, long cnt) { gen_pres(t, cnt, 20, (1ULL<<11)|(1ULL<<13)|(1ULL<<17)|(1ULL<<19)); }


/* base-dependent, continued */
static int fact_digit_sum(u64 n) { int s = 0; for (u64 b = 2; n; b++) { s += (int)(n % b); n /= b; } return s; }
N_FILTER(gen_a049445, 1, n % (u64)popcount_u64(n) == 0)
N_FILTER(gen_a064481, 1, n % (u64)dsum_b(n, 5) == 0)
N_FILTER(gen_a101813, 1, (n & 1) && n % (u64)digitsum(n) == 0)
N_FILTER(gen_a101814, 1, !(n & 1) && n % (u64)digitsum(n) == 0)
N_FILTER(gen_a118363, 1, n % (u64)fact_digit_sum(n) == 0)
N_FILTER(gen_a077436, 0, n < (1ULL << 26) && popcount_u64(n) == popcount_u64(n * n))
N_FILTER(gen_a038366, 1, n % (dprod10(n) + (u64)digitsum(n)) == 0)
N_FILTER(gen_a028839, 1, is_square((u64)digitsum(n)))
static void gen_a016052(u64 *t, long cnt) { gen_digadd_b(t, cnt, 3, 10); }
static void gen_a016096(u64 *t, long cnt) { gen_digadd_b(t, cnt, 9, 10); }
static void gen_a033298(u64 *t, long cnt) { u64 a = 1; for (long i = 0; i < cnt; i++) { t[i] = a; a += (u64)digitsum(a * a); } }
static void gen_a007612(u64 *t, long cnt) { u64 a = 1; for (long i = 0; i < cnt; i++) { t[i] = a; a += 1 + (a - 1) % 9; } }
static void gen_a019551(u64 *t, long cnt) { gen_concat_k(t, cnt, 3, 0, 1); }
static void gen_a019552(u64 *t, long cnt) { gen_concat_k(t, cnt, 4, 0, 1); }

/* ---------------------------- the four hundred added 26 Sep 2026 (1201-1600) */
/* Polynomial prime conditions read from the OEIS names (numbers n with f(n) prime,
   primes of the form f(k), primes p with f(p) prime), semiprime and squarefree
   variants, and more from the verified pools; each checked against its OEIS terms. */

static int big_omega64(u64 v) { if (v < 2) return 0; factor(v); int c = 0; for (int i = 0; i < g_npf; i++) c += g_pe[i]; return c; }
static int squarefree64(u64 v) { if (v < 1) return 0; if (v == 1) return 1; factor(v); for (int i = 0; i < g_npf; i++) if (g_pe[i] > 1) return 0; return 1; }
static int fits64(i128 v) { return v >= 0 && v < ((i128)1 << 64); }
#define NF_PRIME(name, expr) \
static void name(u64 *t, long cnt) { long k = 0; for (u64 n = 0; k < cnt; n++) { i128 N = (i128)n; i128 v = (expr); if (fits64(v) && isprime_s(v)) t[k++] = n; } }
#define PP_PRIME(name, expr) \
static void name(u64 *t, long cnt) \
{ for (u64 lim = 1ULL << 22; ; lim *= 2) { u8 *c = composite_flags(lim); long k = 0; \
    for (u64 p = 2; p <= lim && k < cnt; p++) { if (c[p]) continue; i128 N = (i128)p; i128 v = (expr); if (fits64(v) && isprime_s(v)) t[k++] = p; } \
    free(c); if (k == cnt) return; } }
/* primes f(k), k >= 0, in increasing order: all prime values up to a bound, sorted */
#define PF_PRIME(name, expr) \
static void name(u64 *t, long cnt) \
{ for (i128 bound = (i128)1 << 30; ; bound *= 4) { size_t cap = 1 << 16, m = 0; u64 *v = malloc(cap * sizeof(u64)); \
    for (u64 n = 0; ; n++) { i128 N = (i128)n; i128 x = (expr); if (x > bound && n > 64) break; \
      if (x >= 2 && x <= bound && isprime_s(x)) { if (m == cap) { cap *= 2; v = realloc(v, cap * sizeof(u64)); } v[m++] = (u64)x; } } \
    qsort(v, m, sizeof(u64), cmp_u64); size_t u = 0; for (size_t i = 0; i < m; i++) if (!u || v[i] != v[u - 1]) v[u++] = v[i]; \
    if ((long)u >= cnt) { memcpy(t, v, (size_t)cnt * sizeof(u64)); free(v); return; } free(v); } }
#define NF_SEMI(name, expr) \
static void name(u64 *t, long cnt) { long k = 0; for (u64 n = 0; k < cnt; n++) { i128 N = (i128)n; i128 v = (expr); if (fits64(v) && big_omega64((u64)v) == 2) t[k++] = n; } }
#define PP_TEST(name, expr, test) \
static void name(u64 *t, long cnt) \
{ for (u64 lim = 1ULL << 22; ; lim *= 2) { u8 *c = composite_flags(lim); long k = 0; \
    for (u64 p = 2; p <= lim && k < cnt; p++) { if (c[p]) continue; i128 N = (i128)p; i128 v = (expr); if (fits64(v) && test((u64)v)) t[k++] = p; } \
    free(c); if (k == cnt) return; } }
static int is_semi64(u64 v) { return big_omega64(v) == 2; }
#define SEMI_FORM(name, expr) \
static void name(u64 *t, long cnt) { long k = 0; for (u64 n = 0; k < cnt; n++) { i128 N = (i128)n; i128 v = (expr); if (fits64(v) && v > 3 && big_omega64((u64)v) == 2) t[k++] = (u64)v; } }

PF_PRIME(gen_a002327, (((N*N) - N) - (i128)1))   /* n^2-n-1 */
NF_PRIME(gen_a002328, (((N*N) - N) - (i128)1))   /* n^2-n-1 */
PF_PRIME(gen_a002383, (((N*N) + N) + (i128)1))   /* n^2+n+1 */
NF_PRIME(gen_a002384, (((N*N) + N) + (i128)1))   /* n^2+n+1 */
NF_PRIME(gen_a002970, (((i128)4 * (N*N)) + (i128)9))   /* 4*n^2+9 */
NF_PRIME(gen_a002971, (((i128)4 * (N*N)) + (i128)25))   /* 4*n^2+25 */
PF_PRIME(gen_a005846, (((N*N) + N) + (i128)41))   /* n^2+n+41 */
PF_PRIME(gen_a007635, (((N*N) + N) + (i128)17))   /* n^2+n+17 */
PF_PRIME(gen_a007637, ((((i128)3 * (N*N)) - ((i128)3 * N)) + (i128)23))   /* 3*n^2-3*n+23 */
PF_PRIME(gen_a007639, ((((i128)2 * (N*N)) - ((i128)2 * N)) + (i128)19))   /* 2*n^2-2*n+19 */
PF_PRIME(gen_a007641, (((i128)2 * (N*N)) + (i128)29))   /* 2*n^2+29 */
PP_PRIME(gen_a023219, (((i128)5 * N) + (i128)6))   /* 5*n+6 */
PF_PRIME(gen_a027753, (((N*N) + N) + (i128)3))   /* n^2+n+3 */
PF_PRIME(gen_a027755, (((N*N) + N) + (i128)5))   /* n^2+n+5 */
PF_PRIME(gen_a027758, (((N*N) + N) + (i128)9))   /* n^2+n+9 */
NF_PRIME(gen_a027861, ((N*N) + ((N + (i128)1)*(N + (i128)1))))   /* n^2+(n+1)^2 */
NF_PRIME(gen_a027863, (((N*N) + ((N + (i128)1)*(N + (i128)1))) + ((N + (i128)2)*(N + (i128)2))))   /* n^2+(n+1)^2+(n+2)^2 */
NF_PRIME(gen_a027866, ((((((N*N) + ((N + (i128)1)*(N + (i128)1))) + ((N + (i128)2)*(N + (i128)2))) + ((N + (i128)3)*(N + (i128)3))) + ((N + (i128)4)*(N + (i128)4))) + ((N + (i128)5)*(N + (i128)5))))   /* n^2+(n+1)^2+(n+2)^2+(n+3)^2+(n+4)^2+(n+5)^2 */
PF_PRIME(gen_a027867, ((((((N*N) + ((N + (i128)1)*(N + (i128)1))) + ((N + (i128)2)*(N + (i128)2))) + ((N + (i128)3)*(N + (i128)3))) + ((N + (i128)4)*(N + (i128)4))) + ((N + (i128)5)*(N + (i128)5))))   /* n^2+(n+1)^2+(n+2)^2+(n+3)^2+(n+4)^2+(n+5)^2 */
PF_PRIME(gen_a037029, (((i128)666 * N) + (i128)1))   /* 666*n+1 */
PF_PRIME(gen_a048059, (((N*N) + N) + (i128)11))   /* n^2+n+11 */
PF_PRIME(gen_a048988, ((((i128)4 * (N*N)) + ((i128)4 * N)) + (i128)59))   /* 4*n^2+4*n+59 */
PF_PRIME(gen_a050265, (((i128)2 * (N*N)) + (i128)11))   /* 2*n^2+11 */
PP_PRIME(gen_a051647, (((i128)210 * N) + (i128)1))   /* 210*n+1 */
PP_PRIME(gen_a051653, (((i128)2310 * N) + (i128)1))   /* 2310*n+1 */
PP_PRIME(gen_a051654, (((i128)30030 * N) + (i128)1))   /* 30030*n+1 */
PP_PRIME(gen_a052291, (((i128)4 * (N*N)) + (i128)1))   /* 4*n^2+1 */
PP_PRIME(gen_a053182, (((N*N) + N) + (i128)1))   /* n^2+n+1 */
PP_PRIME(gen_a053184, (((N*N) + N) - (i128)1))   /* n^2+n-1 */
NF_PRIME(gen_a055494, (((N*N) - N) + (i128)1))   /* n^2-n+1 */
NF_PRIME(gen_a056906, (((i128)36 * (N*N)) + (i128)5))   /* 36*n^2+5 */
NF_PRIME(gen_a056908, ((((i128)36 * (N*N)) + ((i128)36 * N)) + (i128)13))   /* 36*n^2+36*n+13 */
PF_PRIME(gen_a057604, (((i128)4 * (N*N)) + (i128)163))   /* 4*n^2+163 */
PF_PRIME(gen_a059425, (((N*N) + ((i128)19 * N)) + (i128)17))   /* n^2+19*n+17 */
PF_PRIME(gen_a060844, ((((i128)6 * (N*N)) + ((i128)6 * N)) + (i128)31))   /* 6*n^2+6*n+31 */
PF_PRIME(gen_a061242, (((i128)9 * N) - (i128)1))   /* 9*n-1 */
PF_PRIME(gen_a062800, (((i128)100 * N) + (i128)1))   /* 100*n+1 */
PF_PRIME(gen_a063472, (((i128)666 * N) - (i128)1))   /* 666*n-1 */
PP_PRIME(gen_a065508, (((N*N) - N) + (i128)1))   /* n^2-n+1 */
NF_PRIME(gen_a066049, (((i128)2 * (N*N)) - (i128)1))   /* 2*n^2-1 */
PF_PRIME(gen_a066436, (((i128)2 * (N*N)) - (i128)1))   /* 2*n^2-1 */
PF_PRIME(gen_a073102, (((i128)210 * N) + (i128)1))   /* 210*n+1 */
PF_PRIME(gen_a076339, (((i128)512 * N) + (i128)1))   /* 512*n+1 */
PF_PRIME(gen_a076727, ((N*N) + ((N + (i128)3)*(N + (i128)3))))   /* n^2+(n+3)^2 */
NF_PRIME(gen_a083022, (((i128)4 * (N*N)) - (i128)3))   /* 4*n^2-3 */
NF_PRIME(gen_a086285, (((i128)1 + ((i128)2 * N)) + ((i128)3 * (N*N))))   /* 1+2*n+3*n^2 */
NF_PRIME(gen_a086298, (((i128)1 - ((i128)2 * N)) + ((i128)3 * (N*N))))   /* 1-2*n+3*n^2 */
NF_PRIME(gen_a086303, (N + (i128)15))   /* n+15 */
NF_PRIME(gen_a086304, (N + (i128)6))   /* n+6 */
NF_PRIME(gen_a088572, (((((i128)2 * N) + (i128)1)*(((i128)2 * N) + (i128)1)) - (i128)2))   /* (2*n+1)^2-2 */
NF_PRIME(gen_a088758, (((((i128)4 * N) + (i128)1)*(((i128)4 * N) + (i128)1)) + ((((i128)4 * N) + (i128)2)*(((i128)4 * N) + (i128)2))))   /* (4*n+1)^2+(4*n+2)^2 */
NF_PRIME(gen_a088759, (((((i128)4 * N) + (i128)3)*(((i128)4 * N) + (i128)3)) + ((((i128)4 * N) + (i128)2)*(((i128)4 * N) + (i128)2))))   /* (4*n+3)^2+(4*n+2)^2 */
PF_PRIME(gen_a088955, (((i128)60 * N) + (i128)1))   /* 60*n+1 */
NF_PRIME(gen_a088967, (N + (i128)9))   /* n+9 */
NF_PRIME(gen_a089001, (((i128)2 * (N*N)) + (i128)1))   /* 2*n^2+1 */
NF_PRIME(gen_a089008, (((i128)18 * (N*N)) + (i128)1))   /* 18*n^2+1 */
NF_PRIME(gen_a089063, (((i128)840 * N) + (i128)175177943))   /* 840*n+175177943 */
NF_PRIME(gen_a089373, (((N*N) - ((i128)7 * N)) + (i128)7))   /* n^2-7*n+7 */
PF_PRIME(gen_a089376, (((N*N) - ((i128)7 * N)) + (i128)7))   /* n^2-7*n+7 */
PP_PRIME(gen_a089438, (((i128)6 * N) + (i128)11))   /* 6*n+11 */
PP_PRIME(gen_a089441, (((i128)16 * N) + (i128)17))   /* 16*n+17 */
NF_PRIME(gen_a089593, (((N*N) + ((i128)2 * N)) + (i128)2))   /* n^2+2*n+2 */
NF_PRIME(gen_a089623, (((N*N) + ((i128)2 * N)) - (i128)1))   /* n^2+2*n-1 */
NF_PRIME(gen_a089681, (((i128)3 * (N*N)) - (i128)1))   /* 3*n^2-1 */
PF_PRIME(gen_a089682, (((i128)3 * (N*N)) - (i128)1))   /* 3*n^2-1 */
NF_PRIME(gen_a089747, (((N*N) - ((i128)2 * N)) + (i128)5))   /* n^2-2*n+5 */
PF_PRIME(gen_a090187, (((i128)11 * N) + (i128)2))   /* 11*n+2 */
PF_PRIME(gen_a090562, ((((i128)5 * (N*N)) + ((i128)5 * N)) + (i128)1))   /* 5*n^2+5*n+1 */
NF_PRIME(gen_a090563, ((((i128)5 * (N*N)) + ((i128)5 * N)) + (i128)1))   /* 5*n^2+5*n+1 */
PF_PRIME(gen_a090684, (((i128)8 * (N*N)) - (i128)1))   /* 8*n^2-1 */
PF_PRIME(gen_a090685, (((i128)8 * (N*N)) + (i128)1))   /* 8*n^2+1 */
PF_PRIME(gen_a090686, (((i128)6 * (N*N)) - (i128)1))   /* 6*n^2-1 */
PF_PRIME(gen_a090687, (((i128)6 * (N*N)) + (i128)1))   /* 6*n^2+1 */
NF_PRIME(gen_a090696, ((N*N) - (i128)11))   /* n^2-11 */
PF_PRIME(gen_a090698, (((i128)2 * (N*N)) + (i128)1))   /* 2*n^2+1 */
NF_PRIME(gen_a091271, (((i128)4 * (N*N)) - (i128)11))   /* 4*n^2-11 */
PP_PRIME(gen_a091567, (((N*N) - N) - (i128)1))   /* n^2-n-1 */
NF_PRIME(gen_a092968, (((i128)2 * (N*N)) + (i128)11))   /* 2*n^2+11 */
PF_PRIME(gen_a093359, (((i128)28 * N) + (i128)1))   /* 28*n+1 */
PF_PRIME(gen_a093838, (((i128)36 * N) + (i128)1))   /* 36*n+1 */
NF_PRIME(gen_a094210, (((N*N) + ((i128)3 * N)) + (i128)1))   /* n^2+3*n+1 */
PF_PRIME(gen_a094407, (((i128)16 * N) + (i128)1))   /* 16*n+1 */
PF_PRIME(gen_a095995, (((i128)100 * N) - (i128)1))   /* 100*n-1 */
NF_PRIME(gen_a096689, ((((i128)2 * (N*N)) + ((i128)3 * N)) + (i128)3))   /* 2*n^2+3*n+3 */
NF_PRIME(gen_a096691, ((((i128)8 * (N*N)) + ((i128)6 * N)) + (i128)3))   /* 8*n^2+6*n+3 */
PF_PRIME(gen_a098828, ((((i128)2 * (N*N)) + ((i128)2 * N)) - (i128)1))   /* 2*n^2+2*n-1 */
PF_PRIME(gen_a099007, ((((i128)6 * (N*N)) - ((i128)2 * N)) - (i128)1))   /* 6*n^2-2*n-1 */
PF_PRIME(gen_a100201, (((i128)23 * N) + (i128)3))   /* 23*n+3 */
PF_PRIME(gen_a100202, (((i128)13 * N) + (i128)3))   /* 13*n+3 */
PF_PRIME(gen_a100203, (((i128)37 * N) + (i128)3))   /* 37*n+3 */
PF_PRIME(gen_a100494, (((i128)47 * N) + (i128)3))   /* 47*n+3 */
PF_PRIME(gen_a100760, (((i128)47 * N) + (i128)5))   /* 47*n+5 */
NF_PRIME(gen_a101444, (((i128)9973 * N) + (i128)10007))   /* (9973*n+10007) */
NF_PRIME(gen_a101567, (((i128)1009 * N) + (i128)10007))   /* 1009*n+10007 */
PF_PRIME(gen_a101780, (((i128)100 * N) + (i128)3))   /* 100*n+3 */
PF_PRIME(gen_a102130, ((((i128)8 * (N*N)) + ((i128)4 * N)) + (i128)1))   /* 8*n^2+4*n+1 */
NF_PRIME(gen_a102166, ((((i128)2 * (N*N)) + ((i128)11 * N)) + (i128)101))   /* 2*n^2+11*n+101 */
NF_PRIME(gen_a102339, ((N * ((i128)10*(i128)10*(i128)10)) + (i128)333))   /* n*10^3+333 */
NF_PRIME(gen_a102343, ((N * ((i128)10*(i128)10*(i128)10)) + (i128)777))   /* n*10^3+777 */
NF_PRIME(gen_a102649, ((((i128)11 * (N*N)) + ((i128)11 * N)) + (i128)3))   /* 11*n^2+11*n+3 */
NF_PRIME(gen_a102657, ((((i128)11 * (N*N)) + ((i128)11 * N)) + (i128)1))   /* 11*n^2+11*n+1 */
PF_PRIME(gen_a102732, (((i128)13 * N) + (i128)5))   /* 13*n+5 */
PF_PRIME(gen_a102734, (((i128)23 * N) + (i128)5))   /* 23*n+5 */
PF_PRIME(gen_a102851, (((i128)19 * N) + (i128)5))   /* 19*n+5 */
PP_PRIME(gen_a103564, (((i128)3 * (N*N)) + (i128)2))   /* 3*n^2+2 */
PP_PRIME(gen_a103776, ((((i128)8 * (N*N)) + ((i128)4 * N)) + (i128)1))   /* 8*n^2+4*n+1 */
NF_PRIME(gen_a105057, (((i128)10000 * N) - (i128)1))   /* 10000*n-1 */
NF_PRIME(gen_a105059, (((i128)100000 * N) - (i128)1))   /* 100000*n-1 */
NF_PRIME(gen_a105107, (((i128)10000 * N) + (i128)1001))   /* 10000*n+1001 */
PF_PRIME(gen_a105126, (((i128)16 * N) + (i128)9))   /* 16*n+9 */
PF_PRIME(gen_a105127, (((i128)32 * N) + (i128)17))   /* 32*n+17 */
PF_PRIME(gen_a105128, (((i128)64 * N) + (i128)33))   /* 64*n+33 */
PF_PRIME(gen_a105129, (((i128)128 * N) + (i128)65))   /* 128*n+65 */
PF_PRIME(gen_a105130, (((i128)256 * N) + (i128)129))   /* 256*n+129 */
PF_PRIME(gen_a105131, (((i128)512 * N) + (i128)257))   /* 512*n+257 */
PF_PRIME(gen_a105132, (((i128)1024 * N) + (i128)513))   /* 1024*n+513 */
NF_PRIME(gen_a105140, (((i128)1024 * N) + (i128)513))   /* 1024*n+513 */
NF_PRIME(gen_a105680, (((i128)1009 * N) + (i128)9973))   /* 1009*n+9973 */
NF_PRIME(gen_a105710, (((i128)9973 * N) + (i128)1009))   /* 9973*n+1009 */
PF_PRIME(gen_a105854, (((i128)20 * N) + (i128)3))   /* 20*n+3 */
PP_PRIME(gen_a106483, (((i128)2 * (N*N)) - (i128)1))   /* 2*n^2-1 */
NF_PRIME(gen_a106699, (((i128)1009 * N) - (i128)9973))   /* 1009*n-9973 */
NF_PRIME(gen_a106700, (((i128)9973 * N) - (i128)1009))   /* 9973*n-1009 */
PF_PRIME(gen_a107003, (((i128)24 * N) + (i128)5))   /* 24*n+5 */
NF_PRIME(gen_a107071, (((i128)1019 * N) + (i128)1021))   /* 1019*n+1021 */
NF_PRIME(gen_a107072, (((i128)1021 * N) + (i128)1019))   /* 1021*n+1019 */
NF_PRIME(gen_a107301, (((i128)10007 * N) + (i128)99991))   /* 10007*n+99991 */
NF_PRIME(gen_a107302, (((i128)99991 * N) + (i128)10007))   /* 99991*n+10007 */
NF_PRIME(gen_a107303, (((i128)3 * N) - (i128)5))   /* (3*n-5) */
NF_PRIME(gen_a107306, (((i128)17 * N) - (i128)19))   /* (17*n-19) */
NF_PRIME(gen_a107308, (((i128)29 * N) - (i128)31))   /* (29*n-31) */
NF_PRIME(gen_a108058, (((i128)179 * N) + (i128)181))   /* 179*n+181 */
NF_PRIME(gen_a108059, (((i128)181 * N) + (i128)179))   /* 181*n+179 */
NF_PRIME(gen_a108060, (((i128)191 * N) + (i128)193))   /* 191*n+193 */
NF_PRIME(gen_a108061, (((i128)193 * N) + (i128)191))   /* 193*n+191 */
NF_PRIME(gen_a108187, (((i128)11 * N) - (i128)5))   /* 11*n-5 */
NF_PRIME(gen_a108232, (((i128)11 * N) - (i128)7))   /* 11*n-7 */
NF_PRIME(gen_a108233, (((i128)11 * N) + (i128)5))   /* 11*n+5 */
NF_PRIME(gen_a108341, (((i128)997 * N) - (i128)1009))   /* 997*n-1009 */
NF_PRIME(gen_a108342, (((i128)1009 * N) - (i128)997))   /* 1009*n-997 */
NF_PRIME(gen_a108584, (((i128)10 * N) - (i128)97))   /* 10*n-97 */
NF_PRIME(gen_a108588, (((i128)10 * N) + (i128)97))   /* 10*n+97 */
NF_PRIME(gen_a108594, (((i128)10 * N) + (i128)101))   /* 10*n+101 */
NF_PRIME(gen_a108595, (((i128)10 * N) + (i128)103))   /* 10*n+103 */
NF_PRIME(gen_a108596, (((i128)911 * N) - (i128)7))   /* 911*n-7 */
NF_PRIME(gen_a108597, (((i128)911 * N) - (i128)11))   /* 911*n-11 */
NF_PRIME(gen_a108724, (((i128)11 * N) + (i128)17))   /* 11*n+17 */
NF_PRIME(gen_a108725, (((i128)11 * N) + (i128)19))   /* 11*n+19 */
NF_PRIME(gen_a108726, (((i128)11 * N) + (i128)29))   /* 11*n+29 */
NF_PRIME(gen_a108727, (((i128)11 * N) + (i128)31))   /* 11*n+31 */
NF_PRIME(gen_a108751, (((i128)11 * N) - (i128)911))   /* 11*n-911 */
NF_PRIME(gen_a108757, (((i128)1000 * N) + (i128)911))   /* 1000*n+911 */
NF_PRIME(gen_a108762, (((i128)911 * N) + (i128)13))   /* 911*n+13 */
NF_PRIME(gen_a108854, (((i128)10 * N) - (i128)127))   /* 10*n-127 */
NF_PRIME(gen_a108855, (((i128)10 * N) + (i128)127))   /* 10*n+127 */
NF_PRIME(gen_a108856, (((i128)10 * N) - (i128)131))   /* 10*n-131 */
NF_PRIME(gen_a108857, (((i128)10 * N) + (i128)131))   /* 10*n+131 */
NF_PRIME(gen_a108874, (((i128)41 * N) + (i128)43))   /* 41*n+43 */
NF_PRIME(gen_a108899, (((i128)11 * N) + (i128)2357))   /* 11*n+2357 */
NF_PRIME(gen_a108900, (((i128)2357 * N) + (i128)11))   /* 2357*n+11 */
NF_PRIME(gen_a108901, (((i128)2357 * N) + (i128)23))   /* 2357*n+23 */
NF_PRIME(gen_a108902, (((i128)23 * N) + (i128)2357))   /* 23*n+2357 */
NF_PRIME(gen_a108936, (((i128)11 * N) + (i128)911))   /* 11*n+911 */
NF_PRIME(gen_a108937, (((i128)911 * N) + (i128)11))   /* 911*n+11 */
NF_PRIME(gen_a108938, (((i128)911 * N) + (i128)7))   /* 911*n+7 */
NF_PRIME(gen_a108969, (((i128)43 * N) + (i128)41))   /* 43*n+41 */
NF_PRIME(gen_a108976, (((i128)17 * N) + (i128)19))   /* 17*n+19 */
NF_PRIME(gen_a108977, (((i128)19 * N) + (i128)17))   /* 19*n+17 */
NF_PRIME(gen_a108978, (((i128)29 * N) + (i128)31))   /* 29*n+31 */
NF_PRIME(gen_a108979, (((i128)31 * N) + (i128)29))   /* 31*n+29 */
NF_PRIME(gen_a109603, (((i128)43 * N) - (i128)41))   /* 43*n-41 */
NF_PRIME(gen_a109604, (((i128)41 * N) - (i128)43))   /* 41*n-43 */
NF_PRIME(gen_a109605, (((i128)100000 * N) + (i128)91111))   /* 100000*n+91111 */
NF_PRIME(gen_a110801, (((i128)12 * N) + (i128)1))   /* 12*n+1 */
NF_PRIME(gen_a110913, (((i128)23 * (N*N)) - (i128)49))   /* 23*n^2-49 */
NF_PRIME(gen_a110959, (((i128)23 * (N*N)) + (i128)1))   /* 23*n^2+1 */
NF_PRIME(gen_a110960, (((i128)23 * (N*N)) + (i128)4))   /* 23*n^2+4 */
NF_PRIME(gen_a110961, (((i128)23 * (N*N)) + (i128)9))   /* 23*n^2+9 */
NF_PRIME(gen_a110964, (((i128)23 * (N*N)) + (i128)16))   /* 23*n^2+16 */
NF_PRIME(gen_a110965, (((i128)23 * (N*N)) + (i128)25))   /* 23*n^2+25 */
NF_PRIME(gen_a110966, (((i128)23 * (N*N)) + (i128)36))   /* 23*n^2+36 */
NF_PRIME(gen_a110967, (((i128)23 * (N*N)) + (i128)49))   /* 23*n^2+49 */
NF_PRIME(gen_a110974, (((i128)23 * (N*N)) - (i128)1))   /* 23*n^2-1 */
NF_PRIME(gen_a110994, (((i128)23 * (N*N)) - (i128)4))   /* 23*n^2-4 */
NF_PRIME(gen_a110998, (((i128)23 * (N*N)) - (i128)9))   /* 23*n^2-9 */
NF_PRIME(gen_a110999, (((i128)23 * (N*N)) - (i128)16))   /* 23*n^2-16 */
NF_PRIME(gen_a111001, (((i128)23 * (N*N)) - (i128)25))   /* 23*n^2-25 */
NF_PRIME(gen_a111040, (((i128)2 * (N*N)) + (i128)9))   /* 2*n^2+9 */
NF_PRIME(gen_a111041, (((i128)2 * (N*N)) + (i128)25))   /* 2*n^2+25 */
NF_PRIME(gen_a111051, (((i128)3 * (N*N)) + (i128)1))   /* 3*n^2+1 */
NF_PRIME(gen_a111052, (((i128)3 * (N*N)) + (i128)4))   /* 3*n^2+4 */
NF_PRIME(gen_a111068, (((i128)3 * (N*N)) + (i128)16))   /* 3*n^2+16 */
NF_PRIME(gen_a111069, (((i128)3 * (N*N)) + (i128)25))   /* 3*n^2+25 */
NF_PRIME(gen_a111082, (((i128)3 * (N*N)) + (i128)49))   /* 3*n^2+49 */
NF_PRIME(gen_a111083, (((i128)3 * (N*N)) + (i128)64))   /* 3*n^2+64 */
NF_PRIME(gen_a111094, (((i128)18 * N) + (i128)1))   /* 18*n+1 */
NF_PRIME(gen_a111147, (((i128)5 * (N*N)) + (i128)1))   /* 5*n^2+1 */
NF_PRIME(gen_a111148, (((i128)5 * (N*N)) + (i128)4))   /* 5*n^2+4 */
NF_PRIME(gen_a111149, (((i128)5 * (N*N)) + (i128)9))   /* 5*n^2+9 */
NF_PRIME(gen_a111174, (((i128)24 * N) + (i128)1))   /* 24*n+1 */
NF_PRIME(gen_a111175, (((i128)30 * N) + (i128)1))   /* 30*n+1 */
NF_PRIME(gen_a111251, ((((i128)3 * (N*N)) + ((i128)3 * N)) + (i128)1))   /* 3*n^2+3*n+1 */
NF_PRIME(gen_a111292, ((((i128)6 * (N*N)) + ((i128)6 * N)) + (i128)1))   /* 6*n^2+6*n+1 */
NF_PRIME(gen_a111294, (((i128)23 * N) + (i128)2))   /* 23*n+2 */
NF_PRIME(gen_a111312, (((i128)11 * N) + (i128)2))   /* 11*n+2 */
NF_PRIME(gen_a111369, (((i128)13 * N) + (i128)11))   /* 13*n+11 */
NF_PRIME(gen_a111455, (((i128)101 * N) + (i128)97))   /* 101*n+97 */
PP_PRIME(gen_a113151, (((i128)19 * N) + (i128)2))   /* 19*n+2 */
NF_PRIME(gen_a113487, (((i128)17 * N) + (i128)2))   /* 17*n+2 */
NF_PRIME(gen_a113488, (((i128)19 * N) + (i128)2))   /* 19*n+2 */
NF_PRIME(gen_a113510, (((i128)29 * N) + (i128)2))   /* 29*n+2 */
PF_PRIME(gen_a117047, (((i128)60 * N) + (i128)11))   /* 60*n+11 */
PF_PRIME(gen_a117049, (((i128)22 * (N*N)) + (i128)1))   /* 22*(n^2)+1 */
NF_PRIME(gen_a119409, (((i128)235 * N) + (i128)1))   /* 235*n+1 */
NF_PRIME(gen_a120344, (((i128)23 * N) + (i128)1))   /* 23*n+1 */
NF_PRIME(gen_a120345, (((i128)2357 * N) + (i128)1))   /* 2357*n+1 */
NF_PRIME(gen_a121068, (((i128)8 * (N*N)) + (i128)7))   /* 8*n^2+7 */
NF_PRIME(gen_a121817, ((i128)23 + (((i128)36 * N) * (N + (i128)1))))   /* 23+36*n*(n+1) */
PF_PRIME(gen_a122114, ((((i128)2 * (N*N)) + ((i128)26 * N)) + (i128)1))   /* 2*n^2+26*n+1 */
PF_PRIME(gen_a122430, (((i128)1 + ((i128)2 * N)) + ((i128)3 * (N*N))))   /* 1+2*n+3*n^2 */
PP_PRIME(gen_a122482, (((i128)1 + ((i128)4 * N)) + ((i128)12 * (N*N))))   /* 1+4*n+12*n^2 */
NF_PRIME(gen_a124127, (((i128)17 * N) + (i128)1))   /* 17*n+1 */
NF_PRIME(gen_a124198, (((i128)21 * N) + (i128)1))   /* 21*n+1 */
NF_PRIME(gen_a124204, (((i128)20 * N) + (i128)1))   /* 20*n+1 */
NF_PRIME(gen_a126332, (((i128)10 * N) + (i128)13))   /* 10*n+13 */
NF_PRIME(gen_a126785, (((i128)10 * N) + (i128)11))   /* 10*n+11 */
PP_PRIME(gen_a126960, ((((i128)3 * N)*((i128)3 * N)) + (i128)2))   /* (3*n)^2+2 */
PP_PRIME(gen_a127435, (((N - (i128)1)*(N - (i128)1)) + (i128)1))   /* (n-1)^2+1 */
NF_PRIME(gen_a127575, (((i128)16 * N) + (i128)15))   /* 16*n+15 */
PF_PRIME(gen_a127576, (((i128)16 * N) + (i128)15))   /* 16*n+15 */
PF_PRIME(gen_a127579, (((i128)64 * N) + (i128)63))   /* 64*n+63 */
NF_PRIME(gen_a127580, (((i128)64 * N) + (i128)63))   /* 64*n+63 */
PF_PRIME(gen_a127589, (((i128)16 * N) + (i128)5))   /* 16*n+5 */
NF_PRIME(gen_a127590, (((i128)16 * N) + (i128)5))   /* 16*n+5 */
NF_PRIME(gen_a127591, (((i128)64 * N) + (i128)21))   /* 64*n+21 */
PF_PRIME(gen_a127592, (((i128)64 * N) + (i128)21))   /* 64*n+21 */
PF_PRIME(gen_a127593, (((i128)256 * N) + (i128)85))   /* 256*n+85 */
NF_PRIME(gen_a127594, (((i128)256 * N) + (i128)85))   /* 256*n+85 */
NF_PRIME(gen_a128829, (((i128)6 * (N*N)) + (i128)17))   /* 6*n^2+17 */
PF_PRIME(gen_a129484, (((i128)17 * N) + (i128)1))   /* 17*n+1 */
PP_TEST(gen_a039787, (N - (i128)1), squarefree64)   /* n-1 */
PP_TEST(gen_a049097, (N + (i128)1), squarefree64)   /* n+1 */
PP_TEST(gen_a049231, (N - (i128)2), squarefree64)   /* n-2 */
PP_TEST(gen_a049233, (N + (i128)2), squarefree64)   /* n+2 */
NF_SEMI(gen_a085722, ((N*N) + (i128)1))   /* n^2+1 */
NF_SEMI(gen_a085746, (((N*N) + N) + (i128)1))   /* n^2+n+1 */
PP_TEST(gen_a092109, (N + (i128)3), is_semi64)   /* n+3 */
SEMI_FORM(gen_a108181, (((i128)4 * N) + (i128)1))   /* 4*n+1 */
NF_SEMI(gen_a108769, ((N*N) + ((N + (i128)1)*(N + (i128)1))))   /* n^2+(n+1)^2 */
PP_TEST(gen_a109953, ((N*N) + (i128)2), is_semi64)   /* n^2+2 */
SEMI_FORM(gen_a112771, (((i128)6 * N) + (i128)1))   /* 6*n+1 */
SEMI_FORM(gen_a112772, (((i128)6 * N) + (i128)2))   /* 6*n+2 */
SEMI_FORM(gen_a112774, (((i128)6 * N) + (i128)4))   /* 6*n+4 */
NF_SEMI(gen_a112775, (((i128)6 * N) + (i128)1))   /* 6*n+1 */
NF_SEMI(gen_a112776, (((i128)6 * N) + (i128)5))   /* 6*n+5 */
NF_SEMI(gen_a112777, (((i128)2 * (N*N)) + (i128)1))   /* 2*n^2+1 */
NF_SEMI(gen_a122488, (((i128)1 + ((i128)2 * N)) + ((i128)3 * (N*N))))   /* 1+2*n+3*n^2 */
static void gen_a033253(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 83, 1); }
static void gen_a033254(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 85, 1); }
static void gen_a033255(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 86, 1); }
static void gen_a033256(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 87, 1); }
static void gen_a033257(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 89, 1); }
static void gen_a033258(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 91, 1); }
static void gen_a106880(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 9, 1); }
static void gen_a106890(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 11, 1); }
static void gen_a106900(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 12, 1); }
static void gen_a106901(u64 *t, long cnt) { gen_bform(t, cnt, 3, 3, 5, 1); }
static void gen_a106903(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 13, 1); }
static void gen_a106905(u64 *t, long cnt) { gen_bform(t, cnt, 2, 2, 7, 1); }
static void gen_a106907(u64 *t, long cnt) { gen_bform(t, cnt, 4, 3, 4, 1); }
static void gen_a106921(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 15, 1); }
static void gen_a106926(u64 *t, long cnt) { gen_bform(t, cnt, 2, 1, 8, 1); }
static void gen_a106929(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 16, 1); }
static void gen_a106931(u64 *t, long cnt) { gen_bform(t, cnt, 4, 4, 5, 1); }
static void gen_a106932(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 17, 1); }
static void gen_a106934(u64 *t, long cnt) { gen_bform(t, cnt, 3, 2, 6, 1); }
static void gen_a106937(u64 *t, long cnt) { gen_bform(t, cnt, 2, 2, 9, 1); }
static void gen_a106939(u64 *t, long cnt) { gen_bform(t, cnt, 4, 3, 5, 1); }
static void gen_a106945(u64 *t, long cnt) { gen_bform(t, cnt, 2, 1, 9, 1); }
static void gen_a106949(u64 *t, long cnt) { gen_bform(t, cnt, 2, 0, 9, 1); }
static void gen_a106950(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 18, 1); }
static void gen_a106951(u64 *t, long cnt) { gen_bform(t, cnt, 3, 3, 7, 1); }
static void gen_a106953(u64 *t, long cnt) { gen_bform(t, cnt, 4, 2, 5, 1); }
static void gen_a106959(u64 *t, long cnt) { gen_bform(t, cnt, 2, 1, 10, 1); }
static void gen_a106964(u64 *t, long cnt) { gen_bform(t, cnt, 3, 2, 7, 1); }
static void gen_a106966(u64 *t, long cnt) { gen_bform(t, cnt, 3, 1, 7, 1); }
static void gen_a106969(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 21, 1); }
static void gen_a106971(u64 *t, long cnt) { gen_bform(t, cnt, 5, 4, 5, 1); }
static void gen_a106974(u64 *t, long cnt) { gen_bform(t, cnt, 2, 2, 11, 1); }
static void gen_a106975(u64 *t, long cnt) { gen_bform(t, cnt, 4, 3, 6, 1); }
static void gen_a106978(u64 *t, long cnt) { gen_bform(t, cnt, 3, 3, 8, 1); }
static void gen_a106980(u64 *t, long cnt) { gen_bform(t, cnt, 2, 1, 11, 1); }
static void gen_a106984(u64 *t, long cnt) { gen_bform(t, cnt, 2, 0, 11, 1); }
static void gen_a106985(u64 *t, long cnt) { gen_bform(t, cnt, 5, 3, 5, 1); }
static void gen_a106988(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 23, 1); }
static void gen_a106990(u64 *t, long cnt) { gen_bform(t, cnt, 5, 5, 6, 1); }
static void gen_a106992(u64 *t, long cnt) { gen_bform(t, cnt, 4, 1, 6, 1); }
static void gen_a106995(u64 *t, long cnt) { gen_bform(t, cnt, 3, 1, 8, 1); }
static void gen_a106998(u64 *t, long cnt) { gen_bform(t, cnt, 2, 1, 12, 1); }
static void gen_a107002(u64 *t, long cnt) { gen_bform(t, cnt, 5, 2, 5, 1); }
static void gen_a107005(u64 *t, long cnt) { gen_bform(t, cnt, 4, 4, 7, 1); }
static void gen_a107007(u64 *t, long cnt) { gen_bform(t, cnt, 3, 0, 8, 1); }
static void gen_a107008(u64 *t, long cnt) { gen_bform(t, cnt, 1, 0, 24, 1); }
static void gen_a107009(u64 *t, long cnt) { gen_bform(t, cnt, 5, 1, 5, 1); }
static void gen_a107012(u64 *t, long cnt) { gen_bform(t, cnt, 1, 1, 25, 1); }
static void gen_a107132(u64 *t, long cnt) { gen_bform(t, cnt, 2, 0, 13, 1); }
static void gen_a107133(u64 *t, long cnt) { gen_bform(t, cnt, 4, 0, 7, 1); }
static void gen_a127736(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)(((N * (((N*N) + ((i128)2 * N)) - (i128)1)) / (i128)2)); } }  /* n*(n^2+2*n-1)/2 */
static void gen_a006000(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N + (i128)1) * (((N*N) + N) + (i128)2)) / (i128)2)); } }  /* (n+1)*(n^2+n+2)/2 */
static void gen_a016061(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N * (N + (i128)1)) * (((i128)4 * N) + (i128)5)) / (i128)6)); } }  /* n*(n+1)*(4*n+5)/6 */
static void gen_a027620(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N + ((N + (i128)1)*(N + (i128)1))) + ((N + (i128)2)*(N + (i128)2)*(N + (i128)2)))); } }  /* n+(n+1)^2+(n+2)^3 */
static void gen_a033994(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((((N * (N + (i128)1)) * (((i128)5 * N) + (i128)1)) / (i128)6)); } }  /* n*(n+1)*(5*n+1)/6 */
static void gen_a035328(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (((i128)2 * N) - (i128)1)) * (((i128)2 * N) + (i128)1))); } }  /* n*(2*n-1)*(2*n+1) */
static void gen_a035329(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (((i128)2 * N) + (i128)5)) * (((i128)2 * N) + (i128)7))); } }  /* n*(2*n+5)*(2*n+7) */
static void gen_a037235(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * ((((i128)2 * (N*N)) - ((i128)3 * N)) + (i128)4)) / (i128)3)); } }  /* n*(2*n^2-3*n+4)/3 */
static void gen_a056578(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((((i128)1 + ((i128)2 * N)) + ((i128)3 * (N*N))) + ((i128)4 * (N*N*N)))); } }  /* 1+2*n+3*n^2+4*n^3 */
static void gen_a071230(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * ((((i128)6 * (N*N)) - ((i128)7 * N)) + (i128)3)) / (i128)2)); } }  /* n*(6*n^2-7*n+3)/2 */
static void gen_a086605(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((((i128)9 * (N*N*N)) - ((i128)18 * (N*N))) + ((i128)10 * N))); } }  /* 9*n^3-18*n^2+10*n */
static void gen_a101853(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)(((N * (((i128)20 + ((i128)15 * N)) + (N*N))) / (i128)6)); } }  /* n*(20+15*n+n^2)/6 */
static void gen_a102094(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)(((((i128)2 * N) - (i128)1) * ((((i128)2 * N) + (i128)1)*(((i128)2 * N) + (i128)1)))); } }  /* (2*n-1)*(2*n+1)^2 */
static void gen_a111144(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N * (N + (i128)13)) * (N + (i128)14)) / (i128)6)); } }  /* n*(n+13)*(n+14)/6 */
static void gen_a115519(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * (((i128)1 + ((i128)3 * N)) + ((i128)6 * (N*N)))) / (i128)2)); } }  /* n*(1+3*n+6*n^2)/2 */
static void gen_a126335(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)(((N * ((((i128)4 * (N*N)) + ((i128)5 * N)) - (i128)3)) / (i128)2)); } }  /* n*(4*n^2+5*n-3)/2 */
static void gen_a074742(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((((N*N*N) + ((i128)6 * (N*N))) - N) + (i128)12) / (i128)6)); } }  /* (n^3+6*n^2-n+12)/6 */
static void gen_a100207(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((((i128)4 + ((i128)8 * N)) + ((i128)10 * (N*N))) + ((i128)4 * (N*N*N)))); } }  /* 4+8*n+10*n^2+4*n^3 */
static void gen_a000297(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + -1; t[i] = (u64)(((((N + (i128)1) * (N + (i128)3)) * (N + (i128)8)) / (i128)6)); } }  /* (n+1)*(n+3)*(n+8)/6 */
static void gen_a027602(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N*N*N) + ((N + (i128)1)*(N + (i128)1)*(N + (i128)1))) + ((N + (i128)2)*(N + (i128)2)*(N + (i128)2)))); } }  /* n^3+(n+1)^3+(n+2)^3 */
static void gen_a027849(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N + (i128)1) * ((((i128)5 * (N*N)) + ((i128)4 * N)) + (i128)1))); } }  /* (n+1)*(5*n^2+4*n+1) */
static void gen_a049480(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((((((i128)2 * N) - (i128)1) * (((N*N) - N) + (i128)6)) / (i128)6)); } }  /* (2*n-1)*(n^2-n+6)/6 */
static void gen_a056520(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N + (i128)2) * ((((i128)2 * (N*N)) - N) + (i128)3)) / (i128)6)); } }  /* (n+2)*(2*n^2-n+3)/6 */
static void gen_a063488(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((((((i128)2 * N) - (i128)1) * (((N*N) - N) + (i128)2)) / (i128)2)); } }  /* (2*n-1)*(n^2-n+2)/2 */
static void gen_a101165(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((((i128)7 * (N*N*N)) + ((i128)6 * (N*N))) + ((i128)5 * N)) / (i128)6)); } }  /* (7*n^3+6*n^2+5*n)/6 */
static void gen_a071229(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((N * ((((i128)14 * (N*N)) - ((i128)21 * N)) + (i128)13)) / (i128)6)); } }  /* n*(14*n^2-21*n+13)/6 */
static void gen_a101860(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((((i128)3 + N) * (((i128)2 + ((i128)33 * N)) + (N*N))) / (i128)6)); } }  /* (3+n)*(2+33*n+n^2)/6 */
static void gen_a114211(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((((((i128)5 * (N*N*N)) + ((i128)12 * (N*N))) + N) + (i128)6) / (i128)6)); } }  /* (5*n^3+12*n^2+n+6)/6 */
static void gen_a011199(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N + (i128)1) * (((i128)2 * N) + (i128)1)) * (((i128)3 * N) + (i128)1))); } }  /* (n+1)*(2*n+1)*(3*n+1) */
static void gen_a079588(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((N + (i128)1) * (((i128)2 * N) + (i128)1)) * (((i128)4 * N) + (i128)1))); } }  /* (n+1)*(2*n+1)*(4*n+1) */
static void gen_a095796(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((i128)1 + ((((((i128)26 * N) + (i128)17) + ((i128)7 * (N*N))) * N) / (i128)2))); } }  /* 1+(26*n+17+7*n^2)*n/2 */
static void gen_a100504(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((((((i128)4 * (N*N*N)) + ((i128)6 * (N*N))) + ((i128)8 * N)) + (i128)6) / (i128)3)); } }  /* (4*n^3+6*n^2+8*n+6)/3 */
static void gen_a034721(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)(((((((i128)10 * (N*N*N)) - ((i128)9 * (N*N))) + ((i128)2 * N)) / (i128)3) + (i128)1)); } }  /* (10*n^3-9*n^2+2*n)/3+1 */
static void gen_a087863(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((((N*N*N) + ((i128)24 * (N*N))) + ((i128)65 * N)) + (i128)36) / (i128)6)); } }  /* (n^3+24*n^2+65*n+36)/6 */
static void gen_a057813(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((((i128)2 * N) + (i128)1) * ((((i128)4 * (N*N)) + ((i128)4 * N)) + (i128)3)) / (i128)3)); } }  /* (2*n+1)*(4*n^2+4*n+3)/3 */
static void gen_a061550(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 0; t[i] = (u64)((((((i128)2 * N) + (i128)1) * (((i128)2 * N) + (i128)3)) * (((i128)2 * N) + (i128)5))); } }  /* (2*n+1)*(2*n+3)*(2*n+5) */
static void gen_a063489(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((((((i128)2 * N) - (i128)1) * ((((i128)5 * (N*N)) - ((i128)5 * N)) + (i128)6)) / (i128)6)); } }  /* (2*n-1)*(5*n^2-5*n+6)/6 */
static void gen_a063490(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((((((i128)2 * N) - (i128)1) * ((((i128)7 * (N*N)) - ((i128)7 * N)) + (i128)6)) / (i128)6)); } }  /* (2*n-1)*(7*n^2-7*n+6)/6 */
static void gen_a063491(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((((((i128)2 * N) - (i128)1) * ((((i128)3 * (N*N)) - ((i128)3 * N)) + (i128)2)) / (i128)2)); } }  /* (2*n-1)*(3*n^2-3*n+2)/2 */
static void gen_a063494(u64 *t, long cnt) { for (long i = 0; i < cnt; i++) { i128 N = (i128)i + 1; t[i] = (u64)((((((i128)2 * N) - (i128)1) * ((((i128)7 * (N*N)) - ((i128)7 * N)) + (i128)3)) / (i128)3)); } }  /* (2*n-1)*(7*n^2-7*n+3)/3 */
static void gen_a047205(u64 *t, long cnt) { gen_res(t, cnt, 5, 25ULL, 0); }
static void gen_a047206(u64 *t, long cnt) { gen_res(t, cnt, 5, 26ULL, 1); }
static void gen_a047207(u64 *t, long cnt) { gen_res(t, cnt, 5, 27ULL, 0); }
static void gen_a047208(u64 *t, long cnt) { gen_res(t, cnt, 5, 17ULL, 0); }
static void gen_a047218(u64 *t, long cnt) { gen_res(t, cnt, 5, 9ULL, 0); }
static void gen_a047219(u64 *t, long cnt) { gen_res(t, cnt, 5, 10ULL, 1); }
static void gen_a047221(u64 *t, long cnt) { gen_res(t, cnt, 5, 12ULL, 2); }
static void gen_a047222(u64 *t, long cnt) { gen_res(t, cnt, 5, 13ULL, 0); }
static void gen_a047223(u64 *t, long cnt) { gen_res(t, cnt, 5, 14ULL, 1); }
static void gen_a047227(u64 *t, long cnt) { gen_res(t, cnt, 6, 30ULL, 1); }
static void gen_a047228(u64 *t, long cnt) { gen_res(t, cnt, 6, 28ULL, 2); }
static void gen_a047230(u64 *t, long cnt) { gen_res(t, cnt, 6, 24ULL, 3); }
static void gen_a047231(u64 *t, long cnt) { gen_res(t, cnt, 6, 25ULL, 0); }
static void gen_a047233(u64 *t, long cnt) { gen_res(t, cnt, 6, 17ULL, 0); }
static void gen_a047234(u64 *t, long cnt) { gen_res(t, cnt, 6, 19ULL, 0); }
static void gen_a047235(u64 *t, long cnt) { gen_res(t, cnt, 6, 20ULL, 2); }
static void gen_a047236(u64 *t, long cnt) { gen_res(t, cnt, 6, 22ULL, 1); }
static void gen_a047237(u64 *t, long cnt) { gen_res(t, cnt, 6, 23ULL, 0); }
static void gen_a047242(u64 *t, long cnt) { gen_res(t, cnt, 6, 11ULL, 0); }
static void gen_a047243(u64 *t, long cnt) { gen_res(t, cnt, 6, 12ULL, 2); }
static void gen_a047244(u64 *t, long cnt) { gen_res(t, cnt, 6, 13ULL, 0); }
static void gen_a047245(u64 *t, long cnt) { gen_res(t, cnt, 6, 14ULL, 1); }
static void gen_a047249(u64 *t, long cnt) { gen_res(t, cnt, 6, 56ULL, 3); }
static void gen_a047252(u64 *t, long cnt) { gen_res(t, cnt, 6, 59ULL, 0); }
static void gen_a047253(u64 *t, long cnt) { gen_res(t, cnt, 6, 62ULL, 1); }
static void gen_a045324(u64 *t, long cnt) { gen_pres(t, cnt, 7, 31ULL); }
static void gen_a045325(u64 *t, long cnt) { gen_pres(t, cnt, 7, 29ULL); }
static void gen_a045328(u64 *t, long cnt) { gen_pres(t, cnt, 7, 79ULL); }
static void gen_a045329(u64 *t, long cnt) { gen_pres(t, cnt, 7, 77ULL); }
static void gen_a045346(u64 *t, long cnt) { gen_pres(t, cnt, 7, 119ULL); }
static void gen_a045347(u64 *t, long cnt) { gen_pres(t, cnt, 7, 117ULL); }
static void gen_a045350(u64 *t, long cnt) { gen_pres(t, cnt, 7, 55ULL); }
static void gen_a045351(u64 *t, long cnt) { gen_pres(t, cnt, 7, 53ULL); }
static void gen_a045352(u64 *t, long cnt) { gen_pres(t, cnt, 8, 166ULL); }
static void gen_a045353(u64 *t, long cnt) { gen_pres(t, cnt, 7, 103ULL); }
static void gen_a045354(u64 *t, long cnt) { gen_pres(t, cnt, 7, 101ULL); }
static void gen_a045358(u64 *t, long cnt) { gen_pres(t, cnt, 7, 39ULL); }
static void gen_a045369(u64 *t, long cnt) { gen_pres(t, cnt, 7, 87ULL); }
static void gen_a045370(u64 *t, long cnt) { gen_pres(t, cnt, 7, 85ULL); }
static void gen_a045376(u64 *t, long cnt) { gen_pres(t, cnt, 7, 71ULL); }
static void gen_a045393(u64 *t, long cnt) { gen_pres(t, cnt, 7, 123ULL); }
static void gen_a045394(u64 *t, long cnt) { gen_pres(t, cnt, 7, 121ULL); }
static void gen_a045396(u64 *t, long cnt) { gen_pres(t, cnt, 7, 59ULL); }
static void gen_a045397(u64 *t, long cnt) { gen_pres(t, cnt, 7, 57ULL); }
static void gen_a045398(u64 *t, long cnt) { gen_pres(t, cnt, 7, 107ULL); }
static void gen_a047254(u64 *t, long cnt) { gen_res(t, cnt, 6, 44ULL, 2); }
static void gen_a047256(u64 *t, long cnt) { gen_res(t, cnt, 6, 47ULL, 0); }
static void gen_a047257(u64 *t, long cnt) { gen_res(t, cnt, 6, 48ULL, 4); }
static void gen_a047258(u64 *t, long cnt) { gen_res(t, cnt, 6, 49ULL, 0); }
static void gen_a047259(u64 *t, long cnt) { gen_res(t, cnt, 6, 50ULL, 1); }
static void gen_a047260(u64 *t, long cnt) { gen_res(t, cnt, 6, 51ULL, 0); }
static void gen_a047262(u64 *t, long cnt) { gen_res(t, cnt, 6, 53ULL, 0); }
static void gen_a047263(u64 *t, long cnt) { gen_res(t, cnt, 6, 55ULL, 0); }

/* ------------------------------------------------------------------ */
/* catalogue                                                           */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *id, *anum, *name, *family, *note;
    long terms; long n0;
    void (*gen)(u64 *, long);
} SeqDef;

static SeqDef defs[] = {
{ "naturals", "A000027", "Natural numbers", "base case",
  "The base case: d = 1 everywhere, so k = spf(a - 1) and the weight sheet is the sieve of Eratosthenes. The level class is the single line L = 1 (a - 1 prime): 9,592 = pi(99,999) terms here. The ties k = L are a - 1 = p^2: 65 = pi(316).",
  100000, 1, gen_naturals },
{ "nonsquares", "A000037", "Non-squares", "complement",
  "Non-squares: gaps of 1, and 2 across each square. The plane is the naturals' almost unchanged (9.65 % level against 9.59 %); the gap-2 terms add a short line at L = 2.",
  100000, 1, gen_nonsquares },
{ "primes", "A000040", "Prime numbers", "primes",
  "weight A117078, level A117563, jump A001223, l(n) A118534. Every weight and every level is odd. A prime > 3 is the lesser of a twin pair iff k = 3. Below 10^8 only 2, 3 and 7 fail to decompose.",
  100000, 1, gen_primes },
{ "wythoff", "A000201", "Lower Wythoff sequence", "Beatty",
  "floor(n phi), a Beatty sequence. The gaps are 1 and 2 only, in the Fibonacci-word pattern, so k is the least divisor of a - d above 1 or 2, as for the naturals and the odd numbers. The level share, 12.20 %, falls between theirs.",
  100000, 1, gen_wythoff },
{ "triangular", "A000217", "Triangular numbers", "polynomial",
  "d = n + 1 and l = (n + 1)(n - 2)/2 = d(n - 2)/2, so l <= d^2 always and L/k < 1/2: every term is level-classified and the cloud stops a factor of two short of the diagonal.",
  100000, 0, gen_triangular },
{ "squares", "A000290", "Squares", "polynomial",
  "d = 2n + 1 and l = (n - 1)^2 - 2 < d^2, so every decomposable term is forced level (100 %). For a ~ c n^2 in general l/d^2 -> 1/(4c), so every quadratic with c > 1/4 ends up entirely level (A002620, c = 1/4, sits on the edge); the nine polynomial-type sequences here differ only in their ray structure.",
  100000, 0, gen_squares },
{ "tetrahedral", "A000292", "Tetrahedral numbers", "polynomial",
  "d = (n + 1)(n + 2)/2 grows like n^2 while l grows like n^3/6, so l/d^2 -> 0: every decomposable term is forced level (100 %), and the line L = 1 is thin (329 terms).",
  100000, 0, gen_tetrahedral },
{ "pentagonal", "A000326", "Pentagonal numbers", "polynomial",
  "d = 3n + 1 against l = (3n^2 - 7n - 2)/2 < d^2: forced level at every decomposable term (100 %).",
  100000, 0, gen_pentagonal },
{ "pyramidal", "A000330", "Square pyramidal numbers", "polynomial",
  "d = (n + 1)^2 against l ~ n^3/3, so l/d^2 -> 0: forced level at every decomposable term (100 %).",
  100000, 0, gen_pyramidal },
{ "hexagonal", "A000384", "Hexagonal numbers", "polynomial",
  "d = 4n + 1 against l = 2n^2 - 5n - 1 < d^2: forced level at every decomposable term (100 %).",
  100000, 0, gen_hexagonal },
{ "cubes", "A000578", "Cubes", "polynomial",
  "d = 3n^2 + 3n + 1 against l ~ n^3, so l/d^2 -> 0: forced level at every decomposable term (100 %). The line L = 2 is slightly fuller than L = 1 (7,847 against 7,814 terms).",
  100000, 0, gen_cubes },
{ "lucky", "A000959", "Lucky numbers", "sieve",
  "A sieve of the same shape as Eratosthenes, with the same density x / log x. As for the primes, every term is odd and every gap even, so l is odd and every weight is odd. The level share is 32.31 %, against 23.00 % for the primes over a comparable range (a up to 1.43e6 and 1.30e6).",
  100000, 1, gen_lucky },
{ "flavius", "A000960", "Flavius Josephus sieve", "sieve",
  "count(x) = 2 sqrt(x/pi), so a(n) ~ pi n^2/4 and the gap outgrows sqrt(l): l <= d^2 on 77 % of terms, and 96.6 % of the sequence is level-classified. It is the spread of the gap, not its mean, that populates the weight class.",
  100000, 1, gen_flavius },
{ "primepowers", "A000961", "Prime powers", "powers",
  "Prime powers, with 1: the primes plus a sparse set of higher powers. The plane is the primes' plane almost exactly (22.97 % level against 23.00 %).",
  100000, 1, gen_primepowers },
{ "primesq", "A001248", "Squares of primes", "powers",
  "a = p^2 and d = q^2 - p^2 >= 4p + 4, so l < p^2 < d^2: every decomposable term is forced level (100 %). a decomposes iff q^2 < 1.5 p^2. That fails only for p = 2, 3, 5, 7, 13, 23; Nagura (a prime in [x, 6x/5] for x >= 25) rules out any other.",
  100000, 1, gen_primesq },
{ "semiprimes", "A001358", "Semiprimes", "multiplicative",
  "Omega(n) = 2. The gaps stay small (1 to 47 here) and l/d^2 is large, so only 5 terms are forced level. The level share, 15.71 %, lies between the naturals' and the primes'. The line L = 1 (l prime) carries 61 % of the level class.",
  100000, 1, gen_semiprimes },
{ "twinlo", "A001359", "Lesser of twin primes", "primes",
  "Past 5 every term is 5 mod 6, so d = 0 and l = 5 (mod 6). No square is 5 mod 6, so k = L can never happen: zero ties, by proof rather than by count. l is prime to 6, so every weight is odd and not a multiple of 3. The level share is 47.77 %.",
  100000, 1, gen_twinlo },
{ "twosquares", "A001481", "Sums of two squares", "multiplicative",
  "n = x^2 + y^2: every prime 3 mod 4 divides n to an even power. Starts at a(1) = 0. Density ~ 0.764 x / sqrt(log x) (Landau-Ramanujan). The gaps stay below 35 here and the level share is 15.95 %.",
  100000, 1, gen_twosq },
{ "non3", "A001651", "Numbers not divisible by 3", "forced divisor",
  "The gaps alternate 1, 2 and 3 divides l at every term. Since d <= 2, the forced divisor 3 is always eligible, so k is 2 or 3 at every decomposable term (24,999 and 74,999 of them) and the cloud collapses to two vertical lines. Only a = 4, 5 and 8 are level-classified: the most weight-dominated sequence in the atlas.",
  100000, 1, gen_non3 },
{ "powerful", "A001694", "Powerful numbers", "multiplicative",
  "p | n implies p^2 | n; each term is a^2 b^3 with b squarefree, uniquely. count(x) ~ (zeta(3/2)/zeta(3)) sqrt(x) = 2.17 sqrt(x), so the sequence thins out like a quadratic one. The jump outgrows sqrt(l) on 42.07 % of terms (forced level) and 71.90 % are level-classified.",
  100000, 1, gen_powerful },
{ "sorting", "A001855", "Sorting numbers (binary insertion)", "summatory",
  "d = ceil(log2 n) exactly, constant on dyadic blocks. l/d^2 grows without bound so the level-forcing criterion never fires: 68 % weight. The gap sets a left wall - no weight column below k = d can exist.",
  100000, 1, gen_sorting },
{ "totientsum", "A002088", "Sum of totients", "summatory",
  "a(n) = sum of phi(k) for k <= n, from a(0) = 0; a ~ 3n^2/pi^2 and d = phi(n+1) <= n, so the sequence is nearly quadratic. l <= d^2 on 45.98 % of terms and 87.35 % are level-classified. For n >= 2 both a and d are even, so 2 | l. Since a > 3d on every decomposable term, l/2 > d, so k <= l/2 and L >= 2: the line L = 1 is empty.",
  100000, 0, gen_totsum },
{ "palindromes", "A002113", "Base-10 palindromes", "digit rule",
  "Starts at a(1) = 0. Only 13 distinct gaps below 10^12, all of the form 10^j or 11*10^j (plus the single value 2). The plane is self-similar with period 2 in the digit length; the level share oscillates 60/94 % with its parity.",
  100000, 1, gen_palindromes },
{ "primes3mod4", "A002145", "Primes congruent to 3 mod 4", "primes",
  "Every term is 3 mod 4 and every gap 0 mod 4, so l = 3 (mod 4). No square is 3 mod 4: zero ties, by proof. l is odd, so every weight is odd. The level share, 28.16 %, is above the primes' 23.00 %.",
  100000, 1, gen_primes3mod4 },
{ "oblong", "A002378", "Oblong numbers", "polynomial",
  "n(n + 1) = 2 T(n): d = 2n + 2 and l = (n - 2)(n + 1) < d^2, forced level (100 %). l is even and l/2 > d for n > 6, so L >= 2: the line L = 1 holds only a = 30 and 42.",
  100000, 0, gen_oblong },
{ "quartersq", "A002620", "Quarter-squares", "polynomial",
  "floor(n^2/4), from a(1) = 0 (a(0) = a(1) = 0). The gap takes each value twice and l <= d^2 throughout: forced level at every decomposable term (100 %).",
  100000, 1, gen_quartersq },
{ "composites", "A002808", "Composite numbers", "complement",
  "Composites: gaps of 1, and 2 across each prime. The naturals' plane, thinned: 10.58 % level against 9.59 %.",
  100000, 1, gen_composites },
{ "ulam", "A002858", "Ulam numbers", "self-referential",
  "Rigid: gap 2 alone carries 37 % of the terms. The signal at alpha = 2.5714475 acts on the classification only through the gap - across phase bins mean gap and level share correlate at r = 0.9957.",
  100000, 1, gen_ulam },
{ "self", "A003052", "Self (Colombian) numbers", "digit rule",
  "Self numbers: not of the form m + digitsum(m). They come in runs spaced 11 apart - gap 11 on 89.8 % of the terms - broken at decade boundaries. The level share is 24.93 %.",
  100000, 1, gen_self },
{ "fibbinary", "A003714", "Fibbinary numbers", "binary rule",
  "No two adjacent 1 bits in binary. The decomposition fails 24 times, each where the sequence jumps to the next power of two (a(n+1) >= 1.5 a(n)). Gap 1 carries 62 % of the terms and the level share is 10.35 %.",
  100000, 0, gen_fibbinary },
{ "skiptake", "A004202", "Skip 1, take 1, skip 2, take 2, ...", "block",
  "Closed form { m^2+j : 1 <= j <= m }. Interior terms (d = 1) follow the naturals rule; block openers m^2+1 give a tie exactly when m is prime; every block closer satisfies l <= d^2 and is level.",
  100000, 1, gen_skiptake },
{ "digitadd", "A004207", "Digit-addition trajectory", "digit rule",
  "a(n) = a(n-1) + digitsum(a(n-1)), so 9 divides l at every term. Above a = 800 every level term lies on a line L in 9Z (a level term off 9Z needs l <= 3d^2, impossible past a = 10^4), and Lemma 3 (L <= d) caps the lines at the gap: with d <= 50 in this range only L = 9, 18, 27, 36 occur. The 12 level terms off 9Z are all below 800.",
  100000, 1, gen_digitadd },
{ "abundant", "A005101", "Abundant numbers", "divisor functions",
  "sigma(n) > 2n. 99.2 % of these terms are even (the first odd one is 945), and so is l on 99.2 % of decomposable terms: the level class sits on even L, and the line L = 1 holds only 32 terms.",
  100000, 1, gen_abundant },
{ "squarefree", "A005117", "Squarefree numbers", "multiplicative",
  "Squarefree numbers, density 6/pi^2. Gaps of 1 to 7; the plane is close to the naturals' (10.72 % level) and L = 1 carries almost all of the level class.",
  100000, 1, gen_squarefree },
{ "practical", "A005153", "Practical numbers", "divisor functions",
  "Practical numbers. Every term > 2 is divisible by 4 or 6, so a and d are even and 2 | l at every decomposable term: the line L = 1 holds only 2 terms and the level class sits on L = 4, 6, 8, ... The level share is 17.50 %.",
  100000, 1, gen_practical },
{ "trisq", "A005214", "Triangular numbers and squares", "polynomial",
  "The union of two sequences that are 100 % level-classified is not: merging triangular numbers and squares shrinks the gaps below sqrt(l) on 58.6 % of terms, and 16.2 % turn weight. No ties occur.",
  100000, 1, gen_trisq },
{ "harshad", "A005349", "Harshad numbers", "digit rule",
  "n divisible by its digit sum. Density by digit sum is 9:3:1 according to gcd(s,9), but that structure does not reach the plane: the weight fence is flat and the level share varies only between 11 and 24 %.",
  100000, 1, gen_harshad },
{ "odd", "A005408", "Odd numbers", "forced divisor",
  "d = 2 and l = a - 2 is odd, so the divisor 2 is never available and k = spf(a - 2). The level class is exactly {a : a - 2 prime}, all on L = 1: 17,982 = pi(199,997) - 1 terms here. The ties are a - 2 = p^2, 85 = pi(447) - 1. Removing the divisor 2 nearly doubles the level share against the naturals: 17.98 % against 9.59 %.",
  100000, 0, gen_odd },
{ "even", "A005843", "Even numbers", "arithmetic progression",
  "d = 2 and l = a - 2 = 2m. The level class is exactly {a : m an odd prime}, all on the line L = 2, plus a = 6 (L = 1) and a = 10: the naturals' level line, moved from L = 1 to L = 2 (9.59 % level, the same share).",
  100000, 0, gen_even },
{ "divsum", "A006218", "Sum of d(k) for k <= n", "summatory",
  "Sum of d(k) for k <= n, ~ n log n; the gap is d(n+1) (1 to 128 here). The level share, 23.41 %, is close to the primes' 23.00 %.",
  100000, 0, gen_divsum },
{ "floorsqrt", "A006446", "Numbers k with floor(sqrt k) | k", "block",
  "Blocks {m^2, m^2 + m, m^2 + 2m} with gaps m, m, 1. The first two terms of each block have l <= d^2 and are forced level - exactly two thirds of the sequence (66.67 %). The third, with d = 1 and l = (m+1)^2 - 2, never a square, behaves like a natural number: no ties.",
  100000, 1, gen_floorsqrt },
{ "twinhi", "A006512", "Greater of twin primes", "primes",
  "Past 5, every term is 1 mod 6, so l = 1 (mod 6): squares are possible and 29 ties occur, where the lesser twins (l = 5 mod 6) have none. The level share, 46.93 %, is close to the lesser twins' 47.77 %.",
  100000, 1, gen_twinhi },
{ "sfsemiprimes", "A006881", "Squarefree semiprimes", "multiplicative",
  "Products of two distinct primes. The plane is the semiprimes' almost exactly (15.72 % level against 15.71 %).",
  100000, 1, gen_sfsemiprimes },
{ "binary", "A007088", "Binary expansion read in decimal", "binary rule",
  "n in binary, read in decimal. The gap is (8*10^j + 1)/9 = 1, 9, 89, 889, ... where j counts the trailing 1 bits of n: sixteen values here, and the plane splits into clean diagonal lines. Capped at 65,535 terms, since a(65,536) = 10^16 passes 2^53.",
  65535, 0, gen_binary },
{ "primesum", "A007504", "Sum of the first n primes", "summatory",
  "Sum of the first n primes: a ~ n^2 log n / 2 and d = p(n+1), so l/d^2 -> 0 like 1/(2 log n): every decomposable term is forced level (100 %).",
  100000, 0, gen_primesum },
{ "pplus1", "A008864", "Primes plus one", "primes",
  "p + 1. The gap is the prime gap and l = p + 1 - g is even, so 2 | l; l/2 > d on all but two terms (a = 3, 6), and the primes' level line L = 1 moves to L = 2, which holds 43 % of the level class.",
  100000, 1, gen_pplus1 },
{ "notpp", "A024619", "Not prime powers", "complement",
  "Numbers with at least two distinct prime factors. Gaps of 1 to 4; the plane is close to the naturals' (10.59 % level).",
  100000, 1, gen_notpp },
{ "sigmasum", "A024916", "Sum of sigma(k) for k <= n", "summatory",
  "Sum of sigma(k) for k <= n: a ~ pi^2 n^2/12 = 0.822 n^2 and d = sigma(n+1) > n: every decomposable term is forced level (100 %).",
  100000, 1, gen_sigmasum },
{ "oddomega", "A026424", "Odd number of prime factors", "multiplicative",
  "Omega(n) odd (Liouville lambda = -1). Gaps of 1 to 16; the level share is 12.82 %.",
  100000, 1, gen_oddomega },
{ "zeroless", "A052382", "Zeroless numbers", "digit rule",
  "No digit 0. Gap 1 almost everywhere, with jumps 2, 12, 112, ... across the zeros; the plane is the naturals' (10.54 % level).",
  100000, 1, gen_zeroless },
{ "odious", "A000069", "Odious numbers", "binary rule",
  "Numbers with an odd number of 1 bits. Gaps are 1, 2 and 3, each on a third of the terms. The level share is 11.01 %, next to the evil numbers' 11.39 %, with 86 ties each. 84 % of the level class sits on L = 1 and the rest on L = 3.",
  100000, 1, gen_odious },
{ "evil", "A001969", "Evil numbers", "binary rule",
  "Numbers with an even number of 1 bits, the complement of the odious numbers. Gaps are 1, 2 and 3, each on a third of the terms. The level share is 11.39 % against 11.01 % for the odious numbers. The level class sits on L = 1 (76 %) and L = 3.",
  100000, 1, gen_evil },
{ "heptagonal", "A000566", "Heptagonal numbers", "polynomial",
  "n(5n - 3)/2. For a ~ c n^2, l/d^2 -> 1/(4c); here c = 5/2, so l <= d^2 and every decomposable term is forced level (100 %), as for the squares and pentagonal numbers.",
  100000, 0, gen_heptagonal },
{ "octagonal", "A000567", "Octagonal numbers", "polynomial",
  "n(3n - 2), a ~ 3 n^2 and d = 6n + 1: l < d^2, so every decomposable term is forced level (100 %). The lines L = 1 and L = 2 are about equally full (6.1 % and 6.0 % of the level class).",
  100000, 0, gen_octagonal },
{ "censquare", "A001844", "Centered square numbers", "polynomial",
  "2n(n + 1) + 1 = n^2 + (n + 1)^2. d = 4(n + 1) and a ~ 2 n^2 < d^2: forced level at every decomposable term (100 %). Every l is odd and the fullest line is L = 3 (10.9 %), ahead of L = 1 (8.4 %).",
  100000, 0, gen_censquare },
{ "cenhex", "A003215", "Centered hexagonal numbers", "polynomial",
  "3n(n + 1) + 1, the hex numbers. d = 6(n + 1) and a ~ 3 n^2 < d^2: forced level at every decomposable term (100 %). L = 1 holds 13.8 % of the class, then L = 5 (6.9 %).",
  100000, 0, gen_cenhex },
{ "centri", "A005448", "Centered triangular numbers", "polynomial",
  "3n(n - 1)/2 + 1. d = 3n and a ~ 1.5 n^2 < d^2: forced level at every decomposable term (100 %).",
  100000, 1, gen_centri },
{ "lazycaterer", "A000124", "Lazy caterer's sequence", "polynomial",
  "n(n + 1)/2 + 1, the triangular numbers moved up by one. With c = 1/2 > 1/4 the terms are all forced level (100 %), like the triangular numbers, but the line L = 1 is much thinner: 1,342 terms against 2,642.",
  100000, 0, gen_lazycaterer },
{ "sqplus1", "A002522", "Squares plus one", "polynomial",
  "n^2 + 1. d = 2n + 1 and l = n^2 - 2n < d^2: every decomposable term is forced level (100 %). The line L = 1 (l prime) is thin: 1,225 terms.",
  100000, 0, gen_sqplus1 },
{ "genpent", "A001318", "Generalized pentagonal numbers", "polynomial",
  "m(3m - 1)/2 for m = 0, 1, -1, 2, -2, ... Two quadratics interleaved: the gaps alternate j and 2j + 1. With a ~ 1.5 j^2 the gap j leaves l > d^2 while the gap 2j + 1 forces level, so exactly half the terms (49,999) are forced level. The other half is split, and the level share ends at 77.61 %.",
  100000, 0, gen_genpent },
{ "pentpyr", "A002411", "Pentagonal pyramidal numbers", "polynomial",
  "n^2 (n + 1)/2, the pentagonal pyramidal numbers. d grows like n^2 while l grows like n^3, so l/d^2 -> 0: every decomposable term is forced level (100 %).",
  100000, 0, gen_pentpyr },
{ "octahedral", "A005900", "Octahedral numbers", "polynomial",
  "n(2n^2 + 1)/3. A cubic sequence: l/d^2 -> 0 and every decomposable term is forced level (100 %). L = 1 holds 9.4 % of the terms and L = 2 half as many.",
  100000, 0, gen_octahedral },
{ "upperwythoff", "A001950", "Upper Wythoff sequence", "Beatty",
  "floor(n phi^2) = floor(n phi) + n, the complement of the lower Wythoff sequence. The gaps are 2 and 3 only (38 % and 62 %). The level share, 15.38 %, is above the lower sequence's 12.20 %; the level class lies on the lines L = 1, 2 and 3.",
  100000, 1, gen_upperwythoff },
{ "beattysqrt2", "A001951", "Beatty sequence floor(n sqrt 2)", "Beatty",
  "floor(n sqrt 2), a Beatty sequence with gaps 1 and 2 only. The level share is 11.37 %, on the lines L = 1 (82 %) and L = 2.",
  100000, 0, gen_beattysqrt2 },
{ "beattysqrt3", "A022838", "Beatty sequence floor(n sqrt 3)", "Beatty",
  "floor(n sqrt 3), a Beatty sequence with gaps 1 and 2 only (27 % and 73 %). The level share is 12.69 %, on the lines L = 1 (72 %) and L = 2.",
  100000, 1, gen_beattysqrt3 },
{ "primes1mod4", "A002144", "Primes congruent to 1 mod 4", "primes",
  "Every term is 1 mod 4 and every gap 0 mod 4, so l = 1 (mod 4) and ties are possible (23 here), where the primes 3 mod 4 have none. The level share, 27.95 %, is close to theirs (28.16 %). The fullest level line is L = 3.",
  100000, 1, gen_primes1mod4 },
{ "primes1mod6", "A002476", "Primes congruent to 1 mod 6", "primes",
  "Every term is 1 mod 6 and every gap 0 mod 6, so l = 1 (mod 6). The level share, 35.24 %, is almost the primes 5 mod 6's 35.16 %, but they have no ties and these have 48.",
  100000, 1, gen_primes1mod6 },
{ "primes5mod6", "A007528", "Primes congruent to 5 mod 6", "primes",
  "Every term is 5 mod 6 and every gap 0 mod 6, so l = 5 (mod 6). No square is 5 mod 6: zero ties, by proof. The level share is 35.16 %, and half the level class sits on L = 1.",
  100000, 1, gen_primes5mod6 },
{ "sophie", "A005384", "Sophie Germain primes", "primes",
  "p with 2p + 1 prime. Past 3 every term is 5 mod 6, so l = 5 (mod 6) and no tie can occur, as for the lesser twin primes. The level share, 47.25 %, is close to theirs (47.77 %).",
  100000, 1, gen_sophie },
{ "safe", "A005385", "Safe primes", "primes",
  "p with (p - 1)/2 prime. Past 7 every term is 11 mod 12 and every gap 0 mod 12, so l = 3 (mod 4): no square, no tie. The level share, 51.96 %, is the highest among the prime subsequences here.",
  100000, 1, gen_safe },
{ "primeidx", "A006450", "Prime-indexed primes", "primes",
  "prime(prime(n)). The gaps are larger and more varied than the primes' (662 distinct values here), and the level share rises to 43.29 %, against 23.00 % for the primes.",
  100000, 1, gen_primeidx },
{ "isolated", "A007510", "Isolated (single) primes", "primes",
  "Primes p with neither p - 2 nor p + 2 prime. Level share 26.11 %, a little above the primes' 23.00 %; L = 1 and L = 3 carry 60 % of the level class.",
  100000, 1, gen_isolated },
{ "twins", "A001097", "Twin primes", "primes",
  "Primes in a twin pair, both members. Past 5 every l is 3 mod 6: after a lesser twin d = 2 and l = p - 2; after a greater twin d = 4 (mod 6). So 3 | l, and the line L = 1 holds only 6 terms. All but 68 level terms sit on odd multiples of 3 (L = 3, 9, 15, ...). The level share is 19.67 %.",
  100000, 1, gen_twins },
{ "cousin", "A023200", "Lesser of cousin primes", "primes",
  "p with p + 4 prime. Past 3 every term is 1 mod 6, so l = 1 (mod 6) and ties can occur (31 here). The level share, 46.77 %, is close to the lesser twins' 47.77 %.",
  100000, 1, gen_cousin },
{ "sexy", "A023201", "Lesser of sexy primes", "primes",
  "p with p + 6 prime. The terms mix 1 and 5 mod 6, and the level share, 34.58 %, lies between the primes' and the twin primes'.",
  100000, 1, gen_sexy },
{ "emirps", "A006567", "Emirps", "primes",
  "Primes whose decimal reversal is a different prime. No emirp begins with 2, 4, 5, 6 or 8 (its reversal would be even or a multiple of 5), so whole blocks are skipped. The largest gap is 3,000,162, and 12 terms fail to decompose. The level share is 31.54 %.",
  100000, 1, gen_emirps },
{ "pminus1", "A006093", "Primes minus one", "primes",
  "p - 1. Past 2 every term and every gap is even, so 2 | l. The primes' level line L = 1 moves to L = 2, which holds 45 % of the level class (compare p + 1). There are no ties, and the line L = 1 holds a single term.",
  100000, 1, gen_primesminus1 },
{ "almost3", "A014612", "3-almost primes", "multiplicative",
  "Omega(n) = 3. The gaps stay between 1 and 34 here. The level share is 15.79 %, beside the semiprimes' 15.71 %, and L = 1 carries 64 % of the level class.",
  100000, 1, gen_almost3 },
{ "almost4", "A014613", "4-almost primes", "multiplicative",
  "Omega(n) = 4. The level share is 16.89 %. The lines L = 1 and L = 2 are nearly equal (36 % and 33 % of the level class); 65 % of the terms are even.",
  100000, 1, gen_almost4 },
{ "sphenic", "A007304", "Sphenic numbers", "multiplicative",
  "Products of three distinct primes. The level share is 17.04 %, against 15.72 % for the squarefree semiprimes.",
  100000, 1, gen_sphenic },
{ "nonsqfree", "A013929", "Non-squarefree numbers", "multiplicative",
  "Numbers divisible by a square > 1. Only four gaps occur here, 1 to 4. The level share is 13.03 %, spread over L = 1, 2, 4 and 3.",
  100000, 1, gen_nonsqfree },
{ "weak", "A052485", "Weak numbers (not powerful)", "multiplicative",
  "Numbers that are not powerful (some prime divides them exactly once). Gap 1 on 99.4 % of the terms. The plane is the naturals' almost unchanged: 9.62 % level against 9.59 %.",
  100000, 1, gen_weak },
{ "deficient", "A005100", "Deficient numbers", "divisor functions",
  "sigma(n) < 2n. Gaps 1 and 2, with 38 gaps of 3. Every level-classified term lies on L = 1, and the level share, 7.79 %, is the lowest here after the numbers not divisible by 3.",
  100000, 1, gen_deficient },
{ "loeschian", "A003136", "Loeschian numbers", "quadratic form",
  "x^2 + xy + y^2: every prime 2 mod 3 divides n to an even power. Starts at a(1) = 0. The level share, 24.23 %, is well above the sums of two squares' 15.95 %.",
  100000, 1, gen_loeschian },
{ "x2p2y2", "A002479", "Numbers x^2 + 2y^2", "quadratic form",
  "x^2 + 2y^2: every prime 5 or 7 mod 8 divides n to an even power. Starts at a(1) = 0. Level share 14.91 %, with half the level class on L = 1.",
  100000, 1, gen_x2p2y2 },
{ "twopossq", "A000404", "Sums of two nonzero squares", "quadratic form",
  "x^2 + y^2 with x, y >= 1. The plane is almost exactly that of all sums of two squares (15.97 % level against 15.95 %).",
  100000, 1, gen_twopossq },
{ "threesq", "A000378", "Sums of three squares", "quadratic form",
  "Sums of three squares, the numbers not of the form 4^a(8b + 7) (Legendre). Density 5/6 and gaps of 1 to 3. The level share is 8.54 %, 99 % of it on L = 1.",
  100000, 1, gen_threesq },
{ "notthreesq", "A004215", "Not sums of three squares", "quadratic form",
  "4^a(8b + 7), the complement of the sums of three squares. Half the gaps are 8. The level share is 23.71 %, with L = 1 holding half the level class.",
  100000, 1, gen_notthreesq },
{ "coprime30", "A007775", "Numbers coprime to 30", "forced divisor",
  "Numbers prime to 2, 3 and 5. The gaps are 2, 4 and 6, repeating with period 8. The level class lies on only three lines: L = 3 (43 %), L = 1 and L = 5.",
  100000, 1, gen_coprime30 },
{ "coprime6", "A007310", "Numbers coprime to 6", "forced divisor",
  "6m + 1 and 6m + 5, with gaps 4 and 2. Then 3 | l at every term, so after a gap of 2 the weight is 3. After a gap of 4, a term is level essentially when l/3 is prime. The level class is the naturals' level line moved to L = 3 (9.59 % level).",
  100000, 1, gen_coprime6 },
{ "mult3", "A008585", "Multiples of 3", "forced divisor",
  "3n: d = 3 and l = 3(n - 1). The level class is essentially {3n : n - 1 prime}, all on the line L = 3, like the even numbers' line L = 2. The level share is 9.60 %.",
  100000, 0, gen_mult3 },
{ "oddnonprimes", "A014076", "Odd nonprimes", "complement",
  "1 and the odd composites. Gaps 2, 4, 6 and 8. The level share, 21.01 %, is above the odd numbers' 17.98 %; 97 % of the level class is on L = 1.",
  100000, 1, gen_oddnonprimes },
{ "happy", "A007770", "Happy numbers", "digit rule",
  "The digit-square-sum iteration reaches 1. The gaps are irregular (58 distinct values up to 73) and the level share is 19.71 %.",
  100000, 1, gen_happy },
{ "haszero", "A011540", "Numbers with a digit 0", "digit rule",
  "Numbers with at least one digit 0; starts at a(1) = 0. Only the gaps 1 and 10 occur here, and the level class lies on the two lines L = 1 and L = 10. The level share is 8.76 %.",
  100000, 1, gen_haszero },
{ "odddigits", "A014261", "Numbers with only odd digits", "digit rule",
  "Every decimal digit odd. The gaps are 2, 12, 112, 1112, ... (80 % are 2). The level share is 17.40 %.",
  100000, 1, gen_odddigits },
{ "ternary", "A007089", "Ternary expansion read in decimal", "digit rule",
  "n in base 3, read in decimal. The gaps are 1, 8, 78, 778, ..., one per count of trailing 2 digits, so the plane splits into clean lines. 21 terms fail to decompose, and the level share is 10.70 %.",
  100000, 0, gen_ternary },
{ "cantor", "A005836", "No digit 2 in base 3", "digit rule",
  "No digit 2 in base 3: the integers of the Cantor set, and the greedy sequence with no three terms in arithmetic progression. The gaps are (3^j + 1)/2 = 1, 2, 5, 14, ... Level share 14.75 %, and no ties.",
  100000, 1, gen_cantor },
{ "moser", "A000695", "Moser-de Bruijn sequence", "binary rule",
  "Sums of distinct powers of 4: binary n read in base 4. The gaps are (2*4^j + 1)/3 = 1, 3, 11, 43, ... The level share is 13.96 %; L = 1 and L = 2 hold 73 % of the level class. No ties.",
  100000, 0, gen_moser },
{ "evenzeros", "A003159", "Even number of trailing 0 bits", "binary rule",
  "n with an even number of trailing 0 bits (the 2-adic valuation is even). Gaps 1 and 2. Every level term lies on L = 1, there are no ties, and the level share is 13.88 %.",
  100000, 1, gen_evenzeros },
{ "binpal", "A006995", "Binary palindromes", "binary rule",
  "Binary palindromes; starts at a(1) = 0. There are about sqrt(x) of them up to x, so the sequence thins out like a quadratic. The gap outgrows sqrt(l) on 35 % of terms (forced level), and 88.40 % are level-classified.",
  100000, 1, gen_binpal },
{ "nnp3", "A000096", "n(n+3)/2", "polynomial",
  "n(n+3)/2, the triangular numbers minus one, shifted: T(n+1) - 1. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_nnp3 },
{ "nnp2", "A005563", "n(n+2)", "polynomial",
  "n(n+2) = (n+1)^2 - 1. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_nnp2 },
{ "twicesq", "A001105", "Twice squares", "polynomial",
  "2n^2. Every term and every gap is even, so 2 | l, and the line L = 1 holds only 2 terms. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_twicesq },
{ "oddsq", "A016754", "Odd squares", "polynomial",
  "(2n+1)^2. d = 8(n+1) and l = (2n+1)^2 - 8(n+1) is odd. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_oddsq },
{ "genoct", "A001082", "Generalized octagonal numbers", "polynomial",
  "m(3m - 2) for m = 0, 1, -1, 2, -2, ... Two quadratics interleaved: the gaps alternate 4j and 2j + 1 while a ~ 3j^2, so both leave l <= d^2. Every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_genoct },
{ "cenpent", "A005891", "Centered pentagonal numbers", "polynomial",
  "(5n^2 + 5n + 2)/2, the centered pentagonal numbers. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_cenpent },
{ "cenocta", "A001845", "Centered octahedral numbers", "polynomial",
  "(2n+1)(2n^2+2n+3)/3, the centered octahedral numbers. A cubic sequence. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_cenocta },
{ "cencube", "A005898", "Centered cube numbers", "polynomial",
  "n^3 + (n+1)^3, the centered cube numbers. A cubic sequence. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_cencube },
{ "stella", "A007588", "Stella octangula numbers", "polynomial",
  "n(2n^2 - 1), the stella octangula numbers. A cubic sequence. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_stella },
{ "oddsqsum", "A000447", "Sums of the first n odd squares", "polynomial",
  "1^2 + 3^2 + ... + (2n-1)^2 = n(2n-1)(2n+1)/3. A cubic sequence. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_oddsqsum },
{ "hexpyr", "A002412", "Hexagonal pyramidal numbers", "polynomial",
  "n(n+1)(4n-1)/6, the hexagonal pyramidal numbers. A cubic sequence. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_hexpyr },
{ "cake", "A000125", "Cake numbers", "polynomial",
  "(n^3 + 5n + 6)/6, the most pieces from n plane cuts of a cube. A cubic sequence with d = T(n) + 1. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_cake },
{ "magic", "A006003", "Magic constants n(n^2+1)/2", "polynomial",
  "n(n^2 + 1)/2, the magic constant of an n x n magic square. A cubic sequence. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_magic },
{ "sechex", "A014105", "Second hexagonal numbers", "polynomial",
  "n(2n+1), the second hexagonal numbers. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_sechex },
{ "secpent", "A005449", "Second pentagonal numbers", "polynomial",
  "n(3n+1)/2, the second pentagonal numbers. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_secpent },
{ "trimatch", "A045943", "Triangular matchstick numbers", "polynomial",
  "3n(n+1)/2 = 3 T(n). 3 | a and 3 | d, so 3 | l, and the line L = 1 holds only 3 terms. Every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_trimatch },
{ "thirdsq", "A000212", "floor(n^2/3)", "polynomial",
  "floor(n^2/3), from a(1) = 0 (a(0) = a(1) = 0). Every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_thirdsq },
{ "halfsq", "A007590", "floor(n^2/2)", "polynomial",
  "floor(n^2/2), from a(1) = 0 (a(0) = a(1) = 0). Terms and gaps are even past the start, so the line L = 1 holds 2 terms. Every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_halfsq },
{ "pow32", "A000093", "floor(n^(3/2))", "polynomial",
  "floor(n^(3/2)). The gap grows like sqrt(n) while l grows like n^(3/2), so l/d^2 grows like sqrt(n) and almost no term is forced level. The level share is 46.64 %.",
  100000, 0, gen_pow32 },
{ "p1mod8", "A007519", "Primes congruent to 1 mod 8", "primes",
  "Primes 1 mod 8. Every gap is 0 mod 8, so l = 1 (mod 8) and ties are possible. The level share is 33.30 %.",
  100000, 1, gen_p1mod8 },
{ "p3mod8", "A007520", "Primes congruent to 3 mod 8", "primes",
  "Primes 3 mod 8. Every gap is 0 mod 8, so l = 3 (mod 8); no square is 3 mod 4, so there are no ties, by proof. The level share is 33.29 %.",
  100000, 1, gen_p3mod8 },
{ "p5mod8", "A007521", "Primes congruent to 5 mod 8", "primes",
  "Primes 5 mod 8. Every gap is 0 mod 8, so l = 5 (mod 8); no square is 5 mod 8, so there are no ties, by proof. The level share is 33.29 %.",
  100000, 1, gen_p5mod8 },
{ "p7mod8", "A007522", "Primes congruent to 7 mod 8", "primes",
  "Primes 7 mod 8. Every gap is 0 mod 8, so l = 7 (mod 8); no square is 3 mod 4, so there are no ties, by proof. The level share is 33.13 %.",
  100000, 1, gen_p7mod8 },
{ "pend1", "A030430", "Primes ending in 1", "primes",
  "Primes ending in 1. Every gap is a multiple of 10, so l ends in 1 as well. The level share is 37.20 %.",
  100000, 1, gen_pend1 },
{ "pend3", "A030431", "Primes ending in 3", "primes",
  "Primes ending in 3. Every gap is a multiple of 10, so l ends in 3; no square ends in 3, so there are no ties, by proof. The level share is 37.28 %.",
  100000, 1, gen_pend3 },
{ "pend7", "A030432", "Primes ending in 7", "primes",
  "Primes ending in 7. Every gap is a multiple of 10, so l ends in 7; no square ends in 7, so there are no ties, by proof. The level share is 37.30 %.",
  100000, 1, gen_pend7 },
{ "pend9", "A030433", "Primes ending in 9", "primes",
  "Primes ending in 9. Every gap is a multiple of 10, so l ends in 9 and ties are possible. The level share is 37.29 %.",
  100000, 1, gen_pend9 },
{ "balanced", "A006562", "Balanced primes", "primes",
  "Primes that are the mean of their two neighbours. They are sparse (the 100,000th is 56,206,697), and the gaps are large and varied. The level share is 52.42 %; 1.0 % of terms are forced level (l <= d^2).",
  100000, 1, gen_balanced },
{ "strong", "A051634", "Strong primes", "primes",
  "Primes above the mean of their two neighbours. The level share is 29.38 %.",
  100000, 1, gen_strong },
{ "weakp", "A051635", "Weak primes", "primes",
  "Primes below the mean of their two neighbours. The level share is 29.87 %.",
  100000, 1, gen_weakp },
{ "p2pm1", "A005382", "Primes p with 2p - 1 prime", "primes",
  "Primes p with 2p - 1 prime. The level share is 47.76 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_p2pm1 },
{ "phalf", "A005383", "Primes p with (p + 1)/2 prime", "primes",
  "Primes p with (p + 1)/2 prime. The level share is 52.16 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_phalf },
{ "pplus8", "A023202", "Primes p with p + 8 prime", "primes",
  "Primes p with p + 8 prime. Past 3 every term is 5 mod 6 (p + 8 must avoid 3), so l = 5 (mod 6): no ties, by proof. The level share is 47.00 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_pplus8 },
{ "gap6", "A031924", "Lower prime of a gap of 6", "primes",
  "Primes whose next prime is p + 6. The level share is 36.03 %.",
  100000, 1, gen_gap6 },
{ "trip26", "A022004", "Prime triplets (p, p+2, p+6)", "primes",
  "p with p, p + 2, p + 6 all prime. Past 5 every term is 5 mod 6, so l = 5 (mod 6): no ties, by proof. The level share is 63.42 %; 3.9 % of terms are forced level (l <= d^2).",
  100000, 1, gen_trip26 },
{ "trip46", "A022005", "Prime triplets (p, p+4, p+6)", "primes",
  "p with p, p + 4, p + 6 all prime. Past 7 every term is 1 mod 6. The level share is 66.20 %; 3.8 % of terms are forced level (l <= d^2); L = 1 holds 31 % of the level class.",
  100000, 1, gen_trip46 },
{ "pprod", "A006094", "Products of two successive primes", "primes",
  "p(n) p(n+1). d = p(n+1)(p(n+2) - p(n)) is at least 2 p(n+1), so d^2 > a. Every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_pprod },
{ "npn", "A033286", "n times the n-th prime", "primes",
  "n p(n). d = n (p(n+1) - p(n)) + p(n+1) > p(n) > n, so d^2 > n p(n) = a. Every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_npn },
{ "pnplusn", "A014688", "n-th prime plus n", "primes",
  "p(n) + n. The gap is the prime gap plus 1. The level share is 25.35 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_pnplusn },
{ "psum2", "A001043", "Sums of two successive primes", "primes",
  "p(n) + p(n+1), all even past the first term. l is even, and the line L = 1 holds a single term. The level share is 24.90 %.",
  100000, 1, gen_psum2 },
{ "twop", "A100484", "Twice the primes", "primes",
  "2p. The primes' decomposition doubled: d = 2g and l = 2(p - g), so the primes' level line L = 1 moves to L = 2. 23,004 terms are level-classified, against 22,999 for the primes. The level share is 23.00 %; L = 2 holds 32 % of the level class; there are no ties.",
  100000, 1, gen_twop },
{ "threep", "A001748", "Three times the primes", "primes",
  "3p. The level share is almost the primes' 23.00 %. The level share is 23.02 %; L = 3 holds 32 % of the level class.",
  100000, 1, gen_threep },
{ "pq2", "A054753", "Numbers p q^2", "multiplicative",
  "p q^2 with distinct primes p, q. The level share is 27.45 %.",
  100000, 1, gen_pq2 },
{ "fourdist", "A046386", "Products of four distinct primes", "multiplicative",
  "Products of four distinct primes. The level share is 21.02 %; L = 2 holds 34 % of the level class.",
  100000, 1, gen_fourdist },
{ "almost5", "A014614", "5-almost primes", "multiplicative",
  "Omega(n) = 5. The level share is 17.97 %.",
  100000, 1, gen_almost5 },
{ "cubefree", "A004709", "Cubefree numbers", "multiplicative",
  "No cube > 1 divides n. Density 1/zeta(3) = 0.83. The level share is 10.16 %; L = 1 holds 99 % of the level class.",
  100000, 1, gen_cubefree },
{ "noncubefree", "A046099", "Non-cubefree numbers", "multiplicative",
  "Divisible by a cube > 1. The level share is 14.06 %; L = 8 holds 37 % of the level class.",
  100000, 1, gen_noncubefree },
{ "omega2", "A007774", "Exactly two distinct prime factors", "multiplicative",
  "Exactly two distinct prime factors, any exponents. The level share is 13.17 %; L = 1 holds 66 % of the level class.",
  100000, 1, gen_omega2 },
{ "omega3", "A033992", "Exactly three distinct prime factors", "multiplicative",
  "Exactly three distinct prime factors, any exponents. The level share is 13.55 %; L = 1 holds 60 % of the level class.",
  100000, 1, gen_omega3 },
{ "omega4", "A033993", "Exactly four distinct prime factors", "multiplicative",
  "Exactly four distinct prime factors, any exponents. The level share is 17.16 %; L = 2 holds 37 % of the level class.",
  100000, 1, gen_omega4 },
{ "mu1", "A030229", "Moebius mu(n) = 1", "multiplicative",
  "Squarefree with an even number of prime factors (Moebius mu = 1). The plane is nearly the mu = -1 plane (14.67 % level against 14.60 %). The level share is 14.67 %; L = 1 holds 71 % of the level class.",
  100000, 1, gen_mu1 },
{ "mum1", "A030059", "Moebius mu(n) = -1", "multiplicative",
  "Squarefree with an odd number of prime factors (Moebius mu = -1). The plane is nearly the mu = 1 plane. The level share is 14.60 %; L = 1 holds 71 % of the level class.",
  100000, 1, gen_mum1 },
{ "evenomega", "A028260", "Even number of prime factors", "multiplicative",
  "Omega(n) even (Liouville lambda = 1), the complement of the odd-Omega numbers (12.82 % level there). The level share is 12.78 %; L = 1 holds 70 % of the level class.",
  100000, 1, gen_evenomega },
{ "arith", "A003601", "Arithmetic numbers", "divisor functions",
  "The mean of the divisors, sigma(n)/d(n), is an integer. The level share is 11.93 %; L = 1 holds 97 % of the level class.",
  100000, 1, gen_arith },
{ "totients", "A002202", "Totients (values of phi)", "divisor functions",
  "Values of Euler's phi. Past 1 every totient is even, so 2 | l and the line L = 1 holds 2 terms. The level share is 12.58 %; L = 4 holds 40 % of the level class.",
  100000, 1, gen_totients },
{ "nontotients", "A007617", "Nontotients", "divisor functions",
  "Not a value of Euler's phi: every odd number > 1 and some even ones. Gaps 1 and 2. The level share is 11.10 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_nontotients },
{ "base4", "A007090", "Base 4 read in decimal", "digit rule",
  "n in base 4, read in decimal. The gaps are 1, 7, 67, 667, ... The level share is 12.65 %; L = 1 holds 59 % of the level class.",
  100000, 0, gen_base4 },
{ "base5", "A007091", "Base 5 read in decimal", "digit rule",
  "n in base 5, read in decimal. The gaps are 1, 6, 56, 556, ... The level share is 11.07 %; L = 1 holds 63 % of the level class.",
  100000, 0, gen_base5 },
{ "base6", "A007092", "Base 6 read in decimal", "digit rule",
  "n in base 6, read in decimal. The gaps are 1, 5, 45, 445, ... The level share is 10.58 %; L = 1 holds 90 % of the level class.",
  100000, 0, gen_base6 },
{ "nozero4", "A023705", "No digit 0 in base 4", "digit rule",
  "No digit 0 in base 4. The level share is 12.21 %; L = 1 holds 85 % of the level class.",
  100000, 1, gen_nozero4 },
{ "nozero3", "A032924", "No digit 0 in base 3", "digit rule",
  "No digit 0 in base 3: the digits are 1 and 2. After the numbers not divisible by 3 it has the lowest level share in the atlas, and the line L = 1 holds only 9 terms. The level share is 4.86 %; L = 3 holds 53 % of the level class; there are no ties.",
  100000, 1, gen_nozero3 },
{ "has1", "A011531", "Numbers with a digit 1", "digit rule",
  "At least one decimal digit 1. Only the gaps 1 to 10 occur. The level share is 11.06 %; L = 1 holds 90 % of the level class.",
  100000, 1, gen_has1 },
{ "has9", "A011539", "Numbers with a digit 9", "digit rule",
  "At least one decimal digit 9. Only the gaps 1 to 10 occur. The level share is 13.62 %; L = 1 holds 81 % of the level class.",
  100000, 1, gen_has9 },
{ "no1", "A052383", "No digit 1", "digit rule",
  "No decimal digit 1. Gaps jump across every block with a 1 (the largest here is 100,001). The level share is 9.31 %; L = 1 holds 79 % of the level class.",
  100000, 1, gen_no1 },
{ "no2", "A052404", "No digit 2", "digit rule",
  "No decimal digit 2. The level share is 10.38 %; L = 1 holds 97 % of the level class.",
  100000, 1, gen_no2 },
{ "evendig", "A014263", "Only even digits", "digit rule",
  "Every decimal digit even. The gaps are 2, 12, 112, ..., as for the only-odd-digit numbers. The level share is 11.08 %; L = 2 holds 63 % of the level class.",
  100000, 1, gen_evendig },
{ "nondec", "A009994", "Nondecreasing digits", "digit rule",
  "Decimal digits never decrease (0, 1, ..., 9, 11, 12, ...). There are only C(d+9, 9) of them with d digits or fewer, so the terms grow fast: the 100,000th has 11 digits. The level share is 19.84 %; 2.7 % of terms are forced level (l <= d^2); L = 1 holds 43 % of the level class.",
  100000, 1, gen_nondec },
{ "noninc", "A009996", "Nonincreasing digits", "digit rule",
  "Decimal digits never increase (0, ..., 9, 10, 11, 20, ...). The level share is 17.19 %; 2.0 % of terms are forced level (l <= d^2).",
  100000, 1, gen_noninc },
{ "dsprime", "A028834", "Prime digit sum", "digit rule",
  "The digit sum is prime. The level share is 10.71 %; L = 3 holds 44 % of the level class.",
  100000, 1, gen_dsprime },
{ "dseven", "A054683", "Even digit sum", "digit rule",
  "The digit sum is even. Gaps 1 to 3. The level share is 10.05 %; L = 1 holds 66 % of the level class.",
  100000, 1, gen_dseven },
{ "unhappy", "A031177", "Unhappy numbers", "digit rule",
  "The digit-square-sum iteration never reaches 1 (it falls into the cycle through 4). The level share is 10.38 %; L = 1 holds 92 % of the level class.",
  100000, 1, gen_unhappy },
{ "nonniven", "A065877", "Non-Harshad numbers", "digit rule",
  "Not divisible by the digit sum: the complement of the Harshad numbers. The level share is 8.25 %; L = 1 holds 99 % of the level class.",
  100000, 1, gen_nonniven },
{ "nonpal", "A029742", "Non-palindromes", "digit rule",
  "Not a decimal palindrome. Gaps 1 and 2 only. The level share is 9.60 %; L = 1 holds 99 % of the level class.",
  100000, 1, gen_nonpal },
{ "binbal", "A031443", "Balanced binary numbers", "binary rule",
  "As many 0 bits as 1 bits. The level share is 11.83 %; L = 1 holds 65 % of the level class.",
  100000, 1, gen_binbal },
{ "oddval", "A036554", "Odd number of trailing 0 bits", "binary rule",
  "An odd number of trailing 0 bits, the complement of A003159. Every term is twice a term of A003159. Both have exactly 13,875 level terms, all on L = 1 there and here almost all on L = 2. The level share is 13.88 %; L = 2 holds 100 % of the level class; there are no ties.",
  100000, 1, gen_oddval },
{ "beatty1s2", "A003151", "floor(n(1 + sqrt 2))", "Beatty",
  "floor(n(1 + sqrt 2)), the complement of floor(n(1 + 1/sqrt 2)). Gaps 2 and 3. The level share is 14.92 %; L = 1 holds 60 % of the level class.",
  100000, 1, gen_beatty_1s2 },
{ "beatty1h2", "A003152", "floor(n(1 + 1/sqrt 2))", "Beatty",
  "floor(n(1 + 1/sqrt 2)), the complement of floor(n(1 + sqrt 2)). Gaps 1 and 2. The level share is 12.55 %; L = 1 holds 73 % of the level class.",
  100000, 1, gen_beatty_1h2 },
{ "beatty2s2", "A001952", "floor(n(2 + sqrt 2))", "Beatty",
  "floor(n(2 + sqrt 2)), the complement of floor(n sqrt 2). Gaps 3 and 4. The level share is 17.24 %; L = 1 holds 50 % of the level class.",
  100000, 1, gen_beatty_2s2 },
{ "beattys5", "A022839", "floor(n sqrt 5)", "Beatty",
  "floor(n sqrt 5). Gaps 2 and 3. The level share is 14.41 %; L = 1 holds 62 % of the level class.",
  100000, 1, gen_beatty_s5 },
{ "beattys6", "A022840", "floor(n sqrt 6)", "Beatty",
  "floor(n sqrt 6). Gaps 2 and 3. The level share is 15.07 %; L = 1 holds 59 % of the level class.",
  100000, 1, gen_beatty_s6 },
{ "beattys7", "A022841", "floor(n sqrt 7)", "Beatty",
  "floor(n sqrt 7). Gaps 2 and 3. The level share is 15.61 %; L = 1 holds 57 % of the level class.",
  100000, 1, gen_beatty_s7 },
{ "beattys8", "A022842", "floor(n sqrt 8)", "Beatty",
  "floor(n sqrt 8). Gaps 2 and 3. The level share is 16.01 %; L = 1 holds 54 % of the level class.",
  100000, 1, gen_beatty_s8 },
{ "twodistsq", "A004431", "Sums of two distinct nonzero squares", "quadratic form",
  "x^2 + y^2 with 0 < x < y. The plane is that of all sums of two squares (15.98 % level against 15.95 %). The level share is 15.98 %; L = 1 holds 39 % of the level class.",
  100000, 1, gen_twodistsq },
{ "threeposq", "A000408", "Sums of three nonzero squares", "quadratic form",
  "x^2 + y^2 + z^2 with x, y, z >= 1. The level share is 8.56 %; L = 1 holds 99 % of the level class.",
  100000, 1, gen_threeposq },
{ "nottwosq", "A022544", "Not sums of two squares", "quadratic form",
  "Not a sum of two squares: some prime 3 mod 4 divides n to an odd power. The level share is 13.38 %; L = 1 holds 79 % of the level class; there are no ties.",
  100000, 1, gen_nottwosq },
{ "nontri", "A014132", "Non-triangular numbers", "complement",
  "Not a triangular number. Gaps 1, and 2 across each triangular number. The level share is 9.51 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_nontri },
{ "noncube", "A007412", "Non-cubes", "complement",
  "Not a cube. Gaps 1, and 2 across each cube: the naturals' plane. The level share is 9.60 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_noncube },
{ "nonfib", "A001690", "Non-Fibonacci numbers", "complement",
  "Not a Fibonacci number. Gaps 1 and 2: the naturals' plane. The level share is 9.60 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_nonfib },
{ "perfpow", "A001597", "Perfect powers", "powers",
  "1 and m^e with e >= 2. About sqrt(x) of them lie below x, so the sequence thins out like the squares, and 97.8 % of terms are forced level. The level share is 99.25 %; 97.8 % of terms are forced level (l <= d^2); there are no ties.",
  100000, 1, gen_perfpow },
{ "nonperfpow", "A007916", "Non-perfect powers", "powers",
  "Not a perfect power. Gaps 1 to 3: the naturals' plane. The level share is 9.65 %; L = 1 holds 99 % of the level class.",
  100000, 1, gen_nonperfpow },
{ "figfig", "A005228", "Hofstadter figure-figure sequence", "self-referential",
  "Hofstadter's figure-figure sequence: a(n+1) = a(n) + (the n-th number not in the sequence). a ~ n^2/2 and every term is forced level (100 %). Every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_figfig },
{ "non5", "A047201", "Numbers not divisible by 5", "forced divisor",
  "Not divisible by 5. Gaps 1 and 2. The level share is 10.38 %; L = 1 holds 85 % of the level class.",
  100000, 1, gen_non5 },
{ "mult4", "A008586", "Multiples of 4", "arithmetic progression",
  "4n: d = 4 and l = 4(n - 1). The level class is essentially {4n : n - 1 an odd prime}, on the line L = 4. The level share is 9.60 %; L = 4 holds 100 % of the level class.",
  100000, 0, gen_mult4 },
{ "mult5", "A008587", "Multiples of 5", "arithmetic progression",
  "5n: d = 5 and l = 5(n - 1). The level class lies on the line L = 5, as for the multiples of 3 on L = 3. The level share is 9.60 %; L = 5 holds 100 % of the level class.",
  100000, 0, gen_mult5 },
{ "mult6", "A008588", "Multiples of 6", "arithmetic progression",
  "6n: d = 6 and l = 6(n - 1). The level class lies on the line L = 6. The level share is 9.60 %; L = 6 holds 100 % of the level class; there are no ties.",
  100000, 0, gen_mult6 },
{ "mult7", "A008589", "Multiples of 7", "arithmetic progression",
  "7n: d = 7 and l = 7(n - 1). The level class lies on the line L = 7. The level share is 9.60 %; L = 7 holds 100 % of the level class.",
  100000, 0, gen_mult7 },
{ "n3p1", "A016777", "3n + 1", "arithmetic progression",
  "3n + 1: d = 3 and l = 3n - 2, prime to 3. The level share is about twice the multiples of 3's (9.60 %). The level share is 19.91 %; L = 1 holds 65 % of the level class.",
  100000, 0, gen_3np1 },
{ "n3p2", "A016789", "3n + 2", "arithmetic progression",
  "3n + 2: d = 3 and l = 3n - 1 = 2 (mod 3), never a square: no ties, by proof. The level share is 19.93 %; L = 1 holds 65 % of the level class.",
  100000, 0, gen_3np2 },
{ "n4p1", "A016813", "4n + 1", "arithmetic progression",
  "4n + 1: d = 4 and l = 4n - 3 = 1 (mod 4), odd. The level share is 23.15 %; L = 1 holds 73 % of the level class.",
  100000, 0, gen_4np1 },
{ "n4p3", "A004767", "4n + 3", "arithmetic progression",
  "4n + 3: d = 4 and l = 4n - 1 = 3 (mod 4), never a square: no ties, by proof. The level share is 23.15 %; L = 1 holds 73 % of the level class.",
  100000, 0, gen_4np3 },
{ "gon9", "A001106", "9-gonal (or enneagonal or nonagonal) numbers: a(n) = n*(7*n-5)/2", "polynomial",
  "Every gap is different, from 1 to 699,994; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_gon9 },
{ "gon10", "A001107", "10-gonal (or decagonal) numbers: a(n) = n*(4*n-3)", "polynomial",
  "Every gap is different, from 1 to 799,993; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_gon10 },
{ "gon11", "A051682", "11-gonal (or hendecagonal) numbers: a(n) = n*(9*n-7)/2", "polynomial",
  "Every gap is different, from 1 to 899,992; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_gon11 },
{ "gon12", "A051624", "12-gonal (or dodecagonal) numbers: a(n) = n*(5*n-4)", "polynomial",
  "Every gap is different, from 1 to 999,991; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_gon12 },
{ "gon13", "A051865", "13-gonal (or tridecagonal) numbers: a(n) = n*(11*n - 9)/2", "polynomial",
  "Every gap is different, from 1 to 1,099,990; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_gon13 },
{ "gon14", "A051866", "14-gonal (or tetradecagonal) numbers: a(n) = n*(6*n-5)", "polynomial",
  "Every gap is different, from 1 to 1,199,989; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_gon14 },
{ "gon15", "A051867", "15-gonal (or pentadecagonal) numbers: n*(13n-11)/2", "polynomial",
  "Every gap is different, from 1 to 1,299,988; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_gon15 },
{ "gon16", "A051868", "16-gonal (or hexadecagonal) numbers: a(n) = n*(7*n-6)", "polynomial",
  "Every gap is different, from 1 to 1,399,987; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_gon16 },
{ "cen9", "A060544", "Centered 9-gonal (also known as nonagonal or enneagonal) numbers. Every third triangular number, starting with a(1)=1", "polynomial",
  "Every gap is different, from 9 to 900,000; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_cen9 },
{ "cen10", "A062786", "Centered 10-gonal numbers", "polynomial",
  "Every gap is different, from 10 to 1,000,000; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_cen10 },
{ "cen11", "A069125", "a(n) = (11*n^2 - 11*n + 2)/2", "polynomial",
  "Every gap is different, from 11 to 1,100,000; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_cen11 },
{ "cen12", "A003154", "Centered 12-gonal numbers, or centered dodecagonal numbers: numbers of the form 6*k*(k-1) + 1", "polynomial",
  "Every gap is different, from 12 to 1,200,000; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_cen12 },
{ "cen7", "A069099", "Centered heptagonal numbers", "polynomial",
  "Every gap is different, from 7 to 700,000; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_cen7 },
{ "p4n13", "A001539", "a(n) = (4*n+1)*(4*n+3)", "polynomial",
  "Every gap is different, from 32 to 3,200,000; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_p4n13 },
{ "n2n2np1", "A002943", "a(n) = 2*n*(2*n+1)", "polynomial",
  "Every gap is different, from 6 to 799,998; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_2n2np1 },
{ "n4nm1", "A033991", "a(n) = n*(4*n-1)", "polynomial",
  "Every gap is different, from 3 to 799,995; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_n4nm1 },
{ "n4np1", "A007742", "a(n) = n*(4*n+1)", "polynomial",
  "Every gap is different, from 5 to 799,997; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_n4np1 },
{ "n2n2nm1", "A002939", "a(n) = 2*n*(2*n-1)", "polynomial",
  "Every gap is different, from 2 to 799,994; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_2n2nm1 },
{ "cpoly", "A002061", "Central polygonal numbers: a(n) = n^2 - n + 1", "polynomial",
  "Every gap is different, from 2 to 200,000; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_cpoly },
{ "n2sq2p1", "A058331", "a(n) = 2*n^2 + 1", "polynomial",
  "Every gap is different, from 2 to 399,998; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_2n2p1 },
{ "n2sq2m1", "A056220", "a(n) = 2*n^2 - 1", "polynomial",
  "Every gap is different, from 6 to 400,002; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_2n2m1 },
{ "n2np2", "A014206", "a(n) = n^2 + n + 2", "polynomial",
  "Every gap is different, from 2 to 200,000; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_n2np2 },
{ "nnp1sq", "A028387", "a(n) = n + (n+1)^2", "polynomial",
  "Every gap is different, from 4 to 200,002; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_nnp1sq },
{ "n2m2", "A008865", "a(n) = n^2 - 2", "polynomial",
  "Every gap is different, from 5 to 200,003; every decomposable term is forced level (l <= d^2).",
  100000, 2, gen_n2m2 },
{ "twopyr", "A006331", "a(n) = n*(n+1)*(2*n+1)/3", "polynomial",
  "Every gap is different, from 2 to 20,000,000,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_2pyr },
{ "pyr7", "A002413", "Heptagonal (or 7-gonal) pyramidal numbers: a(n) = n*(n+1)*(5*n-2)/6", "polynomial",
  "Every gap is different, from 1 to 24,999,850,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_pyr7 },
{ "pyr8", "A002414", "Octagonal pyramidal numbers: a(n) = n*(n+1)*(2*n-1)/2", "polynomial",
  "Every gap is different, from 8 to 30,000,400,001; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_pyr8 },
{ "pyr9", "A007584", "9-gonal (or enneagonal) pyramidal numbers: a(n) = n*(n+1)*(7*n-4)/6", "polynomial",
  "Every gap is different, from 1 to 34,999,750,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_pyr9 },
{ "pyr10", "A007585", "10-gonal (or decagonal) pyramidal numbers: a(n) = n*(n + 1)*(8*n - 5)/6", "polynomial",
  "Every gap is different, from 1 to 39,999,700,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_pyr10 },
{ "centet", "A005894", "Centered tetrahedral numbers", "polynomial",
  "Every gap is different, from 4 to 20,000,000,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_centet },
{ "dodeca", "A006566", "Dodecahedral numbers: a(n) = n*(3*n - 1)*(3*n - 2)/2", "polynomial",
  "Every gap is different, from 1 to 134,997,750,010; every decomposable term is forced level (l <= d^2); 8 terms do not decompose.",
  100000, 0, gen_dodeca },
{ "icosa", "A006564", "Icosahedral numbers: a(n) = n*(5*n^2 - 5*n + 2)/2", "polynomial",
  "Every gap is different, from 11 to 75,000,250,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_icosa },
{ "rhdod", "A005917", "Rhombic dodecahedral numbers: a(n) = n^4 - (n - 1)^4", "polynomial",
  "Every gap is different, from 14 to 120,000,000,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_rhdod },
{ "cubep1", "A001093", "a(n) = n^3 + 1", "polynomial",
  "Every gap is different, from 1 to 29,999,700,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_cubep1 },
{ "cubem1", "A068601", "a(n) = n^3 - 1", "polynomial",
  "Every gap is different, from 7 to 30,000,300,001; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_cubem1 },
{ "cubepn", "A034262", "a(n) = n^3 + n", "polynomial",
  "Every gap is different, from 2 to 29,999,700,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_cubepn },
{ "hexprism", "A005915", "Hexagonal prism numbers: a(n) = (n + 1)*(3*n^2 + 3*n + 1)", "polynomial",
  "Every gap is different, from 13 to 90,000,300,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_hexprism },
{ "cenicosa", "A005902", "Centered icosahedral (or cuboctahedral) numbers, also crystal ball sequence for f.c.c. lattice", "polynomial",
  "Every gap is different, from 12 to 100,000,000,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_cenicosa },
{ "trunctet", "A005906", "Truncated tetrahedral numbers: a(n) = (1/6)*(n+1)*(23*n^2 + 19*n + 6)", "polynomial",
  "Every gap is different, from 15 to 115,000,250,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_trunctet },
{ "secpentodd", "A033568", "Second pentagonal numbers with odd index: a(n) = (2*n-1)*(3*n-1)", "polynomial",
  "Every gap is different, from 1 to 1,199,989; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_secpentodd },
{ "beattye", "A022843", "Beatty sequence for e: a(n) = floor(n*e)", "Beatty",
  "The gaps are 2 and 3; the level share is 15.74 %; L = 1 holds 56 % of the level class.",
  100000, 0, gen_beattye },
{ "beattypi", "A022844", "a(n) = floor(n*Pi)", "Beatty",
  "The gaps are 3 and 4; the level share is 16.83 %; L = 1 holds 52 % of the level class.",
  100000, 0, gen_beattypi },
{ "beattyee1", "A054385", "Beatty sequence for e/(e-1); complement of A022843", "Beatty",
  "The gaps are 1 and 2; the level share is 12.08 %; L = 1 holds 77 % of the level class.",
  100000, 1, gen_beattyee1 },
{ "beatty2s3", "A003512", "A Beatty sequence: floor(n*(sqrt(3) + 2))", "Beatty",
  "The gaps are 3 and 4; the level share is 18.04 %; L = 1 holds 47 % of the level class.",
  100000, 1, gen_beatty2s3 },
{ "halfs2", "A001953", "a(n) = floor((n + 1/2) * sqrt(2))", "Beatty",
  "The gaps are 1 and 2; the level share is 11.32 %; L = 1 holds 82 % of the level class.",
  100000, 0, gen_halfs2 },
{ "halfs22", "A001954", "a(n) = floor((n+1/2)*(2+sqrt(2))); winning positions in the 2-Wythoff game", "Beatty",
  "The gaps are 3 and 4; the level share is 17.27 %; L = 1 holds 49 % of the level class.",
  100000, 0, gen_halfs22 },
{ "beattys55", "A003231", "a(n) = floor(n*(sqrt(5)+5)/2)", "Beatty",
  "The gaps are 3 and 4; the level share is 17.71 %; L = 1 holds 48 % of the level class.",
  100000, 1, gen_beattys55 },
{ "zeckeven", "A022342", "Integers with \"even\" Zeckendorf expansions (do not end with ... + F_2 = ... + 1) (the Fibonacci-even numbers); also, apart from first term, a(n) = Fibonacci successor to n-1", "Beatty",
  "The gaps are 1 and 2; the level share is 12.29 %; L = 1 holds 76 % of the level class.",
  100000, 1, gen_zeckeven },
{ "wythaa", "A003622", "The Wythoff compound sequence AA: a(n) = floor(n*phi^2) - 1, where phi = (1+sqrt(5))/2", "Beatty",
  "The gaps are 2 and 3; the level share is 15.49 %; L = 1 holds 57 % of the level class.",
  100000, 1, gen_wythAA },
{ "near2", "A022846", "Nearest integer to n*sqrt(2)", "Beatty",
  "The gaps are 1 and 2; the level share is 11.32 %; L = 1 holds 82 % of the level class.",
  100000, 0, gen_near2 },
{ "near3", "A022847", "Integer nearest n*sqrt(3)", "Beatty",
  "The gaps are 1 and 2; the level share is 12.67 %; L = 1 holds 72 % of the level class.",
  100000, 0, gen_near3 },
{ "near5", "A022848", "Integer nearest nx, where x = sqrt(5)", "Beatty",
  "The gaps are 2 and 3; the level share is 14.47 %; L = 1 holds 62 % of the level class.",
  100000, 0, gen_near5 },
{ "beattys5m1", "A001961", "A Beatty sequence: floor(n * (sqrt(5) - 1))", "Beatty",
  "The gaps are 1 and 2; the level share is 10.65 %; L = 1 holds 89 % of the level class.",
  100000, 1, gen_beattys5m1 },
{ "p1mod12", "A068228", "Primes congruent to 1 (mod 12)", "primes",
  "40 different gaps occur, from 12 to 576; the level share is 40.17 %; L = 1 holds 42 % of the level class.",
  100000, 1, gen_p1mod12 },
{ "p5mod12", "A040117", "Primes congruent to 5 (mod 12). Also primes p such that x^4 = 9 has no solution mod p", "primes",
  "41 different gaps occur, from 12 to 576; the level share is 40.26 %; L = 1 holds 42 % of the level class; there are no ties.",
  100000, 1, gen_p5mod12 },
{ "p7mod12", "A068229", "Primes congruent to 7 (mod 12)", "primes",
  "40 different gaps occur, from 12 to 492; the level share is 40.19 %; L = 1 holds 42 % of the level class; there are no ties.",
  100000, 1, gen_p7mod12 },
{ "p11mod12", "A068231", "Primes congruent to 11 mod 12", "primes",
  "40 different gaps occur, from 12 to 540; the level share is 40.03 %; L = 1 holds 42 % of the level class; there are no ties.",
  100000, 1, gen_p11mod12 },
{ "p14mod5", "A045468", "Primes congruent to {1, 4} mod 5", "primes",
  "68 different gaps occur, from 2 to 252; the level share is 31.78 %; L = 1 holds 32 % of the level class.",
  100000, 1, gen_p14mod5 },
{ "p23mod5", "A003631", "Primes congruent to 2 or 3 modulo 5", "primes",
  "69 different gaps occur, from 1 to 294; the level share is 32.17 %; L = 3 holds 33 % of the level class.",
  100000, 1, gen_p23mod5 },
{ "p127mod8", "A038873", "Primes p such that 2 is a square mod p; or, primes congruent to {1, 2, 7} mod 8", "primes",
  "85 different gaps occur, from 2 to 262; the level share is 28.26 %.",
  100000, 1, gen_p127mod8 },
{ "ppm3mod8", "A003629", "Primes p == +- 3 (mod 8), or, primes p such that 2 is not a square mod p", "primes",
  "82 different gaps occur, from 2 to 232; the level share is 28.01 %.",
  100000, 1, gen_ppm3mod8 },
{ "ppm1mod8", "A001132", "Primes == +-1 (mod 8)", "primes",
  "84 different gaps occur, from 2 to 262; the level share is 28.26 %.",
  100000, 1, gen_ppm1mod8 },
{ "p13mod8", "A033200", "Primes congruent to {1, 3} (mod 8); or, odd primes of form x^2 + 2*y^2", "primes",
  "83 different gaps occur, from 2 to 226; the level share is 28.46 %.",
  100000, 1, gen_p13mod8 },
{ "pfirst1", "A045707", "Primes with first digit 1", "primes",
  "69 different gaps occur, from 2 to 8,000,026; the level share is 21.45 %; L = 1 holds 32 % of the level class; 6 terms do not decompose.",
  100000, 1, gen_pfirst1 },
{ "p2pp3", "A023204", "Primes p such that 2*p + 3 is also prime", "primes",
  "342 different gaps occur, from 2 to 910; the level share is 36.43 %.",
  100000, 1, gen_p2pp3 },
{ "p3pp2", "A023208", "Primes p such that 3*p + 2 is also prime", "primes",
  "357 different gaps occur, from 2 to 1,132; the level share is 36.36 %.",
  100000, 1, gen_p3pp2 },
{ "p4pp1", "A023212", "Primes p such that 4*p+1 is also prime", "primes",
  "278 different gaps occur, from 4 to 2,184; the level share is 47.59 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_p4pp1 },
{ "chen", "A109611", "Chen primes: primes p such that p + 2 is either a prime or a semiprime", "primes",
  "156 different gaps occur, from 1 to 418; the level share is 34.50 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_chen },
{ "triples", "A007529", "Prime triples: p; p+2 or p+4; p+6 all prime", "primes",
  "1,837 different gaps occur, from 2 to 11,398; the level share is 49.12 %; 1.6 % of terms are forced level (l <= d^2).",
  100000, 1, gen_triples },
{ "root2", "A001122", "Primes with primitive root 2", "primes",
  "114 different gaps occur, from 2 to 410; the level share is 30.54 %.",
  100000, 1, gen_root2 },
{ "root10", "A001913", "Full reptend primes: primes with primitive root 10", "primes",
  "147 different gaps occur, from 2 to 408; the level share is 30.72 %.",
  100000, 1, gen_root10 },
{ "up6", "A031925", "Upper prime of a difference of 6 between consecutive primes", "primes",
  "402 different gaps occur, from 6 to 1,288; the level share is 46.88 %; L = 1 holds 41 % of the level class.",
  100000, 1, gen_up6 },
{ "lo8", "A031926", "Lower prime of a difference of 8 between consecutive primes", "primes",
  "342 different gaps occur, from 12 to 2,988; the level share is 49.17 %; L = 1 holds 34 % of the level class; there are no ties.",
  100000, 1, gen_lo8 },
{ "lo10", "A031928", "Lower prime of a difference of 10 between consecutive primes", "primes",
  "269 different gaps occur, from 12 to 2,742; the level share is 47.71 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_lo10 },
{ "lo12", "A031930", "Lower prime of a difference of 12 between consecutive primes", "primes",
  "580 different gaps occur, from 12 to 1,830; the level share is 40.46 %.",
  100000, 1, gen_lo12 },
{ "ppp", "A038580", "Primes with indices that are primes with prime indices", "primes",
  "7,843 different gaps occur, from 6 to 36,468; the level share is 64.56 %; 6.5 % of terms are forced level (l <= d^2).",
  100000, 1, gen_ppp },
{ "pnonprimeidx", "A007821", "Primes p such that pi(p) is not prime", "primes",
  "60 different gaps occur, from 2 to 132; the level share is 23.82 %; L = 1 holds 32 % of the level class.",
  100000, 1, gen_pnonprimeidx },
{ "pform712", "A106856", "Primes of the form x^2 + xy + 2y^2, with x and y nonnegative", "primes",
  "142 different gaps occur, from 2 to 378; the level share is 27.75 %.",
  100000, 1, gen_pform712 },
{ "pform15", "A033205", "Primes of form x^2 + 5*y^2", "primes",
  "73 different gaps occur, from 8 to 660; the level share is 36.72 %.",
  100000, 1, gen_pform15 },
{ "fourp", "A001749", "Primes multiplied by 4", "primes",
  "54 different gaps occur, from 4 to 456; the level share is 23.02 %; L = 4 holds 32 % of the level class.",
  100000, 1, gen_fourp },
{ "ppm1", "A036689", "Product of a prime and the previous number", "primes",
  "99,993 different gaps occur, from 4 to 238,596,576; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_ppm1 },
{ "ppp1", "A036690", "Product of a prime and the following number", "primes",
  "99,995 different gaps occur, from 6 to 238,596,768; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_ppp1 },
{ "psqm2", "A049001", "a(n) = prime(n)^2 - 2", "primes",
  "94,934 different gaps occur, from 5 to 238,596,672; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_psqm2 },
{ "psum3", "A034961", "Sums of three consecutive primes", "primes",
  "73 different gaps occur, from 5 to 152; the level share is 41.68 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_psum3 },
{ "psum4", "A034963", "Sums of four consecutive primes", "primes",
  "82 different gaps occur, from 9 to 180; the level share is 31.42 %.",
  100000, 1, gen_psum4 },
{ "pavg", "A024675", "Average of two consecutive odd primes", "primes",
  "63 different gaps occur, from 2 to 69; the level share is 24.89 %.",
  100000, 1, gen_pavg },
{ "pminus2", "A040976", "a(n) = prime(n) - 2", "primes",
  "54 different gaps occur, from 1 to 114; the level share is 30.62 %; L = 1 holds 60 % of the level class.",
  100000, 1, gen_pminus2 },
{ "phalfm", "A005097", "(Odd primes - 1)/2", "primes",
  "53 different gaps occur, from 1 to 57; the level share is 22.38 %; L = 1 holds 45 % of the level class.",
  100000, 1, gen_phalfm },
{ "phalfp", "A006254", "Numbers k such that 2k-1 is prime", "primes",
  "53 different gaps occur, from 1 to 57; the level share is 22.30 %; L = 1 holds 43 % of the level class.",
  100000, 1, gen_phalfp },
{ "twin6m", "A002822", "Numbers m such that 6m-1, 6m+1 are twin primes", "primes",
  "246 different gaps occur, from 1 to 365; the level share is 31.70 %.",
  100000, 1, gen_twin6m },
{ "sum2p", "A014091", "Numbers that are the sum of 2 primes", "primes",
  "The gaps are 1 and 2; the level share is 22.64 %; L = 1 holds 68 % of the level class.",
  100000, 1, gen_sum2p },
{ "notsum2p", "A014092", "Numbers that are not the sum of 2 primes", "primes",
  "5 different gaps occur, from 1 to 8; the level share is 4.22 %; L = 1 holds 59 % of the level class.",
  100000, 1, gen_notsum2p },
{ "almost6", "A046306", "Numbers that are divisible by exactly 6 primes with multiplicity", "multiplicative",
  "114 different gaps occur, from 1 to 140; the level share is 18.97 %.",
  100000, 1, gen_almost6 },
{ "almost7", "A046308", "Numbers that are divisible by exactly 7 primes counting multiplicity", "multiplicative",
  "190 different gaps occur, from 1 to 280; the level share is 19.73 %.",
  100000, 1, gen_almost7 },
{ "almost8", "A046310", "Numbers that are divisible by exactly 8 primes counting multiplicity", "multiplicative",
  "327 different gaps occur, from 1 to 560; the level share is 20.51 %.",
  100000, 1, gen_almost8 },
{ "oddsemi", "A046315", "Odd semiprimes: odd numbers divisible by exactly 2 primes (counted with multiplicity)", "multiplicative",
  "25 different gaps occur, from 2 to 50; the level share is 21.00 %; L = 1 holds 58 % of the level class.",
  100000, 1, gen_oddsemi },
{ "oddpq", "A046388", "Odd numbers of the form p*q where p and q are distinct primes", "multiplicative",
  "25 different gaps occur, from 2 to 50; the level share is 21.01 %; L = 1 holds 58 % of the level class.",
  100000, 1, gen_oddpq },
{ "atleast3", "A000977", "Numbers that are divisible by at least three different primes", "multiplicative",
  "12 different gaps occur, from 1 to 18; the level share is 13.76 %; L = 1 holds 52 % of the level class.",
  100000, 1, gen_atleast3 },
{ "bitexp", "A000379", "Numbers where total number of 1-bits in the exponents of their prime factorization is even; a 2-way classification of integers: complement of A000028", "multiplicative",
  "17 different gaps occur, from 1 to 17; the level share is 12.70 %; L = 1 holds 71 % of the level class.",
  100000, 1, gen_bitexp },
{ "biqfree", "A046100", "Biquadratefree numbers: numbers that are not divisible by any 4th power greater than 1", "multiplicative",
  "The gaps are 1, 2, 3 and 4; the level share is 9.93 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_biqfree },
{ "biqful", "A046101", "Biquadrateful numbers", "multiplicative",
  "16 different gaps occur, from 1 to 16; the level share is 13.70 %; L = 16 holds 46 % of the level class.",
  100000, 1, gen_biqful },
{ "sqfpair", "A007674", "Numbers m such that m and m+1 are squarefree", "multiplicative",
  "15 different gaps occur, from 1 to 23; the level share is 14.63 %; L = 1 holds 60 % of the level class.",
  100000, 1, gen_sqfpair },
{ "cyclic", "A003277", "Cyclic numbers: k such that k and phi(k) are relatively prime; also k such that there is just one group of order k, i.e., A000001(k) = 1", "divisor functions",
  "10 different gaps occur, from 1 to 18; the level share is 15.14 %; L = 1 holds 50 % of the level class.",
  100000, 1, gen_cyclic },
{ "rough11", "A008364", "11-rough numbers: not divisible by 2, 3, 5 or 7", "forced divisor",
  "5 different gaps occur, from 2 to 10; the level share is 14.72 %; L = 3 holds 38 % of the level class.",
  100000, 1, gen_rough11 },
{ "rough13", "A008365", "13-rough numbers: positive integers that have no prime factors less than 13", "forced divisor",
  "7 different gaps occur, from 2 to 14; the level share is 15.51 %; L = 3 holds 37 % of the level class.",
  100000, 1, gen_rough13 },
{ "rough17", "A008366", "Smallest prime factor is >= 17", "forced divisor",
  "10 different gaps occur, from 2 to 22; the level share is 16.19 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_rough17 },
{ "smooth11", "A051038", "11-smooth numbers: numbers whose prime divisors are all <= 11", "smooth",
  "82,783 different gaps occur, from 1 to 103,527,495,225; the level share is 98.66 %; 96.6 % of terms are forced level (l <= d^2).",
  100000, 1, gen_smooth11 },
{ "smooth13", "A080197", "13-smooth numbers: numbers whose prime divisors are all <= 13", "smooth",
  "62,590 different gaps occur, from 1 to 115,531,416; the level share is 93.08 %; 87.3 % of terms are forced level (l <= d^2).",
  100000, 1, gen_smooth13 },
{ "smooth17", "A080681", "17-smooth numbers: numbers whose prime divisors are all <= 17", "smooth",
  "39,268 different gaps occur, from 1 to 3,304,800; the level share is 79.21 %; 62.4 % of terms are forced level (l <= d^2).",
  100000, 1, gen_smooth17 },
{ "blum", "A016105", "Blum integers: numbers of the form p * q where p and q are distinct primes congruent to 3 (mod 4)", "multiplicative",
  "44 different gaps occur, from 4 to 244; the level share is 29.94 %; L = 1 holds 43 % of the level class.",
  100000, 1, gen_blum },
{ "tau4", "A030513", "Numbers with 4 divisors", "divisor functions",
  "38 different gaps occur, from 1 to 47; the level share is 15.72 %; L = 1 holds 61 % of the level class.",
  100000, 1, gen_tau4 },
{ "tau6", "A030515", "Numbers with exactly 6 divisors", "divisor functions",
  "189 different gaps occur, from 1 to 240; the level share is 27.45 %.",
  100000, 1, gen_tau6 },
{ "tau8", "A030626", "Numbers with exactly 8 divisors", "divisor functions",
  "37 different gaps occur, from 1 to 44; the level share is 16.53 %; L = 1 holds 62 % of the level class.",
  100000, 1, gen_tau8 },
{ "tau10", "A030628", "1 together with numbers of the form p*q^4 and p^9, where p and q are distinct primes", "divisor functions",
  "661 different gaps occur, from 1 to 1,600; the level share is 29.84 %.",
  100000, 1, gen_tau10 },
{ "tau12", "A030630", "Numbers with 12 divisors", "divisor functions",
  "91 different gaps occur, from 1 to 113; the level share is 21.35 %.",
  100000, 1, gen_tau12 },
{ "sqfsq", "A062503", "Squarefree numbers squared", "multiplicative",
  "95,658 different gaps occur, from 3 to 2,132,781; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_sqfsq },
{ "equidig", "A046758", "Equidigital numbers", "digit rule",
  "47 different gaps occur, from 1 to 58; the level share is 16.23 %; L = 1 holds 60 % of the level class.",
  100000, 1, gen_equidig },
{ "wasteful", "A046760", "Wasteful numbers", "digit rule",
  "8 different gaps occur, from 1 to 8; the level share is 11.49 %; L = 1 holds 69 % of the level class.",
  100000, 1, gen_wasteful },
{ "smith", "A006753", "Smith (or joke) numbers: composite numbers k such that sum of digits of k = sum of digits of prime factors of k (counted with multiplicity)", "digit rule",
  "429 different gaps occur, from 1 to 879; the level share is 29.16 %.",
  100000, 1, gen_smith },
{ "hoax", "A019506", "Hoax numbers: composite numbers whose digit-sum equals the sum of the digit-sums of its distinct prime factors", "digit rule",
  "304 different gaps occur, from 1 to 546; the level share is 27.46 %.",
  100000, 1, gen_hoax },
{ "hati", "A036668", "Hati numbers: of form 2^i*3^j*k, i+j even, (k,6)=1", "multiplicative",
  "The gaps are 1, 2, 3 and 4; the level share is 13.08 %; L = 1 holds 75 % of the level class.",
  100000, 1, gen_hati },
{ "evennontot", "A005277", "Nontotients: even numbers k such that phi(m) = k has no solution", "divisor functions",
  "8 different gaps occur, from 2 to 16; the level share is 10.73 %; L = 2 holds 95 % of the level class; there are no ties.",
  100000, 1, gen_evennontot },
{ "brilliant", "A078972", "Brilliant numbers: semiprimes (products of two primes, A001358) whose prime factors have the same number of decimal digits", "multiplicative",
  "619 different gaps occur, from 1 to 24,072; the level share is 35.07 %.",
  100000, 1, gen_brilliant },
{ "sqfcomp", "A120944", "Composite squarefree numbers", "multiplicative",
  "10 different gaps occur, from 1 to 10; the level share is 11.93 %; L = 1 holds 89 % of the level class.",
  100000, 1, gen_sqfcomp },
{ "sq2sq", "A028982", "Squares and twice squares", "quadratic form",
  "34,684 different gaps occur, from 1 to 117,155; the level share is 66.54 %; 58.6 % of terms are forced level (l <= d^2); there are no ties.",
  100000, 1, gen_sq2sq },
{ "sigmaeven", "A028983", "Numbers whose sum of divisors is even", "divisor functions",
  "The gaps are 1, 2 and 3; the level share is 9.63 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_sigmaeven },
{ "exactly3sq", "A000419", "Numbers that are the sum of 3 but no fewer nonzero squares", "quadratic form",
  "5 different gaps occur, from 1 to 5; the level share is 12.59 %; L = 1 holds 71 % of the level class.",
  100000, 1, gen_exactly3sq },
{ "form16", "A002481", "Numbers of form x^2 + 6y^2", "quadratic form",
  "44 different gaps occur, from 1 to 51; the level share is 22.80 %; L = 1 holds 48 % of the level class.",
  100000, 1, gen_form16 },
{ "form23", "A002480", "Numbers of the form 2x^2 + 3y^2", "quadratic form",
  "44 different gaps occur, from 1 to 45; the level share is 22.29 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_form23 },
{ "form112", "A000401", "Numbers of form x^2 + y^2 + 2*z^2", "quadratic form",
  "The gaps are 1 and 2; the level share is 10.38 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_form112 },
{ "form15", "A020669", "Numbers of form x^2 + 5 y^2", "quadratic form",
  "43 different gaps occur, from 1 to 51; the level share is 22.90 %; L = 1 holds 41 % of the level class.",
  100000, 1, gen_form15 },
{ "twocubes", "A003325", "Numbers that are the sum of 2 positive cubes", "quadratic form",
  "5,379 different gaps occur, from 1 to 16,955; the level share is 50.13 %; 1.3 % of terms are forced level (l <= d^2).",
  100000, 1, gen_twocubes },
{ "twotri", "A020756", "Numbers that are the sum of two triangular numbers", "quadratic form",
  "17 different gaps occur, from 1 to 18; the level share is 14.83 %; L = 1 holds 64 % of the level class.",
  100000, 1, gen_twotri },
{ "nottwotri", "A020757", "Numbers that are not the sum of two triangular numbers", "quadratic form",
  "6 different gaps occur, from 1 to 6; the level share is 12.27 %; L = 1 holds 75 % of the level class.",
  100000, 1, gen_nottwotri },
{ "has2", "A011532", "Numbers that contain a 2", "digit rule",
  "The gaps are 1, 3, 8 and 10; the level share is 11.40 %; L = 1 holds 68 % of the level class.",
  100000, 1, gen_has2 },
{ "has3", "A011533", "Numbers that contain a 3", "digit rule",
  "The gaps are 1, 4, 7 and 10; the level share is 12.68 %; L = 1 holds 82 % of the level class.",
  100000, 1, gen_has3 },
{ "has4", "A011534", "Numbers that contain a 4", "digit rule",
  "The gaps are 1, 5, 6 and 10; the level share is 11.47 %; L = 1 holds 64 % of the level class.",
  100000, 1, gen_has4 },
{ "has5", "A011535", "Numbers that contain a 5", "digit rule",
  "The gaps are 1, 5, 6 and 10; the level share is 10.78 %; L = 1 holds 72 % of the level class.",
  100000, 1, gen_has5 },
{ "has6", "A011536", "Numbers that contain a 6", "digit rule",
  "The gaps are 1, 4, 7 and 10; the level share is 12.01 %; L = 1 holds 61 % of the level class.",
  100000, 1, gen_has6 },
{ "has7", "A011537", "Numbers that contain at least one 7", "digit rule",
  "The gaps are 1, 3, 8 and 10; the level share is 13.62 %; L = 1 holds 81 % of the level class.",
  100000, 1, gen_has7 },
{ "has8", "A011538", "Numbers that contain an 8", "digit rule",
  "The gaps are 1, 2, 9 and 10; the level share is 11.16 %; L = 1 holds 65 % of the level class.",
  100000, 1, gen_has8 },
{ "no3", "A052405", "Numbers without 3 as a digit", "digit rule",
  "6 different gaps occur, from 1 to 10,001; the level share is 8.15 %; L = 1 holds 94 % of the level class.",
  100000, 1, gen_no3 },
{ "no4", "A052406", "Numbers without 4 as a digit", "digit rule",
  "6 different gaps occur, from 1 to 10,001; the level share is 10.51 %; L = 1 holds 97 % of the level class.",
  100000, 1, gen_no4 },
{ "no5", "A052413", "Numbers without 5 as a digit", "digit rule",
  "6 different gaps occur, from 1 to 10,001; the level share is 9.64 %; L = 1 holds 79 % of the level class.",
  100000, 1, gen_no5 },
{ "no6", "A052414", "Numbers without 6 as a digit", "digit rule",
  "6 different gaps occur, from 1 to 10,001; the level share is 13.13 %; L = 1 holds 97 % of the level class.",
  100000, 1, gen_no6 },
{ "no7", "A052419", "Numbers without 7 as a digit", "digit rule",
  "6 different gaps occur, from 1 to 10,001; the level share is 11.92 %; L = 1 holds 85 % of the level class.",
  100000, 1, gen_no7 },
{ "no8", "A052421", "Numbers without 8 as a digit", "digit rule",
  "6 different gaps occur, from 1 to 10,001; the level share is 8.01 %; L = 1 holds 96 % of the level class.",
  100000, 1, gen_no8 },
{ "base7", "A007093", "Numbers in base 7", "digit rule",
  "6 different gaps occur, from 1 to 33,334; the level share is 11.77 %; L = 1 holds 76 % of the level class; 8 terms do not decompose.",
  100000, 0, gen_base7 },
{ "base8", "A007094", "Numbers in base 8", "digit rule",
  "6 different gaps occur, from 1 to 22,223; the level share is 10.02 %; L = 1 holds 82 % of the level class.",
  100000, 0, gen_base8 },
{ "base9", "A007095", "Numbers in base 9", "digit rule",
  "6 different gaps occur, from 1 to 11,112; the level share is 9.22 %; L = 1 holds 83 % of the level class.",
  100000, 0, gen_base9 },
{ "dsodd", "A054684", "Numbers whose sum of digits is odd", "digit rule",
  "The gaps are 1, 2 and 3; the level share is 11.04 %; L = 1 holds 66 % of the level class.",
  100000, 1, gen_dsodd },
{ "dsfib", "A028840", "Numbers k such that sum of digits of k is a Fibonacci number", "digit rule",
  "49 different gaps occur, from 1 to 300; the level share is 24.05 %.",
  100000, 1, gen_dsfib },
{ "pal3", "A014190", "Palindromes in base 3 (written in base 10)", "digit rule",
  "21 different gaps occur, from 1 to 78,732; the level share is 84.19 %; 32.4 % of terms are forced level (l <= d^2).",
  100000, 1, gen_pal3 },
{ "pal4", "A014192", "Palindromes in base 4 (written in base 10)", "digit rule",
  "17 different gaps occur, from 1 to 81,920; the level share is 83.21 %; 28.2 % of terms are forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_pal4 },
{ "pal5", "A029952", "Palindromic in base 5", "digit rule",
  "15 different gaps occur, from 1 to 93,750; the level share is 78.45 %; 18.9 % of terms are forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_pal5 },
{ "pal6", "A029953", "Palindromic in base 6", "digit rule",
  "14 different gaps occur, from 1 to 54,432; the level share is 83.24 %; 40.3 % of terms are forced level (l <= d^2).",
  100000, 1, gen_pal6 },
{ "no3b4", "A023717", "Numbers with no 3's in base-4 expansion", "digit rule",
  "11 different gaps occur, from 1 to 349,526; the level share is 7.81 %; L = 1 holds 66 % of the level class; 12 terms do not decompose.",
  100000, 0, gen_no3b4 },
{ "altpar", "A030141", "Numbers in which parity of the decimal digits alternates", "digit rule",
  "13 different gaps occur, from 1 to 202,021; the level share is 14.28 %; L = 1 holds 49 % of the level class.",
  100000, 1, gen_altpar },
{ "noprimedig", "A084984", "Numbers containing no prime digits", "digit rule",
  "14 different gaps occur, from 1 to 2,000,001; the level share is 15.18 %; L = 2 holds 49 % of the level class; 9 terms do not decompose.",
  100000, 1, gen_noprimedig },
{ "has0b3", "A081605", "Numbers having at least one 0 in their ternary representation", "digit rule",
  "The gaps are 1 and 3; the level share is 9.55 %; L = 1 holds 99 % of the level class.",
  100000, 1, gen_has0b3 },
{ "has11", "A004780", "Binary expansion contains 2 adjacent 1's", "binary rule",
  "The gaps are 1, 3 and 4; the level share is 9.71 %; L = 1 holds 99 % of the level class.",
  100000, 1, gen_has11 },
{ "no00", "A003754", "Numbers with no adjacent 0's in binary expansion", "binary rule",
  "22 different gaps occur, from 1 to 1,398,102; the level share is 11.57 %; L = 1 holds 76 % of the level class.",
  100000, 1, gen_no00 },
{ "rsneg", "A022155", "Values of n at which Golay-Rudin-Shapiro sequence A020985 is negative", "binary rule",
  "5 different gaps occur, from 1 to 5; the level share is 11.47 %; L = 1 holds 80 % of the level class.",
  100000, 1, gen_rsneg },
{ "rspos", "A203463", "Where Golay-Rudin-Shapiro sequence A020985 is positive", "binary rule",
  "5 different gaps occur, from 1 to 5; the level share is 11.09 %; L = 1 holds 80 % of the level class.",
  100000, 1, gen_rspos },
{ "pf1", "A091072", "Positive numbers k such that the Kronecker Symbol (-1 / k) > 0", "binary rule",
  "The gaps are 1, 2, 3 and 4; the level share is 13.06 %; L = 1 holds 69 % of the level class.",
  100000, 1, gen_pf1 },
{ "pf3", "A091067", "Numbers whose odd part is of the form 4*k+3", "binary rule",
  "The gaps are 1, 2, 3 and 4; the level share is 13.08 %; L = 1 holds 69 % of the level class.",
  100000, 1, gen_pf3 },
{ "binself", "A010061", "Binary self or Colombian numbers: numbers that cannot be expressed as the sum of distinct terms of the form 2^k+1 (k>=0), or equivalently, numbers not of form m + sum of binary digits of m", "binary rule",
  "12 different gaps occur, from 2 to 21; the level share is 17.09 %; L = 1 holds 50 % of the level class.",
  100000, 1, gen_binself },
{ "zeckodd", "A020899", "Numbers k with an odd number of terms in their Zeckendorf representation (write k as a sum of non-consecutive distinct Fibonacci numbers)", "digit rule",
  "The gaps are 1, 2, 3 and 4; the level share is 12.98 %; L = 1 holds 69 % of the level class.",
  100000, 1, gen_zeckodd },
{ "ludic", "A003309", "Ludic numbers: apply the same sieve as Eratosthenes, but cross off every k-th remaining number", "sieve",
  "59 different gaps occur, from 1 to 120; the level share is 27.21 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_ludic },
{ "evenlucky", "A045954", "Even-Lucky-Numbers: generated by a sieve process like that for Lucky numbers but starting with even numbers", "sieve",
  "56 different gaps occur, from 2 to 128; the level share is 20.56 %; L = 2 holds 50 % of the level class.",
  100000, 1, gen_evenlucky },
{ "figcomp", "A030124", "Complement (and also first differences) of Hofstadter's sequence A005228", "self-referential",
  "The gaps are 1 and 2; the level share is 9.59 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_figcomp },
{ "ulam13", "A002859", "a(1) = 1, a(2) = 3; for n >= 3, a(n) is smallest number that is uniquely of the form a(j) + a(k) with 1 <= j < k < n", "self-referential",
  "49 different gaps occur, from 1 to 98; the level share is 20.60 %; L = 1 holds 39 % of the level class.",
  100000, 1, gen_ulam13 },
{ "klarner", "A002977", "Klarner-Rado sequence: a(1) = 1; subsequent terms are defined by the rule that if m is present so are 2m+1 and 3m+1", "self-referential",
  "522 different gaps occur, from 1 to 15,737; the level share is 23.55 %; L = 1 holds 67 % of the level class.",
  100000, 1, gen_klarner },
{ "addomega", "A094222", "a(n+1) = a(n) + (number of distinct prime factors of a(n)) for n>1; a(1)=1, a(2)=2", "self-referential",
  "6 different gaps occur, from 1 to 6; the level share is 12.30 %; L = 1 holds 54 % of the level class.",
  100000, 1, gen_addomega },
{ "barriers", "A005236", "Barriers for omega(n): numbers n such that, for all m < n, m + omega(m) <= n", "self-referential",
  "590 different gaps occur, from 1 to 4,566; the level share is 33.13 %.",
  100000, 1, gen_barriers },
{ "bitsum", "A000788", "Total number of 1's in binary expansions of 0, ..., n", "summatory",
  "16 different gaps occur, from 1 to 16; the level share is 24.12 %; L = 1 holds 34 % of the level class.",
  100000, 0, gen_bitsum },
{ "digsum", "A037123", "a(n) = a(n-1) + sum of digits of n", "summatory",
  "45 different gaps occur, from 1 to 45; the level share is 28.31 %; there are no ties.",
  100000, 0, gen_digsum },
{ "omegasum", "A013939", "Partial sums of sequence A001221 (number of distinct primes dividing n)", "summatory",
  "6 different gaps occur, from 1 to 6; the level share is 15.33 %; L = 1 holds 57 % of the level class.",
  100000, 1, gen_omegasum },
{ "bigomegasum", "A022559", "Sum of exponents in prime-power factorization of n!", "summatory",
  "16 different gaps occur, from 1 to 16; the level share is 16.65 %; L = 1 holds 52 % of the level class.",
  100000, 1, gen_bigomegasum },
{ "halfsum", "A005187", "a(n) = a(floor(n/2)) + n; also denominators in expansion of 1/sqrt(1-x) are 2^a(n); also 2n - number of 1's in binary expansion of 2n", "summatory",
  "17 different gaps occur, from 1 to 17; the level share is 12.50 %; L = 1 holds 70 % of the level class.",
  100000, 0, gen_halfsum },
{ "pascalodd", "A006046", "Total number of odd entries in first n rows of Pascal's triangle: a(0) = 0, a(1) = 1, a(2k) = 3*a(k), a(2k+1) = 2*a(k) + a(k+1).  a(n) = Sum_{i=0..n-1} 2^wt(i)", "summatory",
  "17 different gaps occur, from 1 to 65,536; the level share is 55.21 %; 2.2 % of terms are forced level (l <= d^2); L = 1 holds 33 % of the level class.",
  100000, 0, gen_pascalodd },
{ "unitsum", "A064608", "Partial sums of A034444: sum of number of unitary divisors from 1 to n", "summatory",
  "6 different gaps occur, from 2 to 64; the level share is 26.43 %; L = 1 holds 61 % of the level class.",
  100000, 1, gen_unitsum },
{ "r5_0234", "A047203", "Numbers that are congruent to {0, 2, 3, 4} mod 5", "residue class",
  "The gaps are 1 and 2; the level share is 13.32 %; L = 1 holds 88 % of the level class.",
  100000, 1, gen_r5_0234 },
{ "r5_14", "A047209", "Numbers that are congruent to {1, 4} mod 5", "residue class",
  "The gaps are 2 and 3; the level share is 18.90 %; L = 1 holds 58 % of the level class; there are no ties.",
  100000, 1, gen_r5_14 },
{ "r5_013", "A047220", "Numbers that are congruent to {0, 1, 3} mod 5", "residue class",
  "The gaps are 1 and 2; the level share is 15.49 %; L = 1 holds 74 % of the level class.",
  100000, 1, gen_r5_013 },
{ "r6_0234", "A047229", "Numbers that are congruent to {0, 2, 3, 4} mod 6", "residue class",
  "The gaps are 1 and 2; the level share is 14.30 %; L = 2 holds 52 % of the level class.",
  100000, 1, gen_r6_0234 },
{ "r6_02", "A047238", "Numbers that are congruent to {0, 2} mod 6", "residue class",
  "The gaps are 2 and 4; the level share is 17.56 %; L = 2 holds 79 % of the level class.",
  100000, 1, gen_r6_02 },
{ "r6_0123", "A047246", "Numbers that are congruent to {0, 1, 2, 3} mod 6", "residue class",
  "The gaps are 1 and 3; the level share is 13.85 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_r6_0123 },
{ "r6_1235", "A047255", "Numbers that are congruent to {1, 2, 3, 5} mod 6", "residue class",
  "The gaps are 1 and 2; the level share is 13.81 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_r6_1235 },
{ "r6_245", "A047261", "Numbers that are congruent to {2, 4, 5} mod 6", "residue class",
  "The gaps are 1, 2 and 3; the level share is 4.79 %; L = 2 holds 100 % of the level class.",
  100000, 1, gen_r6_245 },
{ "r6_015", "A047266", "Numbers that are congruent to {0, 1, 5} mod 6", "residue class",
  "The gaps are 1 and 4; the level share is 15.64 %; L = 1 holds 58 % of the level class.",
  100000, 1, gen_r6_015 },
{ "r6_0135", "A047273", "Numbers that are congruent to {0, 1, 3, 5} mod 6", "residue class",
  "The gaps are 1 and 2; the level share is 20.79 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_r6_0135 },
{ "r10_1379", "A045572", "Numbers that are odd but not divisible by 5", "residue class",
  "The gaps are 2 and 4; the level share is 18.57 %; L = 1 holds 89 % of the level class.",
  100000, 1, gen_r10_1379 },
{ "cop21", "A160545", "Numbers coprime to 21", "residue class",
  "The gaps are 1, 2 and 3; the level share is 4.05 %; L = 1 holds 65 % of the level class.",
  100000, 1, gen_cop21 },
{ "r4_03", "A014601", "Numbers congruent to 0 or 3 mod 4", "residue class",
  "The gaps are 1 and 3; the level share is 12.31 %; L = 1 holds 73 % of the level class.",
  100000, 0, gen_r4_03 },
{ "r4_12", "A042963", "Numbers congruent to 1 or 2 mod 4", "residue class",
  "The gaps are 1 and 3; the level share is 12.32 %; L = 1 holds 73 % of the level class.",
  100000, 1, gen_r4_12 },
{ "r5_24", "A047211", "Numbers that are congruent to {2, 4} mod 5", "residue class",
  "The gaps are 2 and 3; the level share is 10.49 %; L = 1 holds 52 % of the level class.",
  100000, 1, gen_r5_24 },
{ "r3_02", "A007494", "Numbers that are congruent to 0 or 2 mod 3", "residue class",
  "The gaps are 1 and 2; the level share is 17.52 %; L = 1 holds 79 % of the level class.",
  100000, 0, gen_r3_02 },
{ "r3_01", "A032766", "Numbers that are congruent to 0 or 1 (mod 3)", "residue class",
  "The gaps are 1 and 2; the level share is 17.56 %; L = 1 holds 79 % of the level class; there are no ties.",
  100000, 0, gen_r3_01 },
{ "r5_024", "A047212", "Numbers that are congruent to {0, 2, 4} mod 5", "residue class",
  "The gaps are 1 and 2; the level share is 9.66 %; L = 1 holds 79 % of the level class.",
  100000, 1, gen_r5_024 },
{ "a045372", "A045372", "Primes congruent to {1, 2} mod 5", "primes",
  "66 different gaps occur, from 4 to 246; the level share is 30.27 %.",
  100000, 1, gen_a045372 },
{ "a045429", "A045429", "Primes congruent to {1, 3} mod 5", "primes",
  "63 different gaps occur, from 2 to 258; the level share is 27.08 %.",
  100000, 1, gen_a045429 },
{ "a045378", "A045378", "Primes congruent to {2, 4} mod 5", "primes",
  "69 different gaps occur, from 2 to 262; the level share is 29.00 %.",
  100000, 1, gen_a045378 },
{ "a045435", "A045435", "Primes congruent to {3, 4} mod 5", "primes",
  "67 different gaps occur, from 4 to 250; the level share is 26.03 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a045435 },
{ "a045321", "A045321", "Primes congruent to {1, 2, 3} (mod 5)", "primes",
  "77 different gaps occur, from 1 to 192; the level share is 26.15 %; L = 3 holds 31 % of the level class.",
  100000, 1, gen_a045321 },
{ "a045371", "A045371", "Primes congruent to {1, 2, 4} mod 5", "primes",
  "73 different gaps occur, from 2 to 156; the level share is 26.55 %.",
  100000, 1, gen_a045371 },
{ "a045428", "A045428", "Primes congruent to {1, 3, 4} mod 5", "primes",
  "74 different gaps occur, from 2 to 168; the level share is 24.36 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a045428 },
{ "a045327", "A045327", "Primes congruent to {2, 3, 4} mod 5", "primes",
  "74 different gaps occur, from 1 to 192; the level share is 24.98 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_a045327 },
{ "a045392", "A045392", "Primes congruent to 2 mod 7", "primes",
  "53 different gaps occur, from 14 to 854; the level share is 39.08 %.",
  100000, 1, gen_a045392 },
{ "a045437", "A045437", "Primes congruent to 3 mod 7", "primes",
  "53 different gaps occur, from 14 to 980; the level share is 38.67 %; there are no ties.",
  100000, 1, gen_a045437 },
{ "a045471", "A045471", "Primes congruent to 4 mod 7", "primes",
  "54 different gaps occur, from 14 to 826; the level share is 39.11 %.",
  100000, 1, gen_a045471 },
{ "a045458", "A045458", "Primes congruent to 5 mod 7", "primes",
  "50 different gaps occur, from 14 to 756; the level share is 38.96 %; there are no ties.",
  100000, 1, gen_a045458 },
{ "a045473", "A045473", "Primes congruent to 6 mod 7", "primes",
  "54 different gaps occur, from 14 to 952; the level share is 39.05 %; there are no ties.",
  100000, 1, gen_a045473 },
{ "a045465", "A045465", "Primes congruent to {0, 1} mod 7", "primes",
  "51 different gaps occur, from 14 to 770; the level share is 38.90 %.",
  100000, 1, gen_a045465 },
{ "a045391", "A045391", "Primes congruent to {1, 2} mod 7", "primes",
  "75 different gaps occur, from 6 to 456; the level share is 30.02 %.",
  100000, 1, gen_a045391 },
{ "a045436", "A045436", "Primes congruent to {1, 3} mod 7", "primes",
  "75 different gaps occur, from 2 to 378; the level share is 33.94 %.",
  100000, 1, gen_a045436 },
{ "a045343", "A045343", "Primes congruent to {2, 3} mod 7", "primes",
  "76 different gaps occur, from 1 to 364; the level share is 33.79 %.",
  100000, 1, gen_a045343 },
{ "a045469", "A045469", "Primes congruent to {1, 4} mod 7", "primes",
  "73 different gaps occur, from 4 to 382; the level share is 30.02 %.",
  100000, 1, gen_a045469 },
{ "a045387", "A045387", "Primes congruent to {2, 4} mod 7", "primes",
  "71 different gaps occur, from 2 to 420; the level share is 29.52 %.",
  100000, 1, gen_a045387 },
{ "a045432", "A045432", "Primes congruent to {3, 4} mod 7", "primes",
  "70 different gaps occur, from 6 to 372; the level share is 34.04 %.",
  100000, 1, gen_a045432 },
{ "a045456", "A045456", "Primes congruent to {1, 5} mod 7", "primes",
  "72 different gaps occur, from 4 to 382; the level share is 34.58 %.",
  100000, 1, gen_a045456 },
{ "a045368", "A045368", "Primes congruent to {2, 5} mod 7", "primes",
  "69 different gaps occur, from 3 to 378; the level share is 34.34 %.",
  100000, 1, gen_a045368 },
{ "a045416", "A045416", "Primes congruent to {3, 5} mod 7", "primes",
  "71 different gaps occur, from 2 to 390; the level share is 29.91 %.",
  100000, 1, gen_a045416 },
{ "a045452", "A045452", "Primes congruent to {4, 5} mod 7", "primes",
  "74 different gaps occur, from 6 to 400; the level share is 34.01 %.",
  100000, 1, gen_a045452 },
{ "a045472", "A045472", "Primes congruent to {1, 6} mod 7", "primes",
  "75 different gaps occur, from 2 to 436; the level share is 33.75 %.",
  100000, 1, gen_a045472 },
{ "a045389", "A045389", "Primes congruent to {2, 6} mod 7", "primes",
  "74 different gaps occur, from 4 to 364; the level share is 34.38 %.",
  100000, 1, gen_a045389 },
{ "a045434", "A045434", "Primes congruent to {3, 6} mod 7", "primes",
  "76 different gaps occur, from 4 to 388; the level share is 29.01 %.",
  100000, 1, gen_a045434 },
{ "a045467", "A045467", "Primes congruent to {4, 6} mod 7", "primes",
  "73 different gaps occur, from 2 to 378; the level share is 34.09 %.",
  100000, 1, gen_a045467 },
{ "a045455", "A045455", "Primes congruent to {5, 6} mod 7", "primes",
  "74 different gaps occur, from 6 to 448; the level share is 27.48 %.",
  100000, 1, gen_a045455 },
{ "a045342", "A045342", "Primes congruent to {1, 2, 3} mod 7", "primes",
  "83 different gaps occur, from 1 to 286; the level share is 29.64 %.",
  100000, 1, gen_a045342 },
{ "a045386", "A045386", "Primes congruent to {1, 2, 4} mod 7", "primes",
  "103 different gaps occur, from 2 to 258; the level share is 25.87 %.",
  100000, 1, gen_a045386 },
{ "a023203", "A023203", "Primes p such that p + 10 is also prime", "primes",
  "189 different gaps occur, from 4 to 1,680; the level share is 44.96 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_a023203 },
{ "a046133", "A046133", "Primes p such that p + 12 is also prime", "primes",
  "349 different gaps occur, from 2 to 1,038; the level share is 35.23 %.",
  100000, 1, gen_a046133 },
{ "a049488", "A049488", "Primes p such that p+16 is prime", "primes",
  "251 different gaps occur, from 4 to 2,184; the level share is 47.23 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a049488 },
{ "a049481", "A049481", "Primes p such that p + 30 is also prime", "primes",
  "264 different gaps occur, from 2 to 668; the level share is 33.90 %.",
  100000, 1, gen_a049481 },
{ "a049489", "A049489", "Primes p such that p + 32 is also prime", "primes",
  "247 different gaps occur, from 6 to 1,908; the level share is 47.23 %; L = 1 holds 36 % of the level class; there are no ties.",
  100000, 1, gen_a049489 },
{ "a062284", "A062284", "Primes p such that p + 50 is also prime", "primes",
  "185 different gaps occur, from 6 to 1,476; the level share is 44.98 %; L = 1 holds 37 % of the level class; there are no ties.",
  100000, 1, gen_a062284 },
{ "a049482", "A049482", "Primes p such that p + 210 is also prime", "primes",
  "211 different gaps occur, from 2 to 672; the level share is 32.00 %.",
  100000, 1, gen_a049482 },
{ "a063909", "A063909", "Primes p such that 2*p - 5 is also prime", "primes",
  "201 different gaps occur, from 6 to 1,530; the level share is 45.29 %; L = 1 holds 38 % of the level class; there are no ties.",
  100000, 1, gen_a063909 },
{ "a063910", "A063910", "Primes p such that 2*p - 7 is also prime", "primes",
  "231 different gaps occur, from 2 to 1,770; the level share is 45.83 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_a063910 },
{ "a063911", "A063911", "Primes p such that 2*p - 9 is also prime", "primes",
  "343 different gaps occur, from 2 to 938; the level share is 36.31 %.",
  100000, 1, gen_a063911 },
{ "a063912", "A063912", "Primes p such that 2*p - 11 is also prime", "primes",
  "241 different gaps occur, from 4 to 2,100; the level share is 46.94 %; L = 1 holds 36 % of the level class; there are no ties.",
  100000, 1, gen_a063912 },
{ "a063913", "A063913", "Primes p such that 2*p - 13 is also prime", "primes",
  "247 different gaps occur, from 6 to 2,034; the level share is 47.08 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a063913 },
{ "a023209", "A023209", "Primes p such that 3p + 4 is also prime", "primes",
  "348 different gaps occur, from 2 to 1,030; the level share is 36.37 %.",
  100000, 1, gen_a023209 },
{ "a023210", "A023210", "Primes p such that 3*p + 8 is also prime", "primes",
  "355 different gaps occur, from 2 to 1,114; the level share is 37.24 %.",
  100000, 1, gen_a023210 },
{ "a023211", "A023211", "Primes p such that 3*p + 10 is also prime", "primes",
  "263 different gaps occur, from 2 to 784; the level share is 34.78 %.",
  100000, 1, gen_a023211 },
{ "a062737", "A062737", "Primes p such that 4p-1 is also prime", "primes",
  "270 different gaps occur, from 1 to 2,376; the level share is 47.85 %; L = 1 holds 36 % of the level class; there are no ties; 6 terms do not decompose.",
  100000, 1, gen_a062737 },
{ "a023213", "A023213", "Primes p such that 4p + 3 is prime", "primes",
  "361 different gaps occur, from 2 to 1,178; the level share is 37.18 %.",
  100000, 1, gen_a023213 },
{ "a023214", "A023214", "Primes p such that 4*p + 5 is also prime", "primes",
  "207 different gaps occur, from 1 to 1,536; the level share is 45.38 %; L = 1 holds 37 % of the level class; there are no ties.",
  100000, 1, gen_a023214 },
{ "a023215", "A023215", "Primes p such that 4*p + 7 is also prime", "primes",
  "228 different gaps occur, from 6 to 1,980; the level share is 46.05 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_a023215 },
{ "a023216", "A023216", "Primes p such that 4*p + 9 is also prime", "primes",
  "351 different gaps occur, from 2 to 1,282; the level share is 37.52 %.",
  100000, 1, gen_a023216 },
{ "a023217", "A023217", "Primes p such that 5*p + 2 is also prime", "primes",
  "210 different gaps occur, from 4 to 1,722; the level share is 45.21 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_a023217 },
{ "a023218", "A023218", "Primes p such that 5*p + 4 is also prime", "primes",
  "211 different gaps occur, from 2 to 1,848; the level share is 45.50 %; L = 1 holds 37 % of the level class; there are no ties; 6 terms do not decompose.",
  100000, 1, gen_a023218 },
{ "a023220", "A023220", "Primes p such that 5*p + 8 is also prime", "primes",
  "216 different gaps occur, from 4 to 1,758; the level share is 45.36 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_a023220 },
{ "a007693", "A007693", "Primes p such that 6*p + 1 is also prime", "primes",
  "366 different gaps occur, from 1 to 1,072; the level share is 37.09 %.",
  100000, 1, gen_a007693 },
{ "a023221", "A023221", "Primes p such that 6*p + 5 is also prime", "primes",
  "276 different gaps occur, from 1 to 726; the level share is 34.93 %.",
  100000, 1, gen_a023221 },
{ "a023222", "A023222", "Primes p such that 6*p + 7 is also prime", "primes",
  "318 different gaps occur, from 2 to 830; the level share is 36.21 %.",
  100000, 1, gen_a023222 },
{ "a023223", "A023223", "Primes p such that 7*p + 2 is also prime", "primes",
  "236 different gaps occur, from 2 to 1,680; the level share is 46.73 %; L = 1 holds 36 % of the level class; there are no ties.",
  100000, 1, gen_a023223 },
{ "a023224", "A023224", "Primes p such that 7*p + 4 is also prime", "primes",
  "227 different gaps occur, from 6 to 1,848; the level share is 47.08 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a023224 },
{ "a023225", "A023225", "Primes p such that 7*p + 6 is also prime", "primes",
  "314 different gaps occur, from 2 to 972; the level share is 35.45 %.",
  100000, 1, gen_a023225 },
{ "a023226", "A023226", "Primes p such that 7*p + 8 is also prime", "primes",
  "240 different gaps occur, from 2 to 1,800; the level share is 46.48 %; L = 1 holds 37 % of the level class; there are no ties; 7 terms do not decompose.",
  100000, 1, gen_a023226 },
{ "a023227", "A023227", "Primes p such that 7*p + 10 is also prime", "primes",
  "170 different gaps occur, from 4 to 1,302; the level share is 44.22 %; L = 1 holds 38 % of the level class.",
  100000, 1, gen_a023227 },
{ "a023229", "A023229", "Primes p such that 8*p + 3 is also prime", "primes",
  "379 different gaps occur, from 2 to 1,094; the level share is 37.72 %.",
  100000, 1, gen_a023229 },
{ "a023231", "A023231", "Primes p such that 8*p + 7 is also prime", "primes",
  "239 different gaps occur, from 1 to 1,884; the level share is 46.64 %; L = 1 holds 38 % of the level class; there are no ties.",
  100000, 1, gen_a023231 },
{ "a023232", "A023232", "Primes p such that 8*p + 9 is also prime", "primes",
  "371 different gaps occur, from 2 to 1,216; the level share is 36.70 %.",
  100000, 1, gen_a023232 },
{ "a023233", "A023233", "Primes p such that 9*p + 2 is also prime", "primes",
  "374 different gaps occur, from 2 to 1,102; the level share is 37.37 %.",
  100000, 1, gen_a023233 },
{ "a023234", "A023234", "Primes p such that 9*p + 4 is also prime", "primes",
  "369 different gaps occur, from 2 to 1,092; the level share is 37.66 %.",
  100000, 1, gen_a023234 },
{ "a023235", "A023235", "Primes p such that 9*p + 8 is also prime", "primes",
  "382 different gaps occur, from 2 to 1,198; the level share is 37.70 %.",
  100000, 1, gen_a023235 },
{ "a023236", "A023236", "Primes p such that 9*p + 10 is also prime", "primes",
  "278 different gaps occur, from 2 to 760; the level share is 34.78 %.",
  100000, 1, gen_a023236 },
{ "a023237", "A023237", "Primes p such that 10*p + 1 is also prime", "primes",
  "215 different gaps occur, from 4 to 1,842; the level share is 45.66 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a023237 },
{ "a023238", "A023238", "Primes p such that 10*p + 3 is also prime", "primes",
  "292 different gaps occur, from 2 to 778; the level share is 34.93 %.",
  100000, 1, gen_a023238 },
{ "a023239", "A023239", "Primes p such that 10*p + 7 is also prime", "primes",
  "182 different gaps occur, from 6 to 1,428; the level share is 44.46 %; L = 1 holds 38 % of the level class.",
  100000, 1, gen_a023239 },
{ "a023240", "A023240", "Primes p such that 10*p + 9 is also prime", "primes",
  "292 different gaps occur, from 2 to 1,022; the level share is 35.08 %.",
  100000, 1, gen_a023240 },
{ "a089443", "A089443", "Primes p such that 12*p + 13 is prime", "primes",
  "347 different gaps occur, from 2 to 1,000; the level share is 36.85 %.",
  100000, 1, gen_a089443 },
{ "a113169", "A113169", "Primes p such that 13*p + 2 is also prime", "primes",
  "262 different gaps occur, from 2 to 2,076; the level share is 47.39 %; L = 1 holds 36 % of the level class; there are no ties.",
  100000, 1, gen_a113169 },
{ "a113115", "A113115", "Primes p such that 17*p + 2 is also prime", "primes",
  "278 different gaps occur, from 6 to 2,958; the level share is 47.89 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_a113115 },
{ "a067076", "A067076", "Numbers k such that 2*k + 3 is a prime", "prime values",
  "53 different gaps occur, from 1 to 57; the level share is 15.96 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a067076 },
{ "a098090", "A098090", "Numbers k such that 2k-3 is prime", "prime values",
  "53 different gaps occur, from 1 to 57; the level share is 15.12 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a098090 },
{ "a089253", "A089253", "Numbers n such that 2n - 5 is a prime", "prime values",
  "53 different gaps occur, from 1 to 57; the level share is 19.99 %; L = 1 holds 43 % of the level class.",
  100000, 1, gen_a089253 },
{ "a089192", "A089192", "Numbers n such that 2n - 7 is a prime", "prime values",
  "53 different gaps occur, from 1 to 57; the level share is 21.25 %; L = 1 holds 44 % of the level class.",
  100000, 1, gen_a089192 },
{ "a102733", "A102733", "Numbers n such that 2*n + 101 is prime", "prime values",
  "53 different gaps occur, from 1 to 57; the level share is 22.56 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a102733 },
{ "a024892", "A024892", "Numbers k such that 3*k+1 is prime", "prime values",
  "39 different gaps occur, from 2 to 90; the level share is 17.92 %; L = 2 holds 49 % of the level class; there are no ties.",
  100000, 1, gen_a024892 },
{ "a087370", "A087370", "Numbers n such that 3n - 1 is a prime", "prime values",
  "41 different gaps occur, from 1 to 86; the level share is 18.09 %; L = 2 holds 49 % of the level class.",
  100000, 1, gen_a087370 },
{ "a024893", "A024893", "Numbers k such that 3*k+2 is prime", "prime values",
  "41 different gaps occur, from 1 to 86; the level share is 25.91 %; L = 1 holds 61 % of the level class.",
  100000, 1, gen_a024893 },
{ "a034936", "A034936", "Numbers k such that 3*k + 4 is prime", "prime values",
  "39 different gaps occur, from 2 to 90; the level share is 26.62 %; L = 1 holds 61 % of the level class.",
  100000, 1, gen_a034936 },
{ "a089953", "A089953", "Numbers n such that 3*n+7 is prime", "prime values",
  "39 different gaps occur, from 2 to 90; the level share is 17.15 %; L = 2 holds 47 % of the level class.",
  100000, 1, gen_a089953 },
{ "a005098", "A005098", "Numbers k such that 4k + 1 is prime", "prime values",
  "57 different gaps occur, from 1 to 62; the level share is 22.77 %; L = 1 holds 44 % of the level class.",
  100000, 1, gen_a005098 },
{ "a095278", "A095278", "Numbers k such that 4k + 3 is prime", "prime values",
  "56 different gaps occur, from 1 to 64; the level share is 16.38 %; L = 1 holds 32 % of the level class.",
  100000, 1, gen_a095278 },
{ "a111215", "A111215", "Numbers k such that 4k + 5 is prime", "prime values",
  "57 different gaps occur, from 1 to 62; the level share is 20.68 %; L = 1 holds 42 % of the level class.",
  100000, 1, gen_a111215 },
{ "a111199", "A111199", "Numbers k such that 4k + 9 is prime", "prime values",
  "57 different gaps occur, from 1 to 62; the level share is 16.34 %; L = 1 holds 32 % of the level class.",
  100000, 1, gen_a111199 },
{ "a024894", "A024894", "Numbers k such that 5*k + 1 is prime", "prime values",
  "47 different gaps occur, from 2 to 100; the level share is 21.54 %; L = 2 holds 48 % of the level class.",
  100000, 1, gen_a024894 },
{ "a024896", "A024896", "Numbers k such that 5*k - 2 is prime", "prime values",
  "50 different gaps occur, from 2 to 108; the level share is 30.09 %; L = 1 holds 61 % of the level class.",
  100000, 1, gen_a024896 },
{ "a111223", "A111223", "Numbers n such that 5*n + 2 is prime", "prime values",
  "49 different gaps occur, from 1 to 100; the level share is 30.25 %; L = 1 holds 63 % of the level class.",
  100000, 1, gen_a111223 },
{ "a024895", "A024895", "Numbers k such that 5*k - 3 is prime", "prime values",
  "49 different gaps occur, from 1 to 100; the level share is 14.61 %; L = 2 holds 35 % of the level class; there are no ties.",
  100000, 1, gen_a024895 },
{ "a087505", "A087505", "Numbers k such that 5*k+3 is a prime", "prime values",
  "50 different gaps occur, from 2 to 108; the level share is 15.92 %; L = 2 holds 37 % of the level class.",
  100000, 1, gen_a087505 },
{ "a024897", "A024897", "Numbers k such that 5*k + 4 is prime", "prime values",
  "47 different gaps occur, from 2 to 102; the level share is 29.95 %; L = 1 holds 62 % of the level class.",
  100000, 1, gen_a024897 },
{ "a081759", "A081759", "Numbers k such that 5*k+6 is prime", "prime values",
  "47 different gaps occur, from 2 to 100; the level share is 24.50 %; L = 1 holds 42 % of the level class.",
  100000, 1, gen_a081759 },
{ "a107304", "A107304", "Numbers k such that 5k - 7 is prime", "prime values",
  "50 different gaps occur, from 2 to 108; the level share is 20.27 %; L = 2 holds 46 % of the level class; there are no ties.",
  100000, 1, gen_a107304 },
{ "a111224", "A111224", "Numbers n such that 5*n + 7 is prime", "prime values",
  "48 different gaps occur, from 2 to 100; the level share is 20.22 %; L = 2 holds 46 % of the level class.",
  100000, 1, gen_a111224 },
{ "a111225", "A111225", "Numbers n such that 5*n + 8 is prime", "prime values",
  "50 different gaps occur, from 2 to 108; the level share is 30.58 %; L = 1 holds 62 % of the level class.",
  100000, 1, gen_a111225 },
{ "a111226", "A111226", "Numbers k such that 5*k + 12 is prime", "prime values",
  "48 different gaps occur, from 2 to 100; the level share is 25.10 %; L = 1 holds 41 % of the level class.",
  100000, 1, gen_a111226 },
{ "a111230", "A111230", "Numbers k such that 5*k + 14 is prime", "prime values",
  "47 different gaps occur, from 2 to 102; the level share is 29.06 %; L = 1 holds 60 % of the level class.",
  100000, 1, gen_a111230 },
{ "a024899", "A024899", "Numbers k such that 6*k + 1 is prime", "prime values",
  "39 different gaps occur, from 1 to 45; the level share is 17.92 %; L = 1 holds 49 % of the level class.",
  100000, 1, gen_a024899 },
{ "a059325", "A059325", "Numbers n such that 6n + 5 is prime", "prime values",
  "40 different gaps occur, from 1 to 43; the level share is 16.27 %; L = 1 holds 47 % of the level class.",
  100000, 1, gen_a059325 },
{ "a024905", "A024905", "Numbers k such that 7*k + 1 is prime", "prime values",
  "50 different gaps occur, from 2 to 110; the level share is 22.73 %; L = 2 holds 45 % of the level class; there are no ties.",
  100000, 1, gen_a024905 },
{ "a024901", "A024901", "Numbers k such that 7*k - 2 is prime", "prime values",
  "50 different gaps occur, from 2 to 108; the level share is 30.06 %; L = 1 holds 61 % of the level class.",
  100000, 1, gen_a024901 },
{ "a105772", "A105772", "Numbers k such that 7*k + 2 is prime", "prime values",
  "53 different gaps occur, from 2 to 122; the level share is 32.04 %; L = 1 holds 59 % of the level class.",
  100000, 1, gen_a105772 },
{ "a089033", "A089033", "Numbers n such that 7*n+3 is prime", "prime values",
  "53 different gaps occur, from 2 to 140; the level share is 16.42 %; L = 2 holds 34 % of the level class.",
  100000, 1, gen_a089033 },
{ "a024902", "A024902", "Numbers k such that 7*k + 4 is prime", "prime values",
  "54 different gaps occur, from 2 to 118; the level share is 30.69 %; L = 1 holds 59 % of the level class.",
  100000, 1, gen_a024902 },
{ "a024903", "A024903", "Numbers k such that 7*k - 4 is prime", "prime values",
  "53 different gaps occur, from 2 to 140; the level share is 31.58 %; L = 1 holds 61 % of the level class.",
  100000, 1, gen_a024903 },
{ "a024904", "A024904", "Numbers k such that 7*k - 5 is prime", "prime values",
  "53 different gaps occur, from 2 to 122; the level share is 20.17 %; L = 2 holds 44 % of the level class; there are no ties.",
  100000, 1, gen_a024904 },
{ "a111367", "A111367", "Numbers k such that 7*k + 5 is prime", "prime values",
  "50 different gaps occur, from 2 to 108; the level share is 20.04 %; L = 2 holds 45 % of the level class; there are no ties.",
  100000, 1, gen_a111367 },
{ "a024900", "A024900", "Numbers k such that 7*k + 6 is prime", "prime values",
  "54 different gaps occur, from 2 to 136; the level share is 24.95 %; L = 1 holds 43 % of the level class.",
  100000, 1, gen_a024900 },
{ "a111249", "A111249", "Numbers k such that 7*k + 8 is prime", "prime values",
  "50 different gaps occur, from 2 to 110; the level share is 30.20 %; L = 1 holds 59 % of the level class.",
  100000, 1, gen_a111249 },
{ "a111250", "A111250", "Numbers n such that 7*n + 10 is prime", "prime values",
  "53 different gaps occur, from 2 to 140; the level share is 29.43 %; L = 1 holds 58 % of the level class.",
  100000, 1, gen_a111250 },
{ "a033868", "A033868", "Numbers n such that 7*n-11 is prime", "prime values",
  "53 different gaps occur, from 2 to 140; the level share is 21.25 %; L = 2 holds 44 % of the level class.",
  100000, 1, gen_a033868 },
{ "a089079", "A089079", "Numbers n such that 7*n - 23 is prime", "prime values",
  "50 different gaps occur, from 2 to 108; the level share is 22.20 %; L = 2 holds 45 % of the level class.",
  100000, 1, gen_a089079 },
{ "a108601", "A108601", "Numbers n such that 7*n - 911 is prime", "prime values",
  "54 different gaps occur, from 2 to 136; the level share is 21.89 %; L = 2 holds 45 % of the level class; there are no ties.",
  100000, 1, gen_a108601 },
{ "a108935", "A108935", "Numbers k such that 7*k + 911 is prime", "prime values",
  "50 different gaps occur, from 2 to 110; the level share is 22.17 %; L = 2 holds 45 % of the level class; there are no ties.",
  100000, 1, gen_a108935 },
{ "a005122", "A005122", "Numbers k such that 8k - 1 is prime", "prime values",
  "59 different gaps occur, from 1 to 63; the level share is 22.81 %; L = 1 holds 44 % of the level class.",
  100000, 1, gen_a005122 },
{ "a005123", "A005123", "Numbers k such that 8k + 1 is prime", "prime values",
  "59 different gaps occur, from 1 to 66; the level share is 22.92 %; L = 1 holds 42 % of the level class.",
  100000, 1, gen_a005123 },
{ "a005124", "A005124", "Numbers k such that 8k + 3 is prime", "prime values",
  "60 different gaps occur, from 1 to 75; the level share is 15.83 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a005124 },
{ "a005125", "A005125", "Numbers k such that 8k - 3 is prime", "prime values",
  "59 different gaps occur, from 1 to 77; the level share is 17.44 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a005125 },
{ "a105133", "A105133", "Numbers n such that 8n + 5 is prime", "prime values",
  "59 different gaps occur, from 1 to 77; the level share is 20.81 %; L = 1 holds 43 % of the level class.",
  100000, 1, gen_a105133 },
{ "a024906", "A024906", "Numbers k such that 9*k + 1 is prime", "prime values",
  "43 different gaps occur, from 2 to 100; the level share is 18.15 %; L = 2 holds 44 % of the level class.",
  100000, 1, gen_a024906 },
{ "a024910", "A024910", "Numbers k such that 9*k - 2 is prime", "prime values",
  "44 different gaps occur, from 2 to 88; the level share is 28.03 %; L = 1 holds 58 % of the level class.",
  100000, 1, gen_a024910 },
{ "a024909", "A024909", "Numbers k such that 9*k - 4 is prime", "prime values",
  "43 different gaps occur, from 2 to 102; the level share is 26.11 %; L = 1 holds 59 % of the level class.",
  100000, 1, gen_a024909 },
{ "a024908", "A024908", "Numbers k such that 9*k - 5 is prime", "prime values",
  "43 different gaps occur, from 2 to 90; the level share is 16.84 %; L = 2 holds 45 % of the level class.",
  100000, 1, gen_a024908 },
{ "a024907", "A024907", "Numbers k such that 9*k - 7 is prime", "prime values",
  "45 different gaps occur, from 1 to 92; the level share is 17.65 %; L = 2 holds 47 % of the level class.",
  100000, 1, gen_a024907 },
{ "a024912", "A024912", "Numbers k such that 10*k + 1 is prime", "prime values",
  "47 different gaps occur, from 1 to 50; the level share is 21.54 %; L = 1 holds 48 % of the level class.",
  100000, 1, gen_a024912 },
{ "a105042", "A105042", "Numbers n such that 10n - 1 is prime", "prime values",
  "47 different gaps occur, from 1 to 51; the level share is 21.35 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a105042 },
{ "a024914", "A024914", "Numbers k such that 10*k - 3 is prime", "prime values",
  "48 different gaps occur, from 1 to 50; the level share is 14.60 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_a024914 },
{ "a005574", "A005574", "Numbers k such that k^2 + 1 is prime", "prime values",
  "88 different gaps occur, from 1 to 212; the level share is 23.93 %; L = 2 holds 37 % of the level class; there are no ties.",
  100000, 1, gen_a005574 },
{ "a028870", "A028870", "Numbers k such that k^2 - 2 is prime", "prime values",
  "65 different gaps occur, from 1 to 138; the level share is 30.08 %; L = 1 holds 53 % of the level class.",
  100000, 1, gen_a028870 },
{ "a067201", "A067201", "Numbers k such that k^2 + 2 is prime", "prime values",
  "60 different gaps occur, from 1 to 552; the level share is 28.87 %; L = 3 holds 55 % of the level class.",
  100000, 1, gen_a067201 },
{ "a028873", "A028873", "Numbers k such that k^2 - 3 is prime", "prime values",
  "85 different gaps occur, from 2 to 204; the level share is 18.88 %.",
  100000, 1, gen_a028873 },
{ "a049422", "A049422", "Numbers k such that k^2 + 3 is prime", "prime values",
  "107 different gaps occur, from 2 to 252; the level share is 20.69 %.",
  100000, 1, gen_a049422 },
{ "a007591", "A007591", "Numbers k such that k^2 + 4 is prime", "prime values",
  "82 different gaps occur, from 2 to 178; the level share is 33.26 %; L = 1 holds 51 % of the level class.",
  100000, 1, gen_a007591 },
{ "a028876", "A028876", "Numbers k such that k^2 - 5 is prime", "prime values",
  "66 different gaps occur, from 2 to 150; the level share is 19.27 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a028876 },
{ "a078402", "A078402", "Numbers k such that k^2 + 5 is prime", "prime values",
  "83 different gaps occur, from 6 to 624; the level share is 21.10 %; L = 6 holds 36 % of the level class; there are no ties.",
  100000, 1, gen_a078402 },
{ "a028879", "A028879", "Numbers k such that k^2 - 6 is prime", "prime values",
  "112 different gaps occur, from 2 to 260; the level share is 31.56 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a028879 },
{ "a114269", "A114269", "Numbers k such that k^2 + 6 is prime", "prime values",
  "160 different gaps occur, from 2 to 390; the level share is 35.31 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_a114269 },
{ "a028882", "A028882", "Numbers k such that k^2 - 7 is prime", "prime values",
  "61 different gaps occur, from 3 to 462; the level share is 18.45 %; L = 6 holds 42 % of the level class; there are no ties.",
  100000, 1, gen_a028882 },
{ "a114270", "A114270", "Numbers k such that k^2 + 7 is prime", "prime values",
  "57 different gaps occur, from 2 to 122; the level share is 19.16 %; L = 2 holds 43 % of the level class; there are no ties.",
  100000, 1, gen_a114270 },
{ "a028885", "A028885", "Numbers k such that k^2 - 8 is prime", "prime values",
  "62 different gaps occur, from 2 to 142; the level share is 29.80 %; L = 1 holds 53 % of the level class.",
  100000, 1, gen_a028885 },
{ "a114271", "A114271", "Numbers k such that k^2 + 8 is prime", "prime values",
  "58 different gaps occur, from 6 to 474; the level share is 28.95 %; L = 3 holds 54 % of the level class; there are no ties.",
  100000, 1, gen_a114271 },
{ "a114272", "A114272", "Numbers k such that k^2 + 9 is prime", "prime values",
  "130 different gaps occur, from 2 to 342; the level share is 23.34 %.",
  100000, 1, gen_a114272 },
{ "a114273", "A114273", "Numbers k such that k^2 + 10 is prime", "prime values",
  "107 different gaps occur, from 2 to 350; the level share is 32.69 %; L = 1 holds 44 % of the level class.",
  100000, 1, gen_a114273 },
{ "a114274", "A114274", "Numbers k such that k^2 + 11 is prime", "prime values",
  "86 different gaps occur, from 6 to 630; the level share is 23.33 %; L = 6 holds 38 % of the level class.",
  100000, 1, gen_a114274 },
{ "a114275", "A114275", "Numbers k such that k^2 + 12 is prime", "prime values",
  "105 different gaps occur, from 2 to 294; the level share is 30.11 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_a114275 },
{ "a113536", "A113536", "Numbers k such that k^2 + 13 is prime", "prime values",
  "88 different gaps occur, from 2 to 248; the level share is 22.91 %; L = 2 holds 35 % of the level class.",
  100000, 1, gen_a113536 },
{ "a121250", "A121250", "Numbers n such that n^2 + 14 is prime", "prime values",
  "105 different gaps occur, from 6 to 798; the level share is 34.20 %; L = 3 holds 48 % of the level class.",
  100000, 1, gen_a121250 },
{ "a121982", "A121982", "Numbers k such that k^2 + 15 is prime", "prime values",
  "87 different gaps occur, from 2 to 200; the level share is 18.27 %; there are no ties.",
  100000, 1, gen_a121982 },
{ "a122062", "A122062", "Numbers k such that k^2 + 16 is prime", "prime values",
  "83 different gaps occur, from 2 to 190; the level share is 33.31 %; L = 1 holds 51 % of the level class.",
  100000, 1, gen_a122062 },
{ "a020668", "A020668", "Numbers of the form x^2 + 4*y^2", "quadratic form",
  "32 different gaps occur, from 1 to 48; the level share is 20.23 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a020668 },
{ "a020674", "A020674", "Numbers of the form 2*x^2 + 5*y^2", "quadratic form",
  "51 different gaps occur, from 1 to 60; the level share is 20.66 %; L = 1 holds 38 % of the level class.",
  100000, 1, gen_a020674 },
{ "a020677", "A020677", "Numbers of form 3*x^2 + 4*y^2", "quadratic form",
  "42 different gaps occur, from 1 to 60; the level share is 25.26 %; L = 1 holds 44 % of the level class.",
  100000, 1, gen_a020677 },
{ "a020670", "A020670", "Numbers of form x^2 + 7y^2", "quadratic form",
  "42 different gaps occur, from 1 to 42; the level share is 17.26 %; L = 1 holds 38 % of the level class.",
  100000, 1, gen_a020670 },
{ "a020678", "A020678", "Numbers of form 3 x^2 + 5 y^2", "quadratic form",
  "62 different gaps occur, from 1 to 72; the level share is 30.54 %; L = 1 holds 53 % of the level class.",
  100000, 1, gen_a020678 },
{ "a020671", "A020671", "Numbers of form x^2 + 8 y^2", "quadratic form",
  "42 different gaps occur, from 1 to 72; the level share is 20.38 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_a020671 },
{ "a020675", "A020675", "Numbers of form 2 x^2 + 7 y^2", "quadratic form",
  "54 different gaps occur, from 1 to 57; the level share is 17.64 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_a020675 },
{ "a020682", "A020682", "Numbers of form 4 x^2 + 5 y^2", "quadratic form",
  "55 different gaps occur, from 1 to 87; the level share is 27.74 %; L = 1 holds 42 % of the level class.",
  100000, 1, gen_a020682 },
{ "a020672", "A020672", "Numbers of form x^2 + 9 y^2", "quadratic form",
  "53 different gaps occur, from 1 to 63; the level share is 24.48 %; L = 1 holds 44 % of the level class.",
  100000, 1, gen_a020672 },
{ "a020679", "A020679", "Numbers of form 3*x^2 + 7*y^2", "quadratic form",
  "78 different gaps occur, from 1 to 89; the level share is 24.03 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_a020679 },
{ "a020673", "A020673", "Numbers of form x^2 + 10 y^2", "quadratic form",
  "52 different gaps occur, from 1 to 55; the level share is 21.47 %; L = 1 holds 38 % of the level class.",
  100000, 1, gen_a020673 },
{ "a020676", "A020676", "Numbers of form 2 x^2 + 9 y^2", "quadratic form",
  "65 different gaps occur, from 1 to 72; the level share is 23.05 %; L = 1 holds 42 % of the level class.",
  100000, 1, gen_a020676 },
{ "a020680", "A020680", "Numbers of form 3 x^2 + 8 y^2", "quadratic form",
  "63 different gaps occur, from 1 to 87; the level share is 26.50 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_a020680 },
{ "a020683", "A020683", "Numbers of form 4 x^2 + 7 y^2", "quadratic form",
  "58 different gaps occur, from 1 to 87; the level share is 17.92 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a020683 },
{ "a020685", "A020685", "Numbers of form 5 x^2 + 6 y^2", "quadratic form",
  "85 different gaps occur, from 1 to 91; the level share is 28.77 %; L = 1 holds 41 % of the level class.",
  100000, 1, gen_a020685 },
{ "a020686", "A020686", "Numbers of form 5 x^2 + 7 y^2", "quadratic form",
  "87 different gaps occur, from 1 to 122; the level share is 25.79 %; L = 1 holds 40 % of the level class.",
  100000, 1, gen_a020686 },
{ "a020681", "A020681", "Numbers of form 3 x^2 + 10 y^2", "quadratic form",
  "82 different gaps occur, from 1 to 105; the level share is 28.69 %; L = 1 holds 42 % of the level class.",
  100000, 1, gen_a020681 },
{ "a020684", "A020684", "Numbers of form 4 x^2 + 9 y^2", "quadratic form",
  "67 different gaps occur, from 1 to 119; the level share is 29.90 %; L = 1 holds 42 % of the level class.",
  100000, 1, gen_a020684 },
{ "a020687", "A020687", "Numbers of form 5 x^2 + 8 y^2", "quadratic form",
  "73 different gaps occur, from 1 to 105; the level share is 26.86 %; there are no ties.",
  100000, 1, gen_a020687 },
{ "a020689", "A020689", "Numbers of form 6 x^2 + 7 y^2", "quadratic form",
  "96 different gaps occur, from 1 to 110; the level share is 23.20 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_a020689 },
{ "a020688", "A020688", "Numbers of form 5 x^2 + 9 y^2", "quadratic form",
  "99 different gaps occur, from 1 to 123; the level share is 31.68 %; L = 1 holds 39 % of the level class.",
  100000, 1, gen_a020688 },
{ "a020690", "A020690", "Numbers of form 7 x^2 + 8 y^2", "quadratic form",
  "84 different gaps occur, from 1 to 128; the level share is 22.74 %.",
  100000, 1, gen_a020690 },
{ "a020691", "A020691", "Numbers of form 7 x^2 + 9 y^2", "quadratic form",
  "107 different gaps occur, from 1 to 131; the level share is 27.29 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_a020691 },
{ "a020692", "A020692", "Numbers of form 7 x^2 + 10 y^2", "quadratic form",
  "114 different gaps occur, from 1 to 144; the level share is 22.57 %.",
  100000, 1, gen_a020692 },
{ "a020693", "A020693", "Numbers of the form 8*x^2 + 9*y^2", "quadratic form",
  "98 different gaps occur, from 1 to 160; the level share is 28.67 %; L = 1 holds 32 % of the level class.",
  100000, 1, gen_a020693 },
{ "a020694", "A020694", "Numbers of form 9 x^2 + 10 y^2", "quadratic form",
  "123 different gaps occur, from 1 to 177; the level share is 32.49 %; L = 1 holds 40 % of the level class.",
  100000, 1, gen_a020694 },
{ "a035121", "A035121", "Numbers of the form x^2+82*y^2", "quadratic form",
  "109 different gaps occur, from 1 to 128; the level share is 20.53 %.",
  100000, 1, gen_a035121 },
{ "a084865", "A084865", "Primes of the form 2x^2 + 3y^2", "quadratic form",
  "60 different gaps occur, from 1 to 510; the level share is 39.51 %; L = 1 holds 43 % of the level class; there are no ties.",
  100000, 1, gen_a084865 },
{ "a106857", "A106857", "Primes of the form x^2+xy+3y^2, with x and y nonnegative", "quadratic form",
  "135 different gaps occur, from 2 to 312; the level share is 31.99 %.",
  100000, 1, gen_a106857 },
{ "a106861", "A106861", "Primes of the form x^2+xy+4y^2, with x and y nonnegative", "quadratic form",
  "60 different gaps occur, from 12 to 660; the level share is 44.75 %; L = 1 holds 52 % of the level class.",
  100000, 1, gen_a106861 },
{ "a106866", "A106866", "Primes of the form 2x^2+xy+3y^2, with x and y nonnegative", "quadratic form",
  "405 different gaps occur, from 1 to 1,190; the level share is 36.12 %.",
  100000, 1, gen_a106866 },
{ "a033199", "A033199", "Primes of form x^2+6*y^2", "quadratic form",
  "60 different gaps occur, from 6 to 528; the level share is 39.51 %; L = 1 holds 43 % of the level class.",
  100000, 1, gen_a033199 },
{ "a106862", "A106862", "Primes of the form x^2+xy+5y^2, with x and y nonnegative", "quadratic form",
  "129 different gaps occur, from 2 to 336; the level share is 30.86 %.",
  100000, 1, gen_a106862 },
{ "a106871", "A106871", "Primes of the form 2x^2+xy+4y^2, with x and y nonnegative", "quadratic form",
  "399 different gaps occur, from 2 to 1,206; the level share is 36.56 %.",
  100000, 1, gen_a106871 },
{ "a106877", "A106877", "Primes of the form 3x^2+xy+3y^2, with x and y nonnegative", "quadratic form",
  "314 different gaps occur, from 4 to 1,464; the level share is 39.41 %.",
  100000, 1, gen_a106877 },
{ "a106889", "A106889", "Primes of the form 2x^2 + 5y^2", "quadratic form",
  "106 different gaps occur, from 2 to 584; the level share is 37.15 %.",
  100000, 1, gen_a106889 },
{ "a106869", "A106869", "Primes of the form x^2+xy+6y^2, with x and y nonnegative", "quadratic form",
  "403 different gaps occur, from 2 to 1,110; the level share is 36.23 %.",
  100000, 1, gen_a106869 },
{ "a106875", "A106875", "Primes of the form 3x^2+2xy+3y^2, with x and y nonnegative", "quadratic form",
  "161 different gaps occur, from 8 to 1,760; the level share is 40.24 %; there are no ties.",
  100000, 1, gen_a106875 },
{ "a106885", "A106885", "Primes of the form 2x^2+xy+5y^2, with x and y nonnegative", "quadratic form",
  "190 different gaps occur, from 3 to 1,434; the level share is 46.52 %; L = 1 holds 39 % of the level class; there are no ties.",
  100000, 1, gen_a106885 },
{ "a106894", "A106894", "Primes of the form 3x^2+xy+4y^2, with x and y nonnegative", "quadratic form",
  "639 different gaps occur, from 2 to 2,330; the level share is 39.86 %.",
  100000, 1, gen_a106894 },
{ "a106870", "A106870", "Primes of the form x^2+xy+7y^2, with x and y nonnegative", "quadratic form",
  "45 different gaps occur, from 6 to 360; the level share is 36.05 %; L = 1 holds 50 % of the level class.",
  100000, 1, gen_a106870 },
{ "a106882", "A106882", "Primes of the form 2x^2+2xy+5y^2, with x and y nonnegative", "quadratic form",
  "55 different gaps occur, from 3 to 900; the level share is 41.66 %; L = 1 holds 40 % of the level class; there are no ties.",
  100000, 1, gen_a106882 },
{ "a106892", "A106892", "Primes of the form 3x^2+2xy+4y^2, with x and y nonnegative", "quadratic form",
  "424 different gaps occur, from 2 to 1,172; the level share is 39.30 %.",
  100000, 1, gen_a106892 },
{ "a106897", "A106897", "Primes of the form 2x^2+xy+6y^2, with x and y nonnegative", "quadratic form",
  "633 different gaps occur, from 2 to 2,144; the level share is 40.01 %.",
  100000, 1, gen_a106897 },
{ "a106917", "A106917", "Primes of the form 2x^2 + 7y^2", "quadratic form",
  "343 different gaps occur, from 2 to 1,152; the level share is 35.77 %.",
  100000, 1, gen_a106917 },
{ "a106918", "A106918", "Primes of the form 3x^2+xy+5y^2, with x and y nonnegative", "quadratic form",
  "380 different gaps occur, from 2 to 1,056; the level share is 37.62 %.",
  100000, 1, gen_a106918 },
{ "a106923", "A106923", "Primes of the form 4x^2+xy+4y^2, with x and y nonnegative", "quadratic form",
  "376 different gaps occur, from 6 to 3,558; the level share is 47.14 %.",
  100000, 1, gen_a106923 },
{ "a106963", "A106963", "Primes of the form 4x^2 + 5y^2", "quadratic form",
  "152 different gaps occur, from 8 to 1,380; the level share is 41.59 %.",
  100000, 1, gen_a106963 },
{ "a102271", "A102271", "Primes of the form 3*x^2 + 7*y^2", "quadratic form",
  "88 different gaps occur, from 4 to 1,248; the level share is 41.13 %; there are no ties.",
  100000, 1, gen_a102271 },
{ "a106874", "A106874", "Primes of the form x^2+xy+8y^2, with x and y nonnegative", "quadratic form",
  "394 different gaps occur, from 2 to 1,192; the level share is 36.27 %.",
  100000, 1, gen_a106874 },
{ "a106883", "A106883", "Primes of the form 3x^2+3xy+4y^2, with x and y nonnegative", "quadratic form",
  "243 different gaps occur, from 6 to 1,746; the level share is 47.92 %; L = 1 holds 38 % of the level class.",
  100000, 1, gen_a106883 },
{ "a106910", "A106910", "Primes of the form 2x^2+xy+7y^2, with x and y nonnegative", "quadratic form",
  "313 different gaps occur, from 4 to 1,414; the level share is 44.30 %.",
  100000, 1, gen_a106910 },
{ "a106914", "A106914", "Primes of the form 3x^2+2xy+5y^2, with x and y nonnegative", "quadratic form",
  "416 different gaps occur, from 2 to 1,584; the level share is 36.79 %.",
  100000, 1, gen_a106914 },
{ "a106942", "A106942", "Primes of the form 3x^2+xy+6y^2, with x and y nonnegative", "quadratic form",
  "857 different gaps occur, from 2 to 2,812; the level share is 42.25 %.",
  100000, 1, gen_a106942 },
{ "a106956", "A106956", "Primes of the form 4x^2+xy+5y^2, with x and y nonnegative", "quadratic form",
  "629 different gaps occur, from 2 to 1,738; the level share is 40.00 %.",
  100000, 1, gen_a106956 },
{ "a033201", "A033201", "Primes of the form x^2 + 10*y^2", "quadratic form",
  "101 different gaps occur, from 2 to 552; the level share is 37.01 %.",
  100000, 1, gen_a033201 },
{ "a047215", "A047215", "Numbers that are congruent to {0, 2} mod 5", "residue class",
  "The gaps are 2 and 3; the level share is 18.93 %; L = 1 holds 58 % of the level class.",
  100000, 0, gen_a047215 },
{ "a047216", "A047216", "Numbers that are congruent to {1, 2} mod 5", "residue class",
  "The gaps are 1 and 4; the level share is 12.03 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a047216 },
{ "a047217", "A047217", "Numbers that are congruent to {0, 1, 2} mod 5", "residue class",
  "The gaps are 1 and 3; the level share is 11.05 %; L = 1 holds 69 % of the level class.",
  100000, 1, gen_a047217 },
{ "a047225", "A047225", "Numbers that are congruent to {0, 1} mod 6", "residue class",
  "The gaps are 1 and 5; the level share is 23.64 %; L = 1 holds 55 % of the level class; there are no ties.",
  100000, 1, gen_a047225 },
{ "a047240", "A047240", "Numbers that are congruent to {0, 1, 2} mod 6", "residue class",
  "The gaps are 1 and 4; the level share is 16.36 %; L = 1 holds 55 % of the level class.",
  100000, 1, gen_a047240 },
{ "a047241", "A047241", "Numbers that are congruent to {1, 3} mod 6", "residue class",
  "The gaps are 2 and 4; the level share is 26.05 %; L = 1 holds 100 % of the level class; there are no ties.",
  100000, 1, gen_a047241 },
{ "a047274", "A047274", "Numbers that are congruent to {0, 1} mod 7", "residue class",
  "The gaps are 1 and 6; the level share is 18.04 %; L = 1 holds 55 % of the level class.",
  100000, 1, gen_a047274 },
{ "a047352", "A047352", "Numbers that are congruent to {0, 2} mod 7", "residue class",
  "The gaps are 2 and 5; the level share is 19.71 %; L = 1 holds 51 % of the level class.",
  100000, 1, gen_a047352 },
{ "a047353", "A047353", "Numbers that are congruent to {1, 2} mod 7", "residue class",
  "The gaps are 1 and 6; the level share is 13.06 %; L = 1 holds 38 % of the level class.",
  100000, 1, gen_a047353 },
{ "a008590", "A008590", "Multiples of 8", "residue class",
  "The gap is always 8; the level share is 9.60 %; L = 8 holds 100 % of the level class.",
  100000, 0, gen_a008590 },
{ "a047393", "A047393", "Numbers that are congruent to {0, 1} mod 8", "residue class",
  "The gaps are 1 and 7; the level share is 20.80 %; L = 2 holds 43 % of the level class; there are no ties.",
  100000, 1, gen_a047393 },
{ "a047467", "A047467", "Numbers that are congruent to {0, 2} mod 8", "residue class",
  "The gaps are 2 and 6; the level share is 18.60 %; L = 4 holds 52 % of the level class.",
  100000, 1, gen_a047467 },
{ "a090570", "A090570", "Numbers that are congruent to {0, 1} mod 9", "residue class",
  "The gaps are 1 and 8; the level share is 21.19 %; L = 1 holds 59 % of the level class; there are no ties.",
  100000, 1, gen_a090570 },
{ "a087444", "A087444", "Numbers that are congruent to {1, 4} mod 9", "residue class",
  "The gaps are 3 and 6; the level share is 22.46 %; L = 1 holds 56 % of the level class.",
  100000, 1, gen_a087444 },
{ "a054966", "A054966", "Numbers that are congruent to {0, 1, 8} mod 9", "residue class",
  "The gaps are 1 and 7; the level share is 16.03 %; L = 1 holds 54 % of the level class.",
  100000, 1, gen_a054966 },
{ "a090773", "A090773", "Numbers that are congruent to {4, 6} mod 10", "residue class",
  "The gaps are 2 and 8; the level share is 17.55 %; L = 2 holds 63 % of the level class; there are no ties.",
  100000, 1, gen_a090773 },
{ "a078309", "A078309", "Numbers that are congruent to {1, 4, 7} mod 10", "residue class",
  "The gaps are 3 and 4; the level share is 23.39 %; L = 1 holds 61 % of the level class.",
  100000, 1, gen_a078309 },
{ "a090772", "A090772", "Numbers that are congruent to {2, 8} mod 10", "residue class",
  "The gaps are 4 and 6; the level share is 18.90 %; L = 2 holds 58 % of the level class.",
  100000, 1, gen_a090772 },
{ "a008593", "A008593", "Multiples of 11", "residue class",
  "The gap is always 11; the level share is 9.61 %; L = 11 holds 100 % of the level class.",
  100000, 0, gen_a008593 },
{ "a008594", "A008594", "Multiples of 12", "residue class",
  "The gap is always 12; the level share is 9.60 %; L = 12 holds 100 % of the level class.",
  100000, 0, gen_a008594 },
{ "a083031", "A083031", "Numbers that are congruent to {0, 3, 7} mod 12", "residue class",
  "The gaps are 3, 4 and 5; the level share is 23.72 %; L = 2 holds 38 % of the level class.",
  100000, 1, gen_a083031 },
{ "a083030", "A083030", "Numbers that are congruent to {0, 4, 7} mod 12", "residue class",
  "The gaps are 3, 4 and 5; the level share is 22.22 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a083030 },
{ "a008595", "A008595", "Multiples of 13", "residue class",
  "The gap is always 13; the level share is 9.62 %; L = 13 holds 100 % of the level class.",
  100000, 0, gen_a008595 },
{ "a092476", "A092476", "Numbers that are congruent to {1, 3, 9} mod 13", "residue class",
  "The gaps are 2, 5 and 6; the level share is 19.86 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a092476 },
{ "a008596", "A008596", "Multiples of 14", "residue class",
  "The gap is always 14; the level share is 9.61 %; L = 14 holds 100 % of the level class; there are no ties.",
  100000, 0, gen_a008596 },
{ "a113806", "A113806", "Numbers that are congruent to {6, 8} mod 14", "residue class",
  "The gaps are 2 and 12; the level share is 18.04 %; L = 2 holds 55 % of the level class; there are no ties.",
  100000, 1, gen_a113806 },
{ "a113805", "A113805", "Numbers that are congruent to {5, 9} mod 14", "residue class",
  "The gaps are 4 and 10; the level share is 29.19 %; L = 1 holds 65 % of the level class.",
  100000, 1, gen_a113805 },
{ "a008597", "A008597", "Multiples of 15", "residue class",
  "The gap is always 15; the level share is 9.61 %; L = 15 holds 100 % of the level class; there are no ties.",
  100000, 0, gen_a008597 },
{ "a087446", "A087446", "Numbers that are congruent to {1, 6} mod 15", "residue class",
  "The gaps are 5 and 10; the level share is 29.71 %; L = 1 holds 51 % of the level class; there are no ties.",
  100000, 1, gen_a087446 },
{ "a008598", "A008598", "Multiples of 16", "residue class",
  "The gap is always 16; the level share is 9.61 %; L = 16 holds 100 % of the level class.",
  100000, 0, gen_a008598 },
{ "a106839", "A106839", "Numbers congruent to 11 mod 16", "residue class",
  "The gap is always 16; the level share is 33.07 %; L = 1 holds 46 % of the level class; there are no ties.",
  100000, 0, gen_a106839 },
{ "a008599", "A008599", "Multiples of 17", "residue class",
  "The gap is always 17; the level share is 9.63 %; L = 17 holds 100 % of the level class.",
  100000, 0, gen_a008599 },
{ "a008600", "A008600", "Multiples of 18", "residue class",
  "The gap is always 18; the level share is 9.61 %; L = 18 holds 100 % of the level class.",
  100000, 0, gen_a008600 },
{ "a008601", "A008601", "Multiples of 19", "residue class",
  "The gap is always 19; the level share is 9.64 %; L = 19 holds 100 % of the level class.",
  100000, 0, gen_a008601 },
{ "a008602", "A008602", "Multiples of 20", "residue class",
  "The gap is always 20; the level share is 9.61 %; L = 20 holds 100 % of the level class; there are no ties.",
  100000, 0, gen_a008602 },
{ "a008603", "A008603", "Multiples of 21", "residue class",
  "The gap is always 21; the level share is 9.62 %; L = 21 holds 100 % of the level class; there are no ties.",
  100000, 0, gen_a008603 },
{ "a008604", "A008604", "Multiples of 22", "residue class",
  "The gap is always 22; the level share is 9.62 %; L = 22 holds 100 % of the level class.",
  100000, 0, gen_a008604 },
{ "a008605", "A008605", "Multiples of 23", "residue class",
  "The gap is always 23; the level share is 9.65 %; L = 23 holds 100 % of the level class.",
  100000, 0, gen_a008605 },
{ "a008606", "A008606", "Multiples of 24", "residue class",
  "The gap is always 24; the level share is 9.62 %; L = 24 holds 100 % of the level class; there are no ties.",
  100000, 0, gen_a008606 },
{ "a033430", "A033430", "a(n) = 4*n^3", "polynomial",
  "Every gap is different, from 4 to 119,998,800,004; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a033430 },
{ "a033431", "A033431", "a(n) = 2*n^3", "polynomial",
  "Every gap is different, from 2 to 59,999,400,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a033431 },
{ "a084377", "A084377", "a(n) = n^3 + 7", "polynomial",
  "Every gap is different, from 1 to 29,999,700,001; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a084377 },
{ "a084378", "A084378", "a(n) = n^3 + 3", "polynomial",
  "Every gap is different, from 1 to 29,999,700,001; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a084378 },
{ "a084380", "A084380", "a(n) = n^3 + 2", "polynomial",
  "Every gap is different, from 1 to 29,999,700,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a084380 },
{ "a084381", "A084381", "a(n) = n^3 + 5", "polynomial",
  "Every gap is different, from 1 to 29,999,700,001; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a084381 },
{ "a084382", "A084382", "a(n) = n^3 + 6", "polynomial",
  "Every gap is different, from 1 to 29,999,700,001; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a084382 },
{ "a117642", "A117642", "a(n) = 3*n^3", "polynomial",
  "Every gap is different, from 3 to 89,999,100,003; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a117642 },
{ "a084379", "A084379", "a(n) = n^3 + 17", "polynomial",
  "Every gap is different, from 1 to 29,999,700,001; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a084379 },
{ "a033562", "A033562", "a(n) = 2*n^3 + 1", "polynomial",
  "Every gap is different, from 2 to 59,999,400,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a033562 },
{ "a100214", "A100214", "a(n) = 4*n^3 + 4", "polynomial",
  "Every gap is different, from 4 to 119,998,800,004; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a100214 },
{ "a118465", "A118465", "a(n) = 8*n^3 + n", "polynomial",
  "Every gap is different, from 9 to 239,997,600,009; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a118465 },
{ "a003777", "A003777", "a(n) = n^3 + n^2 - 1", "polynomial",
  "Every gap is different, from 10 to 30,000,500,002; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a003777 },
{ "a005491", "A005491", "a(n) = n^3 + 3*n + 1", "polynomial",
  "Every gap is different, from 4 to 29,999,700,004; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a005491 },
{ "a011379", "A011379", "a(n) = n^2*(n+1)", "polynomial",
  "Every gap is different, from 2 to 29,999,900,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a011379 },
{ "a027444", "A027444", "a(n) = n^3 + n^2 + n", "polynomial",
  "Every gap is different, from 3 to 29,999,900,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a027444 },
{ "a098547", "A098547", "a(n) = n^3 + n^2 + 1", "polynomial",
  "Every gap is different, from 2 to 29,999,900,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a098547 },
{ "a105374", "A105374", "a(n) = 4*n^3 + 4*n", "polynomial",
  "Every gap is different, from 8 to 119,998,800,008; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a105374 },
{ "a114364", "A114364", "a(n) = n*(n+1)^2", "polynomial",
  "Every gap is different, from 14 to 30,000,700,004; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a114364 },
{ "a119536", "A119536", "a(n) = 3*n^3 + 3*n", "polynomial",
  "Every gap is different, from 6 to 89,999,100,006; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a119536 },
{ "a122562", "A122562", "a(n) = n^3 + 114 * n", "polynomial",
  "Every gap is different, from 121 to 30,000,300,115; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a122562 },
{ "a006002", "A006002", "a(n) = n*(n+1)^2/2", "polynomial",
  "Every gap is different, from 2 to 15,000,050,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a006002 },
{ "a006527", "A006527", "a(n) = (n^3 + 2*n)/3", "polynomial",
  "Every gap is different, from 1 to 9,999,900,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a006527 },
{ "a015237", "A015237", "a(n) = (2*n - 1)*n^2", "polynomial",
  "Every gap is different, from 1 to 59,999,200,003; every decomposable term is forced level (l <= d^2); 8 terms do not decompose.",
  100000, 0, gen_a015237 },
{ "a053698", "A053698", "a(n) = n^3 + n^2 + n + 1", "polynomial",
  "Every gap is different, from 3 to 29,999,900,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a053698 },
{ "a084367", "A084367", "a(n) = n*(2*n+1)^2", "polynomial",
  "Every gap is different, from 9 to 119,999,600,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a084367 },
{ "a089207", "A089207", "a(n) = 4*n^3 + 2*n^2", "polynomial",
  "Every gap is different, from 34 to 120,001,600,006; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a089207 },
{ "a099721", "A099721", "a(n) = n^2*(2*n+1)", "polynomial",
  "Every gap is different, from 3 to 59,999,600,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a099721 },
{ "a100109", "A100109", "a(n) = n^3 - 2*n^2 + 2", "polynomial",
  "Every gap is different, from 1 to 29,999,899,999; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a100109 },
{ "a100705", "A100705", "a(n) = n^3 + (n+1)^2", "polynomial",
  "Every gap is different, from 4 to 29,999,900,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a100705 },
{ "a028347", "A028347", "a(n) = n^2 - 4", "polynomial",
  "Every gap is different, from 5 to 200,003; every decomposable term is forced level (l <= d^2).",
  100000, 2, gen_a028347 },
{ "a028872", "A028872", "a(n) = n^2 - 3", "polynomial",
  "Every gap is different, from 5 to 200,003; every decomposable term is forced level (l <= d^2).",
  100000, 2, gen_a028872 },
{ "a028881", "A028881", "a(n) = n^2 - 7", "polynomial",
  "Every gap is different, from 7 to 200,005; every decomposable term is forced level (l <= d^2).",
  100000, 3, gen_a028881 },
{ "a033428", "A033428", "a(n) = 3*n^2", "polynomial",
  "Every gap is different, from 3 to 599,997; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a033428 },
{ "a033429", "A033429", "a(n) = 5*n^2", "polynomial",
  "Every gap is different, from 5 to 999,995; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a033429 },
{ "a033581", "A033581", "a(n) = 6*n^2", "polynomial",
  "Every gap is different, from 6 to 1,199,994; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a033581 },
{ "a033582", "A033582", "a(n) = 7*n^2", "polynomial",
  "Every gap is different, from 7 to 1,399,993; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a033582 },
{ "a059100", "A059100", "a(n) = n^2 + 2", "polynomial",
  "Every gap is different, from 1 to 199,999; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a059100 },
{ "a087475", "A087475", "a(n) = n^2 + 4", "polynomial",
  "Every gap is different, from 1 to 199,999; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a087475 },
{ "a114949", "A114949", "a(n) = n^2 + 6", "polynomial",
  "Every gap is different, from 1 to 199,999; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a114949 },
{ "a117619", "A117619", "a(n) = n^2 + 7", "polynomial",
  "Every gap is different, from 1 to 199,999; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a117619 },
{ "a117950", "A117950", "a(n) = n^2 + 3", "polynomial",
  "Every gap is different, from 1 to 199,999; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a117950 },
{ "a117951", "A117951", "a(n) = n^2 + 5", "polynomial",
  "Every gap is different, from 1 to 199,999; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a117951 },
{ "a033583", "A033583", "a(n) = 10*n^2", "polynomial",
  "Every gap is different, from 10 to 1,999,990; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a033583 },
{ "a033584", "A033584", "a(n) = 11*n^2", "polynomial",
  "Every gap is different, from 11 to 2,199,989; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a033584 },
{ "a064761", "A064761", "a(n) = 15*n^2", "polynomial",
  "Every gap is different, from 15 to 2,999,985; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a064761 },
{ "a064762", "A064762", "a(n) = 21*n^2", "polynomial",
  "Every gap is different, from 21 to 4,199,979; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a064762 },
{ "a064763", "A064763", "a(n) = 28*n^2", "polynomial",
  "Every gap is different, from 28 to 5,599,972; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a064763 },
{ "a114948", "A114948", "a(n) = n^2 + 10", "polynomial",
  "Every gap is different, from 1 to 199,999; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a114948 },
{ "a114962", "A114962", "a(n) = n^2 + 14", "polynomial",
  "Every gap is different, from 1 to 199,999; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a114962 },
{ "a114963", "A114963", "a(n) = n^2 + 22", "polynomial",
  "Every gap is different, from 1 to 199,999; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a114963 },
{ "a114964", "A114964", "a(n) = n^2 + 30", "polynomial",
  "Every gap is different, from 1 to 199,999; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a114964 },
{ "a114965", "A114965", "a(n) = n^2 + 34", "polynomial",
  "Every gap is different, from 1 to 199,999; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a114965 },
{ "a016766", "A016766", "a(n) = (3*n)^2", "polynomial",
  "Every gap is different, from 9 to 1,799,991; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a016766 },
{ "a016802", "A016802", "a(n) = (4*n)^2", "polynomial",
  "Every gap is different, from 16 to 3,199,984; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a016802 },
{ "a016850", "A016850", "a(n) = (5*n)^2", "polynomial",
  "Every gap is different, from 25 to 4,999,975; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a016850 },
{ "a016910", "A016910", "a(n) = (6*n)^2", "polynomial",
  "Every gap is different, from 36 to 7,199,964; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a016910 },
{ "a016982", "A016982", "a(n) = (7*n)^2", "polynomial",
  "Every gap is different, from 49 to 9,799,951; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a016982 },
{ "a017066", "A017066", "a(n) = (8*n)^2", "polynomial",
  "Every gap is different, from 64 to 12,799,936; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a017066 },
{ "a017162", "A017162", "a(n) = (9*n)^2", "polynomial",
  "Every gap is different, from 81 to 16,199,919; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a017162 },
{ "a027688", "A027688", "a(n) = n^2 + n + 3", "polynomial",
  "Every gap is different, from 2 to 200,000; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a027688 },
{ "a027689", "A027689", "a(n) = n^2 + n + 4", "polynomial",
  "Every gap is different, from 2 to 200,000; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a027689 },
{ "a004919", "A004919", "a(n) = floor(n*phi^4), where phi is the golden ratio, A001622", "Beatty",
  "The gaps are 6 and 7; the level share is 22.26 %; L = 1 holds 36 % of the level class.",
  100000, 0, gen_a004919 },
{ "a004920", "A004920", "a(n) = floor(n*phi^5), where phi is the golden ratio, A001622", "Beatty",
  "The gaps are 11 and 12; the level share is 25.81 %; L = 1 holds 30 % of the level class.",
  100000, 0, gen_a004920 },
{ "a004921", "A004921", "a(n) = floor(n*phi^6), phi = golden ratio, A001622", "Beatty",
  "The gaps are 17 and 18; the level share is 29.17 %.",
  100000, 0, gen_a004921 },
{ "a004922", "A004922", "a(n) = floor(n*phi^7), where phi is the golden ratio, A001622", "Beatty",
  "The gaps are 29 and 30; the level share is 32.36 %.",
  100000, 0, gen_a004922 },
{ "a004976", "A004976", "a(n) = floor(n*phi^3), where phi=(1+sqrt(5))/2", "Beatty",
  "The gaps are 4 and 5; the level share is 18.86 %; L = 1 holds 44 % of the level class.",
  100000, 0, gen_a004976 },
{ "a037085", "A037085", "Beatty sequence for Pi^2", "Beatty",
  "The gaps are 9 and 10; the level share is 24.76 %; L = 1 holds 32 % of the level class.",
  100000, 0, gen_a037085 },
{ "a037086", "A037086", "Beatty sequence for sqrt(Pi)", "Beatty",
  "The gaps are 1 and 2; the level share is 12.83 %; L = 1 holds 70 % of the level class.",
  100000, 0, gen_a037086 },
{ "a037087", "A037087", "Beatty sequence for e^(1/e)", "Beatty",
  "The gaps are 1 and 2; the level share is 11.42 %; L = 1 holds 81 % of the level class.",
  100000, 0, gen_a037087 },
{ "a038130", "A038130", "Beatty sequence for 2*Pi", "Beatty",
  "The gaps are 6 and 7; the level share is 21.91 %; L = 1 holds 38 % of the level class.",
  100000, 0, gen_a038130 },
{ "a038152", "A038152", "Beatty sequence for e^Pi", "Beatty",
  "The gaps are 23 and 24; the level share is 30.71 %.",
  100000, 1, gen_a038152 },
{ "a038153", "A038153", "Beatty sequence for Pi^e", "Beatty",
  "The gaps are 22 and 23; the level share is 30.45 %.",
  100000, 1, gen_a038153 },
{ "a054386", "A054386", "Beatty sequence for Pi/(Pi-1); complement of A022844", "Beatty",
  "The gaps are 1 and 2; the level share is 11.64 %; L = 1 holds 80 % of the level class.",
  100000, 1, gen_a054386 },
{ "a054965", "A054965", "Beatty sequence for log_3(10), i.e., for 1/log_10(3); so largest exponent of 3 which produces an n-digit decimal number", "Beatty",
  "The gaps are 2 and 3; the level share is 14.03 %; L = 1 holds 64 % of the level class.",
  100000, 1, gen_a054965 },
{ "a059531", "A059531", "Beatty sequence for 1 + 1/Pi", "Beatty",
  "The gaps are 1 and 2; the level share is 11.03 %; L = 1 holds 86 % of the level class.",
  100000, 1, gen_a059531 },
{ "a059532", "A059532", "Beatty sequence for 1 + Pi", "Beatty",
  "The gaps are 4 and 5; the level share is 18.38 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a059532 },
{ "a059535", "A059535", "Beatty sequence for Pi^2/6, or zeta(2)", "Beatty",
  "The gaps are 1 and 2; the level share is 12.37 %; L = 1 holds 75 % of the level class.",
  100000, 1, gen_a059535 },
{ "a059536", "A059536", "Beatty sequence for zeta(2)/(zeta(2)-1)", "Beatty",
  "The gaps are 2 and 3; the level share is 15.22 %; L = 1 holds 58 % of the level class.",
  100000, 1, gen_a059536 },
{ "a059537", "A059537", "Beatty sequence for zeta(3)", "Beatty",
  "The gaps are 1 and 2; the level share is 10.45 %; L = 1 holds 90 % of the level class.",
  100000, 1, gen_a059537 },
{ "a059538", "A059538", "Beatty sequence for zeta(3)/(zeta(3)-1)", "Beatty",
  "The gaps are 5 and 6; the level share is 21.31 %; L = 1 holds 39 % of the level class.",
  100000, 1, gen_a059538 },
{ "a059539", "A059539", "Beatty sequence for 3^(1/3)", "Beatty",
  "The gaps are 1 and 2; the level share is 11.53 %; L = 1 holds 81 % of the level class.",
  100000, 1, gen_a059539 },
{ "a059540", "A059540", "Beatty sequence for 3^(1/3)/(3^(1/3)-1)", "Beatty",
  "The gaps are 3 and 4; the level share is 17.09 %; L = 1 holds 50 % of the level class.",
  100000, 1, gen_a059540 },
{ "a059541", "A059541", "Beatty sequence for 1 + log(2)", "Beatty",
  "The gaps are 1 and 2; the level share is 12.45 %; L = 1 holds 73 % of the level class.",
  100000, 1, gen_a059541 },
{ "a059542", "A059542", "Beatty sequence for 1 + 1/log(2)", "Beatty",
  "The gaps are 2 and 3; the level share is 15.01 %; L = 1 holds 59 % of the level class.",
  100000, 1, gen_a059542 },
{ "a059543", "A059543", "Beatty sequence for log(3)", "Beatty",
  "The gaps are 1 and 2; the level share is 9.96 %; L = 1 holds 95 % of the level class.",
  100000, 1, gen_a059543 },
{ "a059544", "A059544", "Beatty sequence for log(3)/(log(3)-1)", "Beatty",
  "The gaps are 11 and 12; the level share is 25.73 %; L = 1 holds 30 % of the level class.",
  100000, 1, gen_a059544 },
{ "a059545", "A059545", "Beatty sequence for log(10)", "Beatty",
  "The gaps are 2 and 3; the level share is 14.65 %; L = 1 holds 60 % of the level class.",
  100000, 1, gen_a059545 },
{ "a059546", "A059546", "Beatty sequence for log(10)/(log(10)-1)", "Beatty",
  "The gaps are 1 and 2; the level share is 12.88 %; L = 1 holds 71 % of the level class.",
  100000, 1, gen_a059546 },
{ "a059547", "A059547", "Beatty sequence for 1 + 1/log(3)", "Beatty",
  "The gaps are 1 and 2; the level share is 13.36 %; L = 1 holds 67 % of the level class.",
  100000, 1, gen_a059547 },
{ "a059548", "A059548", "Beatty sequence for 1 + log(3)", "Beatty",
  "The gaps are 2 and 3; the level share is 14.04 %; L = 1 holds 64 % of the level class.",
  100000, 1, gen_a059548 },
{ "a059549", "A059549", "Beatty sequence for 1 + 1/log(10)", "Beatty",
  "The gaps are 1 and 2; the level share is 11.56 %; L = 1 holds 81 % of the level class.",
  100000, 1, gen_a059549 },
{ "a027697", "A027697", "Odious primes: primes with odd number of 1's in binary expansion", "primes",
  "125 different gaps occur, from 2 to 304; the level share is 29.63 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_a027697 },
{ "a027699", "A027699", "Evil primes: primes with even number of 1's in their binary expansion", "primes",
  "130 different gaps occur, from 2 to 352; the level share is 30.29 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a027699 },
{ "a007675", "A007675", "Numbers m such that m, m+1 and m+2 are squarefree", "multiplicative",
  "12 different gaps occur, from 4 to 48; the level share is 28.74 %; L = 1 holds 53 % of the level class.",
  100000, 1, gen_a007675 },
{ "a039955", "A039955", "Squarefree numbers congruent to 1 (mod 4)", "multiplicative",
  "6 different gaps occur, from 4 to 24; the level share is 22.88 %; L = 1 holds 55 % of the level class.",
  100000, 1, gen_a039955 },
{ "a001633", "A001633", "Numbers with an odd number of digits", "digit rule",
  "The gaps are 1, 91, 9001 and 900001; the level share is 9.19 %; L = 1 holds 100 % of the level class; 6 terms do not decompose.",
  100000, 1, gen_a001633 },
{ "a001637", "A001637", "Numbers with an even number of digits", "digit rule",
  "The gaps are 1, 901 and 90001; the level share is 8.73 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a001637 },
{ "a034709", "A034709", "Numbers divisible by their last digit", "digit rule",
  "6 different gaps occur, from 1 to 6; the level share is 15.18 %; L = 1 holds 79 % of the level class.",
  100000, 1, gen_a034709 },
{ "a038770", "A038770", "Numbers divisible by at least one of their digits", "digit rule",
  "6 different gaps occur, from 1 to 6; the level share is 11.07 %; L = 1 holds 87 % of the level class.",
  100000, 1, gen_a038770 },
{ "a064150", "A064150", "Numbers divisible by the sum of their ternary digits", "digit rule",
  "57 different gaps occur, from 1 to 70; the level share is 19.66 %.",
  100000, 1, gen_a064150 },
{ "a036785", "A036785", "Numbers divisible by the squares of two distinct primes", "multiplicative",
  "36 different gaps occur, from 1 to 36; the level share is 19.20 %.",
  100000, 1, gen_a036785 },
{ "a023709", "A023709", "Numbers with no 1's in their base 4 expansion", "digit rule",
  "12 different gaps occur, from 1 to 1,048,577; the level share is 13.81 %; L = 2 holds 56 % of the level class; 12 terms do not decompose.",
  100000, 1, gen_a023709 },
{ "a023713", "A023713", "Numbers with no 2's in their base 4 expansion", "digit rule",
  "11 different gaps occur, from 1 to 262,145; the level share is 14.18 %; L = 1 holds 73 % of the level class; there are no ties; 11 terms do not decompose.",
  100000, 1, gen_a023713 },
{ "a023721", "A023721", "Numbers with no 0's in their base-5 expansion", "digit rule",
  "9 different gaps occur, from 1 to 97,657; the level share is 10.31 %; L = 1 holds 77 % of the level class.",
  100000, 1, gen_a023721 },
{ "a023725", "A023725", "Numbers with no 1's in their base-5 expansion", "digit rule",
  "10 different gaps occur, from 1 to 390,626; the level share is 12.99 %; L = 1 holds 78 % of the level class; 10 terms do not decompose.",
  100000, 1, gen_a023725 },
{ "a023729", "A023729", "Numbers with no 2's in their base-5 expansion", "digit rule",
  "9 different gaps occur, from 1 to 78,126; the level share is 13.34 %; L = 1 holds 78 % of the level class; 9 terms do not decompose.",
  100000, 1, gen_a023729 },
{ "a023733", "A023733", "Numbers with no 3's in base-5 expansion", "digit rule",
  "9 different gaps occur, from 1 to 78,126; the level share is 6.95 %; L = 1 holds 76 % of the level class.",
  100000, 1, gen_a023733 },
{ "a043493", "A043493", "Numbers that contain a single 1", "digit rule",
  "8 different gaps occur, from 1 to 10,001; the level share is 11.33 %; L = 1 holds 75 % of the level class.",
  100000, 1, gen_a043493 },
{ "a023692", "A023692", "Numbers with a single 1 in their ternary expansion", "digit rule",
  "13 different gaps occur, from 2 to 531,442; the level share is 23.50 %; L = 1 holds 86 % of the level class.",
  100000, 1, gen_a023692 },
{ "a020893", "A020893", "Squarefree sums of two squares; or squarefree numbers with no prime factors of the form 4k+3", "quadratic form",
  "41 different gaps occur, from 1 to 60; the level share is 17.28 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_a020893 },
{ "a030634", "A030634", "Numbers with 16 divisors", "divisor functions",
  "68 different gaps occur, from 1 to 97; the level share is 19.89 %; L = 2 holds 36 % of the level class.",
  100000, 1, gen_a030634 },
{ "a003625", "A003625", "Primes congruent to {3, 5, 6} mod 7", "primes",
  "106 different gaps occur, from 2 to 240; the level share is 23.69 %.",
  100000, 1, gen_a003625 },
{ "a030638", "A030638", "Numbers with 20 divisors", "divisor functions",
  "382 different gaps occur, from 1 to 638; the level share is 25.86 %.",
  100000, 1, gen_a030638 },
{ "a080682", "A080682", "19-smooth numbers: numbers whose prime divisors are all <= 19", "smooth",
  "22,333 different gaps occur, from 1 to 333,151; the level share is 55.60 %; 29.9 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a080682 },
{ "a080683", "A080683", "23-smooth numbers: numbers whose prime divisors are all <= 23", "smooth",
  "11,847 different gaps occur, from 1 to 69,696; the level share is 41.50 %; 7.4 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a080683 },
{ "a024913", "A024913", "Numbers k such that 10*k - 7 is prime", "prime values",
  "50 different gaps occur, from 1 to 54; the level share is 20.27 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a024913 },
{ "a037030", "A037030", "Numbers n such that 666*n + 1 is prime", "prime values",
  "52 different gaps occur, from 1 to 54; the level share is 19.47 %; L = 1 holds 43 % of the level class.",
  100000, 1, gen_a037030 },
{ "a073085", "A073085", "Numbers k such that 210*k+1 is prime", "prime values",
  "34 different gaps occur, from 1 to 35; the level share is 16.33 %; L = 1 holds 51 % of the level class.",
  100000, 1, gen_a073085 },
{ "a075745", "A075745", "Numbers n such that 210*n + 13 is prime", "prime values",
  "37 different gaps occur, from 1 to 45; the level share is 15.96 %; L = 1 holds 51 % of the level class.",
  100000, 1, gen_a075745 },
{ "a075746", "A075746", "Numbers n such that 210*n-13 is prime", "prime values",
  "33 different gaps occur, from 1 to 33; the level share is 16.60 %; L = 1 holds 52 % of the level class.",
  100000, 1, gen_a075746 },
{ "a075747", "A075747", "Numbers n such that 210*n + 17 is prime", "prime values",
  "34 different gaps occur, from 1 to 37; the level share is 16.36 %; L = 1 holds 52 % of the level class.",
  100000, 1, gen_a075747 },
{ "a075748", "A075748", "Numbers k such that 210*k-17 is prime", "prime values",
  "33 different gaps occur, from 1 to 34; the level share is 16.65 %; L = 1 holds 53 % of the level class.",
  100000, 1, gen_a075748 },
{ "a076354", "A076354", "Numbers n such that 210*n-1 is prime", "prime values",
  "35 different gaps occur, from 1 to 37; the level share is 16.50 %; L = 1 holds 52 % of the level class.",
  100000, 1, gen_a076354 },
{ "a076355", "A076355", "Numbers n such that 210*n + 11 is prime", "prime values",
  "34 different gaps occur, from 1 to 36; the level share is 16.75 %; L = 1 holds 52 % of the level class.",
  100000, 1, gen_a076355 },
{ "a076356", "A076356", "Numbers n such that 210*n-11 is prime", "prime values",
  "34 different gaps occur, from 1 to 34; the level share is 15.49 %; L = 1 holds 52 % of the level class.",
  100000, 1, gen_a076356 },
{ "a088958", "A088958", "Numbers n such that 60*n+1 is prime", "prime values",
  "36 different gaps occur, from 1 to 36; the level share is 17.03 %; L = 1 holds 49 % of the level class.",
  100000, 1, gen_a088958 },
{ "a090614", "A090614", "Numbers n such that 14n+3 is prime", "prime values",
  "53 different gaps occur, from 1 to 70; the level share is 16.42 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a090614 },
{ "a092022", "A092022", "Numbers k such that 16k + 3 is prime", "prime values",
  "64 different gaps occur, from 1 to 72; the level share is 16.89 %; L = 1 holds 32 % of the level class.",
  100000, 1, gen_a092022 },
{ "a101084", "A101084", "Numbers k such that 97*k + 101 is a prime", "prime values",
  "72 different gaps occur, from 2 to 162; the level share is 24.33 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a101084 },
{ "a101503", "A101503", "Numbers k such that 11*k + 101 is prime", "prime values",
  "56 different gaps occur, from 2 to 182; the level share is 22.49 %; L = 2 holds 45 % of the level class; there are no ties.",
  100000, 1, gen_a101503 },
{ "a101557", "A101557", "Numbers k such that 101*k + 1009 is prime", "prime values",
  "70 different gaps occur, from 2 to 156; the level share is 23.67 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a101557 },
{ "a102148", "A102148", "Numbers k such that 101*k + 11 is prime", "prime values",
  "72 different gaps occur, from 2 to 186; the level share is 23.39 %; L = 2 holds 41 % of the level class.",
  100000, 1, gen_a102148 },
{ "a102338", "A102338", "Numbers k such that 10k+3 is prime", "prime values",
  "50 different gaps occur, from 1 to 54; the level share is 15.91 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_a102338 },
{ "a102342", "A102342", "Numbers k such that 10k + 7 is prime", "prime values",
  "48 different gaps occur, from 1 to 50; the level share is 20.22 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a102342 },
{ "a102656", "A102656", "Numbers k such that 11*k + 1 is prime", "prime values",
  "59 different gaps occur, from 2 to 136; the level share is 23.07 %; L = 2 holds 44 % of the level class; there are no ties.",
  100000, 1, gen_a102656 },
{ "a014752", "A014752", "Primes of the form x^2 + 27y^2", "quadratic form",
  "138 different gaps occur, from 6 to 1,020; the level share is 42.30 %; L = 1 holds 40 % of the level class.",
  100000, 1, gen_a014752 },
{ "a033202", "A033202", "Primes of form x^2+93*y^2", "quadratic form",
  "84 different gaps occur, from 12 to 1,092; the level share is 44.04 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_a033202 },
{ "a033204", "A033204", "Primes of form x^2 + 94*y^2", "quadratic form",
  "692 different gaps occur, from 2 to 3,576; the level share is 42.53 %.",
  100000, 1, gen_a033204 },
{ "a033206", "A033206", "Primes of form x^2+95*y^2", "quadratic form",
  "569 different gaps occur, from 2 to 2,988; the level share is 47.51 %.",
  100000, 1, gen_a033206 },
{ "a033208", "A033208", "Primes of form x^2+97*y^2", "quadratic form",
  "249 different gaps occur, from 4 to 1,320; the level share is 37.98 %.",
  100000, 1, gen_a033208 },
{ "a033209", "A033209", "Primes of form x^2 + 11*y^2", "quadratic form",
  "348 different gaps occur, from 2 to 906; the level share is 37.94 %.",
  100000, 1, gen_a033209 },
{ "a033210", "A033210", "Primes of the form x^2+13*y^2", "quadratic form",
  "118 different gaps occur, from 4 to 764; the level share is 35.16 %.",
  100000, 1, gen_a033210 },
{ "a033211", "A033211", "Primes of form x^2 + 14*y^2", "quadratic form",
  "339 different gaps occur, from 2 to 1,398; the level share is 35.89 %.",
  100000, 1, gen_a033211 },
{ "a033213", "A033213", "Primes of form x^2+17*y^2", "quadratic form",
  "244 different gaps occur, from 4 to 1,212; the level share is 37.02 %.",
  100000, 1, gen_a033213 },
{ "a033214", "A033214", "Primes of form x^2+19*y^2", "quadratic form",
  "356 different gaps occur, from 2 to 1,036; the level share is 37.92 %.",
  100000, 1, gen_a033214 },
{ "a033215", "A033215", "Primes of form x^2+21*y^2", "quadratic form",
  "82 different gaps occur, from 12 to 1,092; the level share is 41.70 %.",
  100000, 1, gen_a033215 },
{ "a033216", "A033216", "Primes of form x^2+22*y^2", "quadratic form",
  "163 different gaps occur, from 2 to 560; the level share is 35.26 %.",
  100000, 1, gen_a033216 },
{ "a051645", "A051645", "Primes p such that 30*p+1 is also prime", "primes",
  "308 different gaps occur, from 2 to 988; the level share is 35.48 %.",
  100000, 1, gen_a051645 },
{ "a105961", "A105961", "Primes p such that 20*p + 3 is prime", "primes",
  "293 different gaps occur, from 2 to 840; the level share is 35.34 %.",
  100000, 1, gen_a105961 },
{ "a112391", "A112391", "Primes p such that 23*p + 2 is also prime", "primes",
  "284 different gaps occur, from 4 to 2,640; the level share is 48.21 %; L = 1 holds 36 % of the level class; 6 terms do not decompose.",
  100000, 1, gen_a112391 },
{ "a002496", "A002496", "Primes of the form k^2 + 1", "primes",
  "99,510 different gaps occur, from 3 to 645,426,204; every decomposable term is forced level (l <= d^2); 8 terms do not decompose.",
  100000, 1, gen_a002496 },
{ "a005473", "A005473", "Primes of form k^2 + 4", "primes",
  "98,965 different gaps occur, from 8 to 611,317,440; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a005473 },
{ "a028871", "A028871", "Primes of the form k^2 - 2", "primes",
  "98,323 different gaps occur, from 5 to 302,990,160; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a028871 },
{ "a028874", "A028874", "Primes of form k^2 - 3", "primes",
  "98,831 different gaps occur, from 36 to 745,165,488; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a028874 },
{ "a028877", "A028877", "Primes of form k^2 - 5", "primes",
  "98,931 different gaps occur, from 20 to 420,235,652; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a028877 },
{ "a028880", "A028880", "Primes of the form n^2 - 6", "primes",
  "98,477 different gaps occur, from 16 to 1,197,443,016; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a028880 },
{ "a028883", "A028883", "Primes of the form k^2 - 7", "primes",
  "99,273 different gaps occur, from 27 to 2,904,756,480; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a028883 },
{ "a028886", "A028886", "Primes of the form k^2 - 8", "primes",
  "98,428 different gaps occur, from 24 to 341,242,296; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a028886 },
{ "a049423", "A049423", "Primes of the form k^2 + 3", "primes",
  "98,756 different gaps occur, from 4 to 1,016,645,616; every decomposable term is forced level (l <= d^2); 9 terms do not decompose.",
  100000, 1, gen_a049423 },
{ "a056899", "A056899", "Primes of the form k^2 + 2", "primes",
  "98,771 different gaps occur, from 1 to 4,027,920,840; every decomposable term is forced level (l <= d^2); 9 terms do not decompose.",
  100000, 1, gen_a056899 },
{ "a056905", "A056905", "Primes of the form k^2 + 5", "primes",
  "98,835 different gaps occur, from 36 to 5,486,941,188; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a056905 },
{ "a056909", "A056909", "Primes of the form k^2+6", "primes",
  "98,620 different gaps occur, from 24 to 2,895,762,480; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a056909 },
{ "a079138", "A079138", "Primes of the form k^2 + 7", "primes",
  "99,212 different gaps occur, from 4 to 302,505,728; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a079138 },
{ "a091272", "A091272", "Primes of the form n^2 - 11", "primes",
  "99,431 different gaps occur, from 36 to 1,078,230,516; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a091272 },
{ "a027752", "A027752", "Numbers k such that k^2 + k + 3 is prime", "prime values",
  "82 different gaps occur, from 1 to 387; the level share is 32.70 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_a027752 },
{ "a027754", "A027754", "Numbers k such that k^2 + k + 5 is prime", "prime values",
  "120 different gaps occur, from 1 to 146; the level share is 24.57 %; L = 1 holds 32 % of the level class.",
  100000, 1, gen_a027754 },
{ "a027756", "A027756", "Numbers k such that k^2 + k + 7 is prime", "prime values",
  "95 different gaps occur, from 1 to 113; the level share is 24.14 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a027756 },
{ "a027757", "A027757", "Numbers k such that k^2 + k + 9 is prime", "prime values",
  "86 different gaps occur, from 3 to 288; the level share is 33.44 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a027757 },
{ "a028823", "A028823", "Numbers k such that k^2 + k + 17 is prime", "prime values",
  "53 different gaps occur, from 1 to 78; the level share is 18.69 %; L = 1 holds 45 % of the level class.",
  100000, 1, gen_a028823 },
{ "a045546", "A045546", "Numbers k such that k^2 + k - 1 is prime", "prime values",
  "62 different gaps occur, from 1 to 64; the level share is 20.21 %; L = 1 holds 42 % of the level class.",
  100000, 1, gen_a045546 },
{ "a048097", "A048097", "Numbers k such that k^2 + k + 11 is prime", "prime values",
  "68 different gaps occur, from 1 to 87; the level share is 19.96 %; L = 1 holds 40 % of the level class.",
  100000, 1, gen_a048097 },
{ "a014312", "A014312", "Numbers with exactly 4 ones in binary expansion", "binary rule",
  "139 different gaps occur, from 1 to 68,719,476,743; the level share is 12.94 %; 9.4 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a014312 },
{ "a014313", "A014313", "Numbers with exactly 5 ones in binary expansion", "binary rule",
  "106 different gaps occur, from 1 to 8,388,623; the level share is 8.89 %; 3.3 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a014313 },
{ "a023688", "A023688", "Numbers with exactly 6 ones in binary expansion", "binary rule",
  "84 different gaps occur, from 1 to 65,567; the level share is 8.93 %; L = 1 holds 48 % of the level class.",
  100000, 1, gen_a023688 },
{ "a023689", "A023689", "Numbers with exactly 7 ones in binary expansion", "binary rule",
  "71 different gaps occur, from 1 to 8,255; the level share is 9.38 %; L = 1 holds 56 % of the level class.",
  100000, 1, gen_a023689 },
{ "a023690", "A023690", "Numbers with exactly 8 ones in binary expansion", "binary rule",
  "61 different gaps occur, from 1 to 2,175; the level share is 11.36 %; L = 1 holds 64 % of the level class.",
  100000, 1, gen_a023690 },
{ "a023691", "A023691", "Numbers with exactly 9 ones in binary expansion", "binary rule",
  "55 different gaps occur, from 1 to 1,279; the level share is 12.25 %; L = 1 holds 63 % of the level class.",
  100000, 1, gen_a023691 },
{ "a005237", "A005237", "Numbers k such that k and k+1 have the same number of divisors", "divisor functions",
  "86 different gaps occur, from 1 to 105; the level share is 21.95 %; L = 1 holds 40 % of the level class.",
  100000, 1, gen_a005237 },
{ "a006049", "A006049", "Numbers k such that k and k+1 have the same number of distinct prime divisors", "multiplicative",
  "39 different gaps occur, from 1 to 45; the level share is 16.68 %; L = 1 holds 55 % of the level class.",
  100000, 1, gen_a006049 },
{ "a045920", "A045920", "Numbers m such that the factorizations of m..m+1 have the same number of primes (including multiplicities)", "multiplicative",
  "64 different gaps occur, from 1 to 71; the level share is 20.06 %; L = 1 holds 44 % of the level class.",
  100000, 1, gen_a045920 },
{ "a063464", "A063464", "Numbers k such that omega(k) = omega(k+2), where omega(k) is the number of distinct prime divisors of k", "multiplicative",
  "31 different gaps occur, from 1 to 37; the level share is 14.70 %; L = 1 holds 54 % of the level class.",
  100000, 1, gen_a063464 },
{ "a069977", "A069977", "Numbers k such that k and k+2 are squarefree", "multiplicative",
  "10 different gaps occur, from 2 to 24; the level share is 19.59 %; L = 1 holds 87 % of the level class.",
  100000, 1, gen_a069977 },
{ "a121495", "A121495", "Numbers k such that k and k+1 are composite and squarefree", "multiplicative",
  "19 different gaps occur, from 1 to 27; the level share is 16.18 %; L = 1 holds 57 % of the level class.",
  100000, 1, gen_a121495 },
{ "a124940", "A124940", "Numbers k such that k and k+3 are 3-almost primes", "multiplicative",
  "164 different gaps occur, from 1 to 211; the level share is 23.02 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_a124940 },
{ "a124941", "A124941", "Numbers k such that k and k+4 are 4-almost primes", "multiplicative",
  "256 different gaps occur, from 1 to 306; the level share is 26.20 %.",
  100000, 1, gen_a124941 },
{ "a070552", "A070552", "Semiprimes k such that k+1 is also a semiprime", "multiplicative",
  "311 different gaps occur, from 1 to 573; the level share is 31.89 %.",
  100000, 1, gen_a070552 },
{ "a031879", "A031879", "Nonprime lucky numbers", "sieve",
  "65 different gaps occur, from 2 to 144; the level share is 32.02 %; L = 1 holds 53 % of the level class.",
  100000, 1, gen_a031879 },
{ "a062324", "A062324", "Primes p such that p^2 + 4 is also prime", "primes",
  "532 different gaps occur, from 2 to 3,614; the level share is 45.85 %.",
  100000, 1, gen_a062324 },
{ "a062326", "A062326", "Primes p such that p^2 - 2 is also prime", "primes",
  "538 different gaps occur, from 1 to 1,618; the level share is 39.78 %.",
  100000, 1, gen_a062326 },
{ "a063637", "A063637", "Primes p such that p+2 is a semiprime", "primes",
  "220 different gaps occur, from 4 to 626; the level share is 36.04 %.",
  100000, 1, gen_a063637 },
{ "a063638", "A063638", "Primes p such that p-2 is a semiprime", "primes",
  "224 different gaps occur, from 4 to 572; the level share is 34.81 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a063638 },
{ "a127333", "A127333", "Numbers that are the sum of 6 consecutive primes", "primes",
  "102 different gaps occur, from 15 to 236; the level share is 34.32 %.",
  100000, 1, gen_a127333 },
{ "a127334", "A127334", "Numbers that are the sum of 7 consecutive primes", "primes",
  "110 different gaps occur, from 17 to 266; the level share is 44.22 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a127334 },
{ "a127336", "A127336", "Numbers that are the sum of 9 consecutive primes", "primes",
  "124 different gaps occur, from 27 to 300; the level share is 45.89 %.",
  100000, 1, gen_a127336 },
{ "a127337", "A127337", "Numbers that are the sum of 10 consecutive primes", "primes",
  "127 different gaps occur, from 29 to 306; the level share is 37.63 %.",
  100000, 1, gen_a127337 },
{ "a127338", "A127338", "Numbers that are the sum of 11 consecutive primes", "primes",
  "134 different gaps occur, from 35 to 318; the level share is 47.37 %.",
  100000, 1, gen_a127338 },
{ "a127339", "A127339", "Numbers that are the sum of 12 consecutive primes", "primes",
  "142 different gaps occur, from 39 to 330; the level share is 38.99 %.",
  100000, 1, gen_a127339 },
{ "a007950", "A007950", "Binary sieve: delete every 2nd number, then every 4th, 8th, etc", "sieve",
  "8 different gaps occur, from 2 to 16; the level share is 21.19 %; L = 1 holds 81 % of the level class.",
  100000, 1, gen_a007950 },
{ "a007951", "A007951", "Ternary sieve: delete every 3rd number, then every 9th, 27th, etc", "sieve",
  "8 different gaps occur, from 1 to 8; the level share is 3.51 %; L = 1 holds 62 % of the level class.",
  100000, 1, gen_a007951 },
{ "a059550", "A059550", "Beatty sequence for 1 + log(10)", "Beatty",
  "The gaps are 3 and 4; the level share is 17.17 %; L = 1 holds 50 % of the level class.",
  100000, 1, gen_a059550 },
{ "a059551", "A059551", "Beatty sequence for Gamma(1/3)", "Beatty",
  "The gaps are 2 and 3; the level share is 15.57 %; L = 1 holds 56 % of the level class.",
  100000, 1, gen_a059551 },
{ "a059552", "A059552", "Beatty sequence for Gamma(1/3)/(Gamma(1/3)-1)", "Beatty",
  "The gaps are 1 and 2; the level share is 12.15 %; L = 1 holds 76 % of the level class.",
  100000, 1, gen_a059552 },
{ "a059553", "A059553", "Beatty sequence for Gamma(2/3)", "Beatty",
  "The gaps are 1 and 2; the level share is 11.00 %; L = 1 holds 84 % of the level class.",
  100000, 1, gen_a059553 },
{ "a059554", "A059554", "Beatty sequence for Gamma(2/3)/(Gamma(2/3)-1)", "Beatty",
  "The gaps are 3 and 4; the level share is 18.11 %; L = 1 holds 47 % of the level class.",
  100000, 1, gen_a059554 },
{ "a059556", "A059556", "Beatty sequence for 1 + 1/gamma", "Beatty",
  "The gaps are 2 and 3; the level share is 15.75 %; L = 1 holds 56 % of the level class.",
  100000, 1, gen_a059556 },
{ "a059557", "A059557", "Beatty sequence for 1 + gamma^2, (gamma is the Euler-Mascheroni constant A001620)", "Beatty",
  "The gaps are 1 and 2; the level share is 11.05 %; L = 1 holds 85 % of the level class.",
  100000, 1, gen_a059557 },
{ "a059559", "A059559", "Beatty sequence for 1 + log(1/gamma), (gamma is the Euler-Mascheroni constant A001620)", "Beatty",
  "The gaps are 1 and 2; the level share is 11.93 %; L = 1 holds 77 % of the level class.",
  100000, 1, gen_a059559 },
{ "a059560", "A059560", "Beatty sequence for 1 - 1/log(gamma)", "Beatty",
  "The gaps are 2 and 3; the level share is 15.92 %; L = 1 holds 54 % of the level class.",
  100000, 1, gen_a059560 },
{ "a059562", "A059562", "Beatty sequence for log(Pi)/(log(Pi)-1)", "Beatty",
  "The gaps are 7 and 8; the level share is 23.26 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_a059562 },
{ "a059563", "A059563", "Beatty sequence for e + 1/e", "Beatty",
  "The gaps are 3 and 4; the level share is 16.78 %; L = 1 holds 52 % of the level class.",
  100000, 1, gen_a059563 },
{ "a059564", "A059564", "Beatty sequence for (e^2 + 1)/(e^2 - e + 1)", "Beatty",
  "The gaps are 1 and 2; the level share is 11.76 %; L = 1 holds 79 % of the level class.",
  100000, 1, gen_a059564 },
{ "a059566", "A059566", "Beatty sequence for e^gamma/(e^gamma-1)", "Beatty",
  "The gaps are 2 and 3; the level share is 14.51 %; L = 1 holds 61 % of the level class.",
  100000, 1, gen_a059566 },
{ "a059567", "A059567", "Beatty sequence for 1 - log(log(2))", "Beatty",
  "The gaps are 1 and 2; the level share is 11.29 %; L = 1 holds 84 % of the level class.",
  100000, 1, gen_a059567 },
{ "a059568", "A059568", "Beatty sequence for 1 - 1/log(log(2))", "Beatty",
  "The gaps are 3 and 4; the level share is 17.91 %; L = 1 holds 48 % of the level class.",
  100000, 1, gen_a059568 },
{ "a066343", "A066343", "Beatty sequence for log_2(10)", "Beatty",
  "The gaps are 3 and 4; the level share is 17.08 %; L = 1 holds 50 % of the level class.",
  100000, 1, gen_a066343 },
{ "a066344", "A066344", "Beatty sequence for log_5(10)", "Beatty",
  "The gaps are 1 and 2; the level share is 11.43 %; L = 1 holds 81 % of the level class.",
  100000, 1, gen_a066344 },
{ "a098005", "A098005", "Beatty sequence for 1/(3 - e): a(n) = floor(n/(3-e))", "Beatty",
  "The gaps are 3 and 4; the level share is 17.54 %; L = 1 holds 49 % of the level class.",
  100000, 1, gen_a098005 },
{ "a108598", "A108598", "a(n) = floor(n*((5+sqrt(5))/4))", "Beatty",
  "The gaps are 1 and 2; the level share is 12.98 %; L = 1 holds 70 % of the level class.",
  100000, 1, gen_a108598 },
{ "a121283", "A121283", "a(n) = floor(n*Pi*e)", "Beatty",
  "The gaps are 8 and 9; the level share is 23.74 %; L = 1 holds 34 % of the level class.",
  100000, 0, gen_a121283 },
{ "a102700", "A102700", "Numbers k such that 10*k + 9 is prime", "prime values",
  "47 different gaps occur, from 1 to 51; the level share is 14.62 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a102700 },
{ "a102703", "A102703", "Numbers k such that 100*k+99 is prime", "prime values",
  "57 different gaps occur, from 1 to 66; the level share is 16.08 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a102703 },
{ "a102711", "A102711", "Numbers k such that 11*k + 7 is prime", "prime values",
  "57 different gaps occur, from 2 to 142; the level share is 21.95 %; L = 2 holds 43 % of the level class; there are no ties.",
  100000, 1, gen_a102711 },
{ "a102721", "A102721", "Numbers n such that 11*n + 13 is prime", "prime values",
  "56 different gaps occur, from 2 to 182; the level share is 22.42 %; L = 2 holds 45 % of the level class.",
  100000, 1, gen_a102721 },
{ "a102731", "A102731", "Numbers k such that 11*k + 23 is prime", "prime values",
  "59 different gaps occur, from 2 to 136; the level share is 23.05 %; L = 2 holds 42 % of the level class.",
  100000, 1, gen_a102731 },
{ "a102768", "A102768", "Numbers k such that 23*k + 11 is prime", "prime values",
  "64 different gaps occur, from 2 to 162; the level share is 22.65 %; L = 2 holds 41 % of the level class; there are no ties.",
  100000, 1, gen_a102768 },
{ "a103118", "A103118", "Numbers k such that 100*k + 57 is prime", "prime values",
  "57 different gaps occur, from 1 to 72; the level share is 16.33 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a103118 },
{ "a103871", "A103871", "Numbers n such that 100n + 69 is prime", "prime values",
  "56 different gaps occur, from 1 to 70; the level share is 15.55 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_a103871 },
{ "a105043", "A105043", "Numbers n such that 100*n - 1 is prime", "prime values",
  "57 different gaps occur, from 1 to 66; the level share is 22.51 %; L = 1 holds 44 % of the level class.",
  100000, 1, gen_a105043 },
{ "a105044", "A105044", "Numbers n such that 1000*n - 1 is prime", "prime values",
  "63 different gaps occur, from 1 to 76; the level share is 22.89 %; L = 1 holds 41 % of the level class.",
  100000, 1, gen_a105044 },
{ "a105134", "A105134", "Numbers n such that 16n+9 is prime", "prime values",
  "65 different gaps occur, from 1 to 87; the level share is 18.59 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a105134 },
{ "a105135", "A105135", "Numbers n such that 32n+17 is prime", "prime values",
  "70 different gaps occur, from 1 to 75; the level share is 23.20 %; L = 1 holds 42 % of the level class.",
  100000, 1, gen_a105135 },
{ "a105136", "A105136", "Numbers n such that 64n+33 is prime", "prime values",
  "69 different gaps occur, from 1 to 76; the level share is 17.30 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a105136 },
{ "a105137", "A105137", "Numbers n such that 128n+65 is prime", "prime values",
  "69 different gaps occur, from 1 to 78; the level share is 21.42 %; L = 1 holds 39 % of the level class.",
  100000, 1, gen_a105137 },
{ "a105138", "A105138", "Numbers n such that 256n+129 is prime", "prime values",
  "77 different gaps occur, from 1 to 81; the level share is 18.98 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a105138 },
{ "a105139", "A105139", "Numbers k such that 512*k+257 is prime", "prime values",
  "77 different gaps occur, from 1 to 94; the level share is 24.31 %; L = 1 holds 39 % of the level class.",
  100000, 1, gen_a105139 },
{ "a105583", "A105583", "Numbers k such that 101*k + 997 is prime", "prime values",
  "75 different gaps occur, from 2 to 164; the level share is 24.04 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a105583 },
{ "a105679", "A105679", "Numbers k such that 997*k + 101 is prime", "prime values",
  "85 different gaps occur, from 2 to 200; the level share is 25.21 %; L = 2 holds 39 % of the level class.",
  100000, 1, gen_a105679 },
{ "a105773", "A105773", "Numbers n such that 11*n + 97 is prime", "prime values",
  "56 different gaps occur, from 2 to 128; the level share is 22.71 %; L = 2 holds 44 % of the level class; there are no ties.",
  100000, 1, gen_a105773 },
{ "a105775", "A105775", "Numbers n such that 97*n + 11 is prime", "prime values",
  "74 different gaps occur, from 2 to 174; the level share is 23.59 %; L = 2 holds 41 % of the level class.",
  100000, 1, gen_a105775 },
{ "a106690", "A106690", "Numbers k such that 11*k - 97 is prime", "prime values",
  "57 different gaps occur, from 1 to 182; the level share is 23.08 %; L = 2 holds 42 % of the level class.",
  100000, 1, gen_a106690 },
{ "a106692", "A106692", "Numbers k such that 97*k - 11 is prime", "prime values",
  "75 different gaps occur, from 2 to 172; the level share is 23.58 %; L = 2 holds 39 % of the level class.",
  100000, 1, gen_a106692 },
{ "a106695", "A106695", "Numbers k such that 101*k - 997 is prime", "prime values",
  "77 different gaps occur, from 2 to 216; the level share is 24.08 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a106695 },
{ "a106697", "A106697", "Numbers k such that 997*k - 101 is prime", "prime values",
  "85 different gaps occur, from 2 to 212; the level share is 25.04 %; L = 2 holds 38 % of the level class.",
  100000, 1, gen_a106697 },
{ "a107305", "A107305", "Numbers k such that 11*k - 13 is prime", "prime values",
  "56 different gaps occur, from 2 to 128; the level share is 22.41 %; L = 2 holds 43 % of the level class.",
  100000, 1, gen_a107305 },
{ "a107366", "A107366", "Numbers k such that 101*k + 103 is prime", "prime values",
  "78 different gaps occur, from 2 to 208; the level share is 24.29 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a107366 },
{ "a107369", "A107369", "Numbers n such that 103*n + 101 is prime", "prime values",
  "73 different gaps occur, from 2 to 180; the level share is 23.86 %; L = 2 holds 39 % of the level class.",
  100000, 1, gen_a107369 },
{ "a107371", "A107371", "Numbers k such that 101*k - 103 is prime", "prime values",
  "76 different gaps occur, from 2 to 192; the level share is 23.96 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a107371 },
{ "a107372", "A107372", "Numbers n such that 103*n - 101 is prime", "prime values",
  "71 different gaps occur, from 2 to 162; the level share is 24.05 %; L = 2 holds 41 % of the level class.",
  100000, 1, gen_a107372 },
{ "a107400", "A107400", "Numbers k such that 107*k + 109 is prime", "prime values",
  "74 different gaps occur, from 2 to 202; the level share is 23.76 %; L = 2 holds 39 % of the level class.",
  100000, 1, gen_a107400 },
{ "a107405", "A107405", "Numbers n such that 109*n + 107 is prime", "prime values",
  "72 different gaps occur, from 2 to 170; the level share is 24.30 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a107405 },
{ "a107406", "A107406", "Numbers n such that 107*n - 109 is prime", "prime values",
  "73 different gaps occur, from 2 to 192; the level share is 24.14 %; L = 2 holds 41 % of the level class; there are no ties.",
  100000, 1, gen_a107406 },
{ "a107407", "A107407", "Numbers n such that 109*n - 107 is prime", "prime values",
  "73 different gaps occur, from 2 to 158; the level share is 24.29 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a107407 },
{ "a107960", "A107960", "Numbers n such that 11*n - 1 is prime", "prime values",
  "59 different gaps occur, from 2 to 140; the level share is 22.57 %; L = 2 holds 43 % of the level class; there are no ties.",
  100000, 1, gen_a107960 },
{ "a107992", "A107992", "Numbers n such that 11*n - 3 is prime", "prime values",
  "57 different gaps occur, from 2 to 132; the level share is 16.12 %; L = 2 holds 32 % of the level class.",
  100000, 1, gen_a107992 },
{ "a107994", "A107994", "Numbers n such that 11*n - 2 is prime", "prime values",
  "56 different gaps occur, from 2 to 128; the level share is 32.73 %; L = 1 holds 59 % of the level class.",
  100000, 1, gen_a107994 },
{ "a108027", "A108027", "Numbers k such that 137*k + 139 is prime", "prime values",
  "76 different gaps occur, from 2 to 168; the level share is 24.16 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a108027 },
{ "a108028", "A108028", "Numbers k such that 139*k + 137 is prime", "prime values",
  "76 different gaps occur, from 2 to 168; the level share is 24.28 %; L = 2 holds 41 % of the level class; there are no ties.",
  100000, 1, gen_a108028 },
{ "a108029", "A108029", "Numbers k such that 149*k + 151 is prime", "prime values",
  "77 different gaps occur, from 2 to 180; the level share is 24.23 %; L = 2 holds 38 % of the level class.",
  100000, 1, gen_a108029 },
{ "a108030", "A108030", "Numbers k such that 151*k + 149 is prime", "prime values",
  "76 different gaps occur, from 2 to 168; the level share is 24.03 %; L = 2 holds 38 % of the level class.",
  100000, 1, gen_a108030 },
{ "a033217", "A033217", "Primes of form x^2 + 23*y^2", "quadratic form",
  "351 different gaps occur, from 2 to 968; the level share is 35.33 %.",
  100000, 1, gen_a033217 },
{ "a033218", "A033218", "Primes of form x^2+26*y^2", "quadratic form",
  "514 different gaps occur, from 2 to 2,328; the level share is 42.46 %.",
  100000, 1, gen_a033218 },
{ "a033219", "A033219", "Primes of form x^2+29*y^2", "quadratic form",
  "370 different gaps occur, from 4 to 2,016; the level share is 41.72 %.",
  100000, 1, gen_a033219 },
{ "a033220", "A033220", "Primes of form x^2+30*y^2", "quadratic form",
  "72 different gaps occur, from 18 to 1,062; the level share is 48.22 %; L = 1 holds 47 % of the level class.",
  100000, 1, gen_a033220 },
{ "a033221", "A033221", "Primes of form x^2+31*y^2", "quadratic form",
  "344 different gaps occur, from 2 to 1,092; the level share is 35.55 %.",
  100000, 1, gen_a033221 },
{ "a033222", "A033222", "Primes of form x^2+33*y^2", "quadratic form",
  "84 different gaps occur, from 12 to 1,188; the level share is 46.29 %; L = 1 holds 40 % of the level class.",
  100000, 1, gen_a033222 },
{ "a033223", "A033223", "Primes of form x^2+34*y^2", "quadratic form",
  "347 different gaps occur, from 2 to 1,320; the level share is 37.35 %.",
  100000, 1, gen_a033223 },
{ "a033224", "A033224", "Primes of form x^2+35*y^2", "quadratic form",
  "427 different gaps occur, from 2 to 2,010; the level share is 42.15 %.",
  100000, 1, gen_a033224 },
{ "a033225", "A033225", "Primes of form x^2+37*y^2", "quadratic form",
  "120 different gaps occur, from 4 to 564; the level share is 34.19 %.",
  100000, 1, gen_a033225 },
{ "a033226", "A033226", "Primes of form x^2+38*y^2", "quadratic form",
  "523 different gaps occur, from 2 to 1,880; the level share is 42.53 %.",
  100000, 1, gen_a033226 },
{ "a033227", "A033227", "Primes of form x^2+39*y^2", "quadratic form",
  "172 different gaps occur, from 6 to 1,296; the level share is 45.15 %; L = 1 holds 40 % of the level class.",
  100000, 1, gen_a033227 },
{ "a033228", "A033228", "Primes of form x^2+41*y^2", "quadratic form",
  "492 different gaps occur, from 4 to 2,824; the level share is 42.50 %.",
  100000, 1, gen_a033228 },
{ "a033229", "A033229", "Primes of form x^2+42*y^2", "quadratic form",
  "121 different gaps occur, from 6 to 1,302; the level share is 42.08 %.",
  100000, 1, gen_a033229 },
{ "a033230", "A033230", "Primes of form x^2+43*y^2", "quadratic form",
  "354 different gaps occur, from 2 to 980; the level share is 37.29 %.",
  100000, 1, gen_a033230 },
{ "a033231", "A033231", "Primes of form x^2+46*y^2", "quadratic form",
  "353 different gaps occur, from 2 to 1,584; the level share is 36.94 %.",
  100000, 1, gen_a033231 },
{ "a033232", "A033232", "Primes of form x^2+47*y^2", "quadratic form",
  "590 different gaps occur, from 2 to 1,720; the level share is 39.46 %.",
  100000, 1, gen_a033232 },
{ "a033233", "A033233", "Primes of form x^2+51*y^2", "quadratic form",
  "249 different gaps occur, from 6 to 2,244; the level share is 45.40 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a033233 },
{ "a033234", "A033234", "Primes of form x^2+53*y^2", "quadratic form",
  "369 different gaps occur, from 4 to 1,992; the level share is 41.47 %.",
  100000, 1, gen_a033234 },
{ "a033235", "A033235", "Primes of the form x^2 + 55*y^2", "quadratic form",
  "284 different gaps occur, from 2 to 1,320; the level share is 43.11 %.",
  100000, 1, gen_a033235 },
{ "a033236", "A033236", "Primes of form x^2+57*y^2", "quadratic form",
  "85 different gaps occur, from 12 to 1,212; the level share is 45.48 %; L = 1 holds 39 % of the level class.",
  100000, 1, gen_a033236 },
{ "a033237", "A033237", "Primes of form x^2+58*y^2", "quadratic form",
  "173 different gaps occur, from 2 to 640; the level share is 34.49 %.",
  100000, 1, gen_a033237 },
{ "a033238", "A033238", "Primes of form x^2+59*y^2", "quadratic form",
  "1,013 different gaps occur, from 2 to 3,426; the level share is 44.49 %.",
  100000, 1, gen_a033238 },
{ "a033239", "A033239", "Primes of form x^2+61*y^2", "quadratic form",
  "372 different gaps occur, from 4 to 2,556; the level share is 41.59 %.",
  100000, 1, gen_a033239 },
{ "a033240", "A033240", "Primes of form x^2+62*y^2", "quadratic form",
  "687 different gaps occur, from 2 to 2,862; the level share is 42.41 %.",
  100000, 1, gen_a033240 },
{ "a033241", "A033241", "Primes of form x^2+65*y^2", "quadratic form",
  "300 different gaps occur, from 8 to 3,172; the level share is 47.85 %.",
  100000, 1, gen_a033241 },
{ "a033242", "A033242", "Primes of form x^2+66*y^2", "quadratic form",
  "254 different gaps occur, from 6 to 2,442; the level share is 50.18 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a033242 },
{ "a033243", "A033243", "Primes of form x^2+67*y^2", "quadratic form",
  "347 different gaps occur, from 2 to 1,118; the level share is 37.05 %.",
  100000, 1, gen_a033243 },
{ "a033244", "A033244", "Primes of form x^2+69*y^2", "quadratic form",
  "183 different gaps occur, from 12 to 2,772; the level share is 47.27 %; L = 1 holds 32 % of the level class.",
  100000, 1, gen_a033244 },
{ "a033245", "A033245", "Primes of form x^2+70*y^2", "quadratic form",
  "202 different gaps occur, from 2 to 1,570; the level share is 39.09 %.",
  100000, 1, gen_a033245 },
{ "a033246", "A033246", "Primes of form x^2+71*y^2", "quadratic form",
  "811 different gaps occur, from 2 to 2,580; the level share is 41.65 %.",
  100000, 1, gen_a033246 },
{ "a033247", "A033247", "Primes of form x^2+73*y^2", "quadratic form",
  "244 different gaps occur, from 4 to 1,280; the level share is 37.84 %.",
  100000, 1, gen_a033247 },
{ "a033248", "A033248", "Primes of the form x^2+74*y^2", "quadratic form",
  "850 different gaps occur, from 2 to 3,822; the level share is 45.26 %.",
  100000, 1, gen_a033248 },
{ "a033249", "A033249", "Primes of form x^2+77*y^2", "quadratic form",
  "464 different gaps occur, from 4 to 2,920; the level share is 42.31 %.",
  100000, 1, gen_a033249 },
{ "a033250", "A033250", "Primes of form x^2+78*y^2", "quadratic form",
  "120 different gaps occur, from 6 to 1,122; the level share is 45.52 %; L = 1 holds 39 % of the level class.",
  100000, 1, gen_a033250 },
{ "a033251", "A033251", "Primes of form x^2+79*y^2", "quadratic form",
  "583 different gaps occur, from 2 to 1,794; the level share is 39.48 %.",
  100000, 1, gen_a033251 },
{ "a033252", "A033252", "Primes of form x^2+82*y^2", "quadratic form",
  "361 different gaps occur, from 2 to 1,442; the level share is 37.84 %.",
  100000, 1, gen_a033252 },
{ "a127989", "A127989", "a(n) = 2*n^3 - 2*n + 9", "polynomial",
  "Every gap is different, from 12 to 60,000,600,000; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a127989 },
{ "a004126", "A004126", "a(n) = n*(7*n^2 - 1)/6", "polynomial",
  "Every gap is different, from 1 to 34,999,650,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a004126 },
{ "a004188", "A004188", "a(n) = n*(3*n^2 - 1)/2", "polynomial",
  "Every gap is different, from 1 to 44,999,550,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a004188 },
{ "a004466", "A004466", "a(n) = n*(5*n^2 - 2)/3", "polynomial",
  "Every gap is different, from 1 to 49,999,500,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a004466 },
{ "a006597", "A006597", "a(n) = n^2*(5*n-3)/2", "polynomial",
  "Every gap is different, from 1 to 74,998,950,004; every decomposable term is forced level (l <= d^2); 8 terms do not decompose.",
  100000, 0, gen_a006597 },
{ "a061804", "A061804", "a(n) = 2*n*(2*n^2 + 1)", "polynomial",
  "Every gap is different, from 6 to 119,998,800,006; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a061804 },
{ "a063521", "A063521", "a(n) = n*(7*n^2-4)/3", "polynomial",
  "Every gap is different, from 1 to 69,999,300,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a063521 },
{ "a063522", "A063522", "a(n) = n*(5*n^2 - 3)/2", "polynomial",
  "Every gap is different, from 1 to 74,999,250,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a063522 },
{ "a063523", "A063523", "a(n) = n*(8*n^2 - 5)/3", "polynomial",
  "Every gap is different, from 1 to 79,999,200,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a063523 },
{ "a067389", "A067389", "a(n) = 3*n^3 + 2*n^2 + n", "polynomial",
  "Every gap is different, from 6 to 89,999,500,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a067389 },
{ "a117560", "A117560", "a(n) = n*(n^2 - 1)/2 - 1", "polynomial",
  "Every gap is different, from 9 to 15,000,450,003; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 2, gen_a117560 },
{ "a004467", "A004467", "a(n) = n*(11*n^2 - 5)/6", "polynomial",
  "Every gap is different, from 1 to 54,999,450,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a004467 },
{ "a060163", "A060163", "a(n) = (n^3 + 5*n + 18)/6", "polynomial",
  "99,998 different gaps occur, from 1 to 4,999,750,004; every decomposable term is forced level (l <= d^2); 8 terms do not decompose.",
  100000, -2, gen_a060163 },
{ "a062025", "A062025", "a(n) = n*(13*n^2 - 7)/6", "polynomial",
  "Every gap is different, from 1 to 64,999,350,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a062025 },
{ "a005586", "A005586", "a(n) = n*(n+4)*(n+5)/6", "polynomial",
  "Every gap is different, from 5 to 5,000,250,002; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a005586 },
{ "a006503", "A006503", "a(n) = n*(n+1)*(n+8)/6", "polynomial",
  "Every gap is different, from 3 to 5,000,250,000; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a006503 },
{ "a027480", "A027480", "a(n) = n*(n+1)*(n+2)/2", "polynomial",
  "Every gap is different, from 3 to 15,000,150,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a027480 },
{ "a027903", "A027903", "a(n) = n*(n + 1)*(3*n + 1)", "polynomial",
  "Every gap is different, from 8 to 89,999,900,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a027903 },
{ "a055112", "A055112", "a(n) = n*(n+1)*(2*n+1)", "polynomial",
  "Every gap is different, from 6 to 60,000,000,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a055112 },
{ "a059722", "A059722", "a(n) = n*(2*n^2 - 2*n + 1)", "polynomial",
  "Every gap is different, from 1 to 59,999,000,005; every decomposable term is forced level (l <= d^2); 8 terms do not decompose.",
  100000, 0, gen_a059722 },
{ "a071233", "A071233", "a(n) = 2*(n-1)*(n^2 + 1)", "polynomial",
  "Every gap is different, from 10 to 60,000,200,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a071233 },
{ "a077414", "A077414", "a(n) = n*(n - 1)*(n + 2)/2", "polynomial",
  "Every gap is different, from 4 to 15,000,250,000; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a077414 },
{ "a077415", "A077415", "a(n) = n*(n+2)*(n-2)/3", "polynomial",
  "Every gap is different, from 5 to 10,000,300,001; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 2, gen_a077415 },
{ "a084990", "A084990", "a(n) = n*(n^2+3*n-1)/3", "polynomial",
  "Every gap is different, from 1 to 10,000,099,999; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a084990 },
{ "a085786", "A085786", "a(n) = n*(2*n^2 + n + 1)/2", "polynomial",
  "Every gap is different, from 9 to 30,000,400,002; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a085786 },
{ "a090197", "A090197", "a(n) = n^3 + 6*n^2 + 6*n + 1", "polynomial",
  "Every gap is different, from 13 to 30,000,900,001; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a090197 },
{ "a094421", "A094421", "a(n) = n * (6*n^2 + 6*n + 1)", "polynomial",
  "Every gap is different, from 61 to 180,003,000,013; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a094421 },
{ "a110451", "A110451", "a(n) = n*(4*n^2 + 2*n + 1)", "polynomial",
  "Every gap is different, from 7 to 119,999,200,003; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a110451 },
{ "a111396", "A111396", "a(n) = n*(n+7)*(n+8)/6", "polynomial",
  "Every gap is different, from 12 to 5,000,450,007; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a111396 },
{ "a125200", "A125200", "a(n) = n*(4*n^2 + n - 1)/2", "polynomial",
  "Every gap is different, from 15 to 60,000,700,002; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a125200 },
{ "a008607", "A008607", "Multiples of 25", "residue class",
  "The gap is always 25; the level share is 9.63 %; L = 25 holds 100 % of the level class.",
  100000, 0, gen_a008607 },
{ "a008854", "A008854", "Numbers that are congruent to {0, 1, 4} mod 5", "residue class",
  "The gaps are 1 and 3; the level share is 14.86 %; L = 1 holds 77 % of the level class.",
  100000, 1, gen_a008854 },
{ "a029739", "A029739", "Numbers that are congruent to {1, 3, 4} mod 6", "residue class",
  "The gaps are 1, 2 and 3; the level share is 17.98 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a029739 },
{ "a032769", "A032769", "Numbers that are congruent to {0, 1, 2, 4} mod 5", "residue class",
  "The gaps are 1 and 2; the level share is 5.87 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a032769 },
{ "a032775", "A032775", "Numbers that are congruent to {0, 1, 2, 3, 5, 6} mod 7", "residue class",
  "The gaps are 1 and 2; the level share is 10.17 %; L = 1 holds 90 % of the level class.",
  100000, 1, gen_a032775 },
{ "a032793", "A032793", "Numbers that are congruent to {1, 2, 4} mod 5", "residue class",
  "The gaps are 1 and 2; the level share is 5.83 %; L = 1 holds 65 % of the level class.",
  100000, 1, gen_a032793 },
{ "a032796", "A032796", "Numbers that are congruent to {1, 2, 3, 5, 6} mod 7", "residue class",
  "The gaps are 1 and 2; the level share is 10.95 %; L = 1 holds 79 % of the level class.",
  100000, 1, gen_a032796 },
{ "a044102", "A044102", "Multiples of 36", "residue class",
  "The gap is always 36; the level share is 9.63 %; L = 36 holds 100 % of the level class.",
  100000, 0, gen_a044102 },
{ "a047202", "A047202", "Numbers that are congruent to {2, 3, 4} mod 5", "residue class",
  "The gaps are 1 and 3; the level share is 14.87 %; L = 1 holds 77 % of the level class.",
  100000, 1, gen_a047202 },
{ "a047204", "A047204", "Numbers that are congruent to {3, 4} mod 5", "residue class",
  "The gaps are 1 and 4; the level share is 5.52 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a047204 },
{ "a003628", "A003628", "Primes congruent to {5, 7} mod 8", "primes",
  "83 different gaps occur, from 2 to 248; the level share is 27.96 %.",
  100000, 1, gen_a003628 },
{ "a033212", "A033212", "Primes congruent to 1 or 19 (mod 30)", "primes",
  "46 different gaps occur, from 12 to 660; the level share is 43.83 %; L = 1 holds 54 % of the level class.",
  100000, 1, gen_a033212 },
{ "a039949", "A039949", "Primes of the form 30n - 13", "primes",
  "33 different gaps occur, from 30 to 1,050; the level share is 48.53 %; L = 1 holds 46 % of the level class; there are no ties.",
  100000, 1, gen_a039949 },
{ "a042987", "A042987", "Primes congruent to {2, 3, 5, 7} mod 8", "primes",
  "74 different gaps occur, from 1 to 164; the level share is 25.07 %.",
  100000, 1, gen_a042987 },
{ "a042988", "A042988", "Primes not congruent to -1 (mod 7)", "primes",
  "66 different gaps occur, from 1 to 140; the level share is 25.73 %; L = 1 holds 31 % of the level class.",
  100000, 0, gen_a042988 },
{ "a042989", "A042989", "Primes congruent to {0, 2, 3, 4, 5} mod 7", "primes",
  "82 different gaps occur, from 1 to 190; the level share is 27.35 %.",
  100000, 1, gen_a042989 },
{ "a042990", "A042990", "Primes not congruent to 4 (mod 7)", "primes",
  "65 different gaps occur, from 1 to 132; the level share is 24.31 %.",
  100000, 0, gen_a042990 },
{ "a042992", "A042992", "Primes congruent to {0, 2, 3, 5, 6} (mod 7)", "primes",
  "82 different gaps occur, from 1 to 198; the level share is 24.76 %; L = 3 holds 30 % of the level class.",
  100000, 1, gen_a042992 },
{ "a042994", "A042994", "Primes congruent to {0, 1, 2, 3, 5} (mod 7)", "primes",
  "82 different gaps occur, from 1 to 200; the level share is 27.68 %.",
  100000, 1, gen_a042994 },
{ "a042995", "A042995", "Primes congruent to {0, 2, 3, 5} (mod 7)", "primes",
  "112 different gaps occur, from 1 to 294; the level share is 29.52 %.",
  100000, 1, gen_a042995 },
{ "a042997", "A042997", "Primes congruent to {2, 3, 4, 5, 6} (mod 7)", "primes",
  "65 different gaps occur, from 1 to 142; the level share is 24.12 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a042997 },
{ "a042998", "A042998", "Primes congruent to {1, 2, 3, 5} (mod 8)", "primes",
  "71 different gaps occur, from 1 to 160; the level share is 25.03 %.",
  100000, 1, gen_a042998 },
{ "a045320", "A045320", "Primes not congruent to 5 (mod 7)", "primes",
  "66 different gaps occur, from 1 to 136; the level share is 24.96 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_a045320 },
{ "a045322", "A045322", "Primes congruent to {0, 2, 3, 4, 6} (mod 7)", "primes",
  "83 different gaps occur, from 1 to 184; the level share is 26.85 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a045322 },
{ "a045323", "A045323", "Primes congruent to {1, 2, 3, 7} (mod 8)", "primes",
  "76 different gaps occur, from 1 to 174; the level share is 25.34 %.",
  100000, 1, gen_a045323 },
{ "a001840", "A001840", "Expansion of g.f. x/((1 - x)^2*(1 - x^3))", "polynomial",
  "33,334 different gaps occur, from 1 to 33,334; the level share is 65.74 %.",
  100000, 0, gen_a001840 },
{ "a001859", "A001859", "Triangular numbers plus quarter-squares: n*(n+1)/2 + floor((n+1)^2/4) (i.e., A000217(n) + A002620(n+1))", "polynomial",
  "Every gap is different, from 2 to 150,000; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a001859 },
{ "a002492", "A002492", "Sum of the first n even squares: a(n) = 2*n*(n+1)*(2*n+1)/3", "polynomial",
  "Every gap is different, from 4 to 40,000,000,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a002492 },
{ "a002623", "A002623", "Expansion of 1/((1-x)^4*(1+x))", "polynomial",
  "Every gap is different, from 2 to 2,500,100,001; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a002623 },
{ "a002717", "A002717", "a(n) = floor(n(n+2)(2n+1)/8)", "polynomial",
  "Every gap is different, from 1 to 7,500,050,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a002717 },
{ "a003600", "A003600", "Maximal number of pieces obtained by slicing a torus (or a bagel) with n cuts: (n^3 + 3*n^2 + 8*n)/6 (n > 0)", "polynomial",
  "Every gap is different, from 1 to 5,000,050,001; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a003600 },
{ "a004006", "A004006", "a(n) = C(n,1) + C(n,2) + C(n,3), or n*(n^2 + 5)/6", "polynomial",
  "Every gap is different, from 1 to 4,999,950,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a004006 },
{ "a005286", "A005286", "a(n) = (n + 3)*(n^2 + 6*n + 2)/6", "polynomial",
  "Every gap is different, from 5 to 5,000,250,002; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a005286 },
{ "a005744", "A005744", "Expansion of x*(1+x-x^2)/((1-x)^4*(1+x))", "polynomial",
  "Every gap is different, from 1 to 2,500,100,000; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a005744 },
{ "a005893", "A005893", "Number of points on surface of tetrahedron; coordination sequence for sodalite net (equals 2*n^2+2 for n > 0)", "polynomial",
  "Every gap is different, from 3 to 399,998; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a005893 },
{ "a005897", "A005897", "a(n) = 6*n^2 + 2 for n > 0, a(0)=1", "polynomial",
  "Every gap is different, from 7 to 1,199,994; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a005897 },
{ "a005899", "A005899", "Number of points on surface of octahedron; also coordination sequence for cubic lattice: a(0) = 1; for n > 0, a(n) = 4n^2 + 2", "polynomial",
  "Every gap is different, from 5 to 799,996; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a005899 },
{ "a005901", "A005901", "Number of points on surface of cuboctahedron (or icosahedron): a(0) = 1; for n > 0, a(n) = 10n^2 + 2. Also coordination sequence for f.c.c. or A_3 or D_3 lattice", "polynomial",
  "Every gap is different, from 11 to 1,999,990; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a005901 },
{ "a005914", "A005914", "Number of points on surface of hexagonal prism: 12*n^2 + 2 for n > 0 (coordination sequence for W(2))", "polynomial",
  "Every gap is different, from 13 to 2,399,988; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a005914 },
{ "a005920", "A005920", "Tricapped prism numbers", "polynomial",
  "Every gap is different, from 8 to 45,000,250,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a005920 },
{ "a005993", "A005993", "Expansion of (1+x^2)/((1-x)^2*(1-x^2)^2)", "polynomial",
  "50,001 different gaps occur, from 1 to 2,500,100,001; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a005993 },
{ "a006004", "A006004", "a(n) = C(n+2,3) + C(n,3) + C(n-1,3)", "polynomial",
  "Every gap is different, from 3 to 14,999,950,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a006004 },
{ "a006918", "A006918", "a(n) = binomial(n+3, 3)/4 for odd n, n*(n+2)*(n+4)/24 for even n", "polynomial",
  "50,000 different gaps occur, from 1 to 1,250,025,000; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a006918 },
{ "a007202", "A007202", "Crystal ball sequence for hexagonal close-packing", "polynomial",
  "Every gap is different, from 12 to 105,000,000,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a007202 },
{ "a007586", "A007586", "11-gonal (or hendecagonal) pyramidal numbers: a(n) = n*(n+1)*(3*n-2)/2", "polynomial",
  "Every gap is different, from 1 to 44,999,650,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a007586 },
{ "a007587", "A007587", "12-gonal (or dodecagonal) pyramidal numbers: a(n) = n*(n+1)*(10*n-7)/6", "polynomial",
  "Every gap is different, from 1 to 49,999,600,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a007587 },
{ "a007904", "A007904", "Crystal ball sequence for diamond", "polynomial",
  "Every gap is different, from 4 to 25,000,000,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a007904 },
{ "a008412", "A008412", "Coordination sequence for 4-dimensional cubic lattice (points on surface of 4-dimensional cross-polytope)", "polynomial",
  "Every gap is different, from 7 to 79,999,200,008; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a008412 },
{ "a008577", "A008577", "Crystal ball sequence for planar net 4.8.8", "polynomial",
  "Every gap is different, from 3 to 266,667; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a008577 },
{ "a008580", "A008580", "Crystal ball sequence for planar net 3.6.3.6", "polynomial",
  "90,000 different gaps occur, from 4 to 499,998; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a008580 },
{ "a008804", "A008804", "Expansion of 1/((1-x)^2*(1-x^2)*(1-x^4))", "polynomial",
  "50,001 different gaps occur, from 1 to 625,050,001; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a008804 },
{ "a008810", "A008810", "a(n) = ceiling(n^2/3)", "polynomial",
  "33,334 different gaps occur, from 1 to 66,667; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a008810 },
{ "a019298", "A019298", "Number of balls in pyramid with base either a regular hexagon or a hexagon with alternate sides differing by 1 (balls in hexagonal pyramid of height n taken from hexagonal close-packing)", "polynomial",
  "Every gap is different, from 1 to 7,500,000,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a019298 },
{ "a024206", "A024206", "Expansion of x^2*(1+x-x^2)/((1-x^2)*(1-x)^2)", "polynomial",
  "50,001 different gaps occur, from 1 to 50,001; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a024206 },
{ "a028552", "A028552", "a(n) = n*(n+3)", "polynomial",
  "Every gap is different, from 4 to 200,002; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a028552 },
{ "a033586", "A033586", "a(n) = 4*n*(2*n + 1)", "polynomial",
  "Every gap is different, from 12 to 1,599,996; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a033586 },
{ "a035005", "A035005", "Number of possible queen moves on an n X n chessboard", "polynomial",
  "Every gap is different, from 12 to 100,000,200,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a035005 },
{ "a035006", "A035006", "Number of possible rook moves on an n X n chessboard", "polynomial",
  "Every gap is different, from 8 to 60,000,200,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a035006 },
{ "a035008", "A035008", "Total number of possible knight moves on an (n+2) X (n+2) chessboard, if the knight is placed anywhere", "polynomial",
  "Every gap is different, from 16 to 1,600,000; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a035008 },
{ "a045944", "A045944", "Rhombic matchstick numbers: a(n) = n*(3*n+2)", "polynomial",
  "Every gap is different, from 5 to 599,999; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a045944 },
{ "a045946", "A045946", "Star of David matchstick numbers: a(n) = 6*n*(3*n+1)", "polynomial",
  "Every gap is different, from 24 to 3,599,988; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a045946 },
{ "a001740", "A001740", "Squares written in base 5", "digit rule",
  "Every gap is different, from 1 to 55,555,569,999,999; every decomposable term is forced level (l <= d^2); 18 terms do not decompose.",
  100000, 0, gen_a001740 },
{ "a001741", "A001741", "Squares written in base 6", "digit rule",
  "Every gap is different, from 1 to 444,447,999,999; every decomposable term is forced level (l <= d^2); 16 terms do not decompose.",
  100000, 0, gen_a001741 },
{ "a002440", "A002440", "Squares written in base 7", "digit rule",
  "99,579 different gaps occur, from 3 to 33,334,673,699; every decomposable term is forced level (l <= d^2); 15 terms do not decompose.",
  100000, 1, gen_a002440 },
{ "a002441", "A002441", "Squares written in base 8", "digit rule",
  "87,767 different gaps occur, from 3 to 22,222,772,023; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a002441 },
{ "a002442", "A002442", "Squares written in base 9", "digit rule",
  "87,369 different gaps occur, from 3 to 1,111,384,093; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a002442 },
{ "a001363", "A001363", "Primes in ternary", "digit rule",
  "508 different gaps occur, from 2 to 777,777,779,909; the level share is 27.92 %; 1.3 % of terms are forced level (l <= d^2); there are no ties; 24 terms do not decompose.",
  100000, 1, gen_a001363 },
{ "a004678", "A004678", "Primes written in base 4", "digit rule",
  "470 different gaps occur, from 1 to 6,666,666,682; the level share is 27.27 %; 1.3 % of terms are forced level (l <= d^2); 20 terms do not decompose.",
  100000, 1, gen_a004678 },
{ "a004679", "A004679", "Primes written in base 5", "digit rule",
  "376 different gaps occur, from 1 to 55,555,731; the level share is 27.61 %; 10 terms do not decompose.",
  100000, 1, gen_a004679 },
{ "a004680", "A004680", "Primes written in base 6", "digit rule",
  "238 different gaps occur, from 1 to 4,444,510; the level share is 30.66 %; L = 1 holds 38 % of the level class; 9 terms do not decompose.",
  100000, 1, gen_a004680 },
{ "a004681", "A004681", "Primes written in base 7", "digit rule",
  "305 different gaps occur, from 1 to 3,333,339; the level share is 23.62 %; 9 terms do not decompose.",
  100000, 1, gen_a004681 },
{ "a004682", "A004682", "Primes written in base 8", "digit rule",
  "201 different gaps occur, from 1 to 222,248; the level share is 30.81 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a004682 },
{ "a004683", "A004683", "Primes written in base 9", "digit rule",
  "228 different gaps occur, from 1 to 111,192; the level share is 24.56 %; L = 1 holds 30 % of the level class.",
  100000, 1, gen_a004683 },
{ "a029954", "A029954", "Palindromic in base 7", "digit rule",
  "12 different gaps occur, from 1 to 19,208; the level share is 77.07 %; 15.1 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a029954 },
{ "a029803", "A029803", "Numbers that are palindromic in base 8", "digit rule",
  "12 different gaps occur, from 1 to 36,864; the level share is 84.85 %; 29.8 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a029803 },
{ "a029955", "A029955", "Palindromic in base 9", "digit rule",
  "11 different gaps occur, from 1 to 65,610; the level share is 79.55 %; 36.6 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a029955 },
{ "a029956", "A029956", "Numbers that are palindromic in base 11", "digit rule",
  "10 different gaps occur, from 1 to 15,972; the level share is 75.41 %; 13.6 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a029956 },
{ "a029957", "A029957", "Numbers that are palindromic in base 12", "digit rule",
  "10 different gaps occur, from 1 to 22,464; the level share is 79.45 %; 19.3 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a029957 },
{ "a029958", "A029958", "Numbers that are palindromic in base 13", "digit rule",
  "10 different gaps occur, from 1 to 30,758; the level share is 80.81 %; 26.7 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a029958 },
{ "a029959", "A029959", "Numbers that are palindromic in base 14", "digit rule",
  "10 different gaps occur, from 1 to 41,160; the level share is 83.44 %; 36.1 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a029959 },
{ "a029960", "A029960", "Numbers that are palindromic in base 15", "digit rule",
  "9 different gaps occur, from 1 to 54,000; the level share is 79.98 %; 46.1 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a029960 },
{ "a029730", "A029730", "Numbers that are palindromic in base 16", "digit rule",
  "9 different gaps occur, from 1 to 69,632; the level share is 76.24 %; 32.3 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a029730 },
{ "a023699", "A023699", "Numbers with a single 2 in their ternary expansion", "digit rule",
  "25 different gaps occur, from 1 to 265,723; the level share is 16.83 %; L = 1 holds 64 % of the level class.",
  100000, 1, gen_a023699 },
{ "a023706", "A023706", "Numbers with a single 0 in their base 4 expansion", "digit rule",
  "17 different gaps occur, from 1 to 21,849; the level share is 12.45 %; L = 1 holds 80 % of the level class.",
  100000, 1, gen_a023706 },
{ "a023710", "A023710", "Numbers with a single 1 in their base 4 expansion", "digit rule",
  "12 different gaps occur, from 1 to 65,537; the level share is 15.88 %; L = 2 holds 53 % of the level class.",
  100000, 1, gen_a023710 },
{ "a023714", "A023714", "Numbers with a single 2 in their base 4 expansion", "digit rule",
  "11 different gaps occur, from 1 to 16,385; the level share is 13.86 %; L = 1 holds 72 % of the level class; there are no ties.",
  100000, 1, gen_a023714 },
{ "a023718", "A023718", "Numbers with a single 3 in their base 4 expansion", "digit rule",
  "17 different gaps occur, from 1 to 21,849; the level share is 9.88 %; L = 1 holds 71 % of the level class.",
  100000, 1, gen_a023718 },
{ "a023722", "A023722", "Numbers with a single 0 in their base 5 expansion", "digit rule",
  "13 different gaps occur, from 1 to 3,911; the level share is 10.58 %; L = 1 holds 70 % of the level class.",
  100000, 1, gen_a023722 },
{ "a023726", "A023726", "Numbers with a single 1 in their base 5 expansion", "digit rule",
  "10 different gaps occur, from 1 to 15,626; the level share is 15.09 %; L = 1 holds 73 % of the level class.",
  100000, 1, gen_a023726 },
{ "a023730", "A023730", "Numbers with a single 2 in their base 5 expansion", "digit rule",
  "10 different gaps occur, from 1 to 15,626; the level share is 14.96 %; L = 1 holds 72 % of the level class.",
  100000, 1, gen_a023730 },
{ "a023734", "A023734", "Numbers with a single 3 in their base-5 expansion", "digit rule",
  "10 different gaps occur, from 1 to 15,626; the level share is 8.15 %; L = 1 holds 70 % of the level class.",
  100000, 1, gen_a023734 },
{ "a023738", "A023738", "Numbers with a single 4 in their base 5 expansion", "digit rule",
  "13 different gaps occur, from 1 to 3,911; the level share is 12.38 %; L = 1 holds 70 % of the level class.",
  100000, 1, gen_a023738 },
{ "a001729", "A001729", "List of numbers whose digits contain no loops (version 1)", "digit rule",
  "13 different gaps occur, from 1 to 333,334; the level share is 14.96 %; L = 1 holds 86 % of the level class.",
  100000, 1, gen_a001729 },
{ "a001742", "A001742", "Numbers whose digits contain no loops (version 2)", "digit rule",
  "15 different gaps occur, from 1 to 3,333,334; the level share is 17.14 %; L = 1 holds 82 % of the level class.",
  100000, 1, gen_a001742 },
{ "a001743", "A001743", "Numbers in which every digit contains at least one loop (version 1)", "digit rule",
  "18 different gaps occur, from 1 to 500,000,001; the level share is 13.48 %; L = 2 holds 43 % of the level class; there are no ties; 9 terms do not decompose.",
  100000, 1, gen_a001743 },
{ "a001744", "A001744", "Numbers n such that every digit contains a loop (version 2)", "digit rule",
  "16 different gaps occur, from 1 to 30,000,001; the level share is 13.73 %; L = 2 holds 48 % of the level class; 9 terms do not decompose.",
  100000, 1, gen_a001744 },
{ "a001745", "A001745", "Numbers such that at least one digit contains a loop (version 2). Also called \"holey\" or \"holy\" numbers", "digit rule",
  "The gaps are 1, 2 and 4; the level share is 9.76 %; L = 1 holds 95 % of the level class.",
  100000, 1, gen_a001745 },
{ "a001746", "A001746", "At least one digit contains a loop (version 1)", "digit rule",
  "The gaps are 1, 2 and 6; the level share is 9.87 %; L = 1 holds 94 % of the level class.",
  100000, 1, gen_a001746 },
{ "a046030", "A046030", "Numbers whose digits are squares", "digit rule",
  "17 different gaps occur, from 1 to 40,000,001; the level share is 15.56 %; L = 1 holds 46 % of the level class; 17 terms do not decompose.",
  100000, 1, gen_a046030 },
{ "a046034", "A046034", "Numbers whose digits are primes", "digit rule",
  "17 different gaps occur, from 1 to 144,444,445; the level share is 17.33 %; L = 1 holds 65 % of the level class; 10 terms do not decompose.",
  100000, 1, gen_a046034 },
{ "a028846", "A028846", "Numbers whose product of digits is a power of 2", "digit rule",
  "25 different gaps occur, from 1 to 32,222,223; the level share is 5.16 %; L = 5 holds 62 % of the level class; there are no ties; 10 terms do not decompose.",
  100000, 1, gen_a028846 },
{ "a000787", "A000787", "Strobogrammatic numbers: the same upside down", "digit rule",
  "209 different gaps occur, from 1 to 40,000,003,333,348; the level share is 93.39 %; 68.8 % of terms are forced level (l <= d^2); there are no ties; 17 terms do not decompose.",
  100000, 1, gen_strobo },
{ "a006507", "A006507", "a(n+1) = a(n) + sum of digits of a(n), with a(1)=7", "digit rule",
  "34 different gaps occur, from 1 to 50; the level share is 14.19 %; L = 9 holds 60 % of the level class.",
  100000, 1, gen_a006507 },
{ "a007618", "A007618", "a(n) = a(n-1) + sum of digits of a(n-1), a(1) = 5", "digit rule",
  "34 different gaps occur, from 1 to 50; the level share is 14.20 %; L = 9 holds 60 % of the level class.",
  100000, 1, gen_a007618 },
{ "a010062", "A010062", "a(0)=1; thereafter a(n+1) = a(n) + number of 1's in binary representation of a(n)", "digit rule",
  "18 different gaps occur, from 1 to 19; the level share is 23.93 %; L = 1 holds 33 % of the level class.",
  100000, 0, gen_a010062 },
{ "a010063", "A010063", "a(n+1) = a(n) + sum of digits in base 3 representation of a(n), with a(0) = 1", "digit rule",
  "12 different gaps occur, from 1 to 22; the level share is 22.31 %; L = 2 holds 40 % of the level class.",
  100000, 0, gen_a010063 },
{ "a010065", "A010065", "a(n+1) = a(n) + sum of digits in base 4 representation of a(n), with a(0) = 1", "digit rule",
  "18 different gaps occur, from 1 to 26; the level share is 18.35 %; L = 3 holds 44 % of the level class.",
  100000, 0, gen_a010065 },
{ "a010066", "A010066", "a(n+1) = a(n) + sum of digits in base 5 representation of a(n)", "digit rule",
  "10 different gaps occur, from 1 to 32; the level share is 19.94 %; L = 4 holds 48 % of the level class.",
  100000, 0, gen_a010066 },
{ "a010068", "A010068", "a(n+1) = a(n) + sum of digits in base 6 representation of a(n)", "digit rule",
  "31 different gaps occur, from 1 to 38; the level share is 16.91 %; L = 5 holds 56 % of the level class.",
  100000, 0, gen_a010068 },
{ "a010069", "A010069", "a(n+1) = a(n) + sum of digits in base 7 representation of a(n)", "digit rule",
  "15 different gaps occur, from 1 to 40; the level share is 15.89 %; L = 6 holds 53 % of the level class.",
  100000, 0, gen_a010069 },
{ "a010071", "A010071", "a(n+1) = a(n) + sum of digits in base 8 representation of a(n)", "digit rule",
  "20 different gaps occur, from 1 to 44; the level share is 15.94 %; L = 7 holds 51 % of the level class.",
  100000, 0, gen_a010071 },
{ "a010072", "A010072", "a(n+1) = a(n) + sum of digits in base 9 representation of a(n)", "digit rule",
  "9 different gaps occur, from 1 to 48; the level share is 17.84 %; L = 8 holds 53 % of the level class.",
  100000, 0, gen_a010072 },
{ "a010064", "A010064", "Base 4 self or Colombian numbers (not of form k + sum of base 4 digits of k)", "digit rule",
  "8 different gaps occur, from 2 to 35; the level share is 19.18 %; L = 1 holds 44 % of the level class.",
  100000, 1, gen_a010064 },
{ "a010067", "A010067", "Base 6 self or Colombian numbers (not of form k + sum of base 6 digits of k)", "digit rule",
  "7 different gaps occur, from 2 to 47; the level share is 20.78 %; L = 1 holds 43 % of the level class.",
  100000, 1, gen_a010067 },
{ "a010070", "A010070", "Base 8 self or Colombian numbers (not of form k + sum of base 8 digits of k)", "digit rule",
  "6 different gaps occur, from 2 to 46; the level share is 22.88 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a010070 },
{ "a001704", "A001704", "a(n) = n concatenated with n + 1", "digit rule",
  "11 different gaps occur, from 11 to 89,999,200,001; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a001704 },
{ "a019550", "A019550", "a(n) is the concatenation of n and 2n", "digit rule",
  "11 different gaps occur, from 12 to 45,000,100,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a019550 },
{ "a019553", "A019553", "a(n) is the concatenation of n and 5n", "digit rule",
  "10 different gaps occur, from 105 to 18,000,100,005; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a019553 },
{ "a009440", "A009440", "a(n) is the concatenation of n and 6n", "digit rule",
  "10 different gaps occur, from 106 to 15,000,400,006; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a009440 },
{ "a009441", "A009441", "a(n) is the concatenation of n and 7n", "digit rule",
  "10 different gaps occur, from 107 to 12,857,500,007; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a009441 },
{ "a009470", "A009470", "a(n) is the concatenation of n and 8n", "digit rule",
  "10 different gaps occur, from 108 to 11,250,100,008; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a009470 },
{ "a009474", "A009474", "a(n) is the concatenation of n and 9n", "digit rule",
  "10 different gaps occur, from 109 to 10,000,900,009; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a009474 },
{ "a015976", "A015976", "One iteration of Reverse and Add is needed to reach a palindrome", "digit rule",
  "388 different gaps occur, from 1 to 992; the level share is 10.71 %; L = 1 holds 76 % of the level class.",
  100000, 1, gen_a015976 },
{ "a015977", "A015977", "Two iterations of Reverse and Add are needed to reach a palindrome", "digit rule",
  "120 different gaps occur, from 1 to 540; the level share is 13.01 %; L = 1 holds 65 % of the level class.",
  100000, 0, gen_a015977 },
{ "a015979", "A015979", "Three iterations of Reverse and Add are needed to reach a palindrome", "digit rule",
  "153 different gaps occur, from 1 to 405; the level share is 16.20 %; L = 1 holds 50 % of the level class.",
  100000, 1, gen_a015979 },
{ "a015980", "A015980", "Four iterations of Reverse and Add are needed to reach a palindrome", "digit rule",
  "249 different gaps occur, from 1 to 790; the level share is 20.37 %; L = 1 holds 40 % of the level class.",
  100000, 0, gen_a015980 },
{ "a015982", "A015982", "Five iterations of Reverse and Add are needed to reach a palindrome", "digit rule",
  "217 different gaps occur, from 1 to 839; the level share is 21.00 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a015982 },
{ "a015984", "A015984", "Six iterations of Reverse and Add are needed to reach a palindrome", "digit rule",
  "272 different gaps occur, from 1 to 1,170; the level share is 23.79 %; L = 1 holds 33 % of the level class.",
  100000, 0, gen_a015984 },
{ "a000966", "A000966", "n! never ends in this many 0's", "digit rule",
  "The gaps are 1 and 6; the level share is 19.42 %; L = 1 holds 43 % of the level class.",
  100000, 1, gen_trailz_gaps },
{ "a001196", "A001196", "Double-bitters: only even length runs in binary expansion", "binary rule",
  "17 different gaps occur, from 3 to 8,589,934,593; the level share is 14.01 %; L = 3 holds 37 % of the level class; there are no ties; 17 terms do not decompose.",
  100000, 0, gen_dblbit },
{ "a003607", "A003607", "Location of 0's when natural numbers are listed in binary", "binary rule",
  "16 different gaps occur, from 1 to 16; the level share is 12.94 %; L = 1 holds 70 % of the level class.",
  100000, 0, gen_bin0pos },
{ "a003219", "A003219", "Self numbers divisible by sum of their digits (or, self numbers which are also Harshad numbers)", "digit rule",
  "216 different gaps occur, from 2 to 810; the level share is 33.17 %; 6 terms do not decompose.",
  100000, 1, gen_selfharshad },
{ "a001101", "A001101", "Moran numbers: k such that k/(sum of digits of k) is prime", "digit rule",
  "814 different gaps occur, from 1 to 1,421; the level share is 35.60 %.",
  100000, 1, gen_a001101 },
{ "a002796", "A002796", "Numbers that are divisible by each nonzero digit", "digit rule",
  "207 different gaps occur, from 1 to 1,512; the level share is 13.86 %.",
  100000, 1, gen_a002796 },
{ "a085370", "A085370", "Niven (or Harshad) numbers that are not divisible by 3", "digit rule",
  "278 different gaps occur, from 1 to 345; the level share is 23.62 %.",
  100000, 1, gen_a085370 },
{ "a085371", "A085371", "Non-Niven (or non-Harshad) numbers that are divisible by 3", "digit rule",
  "7 different gaps occur, from 3 to 24; the level share is 8.27 %; L = 3 holds 97 % of the level class.",
  100000, 1, gen_a085371 },
{ "a070938", "A070938", "Harshad numbers which terminate in their digital sum", "digit rule",
  "48 different gaps occur, from 1 to 909; the level share is 25.26 %.",
  100000, 1, gen_a070938 },
{ "a038367", "A038367", "Numbers n with property that (product of digits of n) is divisible by (sum of digits of n)", "digit rule",
  "10 different gaps occur, from 1 to 10; the level share is 11.84 %; L = 1 holds 64 % of the level class.",
  100000, 1, gen_a038367 },
{ "a038368", "A038368", "n is divisible by |(product of digits) - (sum of digits)|", "digit rule",
  "130 different gaps occur, from 1 to 180; the level share is 19.45 %.",
  100000, 1, gen_a038368 },
{ "a028838", "A028838", "Numbers whose sum of digits is a power of 2", "digit rule",
  "52 different gaps occur, from 1 to 500; the level share is 27.00 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_a028838 },
{ "a059094", "A059094", "Numbers whose sum of digits is a cube", "digit rule",
  "53 different gaps occur, from 1 to 1,000; the level share is 11.08 %; L = 9 holds 75 % of the level class.",
  100000, 1, gen_a059094 },
{ "a085802", "A085802", "Numbers whose sum of digits is a semiprime", "digit rule",
  "11 different gaps occur, from 1 to 12; the level share is 13.66 %; L = 1 holds 67 % of the level class.",
  100000, 1, gen_a085802 },
{ "a062713", "A062713", "Numbers k such that the sum of the digits of k is a prime factor of k", "digit rule",
  "635 different gaps occur, from 1 to 1,030; the level share is 33.54 %.",
  100000, 1, gen_a062713 },
{ "a062996", "A062996", "Numbers whose sum of digits is greater than or equal to its product of digits", "digit rule",
  "8 different gaps occur, from 1 to 10; the level share is 8.87 %; L = 1 holds 82 % of the level class.",
  100000, 1, gen_a062996 },
{ "a062997", "A062997", "Numbers whose sum of digits is strictly greater than its product of digits", "digit rule",
  "7 different gaps occur, from 1 to 10; the level share is 8.86 %; L = 1 holds 82 % of the level class.",
  100000, 1, gen_a062997 },
{ "a062998", "A062998", "Numbers whose sum of digits is less than or equal to its product of digits", "digit rule",
  "26 different gaps occur, from 1 to 11,127; the level share is 10.52 %; L = 1 holds 97 % of the level class.",
  100000, 1, gen_a062998 },
{ "a006364", "A006364", "Numbers k with an even number of 1's in binary, ignoring last bit", "binary rule",
  "The gaps are 1, 3 and 5; the level share is 13.11 %; L = 1 holds 69 % of the level class.",
  100000, 1, gen_a006364 },
{ "a000028", "A000028", "Let k = p_1^e_1 p_2^e_2 p_3^e_3 ... be the prime factorization of n. Sequence gives k such that the sum of the numbers of 1's in the binary expansions of e_1, e_2, e_3, ... is odd", "multiplicative",
  "16 different gaps occur, from 1 to 17; the level share is 12.74 %; L = 1 holds 71 % of the level class.",
  100000, 1, gen_a000028 },
{ "a000415", "A000415", "Numbers that are the sum of 2 but no fewer nonzero squares", "multiplicative",
  "33 different gaps occur, from 1 to 34; the level share is 16.02 %; L = 1 holds 38 % of the level class.",
  100000, 1, gen_a000415 },
{ "a000430", "A000430", "Primes and squares of primes", "multiplicative",
  "54 different gaps occur, from 1 to 114; the level share is 22.98 %; L = 1 holds 32 % of the level class.",
  100000, 1, gen_a000430 },
{ "a002035", "A002035", "Numbers that contain primes to odd powers only", "multiplicative",
  "7 different gaps occur, from 1 to 7; the level share is 10.08 %; L = 1 holds 96 % of the level class.",
  100000, 1, gen_a002035 },
{ "a005238", "A005238", "Numbers k such that k, k+1 and k+2 have the same number of divisors", "multiplicative",
  "598 different gaps occur, from 1 to 992; the level share is 36.72 %.",
  100000, 1, gen_a005238 },
{ "a008846", "A008846", "Hypotenuses of primitive Pythagorean triangles", "multiplicative",
  "18 different gaps occur, from 4 to 72; the level share is 23.76 %; L = 3 holds 42 % of the level class.",
  100000, 1, gen_a008846 },
{ "a014567", "A014567", "Numbers k such that k and sigma(k) are relatively prime, where sigma(k) = sum of divisors of k (A000203)", "multiplicative",
  "12 different gaps occur, from 1 to 14; the level share is 13.55 %; L = 1 holds 58 % of the level class.",
  100000, 1, gen_a014567 },
{ "a030230", "A030230", "Numbers that have an odd number of distinct prime divisors", "multiplicative",
  "17 different gaps occur, from 1 to 20; the level share is 12.70 %; L = 1 holds 71 % of the level class.",
  100000, 1, gen_a030230 },
{ "a030231", "A030231", "Numbers with an even number of distinct prime factors", "multiplicative",
  "18 different gaps occur, from 1 to 19; the level share is 12.71 %; L = 1 holds 71 % of the level class.",
  100000, 1, gen_a030231 },
{ "a033950", "A033950", "Refactorable numbers: number of divisors of k divides k. Also known as tau numbers", "multiplicative",
  "82 different gaps occur, from 1 to 180; the level share is 16.87 %; there are no ties.",
  100000, 1, gen_a033950 },
{ "a037020", "A037020", "Numbers whose sum of proper (or aliquot) divisors is a prime", "multiplicative",
  "91 different gaps occur, from 1 to 152; the level share is 28.35 %; L = 1 holds 49 % of the level class.",
  100000, 1, gen_a037020 },
{ "a038509", "A038509", "Composite numbers congruent to +-1 mod 6", "multiplicative",
  "10 different gaps occur, from 2 to 20; the level share is 15.14 %; L = 3 holds 55 % of the level class.",
  100000, 1, gen_a038509 },
{ "a039956", "A039956", "Even squarefree numbers", "multiplicative",
  "6 different gaps occur, from 4 to 24; the level share is 17.47 %; L = 2 holds 78 % of the level class; there are no ties.",
  100000, 1, gen_a039956 },
{ "a039957", "A039957", "Squarefree numbers congruent to 3 mod 4", "multiplicative",
  "6 different gaps occur, from 4 to 24; the level share is 23.03 %; L = 1 holds 55 % of the level class; there are no ties.",
  100000, 1, gen_a039957 },
{ "a025583", "A025583", "Composite numbers that are not the sum of 2 primes", "multiplicative",
  "9 different gaps occur, from 2 to 22; the level share is 5.40 %; L = 1 holds 57 % of the level class.",
  100000, 1, gen_a025583 },
{ "a031363", "A031363", "Positive numbers of the form x^2 + xy - y^2; or, of the form 5x^2 - y^2", "multiplicative",
  "46 different gaps occur, from 1 to 46; the level share is 23.40 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a031363 },
{ "a007921", "A007921", "Numbers that are not the difference of two primes", "primes",
  "The gaps are 2, 4 and 6; the level share is 18.57 %; L = 1 holds 93 % of the level class.",
  100000, 1, gen_a007921 },
{ "a001912", "A001912", "Numbers k such that 4*k^2 + 1 is prime", "prime values",
  "87 different gaps occur, from 1 to 106; the level share is 23.92 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_a001912 },
{ "a002837", "A002837", "Numbers k such that k^2 - k + 41 is prime", "prime values",
  "35 different gaps occur, from 1 to 39; the level share is 16.07 %; L = 1 holds 56 % of the level class.",
  100000, 1, gen_a002837 },
{ "a002731", "A002731", "Numbers k such that (k^2 + 1)/2 is prime", "prime values",
  "81 different gaps occur, from 2 to 210; the level share is 33.03 %; L = 1 holds 52 % of the level class.",
  100000, 1, gen_a002731 },
{ "a002815", "A002815", "a(n) = n + Sum_{k=1..n} pi(k), where pi() = A000720", "summatory",
  "9,593 different gaps occur, from 1 to 9,593; the level share is 70.24 %.",
  100000, 0, gen_a002815 },
{ "a002821", "A002821", "a(n) = nearest integer to n^(3/2)", "polynomial",
  "475 different gaps occur, from 1 to 475; the level share is 46.95 %.",
  100000, 0, gen_a002821 },
{ "a002984", "A002984", "a(0) = 1; for n > 0, a(n) = a(n-1) + floor(sqrt(a(n-1)))", "self-referential",
  "49,992 different gaps occur, from 1 to 49,992; the level share is 100.00 %; 50.0 % of terms are forced level (l <= d^2); there are no ties.",
  100000, 0, gen_a002984 },
{ "a003072", "A003072", "Numbers that are the sum of 3 positive cubes", "quadratic form",
  "85 different gaps occur, from 1 to 100; the level share is 21.97 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a003072 },
{ "a025395", "A025395", "Numbers that are the sum of 3 positive cubes in exactly 1 way", "quadratic form",
  "98 different gaps occur, from 1 to 124; the level share is 22.89 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_a025395 },
{ "a008917", "A008917", "Numbers that are the sum of 3 positive cubes in more than one way", "quadratic form",
  "509 different gaps occur, from 1 to 812; the level share is 34.53 %.",
  100000, 1, gen_a008917 },
{ "a003623", "A003623", "Wythoff AB-numbers: floor(floor(n*phi^2)*phi), where phi = (1+sqrt(5))/2", "Beatty",
  "The gaps are 3 and 5; the level share is 18.56 %; L = 1 holds 45 % of the level class.",
  100000, 1, gen_a003623 },
{ "a004201", "A004201", "Accept one, reject one, accept two, reject two, ", "block",
  "447 different gaps occur, from 1 to 447; the level share is 9.38 %; L = 1 holds 97 % of the level class; there are no ties.",
  100000, 1, gen_a004201 },
{ "a001974", "A001974", "Numbers that are the sum of 3 distinct squares, i.e., numbers of the form x^2 + y^2 + z^2 with 0 <= x < y < z", "quadratic form",
  "5 different gaps occur, from 1 to 5; the level share is 8.57 %; L = 1 holds 99 % of the level class.",
  100000, 1, gen_a001974 },
{ "a004432", "A004432", "Numbers that are the sum of 3 distinct nonzero squares", "quadratic form",
  "6 different gaps occur, from 1 to 7; the level share is 8.59 %; L = 1 holds 99 % of the level class.",
  100000, 1, gen_a004432 },
{ "a004433", "A004433", "Numbers that are the sum of 4 distinct nonzero squares: of form w^2+x^2+y^2+z^2 with 0<w<x<y<z", "quadratic form",
  "7 different gaps occur, from 1 to 9; the level share is 9.60 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a004433 },
{ "a001983", "A001983", "Numbers that are the sum of 2 distinct squares: of form x^2 + y^2 with 0 <= x < y", "quadratic form",
  "32 different gaps occur, from 1 to 34; the level share is 15.95 %; L = 1 holds 39 % of the level class.",
  100000, 1, gen_a001983 },
{ "a004999", "A004999", "Sums of two nonnegative cubes", "quadratic form",
  "5,360 different gaps occur, from 1 to 16,955; the level share is 49.93 %; 1.3 % of terms are forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a004999 },
{ "a045980", "A045980", "Numbers of the form x^3 + y^3 or x^3 - y^3", "quadratic form",
  "1,559 different gaps occur, from 1 to 3,375; the level share is 40.24 %.",
  100000, 1, gen_a045980 },
{ "a005658", "A005658", "If n appears so do 2n, 3n+2, 6n+3", "self-referential",
  "50 different gaps occur, from 1 to 88; the level share is 13.51 %; L = 1 holds 48 % of the level class.",
  100000, 1, gen_a005658 },
{ "a006489", "A006489", "Numbers k such that k-6, k, and k+6 are primes", "primes",
  "1,781 different gaps occur, from 2 to 11,774; the level share is 60.47 %; 1.5 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a006489 },
{ "a007064", "A007064", "Numbers not of form \"nearest integer to n*tau\", tau = (1+sqrt(5))/2", "Beatty",
  "The gaps are 2 and 3; the level share is 15.43 %; L = 1 holds 57 % of the level class.",
  100000, 1, gen_a007064 },
{ "a007066", "A007066", "a(n) = 1 + ceiling((n-1)*phi^2), phi = (1+sqrt(5))/2", "Beatty",
  "The gaps are 2 and 3; the level share is 15.48 %; L = 1 holds 57 % of the level class.",
  100000, 1, gen_a007066 },
{ "a007491", "A007491", "Smallest prime > n^2", "primes",
  "61,747 different gaps occur, from 3 to 200,016; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a007491 },
{ "a013916", "A013916", "Numbers k such that the sum of the first k primes is prime", "primes",
  "141 different gaps occur, from 1 to 366; the level share is 24.14 %; L = 2 holds 32 % of the level class; there are no ties; 6 terms do not decompose.",
  100000, 1, gen_a013916 },
{ "a014657", "A014657", "Numbers m that divide 2^k + 1 for some nonnegative k", "multiplicative",
  "34 different gaps occur, from 1 to 70; the level share is 23.59 %; L = 1 holds 51 % of the level class.",
  100000, 1, gen_a014657 },
{ "a014661", "A014661", "Numbers that do not divide 2^k + 1 for any k>0", "multiplicative",
  "The gaps are 1 and 2; the level share is 11.18 %; L = 1 holds 77 % of the level class.",
  100000, 1, gen_a014661 },
{ "a022797", "A022797", "a(n) = n-th prime + n-th nonprime", "primes",
  "94 different gaps occur, from 3 to 115; the level share is 26.12 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a022797 },
{ "a025475", "A025475", "1 and the prime powers p^m where m >= 2, thus excluding the primes", "powers",
  "95,090 different gaps occur, from 1 to 238,596,672; the level share is 99.97 %; 99.9 % of terms are forced level (l <= d^2); there are no ties.",
  100000, 1, gen_a025475 },
{ "a026430", "A026430", "a(n) is the sum of first n terms of A001285 (Thue-Morse sequence)", "binary rule",
  "The gaps are 1 and 2; the level share is 17.52 %; L = 1 holds 79 % of the level class.",
  100000, 0, gen_a026430 },
{ "a027862", "A027862", "Primes of the form j^2 + (j+1)^2", "primes",
  "98,990 different gaps occur, from 8 to 269,858,020; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a027862 },
{ "a034707", "A034707", "Numbers that are sums (of a nonempty sequence) of consecutive primes", "primes",
  "16 different gaps occur, from 1 to 16; the level share is 13.14 %; L = 1 holds 69 % of the level class.",
  100000, 1, gen_a034707 },
{ "a035106", "A035106", "1, together with numbers of the form k*(k+1) or k*(k+2), k > 0", "polynomial",
  "50,001 different gaps occur, from 1 to 50,001; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a035106 },
{ "a038550", "A038550", "Products of an odd prime and a power of two (sorted)", "primes",
  "46 different gaps occur, from 1 to 49; the level share is 11.98 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a038550 },
{ "a045636", "A045636", "Numbers of the form p^2 + q^2, with p and q primes", "primes",
  "90 different gaps occur, from 3 to 1,680; the level share is 44.45 %; L = 2 holds 44 % of the level class; there are no ties.",
  100000, 1, gen_a045636 },
{ "a045699", "A045699", "Numbers of the form p^2 + q^3, p,q prime", "primes",
  "3,951 different gaps occur, from 1 to 16,800; the level share is 48.64 %; 1.5 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a045699 },
{ "a001463", "A001463", "Partial sums of A001462; also a(n) is the last occurrence of n in A001462", "self-referential",
  "1,478 different gaps occur, from 2 to 1,479; the level share is 55.89 %.",
  100000, 1, gen_a001463 },
{ "a001768", "A001768", "Sorting numbers: number of comparisons for merge insertion sort of n elements", "summatory",
  "17 different gaps occur, from 1 to 17; the level share is 28.28 %.",
  100000, 1, gen_a001768 },
{ "a005598", "A005598", "a(n) = 1 + Sum_{i=1..n} (n-i+1)*phi(i)", "divisor functions",
  "Every gap is different, from 1 to 3,039,650,754; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a005598 },
{ "a003278", "A003278", "Szekeres's sequence: a(n)-1 in ternary = n-1 in binary; also: a(1) = 1, a(2) = 2, and thereafter a(n) is smallest number k which avoids any 3-term arithmetic progression in a(1), a(2), ..., a(n-1), k", "digit rule",
  "17 different gaps occur, from 1 to 21,523,361; the level share is 0.36 %; 17 terms do not decompose.",
  100000, 1, gen_a003278 },
{ "a006285", "A006285", "Odd numbers not of form p + 2^k (de Polignac numbers)", "primes",
  "115 different gaps occur, from 2 to 264; the level share is 29.68 %; L = 1 holds 39 % of the level class.",
  100000, 1, gen_a006285 },
{ "a040098", "A040098", "Primes p such that x^4 = 2 has a solution mod p", "primes",
  "110 different gaps occur, from 2 to 366; the level share is 30.55 %.",
  100000, 1, gen_a040098 },
{ "a045315", "A045315", "Primes p such that x^8 = 2 has a solution mod p", "primes",
  "125 different gaps occur, from 2 to 504; the level share is 31.67 %.",
  100000, 1, gen_a045315 },
{ "a002407", "A002407", "Cuban primes: primes which are the difference of two consecutive cubes", "primes",
  "98,612 different gaps occur, from 12 to 365,910,360; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a002407 },
{ "a014574", "A014574", "Average of twin prime pairs", "primes",
  "247 different gaps occur, from 2 to 2,190; the level share is 31.79 %; there are no ties.",
  100000, 1, gen_twinavg },
{ "a002081", "A002081", "Numbers congruent to {2, 4, 8, 16} (mod 20)", "residue class",
  "The gaps are 2, 4, 6 and 8; the level share is 0.01 %.",
  100000, 0, gen_a002081 },
{ "a016825", "A016825", "Positive integers congruent to 2 (mod 4): a(n) = 4*n+2, for n >= 0", "residue class",
  "The gap is always 4; the level share is 17.98 %; L = 2 holds 100 % of the level class; there are no ties.",
  100000, 0, gen_a016825 },
{ "a042965", "A042965", "Nonnegative integers not congruent to 2 mod 4", "residue class",
  "The gaps are 1 and 2; the level share is 12.50 %; L = 1 holds 100 % of the level class; there are no ties.",
  100000, 1, gen_a042965 },
{ "a003626", "A003626", "Inert rational primes in Q(sqrt(-5))", "primes",
  "97 different gaps occur, from 2 to 300; the level share is 27.82 %.",
  100000, 1, gen_a003626 },
{ "a049445", "A049445", "Numbers k with the property that the number of 1's in binary expansion of k (see A000120) divides k", "binary rule",
  "75 different gaps occur, from 1 to 103; the level share is 19.01 %.",
  100000, 1, gen_a049445 },
{ "a064481", "A064481", "Numbers which are divisible by the sum of their base-5 digits", "digit rule",
  "61 different gaps occur, from 1 to 63; the level share is 18.59 %.",
  100000, 1, gen_a064481 },
{ "a101813", "A101813", "Odd Niven (or Harshad) numbers: odd numbers that are divisible by the sum of their digits", "digit rule",
  "185 different gaps occur, from 2 to 480; the level share is 33.93 %; L = 1 holds 32 % of the level class.",
  100000, 1, gen_a101813 },
{ "a101814", "A101814", "Even Niven (or Harshad) numbers: even numbers that are divisible by the sum of their digits", "digit rule",
  "45 different gaps occur, from 2 to 180; the level share is 18.82 %.",
  100000, 1, gen_a101814 },
{ "a118363", "A118363", "Factorial base Niven (or Harshad) numbers: numbers that are divisible by the sum of their factorial base digits", "digit rule",
  "173 different gaps occur, from 1 to 210; the level share is 22.45 %.",
  100000, 1, gen_a118363 },
{ "a077436", "A077436", "Let B(n) be the sum of binary digits of n. This sequence contains n such that B(n) = B(n^2)", "binary rule",
  "3,140 different gaps occur, from 1 to 57,366; the level share is 28.47 %; 2.1 % of terms are forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a077436 },
{ "a038366", "A038366", "n is divisible by (product of digits) + (sum of digits)", "digit rule",
  "131 different gaps occur, from 1 to 180; the level share is 19.40 %.",
  100000, 1, gen_a038366 },
{ "a016052", "A016052", "a(1) = 3; for n >= 1, a(n+1) = a(n) + sum of its digits", "digit rule",
  "12 different gaps occur, from 3 to 51; the level share is 14.68 %; L = 9 holds 60 % of the level class.",
  100000, 1, gen_a016052 },
{ "a016096", "A016096", "a(n+1) = a(n) + sum of its digits, with a(1) = 9", "digit rule",
  "6 different gaps occur, from 9 to 54; the level share is 17.41 %; L = 9 holds 55 % of the level class.",
  100000, 1, gen_a016096 },
{ "a033298", "A033298", "a(n+1) = a(n) + sum of digits of a(n)^2, with a(1) = 1", "digit rule",
  "13 different gaps occur, from 1 to 99; the level share is 32.09 %; L = 3 holds 32 % of the level class; there are no ties.",
  100000, 1, gen_a033298 },
{ "a007612", "A007612", "a(n+1) = a(n) + digital root (A010888) of a(n)", "digit rule",
  "6 different gaps occur, from 1 to 8; the level share is 0.01 %; L = 3 holds 38 % of the level class.",
  100000, 1, gen_a007612 },
{ "a019551", "A019551", "a(n) is the concatenation of n and 3n", "digit rule",
  "11 different gaps occur, from 13 to 30,000,700,003; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a019551 },
{ "a019552", "A019552", "a(n) is the concatenation of n and 4n", "digit rule",
  "11 different gaps occur, from 14 to 22,500,100,004; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a019552 },
{ "a028839", "A028839", "Sum of digits of n is a square", "digit rule",
  "27 different gaps occur, from 1 to 49; the level share is 25.18 %; L = 1 holds 37 % of the level class.",
  100000, 1, gen_a028839 },
{ "a002327", "A002327", "Primes of the form k^2 - k - 1", "primes",
  "98,161 different gaps occur, from 6 to 79,615,620; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a002327 },
{ "a002328", "A002328", "Numbers k such that k^2 - k - 1 is prime", "prime values",
  "62 different gaps occur, from 1 to 64; the level share is 20.35 %; L = 1 holds 40 % of the level class.",
  100000, 1, gen_a002328 },
{ "a002383", "A002383", "Primes of form k^2 + k + 1", "primes",
  "97,970 different gaps occur, from 4 to 251,450,628; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a002383 },
{ "a002384", "A002384", "Numbers m such that m^2 + m + 1 is prime", "prime values",
  "101 different gaps occur, from 1 to 132; the level share is 24.72 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a002384 },
{ "a002970", "A002970", "Numbers k such that 4*k^2 + 9 is prime", "prime values",
  "130 different gaps occur, from 1 to 171; the level share is 23.33 %.",
  100000, 1, gen_a002970 },
{ "a002971", "A002971", "Numbers k such that 4*k^2 + 25 is prime", "prime values",
  "64 different gaps occur, from 1 to 90; the level share is 18.76 %; L = 1 holds 40 % of the level class.",
  100000, 1, gen_a002971 },
{ "a005846", "A005846", "Primes of the form k^2 + k + 41", "primes",
  "96,902 different gaps occur, from 2 to 22,664,950; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a005846 },
{ "a007635", "A007635", "Primes of form n^2 + n + 17", "primes",
  "98,426 different gaps occur, from 2 to 61,808,838; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a007635 },
{ "a007637", "A007637", "Primes of form 3*k^2 - 3*k + 23", "primes",
  "97,802 different gaps occur, from 6 to 92,677,830; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a007637 },
{ "a007639", "A007639", "Primes of form 2n^2 - 2n + 19", "primes",
  "98,203 different gaps occur, from 4 to 94,621,500; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a007639 },
{ "a007641", "A007641", "Primes of the form 2*k^2 + 29", "primes",
  "98,737 different gaps occur, from 2 to 53,963,848; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a007641 },
{ "a023219", "A023219", "Primes p such that 5p+6 is a prime", "primes",
  "286 different gaps occur, from 2 to 836; the level share is 34.52 %.",
  100000, 1, gen_a023219 },
{ "a027753", "A027753", "Primes of form n^2 + n + 3", "primes",
  "99,098 different gaps occur, from 2 to 1,357,903,764; every decomposable term is forced level (l <= d^2); 8 terms do not decompose.",
  100000, 1, gen_a027753 },
{ "a027755", "A027755", "Primes of the form k^2 + k + 5", "primes",
  "99,023 different gaps occur, from 2 to 341,654,130; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a027755 },
{ "a027758", "A027758", "Primes of the form k^2 + k + 9", "primes",
  "98,728 different gaps occur, from 18 to 1,652,570,208; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a027758 },
{ "a027861", "A027861", "Numbers k such that k^2 + (k+1)^2 is prime", "prime values",
  "81 different gaps occur, from 1 to 105; the level share is 22.33 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a027861 },
{ "a027863", "A027863", "Numbers k such that k^2 + (k+1)^2 + (k+2)^2 is prime", "prime values",
  "113 different gaps occur, from 2 to 290; the level share is 23.47 %; L = 2 holds 31 % of the level class.",
  100000, 1, gen_a027863 },
{ "a027866", "A027866", "Numbers k such that k^2 + (k+1)^2 + (k+2)^2 + (k+3)^2 + (k+4)^2 + (k+5)^2 is prime", "prime values",
  "102 different gaps occur, from 1 to 121; the level share is 21.37 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a027866 },
{ "a027867", "A027867", "Primes of the form n^2 + (n+1)^2 + (n+2)^2 + (n+3)^2 + (n+4)^2 + (n+5)^2", "primes",
  "99,063 different gaps occur, from 60 to 1,237,739,760; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a027867 },
{ "a037029", "A037029", "Primes of the form 666*n + 1", "primes",
  "52 different gaps occur, from 666 to 35,964; the level share is 69.07 %; 7.0 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a037029 },
{ "a048059", "A048059", "Primes of the form k^2 + k + 11", "primes",
  "98,862 different gaps occur, from 2 to 112,335,618; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a048059 },
{ "a048988", "A048988", "Primes of the form 4*k^2 + 4*k + 59", "primes",
  "97,604 different gaps occur, from 8 to 116,656,960; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a048988 },
{ "a050265", "A050265", "Primes of the form 2*n^2 + 11", "primes",
  "99,333 different gaps occur, from 2 to 166,743,070; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a050265 },
{ "a051647", "A051647", "Primes p such that 210*p + 1 is also prime", "primes",
  "289 different gaps occur, from 1 to 964; the level share is 34.31 %.",
  100000, 1, gen_a051647 },
{ "a051653", "A051653", "Primes p such that 2310*p + 1 is also prime", "primes",
  "280 different gaps occur, from 2 to 728; the level share is 34.77 %.",
  100000, 1, gen_a051653 },
{ "a051654", "A051654", "Primes p such that 30030*p + 1 is also prime", "primes",
  "298 different gaps occur, from 2 to 1,052; the level share is 34.73 %.",
  100000, 1, gen_a051654 },
{ "a052291", "A052291", "Primes p such that 4p^2 + 1 is also prime", "primes",
  "540 different gaps occur, from 1 to 3,630; the level share is 45.97 %; 6 terms do not decompose.",
  100000, 1, gen_a052291 },
{ "a053182", "A053182", "Primes p such that p^2 + p + 1 is prime", "primes",
  "469 different gaps occur, from 1 to 4,392; the level share is 50.42 %; L = 1 holds 31 % of the level class; there are no ties; 6 terms do not decompose.",
  100000, 1, gen_a053182 },
{ "a053184", "A053184", "Primes p such that p^2+p-1 is prime", "primes",
  "562 different gaps occur, from 1 to 2,072; the level share is 39.52 %.",
  100000, 1, gen_a053184 },
{ "a055494", "A055494", "Numbers k such that k^2 - k + 1 is prime", "prime values",
  "101 different gaps occur, from 1 to 132; the level share is 24.96 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a055494 },
{ "a056906", "A056906", "Numbers k such that 36*k^2 + 5 is prime", "prime values",
  "83 different gaps occur, from 1 to 104; the level share is 21.06 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a056906 },
{ "a056908", "A056908", "Numbers k such that 36*k^2 + 36*k + 13 is prime", "prime values",
  "91 different gaps occur, from 1 to 134; the level share is 22.02 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a056908 },
{ "a057604", "A057604", "Primes of the form 4*k^2 + 163", "primes",
  "98,729 different gaps occur, from 4 to 93,888,512; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a057604 },
{ "a059425", "A059425", "Primes of form n^2 + 19n + 17", "primes",
  "98,132 different gaps occur, from 20 to 41,112,432; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a059425 },
{ "a060844", "A060844", "Primes of the form 6*k^2 + 6*k + 31", "primes",
  "96,829 different gaps occur, from 12 to 150,530,940; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a060844 },
{ "a061242", "A061242", "Primes of the form 9*k - 1", "primes",
  "41 different gaps occur, from 18 to 756; the level share is 42.75 %; L = 1 holds 38 % of the level class; there are no ties.",
  100000, 1, gen_a061242 },
{ "a062800", "A062800", "Primes of form 100*k + 1", "primes",
  "58 different gaps occur, from 100 to 6,400; the level share is 53.08 %; 1.0 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a062800 },
{ "a063472", "A063472", "Primes of the form 666*k - 1", "primes",
  "56 different gaps occur, from 666 to 37,962; the level share is 68.86 %; 7.0 % of terms are forced level (l <= d^2); there are no ties.",
  100000, 1, gen_a063472 },
{ "a065508", "A065508", "Primes p such that p^2 - p + 1 is prime", "primes",
  "457 different gaps occur, from 1 to 4,380; the level share is 50.27 %; L = 1 holds 31 % of the level class; 6 terms do not decompose.",
  100000, 1, gen_a065508 },
{ "a066049", "A066049", "Numbers k such that 2*k^2 - 1 is a prime", "prime values",
  "61 different gaps occur, from 1 to 61; the level share is 20.78 %; L = 1 holds 41 % of the level class.",
  100000, 1, gen_a066049 },
{ "a066436", "A066436", "Primes of the form 2*n^2 - 1", "primes",
  "98,979 different gaps occur, from 10 to 142,133,640; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a066436 },
{ "a073102", "A073102", "Primes of the form 210n + 1", "primes",
  "34 different gaps occur, from 210 to 7,350; the level share is 62.77 %; 1.2 % of terms are forced level (l <= d^2); L = 1 holds 39 % of the level class.",
  100000, 1, gen_a073102 },
{ "a076339", "A076339", "Primes of the form 512*k+1", "primes",
  "80 different gaps occur, from 512 to 49,152; the level share is 64.20 %; 8.5 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a076339 },
{ "a076727", "A076727", "Primes of the form x^2 + (x+3)^2", "primes",
  "98,702 different gaps occur, from 12 to 975,329,472; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a076727 },
{ "a083022", "A083022", "Numbers n such that 4*n^2 - 3 is prime", "prime values",
  "85 different gaps occur, from 1 to 102; the level share is 18.87 %.",
  100000, 1, gen_a083022 },
{ "a086285", "A086285", "Numbers k such that 1 + 2k + 3k^2 is prime", "prime values",
  "85 different gaps occur, from 2 to 194; the level share is 24.05 %; L = 2 holds 39 % of the level class.",
  100000, 1, gen_a086285 },
{ "a086298", "A086298", "Numbers n such that 1-2n+3n^2 is prime", "prime values",
  "87 different gaps occur, from 2 to 200; the level share is 23.68 %; L = 2 holds 38 % of the level class.",
  100000, 1, gen_a086298 },
{ "a086303", "A086303", "Numbers n such that n+15 is prime", "prime values",
  "53 different gaps occur, from 2 to 114; the level share is 14.92 %; L = 2 holds 30 % of the level class; there are no ties.",
  100000, 1, gen_a086303 },
{ "a086304", "A086304", "Numbers n such that n+6 is prime", "prime values",
  "53 different gaps occur, from 2 to 114; the level share is 25.65 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_a086304 },
{ "a088572", "A088572", "Numbers n such that (2n+1)^2 - 2 is prime", "prime values",
  "64 different gaps occur, from 1 to 69; the level share is 20.64 %; L = 1 holds 42 % of the level class.",
  100000, 1, gen_a088572 },
{ "a088758", "A088758", "Numbers k such that (4*k + 1)^2 + (4*k + 2)^2 is prime", "prime values",
  "91 different gaps occur, from 1 to 100; the level share is 22.13 %; L = 1 holds 32 % of the level class.",
  100000, 1, gen_a088758 },
{ "a088759", "A088759", "Numbers k such that (4*k+3)^2 + (4*k+2)^2 is prime", "prime values",
  "97 different gaps occur, from 1 to 137; the level share is 23.84 %; L = 1 holds 38 % of the level class.",
  100000, 1, gen_a088759 },
{ "a088955", "A088955", "Primes of the form 60*k + 1", "primes",
  "36 different gaps occur, from 60 to 2,160; the level share is 52.60 %; L = 1 holds 40 % of the level class.",
  100000, 1, gen_a088955 },
{ "a088967", "A088967", "Numbers n such that n+9 is a prime", "prime values",
  "53 different gaps occur, from 2 to 114; the level share is 17.18 %; L = 2 holds 34 % of the level class; there are no ties.",
  100000, 1, gen_a088967 },
{ "a089001", "A089001", "Numbers k such that 2*k^2 + 1 is prime", "prime values",
  "59 different gaps occur, from 2 to 192; the level share is 19.47 %; L = 3 holds 42 % of the level class.",
  100000, 1, gen_a089001 },
{ "a089008", "A089008", "Numbers k such that 18*k^2 + 1 is prime", "prime values",
  "58 different gaps occur, from 1 to 64; the level share is 19.46 %; L = 1 holds 42 % of the level class.",
  100000, 1, gen_a089008 },
{ "a089063", "A089063", "Numbers k such that 840*k + 175177943 is a prime", "prime values",
  "39 different gaps occur, from 1 to 40; the level share is 17.17 %; L = 1 holds 49 % of the level class.",
  100000, 1, gen_a089063 },
{ "a089373", "A089373", "Numbers k such that k^2 - 7*k + 7 is prime", "prime values",
  "120 different gaps occur, from 1 to 168; the level share is 25.41 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_a089373 },
{ "a089376", "A089376", "Primes of the form k^2 - 7*k + 7", "primes",
  "98,615 different gaps occur, from 24 to 464,991,240; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a089376 },
{ "a089438", "A089438", "Primes p such that 6p+11 is also a prime", "primes",
  "337 different gaps occur, from 1 to 898; the level share is 36.85 %.",
  100000, 1, gen_a089438 },
{ "a089441", "A089441", "Primes p such that 16*p+17 is a prime", "primes",
  "275 different gaps occur, from 6 to 2,340; the level share is 48.31 %; L = 1 holds 35 % of the level class; there are no ties.",
  100000, 1, gen_a089441 },
{ "a089593", "A089593", "Numbers k such that k^2 + 2k + 2 is prime", "prime values",
  "88 different gaps occur, from 1 to 212; the level share is 31.44 %; L = 1 holds 48 % of the level class.",
  100000, 1, gen_a089593 },
{ "a089623", "A089623", "Numbers n such that n^2 + 2n - 1 is prime", "prime values",
  "65 different gaps occur, from 1 to 138; the level share is 20.64 %; L = 2 holds 42 % of the level class.",
  100000, 1, gen_a089623 },
{ "a089681", "A089681", "Numbers n such that 3n^2 - 1 is prime", "prime values",
  "58 different gaps occur, from 1 to 120; the level share is 19.74 %; L = 2 holds 43 % of the level class.",
  100000, 1, gen_a089681 },
{ "a089682", "A089682", "Primes of the form 3*m^2 - 1", "primes",
  "99,215 different gaps occur, from 9 to 875,458,032; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a089682 },
{ "a089747", "A089747", "Numbers n such that n^2 - 2n + 5 is prime", "prime values",
  "82 different gaps occur, from 2 to 178; the level share is 21.60 %; L = 2 holds 37 % of the level class; there are no ties.",
  100000, 1, gen_a089747 },
{ "a090187", "A090187", "Primes of the form 11*n+2", "primes",
  "57 different gaps occur, from 11 to 2,002; the level share is 41.80 %; there are no ties.",
  100000, 1, gen_a090187 },
{ "a090562", "A090562", "Primes of the form 5k^2 + 5k + 1", "primes",
  "98,335 different gaps occur, from 20 to 322,474,490; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a090562 },
{ "a090563", "A090563", "Numbers k such that 5*k^2 + 5*k + 1 is prime", "prime values",
  "54 different gaps occur, from 1 to 71; the level share is 19.11 %; L = 1 holds 45 % of the level class.",
  100000, 1, gen_a090563 },
{ "a090684", "A090684", "Primes of the form 8*k^2 - 1", "primes",
  "99,127 different gaps occur, from 24 to 741,263,360; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a090684 },
{ "a090685", "A090685", "Primes of the form 8*k^2 + 1", "primes",
  "99,280 different gaps occur, from 648 to 5,455,465,344; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a090685 },
{ "a090686", "A090686", "Primes of the form 6n^2 - 1", "primes",
  "99,422 different gaps occur, from 18 to 731,095,680; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a090686 },
{ "a090687", "A090687", "Primes of the form 6*k^2 + 1", "primes",
  "99,454 different gaps occur, from 54 to 1,548,693,090; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a090687 },
{ "a090696", "A090696", "Numbers k such that k^2 - 11 is a prime", "prime values",
  "100 different gaps occur, from 2 to 308; the level share is 25.24 %; L = 2 holds 35 % of the level class.",
  100000, 1, gen_a090696 },
{ "a090698", "A090698", "Primes of the form 2*n^2+1", "primes",
  "99,186 different gaps occur, from 16 to 1,438,225,920; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a090698 },
{ "a091271", "A091271", "Numbers k such that 4*k^2-11 is a prime", "prime values",
  "100 different gaps occur, from 1 to 154; the level share is 25.23 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_a091271 },
{ "a091567", "A091567", "Primes p such that p^2-p-1 is prime", "primes",
  "561 different gaps occur, from 2 to 1,542; the level share is 40.35 %.",
  100000, 1, gen_a091567 },
{ "a092968", "A092968", "Numbers n such that 2n^2 + 11 is a prime", "prime values",
  "68 different gaps occur, from 1 to 79; the level share is 19.62 %; L = 1 holds 41 % of the level class.",
  100000, 1, gen_a092968 },
{ "a093359", "A093359", "Primes of the form 28*k + 1", "primes",
  "54 different gaps occur, from 28 to 1,596; the level share is 43.55 %.",
  100000, 1, gen_a093359 },
{ "a093838", "A093838", "Primes of the form 36n + 1", "primes",
  "41 different gaps occur, from 36 to 1,548; the level share is 47.29 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_a093838 },
{ "a094210", "A094210", "Numbers k such that k^2 + 3k + 1 is a prime", "prime values",
  "62 different gaps occur, from 1 to 64; the level share is 20.64 %; L = 1 holds 42 % of the level class.",
  100000, 1, gen_a094210 },
{ "a094407", "A094407", "Primes of the form 16n+1", "primes",
  "59 different gaps occur, from 16 to 1,344; the level share is 38.24 %.",
  100000, 1, gen_a094407 },
{ "a095995", "A095995", "Primes of the form 100n - 1", "primes",
  "57 different gaps occur, from 100 to 6,600; the level share is 53.11 %; 1.0 % of terms are forced level (l <= d^2); there are no ties.",
  100000, 1, gen_a095995 },
{ "a096689", "A096689", "Numbers n such that 2n^2 + 3n + 3 is prime", "prime values",
  "95 different gaps occur, from 2 to 234; the level share is 19.36 %.",
  100000, 1, gen_a096689 },
{ "a096691", "A096691", "Numbers n such that 8n^2 + 6n + 3 is prime", "prime values",
  "95 different gaps occur, from 1 to 117; the level share is 19.35 %.",
  100000, 1, gen_a096691 },
{ "a098828", "A098828", "Primes of the form 2*n^2 + 2*n - 1", "primes",
  "97,983 different gaps occur, from 8 to 343,225,320; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a098828 },
{ "a099007", "A099007", "Primes of the form 6n^2 - 2n - 1", "primes",
  "99,955 different gaps occur, from 16 to 1,052,281,164; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a099007 },
{ "a100201", "A100201", "Primes of the form 23*k+3", "primes",
  "63 different gaps occur, from 46 to 3,220; the level share is 46.46 %.",
  100000, 1, gen_a100201 },
{ "a100202", "A100202", "Primes of the form 13*k + 3", "primes",
  "61 different gaps occur, from 26 to 1,638; the level share is 42.66 %.",
  100000, 1, gen_a100202 },
{ "a100203", "A100203", "Primes of the form 37n+3", "primes",
  "67 different gaps occur, from 74 to 5,402; the level share is 49.18 %.",
  100000, 1, gen_a100203 },
{ "a100494", "A100494", "Primes of the form 47*k + 3", "primes",
  "69 different gaps occur, from 94 to 7,144; the level share is 51.09 %; 1.2 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a100494 },
{ "a100760", "A100760", "Primes of the form 47n+5", "primes",
  "72 different gaps occur, from 94 to 7,050; the level share is 50.99 %; 1.2 % of terms are forced level (l <= d^2); there are no ties.",
  100000, 1, gen_a100760 },
{ "a101444", "A101444", "Numbers k such that (9973*k + 10007) is a prime", "prime values",
  "99 different gaps occur, from 2 to 212; the level share is 25.42 %; L = 2 holds 37 % of the level class; there are no ties.",
  100000, 1, gen_a101444 },
{ "a101567", "A101567", "Numbers n such that 1009*n + 10007 is prime", "prime values",
  "84 different gaps occur, from 2 to 200; the level share is 24.76 %; L = 2 holds 39 % of the level class.",
  100000, 1, gen_a101567 },
{ "a101780", "A101780", "Primes of the form 100*n + 3", "primes",
  "58 different gaps occur, from 100 to 8,700; the level share is 53.06 %; 1.0 % of terms are forced level (l <= d^2); there are no ties.",
  100000, 1, gen_a101780 },
{ "a102130", "A102130", "Primes of the form 8*n^2 + 4*n + 1", "primes",
  "99,898 different gaps occur, from 28 to 1,268,174,400; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a102130 },
{ "a102166", "A102166", "Numbers n such that 2*n^2 + 11*n + 101 is prime", "prime values",
  "81 different gaps occur, from 2 to 214; the level share is 23.38 %; L = 2 holds 40 % of the level class; there are no ties.",
  100000, 1, gen_a102166 },
{ "a102339", "A102339", "Numbers k such that k*10^3 + 333 is prime", "prime values",
  "62 different gaps occur, from 1 to 120; the level share is 16.94 %; L = 1 holds 32 % of the level class.",
  100000, 1, gen_a102339 },
{ "a102343", "A102343", "Numbers k such that k*10^3 + 777 is prime", "prime values",
  "69 different gaps occur, from 1 to 78; the level share is 15.98 %; L = 1 holds 30 % of the level class.",
  100000, 1, gen_a102343 },
{ "a102649", "A102649", "Numbers n such that 11*n^2 + 11*n + 3 is prime", "prime values",
  "79 different gaps occur, from 3 to 270; the level share is 32.67 %; L = 1 holds 34 % of the level class.",
  100000, 1, gen_a102649 },
{ "a102657", "A102657", "Numbers k such that 11*k^2 + 11*k + 1 is prime", "prime values",
  "73 different gaps occur, from 1 to 81; the level share is 20.69 %; L = 1 holds 39 % of the level class.",
  100000, 1, gen_a102657 },
{ "a102732", "A102732", "Primes of the form 13n+5", "primes",
  "55 different gaps occur, from 26 to 2,106; the level share is 42.79 %; there are no ties.",
  100000, 1, gen_a102732 },
{ "a102734", "A102734", "Primes of the form 23n+5", "primes",
  "64 different gaps occur, from 46 to 3,450; the level share is 46.48 %; there are no ties.",
  100000, 1, gen_a102734 },
{ "a102851", "A102851", "Primes of the form 19n + 5", "primes",
  "62 different gaps occur, from 38 to 3,192; the level share is 45.09 %.",
  100000, 1, gen_a102851 },
{ "a103564", "A103564", "Primes p such that 3*p^2 + 2 is prime", "primes",
  "729 different gaps occur, from 4 to 4,020; the level share is 49.16 %.",
  100000, 1, gen_a103564 },
{ "a103776", "A103776", "Primes p such that 8*p^2 + 4*p + 1 is also prime", "primes",
  "547 different gaps occur, from 4 to 3,754; the level share is 43.62 %; 6 terms do not decompose.",
  100000, 1, gen_a103776 },
{ "a105057", "A105057", "Numbers n such that 10000 * n - 1 is prime", "prime values",
  "74 different gaps occur, from 1 to 84; the level share is 23.23 %; L = 1 holds 41 % of the level class.",
  100000, 1, gen_a105057 },
{ "a105059", "A105059", "Numbers n such that 100000n - 1 is prime", "prime values",
  "85 different gaps occur, from 1 to 104; the level share is 23.79 %; L = 1 holds 38 % of the level class.",
  100000, 1, gen_a105059 },
{ "a105107", "A105107", "Numbers n such that 10000n + 1001 is prime", "prime values",
  "75 different gaps occur, from 1 to 88; the level share is 21.57 %; L = 1 holds 40 % of the level class.",
  100000, 1, gen_a105107 },
{ "a105126", "A105126", "Primes of the form 16n+9", "primes",
  "65 different gaps occur, from 16 to 1,392; the level share is 38.10 %.",
  100000, 1, gen_a105126 },
{ "a105127", "A105127", "Primes of the form 32n+17", "primes",
  "70 different gaps occur, from 32 to 2,400; the level share is 43.35 %.",
  100000, 1, gen_a105127 },
{ "a105128", "A105128", "Primes of the form 64n+33", "primes",
  "69 different gaps occur, from 64 to 4,864; the level share is 47.91 %.",
  100000, 1, gen_a105128 },
{ "a105129", "A105129", "Primes of the form 128n+65", "primes",
  "69 different gaps occur, from 128 to 9,984; the level share is 53.03 %; 1.8 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a105129 },
{ "a105130", "A105130", "Primes of the form 256n+129", "primes",
  "77 different gaps occur, from 256 to 20,736; the level share is 58.19 %; 4.1 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a105130 },
{ "a105131", "A105131", "Primes of the form 512n+257", "primes",
  "77 different gaps occur, from 512 to 48,128; the level share is 64.20 %; 8.5 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a105131 },
{ "a105132", "A105132", "Primes of the form 1024n + 513", "primes",
  "81 different gaps occur, from 1,024 to 136,192; the level share is 70.44 %; 16.2 % of terms are forced level (l <= d^2).",
  100000, 1, gen_a105132 },
{ "a105140", "A105140", "Numbers n such that 1024n+513 is prime", "prime values",
  "81 different gaps occur, from 1 to 133; the level share is 20.02 %.",
  100000, 1, gen_a105140 },
{ "a105680", "A105680", "Numbers k such that 1009*k + 9973 is prime", "prime values",
  "86 different gaps occur, from 2 to 222; the level share is 24.77 %; L = 2 holds 38 % of the level class; there are no ties.",
  100000, 1, gen_a105680 },
{ "a105710", "A105710", "Numbers k such that 9973*k + 1009 is prime", "prime values",
  "93 different gaps occur, from 2 to 226; the level share is 25.33 %; L = 2 holds 38 % of the level class.",
  100000, 1, gen_a105710 },
{ "a105854", "A105854", "Primes of the form 20*k + 3", "primes",
  "49 different gaps occur, from 20 to 1,080; the level share is 42.25 %; there are no ties.",
  100000, 1, gen_a105854 },
{ "a106483", "A106483", "Primes p such that 2*p^2 - 1 is also prime", "primes",
  "550 different gaps occur, from 1 to 1,674; the level share is 39.87 %.",
  100000, 1, gen_a106483 },
{ "a106699", "A106699", "Numbers k such that 1009*k - 9973 is prime", "prime values",
  "85 different gaps occur, from 2 to 210; the level share is 24.69 %; L = 2 holds 38 % of the level class.",
  100000, 1, gen_a106699 },
{ "a106700", "A106700", "Numbers k such that 9973*k - 1009 is prime", "prime values",
  "97 different gaps occur, from 2 to 258; the level share is 25.29 %; L = 2 holds 37 % of the level class; there are no ties.",
  100000, 1, gen_a106700 },
{ "a107003", "A107003", "Primes of the form 24*k + 5", "primes",
  "42 different gaps occur, from 24 to 1,056; the level share is 44.94 %; L = 1 holds 36 % of the level class; there are no ties.",
  100000, 1, gen_a107003 },
{ "a107071", "A107071", "Numbers k such that 1019*k + 1021 is prime", "prime values",
  "86 different gaps occur, from 2 to 216; the level share is 24.42 %; L = 2 holds 38 % of the level class; there are no ties.",
  100000, 1, gen_a107071 },
{ "a107072", "A107072", "Numbers k such that 1021*k + 1019 is prime", "prime values",
  "87 different gaps occur, from 2 to 212; the level share is 24.66 %; L = 2 holds 37 % of the level class; there are no ties.",
  100000, 1, gen_a107072 },
{ "a107301", "A107301", "Numbers k such that 10007*k + 99991 is prime", "prime values",
  "94 different gaps occur, from 2 to 220; the level share is 25.22 %; L = 2 holds 37 % of the level class.",
  100000, 1, gen_a107301 },
{ "a107302", "A107302", "Numbers k such that 99991*k + 10007 is prime", "prime values",
  "104 different gaps occur, from 2 to 342; the level share is 25.88 %; L = 2 holds 36 % of the level class.",
  100000, 1, gen_a107302 },
{ "a107303", "A107303", "Numbers k such that (3*k - 5) is prime", "prime values",
  "39 different gaps occur, from 2 to 90; the level share is 16.26 %; L = 2 holds 48 % of the level class.",
  100000, 1, gen_a107303 },
{ "a107306", "A107306", "Numbers k such that (17*k - 19) is prime", "prime values",
  "62 different gaps occur, from 2 to 134; the level share is 23.17 %; L = 2 holds 43 % of the level class.",
  100000, 1, gen_a107306 },
{ "a107308", "A107308", "Numbers k such that (29*k - 31) is prime", "prime values",
  "65 different gaps occur, from 2 to 136; the level share is 23.86 %; L = 2 holds 42 % of the level class.",
  100000, 1, gen_a107308 },
{ "a108058", "A108058", "Numbers k such that 179*k + 181 is prime", "prime values",
  "77 different gaps occur, from 2 to 204; the level share is 24.27 %; L = 2 holds 40 % of the level class; there are no ties.",
  100000, 1, gen_a108058 },
{ "a108059", "A108059", "Numbers k such that 181*k + 179 is prime", "prime values",
  "76 different gaps occur, from 2 to 172; the level share is 24.20 %; L = 2 holds 39 % of the level class.",
  100000, 1, gen_a108059 },
{ "a108060", "A108060", "Numbers k such that 191*k + 193 is prime", "prime values",
  "78 different gaps occur, from 2 to 214; the level share is 24.80 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a108060 },
{ "a108061", "A108061", "Numbers k such that 193*k + 191 is prime", "prime values",
  "75 different gaps occur, from 2 to 184; the level share is 24.38 %; L = 2 holds 38 % of the level class.",
  100000, 1, gen_a108061 },
{ "a108187", "A108187", "Numbers n such that 11*n - 5 is prime", "prime values",
  "60 different gaps occur, from 2 to 150; the level share is 20.78 %; L = 2 holds 42 % of the level class.",
  100000, 1, gen_a108187 },
{ "a108232", "A108232", "Numbers n such that 11*n - 7 is prime", "prime values",
  "57 different gaps occur, from 2 to 118; the level share is 21.62 %; L = 2 holds 43 % of the level class.",
  100000, 1, gen_a108232 },
{ "a108233", "A108233", "Numbers n such that 11*n + 5 is prime", "prime values",
  "58 different gaps occur, from 2 to 130; the level share is 20.45 %; L = 2 holds 42 % of the level class; there are no ties.",
  100000, 1, gen_a108233 },
{ "a108341", "A108341", "Numbers n such that 997*n - 1009 is prime", "prime values",
  "84 different gaps occur, from 2 to 186; the level share is 24.80 %; L = 2 holds 38 % of the level class; there are no ties.",
  100000, 1, gen_a108341 },
{ "a108342", "A108342", "Numbers n such that 1009*n - 997 is prime", "prime values",
  "85 different gaps occur, from 2 to 222; the level share is 24.58 %; L = 2 holds 37 % of the level class.",
  100000, 1, gen_a108342 },
{ "a108584", "A108584", "Numbers k such that 10*k - 97 is prime", "prime values",
  "50 different gaps occur, from 1 to 54; the level share is 21.84 %; L = 1 holds 49 % of the level class.",
  100000, 1, gen_a108584 },
{ "a108588", "A108588", "Numbers k such that 10*k + 97 is prime", "prime values",
  "48 different gaps occur, from 1 to 50; the level share is 20.68 %; L = 1 holds 45 % of the level class.",
  100000, 1, gen_a108588 },
{ "a108594", "A108594", "Numbers k such that 10*k + 101 is prime", "prime values",
  "47 different gaps occur, from 1 to 50; the level share is 22.08 %; L = 1 holds 45 % of the level class.",
  100000, 1, gen_a108594 },
{ "a108595", "A108595", "Numbers k such that 10*k + 103 is prime", "prime values",
  "50 different gaps occur, from 1 to 54; the level share is 21.82 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a108595 },
{ "a108596", "A108596", "Numbers k such that 911*k - 7 is prime", "prime values",
  "85 different gaps occur, from 2 to 210; the level share is 23.94 %; L = 2 holds 38 % of the level class; there are no ties.",
  100000, 1, gen_a108596 },
{ "a108597", "A108597", "Numbers n such that 911*n - 11 is prime", "prime values",
  "84 different gaps occur, from 2 to 180; the level share is 24.27 %; L = 2 holds 38 % of the level class; there are no ties.",
  100000, 1, gen_a108597 },
{ "a108724", "A108724", "Numbers n such that 11*n + 17 is prime", "prime values",
  "60 different gaps occur, from 2 to 150; the level share is 22.37 %; L = 2 holds 43 % of the level class; there are no ties.",
  100000, 1, gen_a108724 },
{ "a108725", "A108725", "Numbers n such that 11*n + 19 is prime", "prime values",
  "57 different gaps occur, from 2 to 132; the level share is 22.58 %; L = 2 holds 42 % of the level class; there are no ties.",
  100000, 1, gen_a108725 },
{ "a108726", "A108726", "Numbers n such that 11*n + 29 is prime", "prime values",
  "57 different gaps occur, from 2 to 142; the level share is 23.44 %; L = 2 holds 42 % of the level class; there are no ties.",
  100000, 1, gen_a108726 },
{ "a108727", "A108727", "Numbers n such that 11*n + 31 is prime", "prime values",
  "56 different gaps occur, from 2 to 128; the level share is 22.91 %; L = 2 holds 44 % of the level class.",
  100000, 1, gen_a108727 },
{ "a108751", "A108751", "Numbers k such that 11*k - 911 is prime", "prime values",
  "57 different gaps occur, from 1 to 182; the level share is 22.91 %; L = 2 holds 43 % of the level class.",
  100000, 1, gen_a108751 },
{ "a108757", "A108757", "Numbers k such that 1000*k + 911 is prime", "prime values",
  "66 different gaps occur, from 1 to 75; the level share is 22.53 %; L = 1 holds 43 % of the level class.",
  100000, 1, gen_a108757 },
{ "a108762", "A108762", "Numbers n such that 911*n + 13 is prime", "prime values",
  "86 different gaps occur, from 2 to 204; the level share is 24.43 %; L = 2 holds 39 % of the level class.",
  100000, 1, gen_a108762 },
{ "a108854", "A108854", "Numbers k such that 10*k - 127 is prime", "prime values",
  "50 different gaps occur, from 1 to 54; the level share is 20.98 %; L = 1 holds 45 % of the level class.",
  100000, 1, gen_a108854 },
{ "a108855", "A108855", "Numbers n such that 10*n + 127 is prime", "prime values",
  "48 different gaps occur, from 1 to 50; the level share is 22.00 %; L = 1 holds 47 % of the level class.",
  100000, 1, gen_a108855 },
{ "a108856", "A108856", "Numbers k such that 10*k - 131 is prime", "prime values",
  "47 different gaps occur, from 1 to 51; the level share is 21.80 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a108856 },
{ "a108857", "A108857", "Numbers n such that 10*n + 131 is prime", "prime values",
  "47 different gaps occur, from 1 to 50; the level share is 21.24 %; L = 1 holds 45 % of the level class.",
  100000, 1, gen_a108857 },
{ "a108874", "A108874", "Numbers k such that 41*k + 43 is prime", "prime values",
  "72 different gaps occur, from 2 to 152; the level share is 23.80 %; L = 2 holds 41 % of the level class.",
  100000, 1, gen_a108874 },
{ "a108899", "A108899", "Numbers k such that 11*k + 2357 is prime", "prime values",
  "60 different gaps occur, from 2 to 126; the level share is 22.71 %; L = 2 holds 44 % of the level class; there are no ties.",
  100000, 1, gen_a108899 },
{ "a108900", "A108900", "Numbers k such that 2357*k + 11 is prime", "prime values",
  "84 different gaps occur, from 2 to 182; the level share is 24.52 %; L = 2 holds 39 % of the level class.",
  100000, 1, gen_a108900 },
{ "a108901", "A108901", "Numbers n such that 2357*n + 23 is prime", "prime values",
  "90 different gaps occur, from 2 to 228; the level share is 24.76 %; L = 2 holds 37 % of the level class; there are no ties.",
  100000, 1, gen_a108901 },
{ "a108902", "A108902", "Numbers k such that 23*k + 2357 is prime", "prime values",
  "64 different gaps occur, from 2 to 162; the level share is 23.62 %; L = 2 holds 41 % of the level class; there are no ties.",
  100000, 1, gen_a108902 },
{ "a108936", "A108936", "Numbers n such that 11*n + 911 is prime", "prime values",
  "56 different gaps occur, from 2 to 128; the level share is 22.74 %; L = 2 holds 44 % of the level class; there are no ties.",
  100000, 1, gen_a108936 },
{ "a108937", "A108937", "Numbers k such that 911*k + 11 is prime", "prime values",
  "85 different gaps occur, from 2 to 238; the level share is 24.58 %; L = 2 holds 39 % of the level class.",
  100000, 1, gen_a108937 },
{ "a108938", "A108938", "Numbers k such that 911*k + 7 is prime", "prime values",
  "81 different gaps occur, from 2 to 180; the level share is 23.61 %; L = 2 holds 38 % of the level class; there are no ties.",
  100000, 1, gen_a108938 },
{ "a108969", "A108969", "Numbers n such that 43*n + 41 is prime", "prime values",
  "66 different gaps occur, from 2 to 156; the level share is 23.86 %; L = 2 holds 41 % of the level class.",
  100000, 1, gen_a108969 },
{ "a108976", "A108976", "Numbers k such that 17*k + 19 is prime", "prime values",
  "58 different gaps occur, from 2 to 124; the level share is 23.15 %; L = 2 holds 42 % of the level class.",
  100000, 1, gen_a108976 },
{ "a108977", "A108977", "Numbers n such that 19*n + 17 is prime", "prime values",
  "67 different gaps occur, from 2 to 142; the level share is 23.61 %; L = 2 holds 42 % of the level class.",
  100000, 1, gen_a108977 },
{ "a108978", "A108978", "Numbers k such that 29*k + 31 is prime", "prime values",
  "65 different gaps occur, from 2 to 156; the level share is 23.07 %; L = 2 holds 41 % of the level class.",
  100000, 1, gen_a108978 },
{ "a108979", "A108979", "Numbers k such that 31*k + 29 is prime", "prime values",
  "66 different gaps occur, from 2 to 154; the level share is 23.05 %; L = 2 holds 41 % of the level class.",
  100000, 1, gen_a108979 },
{ "a109603", "A109603", "Numbers n such that 43*n - 41 is prime", "prime values",
  "71 different gaps occur, from 2 to 156; the level share is 24.24 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a109603 },
{ "a109604", "A109604", "Numbers n such that 41*n - 43 is prime", "prime values",
  "66 different gaps occur, from 2 to 150; the level share is 23.56 %; L = 2 holds 41 % of the level class; there are no ties.",
  100000, 1, gen_a109604 },
{ "a109605", "A109605", "Numbers n such that 100000n + 91111 is prime", "prime values",
  "80 different gaps occur, from 1 to 83; the level share is 23.81 %; L = 1 holds 38 % of the level class.",
  100000, 1, gen_a109605 },
{ "a110801", "A110801", "Numbers k such that 12k + 1 is prime", "prime values",
  "40 different gaps occur, from 1 to 48; the level share is 18.29 %; L = 1 holds 49 % of the level class.",
  100000, 1, gen_a110801 },
{ "a110913", "A110913", "Numbers n such that 23*n^2 - 49 is prime", "prime values",
  "79 different gaps occur, from 2 to 174; the level share is 20.73 %; L = 2 holds 38 % of the level class.",
  100000, 1, gen_a110913 },
{ "a110959", "A110959", "Numbers k such that 23*k^2 + 1 is prime", "prime values",
  "58 different gaps occur, from 6 to 414; the level share is 18.84 %; L = 6 holds 44 % of the level class; there are no ties.",
  100000, 1, gen_a110959 },
{ "a110960", "A110960", "Numbers n such that 23*n^2 + 4 is prime", "prime values",
  "55 different gaps occur, from 6 to 384; the level share is 28.31 %; L = 3 holds 55 % of the level class.",
  100000, 1, gen_a110960 },
{ "a110961", "A110961", "Numbers k such that 23*k^2 + 9 is prime", "prime values",
  "81 different gaps occur, from 2 to 236; the level share is 17.59 %; L = 2 holds 31 % of the level class.",
  100000, 1, gen_a110961 },
{ "a110964", "A110964", "Numbers k such that 23*k^2 + 16 is prime", "prime values",
  "56 different gaps occur, from 6 to 378; the level share is 28.32 %; L = 3 holds 55 % of the level class.",
  100000, 1, gen_a110964 },
{ "a110965", "A110965", "Numbers k such that 23*k^2 + 25 is prime", "prime values",
  "69 different gaps occur, from 6 to 522; the level share is 19.19 %; L = 6 holds 39 % of the level class; there are no ties.",
  100000, 1, gen_a110965 },
{ "a110966", "A110966", "Numbers k such that 23*k^2 + 36 is prime", "prime values",
  "72 different gaps occur, from 2 to 166; the level share is 27.70 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a110966 },
{ "a110967", "A110967", "Numbers k such that 23*k^2 + 49 is prime", "prime values",
  "69 different gaps occur, from 6 to 558; the level share is 19.57 %; L = 6 holds 40 % of the level class; there are no ties.",
  100000, 1, gen_a110967 },
{ "a110974", "A110974", "Numbers n such that 23*n^2 - 1 is prime", "prime values",
  "97 different gaps occur, from 2 to 262; the level share is 23.13 %; L = 2 holds 35 % of the level class.",
  100000, 1, gen_a110974 },
{ "a110994", "A110994", "Numbers n such that 23*n^2 - 4 is prime", "prime values",
  "90 different gaps occur, from 2 to 204; the level share is 32.21 %; L = 1 holds 48 % of the level class.",
  100000, 1, gen_a110994 },
{ "a110998", "A110998", "Numbers n such that 23*n^2 - 9 is prime", "prime values",
  "136 different gaps occur, from 2 to 306; the level share is 22.63 %.",
  100000, 1, gen_a110998 },
{ "a110999", "A110999", "Numbers n such that 23*n^2 - 16 is prime", "prime values",
  "91 different gaps occur, from 2 to 200; the level share is 32.01 %; L = 1 holds 49 % of the level class.",
  100000, 1, gen_a110999 },
{ "a111001", "A111001", "Numbers n such that 23*n^2 - 25 is prime", "prime values",
  "115 different gaps occur, from 2 to 296; the level share is 23.43 %; L = 2 holds 32 % of the level class.",
  100000, 1, gen_a111001 },
{ "a111040", "A111040", "Numbers n such that 2*n^2 + 9 is prime", "prime values",
  "77 different gaps occur, from 1 to 84; the level share is 17.87 %.",
  100000, 1, gen_a111040 },
{ "a111041", "A111041", "Numbers m such that 2*m^2 + 25 is prime", "prime values",
  "73 different gaps occur, from 3 to 237; the level share is 19.77 %; L = 3 holds 38 % of the level class.",
  100000, 1, gen_a111041 },
{ "a111051", "A111051", "Numbers m such that 3*m^2 + 1 is prime", "prime values",
  "73 different gaps occur, from 2 to 234; the level share is 21.51 %; L = 2 holds 38 % of the level class.",
  100000, 1, gen_a111051 },
{ "a111052", "A111052", "Numbers m such that 3*m^2 + 4 is prime", "prime values",
  "71 different gaps occur, from 2 to 186; the level share is 30.79 %; L = 1 holds 51 % of the level class.",
  100000, 1, gen_a111052 },
{ "a111068", "A111068", "Numbers k such that 3*k^2 + 16 is prime", "prime values",
  "73 different gaps occur, from 2 to 170; the level share is 30.78 %; L = 1 holds 51 % of the level class.",
  100000, 1, gen_a111068 },
{ "a111069", "A111069", "Numbers k such that 3*k^2 + 25 is prime", "prime values",
  "92 different gaps occur, from 2 to 218; the level share is 21.86 %; L = 2 holds 35 % of the level class.",
  100000, 1, gen_a111069 },
{ "a111082", "A111082", "Numbers n such that 3*n^2 + 49 is prime", "prime values",
  "65 different gaps occur, from 2 to 150; the level share is 18.80 %; L = 2 holds 42 % of the level class.",
  100000, 1, gen_a111082 },
{ "a111083", "A111083", "Numbers k such that 3*k^2 + 64 is prime", "prime values",
  "75 different gaps occur, from 2 to 208; the level share is 30.72 %; L = 1 holds 51 % of the level class.",
  100000, 1, gen_a111083 },
{ "a111094", "A111094", "Numbers k such that 18*k + 1 is prime", "prime values",
  "43 different gaps occur, from 1 to 50; the level share is 18.15 %; L = 1 holds 44 % of the level class.",
  100000, 1, gen_a111094 },
{ "a111147", "A111147", "Numbers k such that 5*k^2 + 1 is prime", "prime values",
  "70 different gaps occur, from 6 to 450; the level share is 20.72 %; L = 6 holds 40 % of the level class; there are no ties.",
  100000, 1, gen_a111147 },
{ "a111148", "A111148", "Numbers k such that 5*k^2 + 4 is prime", "prime values",
  "69 different gaps occur, from 6 to 510; the level share is 30.21 %; L = 3 holds 53 % of the level class.",
  100000, 1, gen_a111148 },
{ "a111149", "A111149", "Numbers k such that 5*k^2 + 9 is prime", "prime values",
  "96 different gaps occur, from 2 to 256; the level share is 19.84 %.",
  100000, 1, gen_a111149 },
{ "a111174", "A111174", "Numbers k such that 24*k + 1 is prime", "prime values",
  "41 different gaps occur, from 1 to 42; the level share is 18.57 %; L = 1 holds 47 % of the level class.",
  100000, 1, gen_a111174 },
{ "a111175", "A111175", "Numbers k such that 30*k + 1 is prime", "prime values",
  "35 different gaps occur, from 1 to 35; the level share is 16.84 %; L = 1 holds 50 % of the level class.",
  100000, 1, gen_a111175 },
{ "a111251", "A111251", "Numbers k such that 3*k^2 + 3*k + 1 is prime", "prime values",
  "71 different gaps occur, from 1 to 95; the level share is 21.09 %; L = 1 holds 39 % of the level class.",
  100000, 1, gen_a111251 },
{ "a111292", "A111292", "Numbers n such that 6*n^2 + 6*n + 1 is prime", "prime values",
  "60 different gaps occur, from 1 to 68; the level share is 19.75 %; L = 1 holds 43 % of the level class.",
  100000, 1, gen_a111292 },
{ "a111294", "A111294", "Numbers n such that 23*n + 2 is prime", "prime values",
  "63 different gaps occur, from 2 to 148; the level share is 31.41 %; L = 1 holds 57 % of the level class.",
  100000, 1, gen_a111294 },
{ "a111312", "A111312", "Numbers n such that 11*n + 2 is prime", "prime values",
  "57 different gaps occur, from 1 to 182; the level share is 30.66 %; L = 1 holds 59 % of the level class.",
  100000, 1, gen_a111312 },
{ "a111369", "A111369", "Numbers k such that 13*k + 11 is prime", "prime values",
  "59 different gaps occur, from 2 to 138; the level share is 22.61 %; L = 2 holds 43 % of the level class; there are no ties.",
  100000, 1, gen_a111369 },
{ "a111455", "A111455", "Numbers k such that 101*k + 97 is prime", "prime values",
  "74 different gaps occur, from 2 to 170; the level share is 24.14 %; L = 2 holds 40 % of the level class; there are no ties.",
  100000, 1, gen_a111455 },
{ "a113151", "A113151", "Primes p such that 19*p + 2 is also prime", "primes",
  "285 different gaps occur, from 2 to 2,142; the level share is 47.94 %; L = 1 holds 36 % of the level class; there are no ties.",
  100000, 1, gen_a113151 },
{ "a113487", "A113487", "Numbers k such that 17*k + 2 is prime", "prime values",
  "59 different gaps occur, from 1 to 124; the level share is 33.37 %; L = 1 holds 58 % of the level class.",
  100000, 1, gen_a113487 },
{ "a113488", "A113488", "Numbers k such that 19*k + 2 is prime", "prime values",
  "62 different gaps occur, from 2 to 142; the level share is 32.93 %; L = 1 holds 58 % of the level class.",
  100000, 1, gen_a113488 },
{ "a113510", "A113510", "Numbers k such that 29*k + 2 is prime", "prime values",
  "66 different gaps occur, from 1 to 156; the level share is 33.54 %; L = 1 holds 56 % of the level class.",
  100000, 1, gen_a113510 },
{ "a117047", "A117047", "Primes of the form 60*k + 11", "primes",
  "38 different gaps occur, from 60 to 2,340; the level share is 52.66 %; L = 1 holds 40 % of the level class; there are no ties.",
  100000, 1, gen_a117047 },
{ "a117049", "A117049", "Primes of the form 22*(n^2)+1", "primes",
  "99,361 different gaps occur, from 66 to 1,913,043,968; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a117049 },
{ "a119409", "A119409", "Numbers k such that 235*k + 1 is prime", "prime values",
  "62 different gaps occur, from 2 to 174; the level share is 23.06 %; L = 2 holds 42 % of the level class.",
  100000, 1, gen_a119409 },
{ "a120344", "A120344", "Numbers k such that 23*k + 1 is a prime", "prime values",
  "66 different gaps occur, from 2 to 144; the level share is 23.11 %; L = 2 holds 41 % of the level class.",
  100000, 1, gen_a120344 },
{ "a120345", "A120345", "Numbers n such that 2357*n + 1 is prime", "prime values",
  "89 different gaps occur, from 2 to 194; the level share is 25.22 %; L = 2 holds 38 % of the level class.",
  100000, 1, gen_a120345 },
{ "a121068", "A121068", "Numbers k such that 8*k^2 + 7 is prime", "prime values",
  "104 different gaps occur, from 3 to 354; the level share is 24.40 %; L = 3 holds 35 % of the level class.",
  100000, 1, gen_a121068 },
{ "a121817", "A121817", "Numbers m such that 23 + 36*m*(m+1) is prime", "prime values",
  "105 different gaps occur, from 1 to 133; the level share is 23.32 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_a121817 },
{ "a122114", "A122114", "Primes of the form 2n^2 + 26n + 1", "primes",
  "97,916 different gaps occur, from 32 to 68,420,740; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a122114 },
{ "a122430", "A122430", "Primes of the form 1+2*n+3*n^2", "primes",
  "99,938 different gaps occur, from 160 to 1,958,965,080; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a122430 },
{ "a122482", "A122482", "Primes p such that 1 + 4p + 12p^2 is prime", "primes",
  "368 different gaps occur, from 6 to 3,480; the level share is 48.96 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_a122482 },
{ "a124127", "A124127", "Numbers k such that 17k + 1 is prime", "prime values",
  "60 different gaps occur, from 2 to 140; the level share is 23.13 %; L = 2 holds 44 % of the level class; there are no ties.",
  100000, 1, gen_a124127 },
{ "a124198", "A124198", "Numbers k such that 21*k + 1 is prime", "prime values",
  "34 different gaps occur, from 2 to 98; the level share is 18.63 %; L = 2 holds 51 % of the level class.",
  100000, 1, gen_a124198 },
{ "a124204", "A124204", "Numbers k such that 20*k + 1 is prime", "prime values",
  "49 different gaps occur, from 1 to 56; the level share is 21.32 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a124204 },
{ "a126332", "A126332", "Numbers k such that 10k + 13 is prime", "prime values",
  "50 different gaps occur, from 1 to 54; the level share is 20.76 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a126332 },
{ "a126785", "A126785", "Numbers k such that 10*k + 11 is prime", "prime values",
  "47 different gaps occur, from 1 to 50; the level share is 20.99 %; L = 1 holds 47 % of the level class.",
  100000, 1, gen_a126785 },
{ "a126960", "A126960", "Primes p such that (3p)^2 + 2 is prime", "primes",
  "468 different gaps occur, from 2 to 1,694; the level share is 38.52 %.",
  100000, 1, gen_a126960 },
{ "a127435", "A127435", "Primes p such that (p-1)^2 + 1 is prime", "primes",
  "515 different gaps occur, from 1 to 2,480; the level share is 42.79 %; 7 terms do not decompose.",
  100000, 1, gen_a127435 },
{ "a127575", "A127575", "Numbers n such that 16n+15 is prime", "prime values",
  "61 different gaps occur, from 1 to 87; the level share is 16.21 %.",
  100000, 1, gen_a127575 },
{ "a127576", "A127576", "Primes of the form 16n+15", "primes",
  "61 different gaps occur, from 16 to 1,392; the level share is 38.28 %; there are no ties.",
  100000, 1, gen_a127576 },
{ "a127579", "A127579", "Primes of the form 64n+63", "primes",
  "71 different gaps occur, from 64 to 4,992; the level share is 47.76 %; there are no ties.",
  100000, 1, gen_a127579 },
{ "a127580", "A127580", "Numbers k such that 64k+63 is prime", "prime values",
  "71 different gaps occur, from 1 to 78; the level share is 17.57 %; L = 1 holds 30 % of the level class.",
  100000, 1, gen_a127580 },
{ "a127589", "A127589", "Primes of the form 16k + 5", "primes",
  "60 different gaps occur, from 16 to 1,104; the level share is 38.51 %; there are no ties.",
  100000, 1, gen_a127589 },
{ "a127590", "A127590", "Numbers n such that 16n+5 is prime", "prime values",
  "60 different gaps occur, from 1 to 69; the level share is 21.18 %; L = 1 holds 41 % of the level class.",
  100000, 1, gen_a127590 },
{ "a127591", "A127591", "Numbers k such that 64k+21 is prime", "prime values",
  "69 different gaps occur, from 1 to 75; the level share is 17.43 %; L = 1 holds 30 % of the level class.",
  100000, 1, gen_a127591 },
{ "a127592", "A127592", "Primes of the form 64k+21", "primes",
  "69 different gaps occur, from 64 to 4,800; the level share is 47.72 %; there are no ties.",
  100000, 1, gen_a127592 },
{ "a127593", "A127593", "Primes of the form 256 k + 85", "primes",
  "79 different gaps occur, from 256 to 23,552; the level share is 58.29 %; 4.0 % of terms are forced level (l <= d^2); there are no ties.",
  100000, 1, gen_a127593 },
{ "a127594", "A127594", "Numbers k such that 256 k + 85 is prime", "prime values",
  "79 different gaps occur, from 1 to 92; the level share is 22.32 %; L = 1 holds 38 % of the level class.",
  100000, 1, gen_a127594 },
{ "a128829", "A128829", "Numbers k such that 6*k^2 + 17 is prime", "prime values",
  "46 different gaps occur, from 1 to 50; the level share is 17.59 %; L = 1 holds 48 % of the level class.",
  100000, 1, gen_a128829 },
{ "a129484", "A129484", "Primes of the form 17k + 1", "primes",
  "60 different gaps occur, from 34 to 2,380; the level share is 44.48 %.",
  100000, 1, gen_a129484 },
{ "a039787", "A039787", "Primes p such that p-1 is squarefree", "primes",
  "78 different gaps occur, from 1 to 360; the level share is 29.95 %; there are no ties.",
  100000, 1, gen_a039787 },
{ "a049097", "A049097", "Primes p such that p+1 is squarefree", "primes",
  "81 different gaps occur, from 3 to 400; the level share is 31.49 %.",
  100000, 1, gen_a049097 },
{ "a049231", "A049231", "Primes p such that p - 2 is squarefree", "primes",
  "71 different gaps occur, from 2 to 180; the level share is 24.84 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_a049231 },
{ "a049233", "A049233", "Primes p such that p + 2 is squarefree", "primes",
  "72 different gaps occur, from 2 to 196; the level share is 25.84 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a049233 },
{ "a085722", "A085722", "Numbers k such that k^2 + 1 is a semiprime", "multiplicative",
  "37 different gaps occur, from 1 to 40; the level share is 16.73 %; L = 2 holds 36 % of the level class.",
  100000, 1, gen_a085722 },
{ "a085746", "A085746", "Numbers n such that n^2 + n + 1 is a semiprime", "multiplicative",
  "31 different gaps occur, from 1 to 31; the level share is 17.21 %; L = 1 holds 55 % of the level class.",
  100000, 1, gen_a085746 },
{ "a092109", "A092109", "Primes p such that p+3 is a semiprime", "primes",
  "340 different gaps occur, from 4 to 2,328; the level share is 40.93 %; there are no ties.",
  100000, 1, gen_a092109 },
{ "a108181", "A108181", "Semiprimes of the form 4n + 1", "multiplicative",
  "24 different gaps occur, from 4 to 96; the level share is 26.20 %; L = 1 holds 43 % of the level class.",
  100000, 1, gen_a108181 },
{ "a108769", "A108769", "Numbers m such that m^2 + (m+1)^2 is a semiprime", "multiplicative",
  "29 different gaps occur, from 1 to 32; the level share is 15.84 %; L = 1 holds 56 % of the level class.",
  100000, 1, gen_a108769 },
{ "a109953", "A109953", "Primes p such that p^2+2 is a semiprime", "primes",
  "629 different gaps occur, from 2 to 1,916; the level share is 41.16 %.",
  100000, 1, gen_a109953 },
{ "a112771", "A112771", "Semiprimes of the form 6n + 1", "multiplicative",
  "19 different gaps occur, from 6 to 114; the level share is 33.55 %; L = 1 holds 64 % of the level class.",
  100000, 1, gen_a112771 },
{ "a112772", "A112772", "Semiprimes of the form 6n+2", "multiplicative",
  "39 different gaps occur, from 12 to 540; the level share is 35.25 %; L = 2 holds 51 % of the level class; there are no ties.",
  100000, 1, gen_a112772 },
{ "a112774", "A112774", "Semiprimes of the form 6n+4", "multiplicative",
  "41 different gaps occur, from 6 to 516; the level share is 35.17 %; L = 2 holds 51 % of the level class; there are no ties.",
  100000, 1, gen_a112774 },
{ "a112775", "A112775", "Numbers k such that 6k+1 is semiprime", "multiplicative",
  "19 different gaps occur, from 1 to 19; the level share is 14.18 %; L = 1 holds 64 % of the level class.",
  100000, 1, gen_a112775 },
{ "a112776", "A112776", "Numbers k such that 6k+5 is semiprime", "multiplicative",
  "20 different gaps occur, from 1 to 20; the level share is 13.26 %; L = 1 holds 64 % of the level class.",
  100000, 1, gen_a112776 },
{ "a112777", "A112777", "Numbers k such that 2*k^2 + 1 is a semiprime", "multiplicative",
  "40 different gaps occur, from 1 to 45; the level share is 17.90 %; L = 1 holds 41 % of the level class.",
  100000, 1, gen_a112777 },
{ "a122488", "A122488", "Numbers k such that 1 + 2k + 3k^2 is semiprime", "multiplicative",
  "42 different gaps occur, from 1 to 57; the level share is 17.65 %; L = 2 holds 36 % of the level class.",
  100000, 1, gen_a122488 },
{ "a033253", "A033253", "Primes of form x^2+83*y^2", "quadratic form",
  "1,009 different gaps occur, from 2 to 3,378; the level share is 44.13 %.",
  100000, 1, gen_a033253 },
{ "a033254", "A033254", "Primes of form x^2+85*y^2", "quadratic form",
  "140 different gaps occur, from 8 to 1,332; the level share is 40.52 %.",
  100000, 1, gen_a033254 },
{ "a033255", "A033255", "Primes of form x^2+86*y^2", "quadratic form",
  "862 different gaps occur, from 2 to 3,344; the level share is 45.25 %.",
  100000, 1, gen_a033255 },
{ "a033256", "A033256", "Primes of form x^2+87*y^2", "quadratic form",
  "259 different gaps occur, from 6 to 2,028; the level share is 47.38 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a033256 },
{ "a033257", "A033257", "Primes of form x^2+89*y^2", "quadratic form",
  "730 different gaps occur, from 4 to 4,088; the level share is 45.57 %.",
  100000, 1, gen_a033257 },
{ "a033258", "A033258", "Primes of form x^2+91*y^2", "quadratic form",
  "665 different gaps occur, from 2 to 1,848; the level share is 40.30 %.",
  100000, 1, gen_a033258 },
{ "a106880", "A106880", "Primes of the form x^2+xy+9y^2, with x and y nonnegative", "quadratic form",
  "148 different gaps occur, from 2 to 690; the level share is 35.50 %.",
  100000, 1, gen_a106880 },
{ "a106890", "A106890", "Primes of the form x^2 + xy + 11y^2, with x and y nonnegative", "quadratic form",
  "126 different gaps occur, from 2 to 318; the level share is 29.70 %.",
  100000, 1, gen_a106890 },
{ "a106900", "A106900", "Primes of the form x^2+xy+12y^2, with x and y nonnegative", "quadratic form",
  "646 different gaps occur, from 2 to 1,830; the level share is 40.15 %.",
  100000, 1, gen_a106900 },
{ "a106901", "A106901", "Primes of the form 3x^2+3xy+5y^2, with x and y nonnegative", "quadratic form",
  "112 different gaps occur, from 2 to 756; the level share is 40.05 %; L = 1 holds 40 % of the level class; there are no ties.",
  100000, 1, gen_a106901 },
{ "a106903", "A106903", "Primes of the form x^2+xy+13y^2, with x and y nonnegative", "quadratic form",
  "87 different gaps occur, from 6 to 666; the level share is 39.03 %; L = 1 holds 41 % of the level class.",
  100000, 1, gen_a106903 },
{ "a106905", "A106905", "Primes of the form 2x^2+2xy+7y^2, with x and y nonnegative", "quadratic form",
  "148 different gaps occur, from 4 to 780; the level share is 36.70 %; there are no ties; 6 terms do not decompose.",
  100000, 1, gen_a106905 },
{ "a106907", "A106907", "Primes of the form 4x^2+3xy+4y^2, with x and y nonnegative", "quadratic form",
  "741 different gaps occur, from 2 to 4,370; the level share is 49.74 %.",
  100000, 1, gen_a106907 },
{ "a106921", "A106921", "Primes of the form x^2+xy+15y^2, with x and y nonnegative", "quadratic form",
  "390 different gaps occur, from 2 to 1,156; the level share is 37.52 %.",
  100000, 1, gen_a106921 },
{ "a106926", "A106926", "Primes of the form 2x^2+xy+8y^2, with x and y nonnegative", "quadratic form",
  "183 different gaps occur, from 6 to 1,650; the level share is 42.55 %; L = 1 holds 30 % of the level class; there are no ties.",
  100000, 1, gen_a106926 },
{ "a106929", "A106929", "Primes of the form x^2+xy+16y^2, with x and y nonnegative", "quadratic form",
  "192 different gaps occur, from 6 to 1,470; the level share is 42.19 %; L = 1 holds 30 % of the level class.",
  100000, 1, gen_a106929 },
{ "a106931", "A106931", "Primes of the form 4x^2+4xy+5y^2, with x and y nonnegative", "quadratic form",
  "92 different gaps occur, from 8 to 1,024; the level share is 35.85 %; there are no ties.",
  100000, 1, gen_a106931 },
{ "a106932", "A106932", "Primes of the form x^2 + xy + 17y^2, with x and y nonnegative", "quadratic form",
  "124 different gaps occur, from 2 to 300; the level share is 29.48 %.",
  100000, 1, gen_a106932 },
{ "a106934", "A106934", "Primes of the form 3x^2+2xy+6y^2, with x and y nonnegative", "quadratic form",
  "287 different gaps occur, from 4 to 1,428; the level share is 38.06 %; there are no ties.",
  100000, 1, gen_a106934 },
{ "a106937", "A106937", "Primes of the form 2x^2+2xy+9y^2, with x and y nonnegative", "quadratic form",
  "293 different gaps occur, from 4 to 1,408; the level share is 38.18 %.",
  100000, 1, gen_a106937 },
{ "a106939", "A106939", "Primes of the form 4x^2+3xy+5y^2, with x and y nonnegative", "quadratic form",
  "1,023 different gaps occur, from 2 to 3,834; the level share is 43.42 %.",
  100000, 1, gen_a106939 },
{ "a106945", "A106945", "Primes of the form 2x^2+xy+9y^2, with x and y nonnegative", "quadratic form",
  "858 different gaps occur, from 2 to 2,526; the level share is 42.34 %.",
  100000, 1, gen_a106945 },
{ "a106949", "A106949", "Primes of the form 2x^2 + 9y^2", "quadratic form",
  "61 different gaps occur, from 6 to 534; the level share is 39.73 %; L = 1 holds 43 % of the level class; there are no ties.",
  100000, 1, gen_a106949 },
{ "a106950", "A106950", "Primes of the form x^2 + 18y^2", "quadratic form",
  "57 different gaps occur, from 6 to 534; the level share is 39.60 %; L = 1 holds 43 % of the level class.",
  100000, 1, gen_a106950 },
{ "a106951", "A106951", "Primes of the form 3x^2+3xy+7y^2, with x and y nonnegative", "quadratic form",
  "69 different gaps occur, from 4 to 930; the level share is 44.53 %; L = 1 holds 51 % of the level class.",
  100000, 1, gen_a106951 },
{ "a106953", "A106953", "Primes of the form 4x^2+2xy+5y^2, with x and y nonnegative", "quadratic form",
  "400 different gaps occur, from 2 to 1,026; the level share is 38.58 %.",
  100000, 1, gen_a106953 },
{ "a106959", "A106959", "Primes of the form 2x^2+xy+10y^2, with x and y nonnegative", "quadratic form",
  "620 different gaps occur, from 2 to 1,932; the level share is 39.82 %.",
  100000, 1, gen_a106959 },
{ "a106964", "A106964", "Primes of the form 3x^2+2xy+7y^2, with x and y nonnegative", "quadratic form",
  "175 different gaps occur, from 4 to 1,604; the level share is 42.96 %; there are no ties.",
  100000, 1, gen_a106964 },
{ "a106966", "A106966", "Primes of the form 3x^2+xy+7y^2, with x and y nonnegative", "quadratic form",
  "376 different gaps occur, from 2 to 1,238; the level share is 37.32 %; 6 terms do not decompose.",
  100000, 1, gen_a106966 },
{ "a106969", "A106969", "Primes of the form x^2+xy+21y^2, with x and y nonnegative", "quadratic form",
  "384 different gaps occur, from 2 to 1,010; the level share is 37.23 %.",
  100000, 1, gen_a106969 },
{ "a106971", "A106971", "Primes of the form 5x^2+4xy+5y^2, with x and y nonnegative", "quadratic form",
  "245 different gaps occur, from 12 to 4,260; the level share is 48.33 %; there are no ties.",
  100000, 1, gen_a106971 },
{ "a106974", "A106974", "Primes of the form 2x^2+2xy+11y^2, with x and y nonnegative", "quadratic form",
  "99 different gaps occur, from 9 to 1,596; the level share is 43.01 %; there are no ties.",
  100000, 1, gen_a106974 },
{ "a106975", "A106975", "Primes of the form 4x^2+3xy+6y^2, with x and y nonnegative", "quadratic form",
  "316 different gaps occur, from 6 to 2,598; the level share is 48.94 %; L = 1 holds 35 % of the level class.",
  100000, 1, gen_a106975 },
{ "a106978", "A106978", "Primes of the form 3x^2+3xy+8y^2, with x and y nonnegative", "quadratic form",
  "320 different gaps occur, from 6 to 2,838; the level share is 48.92 %; L = 1 holds 35 % of the level class; there are no ties.",
  100000, 1, gen_a106978 },
{ "a106980", "A106980", "Primes of the form 2x^2+xy+11y^2, with x and y nonnegative", "quadratic form",
  "279 different gaps occur, from 6 to 2,070; the level share is 47.85 %; L = 1 holds 36 % of the level class; there are no ties.",
  100000, 1, gen_a106980 },
{ "a106984", "A106984", "Primes of the form 2x^2 + 11y^2", "quadratic form",
  "171 different gaps occur, from 2 to 550; the level share is 35.23 %.",
  100000, 1, gen_a106984 },
{ "a106985", "A106985", "Primes of the form 5x^2+3xy+5y^2, with x and y nonnegative", "quadratic form",
  "566 different gaps occur, from 2 to 1,656; the level share is 38.69 %.",
  100000, 1, gen_a106985 },
{ "a106988", "A106988", "Primes of the form x^2+xy+23y^2, with x and y nonnegative", "quadratic form",
  "229 different gaps occur, from 2 to 714; the level share is 33.27 %.",
  100000, 1, gen_a106988 },
{ "a106990", "A106990", "Primes of the form 5x^2+5xy+6y^2, with x and y nonnegative", "quadratic form",
  "801 different gaps occur, from 2 to 5,800; the level share is 50.24 %.",
  100000, 1, gen_a106990 },
{ "a106992", "A106992", "Primes of the form 4x^2+xy+6y^2, with x and y nonnegative", "quadratic form",
  "602 different gaps occur, from 2 to 3,192; the level share is 47.94 %.",
  100000, 1, gen_a106992 },
{ "a106995", "A106995", "Primes of the form 3x^2+xy+8y^2, with x and y nonnegative", "quadratic form",
  "597 different gaps occur, from 4 to 2,894; the level share is 47.99 %.",
  100000, 1, gen_a106995 },
{ "a106998", "A106998", "Primes of the form 2x^2+xy+12y^2, with x and y nonnegative", "quadratic form",
  "608 different gaps occur, from 4 to 2,940; the level share is 47.83 %.",
  100000, 1, gen_a106998 },
{ "a107002", "A107002", "Primes of the form 5x^2+2xy+5y^2, with x and y nonnegative", "quadratic form",
  "106 different gaps occur, from 24 to 3,168; the level share is 49.73 %; L = 1 holds 32 % of the level class; there are no ties.",
  100000, 1, gen_a107002 },
{ "a107005", "A107005", "Primes of the form 4x^2+4xy+7y^2, with x and y nonnegative", "quadratic form",
  "62 different gaps occur, from 24 to 1,968; the level share is 46.27 %; L = 1 holds 34 % of the level class; there are no ties.",
  100000, 1, gen_a107005 },
{ "a107007", "A107007", "Primes of the form 3*x^2+8*y^2", "quadratic form",
  "42 different gaps occur, from 8 to 1,056; the level share is 44.34 %; L = 1 holds 36 % of the level class; there are no ties.",
  100000, 1, gen_a107007 },
{ "a107008", "A107008", "Primes of the form x^2 + 24*y^2", "quadratic form",
  "41 different gaps occur, from 24 to 1,008; the level share is 44.44 %; L = 1 holds 36 % of the level class.",
  100000, 1, gen_a107008 },
{ "a107009", "A107009", "Primes of the form 5x^2+xy+5y^2, with x and y nonnegative", "quadratic form",
  "182 different gaps occur, from 6 to 1,320; the level share is 46.39 %; L = 1 holds 40 % of the level class; there are no ties.",
  100000, 1, gen_a107009 },
{ "a107012", "A107012", "Primes of the form x^2+xy+25y^2, with x and y nonnegative", "quadratic form",
  "86 different gaps occur, from 6 to 582; the level share is 42.05 %; L = 1 holds 46 % of the level class.",
  100000, 1, gen_a107012 },
{ "a107132", "A107132", "Primes of the form 2x^2 + 13y^2", "quadratic form",
  "510 different gaps occur, from 2 to 2,242; the level share is 42.97 %.",
  100000, 1, gen_a107132 },
{ "a107133", "A107133", "Primes of the form 4x^2 + 7y^2", "quadratic form",
  "114 different gaps occur, from 4 to 580; the level share is 30.41 %; there are no ties.",
  100000, 1, gen_a107133 },
{ "a127736", "A127736", "a(n) = n*(n^2 + 2*n - 1)/2", "polynomial",
  "Every gap is different, from 6 to 15,000,350,001; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a127736 },
{ "a006000", "A006000", "a(n) = (n+1)*(n^2+n+2)/2", "polynomial",
  "Every gap is different, from 3 to 15,000,050,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a006000 },
{ "a016061", "A016061", "a(n) = n*(n+1)*(4*n+5)/6", "polynomial",
  "Every gap is different, from 3 to 20,000,100,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a016061 },
{ "a027620", "A027620", "a(n) = n + (n+1)^2 + (n+2)^3", "polynomial",
  "Every gap is different, from 23 to 30,001,100,009; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a027620 },
{ "a033994", "A033994", "a(n) = n*(n+1)*(5*n+1)/6", "polynomial",
  "Every gap is different, from 9 to 25,000,450,002; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a033994 },
{ "a035328", "A035328", "a(n) = n*(2*n-1)*(2*n+1)", "polynomial",
  "Every gap is different, from 3 to 119,998,800,003; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a035328 },
{ "a035329", "A035329", "a(n) = n*(2*n+5)*(2*n+7)", "polynomial",
  "Every gap is different, from 63 to 120,003,600,015; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a035329 },
{ "a037235", "A037235", "a(n) = n*(2*n^2 - 3*n + 4)/3", "polynomial",
  "Every gap is different, from 1 to 19,999,600,003; every decomposable term is forced level (l <= d^2); 8 terms do not decompose.",
  100000, 0, gen_a037235 },
{ "a056578", "A056578", "a(n) = 1 + 2*n + 3*n^2 + 4*n^3", "polynomial",
  "Every gap is different, from 9 to 119,999,400,003; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a056578 },
{ "a071230", "A071230", "a(n) = n*(6*n^2 - 7*n + 3)/2", "polynomial",
  "Every gap is different, from 1 to 89,998,400,008; every decomposable term is forced level (l <= d^2); 8 terms do not decompose.",
  100000, 0, gen_a071230 },
{ "a086605", "A086605", "a(n) = 9*n^3 - 18*n^2 + 10*n", "polynomial",
  "Every gap is different, from 1 to 269,993,700,037; every decomposable term is forced level (l <= d^2); 8 terms do not decompose.",
  100000, 0, gen_a086605 },
{ "a101853", "A101853", "a(n) = n*(20 + 15*n + n^2)/6", "polynomial",
  "Every gap is different, from 12 to 5,000,550,006; every decomposable term is forced level (l <= d^2).",
  100000, 1, gen_a101853 },
{ "a102094", "A102094", "a(n) = (2*n-1)*(2*n+1)^2", "polynomial",
  "Every gap is different, from 66 to 240,003,200,010; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a102094 },
{ "a111144", "A111144", "a(n) = n*(n+13)*(n+14)/6", "polynomial",
  "Every gap is different, from 35 to 5,000,850,026; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a111144 },
{ "a115519", "A115519", "a(n) = n*(1+3*n+6*n^2)/2", "polynomial",
  "Every gap is different, from 5 to 89,999,400,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a115519 },
{ "a126335", "A126335", "a(n) = n*(4*n^2+5*n-3)/2", "polynomial",
  "Every gap is different, from 20 to 60,001,100,003; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a126335 },
{ "a074742", "A074742", "a(n) = (n^3 + 6n^2 - n + 12)/6", "polynomial",
  "Every gap is different, from 1 to 5,000,149,999; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a074742 },
{ "a100207", "A100207", "a(n) = 4 + 8*n + 10*n^2 + 4*n^3", "polynomial",
  "Every gap is different, from 22 to 120,000,800,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a100207 },
{ "a000297", "A000297", "a(n) = (n+1)*(n+3)*(n+8)/6", "polynomial",
  "Every gap is different, from 4 to 5,000,250,001; every decomposable term is forced level (l <= d^2).",
  100000, -1, gen_a000297 },
{ "a027602", "A027602", "a(n) = n^3 + (n+1)^3 + (n+2)^3", "polynomial",
  "Every gap is different, from 27 to 90,000,900,009; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a027602 },
{ "a027849", "A027849", "a(n) = (n+1)*(5*n^2+4*n+1)", "polynomial",
  "Every gap is different, from 19 to 150,000,300,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a027849 },
{ "a049480", "A049480", "a(n) = (2*n-1)*(n^2 -n +6)/6", "polynomial",
  "Every gap is different, from 3 to 10,000,000,002; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 1, gen_a049480 },
{ "a056520", "A056520", "a(n) = (n + 2)*(2*n^2 - n + 3)/6", "polynomial",
  "Every gap is different, from 1 to 10,000,000,000; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a056520 },
{ "a063488", "A063488", "a(n) = (2*n-1)*(n^2 -n +2)/2", "polynomial",
  "Every gap is different, from 5 to 30,000,000,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a063488 },
{ "a101165", "A101165", "a(n) = (7*n^3 + 6*n^2 + 5*n) / 6", "polynomial",
  "Every gap is different, from 3 to 34,999,850,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a101165 },
{ "a071229", "A071229", "a(n) = n*(14*n^2 - 21*n + 13)/6", "polynomial",
  "Every gap is different, from 1 to 69,998,600,008; every decomposable term is forced level (l <= d^2); 8 terms do not decompose.",
  100000, 0, gen_a071229 },
{ "a101860", "A101860", "a(n) = (3+n)*(2 + 33*n + n^2)/6", "polynomial",
  "Every gap is different, from 23 to 5,001,150,011; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a101860 },
{ "a114211", "A114211", "a(n) = (5*n^3+12*n^2+n+6)/6", "polynomial",
  "Every gap is different, from 3 to 25,000,149,999; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a114211 },
{ "a011199", "A011199", "a(n) = (n+1)*(2*n+1)*(3*n+1)", "polynomial",
  "Every gap is different, from 23 to 180,000,400,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a011199 },
{ "a079588", "A079588", "a(n) = (n+1)*(2*n+1)*(4*n+1)", "polynomial",
  "Every gap is different, from 29 to 240,000,400,001; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a079588 },
{ "a095796", "A095796", "a(n) = 1 + (26*n+17+7*n^2)*n/2", "polynomial",
  "Every gap is different, from 25 to 105,001,549,999; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a095796 },
{ "a100504", "A100504", "a(n) = (4*n^3 + 6*n^2 + 8*n + 6)/3", "polynomial",
  "Every gap is different, from 6 to 40,000,000,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a100504 },
{ "a034721", "A034721", "a(n) = (10*n^3 - 9*n^2 + 2*n)/3 + 1", "polynomial",
  "Every gap is different, from 1 to 99,998,400,007; every decomposable term is forced level (l <= d^2); 8 terms do not decompose.",
  100000, 0, gen_a034721 },
{ "a087863", "A087863", "a(n) = (n^3 + 24*n^2 + 65*n + 36)/6", "polynomial",
  "Every gap is different, from 15 to 5,000,750,007; every decomposable term is forced level (l <= d^2).",
  100000, 0, gen_a087863 },
{ "a057813", "A057813", "a(n) = (2*n+1)*(4*n^2+4*n+3)/3", "polynomial",
  "Every gap is different, from 10 to 80,000,000,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 0, gen_a057813 },
{ "a061550", "A061550", "a(n) = (2*n+1)*(2*n+3)*(2*n+5)", "polynomial",
  "Every gap is different, from 90 to 240,004,800,018; every decomposable term is forced level (l <= d^2); 6 terms do not decompose.",
  100000, 0, gen_a061550 },
{ "a063489", "A063489", "a(n) = (2*n-1)*(5*n^2-5*n+6)/6", "polynomial",
  "Every gap is different, from 7 to 50,000,000,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a063489 },
{ "a063490", "A063490", "a(n) = (2*n - 1)*(7*n^2 - 7*n + 6)/6", "polynomial",
  "Every gap is different, from 9 to 70,000,000,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a063490 },
{ "a063491", "A063491", "a(n) = (2*n - 1)*(3*n^2 - 3*n + 2)/2", "polynomial",
  "Every gap is different, from 11 to 90,000,000,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a063491 },
{ "a063494", "A063494", "a(n) = (2*n - 1)*(7*n^2 - 7*n + 3)/3", "polynomial",
  "Every gap is different, from 16 to 140,000,000,002; every decomposable term is forced level (l <= d^2); 7 terms do not decompose.",
  100000, 1, gen_a063494 },
{ "a047205", "A047205", "Numbers that are congruent to {0, 3, 4} mod 5", "residue class",
  "The gaps are 1 and 3; the level share is 14.85 %; L = 1 holds 77 % of the level class; there are no ties.",
  100000, 1, gen_a047205 },
{ "a047206", "A047206", "Numbers that are congruent to {1, 3, 4} mod 5", "residue class",
  "The gaps are 1 and 2; the level share is 15.48 %; L = 1 holds 74 % of the level class.",
  100000, 1, gen_a047206 },
{ "a047207", "A047207", "Numbers that are congruent to {0, 1, 3, 4} mod 5", "residue class",
  "The gaps are 1 and 2; the level share is 13.30 %; L = 1 holds 88 % of the level class.",
  100000, 1, gen_a047207 },
{ "a047208", "A047208", "Numbers that are congruent to {0, 4} mod 5", "residue class",
  "The gaps are 1 and 4; the level share is 17.57 %; L = 1 holds 63 % of the level class.",
  100000, 1, gen_a047208 },
{ "a047218", "A047218", "Numbers that are congruent to {0, 3} mod 5", "residue class",
  "The gaps are 2 and 3; the level share is 18.92 %; L = 1 holds 58 % of the level class.",
  100000, 1, gen_a047218 },
{ "a047219", "A047219", "Numbers that are congruent to {1, 3} mod 5", "residue class",
  "The gaps are 2 and 3; the level share is 8.46 %; L = 1 holds 65 % of the level class.",
  100000, 1, gen_a047219 },
{ "a047221", "A047221", "Numbers that are congruent to {2, 3} mod 5", "residue class",
  "The gaps are 1 and 4; the level share is 17.55 %; L = 1 holds 63 % of the level class.",
  100000, 1, gen_a047221 },
{ "a047222", "A047222", "Numbers that are congruent to {0, 2, 3} mod 5", "residue class",
  "The gaps are 1 and 2; the level share is 15.50 %; L = 1 holds 74 % of the level class.",
  100000, 1, gen_a047222 },
{ "a047223", "A047223", "Numbers that are congruent to {1, 2, 3} mod 5", "residue class",
  "The gaps are 1 and 3; the level share is 3.81 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a047223 },
{ "a047227", "A047227", "Numbers that are congruent to {1, 2, 3, 4} mod 6", "residue class",
  "The gaps are 1 and 3; the level share is 13.81 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a047227 },
{ "a047228", "A047228", "Numbers that are congruent to {2, 3, 4} mod 6", "residue class",
  "The gaps are 1 and 4; the level share is 8.99 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a047228 },
{ "a047230", "A047230", "Numbers that are congruent to {3, 4} mod 6", "residue class",
  "The gaps are 1 and 5; the level share is 16.04 %; L = 1 holds 81 % of the level class; there are no ties.",
  100000, 1, gen_a047230 },
{ "a047231", "A047231", "Numbers that are congruent to {0, 3, 4} mod 6", "residue class",
  "The gaps are 1, 2 and 3; the level share is 11.43 %; L = 3 holds 58 % of the level class.",
  100000, 1, gen_a047231 },
{ "a047233", "A047233", "Numbers that are congruent to {0, 4} mod 6", "residue class",
  "The gaps are 2 and 4; the level share is 17.53 %; L = 2 holds 79 % of the level class; there are no ties.",
  100000, 1, gen_a047233 },
{ "a047234", "A047234", "Numbers that are congruent to {0, 1, 4} mod 6", "residue class",
  "The gaps are 1, 2 and 3; the level share is 18.59 %; L = 2 holds 52 % of the level class.",
  100000, 1, gen_a047234 },
{ "a047235", "A047235", "Numbers that are congruent to {2, 4} mod 6", "residue class",
  "The gaps are 2 and 4; the level share is 0.01 %; L = 2 holds 40 % of the level class.",
  100000, 1, gen_a047235 },
{ "a047236", "A047236", "Numbers that are congruent to {1, 2, 4} mod 6", "residue class",
  "The gaps are 1, 2 and 3; the level share is 8.99 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a047236 },
{ "a047237", "A047237", "Numbers that are congruent to {0, 1, 2, 4} mod 6", "residue class",
  "The gaps are 1 and 2; the level share is 10.62 %; L = 1 holds 65 % of the level class; there are no ties.",
  100000, 1, gen_a047237 },
{ "a047242", "A047242", "Numbers that are congruent to {0, 1, 3} mod 6", "residue class",
  "The gaps are 1, 2 and 3; the level share is 17.99 %; L = 1 holds 100 % of the level class; there are no ties.",
  100000, 1, gen_a047242 },
{ "a047243", "A047243", "Numbers that are congruent to {2, 3} mod 6", "residue class",
  "The gaps are 1 and 5; the level share is 23.59 %; L = 1 holds 55 % of the level class.",
  100000, 1, gen_a047243 },
{ "a047244", "A047244", "Numbers that are congruent to {0, 2, 3} mod 6", "residue class",
  "The gaps are 1, 2 and 3; the level share is 13.80 %; L = 1 holds 65 % of the level class.",
  100000, 1, gen_a047244 },
{ "a047245", "A047245", "Numbers that are congruent to {1, 2, 3} mod 6", "residue class",
  "The gaps are 1 and 4; the level share is 17.98 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a047245 },
{ "a047249", "A047249", "Numbers that are congruent to {3, 4, 5} mod 6", "residue class",
  "The gaps are 1 and 4; the level share is 8.99 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a047249 },
{ "a047252", "A047252", "Numbers that are congruent to {0, 1, 3, 4, 5} mod 6", "residue class",
  "The gaps are 1 and 2; the level share is 11.33 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a047252 },
{ "a047253", "A047253", "Numbers that are congruent to {1, 2, 3, 4, 5} mod 6", "residue class",
  "The gaps are 1 and 2; the level share is 5.64 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a047253 },
{ "a045324", "A045324", "Primes congruent to {0, 1, 2, 3, 4} (mod 7)", "primes",
  "79 different gaps occur, from 1 to 190; the level share is 26.36 %; L = 1 holds 31 % of the level class.",
  100000, 1, gen_a045324 },
{ "a045325", "A045325", "Primes congruent to {0, 2, 3, 4} (mod 7)", "primes",
  "83 different gaps occur, from 1 to 258; the level share is 28.94 %.",
  100000, 1, gen_a045325 },
{ "a045328", "A045328", "Primes congruent to {0, 1, 2, 3, 6} (mod 7)", "primes",
  "80 different gaps occur, from 1 to 204; the level share is 27.36 %; L = 1 holds 32 % of the level class.",
  100000, 1, gen_a045328 },
{ "a045329", "A045329", "Primes congruent to {0, 2, 3, 6} (mod 7)", "primes",
  "80 different gaps occur, from 1 to 252; the level share is 29.35 %.",
  100000, 1, gen_a045329 },
{ "a045346", "A045346", "Primes congruent to {0, 1, 2, 4, 5, 6} mod 7", "primes",
  "67 different gaps occur, from 2 to 158; the level share is 24.38 %; L = 3 holds 30 % of the level class.",
  100000, 1, gen_a045346 },
{ "a045347", "A045347", "Primes congruent to {0, 2, 4, 5, 6} mod 7", "primes",
  "84 different gaps occur, from 2 to 178; the level share is 26.00 %.",
  100000, 1, gen_a045347 },
{ "a045350", "A045350", "Primes congruent to {0, 1, 2, 4, 5} mod 7", "primes",
  "76 different gaps occur, from 2 to 210; the level share is 26.57 %.",
  100000, 1, gen_a045350 },
{ "a045351", "A045351", "Primes congruent to {0, 2, 4, 5} mod 7", "primes",
  "103 different gaps occur, from 2 to 260; the level share is 29.45 %.",
  100000, 1, gen_a045351 },
{ "a045352", "A045352", "Primes congruent to {1, 2, 5, 7} mod 8", "primes",
  "74 different gaps occur, from 2 to 162; the level share is 25.04 %.",
  100000, 1, gen_a045352 },
{ "a045353", "A045353", "Primes congruent to {0, 1, 2, 5, 6} mod 7", "primes",
  "78 different gaps occur, from 2 to 194; the level share is 26.45 %; L = 3 holds 30 % of the level class.",
  100000, 1, gen_a045353 },
{ "a045354", "A045354", "Primes congruent to {0, 2, 5, 6} mod 7", "primes",
  "81 different gaps occur, from 2 to 276; the level share is 27.84 %; L = 3 holds 31 % of the level class.",
  100000, 1, gen_a045354 },
{ "a045358", "A045358", "Primes congruent to {0, 1, 2, 5} mod 7", "primes",
  "81 different gaps occur, from 2 to 280; the level share is 30.11 %.",
  100000, 1, gen_a045358 },
{ "a045369", "A045369", "Primes congruent to {0, 1, 2, 4, 6} mod 7", "primes",
  "80 different gaps occur, from 2 to 192; the level share is 26.30 %.",
  100000, 1, gen_a045369 },
{ "a045370", "A045370", "Primes congruent to {0, 2, 4, 6} mod 7", "primes",
  "79 different gaps occur, from 2 to 306; the level share is 29.43 %.",
  100000, 1, gen_a045370 },
{ "a045376", "A045376", "Primes congruent to {0, 1, 2, 6} mod 7", "primes",
  "108 different gaps occur, from 2 to 254; the level share is 29.60 %.",
  100000, 1, gen_a045376 },
{ "a045393", "A045393", "Primes congruent to {0, 1, 3, 4, 5, 6} mod 7", "primes",
  "65 different gaps occur, from 2 to 168; the level share is 23.73 %; L = 1 holds 30 % of the level class.",
  100000, 1, gen_a045393 },
{ "a045394", "A045394", "Primes congruent to {0, 3, 4, 5, 6} mod 7", "primes",
  "81 different gaps occur, from 2 to 224; the level share is 24.76 %.",
  100000, 1, gen_a045394 },
{ "a045396", "A045396", "Primes congruent to {0, 1, 3, 4, 5} mod 7", "primes",
  "83 different gaps occur, from 2 to 184; the level share is 27.47 %.",
  100000, 1, gen_a045396 },
{ "a045397", "A045397", "Primes congruent to {0, 3, 4, 5} mod 7", "primes",
  "78 different gaps occur, from 2 to 252; the level share is 29.80 %.",
  100000, 1, gen_a045397 },
{ "a045398", "A045398", "Primes congruent to {0, 1, 3, 5, 6} mod 7", "primes",
  "79 different gaps occur, from 2 to 192; the level share is 24.66 %.",
  100000, 1, gen_a045398 },
{ "a047254", "A047254", "Numbers that are congruent to {2, 3, 5} mod 6", "residue class",
  "The gaps are 1, 2 and 3; the level share is 22.76 %; L = 1 holds 79 % of the level class.",
  100000, 1, gen_a047254 },
{ "a047256", "A047256", "Numbers that are congruent to {0, 1, 2, 3, 5} mod 6", "residue class",
  "The gaps are 1 and 2; the level share is 16.93 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a047256 },
{ "a047257", "A047257", "Numbers that are congruent to {4, 5} mod 6", "residue class",
  "The gaps are 1 and 5; the level share is 0.01 %; L = 1 holds 33 % of the level class.",
  100000, 1, gen_a047257 },
{ "a047258", "A047258", "Numbers that are congruent to {0, 4, 5} mod 6", "residue class",
  "The gaps are 1 and 4; the level share is 7.36 %; L = 2 holds 65 % of the level class.",
  100000, 1, gen_a047258 },
{ "a047259", "A047259", "Numbers that are congruent to {1, 4, 5} mod 6", "residue class",
  "The gaps are 1, 2 and 3; the level share is 4.81 %; L = 2 holds 100 % of the level class.",
  100000, 1, gen_a047259 },
{ "a047260", "A047260", "Numbers that are congruent to {0, 1, 4, 5} mod 6", "residue class",
  "The gaps are 1 and 3; the level share is 10.65 %; L = 1 holds 65 % of the level class.",
  100000, 1, gen_a047260 },
{ "a047262", "A047262", "Numbers that are congruent to {0, 2, 4, 5} mod 6", "residue class",
  "The gaps are 1 and 2; the level share is 3.71 %; L = 2 holds 100 % of the level class.",
  100000, 1, gen_a047262 },
{ "a047263", "A047263", "Numbers that are congruent to {0, 1, 2, 4, 5} mod 6", "residue class",
  "The gaps are 1 and 2; the level share is 5.67 %; L = 1 holds 100 % of the level class.",
  100000, 1, gen_a047263 },
};
static const int NDEF = (int)(sizeof defs / sizeof defs[0]);

/* ------------------------------------------------------------------ */

static void mkdirp(const char *p)
{
    char buf[1024]; snprintf(buf, sizeof buf, "%s", p);
    for (char *q = buf + 1; *q; q++)
        if (*q == '/') { *q = 0; mkdir(buf, 0755); *q = '/'; }
    mkdir(buf, 0755);
}

static char *u64s(char *p, u64 v)
{
    char tmp[24]; int n = 0;
    if (!v) tmp[n++] = '0';
    while (v) { tmp[n++] = '0' + (char)(v % 10); v /= 10; }
    while (n) *p++ = tmp[--n];
    return p;
}

typedef struct {
    long terms, decomposable, level, weight, ties, lone, forced;
    u64 amin, amax, kmax, Lmax, dmin, dmax;
    int chunks;
} Stats;

static Stats run_one(const char *outdir, SeqDef *S, long terms)
{
    Stats st; memset(&st, 0, sizeof st);
    st.terms = terms; st.dmin = ~0ULL; st.amin = ~0ULL;

    u64 *t = malloc((size_t)(terms + 1) * sizeof(u64));
    if (!t) { fprintf(stderr, "term alloc failed\n"); exit(1); }
    clock_t c0 = clock();
    S->gen(t, terms + 1);
    for (long i = 0; i < terms; i++)
        if (t[i + 1] <= t[i]) { fprintf(stderr, "%s: not strictly increasing at %ld\n", S->id, i); exit(1); }
    /* the browser holds every value in a double: exact only below 2^53 */
    if (t[terms] >= (1ULL << 53)) { fprintf(stderr, "%s: a(n) passes 2^53\n", S->id); exit(1); }

    char dir[1024]; snprintf(dir, sizeof dir, "%s/seq/%s", outdir, S->id);
    mkdirp(dir);

    char *buf = malloc(CHUNK_ROWS * 64 + 64);
    long i = 0; int chunk = 0;
    while (i < terms) {
        long end = i + CHUNK_ROWS; if (end > terms) end = terms;
        char *p = buf;
        memcpy(p, "a,d,k,L\n", 8); p += 8;
        for (long j = i; j < end; j++) {
            u64 k, L, d;
            int ok = decomp(t[j], t[j + 1], &k, &L, &d);
            p = u64s(p, t[j]); *p++ = ',';
            p = u64s(p, d);    *p++ = ',';
            p = u64s(p, k);    *p++ = ',';
            p = u64s(p, L);    *p++ = '\n';

            if (t[j] < st.amin) st.amin = t[j];
            if (t[j] > st.amax) st.amax = t[j];
            if (d < st.dmin) st.dmin = d;
            if (d > st.dmax) st.dmax = d;
            if (ok) {
                if (k * L + d != t[j]) { fprintf(stderr, "%s: identity failed at a=%llu\n", S->id, (unsigned long long)t[j]); exit(1); }
                st.decomposable++;
                if (k > L) st.level++; else st.weight++;
                /* forced level: l <= d^2, i.e. alpha = log d / log l >= 1/2, so the
                   divisor window (d, sqrt l] is empty and the term cannot be weight */
                if ((u128)(t[j] - d) <= (u128)d * d) st.forced++;   /* 128-bit: d^2 overflows 64 bits once d > 2^32 */
                if (k == L) st.ties++;
                if (L == 1) st.lone++;
                if (k > st.kmax) st.kmax = k;
                if (L > st.Lmax) st.Lmax = L;
            }
        }
        char path[1200]; snprintf(path, sizeof path, "%s/chunk-%03d.csv", dir, chunk);
        FILE *f = fopen(path, "wb");
        if (!f) { perror(path); exit(1); }
        fwrite(buf, 1, (size_t)(p - buf), f);
        fclose(f);
        chunk++; i = end;
    }
    st.chunks = chunk;
    free(buf); free(t);
    fprintf(stderr, "  %-12s %8ld terms  %2d chunks  a<=%llu  level %.2f%%  %.1fs\n",
            S->id, terms, chunk, (unsigned long long)st.amax,
            100.0 * st.level / (st.decomposable ? st.decomposable : 1),
            (double)(clock() - c0) / CLOCKS_PER_SEC);
    return st;
}

static void csv_quote(FILE *f, const char *s)
{
    fputc('"', f);
    for (const char *p = s; *p; p++) { if (*p == '"') fputc('"', f); fputc(*p, f); }
    fputc('"', f);
}

int main(int argc, char **argv)
{
    const char *outdir = argc > 1 ? argv[1] : "data";
    const char *only   = argc > 2 ? argv[2] : NULL;
    long override      = argc > 3 ? atol(argv[3]) : 0;

    fprintf(stderr, "sieving...\n");
    sieve_init();

    mkdirp(outdir);
    char cpath[1200]; snprintf(cpath, sizeof cpath, "%s/catalog.csv", outdir);

    /* regenerating one sequence: carry every other catalogue row across unchanged */
    char *old = NULL;
    if (only) {
        FILE *f = fopen(cpath, "rb");
        if (f) { fseek(f, 0, SEEK_END); long n = ftell(f); fseek(f, 0, SEEK_SET);
                 old = malloc(n + 1);
                 if (fread(old, 1, n, f) != (size_t)n) { free(old); old = NULL; } else old[n] = 0;
                 fclose(f); }
    }

    FILE *cat = fopen(cpath, "wb");
    if (!cat) { perror(cpath); exit(1); }
    fprintf(cat, "id,anumber,name,family,n0,terms,chunks,chunk_rows,decomposable,"
                 "level,weight,ties,level_one,forced,amin,amax,kmax,Lmax,dmin,dmax,note\n");
    if (old) {
        size_t idlen = strlen(only);
        char *line = strchr(old, '\n');            /* drop the old header */
        while (line && *++line) {
            char *nl = strchr(line, '\n');
            size_t len = nl ? (size_t)(nl - line) : strlen(line);
            if (len && !(strncmp(line, only, idlen) == 0 && line[idlen] == ','))
                fprintf(cat, "%.*s\n", (int)len, line);
            if (!nl) break;
            line = nl;
        }
        free(old);
    }

    for (int i = 0; i < NDEF; i++) {
        if (only && strcmp(only, defs[i].id)) continue;
        long terms = override ? override : defs[i].terms;
        Stats st = run_one(outdir, &defs[i], terms);
        fprintf(cat, "%s,%s,", defs[i].id, defs[i].anum);
        csv_quote(cat, defs[i].name); fputc(',', cat);
        csv_quote(cat, defs[i].family); fputc(',', cat);
        fprintf(cat, "%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,%llu,%llu,%llu,%llu,%llu,%llu,",
                defs[i].n0, st.terms, st.chunks, CHUNK_ROWS, st.decomposable,
                st.level, st.weight, st.ties, st.lone, st.forced,
                (unsigned long long)st.amin, (unsigned long long)st.amax,
                (unsigned long long)st.kmax, (unsigned long long)st.Lmax,
                (unsigned long long)st.dmin, (unsigned long long)st.dmax);
        csv_quote(cat, defs[i].note); fputc('\n', cat);
    }
    fclose(cat);
    fprintf(stderr, "done -> %s\n", outdir);
    return 0;
}
