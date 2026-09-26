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

/* smooth numbers: all of them below 2^53, sorted */
static size_t smooth_m; static u64 *smooth_v;
static void smooth_rec(const u64 *ps, int np, int i, u64 v, u64 lim)
{
    if (i == np) { smooth_v[smooth_m++] = v; return; }
    for (u64 x = v; ; x *= ps[i]) { smooth_rec(ps, np, i + 1, x, lim); if (x > lim / ps[i]) break; }
}
static long smooth_all(u64 *t, long cnt, const u64 *ps, int np)
{
    smooth_v = malloc(4000000 * sizeof(u64)); smooth_m = 0;
    smooth_rec(ps, np, 0, 1, (1ULL << 53) - 1);
    qsort(smooth_v, smooth_m, sizeof(u64), cmp_u64);
    long k = (long)smooth_m < cnt ? (long)smooth_m : cnt;
    memcpy(t, smooth_v, (size_t)k * sizeof(u64));
    free(smooth_v);
    return k;
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
static void gen_res(u64 *t, long cnt, u64 m, u32 mask, u64 from)
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
{ "primepowers", "A000961", "Prime powers", "prime powers",
  "Prime powers, with 1: the primes plus a sparse set of higher powers. The plane is the primes' plane almost exactly (22.97 % level against 23.00 %).",
  100000, 1, gen_primepowers },
{ "primesq", "A001248", "Squares of primes", "prime powers",
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
{ "sorting", "A001855", "Sorting numbers (binary insertion)", "logarithmic gap",
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
{ "abundant", "A005101", "Abundant numbers", "divisor sum",
  "sigma(n) > 2n. 99.2 % of these terms are even (the first odd one is 945), and so is l on 99.2 % of decomposable terms: the level class sits on even L, and the line L = 1 holds only 32 terms.",
  100000, 1, gen_abundant },
{ "squarefree", "A005117", "Squarefree numbers", "multiplicative",
  "Squarefree numbers, density 6/pi^2. Gaps of 1 to 7; the plane is close to the naturals' (10.72 % level) and L = 1 carries almost all of the level class.",
  100000, 1, gen_squarefree },
{ "practical", "A005153", "Practical numbers", "divisor sum",
  "Practical numbers. Every term > 2 is divisible by 4 or 6, so a and d are even and 2 | l at every decomposable term: the line L = 1 holds only 2 terms and the level class sits on L = 4, 6, 8, ... The level share is 17.50 %.",
  100000, 1, gen_practical },
{ "trisq", "A005214", "Triangular numbers and squares", "polynomial union",
  "The union of two sequences that are 100 % level-classified is not: merging triangular numbers and squares shrinks the gaps below sqrt(l) on 58.6 % of terms, and 16.2 % turn weight. No ties occur.",
  100000, 1, gen_trisq },
{ "harshad", "A005349", "Harshad numbers", "digit rule",
  "n divisible by its digit sum. Density by digit sum is 9:3:1 according to gcd(s,9), but that structure does not reach the plane: the weight fence is flat and the level share varies only between 11 and 24 %.",
  100000, 1, gen_harshad },
{ "odd", "A005408", "Odd numbers", "forced divisor",
  "d = 2 and l = a - 2 is odd, so the divisor 2 is never available and k = spf(a - 2). The level class is exactly {a : a - 2 prime}, all on L = 1: 17,982 = pi(199,997) - 1 terms here. The ties are a - 2 = p^2, 85 = pi(447) - 1. Removing the divisor 2 nearly doubles the level share against the naturals: 17.98 % against 9.59 %.",
  100000, 0, gen_odd },
{ "even", "A005843", "Even numbers", "parity",
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
{ "deficient", "A005100", "Deficient numbers", "divisor sum",
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
{ "arith", "A003601", "Arithmetic numbers", "divisor sum",
  "The mean of the divisors, sigma(n)/d(n), is an integer. The level share is 11.93 %; L = 1 holds 97 % of the level class.",
  100000, 1, gen_arith },
{ "totients", "A002202", "Totients (values of phi)", "totient",
  "Values of Euler's phi. Past 1 every totient is even, so 2 | l and the line L = 1 holds 2 terms. The level share is 12.58 %; L = 4 holds 40 % of the level class.",
  100000, 1, gen_totients },
{ "nontotients", "A007617", "Nontotients", "totient",
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
{ "cyclic", "A003277", "Cyclic numbers: k such that k and phi(k) are relatively prime; also k such that there is just one group of order k, i.e., A000001(k) = 1", "totient",
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
{ "tau4", "A030513", "Numbers with 4 divisors", "divisor count",
  "38 different gaps occur, from 1 to 47; the level share is 15.72 %; L = 1 holds 61 % of the level class.",
  100000, 1, gen_tau4 },
{ "tau6", "A030515", "Numbers with exactly 6 divisors", "divisor count",
  "189 different gaps occur, from 1 to 240; the level share is 27.45 %.",
  100000, 1, gen_tau6 },
{ "tau8", "A030626", "Numbers with exactly 8 divisors", "divisor count",
  "37 different gaps occur, from 1 to 44; the level share is 16.53 %; L = 1 holds 62 % of the level class.",
  100000, 1, gen_tau8 },
{ "tau10", "A030628", "1 together with numbers of the form p*q^4 and p^9, where p and q are distinct primes", "divisor count",
  "661 different gaps occur, from 1 to 1,600; the level share is 29.84 %.",
  100000, 1, gen_tau10 },
{ "tau12", "A030630", "Numbers with 12 divisors", "divisor count",
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
{ "evennontot", "A005277", "Nontotients: even numbers k such that phi(m) = k has no solution", "totient",
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
{ "sigmaeven", "A028983", "Numbers whose sum of divisors is even", "divisor sum",
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
{ "zeckodd", "A020899", "Numbers k with an odd number of terms in their Zeckendorf representation (write k as a sum of non-consecutive distinct Fibonacci numbers)", "Fibonacci",
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
