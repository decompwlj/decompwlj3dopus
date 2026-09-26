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
