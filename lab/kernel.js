/* kernel.js — the decomposition, exact for every integer below 2^53.
 *
 *   a(n) = k(n) · L(n) + d(n)       d = a(n+1) − a(n),  ℓ = a − d,
 *                                    k = the least divisor of ℓ above d,  L = ℓ / k
 *
 * The same function as decomp_fact in tools/decompwlj.gp: factor ℓ, walk its
 * divisors, keep the least one above d.  Factoring is trial division by the
 * primes below 2^16, then a deterministic Miller–Rabin test (the seven bases of
 * Jim Sinclair, exact below 2^64) and Brent's variant of Pollard's rho.
 *
 * Arithmetic stays in doubles and stays exact: a product of two residues below
 * m < 2^51 is formed as an unevaluated sum (Dekker's two-product) and reduced
 * with a corrected quotient; at m ≥ 2^51 the code switches to BigInt.  Nothing
 * here rounds: every value is an integer below 2^53 at every step.
 *
 * A classic script (no import/export): the Lab page loads it with <script>, its
 * worker with importScripts(), and the Node tests with vm.  It defines one
 * global, DWLJ.
 */
(function (G) {
  'use strict';

  const TWO26 = 67108864, TWO51 = 2 ** 51, TWO53 = 2 ** 53;

  /* primes below 2^16 */
  const SMALL = (() => {
    const N = 65536, s = new Uint8Array(N), out = [];
    for (let i = 2; i < N; i++) {
      if (s[i]) continue;
      out.push(i);
      for (let j = i * i; j < N; j += i) s[j] = 1;
    }
    return Int32Array.from(out);
  })();
  const TRIAL = 1024;
  let TRIAL_COUNT = 0;
  while (SMALL[TRIAL_COUNT] < TRIAL) TRIAL_COUNT++;

  /* Dekker: x·y = hi + lo exactly, for |x|, |y| < 2^53 */
  const SPLIT = 134217729;              /* 2^27 + 1 */
  function twoProd(x, y) {
    const p = x * y;
    let t = SPLIT * x; const xh = t - (t - x), xl = x - xh;
    t = SPLIT * y;     const yh = t - (t - y), yl = y - yh;
    const e = ((xh * yh - p) + xh * yl + xl * yh) + xl * yl;
    return [p, e];
  }

  /** a·b mod m, for 0 ≤ a, b < m < 2^53 */
  function mulmod(a, b, m) {
    if (m <= TWO26) return (a * b) % m;
    if (m < TWO51) {
      const q = Math.floor((a * b) / m);
      const [ph, pl] = twoProd(a, b), [qh, ql] = twoProd(q, m);
      let r = (ph - qh) + (pl - ql);       /* exact: |r| < 3m < 2^53 */
      while (r < 0) r += m;
      while (r >= m) r -= m;
      return r;
    }
    return Number((BigInt(a) * BigInt(b)) % BigInt(m));
  }

  function powmod(b, e, m) {
    let r = 1 % m; b %= m;
    while (e > 0) {
      if (e % 2 === 1) r = mulmod(r, b, m);
      e = Math.floor(e / 2);
      if (e > 0) b = mulmod(b, b, m);
    }
    return r;
  }

  const MR_BASES = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];

  /** deterministic below 2^64, hence for every n < 2^53 */
  function isPrime(n) {
    if (n < 2) return false;
    for (let i = 0; i < 12; i++) {
      const p = SMALL[i];
      if (n === p) return true;
      if (n % p === 0) return false;
    }
    if (n < 1369) return true;            /* 37² */
    let dd = n - 1, s = 0;
    while (dd % 2 === 0) { dd /= 2; s++; }
    outer:
    for (const B of MR_BASES) {
      const a = B % n;
      if (a === 0) continue;
      let x = powmod(a, dd, n);
      if (x === 1 || x === n - 1) continue;
      for (let r = 1; r < s; r++) {
        x = mulmod(x, x, n);
        if (x === n - 1) continue outer;
      }
      return false;
    }
    return true;
  }

  function gcd(a, b) { while (b) { const t = a % b; a = b; b = t; } return a; }

  /** a nontrivial factor of the odd composite n (Brent's rho) */
  function rho(n) {
    for (let c = 1; ; c++) {
      const f = (x) => { const y = mulmod(x, x, n) + c; return y >= n ? y - n : y; };
      let y = 2, x = 2, g = 1, q = 1, ys = 2, r = 1;
      const M = 64;
      while (g === 1) {
        x = y;
        for (let i = 0; i < r; i++) y = f(y);
        for (let k = 0; k < r && g === 1; k += M) {
          ys = y;
          for (let i = 0; i < Math.min(M, r - k); i++) {
            y = f(y);
            q = mulmod(q, x > y ? x - y : y - x, n);
          }
          g = gcd(q, n);
        }
        r *= 2;
      }
      if (g === n) {
        do { ys = f(ys); g = gcd(x > ys ? x - ys : ys - x, n); } while (g === 1);
      }
      if (g !== n) return g;
    }
  }

  /** prime factorisation of 1 ≤ n < 2^53 as [[p, e], …], p ascending */
  function factor(n) {
    if (!(Number.isInteger(n) && n >= 1 && n < TWO53)) throw new RangeError(`factor: ${n}`);
    const out = [];
    /* trial division: all the way to √n below 2^32, else by the primes below TRIAL */
    const lim = n < 16777216 ? SMALL.length : TRIAL_COUNT;
    let i = 0;
    for (; i < lim; i++) {
      const p = SMALL[i];
      if (p * p > n) break;
      if (n % p === 0) {
        let e = 0;
        do { n /= p; e++; } while (n % p === 0);
        out.push([p, e]);
      }
    }
    if (n > 1) {
      const t = SMALL[i] || 65537;           /* every prime factor left is ≥ t */
      if (t * t > n) out.push([n, 1]);       /* so n itself is prime */
      else {
        const big = [];
        const split = (m) => {
          if (isPrime(m)) { big.push(m); return; }
          const s = Math.sqrt(m);              /* rho cannot split p²: catch squares first */
          const r = Math.round(s);
          if (r * r === m) { split(r); split(r); return; }
          const g = rho(m);
          split(g); split(m / g);
        };
        split(n);
        big.sort((x, y) => x - y);
        for (const p of big) {
          const last = out[out.length - 1];
          if (last && last[0] === p) last[1]++; else out.push([p, 1]);
        }
      }
    }
    return out;
  }

  /** the least divisor of l strictly above d (l > d ≥ 0) */
  function leastDivisorAbove(l, d) {
    if (d < 1) return 1;
    if (l <= d) throw new RangeError('l ≤ d');
    const F = factor(l);
    let best = l;                         /* l itself always qualifies */
    /* depth-first over the divisor lattice, pruning products already ≥ best */
    const walk = (i, v) => {
      if (v > d) { if (v < best) best = v; return; }
      if (i === F.length) return;
      const [p, e] = F[i];
      let w = v;
      for (let j = 0; j <= e; j++) {
        if (w >= best) break;
        walk(i + 1, w);
        w *= p;
      }
    };
    walk(0, 1);
    return best;
  }

  /** decompose a with jump d: [k, L] or null when a ≤ 2d (not decomposable) */
  function decompose(a, d) {
    if (a <= 2 * d) return null;
    const l = a - d, k = leastDivisorAbove(l, d);
    return [k, l / k];
  }

  /** class code: 0 not decomposable, 1 weight (k < L), 2 tie (k = L), 3 level (k > L) */
  function classify(k, L) { return k === 0 ? 0 : k > L ? 3 : k === L ? 2 : 1; }

  /* seeded pseudo-random numbers: xoshiro128** on 32-bit words, seeded by splitmix32 */
  function rng(seed) {
    let s = seed >>> 0;
    const sm = () => { s = (s + 0x9e3779b9) >>> 0; let z = s; z = Math.imul(z ^ (z >>> 16), 0x85ebca6b); z = Math.imul(z ^ (z >>> 13), 0xc2b2ae35); return (z ^ (z >>> 16)) >>> 0; };
    let a = sm(), b = sm(), c = sm(), d = sm();
    const next = () => {
      const r = Math.imul(((Math.imul(b, 5) << 7) | (Math.imul(b, 5) >>> 25)), 9) >>> 0;
      const t = b << 9;
      c ^= a; d ^= b; b ^= c; a ^= d; c ^= t; d = (d << 11) | (d >>> 21);
      return r;
    };
    /* uniform integer in [0, n), n ≤ 2^53, by rejection on 53 random bits */
    const below = (n) => {
      if (n <= 1) return 0;
      const lim = TWO53 - (TWO53 % n);
      for (;;) {
        const v = (next() >>> 11) * 4294967296 + next();   /* 21 + 32 = 53 bits */
        if (v < lim) return v % n;
      }
    };
    return { next, below };
  }

  /** a 32-bit seed from a string (FNV-1a) */
  function seedOf(str) {
    let h = 0x811c9dc5;
    for (let i = 0; i < str.length; i++) { h ^= str.charCodeAt(i); h = Math.imul(h, 0x01000193); }
    return h >>> 0;
  }

  G.DWLJ = { factor, isPrime, mulmod, powmod, leastDivisorAbove, decompose, classify, rng, seedOf, SMALL };
})(typeof self !== 'undefined' ? self : globalThis);
