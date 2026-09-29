/* data.js — loading one sequence.
 *
 * A sequence lives in data/seq/<id>/chunk-NNN.csv, each chunk `chunk_rows`
 * rows of "a,d,k,L", every term present whether it decomposes or not, so
 * row i of the concatenation is always sequence index n = n0 + i.
 *
 * Chunks are fetched with bounded concurrency and parsed straight into the
 * right slice of pre-allocated columns, so they may land in any order.  The
 * contiguous prefix is converted to point attributes as it completes and
 * handed to the renderer, which is why the cloud builds in front of you
 * instead of appearing all at once at the end.
 */

import { parseIntRows } from './csv.js';

const LOG10 = Math.LN10;
const CONCURRENCY = 4;

export const FIELDS = {
  k: { label: 'log₁₀ k  (weight)', short: 'k' },
  L: { label: 'log₁₀ L  (level)',  short: 'L' },
  d: { label: 'log₁₀ d  (jump)',   short: 'd' },
  a: { label: 'log₁₀ a  (term)',   short: 'a' },
  l: { label: 'log₁₀ ℓ  (k·L)',    short: 'l' },
  n: { label: 'log₁₀ n  (index)',  short: 'n' },
};

export class Dataset {
  constructor(rec, baseUrl) {
    this.rec = rec;
    this.base = `${baseUrl}/seq/${rec.id}`;
    const N = rec.terms;

    this.a = new Float64Array(N);
    this.d = new Float64Array(N);
    this.k = new Float64Array(N);
    this.L = new Float64Array(N);

    this.count = 0;                       // decomposable points converted
    this.row   = new Int32Array(N);       // point -> row (n - n0)
    this.cls   = new Uint8Array(N);       // 0 weight, 1 level, 2 tie
    this.one   = new Uint8Array(N);       // 1 when L = 1

    this._cache = new Map();              // field -> Float32Array
    this._built = 0;                      // rows already converted
    this._done = new Uint8Array(rec.chunks);
    this._prefix = 0;
    this.complete = false;
  }

  /** @param {(d: Dataset, frac: number) => void} onGrow */
  async load(onGrow, signal) {
    const { chunks, chunk_rows, terms } = this.rec;
    let next = 0, loaded = 0;

    const worker = async () => {
      for (;;) {
        const c = next++;
        if (c >= chunks) return;
        const url = `${this.base}/chunk-${String(c).padStart(3, '0')}.csv`;
        const res = await fetch(url, { signal });
        if (!res.ok) throw new Error(`${res.status} — ${url}`);
        const bytes = new Uint8Array(await res.arrayBuffer());
        const at = c * chunk_rows;
        const want = Math.min(chunk_rows, terms - at);
        const got = parseIntRows(bytes, [this.a, this.d, this.k, this.L], at);
        if (got < want)
          throw new Error(`${url}: ${got} rows, expected ${want}`);
        this._done[c] = 1;
        loaded++;
        while (this._prefix < chunks && this._done[this._prefix]) this._prefix++;
        this._convert(Math.min(this._prefix * chunk_rows, terms));
        onGrow(this, loaded / chunks);
      }
    };

    await Promise.all(Array.from({ length: Math.min(CONCURRENCY, chunks) }, worker));
    this._convert(terms);
    this.complete = true;
    onGrow(this, 1);
  }

  /** Convert rows [_built, upto) into point attributes. */
  _convert(upto) {
    if (upto <= this._built) return;
    const { k, L } = this;
    let c = this.count;
    for (let i = this._built; i < upto; i++) {
      const kk = k[i];
      if (kk === 0) continue;                     // non-decomposable: no point
      const ll = L[i];
      this.row[c] = i;
      this.cls[c] = kk === ll ? 2 : (kk > ll ? 1 : 0);
      this.one[c] = ll === 1 ? 1 : 0;
      c++;
    }
    this.count = c;
    this._built = upto;
    this._cache.clear();                          // lengths changed
  }

  /** log₁₀ of one field over the decomposable points, cached. */
  field(name) {
    const hit = this._cache.get(name);
    if (hit && hit.length >= this.count) return hit;
    const n = this.count, out = new Float32Array(n);
    const row = this.row, n0 = this.rec.n0;
    if (name === 'n') {
      for (let i = 0; i < n; i++) out[i] = Math.log(row[i] + n0) / LOG10;
    } else if (name === 'l') {
      const k = this.k, L = this.L;
      for (let i = 0; i < n; i++) { const r = row[i]; out[i] = Math.log(k[r] * L[r]) / LOG10; }
    } else {
      const src = this[name];
      for (let i = 0; i < n; i++) out[i] = Math.log(src[row[i]]) / LOG10;
    }
    this._cache.set(name, out);
    return out;
  }

  extent(name) {
    const v = this.field(name);
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < v.length; i++) { const x = v[i]; if (x < lo) lo = x; if (x > hi) hi = x; }
    if (!isFinite(lo)) { lo = 0; hi = 1; }
    if (hi - lo < 1e-9) { lo -= .5; hi += .5; }
    return [lo, hi];
  }

  /** Counts under the live filters — recomputed on the CPU, ~1 ms at 2·10⁵. */
  stats(dLo, dHi, nMax) {
    const logd = this.field('d'), cls = this.cls, one = this.one;
    let weight = 0, level = 0, tie = 0, lone = 0;
    const n = Math.min(this.count, nMax);
    for (let i = 0; i < n; i++) {
      const x = logd[i];
      if (x < dLo || x > dHi) continue;
      const c = cls[i];
      if (c === 1) level++; else { weight++; if (c === 2) tie++; }
      if (one[i]) lone++;
    }
    return { weight, level, tie, one: lone, shown: weight + level };
  }

  /** Everything about one point, for the read-out. */
  point(i) {
    const r = this.row[i];
    const k = this.k[r], L = this.L[r], d = this.d[r], a = this.a[r];
    return { i, n: r + this.rec.n0, a, d, k, L, l: k * L, cls: this.cls[i], one: !!this.one[i] };
  }
}
