/* worker.js — the Lab's arithmetic, off the main thread.
 *
 *   control    the record-level residue control: for each record keep the jump d and
 *              the class of ℓ mod M, redraw ℓ uniformly among the members of that class
 *              in [max(d + 1, ⌈0.75 ℓ⌉), ⌊1.25 ℓ⌋], and decompose the new ℓ' with the
 *              same kernel.  Returns ℓ' and k' for every record; L' = ℓ'/k'.
 *   decompose  a pasted sequence: d, k, L for every term but the last.
 */
importScripts('kernel.js');

const TWO53 = 2 ** 53;
const FROZEN = 16;            /* a class with fewer members than this in the window is "nearly frozen" */

self.onmessage = (e) => {
  const m = e.data;
  try {
    if (m.type === 'control') control(m);
    else if (m.type === 'decompose') decompose(m);
  } catch (err) {
    self.postMessage({ type: 'error', task: m.task, message: String(err && err.message || err) });
  }
};

function control({ task, l, d, M, seed }) {
  const n = l.length, R = DWLJ.rng(seed);
  const lp = new Float64Array(n), k = new Float64Array(n);
  let frozen = 0;
  for (let i = 0; i < n; i++) {
    const li = l[i], di = d[i];
    const lo = Math.max(di + 1, Math.ceil(0.75 * li));
    const hi = Math.min(TWO53 - 1, Math.floor(1.25 * li));
    const c = li % M;
    const first = lo + ((c - (lo % M)) % M + M) % M;        /* least member of the class ≥ lo */
    const count = Math.floor((hi - first) / M) + 1;           /* ≥ 1: ℓ itself is a member */
    if (count < FROZEN) frozen++;
    const x = first + M * R.below(count);
    lp[i] = x;
    k[i] = DWLJ.leastDivisorAbove(x, di);
    if ((i & 2047) === 2047) self.postMessage({ type: 'progress', task, done: i + 1 });
  }
  self.postMessage({ type: 'control', task, lp, k, frozen }, [lp.buffer, k.buffer]);
}

function decompose({ task, a }) {
  const n = a.length - 1;
  const d = new Float64Array(n), k = new Float64Array(n), L = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const di = a[i + 1] - a[i];
    d[i] = di;
    const r = DWLJ.decompose(a[i], di);
    if (r) { k[i] = r[0]; L[i] = r[1]; }
    if ((i & 2047) === 2047) self.postMessage({ type: 'progress', task, done: i + 1 });
  }
  self.postMessage({ type: 'decompose', task, d, k, L }, [d.buffer, k.buffer, L.buffer]);
}
