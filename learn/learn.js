/* decompwlj 3D, written by tools/seo.py: the live example of the "How it works" page.
   The same rule as the generator (tools/decompwlj_gen.c):
     d = a(n+1) - a(n);  if a > 2d:  l = a - d,  k = least divisor of l above d,  L = l / k. */
(() => {
  const $ = (s) => document.querySelector(s);
  const PRESETS = {
    primes: () => { const p = []; for (let n = 2; p.length < 40; n++) if (isPrime(n)) p.push(n); return p; },
    squares: () => Array.from({ length: 40 }, (_, i) => (i + 1) ** 2),
    triangular: () => Array.from({ length: 40 }, (_, i) => (i + 1) * (i + 2) / 2),
    odd: () => Array.from({ length: 40 }, (_, i) => 2 * i + 1),
    fibonacci: () => { const f = [1, 2]; while (f.length < 40) f.push(f[f.length - 1] + f[f.length - 2]); return f; },
  };
  function isPrime(n) { if (n < 2) return false; for (let i = 2; i * i <= n; i++) if (n % i === 0) return false; return true; }
  /** k, L for a term a with jump d, or null when a <= 2d */
  function decomp(a, d) {
    if (a <= 2 * d) return null;
    const l = a - d; let k = l;                      /* l itself always exceeds d */
    for (let i = 1; i * i <= l; i++) {
      if (l % i) continue;
      if (i > d && i < k) k = i;
      const j = l / i; if (j > d && j < k) k = j;
    }
    return { l, k, L: l / k };
  }
  const fmt = (x) => x.toLocaleString('en-US');
  function parse(text) {
    const v = (text.match(/-?\d+/g) || []).map(Number);
    if (v.length < 2) return { err: 'Type at least two terms.' };
    if (v.length > 400) return { err: 'At most 400 terms here; the viewer has 100,000 per sequence.' };
    if (v.some(x => !Number.isSafeInteger(x) || x < 1 || x > 1e12)) return { err: 'Terms must be whole numbers from 1 to 10^12.' };
    for (let i = 1; i < v.length; i++) if (v[i] <= v[i - 1]) return { err: `The terms must increase: ${fmt(v[i - 1])} is followed by ${fmt(v[i])}.` };
    return { v };
  }
  function run() {
    const { v, err } = parse($('#terms').value);
    $('#err').textContent = err || '';
    const body = $('#rows'); body.replaceChildren();
    const pts = [];
    let nw = 0, nl = 0, nt = 0, nn = 0;
    if (!v) { draw(pts); $('#sum').textContent = ''; return; }
    for (let i = 0; i + 1 < v.length; i++) {
      const a = v[i], d = v[i + 1] - a, r = decomp(a, d);
      const tr = document.createElement('tr');
      let cls = 'no', label = `a ≤ 2d: no decomposition`;
      if (r) {
        if (r.k > r.L) { cls = 'lv'; label = 'level (k > L)'; nl++; }
        else { cls = 'wt'; label = r.k === r.L ? 'weight, tie (k = L)' : 'weight (k ≤ L)'; nw++; if (r.k === r.L) nt++; }
        pts.push({ k: r.k, L: r.L, lv: r.k > r.L, i });
      } else nn++;
      const cells = [i + 1, fmt(a), fmt(d), r ? fmt(r.l) : '–', r ? fmt(r.k) : '–', r ? fmt(r.L) : '–',
                     r ? `${fmt(r.k)} × ${fmt(r.L)} + ${fmt(d)}` : '–'];
      for (const c of cells) { const td = document.createElement('td'); td.textContent = c; tr.append(td); }
      const td = document.createElement('td'); td.textContent = label; td.className = cls; tr.append(td);
      tr.dataset.i = i; body.append(tr);
    }
    const dec = nw + nl;
    $('#sum').textContent = `${v.length - 1} terms (the last one has no successor): ${dec} decompose, ` +
      `${nl} in the level class, ${nw} in the weight class` + (nt ? ` (${nt} ties)` : '') + (nn ? `, ${nn} do not decompose.` : '.');
    draw(pts);
  }
  /* the weight-level plate of the example: log k across, log L up, as on the site */
  function draw(pts) {
    const svg = $('#plate'), NS = 'http://www.w3.org/2000/svg';
    svg.replaceChildren();
    const S = 300, M = 24;
    const mx = Math.max(4, ...pts.map(p => Math.max(p.k, p.L)));
    const lg = Math.log(mx);
    const X = (k) => M + Math.log(k) / lg * (S - 2 * M), Y = (L) => S - M - Math.log(L) / lg * (S - 2 * M);
    const el = (n, at) => { const e = document.createElementNS(NS, n); for (const k in at) e.setAttribute(k, at[k]); svg.append(e); return e; };
    el('rect', { x: M, y: M, width: S - 2 * M, height: S - 2 * M, class: 'fr' });
    el('line', { x1: X(1), y1: Y(1), x2: X(mx), y2: Y(mx), class: 'dg' });
    const tk = el('text', { x: S - M, y: S - 6, 'text-anchor': 'end' }); tk.textContent = 'k →';
    const tl = el('text', { x: 6, y: M - 8 }); tl.textContent = '↑ L';
    for (const p of pts) {
      const c = el('circle', { cx: X(p.k), cy: Y(p.L), r: 4, class: p.lv ? 'lv' : 'wt' });
      const t = document.createElementNS(NS, 'title'); t.textContent = `term ${p.i + 1}: k = ${p.k}, L = ${p.L}`; c.append(t);
    }
  }
  document.querySelectorAll('[data-preset]').forEach(b => b.addEventListener('click', () => {
    $('#terms').value = PRESETS[b.dataset.preset]().join(', '); run();
  }));
  $('#terms').addEventListener('input', run);
  $('#terms').value = PRESETS.primes().join(', ');
  run();
})();
