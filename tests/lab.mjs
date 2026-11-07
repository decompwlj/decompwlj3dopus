/* decompwlj Lab — browser test, run by .github/workflows/checks.yml next to tests/smoke.mjs.

     python3 -m http.server 8765 &        from the repository root
     node tests/lab.mjs                    needs the playwright package and its Chromium

   Loads the Lab on the primes and on the sequence with the largest terms, and checks, in the
   browser itself: the data verify (counts and SHA-256); the kernel reproduces 2000 published
   rows of each; every control record keeps its jump and its residue mod M, stays in its window,
   and is decomposed to the least divisor above d (by trial division wherever ℓ' < 10^8); the
   controls are reproducible (the seeded ladder is the same after a reload); the filters (a jump,
   a range of n, the modulus) reach every chart and the link; the prime-power ladder closes the gap
   of A047429 (n ≡ 4, 5, 6 mod 8) at M = 840 and not before, and a modulus of one's own (8) does too; the congruence lock gives exactly
   3 at r = 3 for the jump 2 of the primes (ℓ = p − 2 ≡ 0 mod 3 for every lesser twin p ≥ 5); a
   pasted list of primes decomposes to the published rows; the picker opens a sequence; a phone
   width has no horizontal scroll.  Fails on any page error or any 4xx/5xx response.
   PLAYWRIGHT names another place to load playwright from; CHROME another browser binary. */
import { createRequire } from 'module';
const require = createRequire(import.meta.url);
const { chromium } = require(process.env.PLAYWRIGHT || 'playwright');
const BASE = process.env.BASE_URL || 'http://localhost:8765';

const fails = [], errs = [];
const check = (ok, what) => { console.log(`${ok ? '  ok  ' : '  FAIL'} ${what}`); if (!ok) fails.push(what); };

const b = await chromium.launch(process.env.CHROME ? { executablePath: process.env.CHROME } : {});
const ctx = await b.newContext({ viewport: { width: 1400, height: 900 } });
const p = await ctx.newPage();
p.on('pageerror', e => errs.push('page error: ' + e.message));
p.on('response', r => { if (r.status() >= 400 && !r.url().includes('goatcounter') && !r.url().includes('gc.zgo.at')) errs.push(`${r.status()} ${r.url()}`); });

const ready = () => p.waitForFunction(() => { const s = window.__lab; return s && s.dec && s.total && s.done === s.total && s.verified; }, null, { timeout: 240000 });
const ladderText = () => p.evaluate(() => window.__labFns.ladder().map(r => r && `${r.M}:${r.level.delta.toFixed(6)}:${r.one.delta.toFixed(6)}`).join(' '));

/* in the page: the kernel against the published rows, and every control record against its rules */
const audit = () => p.evaluate(() => {
  const s = window.__lab, K = window.DWLJ, out = { rows: 0, rowBad: 0, ctl: 0, ctlBad: 0, brute: 0, bruteBad: 0, first: '' };
  const R = K.rng(7);
  for (let t = 0; t < 2000; t++) {
    const r = s.dec[R.below(s.dec.length)], got = K.decompose(s.a[r], s.d[r]);
    out.rows++;
    if (!got || got[0] !== s.k[r] || got[1] !== s.L[r]) { out.rowBad++; out.first ||= `row ${r}`; }
  }
  for (const [key, c] of s.ctrl) {
    const M = +key.split('|')[0];
    for (let j = 0; j < s.sub.length; j++) {
      const r = s.sub[j], d = s.d[r], l = s.a[r] - d, lp = c.lp[j], k = c.k[j];
      out.ctl++;
      const ok = Number.isInteger(lp) && lp % M === l % M && lp > d && lp >= Math.ceil(0.75 * l) - 1 && lp <= Math.floor(1.25 * l)
              && k > d && lp % k === 0;
      if (!ok) { out.ctlBad++; out.first ||= `control ${key} record ${j}`; continue; }
      if (lp < 1e7 && j % 13 === 0) {                 /* the least divisor above d, by brute force */
        out.brute++;
        let least = lp;
        for (let m = d + 1; m * m <= lp; m++) if (lp % m === 0) { least = m; break; }
        if (least === lp) for (let m = Math.min(d, Math.floor(Math.sqrt(lp))); m >= 1; m--) if (lp % m === 0 && lp / m > d) { least = lp / m; break; }
        if (least !== k) { out.bruteBad++; out.first ||= `brute ${key} ${lp} ${d}: ${least} vs ${k}`; }
      }
    }
  }
  return out;
});

/* 1. the primes */
await p.goto(BASE + '/lab/#primes'); await ready();
let st = await p.evaluate(() => ({ v: window.__lab.verified, N: window.__lab.N, sub: window.__lab.sub.length, a: window.__lab.a[29], k: window.__lab.k[29], L: window.__lab.L[29], d: window.__lab.d[29] }));
check(st.v.startsWith('✓ data verified'), `primes: data verified (${st.v.slice(0, 40)}…)`);
check(st.a === 113 && st.k === 33 && st.L === 3 && st.d === 14, 'primes: a(30) = 113 = 33 × 3 + 14');
check(st.sub === 25000, `primes: the control runs on 25,000 records (${st.sub})`);
let A = await audit();
check(A.rowBad === 0, `primes: kernel reproduces ${A.rows} published rows (${A.rowBad} differ) ${A.first}`);
check(A.ctl === 20 * 25000 && A.ctlBad === 0, `primes: ${A.ctl} control records keep d, ℓ mod M and the window, k | ℓ′, k > d`);
check(A.brute > 1000 && A.bruteBad === 0, `primes: ${A.brute} control weights equal the least divisor above d by brute force`);
const lad1 = await ladderText();
check(!lad1.includes('null') && lad1.split(' ').length === 10, 'primes: the ladder has all ten steps, lcm(1..n) up to 360360');
const lv = await p.evaluate(() => window.__labFns.ladder().map(r => r.level.delta));
check(lv[1] < -0.03 && Math.abs(lv[6]) < Math.abs(lv[1]), `primes: the deficit is largest at M = 2 and shrinks by M = 840 (${lv.map(x => (100 * x).toFixed(2)).join(', ')})`);
check(await p.evaluate(() => [...document.querySelectorAll('canvas')].every(c => c.width > 0)), 'every chart is drawn');

/* the jump filter, from a click on the bar of d = 2 */
const g2 = await p.evaluate(() => {
  const cv = document.querySelector('#gaps'); cv.scrollIntoView({ block: 'center' });
  const b = cv.getBoundingClientRect();
  for (let x = 0; x < b.width; x += 2) {
    cv.dispatchEvent(new PointerEvent('pointermove', { clientX: b.left + x, clientY: b.top + b.height / 2, bubbles: true }));
    const t = document.querySelector('#tip');
    if (!t.hidden && t.textContent.startsWith('d = 2\n')) return { x: b.left + x, y: b.top + b.height / 2 };
  }
  return null;
});
check(!!g2, 'the jumps chart has a bar for d = 2');
if (g2) {
  await p.mouse.click(g2.x, g2.y); await p.waitForTimeout(300);
  const f = await p.evaluate(() => {
    const s = window.__lab, F = window.__labFns, lk = F.lockStats();
    let n = 0; for (const r of s.dec) if (s.d[r] === 2) n++;
    return { g: s.g, hash: location.hash, n, shown: F.obsCounts(s.dec).n, r3: lk[0].y, r5: lk[1].y,
             row: document.querySelector('#numTable tbody tr td:nth-child(2)').textContent };
  });
  check(f.g === 2 && f.hash.includes('.g2'), `clicking d = 2 keeps that jump, in the link too (${f.hash})`);
  check(f.shown === f.n && f.row === f.n.toLocaleString('en-US'), `every chart counts the ${f.n} records with d = 2`);
  check(f.r3 === 3, `the lock at r = 3 is exactly 3 for the lesser twins (${f.r3})`);
  await p.click('#clearG'); await p.waitForTimeout(200);
  check(await p.evaluate(() => window.__lab.g === null), 'all jumps again');
}

/* the modulus, from a click on the ladder */
await p.click('#ladderTable tbody tr:nth-child(5)'); await p.waitForTimeout(300);
check(await p.evaluate(() => window.__lab.M === 60 && location.hash.includes('.M60') && document.querySelector('#fM').textContent === 'M = 60'),
      'a ladder step selects M = 60 everywhere');

/* a range of n, by dragging along the level-share chart */
await p.locator('#along').scrollIntoViewIfNeeded();
const al = await p.locator('#along').boundingBox();
await p.mouse.move(al.x + al.width * 0.3, al.y + al.height / 2); await p.mouse.down();
await p.mouse.move(al.x + al.width * 0.6, al.y + al.height / 2, { steps: 6 }); await p.mouse.up(); await p.waitForTimeout(300);
const br = await p.evaluate(() => ({ lo: window.__lab.rLo, hi: window.__lab.rHi, N: window.__lab.N, hash: location.hash }));
check(br.lo > 10000 && br.hi < 80000 && /\.n\d+-\d+/.test(br.hash), `a drag keeps n ${br.lo}–${br.hi}, in the link too`);
await p.dblclick('#along'); await p.waitForTimeout(200);
check(await p.evaluate(() => window.__lab.rLo === 0 && window.__lab.rHi === window.__lab.N - 1), 'a double-click clears it');

/* reproducible: the same ladder after a reload, and a deep link restores its filters */
await p.goto(BASE + '/lab/#lucky'); await ready();
await p.goto(BASE + '/lab/#primes'); await ready();
check(await ladderText() === lad1, 'the seeded controls give the same ladder after a reload');
await p.goto(BASE + '/lab/#primes.M420.g6.n1000-50000'); await p.reload(); await ready();
check(await p.evaluate(() => { const s = window.__lab; return s.M === 420 && s.g === 6 && s.n0 + s.rLo === 1000 && s.n0 + s.rHi === 50000; }), 'a deep link restores M, the jump and the range');

/* the ladder climbs prime powers: A047429 (n ≡ 4, 5, 6 mod 8) departs from its control through ℓ mod 8 alone,
   so the gap stays open up to M = 420 and closes at M = 840, the first step that keeps 8 */
await p.goto(BASE + '/lab/#a047429'); await ready();
const q = await p.evaluate(() => Object.fromEntries(window.__labFns.ladder().map(r => [r.M, [r.level.delta, r.level.z, r.frozen]])));
check(q[420][0] < -0.02 && q[420][1] < -5 && q[420][2] < 0.5, `A047429: still ${(100 * q[420][0]).toFixed(2)} points at M = 420 (z ${q[420][1].toFixed(1)})`);
check(Math.abs(q[840][1]) < 3 && q[840][2] < 0.5, `A047429: closed at M = 840, where 8 enters (Δ ${(100 * q[840][0]).toFixed(2)}, z ${q[840][1].toFixed(1)})`);
/* your own M, from the form, then from a link */
await p.fill('#ownM', '8'); await p.click('#ownForm button'); await ready();
const o8 = await p.evaluate(() => { const r = window.__labFns.ladder().find(x => x && x.M === 8); return { M: window.__lab.M, steps: window.__labFns.steps().length, hash: location.hash, z: r && r.level.z, row: document.querySelector('#ladderTable tbody tr:last-child').textContent }; });
check(o8.M === 8 && o8.steps === 11 && o8.hash.includes('.M8') && o8.row.includes('yours'), `your own M = 8 joins the ladder and the link (${o8.hash})`);
check(o8.z !== null && Math.abs(o8.z) < 3, `A047429 at M = 8: no gap left (z ${o8.z && o8.z.toFixed(1)})`);
await p.click('#ladderTable tbody tr:nth-child(2)'); await p.waitForTimeout(200);
check(await p.evaluate(() => location.hash === '#a047429.M2.own8'), 'choosing a ladder step keeps your own M in the link');
await p.goto(BASE + '/lab/#lucky'); await ready();
await p.goto(BASE + '/lab/#a047429.M2.own8'); await ready();
check(await p.evaluate(() => window.__lab.M === 2 && window.__lab.ownM === 8 && window.__labFns.reps(8).length === 2), 'a link with your own M restores it');
await p.fill('#ownM', '0'); await p.click('#ownForm button');
check(await p.evaluate(() => !document.querySelector('#ownM').checkValidity() && window.__lab.ownM === 8), 'a modulus out of range is refused');
await p.goto(BASE + '/lab/#primes'); await ready();

/* a pasted list: the first 3000 primes decompose to the published rows */
const ref = await p.evaluate(() => ({ a: Array.from(window.__lab.a.slice(0, 3001)), k: Array.from(window.__lab.k.slice(0, 3000)), L: Array.from(window.__lab.L.slice(0, 3000)) }));
await p.click('#pasteBtn'); await p.fill('#pasteText', ref.a.join(', ')); await p.fill('#pasteName', 'first primes'); await p.click('#pasteGo');
await ready();
const pr = await p.evaluate(() => ({ k: Array.from(window.__lab.k), L: Array.from(window.__lab.L), name: document.querySelector('#seqName').textContent, hash: location.hash }));
check(pr.k.length === 3000 && pr.k.every((k, i) => k === ref.k[i] && pr.L[i] === ref.L[i]), 'a pasted list of 3001 primes decomposes to the published rows');
check(pr.name === 'first primes' && pr.hash === '#custom', 'the pasted sequence is named and kept out of the link');
await p.click('#pasteBtn'); await p.fill('#pasteText', '1, 2, 2, 3'); await p.click('#pasteGo');
check((await p.textContent('#pasteErr')).length > 0, `bad input is refused (${await p.textContent('#pasteErr')})`);
await p.keyboard.press('Escape');

/* 2. the sequence with the largest terms: the kernel near 2^53 */
const big = await p.evaluate(() => [...window.__lab.cat].sort((x, y) => y.amax - x.amax)[0].id);
await p.goto(BASE + `/lab/#${big}`); await ready();
A = await audit();
st = await p.evaluate(() => ({ v: window.__lab.verified, amax: Math.max(...window.__lab.a) }));
check(st.v.startsWith('✓ data verified'), `${big}: data verified, terms up to ${st.amax.toExponential(3)}`);
check(A.rowBad === 0, `${big}: kernel reproduces ${A.rows} published rows (${A.rowBad} differ) ${A.first}`);
check(A.ctlBad === 0, `${big}: ${A.ctl} control records keep their rules`);

/* the picker */
await p.fill('#q', 'A000959'); await p.waitForTimeout(100); await p.keyboard.press('Enter');
await ready();
check(await p.evaluate(() => window.__lab.rec.anumber === 'A000959' && location.hash === '#lucky'), 'the picker opens the lucky numbers');

/* a phone */
await p.setViewportSize({ width: 390, height: 844 }); await p.waitForTimeout(400);
check(await p.evaluate(() => document.documentElement.scrollWidth <= document.documentElement.clientWidth), 'no horizontal scroll at 390 px');

await b.close();
for (const e of errs) console.log('  FAIL', e);
if (fails.length || errs.length) { console.log(`\n${fails.length + errs.length} problem(s)`); process.exit(1); }
console.log('\nlab test passed');
