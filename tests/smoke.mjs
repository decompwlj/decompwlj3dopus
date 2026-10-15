/* decompwlj 3D — browser smoke test of the published site, run by .github/workflows/checks.yml.

     python3 -m http.server 8765 &        from the repository root
     node tests/smoke.mjs                  needs the playwright package and its Chromium

   Loads the gallery, opens a sequence from a card (the k–L view) and in the 3-D viewer; for a
   sample of 60 sequences checks that the data decode to their catalogue counts, that the graph
   draws every decomposable term exactly at its log k, log L, log d in its class, that the viewer
   finds the published SHA-256 and that the OEIS link is the sequence's own; downloads a CSV from a
   sequence page, and runs the live example of the "How it works" page.  Fails on any page error
   or any 4xx/5xx response.  PLAYWRIGHT names another place to load playwright from. */
import { createRequire } from 'module';
const require = createRequire(import.meta.url);
const { chromium } = require(process.env.PLAYWRIGHT || 'playwright');
const BASE = process.env.BASE_URL || 'http://localhost:8765';

const fails = [], errs = [];
const check = (ok, what) => { console.log(`${ok ? '  ok  ' : '  FAIL'} ${what}`); if (!ok) fails.push(what); };

const b = await chromium.launch();
const ctx = await b.newContext({ viewport: { width: 1400, height: 900 }, acceptDownloads: true });
const p = await ctx.newPage();
p.on('pageerror', e => errs.push('page error: ' + e.message));
p.on('response', r => { if (r.status() >= 400 && !r.url().includes('goatcounter') && !r.url().includes('gc.zgo.at')) errs.push(`${r.status()} ${r.url()}`); });  /* the visit counter is not ours to test */

/* the gallery */
await p.goto(BASE + '/'); await p.waitForSelector('#grid .card');
const total = Number(await p.textContent('#homeTotal'));
check(total >= 1000, `gallery total ${total}`);
check(await p.locator('#grid .card').count() > 10, 'gallery shows cards');
check(/^Today: A\d{6}/.test(await p.textContent('#today')), 'the gallery shows a sequence of the day');

/* a card opens the flat k–L view; a drag turns it into 3-D */
await p.locator('.card .open').first().click(); await p.waitForTimeout(2500);
check(await p.evaluate(() => decompwlj.state.view) === 'xy', 'a card opens the k–L view');
check(await p.evaluate(() => JSON.parse(localStorage.getItem('dwlj-recent') || '[]')[0] === decompwlj.state.id),
      'the opened sequence is remembered as recently viewed');
await p.mouse.move(700, 450); await p.mouse.down(); await p.mouse.move(640, 410, { steps: 8 }); await p.mouse.up();
check(await p.evaluate(() => decompwlj.state.view) === 'iso', 'a drag turns it into 3-D');

/* a view picked during a flight (as when a gallery card opens a second sequence) still lands with
   the graph inside the camera's depth range, and drawn */
for (const [a, z] of [['xy', 'iso'], ['iso', 'xy']]) {
  await p.evaluate(([a, z]) => { const v = decompwlj.viz; v.view(a); v._tick(performance.now() + 200); v.view(z); }, [a, z]);
  await p.waitForTimeout(1500);
  const r = await p.evaluate(() => {
    const v = decompwlj.viz, c = v.camera, g = v.canvas, x = document.createElement('canvas');
    const dist = c.position.distanceTo(v.controls.target), far = c.far;
    x.width = g.width; x.height = g.height; const k = x.getContext('2d'); k.drawImage(g, 0, 0);
    const d = k.getImageData(0, 0, x.width, x.height).data; let n = 0;
    for (let i = 0; i < d.length; i += 16) if (Math.max(d[i], d[i + 1], d[i + 2]) - Math.min(d[i], d[i + 1], d[i + 2]) > 60) n++;
    return { dist: Math.round(dist), far, n };
  });
  check(r.dist < r.far / 2 && r.n > 200, `${a} → ${z} picked mid-flight: camera at ${r.dist} (far ${r.far}), ${r.n} coloured pixels`);
}

check(await p.locator('#simList li').count() === 8, 'eight similar sequences are listed');

/* a sample of sequences, each in the viewer: the data decode to the catalogue's counts; the graph
   holds one point per decomposable term, at exactly log₁₀ k, log₁₀ L, log₁₀ d in the frame the axes
   show, coloured by its class, with k·L + d = a; the viewer's own check finds the published SHA-256;
   and the OEIS link is the sequence's own */
const res = await p.evaluate(async () => {
  const cat = decompwlj.all, step = Math.max(1, Math.floor(cat.length / 60)), bad = [];
  let ok = 0, points = 0;
  const L10 = (v) => Math.log10(v);
  for (let s = 0; s < cat.length; s += step) {
    const id = cat[s].id, why = [];
    await decompwlj.select(id);
    const ds = decompwlj.ds, v = decompwlj.viz, r = ds && ds.rec;
    if (!(ds && ds.complete && r.id === id && ds.count === r.decomposable)) why.push('load');
    else {
      /* the frame of the axes: k and L from 10⁰, d from its smallest decade; centred */
      const lo = [0, 0, Math.floor(L10(r.dmin))];
      const hi = [r.kmax, r.Lmax, r.dmax].map((m, j) => Math.max(lo[j] + 1, Math.ceil(L10(Math.max(1, m)))));
      const off = [0, 1, 2].map(j => lo[j] + (hi[j] - lo[j]) / 2);
      if (v._n !== ds.count) why.push(`${v._n} points drawn for ${ds.count} terms`);
      let worst = 0, cls = 0, sum = 0;
      for (let i = 0; i < v._n && why.length === 0; i++) {
        const row = ds.row[i], a = ds.a[row], k = ds.k[row], L = ds.L[row], d = ds.d[row];
        if (k * L + d !== a) sum++;
        const want = [L10(k) - off[0], L10(L) - off[1], L10(d) - off[2]];
        for (let j = 0; j < 3; j++) worst = Math.max(worst, Math.abs(v._pos[i * 3 + j] - want[j]));
        if (v._meta[i * 4] !== (k === L ? 2 : k > L ? 1 : 0)) cls++;
      }
      if (worst > 1e-4) why.push(`a point is ${worst.toExponential(1)} decades from its term`);
      if (cls) why.push(`${cls} points in the wrong class`);
      if (sum) why.push(`${sum} terms with k·L + d ≠ a`);
      points += v._n;
      if (document.querySelector('#oeisLink').href !== `https://oeis.org/${r.anumber}`) why.push('OEIS link ' + document.querySelector('#oeisLink').href);
      let txt = '';
      for (const t0 = Date.now(); Date.now() - t0 < 30000; await new Promise(f => setTimeout(f, 40))) {
        txt = document.querySelector('#verify').textContent;
        if (/^(✓ data verified|⚠)/.test(txt)) break;
      }
      if (!txt.startsWith('✓ data verified')) why.push('fingerprint: ' + (txt.slice(0, 80) || 'no answer'));
    }
    if (why.length) bad.push(`${id}: ${why.join('; ')}`); else ok++;
  }
  return { ok, bad, points };
});
check(res.bad.length === 0 && res.ok > 30,
      `sample of sequences in the viewer: ${res.ok} decode to their counts, ${res.points.toLocaleString('en')} points drawn exactly at their terms ` +
      `in their class, fingerprints as published, OEIS links their own${res.bad.length ? '; FAILED ' + res.bad.join(' | ') : ''}`);

/* the viewer verifies the numbers it shows against the published fingerprint */
await p.goto(BASE + '/#A000040');
await p.waitForFunction(() => /verified:|failed/.test(document.querySelector('#verify').textContent), null, { timeout: 60000 });
check((await p.textContent('#verify')).startsWith('✓ data verified'), 'viewer: data verified against the published SHA-256');

/* a sequence page and its CSV */
await p.goto(BASE + '/seq/A000040/');
check((await p.getAttribute('meta[property="og:image"]', 'content')).endsWith('/share/primes.jpg'), 'sequence page has its share image');
const img = await p.request.get(BASE + '/share/primes.jpg');
check(img.ok() && (img.headers()['content-type'] || '').includes('jpeg'), 'share image is served');
const [dl] = await Promise.all([p.waitForEvent('download'), p.click('#csv')]);
const csv = await (await import('fs')).promises.readFile(await dl.path(), 'utf8');
const lines = csv.trim().split('\n');
check(lines[0] === 'n;a;weight;level;jump' && lines.length === 100001, `CSV has a header and 100000 rows (${lines.length - 1})`);
check(lines[30] === '30;113;33;3;14', `CSV row of 113 reads 30;113;33;3;14 (${lines[30]})`);

/* how it works: the live example */
await p.goto(BASE + '/learn/'); await p.waitForSelector('#rows tr');
check((await p.textContent('#sum')).includes('36 decompose, 19 in the level class, 17 in the weight class'), 'live example: counts for the first 40 primes');
await p.fill('#terms', '1009, 1013'); await p.waitForTimeout(100);
check((await p.locator('#rows tr').first().innerText()).replace(/\s+/g, ' ').includes('5 × 201 + 4'), 'live example: 1009 = 5 × 201 + 4');

await b.close();
for (const e of errs) console.log('  FAIL', e);
if (fails.length || errs.length) { console.log(`\n${fails.length + errs.length} problem(s)`); process.exit(1); }
console.log('\nsmoke test passed');
