/* decompwlj 3D — browser smoke test of the published site, run by .github/workflows/checks.yml.

     python3 -m http.server 8765 &        from the repository root
     node tests/smoke.mjs                  needs the playwright package and its Chromium

   Loads the gallery, opens a sequence from a card (the k–L view) and in the 3-D viewer, checks
   that the data of a sample of sequences decodes to their catalogue counts, downloads a CSV from a
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
p.on('response', r => { if (r.status() >= 400) errs.push(`${r.status()} ${r.url()}`); });

/* the gallery */
await p.goto(BASE + '/'); await p.waitForSelector('#grid .card');
const total = Number(await p.textContent('#homeTotal'));
check(total >= 1000, `gallery total ${total}`);
check(await p.locator('#grid .card').count() > 10, 'gallery shows cards');

/* a card opens the flat k–L view; a drag turns it into 3-D */
await p.locator('.card .open').first().click(); await p.waitForTimeout(2500);
check(await p.evaluate(() => decompwlj.state.view) === 'xy', 'a card opens the k–L view');
await p.mouse.move(700, 450); await p.mouse.down(); await p.mouse.move(640, 410, { steps: 8 }); await p.mouse.up();
check(await p.evaluate(() => decompwlj.state.view) === 'iso', 'a drag turns it into 3-D');

/* data of a sample of sequences: every chunk decodes to the catalogue's counts */
const res = await p.evaluate(async () => {
  const cat = decompwlj.all, step = Math.max(1, Math.floor(cat.length / 60)), bad = [];
  let ok = 0;
  for (let i = 0; i < cat.length; i += step) {
    await decompwlj.select(cat[i].id); const ds = decompwlj.ds;
    if (ds && ds.complete && ds.rec.id === cat[i].id && ds.count === ds.rec.decomposable) ok++; else bad.push(cat[i].id);
  }
  return { ok, bad };
});
check(res.bad.length === 0 && res.ok > 30, `sample of sequences loads in the viewer (${res.ok} ok${res.bad.length ? ', bad ' + res.bad.join(' ') : ''})`);

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
