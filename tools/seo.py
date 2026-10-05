#!/usr/bin/env python3
"""seo.py — the crawlable side of the site: static pages, sitemap and robots.txt.

    python3 seo.py                     base URL from ../CNAME (https://decompwlj.net)
    python3 seo.py --base https://user.github.io/repo

The app itself is one page driven by the URL hash (#A000040), which search engines
do not index. This script writes plain HTML pages that they do, all from
data/catalog.csv, next to index.html:

    seq/<A-number>/index.html     one page per sequence: OEIS name, plate, counts, note,
                                  links to the 3-D viewer, the OEIS and its neighbours
    seq/index.html                every sequence, by A-number
    learn/index.html, learn.js    how it works: the decomposition step by step, with a live example
    family/<family>/index.html    the sequences of one family
    404.html                      not-found page; /A000040 redirects to /seq/A000040/
    sitemap.xml, robots.txt       for the crawlers
    og.png                        1200 x 630 preview for links shared on social sites
    share/<id>.jpg                600 x 315 preview of each sequence page
                                  (both need Pillow; skipped without it)

Run it again whenever the catalogue changes. The output is deterministic.
"""
import csv, html, json, pathlib, re, sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
SITE = 'decompwlj 3D'
# visit counts: GoatCounter, no cookies and no personal data (see index.html); '' leaves them out
GOATCOUNTER_URL = 'https://decompwlj.goatcounter.com/count'
GOATCOUNTER = (f'<script data-goatcounter="{GOATCOUNTER_URL}" async src="https://gc.zgo.at/count.js"></script>'
               if GOATCOUNTER_URL else '')

def base_url():
    if '--base' in sys.argv:
        return sys.argv[sys.argv.index('--base') + 1].rstrip('/')
    cname = ROOT / 'CNAME'
    if cname.exists() and cname.read_text().strip():
        return 'https://' + cname.read_text().strip().rstrip('/')
    sys.exit('no CNAME: pass --base https://<user>.github.io/<repo>')

BASE = base_url()
rows = list(csv.DictReader(open(ROOT / 'data' / 'catalog.csv', newline='')))
rows.sort(key=lambda r: r['anumber'])
N = len(rows)
e = lambda s: html.escape(str(s), quote=True)
num = lambda v: f'{int(v):,}'
def slug(f): return re.sub(r'[^a-z0-9]+', '-', f.lower()).strip('-')
def pct(a, b): return f'{100 * a / b:.2f} %' if b else '–'
def clip(s, n):
    return s if len(s) <= n else s[:n - 1].rsplit(' ', 1)[0].rstrip(',;:') + '…'

# the published fingerprints (tools/fingerprint.py): the SHA-256 of every sequence's CSV
SHA = {}
if (ROOT / 'data' / 'sha256.txt').exists():
    for l in (ROOT / 'data' / 'sha256.txt').read_text().splitlines():
        if l.strip(): h, f = l.split(); SHA[f[len('decompwlj_'):-len('.csv')]] = h
families = {}
for r in rows: families.setdefault(r['family'], []).append(r)
fam_order = sorted(families, key=lambda f: (-len(families[f]), f.lower()))

# the favicon: files at the site root (favicon.svg, favicon.ico, apple-touch-icon.png), light theme
def icons(up):
    return (f'<link rel="icon" href="{up}favicon.ico" sizes="32x32">\n<link rel="icon" href="{up}favicon.svg" type="image/svg+xml">\n'
            f'<link rel="apple-touch-icon" href="{up}apple-touch-icon.png">')
PM, PX = 0.03 / 1.06, 1.03 / 1.06
PLATE_SVG = (f'<svg viewBox="0 0 1 1" preserveAspectRatio="none" aria-hidden="true">'
             f'<rect x="{PM}" y="{PM}" width="{PX - PM}" height="{PX - PM}"/>'
             f'<line class="ed" x1="{PM}" y1="{PM}" x2="{PX}" y2="{PX}"/>'
             f'<line class="dg" x1="{PM}" y1="{PX}" x2="0.5" y2="0.5"/></svg>')

def page(up, path, title, desc, body, image=None, ld=()):
    """A complete page. up: relative path back to the site root ('../../'). image: a 600 x 315 share image."""
    url = f'{BASE}/{path}'
    img, (iw, ih) = (image, (600, 315)) if image else (f'{BASE}/og.png', (1200, 630))
    lds = ''.join(f'<script type="application/ld+json">{json.dumps(x, ensure_ascii=False, separators=(",", ":"))}</script>\n'
                  for x in ld)
    return f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{e(title)}</title>
<meta name="description" content="{e(desc)}">
<link rel="canonical" href="{e(url)}">
{icons(up)}
<link rel="stylesheet" href="{up}css/pages.css">
{GOATCOUNTER}
<meta property="og:type" content="website">
<meta property="og:site_name" content="{SITE}">
<meta property="og:title" content="{e(title)}">
<meta property="og:description" content="{e(desc)}">
<meta property="og:url" content="{e(url)}">
<meta property="og:image" content="{e(img)}">
<meta property="og:image:width" content="{iw}">
<meta property="og:image:height" content="{ih}">
<meta name="twitter:card" content="summary_large_image">
{lds}</head>
<body>
<header class="top">
  <a class="brand" href="{up}"><span class="brand-mark" aria-hidden="true"></span>decompwlj<span class="brand-dim"> 3D</span></a>
  <nav aria-label="Site"><a href="{up}">Gallery</a><a href="{up}learn/">How it works</a><a href="{up}seq/">All sequences</a><a href="https://decompwlj.com" rel="noopener">decompwlj.com</a></nav>
</header>
<main>
{body}
</main>
<footer class="foot">
  <span>{SITE} — the decomposition a(n) = k·L + d for {N} integer sequences</span>
  <a href="https://arxiv.org/abs/0711.0865" rel="noopener">arXiv:0711.0865</a>
  <a href="https://oeis.org" rel="noopener">OEIS</a>
</footer>
</body>
</html>
'''

def crumbs(up, items):
    """items: [(name, relative href or None)] after Home."""
    parts = [f'<a href="{up}">Home</a>'] + [f'<a href="{h}">{e(n)}</a>' if h else f'<span aria-current="page">{e(n)}</span>'
                                             for n, h in items]
    return '<nav class="crumbs" aria-label="Breadcrumb">' + ' › '.join(parts) + '</nav>'

def crumbs_ld(items):
    lst = [('Home', f'{BASE}/')] + items
    return {'@context': 'https://schema.org', '@type': 'BreadcrumbList',
            'itemListElement': [{'@type': 'ListItem', 'position': i + 1, 'name': n, 'item': u} for i, (n, u) in enumerate(lst)]}

ABOUT = ('<h2>The decomposition</h2>\n<p class="lead">Every term of a strictly increasing sequence is written '
         '<i>a</i>(n) = <i>k</i>(n)·<i>L</i>(n) + <i>d</i>(n): the jump <i>d</i> = a(n+1) − a(n), the weight <i>k</i> the least '
         'divisor of a − d greater than d, the level <i>L</i> = (a − d)/k. A term decomposes when a &gt; 2d; it is in the level '
         'class when k &gt; L and in the weight class when k ≤ L. See <a href="https://decompwlj.com" rel="noopener">decompwlj.com</a> '
         'and <a href="https://arxiv.org/abs/0711.0865" rel="noopener">arXiv:0711.0865</a>, or '
         '<a href="../../learn/">how it works</a>, with worked examples and a live one.</p>')

CSV_JS = r"""/* decompwlj 3D, written by tools/seo.py: the "Download CSV" button of a sequence page.
   Rebuilds n;a;weight;level;jump from the data chunks (the dwj2 format of tools/compact.py),
   with every check the viewer makes, then compares the file's SHA-256 with the one published
   in data/sha256.txt (shown on the page) and saves it only if they are equal. Nothing is stored. */
(function () {
  var btn = document.getElementById('csv');
  if (!btn) return;
  function inflate(bytes) {
    if (bytes.length > 2 && bytes[0] === 0x1f && bytes[1] === 0x8b) {
      var s = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
      return new Response(s).arrayBuffer().then(function (b) { return new Uint8Array(b); });
    }
    return Promise.resolve(bytes);
  }
  function rows(b, url, want, out) {
    if (!(b[0] === 0x64 && b[1] === 0x77 && b[2] === 0x6a && (b[3] === 0x31 || b[3] === 0x32))) throw new Error(url + ': not a dwj1/dwj2 chunk');
    var m = b[3] === 0x32 ? b[4] : 0, i = b[3] === 0x32 ? 5 : 4;
    if (m > 2) throw new Error(url + ': unknown jump coding ' + m);
    function unzz(z) { return z % 2 === 0 ? z / 2 : -(z + 1) / 2; }
    function rd() { var v = 0, m = 1, c; do { if (i >= b.length) throw new Error(url + ': truncated'); c = b[i++]; v += (c & 127) * m; m *= 128; } while (c & 128); return v; }
    var n = rd(), a = rd(), d = new Array(n);
    if (n !== want) throw new Error(url + ': ' + n + ' rows, expected ' + want);
    for (var r = 0, p = 0, pp = 0; r < n; r++) {
      var v = rd(); d[r] = m === 0 ? v : m === 1 ? p + unzz(v) : 2 * p - pp + unzz(v); pp = p; p = d[r];
      if (!(d[r] > 0) || d[r] >= 9007199254740992) throw new Error(url + ': jump out of range at row ' + r);
    }
    for (r = 0; r < n; r++) {
      var s = rd(), k = 0, L = 0;
      if (s) {
        var h = Math.floor(s / 2), q = (a - d[r]) / h;
        k = s % 2 ? q : h; L = s % 2 ? h : q;
        if (!Number.isInteger(q) || q * h !== a - d[r] || k <= d[r] || (s % 2 ? !(L < k) : !(k <= L)))
          throw new Error(url + ': row ' + r + ' does not fit a = ' + a);
      } else if (a > 2 * d[r]) throw new Error(url + ': row ' + r + ' should decompose');
      out.push([a, k, L, d[r]]);
      a += d[r];
    }
    if (i !== b.length) throw new Error(url + ': trailing bytes');
  }
  function hex(buf) { return Array.prototype.map.call(new Uint8Array(buf), function (x) { return ('0' + x.toString(16)).slice(-2); }).join(''); }
  btn.addEventListener('click', function () {
    var id = btn.dataset.id, an = btn.dataset.an, n0 = +btn.dataset.n0, nc = +btn.dataset.chunks;
    var terms = +btn.dataset.terms, per = +btn.dataset.rows, want = btn.dataset.sha256;
    var label = btn.textContent; btn.disabled = true; btn.textContent = 'Building…';
    var jobs = [];
    for (var c = 0; c < nc; c++) {
      (function (c) {
        var url = '../../data/seq/' + id + '/chunk-' + String(c).padStart(3, '0') + '.bin.gz';
        jobs.push(fetch(url).then(function (res) { if (!res.ok) throw new Error(res.status + ' ' + url); return res.arrayBuffer(); })
          .then(function (buf) { return inflate(new Uint8Array(buf)); })
          .then(function (b) { var out = []; rows(b, url, Math.min(per, terms - c * per), out); return out; }));
      })(c);
    }
    var text;
    Promise.all(jobs).then(function (parts) {
      for (var c = 1; c < parts.length; c++) {                /* each chunk starts where the last ended */
        var e = parts[c - 1][parts[c - 1].length - 1];
        if (parts[c][0][0] !== e[0] + e[3]) throw new Error('chunk ' + c + ' does not continue chunk ' + (c - 1));
      }
      text = ['n;a;weight;level;jump\n']; var buf = '', n = n0;
      parts.forEach(function (p) { p.forEach(function (r) {
        buf += n++ + ';' + r[0] + ';' + (r[1] || '') + ';' + (r[1] ? r[2] : '') + ';' + r[3] + '\n';
        if (buf.length > 65536) { text.push(buf); buf = ''; }
      }); });
      text.push(buf);
      if (n - n0 !== terms) throw new Error((n - n0) + ' rows, expected ' + terms);
      if (!want || !window.crypto || !crypto.subtle) return null;
      return new Blob(text).arrayBuffer().then(function (b) { return crypto.subtle.digest('SHA-256', b); });
    }).then(function (digest) {
      /* never save a file that differs from the published one */
      if (digest && hex(digest) !== want)
        throw new Error('the rebuilt file differs from the published one (SHA-256 ' + hex(digest).slice(0, 12) + '… instead of ' + want.slice(0, 12) + '…). Reload the page (Ctrl+F5); if it persists, please report it');
      var url = URL.createObjectURL(new Blob(text, { type: 'text/csv' }));
      var a = document.createElement('a'); a.href = url; a.download = 'decompwlj_' + an + '.csv';
      document.body.appendChild(a); a.click(); a.remove();
      try { window.goatcounter.count({ path: 'csv/' + an, title: 'CSV ' + an, event: true }); } catch (e) {}
      setTimeout(function () { URL.revokeObjectURL(url); }, 30000);
    }).catch(function (e) { alert('Could not build the CSV: ' + e.message); })
      .then(function () { btn.disabled = false; btn.textContent = label; });
  });
})();
"""

written = []
def write(rel, text):
    p = ROOT / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    written.append(rel)

# ── one page per sequence ────────────────────────────────────────────────────
def kw_html(r):
    """the OEIS keywords core and nice, each a link to every sequence the OEIS gives it"""
    ks = (r.get('keywords') or '').split()
    return (' · OEIS ' + ', '.join(f'<a class="kw" href="https://oeis.org/search?q=keyword:{k}" rel="noopener" '
                                    f'title="Every sequence the OEIS marks {k}">{k}</a>' for k in ks)) if ks else ''

for i, r in enumerate(rows):
    A, name, fam = r['anumber'], r['name'], r['family']
    up = '../../'
    dec, lev, wt = int(r['decomposable']), int(r['level']), int(r['weight'])
    terms, n0 = int(r['terms']), int(r['n0'])
    title = f'{clip(name, 64)} ({A}) · weight × level + jump'
    desc = clip(f'{A}, {name}: {num(terms)} terms decomposed as a(n) = k·L + d; {pct(lev, dec)} in the level class (k > L), '
                f'{pct(wt, dec)} in the weight class. {r["note"]}', 300)
    prev_r, next_r = rows[i - 1] if i else None, rows[i + 1] if i + 1 < N else None
    stats = [
        ('Terms', f'{num(terms)} <span class="dim">(n = {n0} … {n0 + terms - 1:,})</span>'),
        ('Decomposable (a &gt; 2d)', num(dec)),
        ('Level class, k &gt; L', f'<span class="lv">{num(lev)} · {pct(lev, dec)}</span>'),
        ('Weight class, k ≤ L', f'<span class="wt">{num(wt)} · {pct(wt, dec)}</span>'),
        ('Ties, k = L', num(r['ties'])),
        ('On the level line L = 1', num(r['level_one'])),
        ('Forced level, l ≤ d²', num(r['forced'])),
        ('Range of a(n)', f'{num(r["amin"])} … {num(r["amax"])}'),
        ('Range of the jump d', f'{num(r["dmin"])} … {num(r["dmax"])}'),
        ('Largest weight k, level L', f'{num(r["kmax"])}, {num(r["Lmax"])}'),
    ]
    body = f'''{crumbs(up, [('All sequences', '../'), (fam, f'{up}family/{slug(fam)}/'), (A, None)])}
<h1>{e(name)}</h1>
<p class="meta"><a class="an" href="https://oeis.org/{A}" rel="noopener">{A}</a> on the OEIS · family <a href="{up}family/{slug(fam)}/">{e(fam)}</a>{f' · also known as {e(r["alias"])}' if r['alias'] else ''}{kw_html(r)}</p>
<div class="seq">
  <figure>
    <a class="plate-link" href="{up}#{A}" title="Open in the 3-D viewer: it starts on this plate; drag it to turn it into 3-D"><div class="plate">{PLATE_SVG}<img src="{up}thumbs/{r['id']}.webp" width="480" height="480" alt="Weight–level plate of {e(name)}"><span class="axl k" aria-hidden="true">k</span><span class="axl L" aria-hidden="true">L</span><span class="go3d" aria-hidden="true">Explore in 3-D →</span></div></a>
    <figcaption>Click the plate to explore it in 3-D. Weight k across, level L up, both on log scales: blue in the weight class (k ≤ L), orange in the level class (k &gt; L). The dashed diagonal is k = L.</figcaption>
  </figure>
  <div>
    <a class="cta" href="{up}#{A}">Open in the 3-D viewer</a><a class="cta alt" href="https://oeis.org/{A}" rel="noopener">{A} on the OEIS</a><button type="button" class="cta alt" id="csv" data-id="{r['id']}" data-an="{A}" data-n0="{n0}" data-terms="{terms}" data-chunks="{r['chunks']}" data-rows="{r['chunk_rows']}" data-sha256="{SHA.get(A, '')}" title="n;a;weight;level;jump, one row per term, rebuilt from the data and checked against its published SHA-256">Download CSV</button>
    <table class="stats">
{''.join(f'      <tr><th scope="row">{k}</th><td>{v}</td></tr>{chr(10)}' for k, v in stats)}    </table>
    <p class="note">{e(r['note'])}</p>{f'{chr(10)}    <p class="sha">SHA-256 of the CSV <code>decompwlj_{A}.csv</code>, checked before it is saved: <code class="h">{SHA[A]}</code> (all of them: <a href="{up}data/sha256.txt">sha256.txt</a>, for <code>sha256sum -c</code>)</p>' if A in SHA else ''}
  </div>
</div>
{ABOUT}
<script src="../csv.js" defer></script>
<nav class="pn" aria-label="Neighbouring sequences">{f'<a class="prev" href="../{prev_r["anumber"]}/" title="{e(prev_r["name"])}">← {prev_r["anumber"]} {e(clip(prev_r["name"], 50))}</a>' if prev_r else ''}{f'<a class="next" href="../{next_r["anumber"]}/" title="{e(next_r["name"])}">{next_r["anumber"]} {e(clip(next_r["name"], 50))} →</a>' if next_r else ''}</nav>'''
    ld = [{'@context': 'https://schema.org', '@type': 'Dataset', 'name': f'{A}: {name} — weight × level + jump',
           'description': desc if len(desc) >= 50 else desc + ' Decomposition into weight × level + jump.',
           'url': f'{BASE}/seq/{A}/', 'identifier': A, 'isBasedOn': f'https://oeis.org/{A}',
           'keywords': ['integer sequence', 'OEIS', A, fam, 'weight', 'level', 'jump', 'decomposition'],
           'image': f'{BASE}/thumbs/{r["id"]}.webp',
           'includedInDataCatalog': {'@type': 'DataCatalog', 'name': SITE, 'url': f'{BASE}/'},
           'distribution': [{'@type': 'DataDownload', 'encodingFormat': 'application/gzip',
                             'contentUrl': f'{BASE}/data/seq/{r["id"]}/chunk-{c:03d}.bin.gz'} for c in range(int(r['chunks']))]},
          crumbs_ld([('All sequences', f'{BASE}/seq/'), (fam, f'{BASE}/family/{slug(fam)}/'), (A, f'{BASE}/seq/{A}/')])]
    write(f'seq/{A}/index.html', page(up, f'seq/{A}/', title, desc, body, image=f'{BASE}/share/{r["id"]}.jpg', ld=ld))

# ── lists: all sequences, and one per family ─────────────────────────────────
def chips(up, current=None):
    cur = ' aria-current="page"'
    return '<ul class="chips">' + ''.join(
        f'<li><a href="{up}family/{slug(f)}/"{cur if f == current else ""}>{e(f)} · {len(families[f])}</a></li>'
        for f in fam_order) + '</ul>'

def table(up, rs, show_fam=True):
    head = '<tr><th>A-number</th><th>Name</th>' + ('<th class="fam">Family</th>' if show_fam else '') + '<th class="num">Level</th></tr>'
    out = []
    for r in rs:
        dec = int(r['decomposable'])
        out.append(f'<tr><td><a href="{up}seq/{r["anumber"]}/">{r["anumber"]}</a></td><td>{e(r["name"])}</td>'
                   + (f'<td class="fam">{e(r["family"])}</td>' if show_fam else '')
                   + f'<td class="num">{100 * int(r["level"]) / dec if dec else 0:.1f} %</td></tr>')
    return f'<table class="list">\n<thead>{head}</thead>\n<tbody>\n' + '\n'.join(out) + '\n</tbody>\n</table>'

body = f'''{crumbs('../', [('All sequences', None)])}
<h1>All {N} sequences</h1>
<p class="lead">Every sequence in {SITE}, by A-number, with the share of its decomposable terms in the level class (k &gt; L).
Each has its own page with its weight–level plate and counts, and opens in the interactive 3-D viewer.</p>
{chips('../')}
{table('../', rows)}'''
write('seq/index.html', page('../', 'seq/', f'All {N} sequences · {SITE}',
      f'The {N} integer sequences of {SITE}, from the OEIS, each decomposed into weight × level + jump over 100,000 terms: '
      f'names, A-numbers, families and level shares.', body,
      ld=[crumbs_ld([('All sequences', f'{BASE}/seq/')])]))

for f in fam_order:
    rs = families[f]
    body = f'''{crumbs('../../', [('All sequences', '../../seq/'), (f, None)])}
<h1>{e(f[:1].upper() + f[1:])} <span class="brand-dim">· {len(rs)} sequence{'s' if len(rs) != 1 else ''}</span></h1>
<p class="lead">The sequences of the family “{e(f)}”, by A-number, with the share of their decomposable terms in the level class (k &gt; L).</p>
{chips('../../', f)}
{table('../../', rs, show_fam=False)}'''
    write(f'family/{slug(f)}/index.html', page('../../', f'family/{slug(f)}/', f'{f[:1].upper() + f[1:]}: {len(rs)} sequences · {SITE}',
          clip(f'{len(rs)} integer sequences of the family “{f}” decomposed into weight × level + jump: '
               + ', '.join(f'{r["anumber"]} {r["name"]}' for r in rs[:4]) + '…', 300), body,
          ld=[crumbs_ld([('All sequences', f'{BASE}/seq/'), (f, f'{BASE}/family/{slug(f)}/')])]))

# ── how it works: the decomposition explained, with a live example ──────────
def _spf(n):
    p = 2
    while p * p <= n:
        if n % p == 0: return p
        p += 1
    return n
_nat = []
for a in range(3, 18):                         # the naturals 3 ... 17: the weight is the sieve's smallest prime factor
    l = a - 1; k = _spf(l); L = l // k
    cls = '<td class="lv">level</td>' if k > L else '<td class="wt">weight, tie</td>' if k == L else '<td class="wt">weight</td>'
    _nat.append(f'<tr><td>{a}</td><td>{l}</td><td>{k}</td><td>{L}</td><td>{k} × {L} + 1</td>{cls}</tr>')
NAT_TABLE = ('<div class="scroll narrow"><table class="list demo"><thead><tr><th>a</th><th>l = a − 1</th><th>k, smallest prime factor of l</th>'
             '<th>L</th><th>k × L + d</th><th>class</th></tr></thead><tbody>' + ''.join(_nat) + '</tbody></table></div>')
LEARN_JS = r"""/* decompwlj 3D, written by tools/seo.py: the live example of the "How it works" page.
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
"""
lbody = f'''{crumbs('../', [('How it works', None)])}
<h1>How it works: weight × level + jump</h1>
<p class="lead">Take a strictly increasing sequence of whole numbers — the primes, the squares, any sequence of the OEIS.
Each term, together with the next one, is written as a product plus a jump:</p>
<p class="formula"><i>a</i>(n) = <i>k</i>(n) · <i>L</i>(n) + <i>d</i>(n)</p>
<p class="lead"><i>k</i> is the <b class="wt">weight</b>, <i>L</i> the <b class="lv">level</b> and <i>d</i> the <b>jump</b>.
Every term that can be decomposed is decomposed in exactly one way, by three steps.</p>

<h2>The three steps</h2>
<ol class="steps">
  <li><b>The jump.</b> <i>d</i> = a(n+1) − a(n), the distance to the next term.</li>
  <li><b>What is left.</b> <i>l</i> = a(n) − d. The term decomposes only if <i>l</i> &gt; <i>d</i>, that is a(n) &gt; 2d;
    otherwise the jump is too large and the term is skipped.</li>
  <li><b>Weight and level.</b> The weight <i>k</i> is the <em>smallest divisor of l that is larger than d</em>, and the
    level is <i>L</i> = l / k. Then a(n) = k·L + d, and k &gt; d: the weight always exceeds the jump.</li>
</ol>

<h2>Three primes, worked by hand</h2>
<div class="cards3">
  <div class="card3"><p class="big">113</p><p>next prime 127, so <i>d</i> = 14<br><i>l</i> = 113 − 14 = 99<br>
    divisors of 99: 1, 3, 9, 11, <b>33</b>, 99<br>the smallest above 14 is <i>k</i> = 33, and <i>L</i> = 99 / 33 = 3</p>
    <p class="res">113 = 33 × 3 + 14 · <span class="lv">level class, k &gt; L</span></p></div>
  <div class="card3"><p class="big">1009</p><p>next prime 1013, so <i>d</i> = 4<br><i>l</i> = 1009 − 4 = 1005<br>
    divisors of 1005: 1, 3, <b>5</b>, 15, 67, …<br>the smallest above 4 is <i>k</i> = 5, and <i>L</i> = 1005 / 5 = 201</p>
    <p class="res">1009 = 5 × 201 + 4 · <span class="wt">weight class, k ≤ L</span></p></div>
  <div class="card3"><p class="big">7</p><p>next prime 11, so <i>d</i> = 4<br><i>l</i> = 7 − 4 = 3<br>
    3 is not larger than 4: no divisor of 3 can exceed the jump</p>
    <p class="res">7 does not decompose (7 ≤ 2 × 4)</p></div>
</div>
<p class="lead">A term is in the <b class="lv">level class</b> when k &gt; L and in the <b class="wt">weight class</b> when
k ≤ L; the ties k = L sit on the diagonal. For the primes the weight is never even, and a prime above 3 is the lesser
of a twin pair exactly when its weight is 3.</p>

<h2>A little theory</h2>

<h3>A Euclidean division</h3>
<p class="lead">Since the weight exceeds the jump, 0 ≤ d &lt; k, and a = k·L + d is exactly the
<a href="https://en.wikipedia.org/wiki/Euclidean_division" rel="noopener">Euclidean division</a> of a by k:
the level L is the quotient and the jump d is the remainder. So the weight can be read another way:
<b>k is the smallest number above d for which dividing a by k leaves the remainder d</b>, the distance to the
next term. For 113: dividing by 15, 16, …, 32 leaves other remainders; 113 = 33 × 3 + 14 is the first division by a
number above 14 whose remainder is 14. It is the same condition as before, since a leaves the remainder d on
division by k exactly when k divides a − d.</p>

<h3>When a term decomposes</h3>
<p class="lead">If a &gt; 2d, then l = a − d &gt; d, and l is itself a divisor of l larger than d: a weight always exists
and is unique, with d &lt; k ≤ l. If a ≤ 2d, then l ≤ d and no divisor of l can exceed d: the term does not decompose.
In other words a term decomposes unless the next term is at least 1.5 times as large, which only fast-growing
sequences do all the time: the powers of 2, or the Fibonacci numbers, whose ratio tends to 1.618.</p>

<h3>Weight or level: the window (d, √l]</h3>
<p class="lead">Because k·L = l, the weight class k ≤ L is the same as k ≤ √l. A term is therefore in the
<b class="wt">weight class</b> exactly when l has a divisor in the window d &lt; k ≤ √l, and in the
<b class="lv">level class</b> when that window holds no divisor of l. Three consequences:</p>
<ul class="steps">
  <li><b>Forced level.</b> When l ≤ d², the window is empty whatever l is made of, and the term is in the level class
    without looking at its divisors. Sequences whose jumps grow faster than the square root of the terms (the
    squares of primes, many polynomials) are almost entirely forced level.</li>
  <li><b>The level line L = 1.</b> L = 1 means k = l: no divisor of l lies strictly between d and l. This always
    happens when l is a prime larger than d.</li>
  <li><b>Ties.</b> k = L means l = k²: the least divisor above d is exactly the square root of l.</li>
</ul>

<h3>The natural numbers and the sieve of Eratosthenes</h3>
<p class="lead">For the natural numbers the jump is always d = 1, so the weight of a is the least divisor of a − 1
above 1: its <b>smallest prime factor</b>. That is what the
<a href="https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes" rel="noopener">sieve of Eratosthenes</a> computes:
it crosses out every composite number first with its smallest prime factor. The weight–level plate of the naturals
is the sieve drawn in the plane: one column k = p for each prime p, holding the terms a for which p is the smallest
prime factor of a − 1.</p>
{NAT_TABLE}
<p class="lead">The window (1, √l] contains a divisor of l exactly when l is composite, so the level class of the
naturals is the set of terms with a − 1 prime, all on the line L = 1: 9,592 of the first 10⁵ terms, one per prime
below 10⁵. The ties are the squares of primes, a − 1 = p² (65 of them). For any other sequence the jump plays the
role of a starting point: the divisors up to d are skipped, as if the sieve began at d + 1.</p>

<h2>Try it</h2>
<p class="lead">Type or paste increasing terms, or start from a sequence below. The table and the plate follow as you type.</p>
<div class="presets">
  <button type="button" class="cta alt" data-preset="primes">Primes</button><button type="button" class="cta alt" data-preset="squares">Squares</button><button type="button" class="cta alt" data-preset="triangular">Triangular numbers</button><button type="button" class="cta alt" data-preset="odd">Odd numbers</button><button type="button" class="cta alt" data-preset="fibonacci">Fibonacci</button>
</div>
<label class="lab" for="terms">Terms, separated by commas or spaces</label>
<textarea id="terms" rows="3" spellcheck="false"></textarea>
<p id="err" class="err" role="alert"></p>
<div class="try">
  <figure class="tryplate"><svg id="plate" viewBox="0 0 300 300" role="img" aria-label="Weight–level plate of the example"></svg>
    <figcaption>The example's plate: weight k across, level L up, on log scales; the dashed line is k = L.
    Hover a point for its values.</figcaption></figure>
  <div class="trytab"><p id="sum" class="note" aria-live="polite"></p>
    <div class="scroll"><table class="list demo"><thead><tr><th>n</th><th>a(n)</th><th>d</th><th>l = a − d</th><th>k</th><th>L</th><th>k × L + d</th><th>class</th></tr></thead>
    <tbody id="rows"></tbody></table></div></div>
</div>
<p class="note">Fibonacci is a sequence where nothing decomposes: each term is less than twice the jump to the next.</p>

<h2>Reading the site</h2>
<ul class="steps">
  <li><b>The plates</b> of the gallery put every decomposable term at (log k, log L), scaled by the largest term:
    <span class="wt">blue</span> for the weight class, <span class="lv">orange</span> for the level class. The dashed
    diagonal is k = L, the edge k·L = a is where the jump is small against the term.</li>
  <li><b>The 3-D viewer</b> adds the jump d as a third axis. Its flat views are the plate (k–L), and the jump against
    the weight (k–d) or the level (L–d); the level line L = 1 holds the terms whose l is itself the weight.</li>
  <li>Each sequence page gives the counts, a note on its shape, and a CSV of all 10⁵ terms with their weight, level
    and jump. Good places to start: the <a href="../seq/A000040/">primes</a>, the
    <a href="../seq/A000027/">natural numbers</a> (d = 1: the weight is the least divisor of a − 1 above 1),
    the <a href="../seq/A000290/">squares</a>.</li>
</ul>
<p class="lead">The decomposition is described in <a href="https://arxiv.org/abs/0711.0865" rel="noopener">arXiv:0711.0865</a>
and on <a href="https://decompwlj.com" rel="noopener">decompwlj.com</a>.</p>
<script src="learn.js" defer></script>'''
write('learn/index.html', page('../', 'learn/', 'How it works: weight × level + jump · decompwlj 3D',
      'The decomposition a(n) = k·L + d of an integer sequence, step by step: the jump, the weight and the level, three worked '
      'primes, and a live example where you type your own terms.', lbody,
      ld=[crumbs_ld([('How it works', f'{BASE}/learn/')])]))
write('learn/learn.js', LEARN_JS)

# ── 404: short URLs such as /A000040 go to their page ───────────────────────
write('404.html', f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Page not found · {SITE}</title>
<meta name="robots" content="noindex">
{icons('/')}
<link rel="stylesheet" href="/css/pages.css">
{GOATCOUNTER}
<script>
  /* /A000040, /a000040/ or /seq/a000040 -> /seq/A000040/ ; the site may live under a sub-path */
  (function () {{
    var m = location.pathname.match(/^(.*?)\\/(?:seq\\/)?(a\\d{{6}})\\/?$/i);
    if (m) location.replace(m[1] + '/seq/' + m[2].toUpperCase() + '/');
  }})();
</script>
</head>
<body>
<header class="top"><a class="brand" href="/"><span class="brand-mark" aria-hidden="true"></span>decompwlj<span class="brand-dim"> 3D</span></a>
  <nav aria-label="Site"><a href="/">Gallery</a><a href="/seq/">All sequences</a></nav></header>
<main>
<h1>Page not found</h1>
<p class="lead">There is nothing at this address. Try the <a href="/">gallery</a> or the <a href="/seq/">list of all {N} sequences</a>;
a sequence also opens by its A-number, as in <a href="/seq/A000040/">/seq/A000040/</a>.</p>
</main>
</body>
</html>
''')

# ── sitemap and robots ──────────────────────────────────────────────────────
urls = [(f'{BASE}/', f'{BASE}/og.png'), (f'{BASE}/learn/', None), (f'{BASE}/seq/', None)]
urls += [(f'{BASE}/family/{slug(f)}/', None) for f in fam_order]
urls += [(f'{BASE}/seq/{r["anumber"]}/', f'{BASE}/thumbs/{r["id"]}.webp') for r in rows]
sm = ['<?xml version="1.0" encoding="UTF-8"?>',
      '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9" xmlns:image="http://www.google.com/schemas/sitemap-image/1.1">']
for u, im in urls:
    sm.append(f'<url><loc>{e(u)}</loc>' + (f'<image:image><image:loc>{e(im)}</image:loc></image:image>' if im else '') + '</url>')
sm.append('</urlset>')
write('sitemap.xml', '\n'.join(sm) + '\n')
write('seq/csv.js', CSV_JS)
write('robots.txt', f'# {SITE}: everything may be crawled.\nUser-agent: *\nAllow: /\n\nSitemap: {BASE}/sitemap.xml\n')

# ── social preview: a mosaic of plates ──────────────────────────────────────
try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:
    print('Pillow missing: og.png not written')
else:
    # the link previews always use the light palette of the pages (css/pages.css), whatever the
    # theme a visitor picks in the app: they are shown on social sites and in chats, mostly light
    BG, CARD, FG, FG2, FG3 = (246, 247, 250), (255, 255, 255), (17, 21, 31), (76, 86, 112), (107, 118, 141)
    ACCENT, LEVEL = (47, 95, 208), (226, 102, 26)
    W, H, T = 1200, 630, 150
    img = Image.new('RGB', (W, H), BG)
    cols, rws = W // T + 1, H // T + 1
    step = max(1, N // (cols * rws))
    pick = [rows[(j * step) % N] for j in range(cols * rws)]
    for j, r in enumerate(pick):
        p = ROOT / 'thumbs' / f'{r["id"]}.webp'
        if not p.exists(): continue
        th = Image.open(p).convert('RGBA').resize((T - 6, T - 6), Image.LANCZOS)
        tile = Image.new('RGBA', (T - 6, T - 6), CARD + (255,)); tile.alpha_composite(th)
        img.paste(tile.convert('RGB'), ((j % cols) * T + 3 - 30, (j // cols) * T + 3 - 30))
    shade = Image.new('RGBA', (W, H), (0, 0, 0, 0)); d = ImageDraw.Draw(shade)
    d.rectangle([0, 190, W, 440], fill=BG + (232,))
    img = Image.alpha_composite(img.convert('RGBA'), shade)
    d = ImageDraw.Draw(img)
    def font(bold, size):
        for f in (f'/usr/share/fonts/truetype/dejavu/DejaVuSans{"-Bold" if bold else ""}.ttf',
                  f'/Library/Fonts/Arial{" Bold" if bold else ""}.ttf', 'C:/Windows/Fonts/arial.ttf'):
            try: return ImageFont.truetype(f, size)
            except OSError: pass
        return ImageFont.load_default()
    d.text((W // 2, 262), 'decompwlj 3D', font=font(True, 72), fill=FG, anchor='mm')
    d.text((W // 2, 345), f'weight × level + jump · {N} integer sequences from the OEIS', font=font(False, 32),
           fill=FG2, anchor='mm')
    d.text((W // 2, 398), 'a(n) = k·L + d', font=font(False, 28), fill=LEVEL, anchor='mm')
    img.convert('RGB').save(ROOT / 'og.png', optimize=True)
    written.append('og.png')

    # one share image per sequence (share/<id>.jpg, 600 x 315): the plate and the A-number; the name and
    # the counts are in the card's title and description, since text costs as much as the plate in a
    # JPEG.  Rewritten only when its bytes change; the 5000 take about 34 MB.
    import io
    SW, SH, P = 600, 315, 271
    f_an, f_sm = font(True, 20), font(False, 14)
    nshare = 0
    for r in rows:
        src = ROOT / 'thumbs' / f'{r["id"]}.webp'
        if not src.exists(): continue
        im = Image.new('RGB', (SW, SH), CARD); d = ImageDraw.Draw(im)
        x0, y0 = (SW - P) // 2, (SH - P) // 2
        th = Image.open(src).convert('RGBA').resize((P, P), Image.LANCZOS)
        im.paste(th, (x0, y0), th)
        d.text((18, 16), r['anumber'], font=f_an, fill=ACCENT)
        d.text((SW - 18, SH - 14), 'decompwlj 3D · k·L + d', font=f_sm, fill=FG3, anchor='rs')
        b = io.BytesIO(); im.save(b, 'JPEG', quality=60, optimize=True, progressive=True)
        out = ROOT / 'share' / f'{r["id"]}.jpg'; out.parent.mkdir(exist_ok=True)
        if not out.exists() or out.read_bytes() != b.getvalue(): out.write_bytes(b.getvalue()); nshare += 1
    print(f'{nshare} share images written or updated')

print(f'{len(written)} files written for {BASE}: {N} sequence pages, {len(fam_order)} family pages, '
      f'the index, 404.html, sitemap.xml ({len(urls)} URLs), robots.txt')
