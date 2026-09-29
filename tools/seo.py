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
    family/<family>/index.html    the sequences of one family
    404.html                      not-found page; /A000040 redirects to /seq/A000040/
    sitemap.xml, robots.txt       for the crawlers
    og.png                        1200 x 630 preview for links shared on social sites
                                  (needs Pillow; skipped without it)

Run it again whenever the catalogue changes. The output is deterministic.
"""
import csv, html, json, pathlib, re, sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
SITE = 'decompwlj 3D'

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

families = {}
for r in rows: families.setdefault(r['family'], []).append(r)
fam_order = sorted(families, key=lambda f: (-len(families[f]), f.lower()))

FAVICON = ("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'%3E%3Crect width='32' height='32' rx='6' "
           "fill='%230b0e14'/%3E%3Ccircle cx='10' cy='21' r='3' fill='%234c8dff'/%3E%3Ccircle cx='19' cy='13' r='3' fill='%23ff8c42'/%3E"
           "%3Ccircle cx='25' cy='7' r='2' fill='%23ffe45c'/%3E%3C/svg%3E")
PM, PX = 0.03 / 1.06, 1.03 / 1.06
PLATE_SVG = (f'<svg viewBox="0 0 1 1" preserveAspectRatio="none" aria-hidden="true">'
             f'<rect x="{PM}" y="{PM}" width="{PX - PM}" height="{PX - PM}"/>'
             f'<line class="ed" x1="{PM}" y1="{PM}" x2="{PX}" y2="{PX}"/>'
             f'<line class="dg" x1="{PM}" y1="{PX}" x2="0.5" y2="0.5"/></svg>')

def page(up, path, title, desc, body, image=None, ld=()):
    """A complete page. up: relative path back to the site root ('../../')."""
    url = f'{BASE}/{path}'
    img = image or f'{BASE}/og.png'
    card = 'summary' if image else 'summary_large_image'
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
<link rel="icon" href="{FAVICON}">
<link rel="stylesheet" href="{up}css/pages.css">
<meta property="og:type" content="website">
<meta property="og:site_name" content="{SITE}">
<meta property="og:title" content="{e(title)}">
<meta property="og:description" content="{e(desc)}">
<meta property="og:url" content="{e(url)}">
<meta property="og:image" content="{e(img)}">
<meta name="twitter:card" content="{card}">
{lds}</head>
<body>
<header class="top">
  <a class="brand" href="{up}"><span class="brand-mark" aria-hidden="true"></span>decompwlj<span class="brand-dim"> 3D</span></a>
  <nav aria-label="Site"><a href="{up}">Gallery</a><a href="{up}seq/">All sequences</a><a href="https://decompwlj.com" rel="noopener">decompwlj.com</a></nav>
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
         'and <a href="https://arxiv.org/abs/0711.0865" rel="noopener">arXiv:0711.0865</a>.</p>')

CSV_JS = r"""/* decompwlj 3D, written by tools/seo.py: the "Download CSV" button of a sequence page.
   Rebuilds n;a;weight;level;jump from the data chunks (the dwj2 format of tools/compact.py),
   checking every division as the viewer does; nothing is stored. */
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
  function rows(b, url, out) {
    if (!(b[0] === 0x64 && b[1] === 0x77 && b[2] === 0x6a && (b[3] === 0x31 || b[3] === 0x32))) throw new Error(url + ': not a dwj1/dwj2 chunk');
    var m = b[3] === 0x32 ? b[4] : 0, i = b[3] === 0x32 ? 5 : 4;
    function unzz(z) { return z % 2 === 0 ? z / 2 : -(z + 1) / 2; }
    function rd() { var v = 0, m = 1, c; do { if (i >= b.length) throw new Error(url + ': truncated'); c = b[i++]; v += (c & 127) * m; m *= 128; } while (c & 128); return v; }
    var n = rd(), a = rd(), d = new Array(n);
    for (var r = 0, p = 0, pp = 0; r < n; r++) { var v = rd(); d[r] = m === 0 ? v : m === 1 ? p + unzz(v) : 2 * p - pp + unzz(v); pp = p; p = d[r]; }
    for (r = 0; r < n; r++) {
      var s = rd(), k = 0, L = 0;
      if (s) {
        var h = Math.floor(s / 2), q = (a - d[r]) / h;
        k = s % 2 ? q : h; L = s % 2 ? h : q;
        if (!Number.isInteger(q) || q * h !== a - d[r] || k <= d[r]) throw new Error(url + ': row ' + r + ' does not fit');
      }
      out.push([a, k, L, d[r]]);
      a += d[r];
    }
  }
  btn.addEventListener('click', function () {
    var id = btn.dataset.id, an = btn.dataset.an, n0 = +btn.dataset.n0, nc = +btn.dataset.chunks;
    var label = btn.textContent; btn.disabled = true; btn.textContent = 'Building…';
    var jobs = [];
    for (var c = 0; c < nc; c++) {
      (function (c) {
        var url = '../../data/seq/' + id + '/chunk-' + String(c).padStart(3, '0') + '.bin.gz';
        jobs.push(fetch(url).then(function (res) { if (!res.ok) throw new Error(res.status + ' ' + url); return res.arrayBuffer(); })
          .then(function (buf) { return inflate(new Uint8Array(buf)); })
          .then(function (b) { var out = []; rows(b, url, out); return out; }));
      })(c);
    }
    Promise.all(jobs).then(function (parts) {
      var text = ['n;a;weight;level;jump\n'], buf = '', n = n0;
      parts.forEach(function (p) { p.forEach(function (r) {
        buf += n++ + ';' + r[0] + ';' + (r[1] || '') + ';' + (r[1] ? r[2] : '') + ';' + r[3] + '\n';
        if (buf.length > 65536) { text.push(buf); buf = ''; }
      }); });
      text.push(buf);
      var url = URL.createObjectURL(new Blob(text, { type: 'text/csv' }));
      var a = document.createElement('a'); a.href = url; a.download = 'decompwlj_' + an + '.csv';
      document.body.appendChild(a); a.click(); a.remove();
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
<p class="meta"><a class="an" href="https://oeis.org/{A}" rel="noopener">{A}</a> on the OEIS · family <a href="{up}family/{slug(fam)}/">{e(fam)}</a>{f' · also known as {e(r["alias"])}' if r['alias'] else ''}</p>
<div class="seq">
  <figure>
    <div class="plate">{PLATE_SVG}<img src="{up}thumbs/{r['id']}.webp" width="480" height="480" alt="Weight–level plate of {e(name)}"><span class="axl k" aria-hidden="true">k</span><span class="axl L" aria-hidden="true">L</span></div>
    <figcaption>Weight k across, level L up, both on log scales: blue in the weight class (k ≤ L), orange in the level class (k &gt; L). The dashed diagonal is k = L.</figcaption>
  </figure>
  <div>
    <a class="cta" href="{up}#{A}">Open in the 3-D viewer</a><a class="cta alt" href="https://oeis.org/{A}" rel="noopener">{A} on the OEIS</a><button type="button" class="cta alt" id="csv" data-id="{r['id']}" data-an="{A}" data-n0="{n0}" data-terms="{terms}" data-chunks="{r['chunks']}" title="n;a;weight;level;jump, one row per term, rebuilt from the data">Download CSV</button>
    <table class="stats">
{''.join(f'      <tr><th scope="row">{k}</th><td>{v}</td></tr>{chr(10)}' for k, v in stats)}    </table>
    <p class="note">{e(r['note'])}</p>
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
    write(f'seq/{A}/index.html', page(up, f'seq/{A}/', title, desc, body, image=f'{BASE}/thumbs/{r["id"]}.webp', ld=ld))

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

# ── 404: short URLs such as /A000040 go to their page ───────────────────────
write('404.html', f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Page not found · {SITE}</title>
<meta name="robots" content="noindex">
<link rel="icon" href="{FAVICON}">
<link rel="stylesheet" href="/css/pages.css">
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
urls = [(f'{BASE}/', f'{BASE}/og.png'), (f'{BASE}/seq/', None)]
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
    W, H, T = 1200, 630, 150
    img = Image.new('RGB', (W, H), (11, 14, 20))
    cols, rws = W // T + 1, H // T + 1
    step = max(1, N // (cols * rws))
    pick = [rows[(j * step) % N] for j in range(cols * rws)]
    for j, r in enumerate(pick):
        p = ROOT / 'thumbs' / f'{r["id"]}.webp'
        if not p.exists(): continue
        th = Image.open(p).convert('RGBA').resize((T - 6, T - 6), Image.LANCZOS)
        tile = Image.new('RGBA', (T - 6, T - 6), (17, 21, 31, 255)); tile.alpha_composite(th)
        img.paste(tile.convert('RGB'), ((j % cols) * T + 3 - 30, (j // cols) * T + 3 - 30))
    shade = Image.new('RGBA', (W, H), (0, 0, 0, 0)); d = ImageDraw.Draw(shade)
    d.rectangle([0, 190, W, 440], fill=(11, 14, 20, 225))
    img = Image.alpha_composite(img.convert('RGBA'), shade)
    d = ImageDraw.Draw(img)
    def font(bold, size):
        for f in (f'/usr/share/fonts/truetype/dejavu/DejaVuSans{"-Bold" if bold else ""}.ttf',
                  f'/Library/Fonts/Arial{" Bold" if bold else ""}.ttf', 'C:/Windows/Fonts/arial.ttf'):
            try: return ImageFont.truetype(f, size)
            except OSError: pass
        return ImageFont.load_default()
    d.text((W // 2, 262), 'decompwlj 3D', font=font(True, 72), fill=(230, 235, 245), anchor='mm')
    d.text((W // 2, 345), f'weight × level + jump · {N} integer sequences from the OEIS', font=font(False, 32),
           fill=(154, 166, 189), anchor='mm')
    d.text((W // 2, 398), 'a(n) = k·L + d', font=font(False, 28), fill=(255, 140, 66), anchor='mm')
    img.convert('RGB').save(ROOT / 'og.png', optimize=True)
    written.append('og.png')

print(f'{len(written)} files written for {BASE}: {N} sequence pages, {len(fam_order)} family pages, '
      f'the index, 404.html, sitemap.xml ({len(urls)} URLs), robots.txt')
