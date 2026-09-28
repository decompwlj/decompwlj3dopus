# decompwlj 3D

An interactive atlas of Rémi Eismann's **decomposition into weight × level + jump**, built as a
static website: 2000 integer sequences from the OEIS, 10⁵ terms each, shown as a browsable gallery
of 2-D plates and an explorable 3-D point cloud (three.js / WebGL).

**Live site:** <https://decompwlj.net/> (also <https://decompwlj.github.io/decompwlj3dopus/>)

Every term of a strictly increasing sequence is written

    a(n) = k(n) · L(n) + d(n)

where the **jump** d(n) = a(n+1) − a(n), the **weight** k(n) is the least divisor of a(n) − d(n)
greater than d(n), and the **level** L(n) = (a(n) − d(n)) / k(n). A term is level-classified when
k > L and weight-classified otherwise. Background: [decompwlj.com](https://decompwlj.com) and
[arXiv:0711.0865](https://arxiv.org/abs/0711.0865).

The site is plain files: one HTML page, compressed data files and WebP previews. There is no server
code, database, build step or third-party request, so it runs on any static host.

---

## Contents

- [Features](#features)
- [Quick start](#quick-start)
- [Deploying](#deploying)
  - [GitHub Pages](#github-pages)
  - [A custom domain with HTTPS](#a-custom-domain-with-https)
  - [Any other static host](#any-other-static-host)
  - [Search engines and link previews](#search-engines-and-link-previews)
  - [Updating a deployment](#updating-a-deployment)
- [Requirements](#requirements)
- [Repository layout](#repository-layout)
- [Using the site](#using-the-site)
- [Data format](#data-format)
- [Rebuilding the data](#rebuilding-the-data)
- [Troubleshooting](#troubleshooting)
- [Credits and licences](#credits-and-licences)

---

## Features

- **Gallery** of 2000 weight–level plates, sorted by A-number, with full-text search, a
  **family selector** (primes, polynomial, digit rule, Beatty, sieve, …) and a
  **Random sequence** button.
- **Progressive loading**: the gallery builds cards thirty at a time as you scroll, and each preview
  image is requested only when its card nears the screen.
- **3-D viewer** of log k, log L, log d with preset views (3D, k–L, k–d, L–d, edge), point size and
  opacity, a density mode, the L = 1 highlight, filters on n and d, and an animated sweep along n.
- **OEIS integration**: each sequence is shown under its official OEIS name, and every A-number
  links to its OEIS entry. `#A000040`-style links open a sequence directly.
- **Shareable views**: the URL fragment records the sequence, view and point mode.
- Light and dark themes, keyboard shortcuts, PNG snapshots.

## Quick start

To preview the site on your own computer, serve the repository root over HTTP:

```sh
git clone https://github.com/decompwlj/decompwlj3dopus.git
cd decompwlj3dopus
python3 -m http.server 8000
```

Then open <http://localhost:8000/>. Any static server works (`npx serve`, `php -S`, …).

> Opening `index.html` directly from disk (`file://…`) does not work: browsers block ES modules
> and `fetch()` for local files. The page detects this and says so.

## Deploying

The deployable site is the repository root: `index.html`, `data/`, `thumbs/`, `vendor/` and the
empty `.nojekyll` file. Paths are relative, so the site works at a domain root or in any
subfolder.

### GitHub Pages

1. Push the repository to GitHub (or fork this one).
2. Open **Settings → Pages**.
3. Under **Build and deployment**, set **Source** to **Deploy from a branch**, choose **main** and
   **/ (root)**, and click **Save**.
4. After one or two minutes the site is live at `https://<user>.github.io/<repository>/`. The
   **Actions** tab shows each publication.

Every push to `main` republishes the site. `.nojekyll` tells GitHub Pages to serve the files as
they are, without running Jekyll.

The site is about 450 MB, within GitHub Pages' limits: 1 GB per published site, with a soft
bandwidth limit of 100 GB per month.

### A custom domain with HTTPS

This repository is configured for **decompwlj.net** through the `CNAME` file. To use your own
domain, replace the name in `CNAME` (or set it under **Settings → Pages → Custom domain**, which
rewrites the file), then create these DNS records at your registrar:

| Type | Name | Value |
|---|---|---|
| A | @ (apex) | `185.199.108.153` |
| A | @ | `185.199.109.153` |
| A | @ | `185.199.110.153` |
| A | @ | `185.199.111.153` |
| AAAA | @ | `2606:50c0:8000::153` |
| AAAA | @ | `2606:50c0:8001::153` |
| AAAA | @ | `2606:50c0:8002::153` |
| AAAA | @ | `2606:50c0:8003::153` |
| CNAME | www | `<user>.github.io.` |

Also:

- Remove any other A or AAAA records on the apex. Registrars often add a default parking or hosting
  record, and also a web "redirection", which prevents HTTPS.
- If the zone has a CAA record, it must allow `letsencrypt.org`.

Once the DNS check in **Settings → Pages** passes, GitHub requests a Let's Encrypt certificate.
This usually takes minutes, sometimes up to 24 hours. Then tick **Enforce HTTPS**.

It is also recommended to verify the domain under your account's **Settings → Pages**. This adds a
TXT record and prevents anyone else from using the domain with GitHub Pages.

### Any other static host

Upload the repository root, or at least `index.html`, `data/`, `thumbs/` and `vendor/`, to any
host that serves static files over HTTP(S): Apache, nginx, Caddy, Netlify, Cloudflare Pages, an S3
bucket with a CDN, and so on.

Server requirements:

- Serve `.js` as JavaScript. Standard MIME tables already do.
- Serve `data/seq/**/chunk-NNN.bin.gz` as the raw gzip file (any content type), **or** with
  `Content-Encoding: gzip`. The page handles both.
- Do not rewrite missing files to `index.html`. A missing chunk must return 404.
- No authentication in front of the data files.

Ready-made examples are in [`deploy/`](deploy/): `apache-vhost.conf` covers a subfolder or a
virtual host, and `nginx.conf` a subfolder or a server block. Long cache lifetimes for `data/`,
`thumbs/` and `vendor/` are safe. Chunk files do not change when other sequences are added.

### Search engines and link previews

The gallery and viewer run from one page and a URL fragment (`#A000040`), which search engines do
not index. `tools/seo.py` therefore writes plain HTML pages that they can index, all from
`data/catalog.csv`:

| Path | Contents |
|---|---|
| `seq/A000040/` | One page per sequence: OEIS name, plate, counts, note, links to the 3-D viewer, the OEIS and the neighbouring sequences |
| `seq/` | The list of all sequences, by A-number |
| `family/<family>/` | The sequences of one family |
| `404.html` | Not-found page; short URLs such as `/A000040` redirect to `/seq/A000040/` |
| `sitemap.xml`, `robots.txt` | Every page for the crawlers (with the plates as images); everything may be crawled |
| `og.png` | The 1200 × 630 preview shown when a link is shared |

Every page carries a canonical URL, a description, Open Graph tags and schema.org data
(`Dataset` and `BreadcrumbList`). The base URL comes from `CNAME`; a site without a custom domain
passes it: `python3 seo.py --base https://<user>.github.io/<repo>`. The 404 page's links assume the
site is at the root of its domain.

To get the site indexed:

1. **Google Search Console** (search.google.com/search-console): add a *Domain* property for the
   domain and verify it with the TXT record it gives, added at your DNS provider (for OVH:
   *Web Cloud → Domain names → DNS zone → Add an entry → TXT*). Then open **Sitemaps** and submit
   `https://<domain>/sitemap.xml`.
2. **Bing Webmaster Tools** (bing.com/webmasters): sign in and import the site from Search
   Console, or verify it the same way and submit the sitemap. Bing also feeds DuckDuckGo and
   Yahoo.
3. On GitHub, fill in the repository's **About** box (gear icon on the repository page): a
   description, the website URL and topics such as `oeis`, `integer-sequences`, `number-theory`,
   `prime-numbers`, `mathematics`, `visualization`, `threejs` and `webgl`. Under **Settings →
   General → Social preview**, upload `og.png`.

Indexing takes days to weeks. Search Console's **Pages** report shows progress.

### Updating a deployment

Replace the files and ask visitors to reload (Ctrl+F5, or Cmd+Shift+R on a Mac) if their browser
cached the previous `index.html` or `data/catalog.csv`. For GitHub Pages, merging into `main` is
enough.

## Requirements

**Visitors** need a current browser with WebGL, import maps and `DecompressionStream`:
Chrome/Edge 111+, Firefox 113+ or Safari 16.4+. The gallery works without WebGL; the 3-D viewer
needs it.

**Rebuilding the data** (optional) needs a C compiler (gcc or clang), Python 3.9+ with `numpy` and
`Pillow`, about 2 GB of RAM, 7 GB of free disk space and about 55 minutes of CPU time for the full pipeline.

## Repository layout

| Path | Contents |
|---|---|
| `index.html` | The whole application: gallery, viewer, CSS and JavaScript (~110 kB) |
| `data/catalog.csv` | One row per sequence: id, A-number, OEIS name, family, index range, counts, ranges, note |
| `data/seq/<id>/chunk-000.bin.gz`, `chunk-001.bin.gz` | The sequence data, 50,000 terms per chunk, compact binary, gzip-compressed |
| `thumbs/<id>.webp` | Gallery previews (480 × 480, transparent, lossless WebP) |
| `vendor/` | three.js r169 and OrbitControls, unmodified (MIT licence included) |
| `deploy/` | Example Apache and nginx configurations |
| `tools/` | Data generator, OEIS metadata and verification scripts (not needed at runtime) |
| `docs/SEQUENCES.md` | All 2000 sequences by family, and how the data was verified |
| `CNAME`, `.nojekyll` | GitHub Pages settings: custom domain; serve files as they are |
| `seq/`, `family/`, `404.html` | Static pages for search engines, written by `tools/seo.py` |
| `sitemap.xml`, `robots.txt`, `og.png` | Sitemap, crawler rules and link preview image |
| `css/pages.css` | The static pages' stylesheet |
| `js/`, `css/app.css` | Earlier modular sources, kept for reference; the page does not load them |

## Using the site

- `…/` or `…/#home` opens the gallery. `…/#primes` or `…/#A000040` opens a sequence.
- `…/seq/A000040/` is the sequence's own page (also reached from `…/A000040`), `…/seq/` lists
  all sequences, and `…/family/primes/` lists one family.
- The full fragment `#primes.kLd.iso.solid.one` also sets the view (`iso`, `xy`, `xz`, `yz`,
  `edge`), the point mode (`solid`, `density`) and the L = 1 highlight. **Copy link** in the viewer
  copies it.
- **Gallery keys:** `/` focuses the search; `Home` (or the **Top** button that appears while
  scrolling) returns to the top; `R` opens a random sequence (within the selected family and
  search); `T` switches the theme.
- **Compare** (panel on the right): pick a second sequence by A-number or name, or **random**.
  **overlay** draws both clouds in one scene, the second in violet and green; **side by side**
  shows two views turned by one camera (stacked on a tall screen). Both share one box, so their
  scales match. **swap** exchanges the two, **clear** ends the comparison. The link carries it:
  `#primes.kLd.iso.solid.vs-A001359` (overlay) or `….vs-A001359.split` (side by side).
- **Viewer keys:** `G` or `Esc` returns to the gallery; `↑`/`↓` moves to the neighbouring
  sequence; `R` opens a random one; `1`–`5` switch views; `V` switches a comparison between overlay and side by side; `C` and `P` toggle the side panels; `T`
  switches the theme; `Space` sweeps along n.

## Data format

`data/catalog.csv` is RFC 4180 CSV with a header row. The main columns:

| Column | Meaning |
|---|---|
| `id` | Short identifier, also the folder name under `data/seq/` and the URL fragment |
| `anumber`, `name`, `alias` | OEIS A-number, OEIS name, and an older short name (searchable) |
| `family` | Family label used by the gallery's selector |
| `n0`, `terms`, `chunks`, `chunk_rows` | Index of the first term (the OEIS offset), number of terms, chunk layout |
| `decomposable`, `level`, `weight`, `ties`, `level_one`, `forced` | Classification counts |
| `amin` … `dmax` | Ranges of a, k, L and d |
| `note` | The text shown under the statistics |

Each chunk is a small binary file, gzip-compressed. Every number is an unsigned LEB128 varint
(7 bits per byte, low bits first):

| Part | Contents |
|---|---|
| `dwj1` | 4-byte format tag |
| `n`, `a0` | Number of rows in the chunk, first term |
| `d[0]` … `d[n−1]` | The jumps |
| `s[0]` … `s[n−1]` | The smaller factor of a − d = k·L: `0` if the term does not decompose, `2k` if k ≤ L, `2L + 1` if L < k |

Every term is present, decomposable or not, so row i is index n = n0 + i. The page rebuilds a(n)
as a running sum of the jumps, and the factor not stored as (a − d) divided by the stored one.
Storing the smaller factor, at most √(a − d), makes the files about 40 % smaller than storing k.
The page rejects a chunk if a division is not exact, if k ≤ d, if the stored factor is not the
smaller one, or if a chunk does not continue the previous one. All values are below 2⁵³, so
JavaScript numbers hold them exactly. `tools/compact.py` has an `encode()` and a `decode()` for
reading the files from Python.

A sequence costs from under 1 kB to 0.6 MB to download (median 0.16 MB). The gallery previews
total 13 MB, but only the visible ones are fetched.

## Rebuilding the data

The published data is reproducible. The generator writes identical files on every run, so
unchanged sequences produce no diff.

```sh
cd tools
cc -O2 -o decompwlj_gen decompwlj_gen.c -lm
./decompwlj_gen raw                 # all sequences: raw a,d,k,L chunks + catalog.csv (~25 min)
./decompwlj_gen raw primes 200000   # or a single sequence, at any length below 2^53
python3 names.py raw                # OEIS names into the catalogue (from oeis.json)
python3 audit.py raw                # independent verification of every row and of the OEIS terms
python3 compact.py raw ../data      # the gzip chunks the site loads
python3 thumbs.py  raw ../thumbs    # the gallery previews
python3 seo.py                      # the static pages, sitemap.xml and robots.txt
```

To **add a sequence**:

1. Write a generator function in `tools/decompwlj_gen.c`. It fills `t[0..cnt-1]` with strictly
   increasing terms.
2. Add an entry to the `defs[]` table: id, A-number, name, family, note, number of terms, first
   index, function.
3. Run `python3 fetch_oeis.py A123456` to record its OEIS name, offset and first terms in
   `oeis.json`.
4. Run the pipeline above; `audit.py` must report "all clear".

[`docs/SEQUENCES.md`](docs/SEQUENCES.md) lists every sequence by family and describes the
verification.

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| Blank page, or "Open this page through a web server" | The page was opened from disk. Serve it over HTTP (see [Quick start](#quick-start)). |
| An old version is shown after an update | Browser cache. Reload with Ctrl+F5 (Cmd+Shift+R on a Mac). GitHub Pages can take a few minutes to publish. |
| "Enforce HTTPS" stays greyed out | DNS is not yet pointing at GitHub, an old A record or registrar redirection remains, or a CAA record blocks Let's Encrypt. Fix the records, then remove and re-add the custom domain to retry. |
| A sequence fails to load (404 on a chunk) | The server rewrites missing files, or `data/seq/` was not fully uploaded. Upload the whole `data/` folder. |
| The gallery works but the 3-D view does not | WebGL is disabled or unavailable (some remote desktops and VMs). Enable hardware acceleration in the browser. |

## Credits and licences

- The decomposition into weight × level + jump: Rémi Eismann,
  [arXiv:0711.0865](https://arxiv.org/abs/0711.0865), [decompwlj.com](https://decompwlj.com).
- Sequence names, offsets and reference terms: [The On-Line Encyclopedia of Integer
  Sequences](https://oeis.org), licensed under
  [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/); read from
  [github.com/oeis/oeisdata](https://github.com/oeis/oeisdata).
- [three.js](https://threejs.org) r169, © 2010–2024 three.js authors, MIT licence
  (`vendor/LICENSE-three.js.txt`).
