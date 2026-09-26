# decompwlj 3D — static site

The gallery and the interactive 3-D atlas of the decomposition into weight × level + jump:
fifty integer sequences, 10⁵ terms each (A007088: 65,535, the most that fit below 2⁵³).
It is plain files, with no server code, no database and no build step.

## Put it online

1. Upload this whole folder to your web server, for example as `/var/www/decompwlj-3d`,
   or straight under your document root as `3D/`.
2. Open it over HTTP(S), for example `https://decompwlj.com/3D/`. Opening `index.html` as a
   file (`file://`) does not work: browsers block ES modules and `fetch()` there.
3. Optional: let Apache read `.htaccess` (`AllowOverride All` on the folder) for gzip and
   cache headers; `deploy/apache-vhost.conf` has a subfolder Alias and a separate-vhost form.
   On nginx use `deploy/nginx.conf`. If `.htaccess` ever gives a 500, delete it — nothing in
   it is required.

Every path in the page is relative, so any folder name or prefix works.

## What is in the folder

| path | what |
|---|---|
| `index.html` | the whole page: gallery, 3-D viewer, all CSS and JavaScript (100 kB) |
| `data/catalog.csv` | one row per sequence: counts, ranges, the note shown under the stats |
| `data/seq/<id>/chunk-000.csv`, `chunk-001.csv` | the sequence, 50,000 terms per chunk |
| `thumbs/<id>.png` | the gallery plates |
| `vendor/` | three.js r169 and OrbitControls, unmodified (MIT, licence included) |
| `deploy/` | Apache and nginx configuration |
| `tools/` | the generator and the scripts that rebuild `data/` and `thumbs/` |

Size on disk is 46 MB. With gzip on, a visitor downloads about 0.2 MB for the primes and 1.1 MB
for the cubes, the largest sequence; the gallery alone is under 1 MB. Nothing is loaded from
another site.

## Links

- `…/3D/` or `…/3D/#home` opens the gallery.
- `…/3D/#primes` opens a sequence. The full form, `#primes.kLd.iso.solid.one`, also sets the axes
  (X, Y, Z from k L d a l n), the view (iso, xy, xz, yz, edge), the point mode (solid, density)
  and the L = 1 highlight. The viewer's **copy link** button copies exactly that.

Keys in the viewer: G or Esc back to the gallery, ↑/↓ the next sequence, 1–5 the views,
C and P the side panels, T the theme, Space the sweep along n. In the gallery, / focuses the search.

## Browser support

Any current browser with WebGL and import maps: Chrome and Edge 111+, Firefox 113+, Safari 16.4+.

## Data format

`chunk-NNN.csv` starts with `#a0=<first term>` and `d,k`, then one `d,k` row per term. Every
term is written, decomposable or not (non-decomposable as k = 0), so row i is index n = n0 + i
(n0 is in the catalogue; it follows the OEIS offset). The page rebuilds a(n) as a running sum
of the gaps and L = (a − d)/k, and refuses a chunk if a division is not exact, if k ≤ d, or if
a chunk does not continue the previous one. It also reads the generator's plain `a,d,k,L`
chunks, sniffed on the first byte.

## Regenerate

```
cd tools
cc -O2 -o decompwlj_gen decompwlj_gen.c -lm
./decompwlj_gen raw                 # all fifty, raw a,d,k,L chunks + catalog.csv  (about 40 s)
./decompwlj_gen raw primes 200000   # or one sequence, at any size below 2^53
python3 audit.py raw                # independent check of every row
python3 compact.py raw ../data      # the chunks the page loads
python3 thumbs.py  raw ../thumbs    # the gallery plates (numpy, Pillow)
```

The PARI/GP cross-check (`decompwlj_3d_check.gp`) and the session records are in the decompwlj
project. This build was verified there: 4,965,535 rows byte-identical between the C engine and
PARI/GP 2.15.4, the OEIS offsets and first terms checked, and the 28 Aug 50-sequence census
reproduced exactly.

## Credits

The decomposition: Rémi Eismann, arXiv:0711.0865, https://decompwlj.com.
three.js © 2010–2024 three.js authors, MIT licence (`vendor/LICENSE-three.js.txt`).
