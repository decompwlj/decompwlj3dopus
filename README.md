# decompwlj 3D — static site

The gallery and the interactive 3-D atlas of the decomposition into weight × level + jump:
one hundred integer sequences, 10⁵ terms each (A007088: 65,535, the most that fit below 2⁵³).
It is plain files, with no server code, no database and no build step.

## Put it online with GitHub Pages

The site is published straight from this repository, at
**https://decompwlj.github.io/decompwlj3dopus/**.

1. On GitHub, open the repository, then **Settings → Pages** (in the left column).
2. Under **Build and deployment**, set **Source** to **Deploy from a branch**, then choose
   the branch **main** and the folder **/ (root)**, and press **Save**.
3. Wait one or two minutes. The address appears at the top of the Pages settings;
   the **Actions** tab shows the publication running (a green tick when it is done).

Every later commit to `main` republishes the site by itself. The empty file `.nojekyll` tells
GitHub to serve the files as they are, without running its Jekyll site builder.

To look at the site on your own computer first, run `python3 -m http.server 8000` in this
folder and open `http://localhost:8000/`. Opening `index.html` as a file (`file://`) does not
work: browsers block ES modules and `fetch()` there.

## Put it on your own server

1. Upload this whole folder to your web server, for example as `/var/www/decompwlj-3d`,
   or straight under your document root as `3D/`.
2. Open it over HTTP(S), for example `https://decompwlj.com/3D/`.
3. Optional: `deploy/apache-vhost.conf` has a subfolder Alias and a separate-vhost form for
   Apache; on nginx use `deploy/nginx.conf`.

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
| `.nojekyll` | tells GitHub Pages to serve the files as they are |

Size on disk is 86 MB. With gzip on, a visitor downloads about 0.2 MB for the primes and 1.1 MB
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
./decompwlj_gen raw                 # all one hundred, raw a,d,k,L chunks + catalog.csv  (about 90 s)
./decompwlj_gen raw primes 200000   # or one sequence, at any size below 2^53
python3 audit.py raw                # independent check of every row
python3 compact.py raw ../data      # the chunks the page loads
python3 thumbs.py  raw ../thumbs    # the gallery plates (numpy, Pillow)
```

The PARI/GP cross-check (`decompwlj_3d_check.gp`) and the session records are in the decompwlj
project. The first fifty were verified there: 4,965,535 rows byte-identical between the C engine and
PARI/GP 2.15.4, the OEIS offsets and first terms checked, and the 28 Aug 50-sequence census
reproduced exactly.

The second fifty (26 Sep 2026) were checked three ways. `audit.py` passes on all 9,965,535
rows, re-deriving 40,000 weights by exhaustive search. The 24 sequences shared with the
separate WLJ Atlas project (`decompwlj/wlj-atlas`) agree with its data row for row (a, k, L, d).
Rebuilding the first fifty reproduces the published chunks, the catalogue and the 50 gallery
plates byte for byte. `audit.py` compares the new fifty's first terms with values written out
from each definition and OEIS offset, because oeis.org could not be reached from the build
machine.

## The one hundred sequences

The first fifty: A000027, A000037, A000040, A000201, A000217, A000290, A000292, A000326,
A000330, A000384, A000578, A000959, A000960, A000961, A001248, A001358, A001359, A001481,
A001651, A001694, A001855, A002088, A002113, A002145, A002378, A002620, A002808, A002858,
A003052, A003714, A004202, A004207, A005101, A005117, A005153, A005214, A005349, A005408,
A005843, A006218, A006446, A006512, A006881, A007088, A007504, A008864, A024619, A024916,
A026424, A052382.

The second fifty:

| family | sequences |
|---|---|
| binary rule | odious A000069, evil A001969, Moser–de Bruijn A000695, even trailing zeros A003159, binary palindromes A006995 |
| polynomial | heptagonal A000566, octagonal A000567, centered square A001844, centered hexagonal A003215, centered triangular A005448, lazy caterer A000124, n² + 1 A002522, generalized pentagonal A001318, pentagonal pyramidal A002411, octahedral A005900 |
| Beatty | upper Wythoff A001950, ⌊n√2⌋ A001951, ⌊n√3⌋ A022838 |
| primes | 1 mod 4 A002144, 1 mod 6 A002476, 5 mod 6 A007528, Sophie Germain A005384, safe A005385, prime-indexed A006450, isolated A007510, twin A001097, cousin A023200, sexy A023201, emirps A006567, p − 1 A006093 |
| multiplicative | 3-almost primes A014612, 4-almost primes A014613, sphenic A007304, non-squarefree A013929, weak A052485 |
| divisor sum | deficient A005100 |
| quadratic form | Loeschian A003136, x² + 2y² A002479, two nonzero squares A000404, three squares A000378, not three squares A004215 |
| forced divisor | coprime to 30 A007775, coprime to 6 A007310, multiples of 3 A008585 |
| complement | odd nonprimes A014076 |
| digit rule | happy A007770, with a digit 0 A011540, only odd digits A014261, base 3 read in decimal A007089, no digit 2 in base 3 A005836 |

## Credits

The decomposition: Rémi Eismann, arXiv:0711.0865, https://decompwlj.com.
three.js © 2010–2024 three.js authors, MIT licence (`vendor/LICENSE-three.js.txt`).
