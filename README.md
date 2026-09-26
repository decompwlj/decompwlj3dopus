# decompwlj 3D — static site

The gallery and the interactive 3-D atlas of the decomposition into weight × level + jump:
four hundred integer sequences, 10⁵ terms each (A007088: 65,535, the most that fit below 2⁵³).
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

Size on disk is 345 MB. With gzip on, a visitor downloads about 0.2 MB for the primes and 1.1 MB
for the stella octangula numbers, the largest sequence; the 400 gallery plates total 6.4 MB and load as you scroll. Nothing is loaded from
another site.

## Links

- `…/3D/` or `…/3D/#home` opens the gallery.
- `…/3D/#A000040` (any A-number in the atlas) opens that sequence too.
- `…/3D/#primes` opens a sequence. The full form, `#primes.kLd.iso.solid.one`, also sets the axes
  (always kLd now; other letters in older links are ignored), the view (iso, xy, xz, yz, edge),
  the point mode (solid, density) and the L = 1 highlight. The viewer's **copy link** button copies exactly that.

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
./decompwlj_gen raw                 # all four hundred, raw a,d,k,L chunks + catalog.csv  (about 5 min)
./decompwlj_gen raw primes 200000   # or one sequence, at any size below 2^53
python3 fetch_oeis.py --catalog raw # OEIS names, offsets and first terms -> oeis.json (already there)
python3 names.py raw                # the OEIS names into the catalogue
python3 audit.py raw                # independent check of every row, and of the OEIS terms
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

The last hundred (also 26 Sep 2026) was checked the same way: `audit.py` passes on all
19,965,535 rows of the 200 sequences (80,000 weights re-derived), the 12 new sequences also in
WLJ Atlas agree with it row for row, and the first hundred and their plates rebuild byte for byte.

For the last two hundred, the OEIS records became reachable through the OEIS's own data
repository on GitHub (github.com/oeis/oeisdata). `fetch_oeis.py` keeps each sequence's name,
offset and first terms in `tools/oeis.json`, and `audit.py` now checks the first rows of all
four hundred against those real terms (up to 40 per sequence, from the offset). That check found
three wrong offsets among the earlier sequences, now corrected: A000069, A001969 and A001082
start at n = 1, not n = 0. Their terms were right. `audit.py` passes on all 39,965,535 rows
(160,000 weights re-derived by exhaustive search).

## The four hundred sequences

The page shows each sequence under its OEIS name, and its A-number links to the OEIS entry.
`data/catalog.csv` lists all four hundred.

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

The last hundred:

| family | sequences |
|---|---|
| polynomial | n(n+3)/2 A000096, n(n+2) A005563, 2n² A001105, odd squares A016754, generalized octagonal A001082, centered pentagonal A005891, centered octahedral A001845, centered cube A005898, stella octangula A007588, sums of odd squares A000447, hexagonal pyramidal A002412, cake A000125, magic constants A006003, second hexagonal A014105, second pentagonal A005449, triangular matchstick A045943, ⌊n²/3⌋ A000212, ⌊n²/2⌋ A007590, ⌊n^(3/2)⌋ A000093 |
| primes | 1, 3, 5, 7 mod 8 A007519–A007522, ending in 1, 3, 7, 9 A030430–A030433, balanced A006562, strong A051634, weak A051635, 2p − 1 prime A005382, (p + 1)/2 prime A005383, p + 8 prime A023202, gap of 6 A031924, triplets A022004 and A022005, p(n)p(n+1) A006094, n·p(n) A033286, p(n) + n A014688, p(n) + p(n+1) A001043, 2p A100484, 3p A001748 |
| multiplicative | pq² A054753, four distinct primes A046386, 5-almost primes A014614, cubefree A004709, non-cubefree A046099, ω = 2, 3, 4 A007774 A033992 A033993, μ = 1 A030229, μ = −1 A030059, Ω even A028260 |
| divisor sum, totient | arithmetic numbers A003601, totients A002202, nontotients A007617 |
| digit rule | bases 4, 5, 6 read in decimal A007090–A007092, no 0 in base 4 A023705, no 0 in base 3 A032924, with a 1 A011531, with a 9 A011539, no 1 A052383, no 2 A052404, only even digits A014263, nondecreasing digits A009994, nonincreasing digits A009996, prime digit sum A028834, even digit sum A054683, unhappy A031177, non-Harshad A065877, non-palindromes A029742 |
| binary rule | balanced binary A031443, odd trailing zeros A036554 |
| Beatty | ⌊n(1+√2)⌋ A003151, ⌊n(1+1/√2)⌋ A003152, ⌊n(2+√2)⌋ A001952, ⌊n√5⌋ ⌊n√6⌋ ⌊n√7⌋ ⌊n√8⌋ A022839–A022842 |
| quadratic form | two distinct nonzero squares A004431, three nonzero squares A000408, not two squares A022544 |
| complement, powers | non-triangular A014132, non-cubes A007412, non-Fibonacci A001690, perfect powers A001597, non-perfect powers A007916 |
| other | figure-figure A005228, not divisible by 5 A047201, 4n 5n 6n 7n A008586–A008589, 3n+1 A016777, 3n+2 A016789, 4n+1 A016813, 4n+3 A004767 |

The last two hundred:

| family | sequences |
|---|---|
| polynomial | A001106, A001107, A051682, A051624, A051865, A051866, A051867, A051868, A060544, A062786, A069125, A003154, A069099, A001539, A002943, A033991, A007742, A002939, A002061, A058331, A056220, A014206, A028387, A008865, A006331, A002413, A002414, A007584, A007585, A005894, A006566, A006564, A005917, A001093, A068601, A034262, A005915, A005902, A005906, A033568 |
| Beatty | A022843, A022844, A054385, A003512, A001953, A001954, A003231, A022342, A003622, A022846, A022847, A022848, A001961 |
| primes | A068228, A040117, A068229, A068231, A045468, A003631, A038873, A003629, A001132, A033200, A045707, A023204, A023208, A023212, A109611, A007529, A001122, A001913, A031925, A031926, A031928, A031930, A038580, A007821, A106856, A033205, A001749, A036689, A036690, A049001, A034961, A034963, A024675, A040976, A005097, A006254, A002822, A014091, A014092 |
| multiplicative | A046306, A046308, A046310, A046315, A046388, A000977, A000379, A046100, A046101, A007674, A016105, A062503, A036668, A078972, A120944 |
| totient | A003277, A005277 |
| forced divisor | A008364, A008365, A008366 |
| smooth | A051038, A080197, A080681 |
| divisor count | A030513, A030515, A030626, A030628, A030630 |
| digit rule | A046758, A046760, A006753, A019506, A011532, A011533, A011534, A011535, A011536, A011537, A011538, A052405, A052406, A052413, A052414, A052419, A052421, A007093, A007094, A007095, A054684, A028840, A014190, A014192, A029952, A029953, A023717, A030141, A084984, A081605 |
| quadratic form | A028982, A000419, A002481, A002480, A000401, A020669, A003325, A020756, A020757 |
| divisor sum | A028983 |
| binary rule | A004780, A003754, A022155, A203463, A091072, A091067, A010061 |
| Fibonacci | A020899 |
| sieve | A003309, A045954 |
| self-referential | A030124, A002859, A002977, A094222, A005236 |
| summatory | A000788, A037123, A013939, A022559, A005187, A006046, A064608 |
| residue class | A047203, A047209, A047220, A047229, A047238, A047246, A047255, A047261, A047266, A047273, A045572, A160545, A014601, A042963, A047211, A007494, A032766, A047212 |

## Credits

The decomposition: Rémi Eismann, arXiv:0711.0865, https://decompwlj.com.
three.js © 2010–2024 three.js authors, MIT licence (`vendor/LICENSE-three.js.txt`).
