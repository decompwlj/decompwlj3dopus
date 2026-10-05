three.js r169 (0.169.0), MIT licence — see LICENSE-three.js.txt.

  three.module.min.js                 the minified ESM build (three/build/three.module.min.js)
  addons/controls/OrbitControls.js    three/examples/jsm/controls/OrbitControls.js

Vendored rather than loaded from a CDN, so the app works offline, behind a proxy and from a
plain copy of the files. The import map at the top of index.html resolves "three" and
"three/addons/" to these files; nothing else needs configuring.

To upgrade: npm pack three@<version> and replace these two files.

mp4-muxer 5.2.2 (Vanilagy/mp4-muxer), MIT licence — see LICENSE-mp4-muxer.txt.

  mp4-muxer.min.mjs                   build/mp4-muxer.mjs, minified with esbuild

Loaded only when a video is recorded: it packs the browser's own H.264 frames (WebCodecs)
into the .mp4 file. To upgrade: npm pack mp4-muxer@<version>, then
  npx esbuild package/build/mp4-muxer.mjs --minify --format=esm --outfile=vendor/mp4-muxer.min.mjs
