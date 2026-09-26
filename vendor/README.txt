three.js r169 (0.169.0), MIT licence — see three-LICENSE.txt.

  three.module.js      the minified ESM build (three/build/three.module.min.js)
  OrbitControls.js     three/examples/jsm/controls/OrbitControls.js
  TrackballControls.js three/examples/jsm/controls/TrackballControls.js  (unused;
                       kept for anyone who prefers the free-rotation controls of
                       the original decompwlj.com pages)

Vendored rather than loaded from a CDN so the app works offline, behind a proxy
and from a plain file copy. The bare specifier "three" is resolved by the import
map at the top of index.html; nothing else needs configuring.

To upgrade:  npm pack three@<version>  and replace these three files.
