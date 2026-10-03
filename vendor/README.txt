three.js r169 (0.169.0), MIT licence — see LICENSE-three.js.txt.

  three.module.min.js                 the minified ESM build (three/build/three.module.min.js)
  addons/controls/OrbitControls.js    three/examples/jsm/controls/OrbitControls.js

Vendored rather than loaded from a CDN, so the app works offline, behind a proxy and from a
plain copy of the files. The import map at the top of index.html resolves "three" and
"three/addons/" to these files; nothing else needs configuring.

To upgrade: npm pack three@<version> and replace these two files.
