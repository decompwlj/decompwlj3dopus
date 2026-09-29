/* viz.js — the WebGL scene.
 *
 * One unit of world space is one decade on every axis, so the aspect is
 * equal by construction: the k = L plane really is at 45°, and how far a
 * point sits from it is exactly the level/weight margin log k − log L.
 *
 * Every point of the sequence is uploaded once.  Class visibility, the jump
 * filter and the sweep along n are uniforms applied in the vertex shader —
 * filtering never touches the buffers, so it costs nothing and stays smooth
 * whether the cloud holds 10³ or 10⁷ points.
 */

import * as THREE from 'three';
import { OrbitControls } from '../vendor/OrbitControls.js';

const SUP = ['⁰','¹','²','³','⁴','⁵','⁶','⁷','⁸','⁹'];
const sup = (n) => (n < 0 ? '⁻' : '') + String(Math.abs(n)).split('').map(d => SUP[+d]).join('');

const VERT = /* glsl */`
  attribute vec4 meta;          // x class(0 w,1 l,2 tie)  y log10 d  z index  w L==1
  uniform float uSize, uDpr, uAtten, uNMax, uHiOne, uRef, uOrthoAtt;
  uniform vec2  uD;
  uniform vec3  uShow;          // weight, level, tie
  uniform vec3  cWeight, cLevel, cTie, cOne;
  varying vec3 vColor;
  void main() {
    float cls = meta.x;
    float show = (cls > 0.5 && cls < 1.5) ? uShow.y : (cls > 1.5 ? uShow.z : uShow.x);
    if (meta.z > uNMax) show = 0.0;
    if (meta.y < uD.x || meta.y > uD.y) show = 0.0;

    vec4 mv = modelViewMatrix * vec4(position, 1.0);
    gl_Position = projectionMatrix * mv;

    float boost = 1.0;
    vColor = (cls > 1.5) ? cTie : ((cls > 0.5) ? cLevel : cWeight);
    if (uHiOne > 0.5 && meta.w > 0.5) { vColor = cOne; boost = 1.7; }

    /* uSize is the pixel size at the framing distance; nearer points grow. */
    float att = (uOrthoAtt > 0.0) ? uOrthoAtt : clamp(uRef / max(-mv.z, 0.001), 0.2, 8.0);
    gl_PointSize = show * uSize * uDpr * mix(1.0, att, uAtten) * boost;
  }`;

const FRAG = /* glsl */`
  uniform float uOpacity, uSoft;
  varying vec3 vColor;
  void main() {
    vec2 c = gl_PointCoord - 0.5;
    float r2 = dot(c, c);
    if (r2 > 0.25) discard;
    float a = uOpacity * mix(1.0, smoothstep(0.25, 0.02, r2), uSoft);
    gl_FragColor = vec4(vColor, a);
    #include <colorspace_fragment>
  }`;

export class Viz {
  constructor(canvas, overlay) {
    this.canvas = canvas;
    this.overlay = overlay;

    this.renderer = new THREE.WebGLRenderer({
      canvas, antialias: true, alpha: false, powerPreference: 'high-performance',
      preserveDrawingBuffer: true,
    });
    this.renderer.setPixelRatio(Math.min(devicePixelRatio || 1, 2));

    this.scene = new THREE.Scene();
    this.persp = new THREE.PerspectiveCamera(32, 1, 0.05, 4000);
    this.ortho = new THREE.OrthographicCamera(-1, 1, 1, -1, -5000, 5000);
    this.camera = this.persp;
    this.projection = 'persp';
    this.camera.position.set(14, 11, 18);
    this._orthoH = 8;
    this._autoRotate = false;

    this._makeControls();

    this.box = new THREE.Vector3(10, 10, 10);
    this.lo = [0, 0, 0];
    this.labels = [];
    this.dirty = true;
    this.tween = null;
    this.onPick = () => {};

    this.uniforms = {
      uSize:    { value: 2.4 },
      uDpr:     { value: this.renderer.getPixelRatio() },
      uAtten:   { value: 1 },
      uRef:     { value: 16 },
      uOrthoAtt:{ value: 0 },
      uNMax:    { value: 1e9 },
      uHiOne:   { value: 0 },
      uD:       { value: new THREE.Vector2(-9, 99) },
      uShow:    { value: new THREE.Vector3(1, 1, 1) },
      uOpacity: { value: 0.85 },
      uSoft:    { value: 1 },
      cWeight:  { value: new THREE.Color(0x4c8dff) },
      cLevel:   { value: new THREE.Color(0xff8c42) },
      cTie:     { value: new THREE.Color(0xffe45c) },
      cOne:     { value: new THREE.Color(0x55e0c0) },
    };

    this.material = new THREE.ShaderMaterial({
      uniforms: this.uniforms, vertexShader: VERT, fragmentShader: FRAG,
      transparent: true, depthWrite: true, depthTest: true,
    });

    this.geometry = new THREE.BufferGeometry();
    this.points = new THREE.Points(this.geometry, this.material);
    this.points.frustumCulled = false;
    this.scene.add(this.points);

    this.frame = new THREE.Group();
    this.scene.add(this.frame);

    new ResizeObserver(() => this.resize()).observe(canvas.parentElement);
    this.resize();

    canvas.addEventListener('pointermove', (e) => this._hover(e));
    canvas.addEventListener('pointerleave', () => this.onPick(null));
    this._raf = requestAnimationFrame(this._loop.bind(this));
  }

  /* OrbitControls caches a quaternion built from camera.up, so it is rebuilt
     whenever the camera object or its up vector changes. */
  _makeControls() {
    const target = this.controls ? this.controls.target.clone() : new THREE.Vector3();
    if (this.controls) this.controls.dispose();
    const c = new OrbitControls(this.camera, this.canvas);
    c.enableDamping = true;
    c.dampingFactor = 0.075;
    c.rotateSpeed = 0.85;
    c.zoomSpeed = 0.9;
    c.autoRotateSpeed = 0.55;
    c.autoRotate = this._autoRotate;
    c.target.copy(target);
    c.addEventListener('change', () => { this.dirty = true; });
    this.controls = c;
  }

  setProjection(kind) {
    if (kind === this.projection) return;
    const from = this.camera, to = kind === 'ortho' ? this.ortho : this.persp;
    to.position.copy(from.position);
    to.up.copy(from.up);
    to.quaternion.copy(from.quaternion);
    this.camera = to;
    this.projection = kind;
    this._makeControls();
    this.resize();
    this.dirty = true;
  }

  /* ── theme ─────────────────────────────────────────────── */

  applyTheme() {
    const cs = getComputedStyle(document.documentElement);
    const c = (n) => cs.getPropertyValue(n).trim() || '#888';
    this.scene.background = new THREE.Color().setStyle(c('--bg'));
    this.uniforms.cWeight.value.setStyle(c('--c-weight'));
    this.uniforms.cLevel.value.setStyle(c('--c-level'));
    this.uniforms.cTie.value.setStyle(c('--c-tie'));
    this.uniforms.cOne.value.setStyle(c('--c-one'));
    this.axCol = [c('--ax-x'), c('--ax-y'), c('--ax-z')];
    this.gridCol = c('--line-2');
    this.planeCol = c('--fg-3');
    if (this._frameArgs) this.buildFrame(...this._frameArgs);
    this.dirty = true;
  }

  resize() {
    const el = this.canvas.parentElement;
    const w = el.clientWidth || 1, h = el.clientHeight || 1;
    const aspect = w / h;
    this.persp.aspect = aspect;
    this.persp.updateProjectionMatrix();
    const H = this._orthoH;
    this.ortho.top = H; this.ortho.bottom = -H;
    this.ortho.left = -H * aspect; this.ortho.right = H * aspect;
    this.ortho.updateProjectionMatrix();
    this.renderer.setSize(w, h, false);
    this.uniforms.uDpr.value = this.renderer.getPixelRatio();
    const tight = w < 720;
    if (this._tight !== undefined && tight !== this._tight && this._frameArgs) this.buildFrame(...this._frameArgs);
    this._tight = tight;
    this.dirty = true;
  }

  /* ── geometry upload ───────────────────────────────────── */

  /** Bind the point cloud.  `fields` are three Float32Arrays already in
   *  log₁₀; `lo` their per-axis decade floor; `box` the decade span. */
  setPoints(dataset, fx, fy, fz, lo, box) {
    const n = dataset.count;
    const pos = new Float32Array(n * 3);
    const meta = new Float32Array(n * 4);
    const hx = box[0] / 2, hy = box[1] / 2, hz = box[2] / 2;
    const logd = dataset.field('d');
    const cls = dataset.cls, one = dataset.one;
    for (let i = 0; i < n; i++) {
      pos[i * 3]     = fx[i] - lo[0] - hx;
      pos[i * 3 + 1] = fy[i] - lo[1] - hy;
      pos[i * 3 + 2] = fz[i] - lo[2] - hz;
      meta[i * 4]     = cls[i];
      meta[i * 4 + 1] = logd[i];
      meta[i * 4 + 2] = i;
      meta[i * 4 + 3] = one[i];
    }
    this.geometry.dispose();
    this.geometry = new THREE.BufferGeometry();
    this.geometry.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    this.geometry.setAttribute('meta', new THREE.BufferAttribute(meta, 4));
    this.points.geometry = this.geometry;
    this.points.frustumCulled = false;
    this._pos = pos;
    this._meta = meta;
    this._n = n;
    this.dirty = true;
  }

  /* ── axes, ticks, decade box, k = L plane ──────────────── */

  buildFrame(names, lo, hi, box, opts) {
    this._frameArgs = [names, lo, hi, box, opts];
    for (const o of [...this.frame.children]) {
      this.frame.remove(o);
      o.geometry && o.geometry.dispose();
      o.material && o.material.dispose();
    }
    for (const l of this.labels) l.remove();
    this.labels = [];

    this.lo = lo; this.box = new THREE.Vector3(...box);
    const h = box.map(v => v / 2);
    const O = [-h[0], -h[1], -h[2]];

    const line = (pts, color, opacity = 1) => {
      const g = new THREE.BufferGeometry().setFromPoints(pts.map(p => new THREE.Vector3(...p)));
      const m = new THREE.LineBasicMaterial({ color: new THREE.Color().setStyle(color), transparent: opacity < 1, opacity });
      const o = new THREE.Line(g, m);
      this.frame.add(o);
      return o;
    };
    const seg = (pts, color, opacity = 1) => {
      const g = new THREE.BufferGeometry().setFromPoints(pts.map(p => new THREE.Vector3(...p)));
      const m = new THREE.LineBasicMaterial({ color: new THREE.Color().setStyle(color), transparent: opacity < 1, opacity });
      const o = new THREE.LineSegments(g, m);
      this.frame.add(o);
      return o;
    };

    /* the decade box: twelve faint edges */
    if (opts.box) {
      const [X, Y, Z] = box, e = [];
      const c = [[0,0,0],[X,0,0],[X,Y,0],[0,Y,0],[0,0,Z],[X,0,Z],[X,Y,Z],[0,Y,Z]]
                .map(p => [p[0] + O[0], p[1] + O[1], p[2] + O[2]]);
      for (const [i, j] of [[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]])
        e.push(c[i], c[j]);
      seg(e, this.gridCol, 0.34);
    }

    /* the three axes, from the box corner */
    const dirs = [[1,0,0],[0,1,0],[0,0,1]];
    for (let ax = 0; ax < 3; ax++) {
      const end = O.map((v, i) => v + dirs[ax][i] * box[ax]);
      line([O, end], this.axCol[ax]);
      this._label(names[ax].t, end.map((v, i) => v + dirs[ax][i] * 0.55), 'name ' + 'xyz'[ax]);

      if (!opts.box) continue;
      const span = hi[ax] - lo[ax];
      const tight = this.canvas.clientWidth < 720;
      const step = span > 26 ? 5 : span > 12 ? 2 : (tight ? 2 : 1);
      const ticks = [];
      const start = Math.ceil(lo[ax]);
      const off = (ax + 1) % 3;
      for (let v = start; v <= hi[ax] + 1e-9; v += step) {
        const t = v - lo[ax];
        const p = O.slice(); p[ax] += t;
        const q = p.slice(); q[off] -= box[off] * 0.02 + 0.1;
        ticks.push(p, q);
        /* the three axes share a corner; only X labels it, or the three
           minimum labels pile up on top of one another */
        if (v === start && ax > 0) continue;
        this._label('10' + sup(v), q.map((c2, i) => c2 + (i === off ? -0.3 : 0)), '');
      }
      seg(ticks, this.axCol[ax], 0.8);
    }

    /* the k = L plane, wherever k and L are actually bound */
    const ik = names.findIndex(s => s.f === 'k'), iL = names.findIndex(s => s.f === 'L');
    this.plane = null;
    if (opts.plane && ik >= 0 && iL >= 0) {
      const t0 = Math.max(lo[ik], lo[iL]), t1 = Math.min(hi[ik], hi[iL]);
      if (t1 > t0) {
        const other = 3 - ik - iL;
        const corner = (t, s) => {
          const p = [0, 0, 0];
          p[ik] = t - lo[ik] + O[ik];
          p[iL] = t - lo[iL] + O[iL];
          p[other] = O[other] + s * box[other];
          return p;
        };
        const A = corner(t0, 0), B = corner(t1, 0), C = corner(t1, 1), D = corner(t0, 1);
        const g = new THREE.BufferGeometry();
        g.setAttribute('position', new THREE.BufferAttribute(new Float32Array([
          ...A, ...B, ...C, ...A, ...C, ...D]), 3));
        const m = new THREE.MeshBasicMaterial({
          color: new THREE.Color().setStyle(this.planeCol),
          transparent: true, opacity: 0.075, side: THREE.DoubleSide, depthWrite: false,
        });
        const mesh = new THREE.Mesh(g, m);
        this.frame.add(mesh);
        this.plane = mesh;
        line([corner(t0, 0.5), corner(t1, 0.5)], this.planeCol, 0.5);
        this._label('k = L', corner(t1, 0.5), '');
      }
    }
    this.dirty = true;
  }

  _label(text, pos, cls) {
    const el = document.createElement('div');
    el.className = 'lbl ' + cls;
    el.textContent = typeof text === 'string' ? text : text.t;
    el._p = new THREE.Vector3(...pos);
    this.overlay.appendChild(el);
    this.labels.push(el);
  }

  _drawLabels() {
    const w = this.canvas.clientWidth, h = this.canvas.clientHeight;
    const v = new THREE.Vector3();
    for (const el of this.labels) {
      v.copy(el._p).project(this.camera);
      if (v.z > 1) { el.style.display = 'none'; continue; }
      el.style.display = '';
      /* keep labels clear of the frame edges and of the top and status bars */
      const name = el.classList.contains('name');
      const mx = name ? 78 : 16, my0 = name ? 62 : 14, my1 = name ? 46 : 26;
      const x = Math.min(w - mx, Math.max(mx, (v.x * 0.5 + 0.5) * w));
      const y = Math.min(h - my1, Math.max(my0, (-v.y * 0.5 + 0.5) * h));
      el.style.transform = `translate(-50%,-50%) translate(${x}px,${y}px)`;
    }
  }

  /* ── camera presets ────────────────────────────────────── */

  /** Camera basis for a view direction (camera sits at dir·D looking at the origin). */
  _basis(dir, up0) {
    const fwd = dir.clone().negate();
    const right = new THREE.Vector3().crossVectors(fwd, up0).normalize();
    const up = new THREE.Vector3().crossVectors(right, fwd).normalize();
    return { fwd, right, up };
  }

  /** Distance at which the decade box just fills the perspective frustum,
   *  and the half-height that just fits it under orthographic projection. */
  _fit(dir, up0) {
    const b = this.box, aspect = this.persp.aspect || 1;
    const tv = Math.tan(THREE.MathUtils.degToRad(this.persp.fov) / 2), th = tv * aspect;
    const { right, up } = this._basis(dir, up0);
    const c = new THREE.Vector3();
    let dist = 0, hh = 0;
    for (const sx of [-.5, .5]) for (const sy of [-.5, .5]) for (const sz of [-.5, .5]) {
      c.set(sx * b.x, sy * b.y, sz * b.z);
      const cz = c.dot(dir), cr = Math.abs(c.dot(right)), cu = Math.abs(c.dot(up));
      dist = Math.max(dist, cz + cr / th, cz + cu / tv);
      hh = Math.max(hh, cu, cr / aspect);
    }
    return { dist: dist * 1.1, half: hh * 1.14 };
  }

  view(kind) {
    const dirs = {
      iso:  new THREE.Vector3(1, 0.58, 1.18),
      xy:   new THREE.Vector3(0, 0, 1),
      xz:   new THREE.Vector3(0, 1, 0.0001),
      yz:   new THREE.Vector3(1, 0, 0.0001),
      edge: new THREE.Vector3(1, 1, 0.0001),
    };
    const dir = (dirs[kind] || dirs.iso).normalize();
    const up = kind === 'xz' ? new THREE.Vector3(0, 0, -1) : new THREE.Vector3(0, 1, 0);
    const fit = this._fit(dir, up);
    this._orthoH = fit.half;
    this.uniforms.uRef.value = fit.dist;
    this.resize();
    const to = dir.clone().multiplyScalar(this.camera.isOrthographicCamera
      ? Math.max(fit.dist, 40) : fit.dist);
    this.tween = {
      t: 0, t0: performance.now(), from: this.camera.position.clone(), to,
      fromUp: this.camera.up.clone(), toUp: up,
    };
    this.controls.target.set(0, 0, 0);
    this.dirty = true;
  }

  frameAll() { this.view('iso'); }

  /* ── picking ───────────────────────────────────────────── */

  _hover(e) {
    const r = this.canvas.getBoundingClientRect();
    this._ptr = { x: e.clientX - r.left, y: e.clientY - r.top, cx: e.clientX, cy: e.clientY };
    this._pickDirty = true;
  }

  _doPick() {
    this._pickDirty = false;
    const p = this._ptr;
    if (!p || !this._pos || !this._n) return this.onPick(null);
    const w = this.canvas.clientWidth, h = this.canvas.clientHeight;
    const nx = (p.x / w) * 2 - 1, ny = -(p.y / h) * 2 + 1;
    const ortho = this.camera.isOrthographicCamera;
    let ox, oy, oz, ux, uy, uz, slack;
    if (ortho) {
      const o = new THREE.Vector3(nx, ny, -1).unproject(this.camera);
      const d = new THREE.Vector3(0, 0, -1).applyQuaternion(this.camera.quaternion);
      ox = o.x; oy = o.y; oz = o.z; ux = d.x; uy = d.y; uz = d.z;
      slack = (2 * this._orthoH / h) * 9;                 /* 9 px, in world units */
    } else {
      const q = new THREE.Vector3(nx, ny, 0.5).unproject(this.camera);
      const o = this.camera.position;
      const dx = q.x - o.x, dy = q.y - o.y, dz = q.z - o.z;
      const inv = 1 / Math.hypot(dx, dy, dz);
      ox = o.x; oy = o.y; oz = o.z; ux = dx * inv; uy = dy * inv; uz = dz * inv;
      slack = 2 * Math.tan(THREE.MathUtils.degToRad(this.camera.fov) / 2) / h * 9;
    }

    const pos = this._pos, meta = this._meta, u = this.uniforms;
    const nMax = u.uNMax.value, dLo = u.uD.value.x, dHi = u.uD.value.y;
    const show = u.uShow.value;
    let best = -1, bestScore = Infinity;
    for (let i = 0; i < this._n; i++) {
      const c = meta[i * 4];
      const vis = (c > 0.5 && c < 1.5) ? show.y : (c > 1.5 ? show.z : show.x);
      if (vis < 0.5) continue;
      if (meta[i * 4 + 2] > nMax) break;
      const md = meta[i * 4 + 1];
      if (md < dLo || md > dHi) continue;
      const ax = pos[i * 3] - ox, ay = pos[i * 3 + 1] - oy, az = pos[i * 3 + 2] - oz;
      const t = ax * ux + ay * uy + az * uz;
      if (t <= 0.01) continue;
      const px = ax - ux * t, py = ay - uy * t, pz = az - uz * t;
      const perp2 = px * px + py * py + pz * pz;
      const lim = ortho ? slack : slack * t;
      if (perp2 > lim * lim) continue;
      const score = Math.sqrt(perp2) / (ortho ? 1 : t) + t * 1e-4;
      if (score < bestScore) { bestScore = score; best = i; }
    }
    this.onPick(best < 0 ? null : { i: best, cx: this._ptr.cx, cy: this._ptr.cy });
  }

  /* ── loop ──────────────────────────────────────────────── */

  set autoRotate(v) { this._autoRotate = v; this.controls.autoRotate = v; this.dirty = this.dirty || v; }

  _loop() {
    this._raf = requestAnimationFrame(this._loop.bind(this));
    const now = performance.now();
    if (this.tween) {
      const tw = this.tween;
      tw.t = Math.min(1, (now - tw.t0) / 520);
      const e = tw.t < 0.5 ? 4 * tw.t ** 3 : 1 - (-2 * tw.t + 2) ** 3 / 2;
      this.camera.position.lerpVectors(tw.from, tw.to, e);
      this.camera.up.lerpVectors(tw.fromUp, tw.toUp, e).normalize();
      if (tw.t >= 1) { this.tween = null; this._makeControls(); }
      this.dirty = true;
    }
    if (this.controls.autoRotate) this.dirty = true;
    if (this.controls.update()) this.dirty = true;
    this.uniforms.uOrthoAtt.value =
      this.camera.isOrthographicCamera ? THREE.MathUtils.clamp(this.camera.zoom, 0.15, 8) : 0;
    if (this._pickDirty) this._doPick();
    if (!this.dirty) return;
    this.dirty = false;
    this.renderer.render(this.scene, this.camera);
    this._drawLabels();
  }

  stop() { cancelAnimationFrame(this._raf); this._raf = 0; }
  start() { if (!this._raf) { this.dirty = true; this._raf = requestAnimationFrame(this._loop.bind(this)); } }

  snapshot() {
    this.renderer.render(this.scene, this.camera);
    return this.canvas.toDataURL('image/png');
  }
}
