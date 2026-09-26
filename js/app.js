/* app.js — wiring: catalogue → dataset → scene → controls. */

import { NormalBlending, AdditiveBlending } from 'three';
import { loadCatalog, filterCatalog, VirtualList } from './catalog.js';
import { Dataset, FIELDS } from './data.js';
import { Viz } from './viz.js';

const DATA = 'data';
const $ = (s) => document.querySelector(s);
const fmt = new Intl.NumberFormat('en-US');
const pct = (x) => (100 * x).toFixed(2) + '%';

const state = {
  id: null, x: 'k', y: 'L', z: 'd', view: 'iso', mode: 'solid',
  show: { weight: true, level: true, tie: true }, hiOne: false,
  reveal: 1, dLo: 0, dHi: 1,
};

let all = [], view = [], rec = null, ds = null, viz = null, list = null, abort = null;

/* console hook — `decompwlj.ds`, `.rec`, `.viz`, `.state` for poking around */
globalThis.decompwlj = {
  get all() { return all; }, get rec() { return rec; }, get ds() { return ds; },
  get viz() { return viz; }, get state() { return state; }, select: (id) => select(id),
};

/* ── URL state ─────────────────────────────────────────────── */

function readHash() {
  const p = new URLSearchParams(location.hash.slice(1));
  for (const [k, v] of [['s', 'id'], ['x', 'x'], ['y', 'y'], ['z', 'z'], ['v', 'view'], ['m', 'mode']])
    if (p.get(k)) state[v] = p.get(k);
  if (p.get('h')) state.hiOne = p.get('h') === '1';
}
function writeHash() {
  const p = new URLSearchParams();
  p.set('s', state.id); p.set('x', state.x); p.set('y', state.y); p.set('z', state.z);
  p.set('v', state.view); p.set('m', state.mode);
  if (state.hiOne) p.set('h', '1');
  history.replaceState(null, '', '#' + p.toString());
}

/* ── boot ──────────────────────────────────────────────────── */

async function boot() {
  readHash();
  viz = new Viz($('#gl'), $('#overlay'));
  viz.applyTheme();
  viz.onPick = showTip;

  list = new VirtualList($('#catList'), $('#catSpacer'), $('#catRows'), (id) => select(id));

  for (const sel of ['#axX', '#axY', '#axZ']) {
    const el = $(sel);
    for (const [key, f] of Object.entries(FIELDS)) {
      const o = document.createElement('option');
      o.value = key; o.textContent = f.label;
      el.appendChild(o);
    }
  }

  try {
    all = await loadCatalog(`${DATA}/catalog.csv`, (_, sofar) => {
      view = filterCatalog(sofar, $('#search').value);
      list.setItems(view);
      $('#searchCount').textContent = fmt.format(view.length);
      $('#search').placeholder = `Search ${fmt.format(sofar.length)} sequences…`;
    });
  } catch (e) {
    return fail(`Could not read ${DATA}/catalog.csv — ${e.message}.\nServe this folder over HTTP (ES modules and fetch do not work from file://).`);
  }
  if (!all.length) return fail('The catalogue is empty.');

  view = all;
  list.setItems(view);
  $('#searchCount').textContent = fmt.format(all.length);
  $('#search').placeholder = `Search ${fmt.format(all.length)} sequences…`;

  bindUI();
  const start = all.find(r => r.id === state.id) || all[0];
  await select(start.id);
  const boot = $('#boot');
  boot.classList.add('gone');
  setTimeout(() => boot.remove(), 500);
}

function fail(msg) {
  const b = $('#boot');
  b.classList.add('err');
  $('#bootMsg').textContent = msg;
  $('#bootMsg').style.whiteSpace = 'pre-line';
  $('#bootMsg').style.maxWidth = '46ch';
}

/* ── selecting and loading a sequence ──────────────────────── */

async function select(id) {
  const r = all.find(s => s.id === id);
  if (!r || (rec && rec.id === id)) return;
  if (abort) abort.abort();
  abort = new AbortController();
  const signal = abort.signal;

  rec = r; state.id = id;
  list.select(id); list.scrollTo(id);
  writeHash();

  $('#seqName').textContent = r.name;
  $('#seqSub').textContent = `${r.anumber} · ${r.family} · n = ${fmt.format(r.n0)} … ${fmt.format(r.n0 + r.terms - 1)}`;
  $('#note').textContent = r.note || '';
  $('#oeisLink').href = `https://oeis.org/${r.anumber}`;
  $('#oeisLink').textContent = `OEIS ${r.anumber}`;
  $('#loadState').textContent = 'loading…';
  $('#progressBar').style.opacity = '1';
  $('#progressBar').style.width = '2%';

  ds = new Dataset(r, DATA);
  rebuildFrame();
  viz.view(state.view);

  try {
    await ds.load((d, frac) => {
      if (signal.aborted) return;
      $('#progressBar').style.width = (frac * 100).toFixed(1) + '%';
      $('#loadState').textContent = frac < 1 ? `${fmt.format(d.count)} pts…` : '';
      upload();
      refreshStats();
    }, signal);
  } catch (e) {
    if (e.name === 'AbortError') return;
    $('#loadState').textContent = 'load failed';
    return fail(`Could not read the data for ${r.name}: ${e.message}`);
  }
  if (signal.aborted) return;

  $('#progressBar').style.opacity = '0';
  $('#loadState').textContent = '';
  applyFilters();
}

/* ── frame + upload ────────────────────────────────────────── */

/** Axis bounds come from the catalogue, not the loaded rows, so the box
 *  is fixed before the first chunk lands and nothing shifts while loading. */
function bounds(f) {
  const L10 = (v) => Math.log10(Math.max(1, v));
  switch (f) {
    case 'k': return [0, L10(rec.kmax)];
    case 'L': return [0, L10(rec.Lmax)];
    case 'd': return [L10(rec.dmin), L10(rec.dmax)];
    case 'a': return [L10(rec.amin), L10(rec.amax)];
    case 'l': return [0, L10(rec.amax)];
    case 'n': return [L10(rec.n0), L10(rec.n0 + rec.terms - 1)];
  }
  return [0, 1];
}

let frame = null;

function rebuildFrame() {
  const fs = [state.x, state.y, state.z];
  const lo = [], hi = [];
  for (const f of fs) {
    const [a, b] = bounds(f);
    lo.push(Math.floor(a));
    hi.push(Math.max(Math.floor(a) + 1, Math.ceil(b)));
  }
  const box = [0, 1, 2].map(i => hi[i] - lo[i]);
  frame = { fs, lo, hi, box };
  viz.buildFrame(
    fs.map(f => ({ f, t: FIELDS[f].label })),
    lo, hi, box,
    { box: $('#showBox').checked, plane: $('#showPlane').checked }
  );
  const [dl, dh] = bounds('d');
  $('#fDLo').value = 0; $('#fDHi').value = 1;
  state.dLo = 0; state.dHi = 1;
  viz.uniforms.uD.value.set(dl - 0.001, dh + 0.001);
  $('#fDOut').textContent = 'all';
  $('#stRange').textContent =
    `a ≤ ${fmt.format(rec.amax)} · k ≤ ${fmt.format(rec.kmax)} · L ≤ ${fmt.format(rec.Lmax)} · d ≤ ${fmt.format(rec.dmax)}`;
}

function upload() {
  if (!ds || !frame) return;
  viz.setPoints(ds, ds.field(frame.fs[0]), ds.field(frame.fs[1]), ds.field(frame.fs[2]), frame.lo, frame.box);
}

/* ── filters and read-outs ─────────────────────────────────── */

function applyFilters() {
  if (!ds) return;
  const u = viz.uniforms;
  u.uShow.value.set(state.show.weight ? 1 : 0, state.show.level ? 1 : 0, state.show.tie ? 1 : 0);
  u.uNMax.value = Math.round(state.reveal * ds.count);
  u.uHiOne.value = state.hiOne ? 1 : 0;
  const [dl, dh] = bounds('d');
  const span = dh - dl || 1;
  u.uD.value.set(dl + state.dLo * span - 0.001, dl + state.dHi * span + 0.001);
  viz.dirty = true;
  refreshStats();
}

let statTimer = 0;
function refreshStats() {
  clearTimeout(statTimer);
  statTimer = setTimeout(() => {
    if (!ds) return;
    const u = viz.uniforms;
    const s = ds.stats(u.uD.value.x, u.uD.value.y, u.uNMax.value);
    $('#lgWeight').textContent = fmt.format(s.weight);
    $('#lgLevel').textContent = fmt.format(s.level);
    $('#lgTie').textContent = fmt.format(s.tie);
    $('#lgOne').textContent = fmt.format(s.one);
    $('#stDots').textContent = `${fmt.format(s.shown)} dots` +
      (ds.complete ? '' : ` of ${fmt.format(rec.terms)} terms loading`);

    const dl = document.createDocumentFragment();
    const put = (k, v, cls) => {
      const dt = document.createElement('dt'); dt.textContent = k;
      const dd = document.createElement('dd'); dd.textContent = v; if (cls) dd.className = cls;
      dl.append(dt, dd);
    };
    put('terms', fmt.format(rec.terms));
    put('decomposable', fmt.format(rec.decomposable));
    put('non-decomposable', fmt.format(rec.terms - rec.decomposable));
    put('level  k > L', `${fmt.format(rec.level)}  ${pct(rec.level / rec.decomposable)}`, 'c-level');
    put('weight  k ≤ L', `${fmt.format(rec.weight)}  ${pct(rec.weight / rec.decomposable)}`, 'c-weight');
    put('ties  k = L', fmt.format(rec.ties));
    put('level line L = 1', fmt.format(rec.level_one));
    $('#stats').replaceChildren(dl);
  }, 40);
}

/* ── point read-out ────────────────────────────────────────── */

function showTip(hit) {
  const tip = $('#tip');
  if (!hit || !ds) { tip.hidden = true; return; }
  const p = ds.point(hit.i);
  const cname = p.cls === 1 ? 'level' : p.cls === 2 ? 'tie' : 'weight';
  const cclass = p.cls === 1 ? 'cl' : p.cls === 2 ? 'ct' : 'cw';
  tip.innerHTML =
    `n = <b>${fmt.format(p.n)}</b>   <span class="${cclass}">${cname}${p.one ? ' · L = 1' : ''}</span>\n` +
    `a = ${fmt.format(p.a)}\n` +
    `k = ${fmt.format(p.k)}\n` +
    `L = ${fmt.format(p.L)}\n` +
    `d = ${fmt.format(p.d)}\n` +
    `ℓ = ${fmt.format(p.l)} = k·L\n` +
    `${fmt.format(p.k)}·${fmt.format(p.L)} + ${fmt.format(p.d)} = ${fmt.format(p.k * p.L + p.d)}`;
  tip.hidden = false;
  const r = $('#stage').getBoundingClientRect();
  const w = tip.offsetWidth, h = tip.offsetHeight;
  let x = hit.cx - r.left + 16, y = hit.cy - r.top + 16;
  if (x + w > r.width - 8) x = hit.cx - r.left - w - 16;
  if (y + h > r.height - 8) y = hit.cy - r.top - h - 16;
  tip.style.left = x + 'px';
  tip.style.top = y + 'px';
}

/* ── UI ────────────────────────────────────────────────────── */

function bindUI() {
  $('#search').addEventListener('input', (e) => {
    view = filterCatalog(all, e.target.value);
    list.setItems(view);
    $('#searchCount').textContent = fmt.format(view.length);
  });

  for (const [sel, key] of [['#axX', 'x'], ['#axY', 'y'], ['#axZ', 'z']]) {
    const el = $(sel);
    el.value = state[key];
    el.addEventListener('change', () => {
      state[key] = el.value;
      rebuildFrame(); upload(); applyFilters(); writeHash();
    });
  }

  document.querySelectorAll('[data-view]').forEach(b => b.addEventListener('click', () => {
    document.querySelectorAll('[data-view]').forEach(o => o.classList.toggle('on', o === b));
    state.view = b.dataset.view;
    /* a flat view is only a true plate without foreshortening */
    setProjection(state.view === 'iso' ? 'persp' : 'ortho');
    viz.view(state.view);
    writeHash();
  }));

  document.querySelectorAll('[data-proj]').forEach(b => b.addEventListener('click', () => {
    setProjection(b.dataset.proj);
    viz.view(state.view);
  }));

  document.querySelectorAll('[data-mode]').forEach(b => b.addEventListener('click', () => {
    document.querySelectorAll('[data-mode]').forEach(o => o.classList.toggle('on', o === b));
    setMode(b.dataset.mode);
  }));

  document.querySelectorAll('.lg').forEach(b => b.addEventListener('click', () => {
    const c = b.dataset.class;
    if (c === 'one') { state.hiOne = !state.hiOne; $('#hiOne').checked = state.hiOne; }
    else state.show[c] = !state.show[c];
    b.classList.toggle('off', c === 'one' ? !state.hiOne : !state.show[c]);
    applyFilters(); writeHash();
  }));

  const slider = (sel, out, fn, fmtOut) => {
    const el = $(sel);
    el.addEventListener('input', () => {
      fn(parseFloat(el.value));
      if (out) $(out).textContent = fmtOut ? fmtOut(parseFloat(el.value)) : el.value;
    });
  };
  slider('#ptSize', '#ptSizeOut', v => { viz.uniforms.uSize.value = v; viz.dirty = true; });
  slider('#ptOpacity', '#ptOpacityOut', v => { viz.uniforms.uOpacity.value = v; viz.dirty = true; });
  slider('#fReveal', '#fRevealOut', v => { state.reveal = v; applyFilters(); },
         v => v >= 1 ? 'all' : fmt.format(Math.round(v * (ds ? ds.count : 0))));

  const dOut = () => {
    const [dl, dh] = bounds('d'); const span = dh - dl || 1;
    const a = Math.round(10 ** (dl + state.dLo * span)), b = Math.round(10 ** (dl + state.dHi * span));
    $('#fDOut').textContent = (state.dLo <= 0 && state.dHi >= 1) ? 'all' : `${fmt.format(a)} … ${fmt.format(b)}`;
  };
  $('#fDLo').addEventListener('input', () => {
    state.dLo = Math.min(parseFloat($('#fDLo').value), state.dHi);
    $('#fDLo').value = state.dLo; applyFilters(); dOut();
  });
  $('#fDHi').addEventListener('input', () => {
    state.dHi = Math.max(parseFloat($('#fDHi').value), state.dLo);
    $('#fDHi').value = state.dHi; applyFilters(); dOut();
  });
  $('#resetFilters').addEventListener('click', () => {
    state.reveal = 1; state.dLo = 0; state.dHi = 1;
    state.show = { weight: true, level: true, tie: true };
    $('#fReveal').value = 1; $('#fRevealOut').textContent = 'all';
    $('#fDLo').value = 0; $('#fDHi').value = 1; $('#fDOut').textContent = 'all';
    document.querySelectorAll('.lg').forEach(b => b.classList.toggle('off', b.dataset.class === 'one' ? !state.hiOne : false));
    applyFilters();
  });

  $('#showPlane').addEventListener('change', rebuildAndUpload);
  $('#showBox').addEventListener('change', rebuildAndUpload);
  $('#attenuate').addEventListener('change', (e) => { viz.uniforms.uAtten.value = e.target.checked ? 1 : 0; viz.dirty = true; });
  $('#hiOne').addEventListener('change', (e) => {
    state.hiOne = e.target.checked;
    document.querySelector('.lg[data-class="one"]').classList.toggle('off', !state.hiOne);
    applyFilters(); writeHash();
  });
  $('#spin').addEventListener('change', (e) => { viz.autoRotate = e.target.checked; });

  let playing = 0;
  $('#play').addEventListener('click', () => {
    if (playing) { cancelAnimationFrame(playing); playing = 0; $('#play').textContent = '▶ sweep along the sequence'; return; }
    $('#play').textContent = '❚❚ pause';
    let t = state.reveal >= 1 ? 0 : state.reveal, last = performance.now();
    const step = (now) => {
      t += (now - last) / 9000; last = now;
      if (t >= 1) { t = 1; }
      state.reveal = t;
      $('#fReveal').value = t;
      $('#fRevealOut').textContent = t >= 1 ? 'all' : fmt.format(Math.round(t * ds.count));
      applyFilters();
      if (t >= 1) { playing = 0; $('#play').textContent = '▶ sweep along the sequence'; return; }
      playing = requestAnimationFrame(step);
    };
    playing = requestAnimationFrame(step);
  });

  $('#snap').addEventListener('click', () => {
    const a = document.createElement('a');
    a.href = viz.snapshot();
    a.download = `decompwlj_${rec.id}_${state.x}${state.y}${state.z}.png`;
    a.click();
  });
  $('#share').addEventListener('click', async () => {
    try { await navigator.clipboard.writeText(location.href); $('#share').textContent = 'copied'; }
    catch { $('#share').textContent = location.hash; }
    setTimeout(() => { $('#share').textContent = 'copy link'; }, 1400);
  });

  $('#themeBtn').addEventListener('click', toggleTheme);
  $('#catToggle').addEventListener('click', () => toggle('cat'));
  $('#catShow').addEventListener('click', () => toggle('cat'));
  $('#panelToggle').addEventListener('click', () => toggle('panel'));

  document.addEventListener('visibilitychange', () => {
    if (document.hidden) viz.stop(); else viz.start();
  });

  addEventListener('keydown', (e) => {
    if (e.target.matches('input, select, textarea')) return;
    const k = e.key.toLowerCase();
    if (k === 'c') toggle('cat');
    else if (k === 'p') toggle('panel');
    else if (k === 't') toggleTheme();
    else if (k === ' ') { e.preventDefault(); $('#play').click(); }
    else if ('12345'.includes(k)) document.querySelectorAll('[data-view]')[+k - 1].click();
    else if (k === 'arrowdown' || k === 'arrowup') {
      e.preventDefault();
      const i = view.findIndex(r => r.id === state.id);
      const j = Math.max(0, Math.min(view.length - 1, i + (k === 'arrowdown' ? 1 : -1)));
      if (view[j]) select(view[j].id);
    }
  });

  setMode(state.mode);
  document.querySelectorAll('[data-view]').forEach(o => o.classList.toggle('on', o.dataset.view === state.view));
  document.querySelectorAll('[data-mode]').forEach(o => o.classList.toggle('on', o.dataset.mode === state.mode));
  $('#hiOne').checked = state.hiOne;
  document.querySelector('.lg[data-class="one"]').classList.toggle('off', !state.hiOne);
}

function rebuildAndUpload() { rebuildFrame(); upload(); applyFilters(); }

function setProjection(kind) {
  state.proj = kind;
  document.querySelectorAll('[data-proj]').forEach(o => o.classList.toggle('on', o.dataset.proj === kind));
  viz.setProjection(kind);
}

function setMode(mode) {
  state.mode = mode;
  const m = viz.material;
  if (mode === 'density') {
    m.blending = AdditiveBlending; m.depthWrite = false;
    viz.uniforms.uOpacity.value = Math.min(0.3, parseFloat($('#ptOpacity').value));
  } else {
    m.blending = NormalBlending; m.depthWrite = true;
    viz.uniforms.uOpacity.value = parseFloat($('#ptOpacity').value);
  }
  $('#ptOpacityOut').textContent = viz.uniforms.uOpacity.value.toFixed(2);
  m.needsUpdate = true;
  viz.dirty = true;
  writeHash();
}

function toggle(what) {
  const app = $('#app');
  /* matchMedia, not innerWidth — they disagree under mobile emulation and zoom */
  if (what === 'cat' && matchMedia('(max-width: 720px)').matches) { app.classList.toggle('show-cat'); return; }
  if (what === 'panel' && matchMedia('(max-width: 1080px)').matches) { app.classList.toggle('show-panel'); return; }
  app.classList.toggle(what === 'cat' ? 'no-cat' : 'no-panel');
  requestAnimationFrame(() => viz.resize());
  setTimeout(() => viz.resize(), 320);
}

function toggleTheme() {
  const root = document.documentElement;
  const light = root.dataset.theme === 'light';
  root.dataset.theme = light ? 'dark' : 'light';
  $('#themeBtn').textContent = light ? 'light theme' : 'dark theme';
  try { localStorage.setItem('decompwlj.theme', root.dataset.theme); } catch {}
  viz.applyTheme();
}

try {
  const t = localStorage.getItem('decompwlj.theme');
  if (t) document.documentElement.dataset.theme = t;
  if (t === 'light') addEventListener('DOMContentLoaded', () => { $('#themeBtn').textContent = 'dark theme'; });
} catch {}

boot().catch(e => fail(e.message || String(e)));
