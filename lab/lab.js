/* lab.js — decompwlj Lab: each sequence against its residue controls.
 *
 * The page loads one sequence from the published chunks of decompwlj 3D (the same files, the
 * same checks: exact divisions, continuity, the catalogue's counts, the SHA-256 fingerprint),
 * draws a random subsample of its decomposable records, and asks worker threads for the
 * record-level control at M = lcm(1..n) = 1, 2, 6, 12, …, 360360, or a modulus of one's own (see worker.js).  Every chart reads
 * the same three filters: a range of n, one jump d, and the control's modulus M.
 */
(() => {
'use strict';

const $ = (s) => document.querySelector(s);
const DATA = '../data';
/* lcm(1, …, n): each step adds exactly one prime power, named in ADDS, so a change in Δ at a step is
 * attributable to that prime power alone (2, 3, 4 = 2², 5, 7, 8 = 2³, 9 = 3², 11, 13) */
const LADDER = [1, 2, 6, 12, 60, 420, 840, 2520, 27720, 360360];
const ADDS = { 2: '2', 6: '3', 12: '4 = 2²', 60: '5', 420: '7', 840: '8 = 2³', 2520: '9 = 3²', 27720: '11', 360360: '13' };
const OWN_MAX = 10000000;        /* a modulus of your own: 2 … 10^7 */
const FROZEN_MAX = 0.5;          /* a step where more than half the records are frozen says little */
const REPS = 2;
const BINS = 30;
const LOCK_R = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113];
const MAXC = 60;                 /* columns k = 2..MAXC+1, lines L = 1..MAXC */
const FORCED_D = 94906265;       /* above this d² passes 2^53, beyond every ℓ */

const st = {
  cat: [], byId: new Map(), byA: new Map(),
  rec: null,                     /* catalogue row, or the pasted sequence's own */
  N: 0, n0: 1, a: null, d: null, k: null, L: null,
  dec: null,                     /* rows that decompose */
  sub: null,                     /* the control's records (rows), ascending */
  ctrl: new Map(),               /* `${M}|${rep}` → { lp, k, frozen } */
  M: 1, rLo: 0, rHi: -1, g: null, colMode: 'k',
  gen: 0, cgen: 0, workers: [], pending: 0, done: 0, total: 0,
  sample: 25000,
};

/* ═══════════════ small things ═══════════════ */

const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const fmtInt = (x) => x.toLocaleString('en-US');
const fmtNum = (x) => x < 1e6 ? fmtInt(x) : x.toExponential(2).replace('e+', 'e');
const pct = (x, dp = 2) => (100 * x).toFixed(dp) + ' %';
const signed = (v, dp) => { const t = Math.abs(v).toFixed(dp); return +t === 0 ? '0' : (v < 0 ? '−' : '+') + t; };
const fz = (z) => (Math.abs(z) < 0.05 ? '0.0' : (z < 0 ? '−' : '') + Math.abs(z).toFixed(1));
const pts = (x, dp = 2) => { const t = Math.abs(100 * x).toFixed(dp); return +t === 0 ? t : (x > 0 ? '+' : '−') + t; };
const setStatus = (t, cls = '') => { const s = $('#status'); s.textContent = t; s.className = 'status ' + cls; };
const setProgress = (f) => { $('#progress i').style.width = `${Math.round(100 * Math.min(1, f))}%`; };

function parseCSV(text) {
  const rows = []; let row = [], f = '', q = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (q) {
      if (c === '"') { if (text[i + 1] === '"') { f += '"'; i++; } else q = false; }
      else f += c;
    } else if (c === '"') q = true;
    else if (c === ',') { row.push(f); f = ''; }
    else if (c === '\n') { row.push(f); rows.push(row); row = []; f = ''; }
    else if (c !== '\r') f += c;
  }
  if (f || row.length) { row.push(f); rows.push(row); }
  const head = rows.shift();
  const NUM = new Set(['n0', 'terms', 'chunks', 'chunk_rows', 'decomposable', 'level', 'weight', 'ties', 'level_one', 'forced',
                       'amin', 'amax', 'kmax', 'Lmax', 'dmin', 'dmax']);
  return rows.filter(r => r.length === head.length).map(r => {
    const o = {};
    head.forEach((h, j) => { o[h] = NUM.has(h) ? Number(r[j]) : r[j]; });
    return o;
  });
}

async function fetchBytes(url) {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`${res.status}: ${url}`);
  let b = new Uint8Array(await res.arrayBuffer());
  if (b.length > 2 && b[0] === 0x1f && b[1] === 0x8b) {
    if (typeof DecompressionStream === 'undefined') throw new Error('this browser cannot inflate gzip (no DecompressionStream)');
    b = new Uint8Array(await new Response(new Blob([b]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer());
  }
  return b;
}

/** A 'dwj1' / 'dwj2' chunk (tools/compact.py) into rows [at, at + n) of the columns. Every
 *  division is checked exact, and k > d, and the stored factor is the smaller one. */
function readChunk(bytes, at, url) {
  if (!(bytes.length > 4 && bytes[0] === 0x64 && bytes[1] === 0x77 && bytes[2] === 0x6a && (bytes[3] === 0x31 || bytes[3] === 0x32)))
    throw new Error(`${url}: not a dwj chunk`);
  const two = bytes[3] === 0x32, m = two ? bytes[4] : 0;
  let i = two ? 5 : 4;
  if (m > 2) throw new Error(`${url}: unknown jump coding ${m}`);
  const rd = () => {
    let v = 0, s = 1, c;
    do { if (i >= bytes.length) throw new Error(`${url}: truncated`); c = bytes[i++]; v += (c & 127) * s; s *= 128; } while (c & 128);
    return v;
  };
  const unzz = (z) => (z % 2 === 0 ? z / 2 : -(z + 1) / 2);
  const n = rd(); let a = rd();
  if (at + n > st.N) throw new Error(`${url}: ${n} rows do not fit`);
  let p = 0, pp = 0;
  for (let r = at; r < at + n; r++) {
    const v = rd(), d = m === 0 ? v : m === 1 ? p + unzz(v) : 2 * p - pp + unzz(v);
    if (!(d >= 1) || d >= 2 ** 53) throw new Error(`${url}: jump out of range at row ${r}`);
    st.d[r] = d; pp = p; p = d;
  }
  for (let r = at; r < at + n; r++) {
    const d = st.d[r], s = rd();
    st.a[r] = a;
    if (s === 0) {
      if (a > 2 * d) throw new Error(`${url}: row ${r} decomposes but is stored as not`);
      st.k[r] = 0; st.L[r] = 0;
    } else {
      const h = Math.floor(s / 2), q = (a - d) / h;
      const k = s % 2 ? q : h, L = s % 2 ? h : q;
      if (!Number.isInteger(q) || q * h !== a - d || k <= d || (s % 2 ? L >= k : k > L))
        throw new Error(`${url}: row ${r} does not fit a = ${a}, d = ${d}`);
      st.k[r] = k; st.L[r] = L;
    }
    a += d;
  }
  if (i !== bytes.length) throw new Error(`${url}: trailing bytes`);
  return n;
}

/* ═══════════════ catalogue and picker ═══════════════ */

async function loadCatalog() {
  const res = await fetch(`${DATA}/catalog.csv`, { cache: 'no-cache' });
  if (!res.ok) throw new Error(`catalog.csv: ${res.status}`);
  st.cat = parseCSV(await res.text());
  for (const r of st.cat) { st.byId.set(r.id, r); st.byA.set(r.anumber, r); r._s = `${r.anumber} ${r.id} ${r.name} ${r.alias} ${r.family}`.toLowerCase(); }
  $('#q').placeholder = `Search ${fmtInt(st.cat.length)} sequences: name, A-number or family`;
}

function initPicker() {
  const q = $('#q'), ul = $('#hits');
  let hits = [], sel = -1;
  const close = () => { ul.hidden = true; q.setAttribute('aria-expanded', 'false'); };
  const show = () => {
    const t = q.value.trim().toLowerCase();
    if (!t) { close(); return; }
    const words = t.split(/\s+/);
    hits = st.cat.filter(r => words.every(w => r._s.includes(w))).slice(0, 40);
    ul.replaceChildren(...hits.map((r, j) => {
      const li = document.createElement('li');
      li.setAttribute('role', 'option'); li.id = `hit${j}`;
      const an = document.createElement('span'); an.className = 'an'; an.textContent = r.anumber;
      const nm = document.createElement('span'); nm.className = 'nm'; nm.textContent = r.alias || r.name;
      li.append(an, nm);
      li.addEventListener('mousedown', (e) => { e.preventDefault(); pick(r); });
      return li;
    }));
    if (!hits.length) { const li = document.createElement('li'); li.textContent = 'no sequence matches'; ul.replaceChildren(li); }
    sel = hits.length ? 0 : -1; mark();
    ul.hidden = false; q.setAttribute('aria-expanded', 'true');
  };
  const mark = () => [...ul.children].forEach((li, j) => li.setAttribute('aria-selected', String(j === sel)));
  const pick = (r) => { close(); q.value = ''; q.blur(); go(r.id); };
  q.addEventListener('input', show);
  q.addEventListener('blur', () => setTimeout(close, 100));
  q.addEventListener('keydown', (e) => {
    if (e.key === 'ArrowDown' && hits.length) { sel = (sel + 1) % hits.length; mark(); e.preventDefault(); }
    else if (e.key === 'ArrowUp' && hits.length) { sel = (sel - 1 + hits.length) % hits.length; mark(); e.preventDefault(); }
    else if (e.key === 'Enter' && sel >= 0 && hits[sel]) pick(hits[sel]);
    else if (e.key === 'Escape') close();
  });
}

/* ═══════════════ routing: #<id>[.M<m>][.g<d>][.n<lo>-<hi>] ═══════════════ */

function parseHash() {
  const parts = decodeURIComponent(location.hash.slice(1)).split('.');
  const o = { id: parts[0] || 'primes' };
  for (const p of parts.slice(1)) {
    let m;
    if ((m = /^M(\d+)$/.exec(p)) && +m[1] >= 1 && +m[1] <= OWN_MAX) o.M = +m[1];
    else if ((m = /^g(\d+)$/.exec(p))) o.g = +m[1];
    else if ((m = /^own(\d+)$/.exec(p)) && +m[1] >= 1 && +m[1] <= OWN_MAX && !LADDER.includes(+m[1])) o.own = +m[1];
    else if ((m = /^n(\d+)-(\d+)$/.exec(p))) o.n = [+m[1], +m[2]];
  }
  return o;
}
function writeHash() {
  if (!st.rec || st.rec.custom) return;
  let h = `#${st.rec.id}`;
  if (st.M !== 1) h += `.M${st.M}`;
  if (st.ownM && st.ownM !== st.M) h += `.own${st.ownM}`;
  if (st.g !== null) h += `.g${st.g}`;
  if (st.rLo > 0 || st.rHi < st.N - 1) h += `.n${st.n0 + st.rLo}-${st.n0 + st.rHi}`;
  if (h !== location.hash) history.replaceState(null, '', h);
}
function go(id) { location.hash = `#${id}`; }

/* ═══════════════ loading a sequence ═══════════════ */

async function loadSequence(o) {
  const r = st.byId.get(o.id) || st.byA.get(o.id.toUpperCase());
  if (!r) { setStatus(`no sequence “${o.id}” in the catalogue`, 'bad'); return; }
  const gen = ++st.gen;
  stopWorkers();
  st.dec = null; st.sub = null;
  st.rec = r; st.N = r.terms; st.n0 = r.n0;
  st.a = new Float64Array(st.N); st.d = new Float64Array(st.N); st.k = new Float64Array(st.N); st.L = new Float64Array(st.N);
  header();
  setStatus('loading…'); setProgress(0);
  let next = 0, got = 0;
  const one = async () => {
    for (;;) {
      const c = next++;
      if (c >= r.chunks) return;
      const url = `${DATA}/seq/${r.id}/chunk-${String(c).padStart(3, '0')}.bin.gz`;
      const bytes = await fetchBytes(url);
      if (gen !== st.gen) return;
      const at = c * r.chunk_rows, want = Math.min(r.chunk_rows, r.terms - at);
      if (readChunk(bytes, at, url) !== want) throw new Error(`${url}: wrong row count`);
      setProgress(0.15 * (++got / r.chunks));
    }
  };
  try {
    await Promise.all(Array.from({ length: Math.min(4, r.chunks) }, one));
    if (gen !== st.gen) return;
    for (let c = 1; c < r.chunks; c++) {
      const i = c * r.chunk_rows;
      if (st.a[i] !== st.a[i - 1] + st.d[i - 1]) throw new Error(`chunk ${c} does not continue chunk ${c - 1}`);
    }
    const bad = checkCounts(r);
    if (bad) throw new Error(`the data disagree with the catalogue: ${bad}`);
  } catch (e) {
    if (gen === st.gen) setStatus(`could not load ${r.anumber}: ${e.message}`, 'bad');
    return;
  }
  setStatus(`${fmtInt(r.terms)} terms, counts verified; checking the fingerprint…`, 'ok');
  afterLoad(o);
  verifyFingerprint(r, gen);
}

function checkCounts(r) {
  const c = { decomposable: 0, level: 0, ties: 0, level_one: 0, forced: 0 };
  for (let i = 0; i < st.N; i++) {
    const k = st.k[i]; if (!k) continue;
    const L = st.L[i], d = st.d[i];
    c.decomposable++; if (k > L) c.level++; if (k === L) c.ties++; if (L === 1) c.level_one++;
    if (d > FORCED_D || st.a[i] - d <= d * d) c.forced++;
  }
  const off = Object.keys(c).filter(key => r[key] !== c[key]);
  return off.length ? off.map(key => `${key} ${c[key]} ≠ ${r[key]}`).join(', ') : '';
}

let sums = null;
async function verifyFingerprint(r, gen) {
  try {
    if (!sums) sums = fetch(`${DATA}/sha256.txt`, { cache: 'no-cache' }).then(x => { if (!x.ok) throw new Error('sha256.txt'); return x.text(); })
      .then(t => new Map(t.split('\n').filter(Boolean).map(l => { const [h, f] = l.trim().split(/\s+/); return [f.slice(10, -4), h]; })));
    const want = (await sums).get(r.anumber);
    await new Promise(res => (window.requestIdleCallback || setTimeout)(res, { timeout: 1500 }));
    if (gen !== st.gen) return;
    const parts = ['n;a;weight;level;jump\n']; let buf = '';
    for (let i = 0; i < st.N; i++) {
      const k = st.k[i];
      buf += `${st.n0 + i};${st.a[i]};${k ? k : ''};${k ? st.L[i] : ''};${st.d[i]}\n`;
      if (buf.length > 65536) { parts.push(buf); buf = ''; }
    }
    parts.push(buf);
    const h = [...new Uint8Array(await crypto.subtle.digest('SHA-256', await new Blob(parts).arrayBuffer()))]
      .map(b => b.toString(16).padStart(2, '0')).join('');
    if (gen !== st.gen) return;
    st.verified = h === want ? `✓ data verified: SHA-256 ${h.slice(0, 16)}…, as published` : `⚠ the data differ from the published fingerprint`;
    st.verifiedOk = h === want;
  } catch {
    st.verified = 'counts verified (the fingerprint could not be checked here)'; st.verifiedOk = true;
  }
  if (gen === st.gen) refreshStatus();
}

function refreshStatus() {
  const left = st.total ? st.total - st.done : 0;
  const v = st.verified || 'counts verified';
  const c = left > 0 ? ` · controls ${Math.floor(100 * st.done / st.total)} %` : st.total ? ` · ${st.ctrl.size} controls on ${fmtInt(st.sub.length)} records` : '';
  setStatus(v + c, st.verifiedOk === false ? 'bad' : 'ok');
}

function header() {
  const r = st.rec;
  $('#seqName').textContent = r.alias || r.name;
  const sub = $('#seqSub'); sub.replaceChildren();
  if (r.custom) {
    sub.textContent = `your sequence · ${fmtInt(r.terms)} terms · not stored anywhere`;
    $('#viewerLink').hidden = $('#pageLink').hidden = true;
  } else {
    const a = document.createElement('a'); a.href = `https://oeis.org/${r.anumber}`; a.target = '_blank'; a.rel = 'noopener'; a.textContent = r.anumber;
    sub.append(a, ` · ${r.family} · ${fmtInt(r.terms)} terms from n = ${r.n0} · a ≤ ${fmtNum(r.amax)}`);
    $('#viewerLink').hidden = $('#pageLink').hidden = false;
    $('#viewerLink').href = `../#${r.id}`;
    $('#pageLink').href = `../seq/${r.anumber}/`;
  }
  document.title = `${r.alias || r.name} · decompwlj Lab`;
}

/** After the columns are in place: the records, the subsample, the filters, the controls. */
function afterLoad(o) {
  const dec = [];
  for (let i = 0; i < st.N; i++) if (st.k[i]) dec.push(i);
  st.dec = Int32Array.from(dec);
  st.M = o.M || 1;
  st.ownM = o.M && !LADDER.includes(o.M) ? o.M : o.own || null;
  st.g = o.g !== undefined ? o.g : null;
  st.rLo = 0; st.rHi = st.N - 1;
  if (o.n) { st.rLo = Math.max(0, o.n[0] - st.n0); st.rHi = Math.min(st.N - 1, o.n[1] - st.n0); if (st.rHi < st.rLo) { st.rLo = 0; st.rHi = st.N - 1; } }
  st.verified = null; st.verifiedOk = null;
  startControls();
}

function drawSubsample() {
  const n = st.dec.length, S = Math.min(st.sample, n);
  if (S === n) return st.dec.slice();
  const R = DWLJ.rng(DWLJ.seedOf(`${st.rec.id}|sample|${S}`));
  const idx = st.dec.slice();
  for (let i = 0; i < S; i++) { const j = i + R.below(n - i); const t = idx[i]; idx[i] = idx[j]; idx[j] = t; }
  return idx.subarray(0, S).sort();
}

/* ═══════════════ workers ═══════════════ */

function stopWorkers() {
  st.cgen++;
  for (const w of st.workers) w.terminate();
  st.workers = []; st.ctrl.clear(); st.done = 0; st.total = 0;
  st.queue = []; st.queued = new Set(); st.progress = new Map();
}

/** the moduli shown: the ladder, then your own M if it is not on it */
const steps = () => (st.ownM && !LADDER.includes(st.ownM) ? [...LADDER, st.ownM] : LADDER);

function startControls() {
  stopWorkers();
  st.sub = drawSubsample();
  const S = st.sub.length;
  st.inL = new Float64Array(S); st.inD = new Float64Array(S);
  for (let j = 0; j < S; j++) { const r = st.sub[j]; st.inD[j] = st.d[r]; st.inL[j] = st.a[r] - st.d[r]; }
  st.queue = []; st.queued = new Set(); st.progress = new Map();
  queueControls(steps());
  setProgress(0.15);
  schedule();
}

/** run the controls at these moduli (both replicates) unless already done or queued;
 *  the selected M goes first, so its charts fill in first */
function queueControls(Ms) {
  const tasks = [];
  for (let rep = 0; rep < REPS; rep++) for (const M of Ms) {
    const key = `${M}|${rep}`;
    if (!st.queued.has(key)) { st.queued.add(key); tasks.push({ M, rep, key }); }
  }
  if (!tasks.length) return;
  tasks.sort((x, y) => (x.rep - y.rep) || ((y.M === st.M) - (x.M === st.M)));
  st.queue.push(...tasks);
  st.total += tasks.length * st.inL.length;
  const W = Math.max(1, Math.min(4, (navigator.hardwareConcurrency || 2) - 1));
  while (st.workers.length < W && st.queue.length) spawn();
}

function spawn() {
  const gen = st.cgen, S = st.inL.length;
  const w = new Worker('worker.js');
  const feed = () => {
    if (!st.queue.length) { w.terminate(); st.workers = st.workers.filter(x => x !== w); return; }
    const t = st.queue.shift();
    st.progress.set(t.key, 0);
    w.postMessage({ type: 'control', task: t.key, l: st.inL, d: st.inD, M: t.M, seed: DWLJ.seedOf(`${st.rec.id}|M${t.M}|r${t.rep}|${S}`) });
  };
  w.onmessage = (e) => {
    if (gen !== st.cgen) return;
    const m = e.data;
    if (m.type === 'progress') {
      st.done += m.done - st.progress.get(m.task); st.progress.set(m.task, m.done);
      setProgress(0.15 + 0.85 * st.done / st.total); refreshStatus();
    } else if (m.type === 'control') {
      st.done += S - st.progress.get(m.task); st.progress.set(m.task, S);
      st.ctrl.set(m.task, { lp: m.lp, k: m.k, frozen: m.frozen });
      setProgress(0.15 + 0.85 * st.done / st.total); refreshStatus();
      schedule();
      feed();
    } else if (m.type === 'error') {
      setStatus(`control failed: ${m.message}`, 'bad');
    }
  };
  w.onerror = (e) => { if (gen === st.cgen) setStatus(`a worker failed: ${e.message || 'unknown error'}`, 'bad'); };
  st.workers.push(w);
  feed();
}

/** your own modulus: shown after the ladder, selected, computed on demand */
function setOwnM(M) {
  if (!(Number.isInteger(M) && M >= 1 && M <= OWN_MAX)) return false;
  if (!LADDER.includes(M)) st.ownM = M;
  st.M = M;
  if (st.sub) queueControls([M]);
  schedule();
  return true;
}

/* ═══════════════ statistics ═══════════════ */

const inN = (r) => r >= st.rLo && r <= st.rHi;
const inG = (r) => st.g === null || st.d[r] === st.g;

/** class counts of the observed records `rows` under both filters */
function obsCounts(rows) {
  const c = { n: 0, level: 0, weight: 0, ties: 0, one: 0, forced: 0 };
  for (let j = 0; j < rows.length; j++) {
    const r = rows[j];
    if (!inN(r) || !inG(r)) continue;
    const k = st.k[r], L = st.L[r], d = st.d[r];
    c.n++; if (k > L) c.level++; else c.weight++;
    if (k === L) c.ties++; if (L === 1) c.one++;
    if (d > FORCED_D || st.a[r] - d <= d * d) c.forced++;
  }
  return c;
}

/** the same counts for one control replicate, over the subsample */
function ctrlCounts(ct) {
  const c = { n: 0, level: 0, weight: 0, ties: 0, one: 0, forced: 0 };
  for (let j = 0; j < st.sub.length; j++) {
    const r = st.sub[j];
    if (!inN(r) || !inG(r)) continue;
    const lp = ct.lp[j], k = ct.k[j], L = lp / k, d = st.d[r];
    c.n++; if (k > L) c.level++; else c.weight++;
    if (k === L) c.ties++; if (L === 1) c.one++;
    if (d > FORCED_D || lp <= d * d) c.forced++;
  }
  return c;
}

const reps = (M) => { const o = []; for (let rep = 0; rep < REPS; rep++) { const c = st.ctrl.get(`${M}|${rep}`); if (c) o.push(c); } return o; };

/** Δ = observed − mean control on the paired records, with σ of a single replicate's share
 *  (binomial bound, or the replicates' spread if larger), inflated by √(1 + 1/R). */
function ladderRow(M, obsSub) {
  const rs = reps(M);
  if (!rs.length || !obsSub.n) return null;
  const cs = rs.map(ctrlCounts);
  const row = { M, R: rs.length, frozen: rs[0].frozen / st.sub.length };
  for (const key of ['level', 'one']) {
    const ps = cs.map(c => c[key] / c.n), p = ps.reduce((s, x) => s + x, 0) / ps.length;
    let sd = Math.sqrt(p * (1 - p) / cs[0].n);
    if (ps.length > 1) { const v = ps.reduce((s, x) => s + (x - p) ** 2, 0) / (ps.length - 1); sd = Math.max(sd, Math.sqrt(v)); }
    const s = sd * Math.sqrt(1 + 1 / ps.length), dlt = obsSub[key] / obsSub.n - p;
    row[key] = { p, delta: dlt, sigma: s, z: s > 0 ? dlt / s : (dlt === 0 ? 0 : Infinity) };
  }
  return row;
}

/* ═══════════════ drawing ═══════════════ */

function setup(cv) {
  const box = cv.getBoundingClientRect(), dpr = Math.min(2, window.devicePixelRatio || 1);
  const w = Math.max(20, Math.round(box.width)), h = Math.max(20, Math.round(box.height));
  if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
  const ctx = cv.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, w, h);
  ctx.font = `11px ${css('--mono') || 'monospace'}`; ctx.textBaseline = 'middle';
  return { ctx, w, h };
}

function ticks(lo, hi, n = 5) {
  const span = hi - lo || 1, step0 = span / n, mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 2.5, 5, 10].map(f => f * mag).find(s => s >= step0) || 10 * mag;
  const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toPrecision(12));
  return out;
}

/** frame: axes box with recessive grid; returns scale functions */
function frame(ctx, box, xr, yr, o = {}) {
  const { x0, y0, x1, y1 } = box;
  const X = (v) => x0 + (v - xr[0]) / (xr[1] - xr[0]) * (x1 - x0);
  const Y = (v) => y1 - (v - yr[0]) / (yr[1] - yr[0]) * (y1 - y0);
  ctx.lineWidth = 1; ctx.strokeStyle = css('--line'); ctx.fillStyle = css('--fg-3');
  for (const v of o.yt || ticks(yr[0], yr[1])) {
    const y = Math.round(Y(v)) + 0.5; ctx.beginPath(); ctx.moveTo(x0, y); ctx.lineTo(x1, y); ctx.stroke();
    ctx.textAlign = 'right'; ctx.fillText(o.yf ? o.yf(v) : String(v), x0 - 6, y);
  }
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  let right = -Infinity;
  for (const v of o.xt || ticks(xr[0], xr[1])) {
    const x = X(v); if (o.xgrid) { ctx.beginPath(); ctx.moveTo(Math.round(x) + 0.5, y0); ctx.lineTo(Math.round(x) + 0.5, y1); ctx.stroke(); }
    const t = o.xf ? o.xf(v) : String(v), tw = ctx.measureText(t).width;
    const cx = Math.min(Math.max(x, x0 + tw / 2), x1 - tw / 2);       /* kept inside the frame */
    if (cx - tw / 2 < right + 5) continue;                             /* and never on top of the last one */
    ctx.fillText(t, cx, y1 + 6); right = cx + tw / 2;
  }
  ctx.textBaseline = 'middle';
  if (o.xl) { ctx.textAlign = 'right'; ctx.textBaseline = 'bottom'; ctx.fillText(o.xl, x1, box.h); ctx.textBaseline = 'middle'; }
  if (o.yl) { ctx.save(); ctx.translate(11, y0 + (y1 - y0) / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'center'; ctx.fillText(o.yl, 0, 0); ctx.restore(); }
  return { X, Y };
}

function message(ctx, w, h, t) {
  ctx.fillStyle = css('--fg-3'); ctx.textAlign = 'center'; ctx.fillText(t, w / 2, h / 2);
}

/* hover targets per canvas: [{x0,y0,x1,y1,text}] or a function (x, y) → text */
const hits = new Map();

/* ── the two plates ── */

function plateRange() {
  let m = 1;
  for (let j = 0; j < st.sub.length; j++) { const r = st.sub[j]; m = Math.max(m, st.k[r], st.L[r]); }
  for (const c of reps(st.M).slice(0, 1)) for (let j = 0; j < c.k.length; j++) m = Math.max(m, c.k[j], c.lp[j] / c.k[j]);
  return Math.max(1, Math.ceil(Math.log10(m) * 4) / 4);
}

function drawPlate(cv, which) {
  const { ctx, w, h } = setup(cv);
  if (!st.sub) return;
  const T = plateRange();
  const ml = 40, mb = 34, mt = 8, mr = 8, s = Math.min(w - ml - mr, h - mt - mb);
  const box = { x0: ml, y0: mt, x1: ml + s, y1: mt + s, h };
  const tk = ticks(0, T, 5);
  const pad = 0.025 * T;                     /* the same on both axes: the aspect stays equal */
  const { X, Y } = frame(ctx, box, [-pad, T], [-pad, T], { xt: tk, yt: tk, xgrid: true, xl: 'log₁₀ k', yl: 'log₁₀ L' });
  const ct = which === 'B' ? reps(st.M)[0] : null;
  if (which === 'B' && !ct) { ctx.fillStyle = css('--fg-3'); ctx.textAlign = 'center'; ctx.fillText(`control M = ${st.M}: computing…`, ml + s / 2, mt + s / 2); return; }
  /* guides: k·L = 10^j, then the diagonal */
  ctx.save(); ctx.beginPath(); ctx.rect(box.x0, box.y0, s, s); ctx.clip();
  ctx.setLineDash([3, 4]); ctx.strokeStyle = css('--line-2');
  for (let j = 1; j <= 2 * T; j++) { ctx.beginPath(); ctx.moveTo(X(-pad), Y(j + pad)); ctx.lineTo(X(j + pad), Y(-pad)); ctx.stroke(); }
  ctx.setLineDash([]); ctx.strokeStyle = css('--fg-3'); ctx.beginPath(); ctx.moveTo(X(-pad), Y(-pad)); ctx.lineTo(X(T), Y(T)); ctx.stroke();
  const cw = css('--c-weight'), cl = css('--c-level'), cf = css('--fg-3');
  const n = st.sub.length, ps = n > 40000 ? 1.4 : n > 15000 ? 1.8 : 2.4;
  const pass = (filtered) => {
    for (let j = 0; j < n; j++) {
      const r = st.sub[j], on = inN(r) && inG(r);
      if (on !== filtered) continue;
      const k = ct ? ct.k[j] : st.k[r], L = ct ? ct.lp[j] / ct.k[j] : st.L[r];
      ctx.fillStyle = !on ? cf : k > L ? cl : cw;
      ctx.fillRect(X(Math.log10(k)) - ps / 2, Y(Math.log10(L)) - ps / 2, ps, ps);
    }
  };
  ctx.globalAlpha = 0.12; pass(false);
  ctx.globalAlpha = n > 40000 ? 0.45 : 0.65; pass(true);
  ctx.globalAlpha = 1; ctx.restore();
  /* hover: the nearest record within 8 px, through a coarse grid */
  const G = 48, cell = new Map();
  for (let j = 0; j < n; j++) {
    const r = st.sub[j]; if (!inN(r) || !inG(r)) continue;
    const k = ct ? ct.k[j] : st.k[r], L = ct ? ct.lp[j] / ct.k[j] : st.L[r];
    const x = X(Math.log10(k)), y = Y(Math.log10(L)), key = Math.floor(x / G) * 1e5 + Math.floor(y / G);
    let a = cell.get(key); if (!a) cell.set(key, a = []); a.push([x, y, j]);
  }
  hits.set(cv, (px, py) => {
    let best = null, bd = 64;
    for (let gx = -1; gx <= 1; gx++) for (let gy = -1; gy <= 1; gy++) {
      for (const [x, y, j] of cell.get((Math.floor(px / G) + gx) * 1e5 + Math.floor(py / G) + gy) || []) {
        const dd = (x - px) ** 2 + (y - py) ** 2; if (dd < bd) { bd = dd; best = j; }
      }
    }
    if (best === null) return null;
    const r = st.sub[best], d = st.d[r];
    if (ct) {
      const lp = ct.lp[best], k = ct.k[best], L = lp / k;
      return [`<b>control, M = ${st.M}</b>  record n = ${st.n0 + r}`, `d = ${fmtInt(d)}   ℓ = ${fmtInt(st.a[r] - d)} → ℓ′ = ${fmtInt(lp)}`,
              `k = ${fmtInt(k)}   L = ${fmtInt(L)}   <span class="m">${k > L ? 'level' : k === L ? 'tie (weight)' : 'weight'}</span>`];
    }
    const k = st.k[r], L = st.L[r];
    return [`<b>n = ${st.n0 + r}</b>   a = ${fmtInt(st.a[r])}`, `${fmtInt(st.a[r])} = ${fmtInt(k)} × ${fmtInt(L)} + ${fmtInt(d)}`,
            `<span class="m">${k > L ? 'level' : k === L ? 'tie (weight)' : 'weight'}</span>`];
  });
}

/* ── the ladder ── */

let ladderRows = [];
function computeLadder() {
  const obsSub = obsCounts(st.sub || []);
  ladderRows = steps().map(M => ladderRow(M, obsSub));
}

function drawLadder() {
  const cv = $('#ladder'), { ctx, w, h } = setup(cv);
  const Ms = steps(), nL = LADDER.length;
  const rows = ladderRows.filter(Boolean);
  const weak = (r) => r && r.frozen > FROZEN_MAX;
  const tb = $('#ladderTable tbody');
  tb.replaceChildren(...Ms.map((M, i) => {
    const r = ladderRows[i], tr = document.createElement('tr');
    if (M === st.M) tr.className = 'sel';
    if (weak(r)) tr.classList.add('weak');
    const adds = i >= nL ? 'yours' : i === 0 ? '' : `+${ADDS[M]}`;
    const cells = r ? [fmtInt(M), adds, pts(r.level.delta), fz(r.level.z), pts(r.one.delta), fz(r.one.z), pct(r.frozen, 1)]
                    : [fmtInt(M), adds, '…', '', '…', '', ''];
    cells.forEach((t, j) => { const td = document.createElement('td'); td.textContent = t; if (j === 3 || j === 5) td.className = 'z'; tr.append(td); });
    tr.addEventListener('click', () => setM(M));
    return tr;
  }));
  if (!rows.length) { message(ctx, w, h, 'controls: computing…'); hits.delete(cv); return; }
  let lo = 0, hi = 0;
  for (const r of rows) for (const key of ['level', 'one']) { lo = Math.min(lo, r[key].delta - 2 * r[key].sigma); hi = Math.max(hi, r[key].delta + 2 * r[key].sigma); }
  const pad = (hi - lo) * 0.08 || 0.005; lo -= pad; hi += pad;
  const box = { x0: 58, y0: 8, x1: w - 10, y1: h - 34, h };
  /* your own M sits after a gap */
  const xs = (i) => (i < nL ? i : i + 0.6);
  const xmax = xs(Ms.length - 1) + 0.5;
  const crowded = (box.x1 - box.x0) / (xmax + 0.5) < 46;
  const { X, Y } = frame(ctx, box, [-0.5, xmax], [lo, hi], {
    xt: Ms.map((_, i) => xs(i)), xf: (x) => { const M = Ms[x < nL ? x : Math.round(x - 0.6)]; return !crowded || M < 1000 ? fmtInt(M) : M < 10000 ? `${(M / 1000).toFixed(1)}k` : `${Math.round(M / 1000)}k`; },
    yf: (v) => signed(100 * v, Math.abs(hi - lo) < 0.02 ? 2 : 1), xl: 'M', yl: 'Δ, points' });
  ctx.strokeStyle = css('--fg-3'); ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(box.x0, Math.round(Y(0)) + 0.5); ctx.lineTo(box.x1, Math.round(Y(0)) + 0.5); ctx.stroke();
  if (Ms.length > nL) { const xg = X(nL - 0.5 + 0.3); ctx.strokeStyle = css('--line-2'); ctx.setLineDash([2, 3]); ctx.beginPath(); ctx.moveTo(xg, box.y0); ctx.lineTo(xg, box.y1); ctx.stroke(); ctx.setLineDash([]); }
  const sel = Ms.indexOf(st.M);
  if (sel >= 0) { ctx.fillStyle = css('--accent'); ctx.globalAlpha = 0.1; ctx.fillRect(X(xs(sel) - 0.5), box.y0, X(xs(sel) + 0.5) - X(xs(sel) - 0.5), box.y1 - box.y0); ctx.globalAlpha = 1; }
  const targets = [];
  for (const [key, col, off] of [['level', css('--c-level'), -0.12], ['one', css('--c-one'), 0.12]]) {
    /* the line joins the informative ladder steps only */
    ctx.strokeStyle = col; ctx.lineWidth = 2; ctx.beginPath();
    let first = true;
    for (let i = 0; i < nL; i++) { const r = ladderRows[i]; if (!r || weak(r)) continue; const x = X(xs(i) + off), y = Y(r[key].delta); if (first) ctx.moveTo(x, y); else ctx.lineTo(x, y); first = false; }
    ctx.globalAlpha = 0.5; ctx.stroke(); ctx.globalAlpha = 1;
    Ms.forEach((M, i) => {
      const r = ladderRows[i]; if (!r) return;
      const x = X(xs(i) + off), v = r[key];
      ctx.globalAlpha = weak(r) ? 0.45 : 1;
      ctx.strokeStyle = col; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(x, Y(v.delta - 2 * v.sigma)); ctx.lineTo(x, Y(v.delta + 2 * v.sigma)); ctx.stroke();
      ctx.beginPath(); ctx.arc(x, Y(v.delta), 4, 0, 2 * Math.PI);
      if (weak(r)) { ctx.fillStyle = css('--bg-2'); ctx.fill(); ctx.strokeStyle = col; ctx.lineWidth = 1.5; ctx.stroke(); }   /* hollow: mostly frozen */
      else { ctx.fillStyle = col; ctx.fill(); ctx.strokeStyle = css('--bg-2'); ctx.lineWidth = 2; ctx.stroke(); }
      ctx.globalAlpha = 1;
    });
  }
  Ms.forEach((M, i) => {
    const r = ladderRows[i];
    const head = `<b>M = ${fmtInt(M)}</b>  <span class="m">${i >= nL ? 'your own modulus' : i ? 'adds ' + ADDS[M] : 'keeps no residue'}</span>`;
    targets.push({ x0: X(xs(i) - 0.5), x1: X(xs(i) + 0.5), y0: box.y0, y1: box.y1, M,
      text: r ? [head,
                 `level  Δ ${pts(r.level.delta)} pts   z ${fz(r.level.z)}`, `       control ${pct(r.level.p)}`,
                 `L = 1  Δ ${pts(r.one.delta)} pts   z ${fz(r.one.z)}`, `       control ${pct(r.one.p)}`,
                 `frozen ${pct(r.frozen, 1)} of the records${weak(r) ? ': Δ is pushed toward 0,\n       this step says little' : ''}`,
                 `<span class="m">${r.R} replicate${r.R > 1 ? 's' : ''} · click: use this control</span>`] : [head, 'computing…'] });
  });
  hits.set(cv, targets);
}

/* ── level share along n ── */

function drawAlong() {
  const cv = $('#along'), { ctx, w, h } = setup(cv);
  if (!st.dec) return;
  const lo = 0, hi = st.N, bw = (hi - lo) / BINS;
  const obs = Array.from({ length: BINS }, () => [0, 0]), ctl = Array.from({ length: BINS }, () => [0, 0]);
  for (let j = 0; j < st.dec.length; j++) { const r = st.dec[j]; if (!inG(r)) continue; const b = Math.min(BINS - 1, Math.floor((r - lo) / bw)); obs[b][0]++; if (st.k[r] > st.L[r]) obs[b][1]++; }
  for (const c of reps(st.M)) for (let j = 0; j < st.sub.length; j++) {
    const r = st.sub[j]; if (!inG(r)) continue; const b = Math.min(BINS - 1, Math.floor((r - lo) / bw));
    ctl[b][0]++; if (c.k[j] * c.k[j] > c.lp[j]) ctl[b][1]++;
  }
  let ymax = 0.05;
  const P = obs.map(([n, l]) => n ? l / n : null), C = ctl.map(([n, l]) => n ? { p: l / n, s: Math.sqrt((l / n) * (1 - l / n) / (n / Math.max(1, reps(st.M).length))) } : null);
  for (const p of P) if (p !== null) ymax = Math.max(ymax, p);
  for (const c of C) if (c) ymax = Math.max(ymax, c.p + 2 * c.s);
  ymax = Math.min(1, ymax * 1.08);
  const box = { x0: 56, y0: 8, x1: w - 10, y1: h - 34, h };
  const nx = (r) => st.n0 + r;
  const { X, Y } = frame(ctx, box, [nx(lo), nx(hi)], [0, ymax], { yf: (v) => `${Math.round(100 * v)} %`, xf: (v) => fmtNum(v), xl: 'n', yl: 'level share' });
  /* the kept range of n */
  if (st.rLo > 0 || st.rHi < st.N - 1) {
    ctx.fillStyle = css('--accent'); ctx.globalAlpha = 0.1;
    ctx.fillRect(X(nx(st.rLo)), box.y0, X(nx(st.rHi + 1)) - X(nx(st.rLo)), box.y1 - box.y0); ctx.globalAlpha = 1;
  }
  const xc = (b) => X(nx(lo + (b + 0.5) * bw));
  /* control band */
  ctx.fillStyle = css('--c-ctrl'); ctx.globalAlpha = 0.3; ctx.beginPath();
  let started = false;
  for (let b = 0; b < BINS; b++) if (C[b]) { const y = Y(Math.min(1, C[b].p + 2 * C[b].s)); if (!started) { ctx.moveTo(xc(b), y); started = true; } else ctx.lineTo(xc(b), y); }
  for (let b = BINS - 1; b >= 0; b--) if (C[b]) ctx.lineTo(xc(b), Y(Math.max(0, C[b].p - 2 * C[b].s)));
  if (started) ctx.fill();
  ctx.globalAlpha = 1;
  if (started) { ctx.strokeStyle = css('--c-ctrl'); ctx.lineWidth = 1; ctx.beginPath(); let f = true; for (let b = 0; b < BINS; b++) if (C[b]) { if (f) ctx.moveTo(xc(b), Y(C[b].p)); else ctx.lineTo(xc(b), Y(C[b].p)); f = false; } ctx.stroke(); }
  else { ctx.fillStyle = css('--fg-3'); ctx.textAlign = 'center'; ctx.fillText('control: computing…', (box.x0 + box.x1) / 2, box.y0 + 12); }
  /* the sequence */
  ctx.strokeStyle = css('--c-level'); ctx.lineWidth = 2; ctx.beginPath();
  let f = true; for (let b = 0; b < BINS; b++) if (P[b] !== null) { if (f) ctx.moveTo(xc(b), Y(P[b])); else ctx.lineTo(xc(b), Y(P[b])); f = false; }
  ctx.stroke();
  ctx.fillStyle = css('--c-level');
  for (let b = 0; b < BINS; b++) if (P[b] !== null) { ctx.beginPath(); ctx.arc(xc(b), Y(P[b]), 2.5, 0, 2 * Math.PI); ctx.fill(); }
  hits.set(cv, (px, py) => {
    if (px < box.x0 || px > box.x1) return null;
    const b = Math.max(0, Math.min(BINS - 1, Math.floor(((px - box.x0) / (box.x1 - box.x0)) * BINS)));
    const a0 = nx(Math.floor(lo + b * bw)), a1 = nx(Math.ceil(lo + (b + 1) * bw) - 1);
    return [`<b>n ${fmtInt(a0)} – ${fmtInt(a1)}</b>`,
            `sequence  ${P[b] === null ? '—' : pct(P[b])}   <span class="m">${fmtInt(obs[b][0])} records</span>`,
            `control   ${C[b] ? pct(C[b].p) + ' ± ' + (200 * C[b].s).toFixed(2) : '—'}`,
            `<span class="m">drag to keep a range of n</span>`];
  });
  cv._map = (px) => { const f = (px - box.x0) / (box.x1 - box.x0); return Math.round(lo + Math.max(0, Math.min(1, f)) * (hi - lo)); };
  cv._box = box;
}

/* ── jumps ── */

function drawGaps() {
  const cv = $('#gaps'), { ctx, w, h } = setup(cv);
  if (!st.dec) return;
  const cnt = new Map();
  for (let j = 0; j < st.dec.length; j++) { const r = st.dec[j]; if (!inN(r)) continue; const d = st.d[r]; cnt.set(d, (cnt.get(d) || 0) + 1); }
  let top = [...cnt.entries()].sort((x, y) => y[1] - x[1]).slice(0, 24);
  if (st.g !== null && !top.some(([d]) => d === st.g)) top.push([st.g, cnt.get(st.g) || 0]);
  top.sort((x, y) => x[0] - y[0]);
  const total = [...cnt.values()].reduce((s, x) => s + x, 0);
  if (!top.length) { message(ctx, w, h, 'no records in this range'); return; }
  const ymax = Math.max(...top.map(t => t[1] / total)) * 1.1;
  const box = { x0: 56, y0: 8, x1: w - 8, y1: h - 34, h };
  const n = top.length, bw = (box.x1 - box.x0) / n;
  const { Y } = frame(ctx, box, [0, n], [0, ymax], { xt: [], yf: (v) => `${Math.round(100 * v)} %`, yl: 'share' });
  ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillStyle = css('--fg-3');
  const every = Math.ceil(n / Math.max(1, Math.floor((box.x1 - box.x0) / 34)));
  top.forEach(([d], i) => { if (i % every === 0) ctx.fillText(fmtNum(d), box.x0 + (i + 0.5) * bw, box.y1 + 6); });
  ctx.textBaseline = 'bottom'; ctx.textAlign = 'right'; ctx.fillText('jump d', box.x1, h); ctx.textBaseline = 'middle';
  const targets = [];
  top.forEach(([d, c], i) => {
    const x = box.x0 + i * bw + 1, bwid = Math.max(1, bw - 2), y = Y(c / total);
    ctx.fillStyle = st.g === d ? css('--accent') : css('--fg-2'); ctx.globalAlpha = st.g === null || st.g === d ? 0.85 : 0.35;
    roundTop(ctx, x, y, bwid, box.y1 - y, Math.min(4, bwid / 2)); ctx.globalAlpha = 1;
    targets.push({ x0: box.x0 + i * bw, x1: box.x0 + (i + 1) * bw, y0: box.y0, y1: box.y1, d,
      text: [`<b>d = ${fmtInt(d)}</b>`, `${fmtInt(c)} records   ${pct(c / total)}`, `<span class="m">${st.g === d ? 'click: all jumps' : 'click: keep this jump'}</span>`] });
  });
  hits.set(cv, targets);
}

function roundTop(ctx, x, y, w, h, r) {
  if (h <= 0) return; r = Math.min(r, h);
  ctx.beginPath(); ctx.moveTo(x, y + h); ctx.lineTo(x, y + r); ctx.quadraticCurveTo(x, y, x + r, y);
  ctx.lineTo(x + w - r, y); ctx.quadraticCurveTo(x + w, y, x + w, y + r); ctx.lineTo(x + w, y + h); ctx.closePath(); ctx.fill();
}

/* ── the congruence lock ── */

function lockStats() {
  const R = LOCK_R.length, on = new Float64Array(R), hit = new Float64Array(R), con = new Float64Array(R), chit = new Float64Array(R);
  for (let j = 0; j < st.dec.length; j++) {
    const r = st.dec[j]; if (!inN(r) || !inG(r)) continue;
    const d = st.d[r], l = st.a[r] - d;
    for (let i = 0; i < R; i++) { const p = LOCK_R[i]; if (d < p) { on[i]++; if (l % p === 0) hit[i]++; } }
  }
  for (const c of reps(st.M)) for (let j = 0; j < st.sub.length; j++) {
    const r = st.sub[j]; if (!inN(r) || !inG(r)) continue;
    const d = st.d[r], l = c.lp[j];
    for (let i = 0; i < R; i++) { const p = LOCK_R[i]; if (d < p) { con[i]++; if (l % p === 0) chit[i]++; } }
  }
  return LOCK_R.map((p, i) => ({
    r: p, n: on[i], y: on[i] ? p * hit[i] / on[i] : null, s: on[i] ? p * Math.sqrt((hit[i] / on[i]) * (1 - hit[i] / on[i]) / on[i]) : 0,
    cn: con[i], cy: con[i] ? p * chit[i] / con[i] : null, cs: con[i] ? p * Math.sqrt((chit[i] / con[i]) * (1 - chit[i] / con[i]) / con[i]) : 0,
  }));
}

function drawLock() {
  const cv = $('#lock'), { ctx, w, h } = setup(cv);
  if (!st.dec) return;
  const S = lockStats().filter(s => s.n >= 50);
  if (!S.length) { message(ctx, w, h, 'no prime r ≤ 113 exceeds the jumps here (d < r is needed)'); hits.delete(cv); return; }
  let lo = 0.5, hi = 2;
  for (const s of S) { lo = Math.min(lo, s.y - 2 * s.s); hi = Math.max(hi, s.y + 2 * s.s); if (s.cy !== null && s.cn >= 50) { lo = Math.min(lo, s.cy - 2 * s.cs); hi = Math.max(hi, s.cy + 2 * s.cs); } }
  lo = Math.max(0, lo);
  const box = { x0: 52, y0: 8, x1: w - 10, y1: h - 34, h };
  const { X, Y } = frame(ctx, box, [0, 116], [lo, hi * 1.04], { xt: [3, 13, 23, 37, 53, 71, 89, 113], xl: 'prime r', yl: 'r · Pr[r | ℓ]' });
  ctx.save(); ctx.beginPath(); ctx.rect(box.x0, box.y0, box.x1 - box.x0, box.y1 - box.y0); ctx.clip();
  ctx.strokeStyle = css('--fg-3'); ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(box.x0, Y(1)); ctx.lineTo(box.x1, Y(1)); ctx.stroke();
  ctx.strokeStyle = css('--c-ref'); ctx.lineWidth = 1.5; ctx.setLineDash([4, 3]); ctx.beginPath();
  for (let r = 3; r <= 116; r += 0.5) { const y = Y(r / (r - 2)); if (r === 3) ctx.moveTo(X(r), y); else ctx.lineTo(X(r), y); }
  ctx.stroke(); ctx.setLineDash([]);
  const targets = [];
  for (const s of S) {
    if (s.cy !== null && s.cn >= 50) {
      const x = X(s.r) + 3; ctx.strokeStyle = css('--c-ctrl'); ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.moveTo(x, Y(s.cy - 2 * s.cs)); ctx.lineTo(x, Y(s.cy + 2 * s.cs)); ctx.stroke();
      ctx.fillStyle = css('--c-ctrl'); ctx.beginPath(); ctx.arc(x, Y(s.cy), 3, 0, 2 * Math.PI); ctx.fill();
    }
    const x = X(s.r) - 2; ctx.strokeStyle = css('--accent'); ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(x, Y(s.y - 2 * s.s)); ctx.lineTo(x, Y(s.y + 2 * s.s)); ctx.stroke();
    ctx.fillStyle = css('--accent'); ctx.beginPath(); ctx.arc(x, Y(s.y), 4, 0, 2 * Math.PI); ctx.fill();
    ctx.strokeStyle = css('--bg-2'); ctx.lineWidth = 1.5; ctx.stroke();
    targets.push({ x0: X(s.r) - 9, x1: X(s.r) + 9, y0: box.y0, y1: box.y1,
      text: [`<b>r = ${s.r}</b>   <span class="m">${fmtInt(s.n)} records with d < r</span>`,
             `sequence  ${s.y.toFixed(3)} ± ${(2 * s.s).toFixed(3)}`,
             `control   ${s.cy === null ? '—' : s.cy.toFixed(3) + ' ± ' + (2 * s.cs).toFixed(3)}${st.M % s.r === 0 ? '  (kept: r | M)' : ''}`,
             `r/(r − 2) ${(s.r / (s.r - 2)).toFixed(3)}   chance 1`] });
  }
  ctx.restore();
  hits.set(cv, targets);
}

/* ── columns and lines ── */

function drawCols() {
  const cv = $('#cols'), { ctx, w, h } = setup(cv);
  if (!st.dec) return;
  const byK = st.colMode === 'k', off = byK ? 2 : 1;
  const o = new Float64Array(MAXC), c = new Float64Array(MAXC);
  let on = 0, cn = 0;
  for (let j = 0; j < st.dec.length; j++) {
    const r = st.dec[j]; if (!inN(r) || !inG(r)) continue;
    on++; const v = (byK ? st.k[r] : st.L[r]) - off; if (v >= 0 && v < MAXC) o[v]++;
  }
  for (const ct of reps(st.M)) for (let j = 0; j < st.sub.length; j++) {
    const r = st.sub[j]; if (!inN(r) || !inG(r)) continue;
    cn++; const v = (byK ? ct.k[j] : ct.lp[j] / ct.k[j]) - off; if (v >= 0 && v < MAXC) c[v]++;
  }
  if (!on) { message(ctx, w, h, 'no records under these filters'); return; }
  let ymax = 0;
  for (let i = 0; i < MAXC; i++) { ymax = Math.max(ymax, o[i] / on); if (cn) ymax = Math.max(ymax, c[i] / cn); }
  if (ymax === 0) { message(ctx, w, h, `no ${byK ? 'weight' : 'level'} below ${MAXC + off} here`); hits.delete(cv); return; }
  ymax *= 1.08;
  const box = { x0: 60, y0: 8, x1: w - 10, y1: h - 34, h };
  const { X, Y } = frame(ctx, box, [off - 0.5, off + MAXC - 0.5], [0, ymax],
    { xt: ticks(off, off + MAXC - 1, 8).filter(v => Number.isInteger(v)), yf: (v) => `${+(100 * v).toFixed(1)} %`, xl: byK ? 'weight k' : 'level L', yl: 'share' });
  const bw = X(1) - X(0), col = byK ? css('--c-weight') : css('--c-level');
  const targets = [];
  for (let i = 0; i < MAXC; i++) {
    const v = i + off, x = X(v) - bw / 2 + 1, y = Y(o[i] / on);
    ctx.fillStyle = col; ctx.globalAlpha = 0.85; roundTop(ctx, x, y, Math.max(1, bw - 2), box.y1 - y, Math.min(3, (bw - 2) / 2)); ctx.globalAlpha = 1;
    if (cn) { const yc = Y(c[i] / cn); ctx.strokeStyle = css('--c-ctrl'); ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(X(v) - bw / 2, yc); ctx.lineTo(X(v) + bw / 2, yc); ctx.stroke(); }
    targets.push({ x0: X(v) - bw / 2, x1: X(v) + bw / 2, y0: box.y0, y1: box.y1,
      text: [`<b>${byK ? 'k' : 'L'} = ${v}</b>`, `sequence ${pct(o[i] / on, 3)}   <span class="m">${fmtInt(o[i])}</span>`, `control  ${cn ? pct(c[i] / cn, 3) : '—'}`] });
  }
  hits.set(cv, targets);
}

/* ── the numbers ── */

function drawNumbers() {
  const tb = $('#numTable tbody');
  if (!st.dec) return;
  const rows = [];
  const add = (label, c, sel) => {
    const tr = document.createElement('tr'); if (sel) tr.className = 'sel';
    const v = c ? [fmtInt(c.n), pct(c.level / c.n), pct(c.weight / c.n), fmtInt(c.ties), pct(c.one / c.n), pct(c.forced / c.n)] : ['…', '', '', '', '', ''];
    [label, ...v].forEach(t => { const td = document.createElement('td'); td.textContent = t; tr.append(td); });
    rows.push(tr);
  };
  const all = obsCounts(st.dec);
  if (!all.n) { tb.replaceChildren(); return; }
  add('sequence, every record', all);
  add('sequence, the control’s records', obsCounts(st.sub));
  for (const M of steps()) {
    const rs = reps(M);
    if (!rs.length) { add(`control M = ${fmtInt(M)}`, null, M === st.M); continue; }
    const cs = rs.map(ctrlCounts), m = {};
    for (const key of ['n', 'level', 'weight', 'ties', 'one', 'forced']) m[key] = cs.reduce((s, c) => s + c[key], 0) / cs.length;
    m.ties = Math.round(m.ties * 10) / 10;
    add(`control M = ${fmtInt(M)}  (${rs.length} rep.)`, m, M === st.M);
  }
  tb.replaceChildren(...rows);
}

/* ── one redraw for everything, at most once a frame ── */

let raf = 0;
function schedule() { if (!raf) raf = requestAnimationFrame(() => { raf = 0; drawAll(); }); }
function drawAll() {
  if (!st.dec) return;
  computeLadder();
  drawPlate($('#plateA'), 'A'); drawPlate($('#plateB'), 'B');
  drawLadder(); drawAlong(); drawGaps(); drawLock(); drawCols(); drawNumbers();
  $('#fN').textContent = st.rLo > 0 || st.rHi < st.N - 1 ? `n = ${fmtInt(st.n0 + st.rLo)} – ${fmtInt(st.n0 + st.rHi)}` : 'all';
  $('#clearN').hidden = !(st.rLo > 0 || st.rHi < st.N - 1);
  $('#fG').textContent = st.g === null ? 'all' : `d = ${fmtInt(st.g)}`;
  $('#clearG').hidden = st.g === null;
  $('#fM').textContent = `M = ${fmtInt(st.M)}`;
  if (document.activeElement !== $('#ownM')) $('#ownM').value = st.ownM || '';
  $('#hPlateB').textContent = `The same records, controlled at M = ${fmtInt(st.M)}`;
  writeHash();
}

function setM(M) { if (M !== st.M) { st.M = M; schedule(); } }

/* ═══════════════ interaction ═══════════════ */

function initHover() {
  const tip = $('#tip');
  const hide = () => { tip.hidden = true; };
  const find = (cv, e) => {
    const h = hits.get(cv); if (!h) return null;
    const b = cv.getBoundingClientRect(), x = e.clientX - b.left, y = e.clientY - b.top;
    if (typeof h === 'function') return { text: h(x, y) };
    return h.find(t => x >= t.x0 && x < t.x1 && y >= t.y0 && y <= t.y1) || null;
  };
  for (const cv of document.querySelectorAll('canvas')) {
    cv.addEventListener('pointermove', (e) => {
      if (cv._drag) return;
      const t = find(cv, e);
      if (!t || !t.text) { hide(); cv.style.cursor = ''; return; }
      tip.innerHTML = t.text.join('\n');      /* built from numbers and fixed strings only */
      tip.hidden = false;
      const tw = tip.offsetWidth, th = tip.offsetHeight;
      let x = e.clientX + 14, y = e.clientY + 14;
      if (x + tw > innerWidth - 8) x = e.clientX - tw - 14;
      if (y + th > innerHeight - 8) y = e.clientY - th - 14;
      tip.style.left = `${Math.max(8, x)}px`; tip.style.top = `${Math.max(8, y)}px`;
      cv.style.cursor = t.M !== undefined || t.d !== undefined ? 'pointer' : '';
    });
    cv.addEventListener('pointerleave', hide);
  }
  $('#ladder').addEventListener('click', (e) => { const t = find($('#ladder'), e); if (t && t.M !== undefined) setM(t.M); });
  $('#gaps').addEventListener('click', (e) => { const t = find($('#gaps'), e); if (t && t.d !== undefined) { st.g = st.g === t.d ? null : t.d; schedule(); } });
  /* brushing along n */
  const al = $('#along');
  al.addEventListener('pointerdown', (e) => {
    if (!al._map) return;
    const b = al.getBoundingClientRect();
    al._drag = { x: e.clientX - b.left }; al.setPointerCapture(e.pointerId); hide();
  });
  al.addEventListener('pointermove', (e) => {
    if (!al._drag) return;
    const b = al.getBoundingClientRect(), x = e.clientX - b.left;
    if (Math.abs(x - al._drag.x) < 4) return;
    const r0 = al._map(Math.min(x, al._drag.x)), r1 = al._map(Math.max(x, al._drag.x)) - 1;
    if (r1 > r0) { st.rLo = r0; st.rHi = Math.min(st.N - 1, r1); al._moved = true; schedule(); }
  });
  al.addEventListener('pointerup', () => { al._drag = null; al._moved = false; });
  al.addEventListener('dblclick', () => { st.rLo = 0; st.rHi = st.N - 1; schedule(); });
  $('#clearN').addEventListener('click', () => { st.rLo = 0; st.rHi = st.N - 1; schedule(); });
  $('#clearG').addEventListener('click', () => { st.g = null; schedule(); });
  $('#colK').addEventListener('click', () => { st.colMode = 'k'; colButtons(); schedule(); });
  $('#colL').addEventListener('click', () => { st.colMode = 'L'; colButtons(); schedule(); });
  $('#ownForm').addEventListener('submit', (e) => {
    e.preventDefault();
    const v = $('#ownM').value.trim();
    if (!v) { if (st.ownM) { st.ownM = null; if (!LADDER.includes(st.M)) st.M = 1; schedule(); } return; }
    if (!/^\d+$/.test(v) || !setOwnM(+v)) { $('#ownM').setCustomValidity(`a whole number from 1 to ${fmtInt(OWN_MAX)}`); $('#ownM').reportValidity(); return; }
    $('#ownM').setCustomValidity('');
  });
  $('#ownM').addEventListener('input', () => $('#ownM').setCustomValidity(''));
  $('#sampleSel').addEventListener('change', (e) => { st.sample = +e.target.value; if (st.dec) startControls(); });
  $('#themeBtn').addEventListener('click', () => {
    const light = document.documentElement.dataset.theme ? document.documentElement.dataset.theme === 'light' : matchMedia('(prefers-color-scheme: light)').matches;
    document.documentElement.dataset.theme = light ? 'dark' : 'light'; schedule();
  });
  matchMedia('(prefers-color-scheme: light)').addEventListener('change', schedule);
  const ro = new ResizeObserver(schedule);
  for (const p of document.querySelectorAll('.plot')) ro.observe(p);
}
function colButtons() {
  $('#colK').setAttribute('aria-pressed', String(st.colMode === 'k'));
  $('#colL').setAttribute('aria-pressed', String(st.colMode === 'L'));
  $('#colsKey').style.setProperty('--c', st.colMode === 'k' ? 'var(--c-weight)' : 'var(--c-level)');
}

/* ═══════════════ your own sequence ═══════════════ */

function initPaste() {
  const dlg = $('#paste');
  $('#pasteBtn').addEventListener('click', () => { $('#pasteErr').textContent = ''; dlg.showModal(); });
  $('#pasteForm').addEventListener('submit', (e) => {
    if (e.submitter && e.submitter.value !== 'go') return;
    e.preventDefault();
    const res = parsePasted($('#pasteText').value);
    if (typeof res === 'string') { $('#pasteErr').textContent = res; return; }
    dlg.close();
    runPasted(res, $('#pasteName').value.trim() || 'Your sequence');
  });
}

/** b-file ("n a(n)" per line) or a plain list; → { n0, a } or an error message */
function parsePasted(text) {
  const lines = text.split(/\r?\n/).map(l => l.replace(/#.*/, '').trim()).filter(Boolean);
  const bfile = lines.length >= 3 && lines.every(l => /^-?\d+\s+-?\d+$/.test(l));
  let n0 = 1, nums;
  if (bfile) { const pr = lines.map(l => l.split(/\s+/)); n0 = +pr[0][0]; nums = pr.map(p => p[1]); }
  else nums = text.replace(/#.*$/gm, '').split(/[\s,;]+/).filter(Boolean);
  if (nums.some(s => !/^\d+$/.test(s))) return 'only non-negative integers, please';
  if (nums.length < 51) return `at least 51 terms (50 decomposed), not ${nums.length}`;
  if (nums.length > 100001) return 'at most 100,001 terms';
  const a = new Float64Array(nums.length);
  for (let i = 0; i < nums.length; i++) {
    const v = Number(nums[i]);
    if (!(v < 2 ** 53) || String(v) !== nums[i].replace(/^0+(?=\d)/, '')) return `term ${i + 1} is not below 2^53`;
    if (i && v <= a[i - 1]) return `not strictly increasing at term ${i + 1}`;
    a[i] = v;
  }
  return { n0, a };
}

function runPasted({ n0, a }, name) {
  const gen = ++st.gen; stopWorkers();
  history.replaceState(null, '', '#custom');
  const N = a.length - 1;
  st.rec = { custom: true, id: `custom-${DWLJ.seedOf(Array.from(a.subarray(0, 64)).join(','))}`, name, alias: name, terms: N, n0 };
  st.N = N; st.n0 = n0; st.a = a.slice(0, N); st.dec = null;
  header(); setStatus('decomposing…'); setProgress(0);
  const w = new Worker('worker.js'); st.workers.push(w);
  w.onmessage = (e) => {
    if (gen !== st.gen) return;
    const m = e.data;
    if (m.type === 'progress') setProgress(0.15 * m.done / N);
    else if (m.type === 'decompose') {
      w.terminate(); st.workers = [];
      st.d = m.d; st.k = m.k; st.L = m.L;
      let dec = 0; for (let i = 0; i < N; i++) if (st.k[i]) dec++;
      if (!dec) { setStatus('no term decomposes (a(n) > 2d never holds): nothing to show', 'bad'); return; }
      afterLoad({});
      st.verified = `${fmtInt(N)} terms decomposed here, ${fmtInt(dec)} decomposable`; st.verifiedOk = true;
      refreshStatus();
    } else if (m.type === 'error') setStatus(`could not decompose: ${m.message}`, 'bad');
  };
  w.postMessage({ type: 'decompose', task: 'paste', a });
}

/* ═══════════════ start ═══════════════ */

async function route() {
  const o = parseHash();
  if (o.id === 'custom' && st.rec && st.rec.custom) return;
  if (o.id === 'custom') { history.replaceState(null, '', '#primes'); o.id = 'primes'; }
  if (st.rec && !st.rec.custom && (o.id === st.rec.id || o.id.toUpperCase() === st.rec.anumber)) {
    /* same sequence, new filters (a pasted link) */
    if (o.own) setOwnM(o.own);
    if (o.M) setOwnM(o.M); st.g = o.g !== undefined ? o.g : null;
    if (o.n) { st.rLo = Math.max(0, o.n[0] - st.n0); st.rHi = Math.min(st.N - 1, o.n[1] - st.n0); } else { st.rLo = 0; st.rHi = st.N - 1; }
    schedule(); return;
  }
  await loadSequence(o);
}

async function start() {
  if (!window.Worker || !window.fetch) { setStatus('this browser lacks web workers or fetch', 'bad'); return; }
  if (location.protocol === 'file:') { setStatus('open this page through a web server: browsers block fetch() on file://', 'bad'); return; }
  initPicker(); initHover(); initPaste(); colButtons();
  try { await loadCatalog(); } catch (e) { setStatus(`could not load the catalogue: ${e.message}`, 'bad'); return; }
  addEventListener('hashchange', route);
  route();
}
start();

/* for the tests (tests/lab.mjs) */
window.__lab = st;
window.__labFns = { lockStats, obsCounts, ctrlCounts, reps, steps, setOwnM, ladder: () => ladderRows };
})();
