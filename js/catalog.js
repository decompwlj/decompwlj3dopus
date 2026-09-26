/* catalog.js — the sequence catalogue and its virtualised list.
 *
 * The catalogue is one CSV row per sequence.  It is streamed and parsed
 * incrementally, and the list renders only the rows inside the viewport,
 * so the cost of showing it is independent of how many sequences exist:
 * ten today, tens of thousands later, same DOM, same frame budget.
 */

import { streamRecords } from './csv.js';

const NUM = new Set(['n0', 'terms', 'chunks', 'chunk_rows', 'decomposable',
                     'level', 'weight', 'ties', 'level_one',
                     'amin', 'amax', 'kmax', 'Lmax', 'dmin', 'dmax']);

/** Stream data/catalog.csv.  onBatch(records) fires as rows arrive. */
export async function loadCatalog(url, onBatch, signal) {
  let header = null;
  const all = [];
  let batch = [];
  await streamRecords(url, (row) => {
    if (!header) { header = row.map(s => s.trim()); return; }
    const r = {};
    for (let i = 0; i < header.length; i++) {
      const key = header[i];
      const v = row[i] === undefined ? '' : row[i];
      r[key] = NUM.has(key) ? Number(v) : v;
    }
    if (!r.id) return;
    r._key = `${r.id} ${r.anumber} ${r.name} ${r.family}`.toLowerCase();
    r.levelShare = r.decomposable ? r.level / r.decomposable : 0;
    all.push(r);
    batch.push(r);
    if (batch.length >= 512) { onBatch(batch, all); batch = []; }
  }, signal);
  if (batch.length) onBatch(batch, all);
  return all;
}

/** Case-insensitive token search: every token must appear in the key. */
export function filterCatalog(all, query) {
  const q = query.trim().toLowerCase();
  if (!q) return all;
  const toks = q.split(/\s+/);
  const out = [];
  for (let i = 0; i < all.length; i++) {
    const k = all[i]._key;
    let ok = true;
    for (let t = 0; t < toks.length; t++) if (k.indexOf(toks[t]) === -1) { ok = false; break; }
    if (ok) out.push(all[i]);
  }
  return out;
}

const fmt = new Intl.NumberFormat('en-US');

export class VirtualList {
  /** @param {HTMLElement} scroller  the scrolling container
   *  @param {HTMLElement} spacer    child sized to the full list height
   *  @param {HTMLElement} rows      absolutely-positioned row host */
  constructor(scroller, spacer, rows, onPick) {
    this.scroller = scroller;
    this.spacer = spacer;
    this.rows = rows;
    this.onPick = onPick;
    this.items = [];
    this.selected = null;
    this.rowH = 58;
    this.pool = [];
    this.overscan = 6;

    scroller.addEventListener('scroll', () => this.render(), { passive: true });
    rows.addEventListener('click', (e) => {
      const el = e.target.closest('.row');
      if (el && el.dataset.id) this.onPick(el.dataset.id);
    });
    new ResizeObserver(() => this.render()).observe(scroller);
  }

  setItems(items) {
    this.items = items;
    this.spacer.style.height = (items.length * this.rowH) + 'px';
    this.render();
  }

  select(id) {
    this.selected = id;
    for (const el of this.pool) el.classList.toggle('sel', el.dataset.id === id);
  }

  scrollTo(id) {
    const i = this.items.findIndex(r => r.id === id);
    if (i < 0) return;
    const top = i * this.rowH, h = this.scroller.clientHeight;
    if (top < this.scroller.scrollTop || top + this.rowH > this.scroller.scrollTop + h)
      this.scroller.scrollTop = Math.max(0, top - h / 2 + this.rowH / 2);
  }

  render() {
    const { scroller, items, rowH } = this;
    const top = scroller.scrollTop;
    const h = scroller.clientHeight || 600;
    let first = Math.max(0, Math.floor(top / rowH) - this.overscan);
    let last = Math.min(items.length, Math.ceil((top + h) / rowH) + this.overscan);
    const need = last - first;

    while (this.pool.length < need) {
      const el = document.createElement('div');
      el.className = 'row';
      el.setAttribute('role', 'option');
      el.innerHTML = '<div class="r1"><span class="nm"></span><span class="an"></span></div>' +
                     '<div class="r2"><span class="ct"></span><span class="bar"><i></i></span><span class="pct"></span></div>';
      el._nm = el.querySelector('.nm');
      el._an = el.querySelector('.an');
      el._ct = el.querySelector('.ct');
      el._bar = el.querySelector('.bar > i');
      el._pct = el.querySelector('.pct');
      this.rows.appendChild(el);
      this.pool.push(el);
    }
    for (let i = 0; i < this.pool.length; i++) {
      const el = this.pool[i];
      if (i >= need) { el.style.display = 'none'; continue; }
      const idx = first + i, r = items[idx];
      if (!r) { el.style.display = 'none'; continue; }
      el.style.display = '';
      el.style.top = (idx * rowH) + 'px';
      el.dataset.id = r.id;
      el.classList.toggle('sel', r.id === this.selected);
      el._nm.textContent = r.name;
      el._an.textContent = r.anumber;
      el._ct.textContent = fmt.format(r.terms) + ' terms';
      el._bar.style.width = (r.levelShare * 100).toFixed(2) + '%';
      el._pct.textContent = (r.levelShare * 100).toFixed(1) + '% lvl';
      el.title = `${r.name} (${r.anumber}) — ${r.family}\n${fmt.format(r.decomposable)} decomposable, ` +
                 `${fmt.format(r.level)} level-classified`;
    }
  }
}
