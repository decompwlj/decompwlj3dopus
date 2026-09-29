/* csv.js — CSV readers.
 *
 * Two readers, because the two files have opposite shapes:
 *
 *   parseIntRows()  bulk numeric data.  Reads bytes straight out of an
 *                   ArrayBuffer into pre-allocated Float64Array columns —
 *                   no strings, no per-row allocation, no split().  ~1e8
 *                   rows/minute in a warm JIT, and it scales linearly.
 *
 *   parseRecords()  the catalogue.  Full RFC-4180 handling (quotes, escaped
 *                   quotes, embedded commas and newlines) because the note
 *                   column is prose.  Fed incrementally so a large catalogue
 *                   can be shown while it is still downloading.
 */

/** Read `cols.length` unsigned integers per line from `bytes` into the
 *  Float64Array columns, starting at row `row0`.  Any line whose first
 *  byte is not a digit (a header, a comment, a blank) is skipped.
 *  Values are exact to 2^53, i.e. well past any 64-bit sequence term
 *  the decomposition produces in practice.
 *  @returns {number} rows written */
export function parseIntRows(bytes, cols, row0) {
  const n = bytes.length, ncol = cols.length;
  let i = 0, row = row0;
  while (i < n) {
    let c = bytes[i];
    if (c < 48 || c > 57) {                    // not a digit: skip the line
      while (i < n && bytes[i] !== 10) i++;
      i++;
      continue;
    }
    for (let col = 0; ; col++) {
      let v = 0;
      while (i < n) {
        const ch = bytes[i];
        if (ch < 48 || ch > 57) break;
        v = v * 10 + (ch - 48);
        i++;
      }
      if (col < ncol) cols[col][row] = v;
      if (i < n && bytes[i] === 44) { i++; continue; }   // ','
      while (i < n && bytes[i] !== 10) i++;              // to end of line
      i++;
      break;
    }
    row++;
  }
  return row - row0;
}

/** Incremental RFC-4180 parser.  Feed it text with push(); it calls
 *  onRow(array of strings) for every complete row.  Call end() at EOF. */
export class RecordParser {
  constructor(onRow) {
    this.onRow = onRow;
    this.buf = '';
    this.field = '';
    this.row = [];
    this.inQuote = false;
    this.wasQuote = false;
  }
  push(text) {
    const s = this.buf + text;
    this.buf = '';
    let i = 0;
    const n = s.length;
    while (i < n) {
      const ch = s[i];
      if (this.inQuote) {
        if (ch === '"') {
          if (i + 1 < n) {
            if (s[i + 1] === '"') { this.field += '"'; i += 2; continue; }
            this.inQuote = false; i++; continue;
          }
          this.buf = s.slice(i);                 // '"' at the very end: wait
          return;
        }
        this.field += ch; i++; continue;
      }
      if (ch === '"' && this.field === '') { this.inQuote = true; i++; continue; }
      if (ch === ',') { this.row.push(this.field); this.field = ''; i++; continue; }
      if (ch === '\n' || ch === '\r') {
        if (ch === '\r' && i + 1 === n) { this.buf = s.slice(i); return; }
        if (ch === '\r' && s[i + 1] === '\n') i++;
        this.row.push(this.field); this.field = '';
        if (this.row.length > 1 || this.row[0] !== '') this.onRow(this.row);
        this.row = []; i++; continue;
      }
      this.field += ch; i++;
    }
  }
  end() {
    if (this.field !== '' || this.row.length) {
      this.row.push(this.field);
      if (this.row.length > 1 || this.row[0] !== '') this.onRow(this.row);
    }
    this.field = ''; this.row = [];
  }
}

/** Stream a CSV over the network, emitting rows as they arrive.
 *  Falls back to a single text() read where streaming is unavailable. */
export async function streamRecords(url, onRow, signal) {
  const res = await fetch(url, { signal });
  if (!res.ok) throw new Error(`${res.status} ${res.statusText} — ${url}`);
  const parser = new RecordParser(onRow);
  if (!res.body || !res.body.getReader) {
    parser.push(await res.text());
    parser.end();
    return;
  }
  const reader = res.body.getReader();
  const dec = new TextDecoder('utf-8');
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    parser.push(dec.decode(value, { stream: true }));
  }
  parser.push(dec.decode());
  parser.end();
}
