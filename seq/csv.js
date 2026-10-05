/* decompwlj 3D, written by tools/seo.py: the "Download CSV" button of a sequence page.
   Rebuilds n;a;weight;level;jump from the data chunks (the dwj2 format of tools/compact.py),
   with every check the viewer makes, then compares the file's SHA-256 with the one published
   in data/sha256.txt (shown on the page) and saves it only if they are equal. Nothing is stored. */
(function () {
  var btn = document.getElementById('csv');
  if (!btn) return;
  function inflate(bytes) {
    if (bytes.length > 2 && bytes[0] === 0x1f && bytes[1] === 0x8b) {
      var s = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
      return new Response(s).arrayBuffer().then(function (b) { return new Uint8Array(b); });
    }
    return Promise.resolve(bytes);
  }
  function rows(b, url, want, out) {
    if (!(b[0] === 0x64 && b[1] === 0x77 && b[2] === 0x6a && (b[3] === 0x31 || b[3] === 0x32))) throw new Error(url + ': not a dwj1/dwj2 chunk');
    var m = b[3] === 0x32 ? b[4] : 0, i = b[3] === 0x32 ? 5 : 4;
    if (m > 2) throw new Error(url + ': unknown jump coding ' + m);
    function unzz(z) { return z % 2 === 0 ? z / 2 : -(z + 1) / 2; }
    function rd() { var v = 0, m = 1, c; do { if (i >= b.length) throw new Error(url + ': truncated'); c = b[i++]; v += (c & 127) * m; m *= 128; } while (c & 128); return v; }
    var n = rd(), a = rd(), d = new Array(n);
    if (n !== want) throw new Error(url + ': ' + n + ' rows, expected ' + want);
    for (var r = 0, p = 0, pp = 0; r < n; r++) {
      var v = rd(); d[r] = m === 0 ? v : m === 1 ? p + unzz(v) : 2 * p - pp + unzz(v); pp = p; p = d[r];
      if (!(d[r] > 0) || d[r] >= 9007199254740992) throw new Error(url + ': jump out of range at row ' + r);
    }
    for (r = 0; r < n; r++) {
      var s = rd(), k = 0, L = 0;
      if (s) {
        var h = Math.floor(s / 2), q = (a - d[r]) / h;
        k = s % 2 ? q : h; L = s % 2 ? h : q;
        if (!Number.isInteger(q) || q * h !== a - d[r] || k <= d[r] || (s % 2 ? !(L < k) : !(k <= L)))
          throw new Error(url + ': row ' + r + ' does not fit a = ' + a);
      } else if (a > 2 * d[r]) throw new Error(url + ': row ' + r + ' should decompose');
      out.push([a, k, L, d[r]]);
      a += d[r];
    }
    if (i !== b.length) throw new Error(url + ': trailing bytes');
  }
  function hex(buf) { return Array.prototype.map.call(new Uint8Array(buf), function (x) { return ('0' + x.toString(16)).slice(-2); }).join(''); }
  btn.addEventListener('click', function () {
    var id = btn.dataset.id, an = btn.dataset.an, n0 = +btn.dataset.n0, nc = +btn.dataset.chunks;
    var terms = +btn.dataset.terms, per = +btn.dataset.rows, want = btn.dataset.sha256;
    var label = btn.textContent; btn.disabled = true; btn.textContent = 'Building…';
    var jobs = [];
    for (var c = 0; c < nc; c++) {
      (function (c) {
        var url = '../../data/seq/' + id + '/chunk-' + String(c).padStart(3, '0') + '.bin.gz';
        jobs.push(fetch(url).then(function (res) { if (!res.ok) throw new Error(res.status + ' ' + url); return res.arrayBuffer(); })
          .then(function (buf) { return inflate(new Uint8Array(buf)); })
          .then(function (b) { var out = []; rows(b, url, Math.min(per, terms - c * per), out); return out; }));
      })(c);
    }
    var text;
    Promise.all(jobs).then(function (parts) {
      for (var c = 1; c < parts.length; c++) {                /* each chunk starts where the last ended */
        var e = parts[c - 1][parts[c - 1].length - 1];
        if (parts[c][0][0] !== e[0] + e[3]) throw new Error('chunk ' + c + ' does not continue chunk ' + (c - 1));
      }
      text = ['n;a;weight;level;jump\n']; var buf = '', n = n0;
      parts.forEach(function (p) { p.forEach(function (r) {
        buf += n++ + ';' + r[0] + ';' + (r[1] || '') + ';' + (r[1] ? r[2] : '') + ';' + r[3] + '\n';
        if (buf.length > 65536) { text.push(buf); buf = ''; }
      }); });
      text.push(buf);
      if (n - n0 !== terms) throw new Error((n - n0) + ' rows, expected ' + terms);
      if (!want || !window.crypto || !crypto.subtle) return null;
      return new Blob(text).arrayBuffer().then(function (b) { return crypto.subtle.digest('SHA-256', b); });
    }).then(function (digest) {
      /* never save a file that differs from the published one */
      if (digest && hex(digest) !== want)
        throw new Error('the rebuilt file differs from the published one (SHA-256 ' + hex(digest).slice(0, 12) + '… instead of ' + want.slice(0, 12) + '…). Reload the page (Ctrl+F5); if it persists, please report it');
      var url = URL.createObjectURL(new Blob(text, { type: 'text/csv' }));
      var a = document.createElement('a'); a.href = url; a.download = 'decompwlj_' + an + '.csv';
      document.body.appendChild(a); a.click(); a.remove();
      try { window.goatcounter.count({ path: 'csv/' + an, title: 'CSV ' + an, event: true }); } catch (e) {}
      setTimeout(function () { URL.revokeObjectURL(url); }, 30000);
    }).catch(function (e) { alert('Could not build the CSV: ' + e.message); })
      .then(function () { btn.disabled = false; btn.textContent = label; });
  });
})();
