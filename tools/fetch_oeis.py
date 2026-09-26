#!/usr/bin/env python3
"""fetch_oeis.py — the OEIS records the site needs: name, offset and first terms.

    python3 fetch_oeis.py A000040 A000027 ...     add or refresh these entries
    python3 fetch_oeis.py --catalog raw           every A-number in raw/catalog.csv

Reads the OEIS's own data repository (github.com/oeis/oeisdata, one .seq file per
sequence) and keeps what it needs in oeis.json beside this script:
    { "A000040": { "name": "The prime numbers.", "offset": 1, "terms": [2, 3, ...] }, ... }
audit.py checks every sequence's first rows against these terms, and the
catalogue takes its names from here.
"""
import concurrent.futures, csv, json, pathlib, sys, urllib.request

HERE = pathlib.Path(__file__).resolve().parent
STORE = HERE / 'oeis.json'
URL = 'https://raw.githubusercontent.com/oeis/oeisdata/main/seq/{d}/{a}.seq'

def fetch(anum):
    url = URL.format(d=anum[:4], a=anum)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                text = r.read().decode('utf-8')
            break
        except Exception as e:
            if attempt == 3: raise RuntimeError(f'{anum}: {e}')
    rec = {'name': '', 'offset': None, 'terms': ''}
    for line in text.splitlines():
        tag, rest = line[:2], line[11:] if len(line) > 11 else ''
        if tag == '%N': rec['name'] = rest.strip()
        elif tag == '%O': rec['offset'] = int(rest.split(',')[0])
        elif tag in ('%S', '%T', '%U'): rec['terms'] += rest.strip()
    rec['terms'] = [int(x) for x in rec['terms'].split(',') if x]
    if not rec['name'] or rec['offset'] is None or not rec['terms']:
        raise RuntimeError(f'{anum}: incomplete record')
    return anum, rec

def main():
    args = sys.argv[1:]
    if args and args[0] == '--catalog':
        args = [r['anumber'] for r in csv.DictReader(open(pathlib.Path(args[1]) / 'catalog.csv'))]
    store = json.loads(STORE.read_text()) if STORE.exists() else {}
    with concurrent.futures.ThreadPoolExecutor(8) as ex:
        for anum, rec in ex.map(fetch, args):
            store[anum] = rec
    STORE.write_text(json.dumps(dict(sorted(store.items())), indent=0, ensure_ascii=False) + '\n')
    print(f'{len(args)} fetched, {len(store)} in {STORE.name}')

if __name__ == '__main__':
    main()
