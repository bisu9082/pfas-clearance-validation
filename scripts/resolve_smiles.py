"""Resolve ECOTOX CAS numbers to SMILES via PubChem PUG REST, with resume."""
import sys, os, time, json
import pandas as pd, requests

SRC, OUT = sys.argv[1], sys.argv[2]
CACHE = f"{OUT}/cas_smiles_cache.json"
RATE = 0.21          # ~5 requests/second, PubChem's published limit

def dashed(cas):
    s = "".join(ch for ch in str(cas) if ch.isdigit())
    if len(s) < 5: return None
    return f"{s[:-3]}-{s[-3:-1]}-{s[-1]}"

cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
d = pd.read_csv(SRC, dtype={"test_cas": str})
todo = [c for c in d.test_cas.astype(str) if c not in cache]
print(f"{len(d)} chemicals | cached {len(cache)} | to resolve {len(todo)}", flush=True)

sess = requests.Session()
sess.headers.update({"User-Agent": "research-recuration/1.0"})
BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/SMILES,CanonicalSMILES/JSON"

for i, c in enumerate(todo):
    dc = dashed(c)
    smi = None
    if dc:
        try:
            r = sess.get(BASE.format(dc), timeout=20)
            if r.status_code == 200:
                props = r.json()["PropertyTable"]["Properties"][0]
                smi = props.get("SMILES") or props.get("CanonicalSMILES")
        except Exception:
            smi = None
    cache[c] = smi
    time.sleep(RATE)
    if (i + 1) % 250 == 0:
        json.dump(cache, open(CACHE, "w"))
        hit = sum(1 for v in cache.values() if v)
        print(f"  {i+1}/{len(todo)}  resolved so far: {hit}/{len(cache)}", flush=True)

json.dump(cache, open(CACHE, "w"))
d["smiles"] = d.test_cas.astype(str).map(cache)
d.to_csv(f"{OUT}/ecotox_acute_smiles.csv", index=False)
ok = d.smiles.notna().sum()
print(f"\nDONE  resolved {ok}/{len(d)} ({ok/len(d):.1%})", flush=True)
