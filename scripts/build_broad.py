"""Build the acute-median training candidate set from ECOTOX (all chemicals)."""
import sys, numpy as np, pandas as pd
RAW, OUT = sys.argv[1], sys.argv[2]
ACUTE = {"LC50","EC50","IC50"}
UNIT={"mg/L":1.0,"ug/L":1e-3,"ng/L":1e-6,"g/L":1e3,"pg/L":1e-9,"ppm":1.0,"ppb":1e-3,
      "ppt":1e-6,"mg/dm3":1.0,"ug/dm3":1e-3,"AI mg/L":1.0,"AI ug/L":1e-3,
      "AE mg/L":1.0,"AE ug/L":1e-3,"mg/l":1.0,"ug/l":1e-3,"ng/l":1e-6}

t = pd.read_csv(f"{RAW}/tests.txt", sep="|", dtype=str, low_memory=False, on_bad_lines="skip",
                usecols=["test_id","test_cas","organism_habitat","media_type","species_number"])
t = t[t.organism_habitat.eq("Water") & t.media_type.fillna("").str.upper().str.startswith("FW")]
r = pd.read_csv(f"{RAW}/results.txt", sep="|", dtype=str, low_memory=False, on_bad_lines="skip",
                usecols=["test_id","endpoint","conc1_mean_op","conc1_mean","conc1_unit"])
m = r.merge(t, on="test_id")
m["ep"] = m.endpoint.astype(str).str.strip().str.upper().str.replace("*","",regex=False)
m = m[m.ep.isin(ACUTE)]
m["mult"] = m.conc1_unit.astype(str).str.strip().map(UNIT)
m = m[m.mult.notna()]
m["val"] = pd.to_numeric(m.conc1_mean, errors="coerce")
m = m[m.val.notna() & (m.val > 0)]
m = m[~m.conc1_mean_op.astype(str).str.strip().isin([">","<",">=","<=","~"])]
m["logTox"] = np.log10(m.val * m.mult)

g = m.groupby("test_cas").logTox
d = pd.DataFrame({"n_rec": g.size(), "acute_median": g.median(),
                  "acute_min": g.min(), "acute_spread": g.max()-g.min()}).reset_index()
d.to_csv(f"{OUT}/ecotox_acute_all.csv", index=False)
print(f"acute freshwater records: {len(m)}")
print(f"unique chemicals with acute median: {len(d)}")
for k in (1,2,3,5,10):
    print(f"   with >= {k:2d} records: {(d.n_rec>=k).sum()}")
from scipy.stats import spearmanr
print(f"\neffort confound on acute median: rho = {spearmanr(np.log10(d.n_rec), d.acute_median).statistic:+.3f}")
print(f"effort confound on acute minimum: rho = {spearmanr(np.log10(d.n_rec), d.acute_min).statistic:+.3f}")
