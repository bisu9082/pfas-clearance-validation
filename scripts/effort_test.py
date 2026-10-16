"""Decisive test: is the protective minimum driven by sampling depth?

If each substance is subsampled to the SAME number of records, an effort-driven
artefact disappears while a real chemistry signal survives.
"""
import sys, numpy as np, pandas as pd
from scipy.stats import spearmanr
RAW, REPO, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
SEED = 20260620

inv = pd.read_csv(f"{REPO}/untested_pfas_structure_locked.csv", usecols=["cas_norm"])
mea = pd.read_csv(f"{REPO}/measured_pfas_locked.csv", usecols=["cas"])
mea["cas_norm"] = mea.cas.astype(str).str.replace(r"\D","",regex=True)
cas = set(inv.cas_norm.astype(str)) | set(mea.cas_norm.astype(str))

tests = pd.read_csv(f"{RAW}/tests.txt", sep="|", dtype=str, low_memory=False,
                    on_bad_lines="skip",
                    usecols=["test_id","test_cas","organism_habitat","media_type"])
t = tests[tests.test_cas.astype(str).str.strip().isin(cas) &
          tests.organism_habitat.eq("Water") &
          tests.media_type.fillna("").str.upper().str.startswith("FW")]
res = pd.read_csv(f"{RAW}/results.txt", sep="|", dtype=str, low_memory=False,
                  on_bad_lines="skip",
                  usecols=["test_id","endpoint","conc1_mean_op","conc1_mean","conc1_unit"])
m = res.merge(t, on="test_id")
ACUTE={"LC50","EC50","IC50"}; CHRON={"NOEC","LOEC","MATC","NOEL","LOEL"}
m["ep"]=m.endpoint.astype(str).str.strip().str.upper().str.replace("*","",regex=False)
m=m[m.ep.isin(ACUTE|CHRON)]
UNIT={"mg/L":1.0,"ug/L":1e-3,"ng/L":1e-6,"g/L":1e3,"pg/L":1e-9,"ppm":1.0,"ppb":1e-3,
      "ppt":1e-6,"mg/dm3":1.0,"ug/dm3":1e-3,"AI mg/L":1.0,"AI ug/L":1e-3,
      "AE mg/L":1.0,"AE ug/L":1e-3,"mg/l":1.0,"ug/l":1e-3,"ng/l":1e-6}
m["mult"]=m.conc1_unit.astype(str).str.strip().map(UNIT)
m=m[m.mult.notna()]
m["val"]=pd.to_numeric(m.conc1_mean,errors="coerce")
m=m[m.val.notna()&(m.val>0)]
m=m[~m.conc1_mean_op.astype(str).str.strip().isin([">","<",">=","<=","~"])]
m["logTox"]=np.log10(m.val*m.mult)
m.to_csv(f"{OUT}/ecotox_pfas_records.csv", index=False)
print(f"usable records: {len(m)}  across {m.test_cas.nunique()} substances", flush=True)

g = m.groupby("test_cas").logTox
full = pd.DataFrame({"n": g.size(), "minimum": g.min(), "median": g.median()})
full["log_n"] = np.log10(full.n)
r_full = spearmanr(full.log_n, full.minimum)
r_med  = spearmanr(full.log_n, full["median"])
print(f"\nAll records:   minimum vs log(n): rho={r_full.statistic:+.3f} (p={r_full.pvalue:.1e})")
print(f"               median  vs log(n): rho={r_med.statistic:+.3f} (p={r_med.pvalue:.1e})")

print("\n=== subsample every substance to the SAME depth k ===")
rng = np.random.default_rng(SEED)
rows=[]
for k in (3, 5, 10, 20):
    elig = full[full.n >= k].index
    rr=[]
    for _ in range(200):
        mins = [rng.choice(m.loc[m.test_cas==c,"logTox"].to_numpy(), k, replace=False).min()
                for c in elig]
        rr.append(spearmanr(np.log10(full.loc[elig,"n"]), mins).statistic)
    rows.append((k, len(elig), float(np.mean(rr)), float(np.percentile(rr,2.5)), float(np.percentile(rr,97.5))))
    print(f"  k={k:2d}  substances={len(elig):3d}   rho(min, log n) = {np.mean(rr):+.3f} "
          f"[{np.percentile(rr,2.5):+.3f}, {np.percentile(rr,97.5):+.3f}]", flush=True)
pd.DataFrame(rows, columns=["k","n_sub","rho_mean","lo","hi"]).to_csv(f"{OUT}/effort_subsample.csv", index=False)
