"""Rebuild the 4,757-chemical broad acute set at record level, gate it, then
   test the effort dependence of the minimum against the median at equal depth."""
import numpy as np, pandas as pd, json
from scipy.stats import spearmanr
RAW="/home/claude/ecotox_raw/ecotox_ascii_09_15_2026"
ACUTE={"LC50","EC50","IC50"}
UNIT={"mg/L":1.0,"ug/L":1e-3,"ng/L":1e-6,"g/L":1e3,"pg/L":1e-9,"ppm":1.0,"ppb":1e-3,
      "ppt":1e-6,"mg/dm3":1.0,"ug/dm3":1e-3,"AI mg/L":1.0,"AI ug/L":1e-3,
      "AE mg/L":1.0,"AE ug/L":1e-3,"mg/l":1.0,"ug/l":1e-3,"ng/l":1e-6}
t=pd.read_csv(f"{RAW}/tests.txt",sep="|",dtype=str,low_memory=False,on_bad_lines="skip",
              usecols=["test_id","test_cas","organism_habitat","media_type","species_number"])
t=t[t.organism_habitat.eq("Water") & t.media_type.fillna("").str.upper().str.startswith("FW")]
r=pd.read_csv(f"{RAW}/results.txt",sep="|",dtype=str,low_memory=False,on_bad_lines="skip",
              usecols=["test_id","endpoint","conc1_mean_op","conc1_mean","conc1_unit"])
m=r.merge(t,on="test_id")
m["ep"]=m.endpoint.astype(str).str.strip().str.upper().str.replace("*","",regex=False)
m=m[m.ep.isin(ACUTE)]
m["mult"]=m.conc1_unit.astype(str).str.strip().map(UNIT)
m=m[m["mult"].notna()]
m["val"]=pd.to_numeric(m.conc1_mean,errors="coerce")
m=m[m.val.notna()&(m.val>0)]
m=m[~m.conc1_mean_op.astype(str).str.strip().isin([">","<",">=","<=","~"])]
m["logTox"]=np.log10(m.val*m["mult"])
print("acute freshwater records:",len(m))
g=m.groupby("test_cas").logTox
d=pd.DataFrame({"n_rec":g.size(),"acute_median":g.median(),"acute_min":g.min(),
                "acute_spread":g.max()-g.min()}).reset_index()
print("chemicals:",len(d))

# ---------- GATE against the deposited file ----------
dep=pd.read_csv("/home/claude/pfasv32/canonical_results/ecotox_acute_smiles.csv")
dep["test_cas"]=dep.test_cas.astype(str); d["test_cas"]=d.test_cas.astype(str)
j=d.merge(dep,on="test_cas",suffixes=("_new","_dep"))
ok=(len(d)==len(dep)==len(j)
    and (j.n_rec_new==j.n_rec_dep).all()
    and np.allclose(j.acute_min_new,j.acute_min_dep)
    and np.allclose(j.acute_median_new,j.acute_median_dep))
print("GATE reproduce the deposited 4,757-chemical table:","PASS" if ok else "FAIL",
      f"({len(d)} vs {len(dep)}, matched {len(j)})")
assert ok
m.to_pickle("broad_records.pkl"); d.to_csv("broad_summary.csv",index=False)

# ---------- equal-depth subsampling ----------
pools={c:v.logTox.values for c,v in m.groupby("test_cas")}
nrec={c:len(v) for c,v in pools.items()}
rng=np.random.default_rng(20261121)
print("\n"+"="*104)
print("EQUAL-DEPTH SUBSAMPLING, 4,757 ECOTOX CHEMICALS OF ALL CLASSES")
print("="*104)
print(f"{'k':>3} {'eligible':>9} | {'MIN baseline':>12} {'MIN at depth k':>16} {'depth effect':>12}"
      f" | {'MEDIAN baseline':>15} {'MEDIAN at depth k':>18}")
rows=[]
for k in (3,5,10,20,50):
    el=[c for c,v in pools.items() if len(v)>=k]
    if len(el)<30: continue
    ln=np.log10([nrec[c] for c in el])
    bmin=spearmanr([pools[c].min() for c in el],ln)[0]
    bmed=spearmanr([np.median(pools[c]) for c in el],ln)[0]
    rmin=[];rmed=[]
    for _ in range(200):
        idx=[rng.choice(pools[c],k,replace=False) for c in el]
        rmin.append(spearmanr([v.min() for v in idx],ln)[0])
        rmed.append(spearmanr([np.median(v) for v in idx],ln)[0])
    smin,sdmin=float(np.mean(rmin)),float(np.std(rmin))
    smed=float(np.mean(rmed))
    print(f"{k:>3} {len(el):>9} | {bmin:>12.3f} {smin:>11.3f}±{sdmin:.3f} {bmin-smin:>12.3f}"
          f" | {bmed:>15.3f} {smed:>18.3f}")
    rows.append(dict(k=k,n_eligible=len(el),min_baseline=float(bmin),min_depth_equal=smin,
                     min_sd=sdmin,depth_effect=float(bmin-smin),
                     med_baseline=float(bmed),med_depth_equal=smed))
json.dump(rows,open("broad_subsample.json","w"),indent=1)
