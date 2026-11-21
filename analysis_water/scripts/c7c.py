import pandas as pd, numpy as np
from scipy import stats
R=pd.read_csv("results_parent.csv.gz")
DEP_EP={"EC50","IC50","LC50","LOEC","LOEL","MATC","NOEC","NOEL"}
G=R[R.ep.isin(DEP_EP)].copy()
dep=pd.read_csv("/home/claude/pfasv32/canonical_results/ecotox_pfas_measured_recurated.csv")
mn=G.drop_duplicates(subset=["test_id","ep","conc1_mean","conc1_unit"]).groupby("parent").logTox.min()
d=dep.set_index("parent").logTox
j=pd.DataFrame({"mine":mn,"dep":d}).dropna()
ok=int(((j.mine-j.dep).abs()<1e-9).sum())
print("GATE: deposited logTox min reproduced for %d / %d parents"%(ok,len(j)))
assert ok==len(j)==87, "gate failed"

APICAL={"MOR","GRO","REP","DVP","POP"}
C7EP={"NOEC","NOEL","MATC","EC10"}
sel=R[R.ep.isin(C7EP) & R.effect.isin(APICAL) & (R.conc1_type=="A")].copy()
print("\nC7 screen: chronic apical, active ingredient -> %d results, %d parents, %d species"
      %(len(sel),sel.parent.nunique(),sel.latin_name.nunique()))
print("  endpoint mix:",sel.ep.value_counts().to_dict())
print("  effect mix:  ",sel.effect.value_counts().to_dict())

# Species Mean Chronic Value = geometric mean within species
smcv=sel.groupby(["parent","latin_name"]).logTox.mean().reset_index()
rows=[]
for p,g in smcv.groupby("parent"):
    v=np.sort(g.logTox.values); n=len(v)
    mu=v.mean(); sd=v.std(ddof=1) if n>1 else np.nan
    hc5=mu+stats.norm.ppf(.05)*sd if (n>=5 and sd==sd) else np.nan
    rows.append(dict(parent=p,n_sp_C7=n,C7_min=v[0],C7_gm=mu,C7_sd=sd,C7_HC5=hc5))
C=pd.DataFrame(rows)
C.to_csv("c7_ssd.csv",index=False)
print("\nparents with >=5 species (SSD-able): %d ; >=3: %d"%(C.n_sp_C7.ge(5).sum(),C.n_sp_C7.ge(3).sum()))

NAME={"O=S(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFOS",
      "O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFOA"}
OFF={"PFOS":{"EPA 2024 CCC":250,"EU direct-tox AA-EQS":230,"AU DGV 95%":480,
             "AU DGV 99%":9.1,"CA FEQG":6800},
     "PFOA":{"EPA 2024 CCC":100000}}
print("\n"+"="*92)
print("C7 (regulatory-style screen) vs OFFICIAL FRESHWATER AQUATIC-LIFE CRITERIA, ng/L")
print("="*92)
for smi,nm in NAME.items():
    r=C[C.parent==smi]
    if not len(r): print(nm,"absent"); continue
    r=r.iloc[0]
    gm=10**np.mean(np.log10(list(OFF[nm].values())))
    print(f"\n{nm}: {int(r.n_sp_C7)} species after screen")
    for lbl,col,af in [("C7 HC5 (SSD 5th pct), no AF","C7_HC5",1.0),
                       ("C7 lowest SMCV / AF 10","C7_min",10.0),
                       ("C7 geometric mean SMCV","C7_gm",1.0)]:
        if pd.isna(r[col]): print(f"  {lbl:32s} n/a"); continue
        v=(10**r[col])*1e6/af      # mg/L -> ng/L
        print(f"  {lbl:32s} {v:>12,.4g} ng/L   {v/gm:>10,.3g}x official gm ({np.log10(v/gm):+.2f} log)")
    print(f"  official: "+", ".join(f"{k} {v:,.4g}" for k,v in OFF[nm].items())+f"  | gm {gm:,.4g}")
