import pandas as pd, numpy as np
from scipy import stats
C=pd.read_csv("c7_ssd.csv"); L=pd.read_csv("wqp/ladder.csv")
m=C.merge(L[["parent","n_records","C1_min_all","C6_min_noec","AF_C6","C5_spmed_acute"]],on="parent")
rng=np.random.default_rng(20261121)
def rho(x,y,lbl):
    k=(~pd.isna(x))&(~pd.isna(y)); x,y=np.asarray(x)[k],np.asarray(y)[k]
    r,p=stats.spearmanr(x,y)
    bs=[stats.spearmanr(x[i],y[i])[0] for i in (rng.integers(0,len(y),len(y)) for _ in range(4000))]
    lo,hi=np.nanpercentile(bs,[2.5,97.5])
    print(f"  {lbl:46s} n={len(y):3d}  rho={r:+.3f}  95%CI [{lo:+.3f},{hi:+.3f}]  p={p:.2g}")
    return r
print("="*98)
print("DOES THE REGULATORY-STYLE SCREEN REMOVE THE TESTING-EFFORT DEPENDENCE?")
print("="*98)
x=np.log10(m.n_records)
rho(x,m.C1_min_all,          "C1  lowest, all endpoints (as published)")
rho(x,np.log10((10**m.C6_min_noec)/m.AF_C6), "C6  lowest chronic NOEC / EU-TGD AF")
rho(x,m.C5_spmed_acute,      "C5  median acute, species-resolved")
rho(x,m.C7_min,              "C7  lowest apical chronic SMCV")
rho(x,m.C7_gm,               "C7  geometric mean apical chronic SMCV")
s=m[m.C7_HC5.notna()]
rho(np.log10(s.n_records),s.C7_HC5, "C7  SSD HC5 (5th pct), >=5 species")
print("\n  n species surviving the C7 screen vs effort:")
rho(x,m.n_sp_C7,             "      (species count itself)")
print("\n" + "="*98)
print("WHY: does the screen remove MORE records from better-studied substances?")
R=pd.read_csv("results_parent.csv.gz")
DEP={"EC50","IC50","LC50","LOEC","LOEL","MATC","NOEC","NOEL"}
g=R[R.ep.isin(DEP)].groupby("parent").agg(n_all=("logTox","size"),
      n_apical=("effect",lambda s:s.isin({"MOR","GRO","REP","DVP","POP"}).sum()))
g["frac_apical"]=g.n_apical/g.n_all
g=g[g.n_all>=10]
rho(np.log10(g.n_all),g.frac_apical,"  apical fraction vs n records (>=10 recs)")
print(f"\n  apical fraction: substances with <50 records {g[g.n_all<50].frac_apical.median():.3f}"
      f" ; >=200 records {g[g.n_all>=200].frac_apical.median():.3f}")
