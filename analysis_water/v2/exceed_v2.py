import pandas as pd, numpy as np
from scipy import stats
rng=np.random.default_rng(20261122)
out=pd.read_csv("ladder_v2.csv")
lad=pd.read_csv("/home/claude/wq/wqp/ladder.csv")
w=pd.read_pickle("/home/claude/wq/wqp/m_ladder.pkl")
w["ugL"]=w.ngL/1000.0
L=lad.copy()

def pnec_map(col,AF,src=None):
    s=(out if src is None else src).dropna(subset=[col]).set_index("parent")
    a=s[AF] if isinstance(AF,str) else AF
    if isinstance(AF,str): s=s.dropna(subset=[AF]); a=s[AF]
    return (10**s[col])*1e6/a      # ng/L

CONV=[("C1 lowest, all endpoints",     pnec_map("C1_min_all",1000.0,L)),
      ("C3 lowest acute",              pnec_map("C3_min_acute",1000.0,L)),
      ("C5 median acute, species",     pnec_map("C5_spmed_acute",1000.0,L)),
      ("C6 EU-TGD (as published)",     pnec_map("C6_pub","AF")),
      ("C6' EU-TGD (apical+duration)", pnec_map("C6_prime","AF")),
      ("C7a lowest apical SMCV /10",   pnec_map("minSMCV_d",10.0)),
      ("C7b SSD HC5 (duration-filt.)", pnec_map("HC5_d",1.0))]

def tab(restrict,title):
    print("\n"+"="*94); print(title); print("="*94)
    print(f"{'convention':32s} {'n_meas':>7s} {'n_sub':>6s} {'medPNEC ng/L':>13s} {'%samples':>9s} {'%stations':>10s}")
    res={}
    for nm,pm in CONV:
        x=w[w.parent.isin(pm.index)].copy()
        if restrict is not None: x=x[x.parent.isin(restrict)]
        if not len(x): continue
        x["RQ"]=x.ugL/(x.parent.map(pm)/1000.0)
        st=x.groupby("MonitoringLocationIdentifier").RQ.max()
        med=pm[pm.index.isin(x.parent.unique())].median()
        print(f"{nm:32s} {len(x):7d} {x.parent.nunique():6d} {med:13,.3g} "
              f"{100*(x.RQ>1).mean():9.2f} {100*(st>1).mean():10.2f}")
        res[nm]=x
    return res

A=tab(None,"EACH CONVENTION ON ITS OWN SUBSTANCE SET")
common=set.intersection(*[set(pm.index) for _,pm in CONV])
common={c for c in common if c in set(w.parent)}
B=tab(common,f"COMMON SUBSET: {len(common)} substances, identical measurements")

# matched pair under C6'
PFOA="O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
PFHpA="O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
print("\n"+"="*94); print("MATCHED PAIR under C6' (apical + duration)"); print("="*94)
pm=dict(CONV)["C6' EU-TGD (apical+duration)"]
st_a=set(w[w.parent==PFOA].MonitoringLocationIdentifier); st_b=set(w[w.parent==PFHpA].MonitoringLocationIdentifier)
sh=st_a&st_b
print(f"  shared stations: {len(sh)}")
for p,nm in [(PFOA,"PFOA"),(PFHpA,"PFHpA")]:
    x=w[(w.parent==p)&w.MonitoringLocationIdentifier.isin(sh)]
    cr=pm.get(p,np.nan)
    nrec=int(out[out.parent==p].n_records.iloc[0])
    nap=int(( pd.read_csv("/home/claude/wq/results_parent.csv.gz",low_memory=False)
              .pipe(lambda z: z[(z.parent==p)&(z.conc1_type=="A")&z.apical.astype(bool)
                   &z.ep.isin(["NOEC","NOEL"])]) ).shape[0])
    print(f"  {nm:6s} criterion {cr:10,.1f} ng/L   median {x.ngL.median():7.2f} ng/L   "
          f"exceeding {100*(x.ngL>cr).mean():6.2f}%   records {nrec:5d} (apical chronic {nap})")
print(f"  criterion ratio: {pm[PFHpA]/pm[PFOA]:,.0f}x")
