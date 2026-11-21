"""How much of the required compensation does the EU-TGD assessment factor deliver?"""
import pandas as pd, numpy as np, json
from scipy import stats
L=pd.read_csv("ladder.csv")
A=L[L.C6_min_noec.notna()&L.AF_C6.notna()].copy()
A["x"]=np.log10(A.n_records); A["laf"]=np.log10(A.AF_C6)
A["crit"]=A.C6_min_noec-A.laf
sn=stats.linregress(A.x,A.C6_min_noec); sa=stats.linregress(A.x,A.laf); sc=stats.linregress(A.x,A.crit)
out=dict(n=len(A),
  slope_noec=sn.slope,p_noec=sn.pvalue,
  slope_af=sa.slope,p_af=sa.pvalue,
  slope_criterion=sc.slope,p_criterion=sc.pvalue,
  pct_compensation_delivered=100*abs(sa.slope)/abs(sn.slope),
  pct_passthrough=100*abs(sc.slope)/abs(sn.slope),
  af_range_log=float(np.log10(A.AF_C6.max()/A.AF_C6.min())),
  noec_range_log=float(A.C6_min_noec.max()-A.C6_min_noec.min()),
  af_tiers={str(int(k)):int(v) for k,v in A.AF_C6.value_counts().sort_index().items()})
for k,v in out.items(): print(f"  {k:32s} {v}")
json.dump(out,open("af_compensation.json","w"),indent=1)
