import pandas as pd, numpy as np
from collections import Counter
m=pd.read_pickle("m_ladder.pkl")
v=m.ngL.values
print("detections:",len(v))
# how many distinct values, and concentration of mass on the commonest
u=pd.Series(v).value_counts()
print("distinct values:",len(u))
print("top 15 values and their share:")
for val,c in u.head(15).items(): print(f"   {val:>10.4g}  n={c:5d}  {100*c/len(v):5.2f}%")
print(f"  top 10 values cover {100*u.head(10).sum()/len(v):.1f}% of detections")
# significant figures
def sig(x):
    s=f"{x:.10g}"
    s=s.replace("-","").replace(".","").lstrip("0").rstrip("0")
    return len(s) if s else 1
sg=pd.Series([sig(x) for x in v])
print("\nsignificant figures in the reported value:")
print((100*sg.value_counts(normalize=True).sort_index()).round(1).to_string())
# grid membership
def ongrid(x,step): 
    r=x/step
    return abs(r-round(r))<1e-9
for step,lbl in [(0.1,"0.1"),(0.5,"0.5"),(1.0,"1.0")]:
    print(f"  on a {lbl} ng/L grid: {100*np.mean([ongrid(x,step) for x in v]):.1f}%")
print("\nPFOS vs PFOA medians:")
for cas,n in [("1763-23-1","PFOS"),("335-67-1","PFOA")]:
    g=m[m.cas==cas]
    print(f"  {n}: n={len(g)} median={g.ngL.median():.4f}  "
          f"q25={g.ngL.quantile(.25):.3f} q75={g.ngL.quantile(.75):.3f}  "
          f"n at exactly the median={int((g.ngL==g.ngL.median()).sum())}")
# organisations
print("\ntop reporting organisations:")
print(m.OrganizationFormalName.value_counts().head(6).to_string())
print(f"  top organisation share: {100*m.OrganizationFormalName.value_counts().iloc[0]/len(m):.1f}%")
# does the grid change any conclusion? recompute C1 exceedance on a jittered set
rng=np.random.default_rng(7)
mm=m.copy(); mm["c1"]=(10**mm.C1_min_all)*1e6/1000
base=100*(mm.ngL>mm.c1).mean()
jit=[]
for _ in range(200):
    j=mm.ngL*rng.uniform(0.9,1.1,len(mm))
    jit.append(100*(j>mm.c1).mean())
print(f"\nC1 exceedance: observed {base:.2f}%  |  with +/-10% multiplicative jitter: "
      f"{np.mean(jit):.2f}% [{np.percentile(jit,2.5):.2f}, {np.percentile(jit,97.5):.2f}]")
