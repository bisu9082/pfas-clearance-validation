import pandas as pd, numpy as np, json
from scipy import stats
rng=np.random.default_rng(20261122)
out=pd.read_csv("ladder_v2.csv")
d=pd.read_csv("/home/claude/wq/results_parent.csv.gz",low_memory=False)

def boot_rho(x,y,B=5000,seed=1):
    r=np.random.default_rng(seed); n=len(x); acc=[]
    x=np.asarray(x); y=np.asarray(y)
    for _ in range(B):
        j=r.integers(0,n,n)
        if len(set(x[j]))<3: continue
        acc.append(stats.spearmanr(x[j],y[j])[0])
    return np.nanpercentile(acc,[2.5,97.5])

print("="*86); print("ROBUSTNESS — every check the panel demanded"); print("="*86)
CS=[("C6_prime","C6' EU-TGD apical+duration"),("HC5_d","C7b SSD HC5 duration-filt."),
    ("minSMCV_d","C7a lowest apical SMCV")]

# 1. chain length control
from rdkit import Chem
def nC(s):
    m=Chem.MolFromSmiles(s); return sum(1 for a in m.GetAtoms() if a.GetSymbol()=="C") if m else np.nan
out["nC"]=out.parent.map(nC)
def partial(x,y,z):
    rx,ry,rz=(stats.rankdata(v) for v in (x,y,z))
    res=lambda a,b: a-np.polyval(np.polyfit(b,a,1),b)
    return np.corrcoef(res(rx,rz),res(ry,rz))[0,1]
print("\n[1] chain-length confound")
for c,nm in CS:
    s=out.dropna(subset=[c,"nC"]); ln=np.log10(s.n_records)
    print(f"  {nm:30s} rho={stats.spearmanr(ln,s[c])[0]:+.3f}  "
          f"partial|nC={partial(ln,s[c],s.nC):+.3f}  (rho(logn,nC)={stats.spearmanr(ln,s.nC)[0]:+.2f})")

print("\n[2] effort variable: result rows vs distinct tests")
for c,nm in CS:
    s=out.dropna(subset=[c])
    print(f"  {nm:30s} records={stats.spearmanr(np.log10(s.n_records),s[c])[0]:+.3f}  "
          f"tests={stats.spearmanr(np.log10(s.n_tests),s[c])[0]:+.3f}")

print("\n[3] drop the 3 non-OECD-PFAS parents")
bad=["O=Cc1c(F)c(F)c(F)c(F)c1F","Nc1c(F)c(F)c(F)c(F)c1F","O=C(O)C(F)F"]
for c,nm in CS:
    s=out.dropna(subset=[c]); s2=s[~s.parent.isin(bad)]
    print(f"  {nm:30s} all={stats.spearmanr(np.log10(s.n_records),s[c])[0]:+.3f}  "
          f"PFAS-only={stats.spearmanr(np.log10(s2.n_records),s2[c])[0]:+.3f}  (n {len(s)}->{len(s2)})")

print("\n[4] bootstrap CIs, fixed seed, 5000 reps")
for i,(c,nm) in enumerate(CS):
    s=out.dropna(subset=[c]); lo,hi=boot_rho(np.log10(s.n_records),s[c],seed=100+i)
    print(f"  {nm:30s} rho={stats.spearmanr(np.log10(s.n_records),s[c])[0]:+.3f}  CI=[{lo:+.2f}, {hi:+.2f}]  n={len(s)}")

print("\n[5] null simulation: is rho=-0.68 at n=13 reachable by chance?")
s=out.dropna(subset=["HC5_d"]); n=len(s)
r0=np.random.default_rng(9); nullr=[]
for _ in range(20000):
    nullr.append(stats.spearmanr(np.log10(s.n_records),r0.permutation(s.HC5_d.values))[0])
nullr=np.array(nullr); obs=stats.spearmanr(np.log10(s.n_records),s.HC5_d)[0]
print(f"  observed {obs:+.3f};  null mean {nullr.mean():+.3f} sd {nullr.std():.3f};  "
      f"95% band [{np.percentile(nullr,2.5):+.2f},{np.percentile(nullr,97.5):+.2f}];  "
      f"P(rho<=obs)={np.mean(nullr<=obs):.4f}")

print("\n[6] leave-one-organisation-out on the C1 exceedance")
w=pd.read_pickle("/home/claude/wq/wqp/m_ladder.pkl"); w["ugL"]=w.ngL/1000.0
L=pd.read_csv("/home/claude/wq/wqp/ladder.csv").dropna(subset=["C1_min_all"]).set_index("parent")
pm=(10**L.C1_min_all)*1e6
x=w[w.parent.isin(pm.index)].copy(); x["RQ"]=x.ngL/x.parent.map(pm)
full=100*(x.RQ>1).mean()
oc=[c for c in x.columns if "Organization" in c and "Formal" in c] or [c for c in x.columns if "Organization" in c]
if oc:
    o=oc[0]; top=x[o].value_counts().head(8)
    vals=[100*(x[x[o]!=k].RQ>1).mean() for k in top.index]
    print(f"  full {full:.2f}% ;  leave-one-out range [{min(vals):.2f}, {max(vals):.2f}]  (top org {100*top.iloc[0]/len(x):.1f}%)")
print("\n[7] one value per station x substance")
m1=x.groupby(["MonitoringLocationIdentifier","parent"]).RQ.median()
print(f"  station x substance median weighting: {100*(m1>1).mean():.2f}%  (vs {full:.2f}%)")
