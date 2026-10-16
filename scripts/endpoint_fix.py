"""Does a depth-controlled endpoint improve external transfer?

The broad non-PFAS model is held FIXED (as published). Only the external PFAS
target changes. Any improvement is attributable to endpoint construction.
"""
import sys, numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')
BROAD, RECS, PFASMAP, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
SEED = 20260620

def feats(sm):
    X,k=[],[]
    for i,s in enumerate(sm):
        mm=Chem.MolFromSmiles(str(s))
        if mm is not None:
            X.append(np.array(AllChem.GetMorganFingerprintAsBitVect(mm,2,2048),dtype=np.int8)); k.append(i)
    return np.vstack(X), np.array(k)
def met(y,p):
    ss=float(np.sum((y-np.mean(y))**2))
    return (1-float(np.sum((y-p)**2))/ss if ss>0 else np.nan, float(spearmanr(y,p).statistic))
def ci(y,p,B=3000,rng=None):
    rng=rng or np.random.default_rng(SEED); n=len(y); o=[]
    for _ in range(B):
        i=rng.integers(0,n,n)
        if len(np.unique(y[i]))>2: o.append(spearmanr(y[i],p[i]).statistic)
    return float(np.nanpercentile(o,2.5)), float(np.nanpercentile(o,97.5))

broad=pd.read_csv(BROAD); Xb,kb=feats(broad.canonical_smiles); yb=broad.logTox.to_numpy()[kb]
full=RandomForestRegressor(n_estimators=300,max_features="sqrt",random_state=SEED,n_jobs=-1).fit(Xb,yb)
print("broad model fitted (unchanged from the manuscript)\n", flush=True)

rec=pd.read_csv(RECS, dtype={"test_cas":str})
smap=pd.read_csv(PFASMAP, dtype={"cas_norm":str})[["cas_norm","canonical_smiles","name"]].drop_duplicates("cas_norm")
ACUTE={"LC50","EC50","IC50"}
rng=np.random.default_rng(SEED)

def build(kind, k=5, nrep=50):
    g=rec.groupby("test_cas")
    rows=[]
    for c,sub in g:
        v=sub.logTox.to_numpy()
        if kind=="minimum":        val=v.min()
        elif kind=="median":       val=float(np.median(v))
        elif kind=="acute_median":
            a=sub[sub.ep.isin(ACUTE)].logTox.to_numpy()
            if len(a)==0: continue
            val=float(np.median(a))
        elif kind=="depth_min":
            if len(v)<k: continue
            val=float(np.mean([rng.choice(v,k,replace=False).min() for _ in range(nrep)]))
        rows.append((c,val,len(v)))
    d=pd.DataFrame(rows,columns=["cas_norm","y","n_rec"]).merge(smap,on="cas_norm",how="inner")
    d=d[d.canonical_smiles.notna()]
    return d

print(f"{'endpoint':16s} {'n':>4s} {'R2':>8s} {'rho':>8s}   {'95% CI':>18s}   rho(y, log n)")
res=[]
for kind in ("minimum","median","acute_median","depth_min"):
    d=build(kind)
    X,kk=feats(d.canonical_smiles); d=d.iloc[kk].reset_index(drop=True)
    y=d.y.to_numpy(); p=full.predict(X)
    r2,rho=met(y,p); lo,hi=ci(y,p)
    reff=spearmanr(np.log10(d.n_rec), y).statistic
    print(f"{kind:16s} {len(y):4d} {r2:+8.3f} {rho:+8.3f}   [{lo:+.2f}, {hi:+.2f}]   {reff:+.3f}")
    res.append(dict(endpoint=kind,n=len(y),r2=r2,rho=rho,lo=lo,hi=hi,rho_effort=reff))
pd.DataFrame(res).to_csv(f"{OUT}/endpoint_fix_results.csv", index=False)
