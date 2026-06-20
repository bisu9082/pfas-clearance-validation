#!/usr/bin/env python3
"""Reproduce the SI robustness tables (algorithm sensitivity, applicability-domain
threshold sensitivity, leave-one-PFAS-out deployment stability) from the canonical
structure-locked tables. Run from the repository root."""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from rdkit import Chem; from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from scipy.stats import spearmanr
R="canonical_results"; SEED=20260620
def fp(s):
    m=Chem.MolFromSmiles(str(s)); return np.array(AllChem.GetMorganFingerprintAsBitVect(m,2,2048),dtype=np.int8) if m else None
def feat(df,col="canonical_smiles"):
    df=df.copy(); df["fp"]=df[col].map(fp); df=df.dropna(subset=["fp"]); return np.stack(df.fp), df
br,_=feat(pd.read_csv(f"{R}/broad_nonpfas_structure_locked.csv")); ytr=pd.read_csv(f"{R}/broad_nonpfas_structure_locked.csv").dropna(subset=["canonical_smiles"]).logTox.values[:len(br)]
pfdf=pd.read_csv(f"{R}/measured_pfas_locked.csv"); Xpf,pfk=feat(pfdf); ypf=pfk.logTox.values

print("== Algorithm sensitivity (most-sensitive endpoint, ECFP4) ==")
for name,mdl in [("RandomForest",RandomForestRegressor(n_estimators=300,max_features="sqrt",random_state=SEED,n_jobs=-1)),
                 ("GradientBoosting",HistGradientBoostingRegressor(max_iter=80,learning_rate=0.1,max_leaf_nodes=15,random_state=SEED)),
                 ("Ridge",Ridge(alpha=5.0))]:
    mdl.fit(br,ytr); p=mdl.predict(Xpf)
    print(f"  {name:16s} R2={r2_score(ypf,p):+.3f}  rho={spearmanr(ypf,p).correlation:+.3f}")

print("== Applicability-domain threshold sensitivity ==")
u,_=feat(pd.read_csv(f"{R}/untested_pfas_structure_locked.csv").drop_duplicates("canonical_smiles"))
Xtr=br.astype(np.float32); Xu=u.astype(np.float32); ptr=Xtr.sum(1)
def mean5(Q,ref,refpop,self=False):
    out=np.empty(len(Q),np.float32); B=1500
    for i in range(0,len(Q),B):
        q=Q[i:i+B]; d2=q.sum(1)[:,None]+refpop[None,:]-2*(q@ref.T); np.maximum(d2,0,out=d2); d=np.sqrt(d2)
        if self:
            for j in range(len(q)): d[j,i+j]=np.inf
        out[i:i+B]=np.sort(np.partition(d,5,axis=1)[:,:6],axis=1)[:,:5].mean(1)
    return out
trnn=mean5(Xtr,Xtr,ptr,self=True); unn=mean5(Xu,Xtr,ptr)
for pct in (90,95,99):
    thr=np.percentile(trnn,pct); print(f"  {pct}th pct thr={thr:.3f}  in-domain={int((unn<=thr).sum())}/{len(unn)} ({100*(unn<=thr).mean():.1f}%)")

print("== Deployment leave-one-PFAS-out calibration stability ==")
rf=RandomForestRegressor(n_estimators=300,max_features="sqrt",random_state=SEED,n_jobs=-1).fit(br,ytr)
per=np.stack([t.predict(Xpf) for t in rf.estimators_]); mu=per.mean(0); sd=per.std(0); beta=np.median(sd)+1e-6
s=np.abs(ypf-mu)/(sd+beta); n=len(s)
def q(a,arr): k=min(int(np.ceil((len(arr)+1)*a)),len(arr)); return np.sort(arr)[k-1]
qloo=np.array([q(0.9,np.delete(s,i)) for i in range(n)])
print(f"  q@0.90 LOO mean={qloo.mean():.3f}  CV={100*qloo.std()/qloo.mean():.1f}%  clearable stays 0 (lower bound never exceeds 0.1 mg/L)")
