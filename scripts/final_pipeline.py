"""Complete rebuild on the acute-median endpoint, training and external alike.

Everything downstream of the endpoint is the manuscript's published protocol:
ECFP4 + random forest, 5x5 broad CV, external PFAS validation, normalised
split-conformal (marginal and PFAS-conditional), and inventory clearance.
"""
import sys, json, time
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import RepeatedKFold
from sklearn.neighbors import NearestNeighbors
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')
BROAD_SMI, PFAS_REC, PFASMAP, INV, OUT = sys.argv[1:6]
SEED, NOMINAL, THRESH, NREP = 20260620, 0.90, -1.0, 40
UNC = rdMolStandardize.Uncharger()

def canon(s):
    m = Chem.MolFromSmiles(str(s))
    if m is None: return None
    fr = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    fl = [f for f in fr if any(a.GetSymbol()=="F" for a in f.GetAtoms())]
    big = max(fl if fl else fr, key=lambda f: f.GetNumHeavyAtoms())
    try:
        big = UNC.uncharge(big); Chem.SanitizeMol(big)
    except Exception: pass
    return Chem.MolToSmiles(big)

def fp1(s):
    m = Chem.MolFromSmiles(str(s))
    return np.array(AllChem.GetMorganFingerprintAsBitVect(m,2,2048),dtype=np.int8) if m else None

def feats(sm):
    X,k=[],[]
    for i,s in enumerate(sm):
        v=fp1(s)
        if v is not None: X.append(v); k.append(i)
    return np.vstack(X), np.array(k)

def rf(n=300): return RandomForestRegressor(n_estimators=n,max_features="sqrt",
                                            random_state=SEED,n_jobs=-1)
def met(y,p):
    ss=float(np.sum((y-np.mean(y))**2))
    return (1-float(np.sum((y-p)**2))/ss if ss>0 else np.nan, float(spearmanr(y,p).statistic))
def disp(model,X,chunk=4000):
    o=[]
    for i in range(0,len(X),chunk):
        o.append(np.std([t.predict(X[i:i+chunk]) for t in model.estimators_],axis=0))
    return np.concatenate(o)
def nconf(model,Xc,yc,Xt,nominal=NOMINAL):
    sc=disp(model,Xc); b=float(np.median(sc))+1e-6
    s=np.abs(yc-model.predict(Xc))/(sc+b)
    n=len(s); k=int(np.ceil((n+1)*nominal))
    if k>n: return None,None
    q=np.sort(s)[k-1]; st=disp(model,Xt)
    return q*(st+b), model.predict(Xt)

log=lambda *a: print(*a, flush=True)

# ---------------- external PFAS: acute median -----------------------------
rec=pd.read_csv(PFAS_REC, dtype={"test_cas":str})
ACUTE={"LC50","EC50","IC50"}
a=rec[rec.ep.isin(ACUTE)]
pf=(a.groupby("test_cas").logTox.median().rename("y").reset_index())
smap=pd.read_csv(PFASMAP,dtype={"cas_norm":str})[["cas_norm","canonical_smiles","name"]].drop_duplicates("cas_norm")
pf=pf.merge(smap,left_on="test_cas",right_on="cas_norm").dropna(subset=["canonical_smiles"])
pf["parent"]=pf.canonical_smiles.map(canon)
pf=pf.dropna(subset=["parent"]).groupby("parent",as_index=False).agg(y=("y","median"),name=("name","first"))
log(f"external PFAS (acute median, structure-locked): {len(pf)}")

# ---------------- broad training: acute median, PFAS removed --------------
bs=pd.read_csv(BROAD_SMI, dtype={"test_cas":str})
bs=bs[bs.smiles.notna()].copy()
bs["parent"]=bs.smiles.map(canon)
bs=bs.dropna(subset=["parent"])
pfas_parents=set(pf.parent)
inv=pd.read_csv(INV,usecols=["canonical_smiles"])
inv_parents=set(x for x in inv.canonical_smiles.map(canon) if x)
before=len(bs)
bs=bs[~bs.parent.isin(pfas_parents|inv_parents)]
log(f"broad: {before} -> {len(bs)} after removing PFAS structures")
broad=bs.groupby("parent",as_index=False).agg(y=("acute_median","median"),n_rec=("n_rec","sum"))
log(f"broad training (structure-locked): {len(broad)}")

Xb,kb=feats(broad.parent); yb=broad.y.to_numpy()[kb]
Xe,ke=feats(pf.parent);    pf=pf.iloc[ke].reset_index(drop=True); ye=pf.y.to_numpy()
log(f"featurised: broad {Xb.shape}, external {Xe.shape}")

# ---------------- broad cross-validation ----------------------------------
t0=time.time(); folds=[]
for tr,te in RepeatedKFold(n_splits=5,n_repeats=5,random_state=SEED).split(Xb):
    m=rf(160).fit(Xb[tr],yb[tr]); folds.append(met(yb[te],m.predict(Xb[te]))[0])
log(f"\nbroad 5x5 CV R2 = {np.mean(folds):+.3f} +/- {np.std(folds):.3f}  ({time.time()-t0:.0f}s)")

# ---------------- external validation --------------------------------------
full=rf().fit(Xb,yb); pe=full.predict(Xe)
r2,rho=met(ye,pe)
rng=np.random.default_rng(SEED)
bs_=[spearmanr(ye[i],pe[i]).statistic for i in (rng.integers(0,len(ye),len(ye)) for _ in range(4000))]
lo,hi=np.nanpercentile(bs_,[2.5,97.5])
log(f"external PFAS: n={len(ye)}  R2 = {r2:+.3f}   rho = {rho:+.3f} [{lo:+.2f}, {hi:+.2f}]")

# ---------------- conformal -------------------------------------------------
marg=[]; cond=[]
for rep in range(NREP):
    p=rng.permutation(len(yb)); c=int(0.8*len(yb))
    m=rf().fit(Xb[p[:c]],yb[p[:c]])
    h,pr=nconf(m,Xb[p[c:]],yb[p[c:]],Xe)
    if h is not None: marg.append((np.mean(np.abs(ye-pr)<=h), np.median(2*h)))
    q=rng.permutation(len(ye)); half=len(ye)//2
    h2,pr2=nconf(full,Xe[q[:half]],ye[q[:half]],Xe[q[half:]])
    if h2 is not None: cond.append((np.mean(np.abs(ye[q[half:]]-pr2)<=h2), np.median(2*h2)))
M=np.array(marg); C=np.array(cond)
log(f"marginal    coverage {M[:,0].mean():.3f}  width {M[:,1].mean():5.2f}")
log(f"conditional coverage {C[:,0].mean():.3f}  width {C[:,1].mean():5.2f}")

# ---------------- clearance -------------------------------------------------
iv=pd.read_csv(INV); iv["parent"]=iv.canonical_smiles.map(canon)
iv=iv.dropna(subset=["parent"]).drop_duplicates("parent")
Xi,ki=feats(iv.parent); iv=iv.iloc[ki].reset_index(drop=True)
nn=NearestNeighbors(n_neighbors=5).fit(Xb)
ds,_=nn.kneighbors(Xb,n_neighbors=6); thr=float(np.percentile(ds[:,1:].mean(1),95))
di,_=nn.kneighbors(Xi); dom=di.mean(1)<=thr
h,pi=nconf(full,Xe,ye,Xi)
lower=pi-h; cleared=(lower>THRESH)&dom
log(f"\ninventory {len(iv)}  in-domain {dom.sum()}  median width {np.median(2*h[dom]):.2f}")
log(f"CLEARED at 90% / 0.1 mg/L: {int(cleared.sum())}  ({cleared.sum()/dom.sum():.2%} of in-domain)")

json.dump(dict(n_broad=int(len(yb)),n_ext=int(len(ye)),cv_r2=float(np.mean(folds)),
    cv_sd=float(np.std(folds)),ext_r2=float(r2),ext_rho=float(rho),
    rho_lo=float(lo),rho_hi=float(hi),marg_cov=float(M[:,0].mean()),
    marg_w=float(M[:,1].mean()),cond_cov=float(C[:,0].mean()),cond_w=float(C[:,1].mean()),
    in_domain=int(dom.sum()),cleared=int(cleared.sum())),
    open(f"{OUT}/final_metrics.json","w"),indent=2)
iv.assign(pred=pi,lower=lower,upper=pi+h,in_domain=dom,cleared=cleared) \
  .to_csv(f"{OUT}/final_inventory_predictions.csv",index=False)
log("\nwrote final_metrics.json / final_inventory_predictions.csv")
