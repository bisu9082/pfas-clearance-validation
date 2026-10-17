"""Endpoint-swap experiment on the parent-collapsed, de-duplicated data.

Adds what reviewers asked for: a matched-sample comparison, bias-corrected R2,
variance-compression diagnostics, and a trivial baseline.
"""
import numpy as np, pandas as pd, json
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*'); UNC=rdMolStandardize.Uncharger()
SEED=20260620; rng=np.random.default_rng(SEED)

def canon(s):
    m=Chem.MolFromSmiles(str(s))
    if m is None: return None
    fr=Chem.GetMolFrags(m,asMols=True,sanitizeFrags=False)
    fl=[f for f in fr if any(a.GetSymbol()=="F" for a in f.GetAtoms())]
    big=max(fl if fl else fr,key=lambda f:f.GetNumHeavyAtoms())
    try: big=UNC.uncharge(big); Chem.SanitizeMol(big)
    except Exception: pass
    return Chem.MolToSmiles(big)
def fp(s):
    m=Chem.MolFromSmiles(str(s))
    return np.array(AllChem.GetMorganFingerprintAsBitVect(m,2,2048),dtype=np.int8) if m else None
def feats(sm):
    X,k=[],[]
    for i,s in enumerate(sm):
        v=fp(s)
        if v is not None: X.append(v); k.append(i)
    return np.vstack(X), np.array(k)

broad=pd.read_csv('cr/broad_nonpfas_structure_locked.csv')
Xb,kb=feats(broad.canonical_smiles); yb=broad.logTox.to_numpy()[kb]
print(f'broad training structures: {len(yb)}')
model=RandomForestRegressor(n_estimators=300,max_features='sqrt',random_state=SEED,n_jobs=-1).fit(Xb,yb)
ybar=float(yb.mean()); print(f'training target mean = {ybar:+.3f}\n')

rec=pd.read_csv('cr/ecotox_pfas_records.csv.gz').drop_duplicates(
      subset=['test_id','ep','conc1_mean','conc1_unit'])
src=pd.read_csv('cr/ecotox_pfas_measured.csv'); src['parent']=[canon(s) for s in src.canonical_smiles]
rec=rec.assign(test_cas=rec.test_cas.astype(str)).merge(
      src.assign(test_cas=src.test_cas.astype(str))[['test_cas','parent']],on='test_cas',how='inner')
ACUTE={'LC50','EC50','IC50'}

def build(kind,k=5,nrep=50):
    rows=[]
    for p,g in rec.groupby('parent'):
        v=g.logTox.to_numpy()
        if kind=='minimum': val=v.min()
        elif kind=='median': val=float(np.median(v))
        elif kind=='acute_median':
            a=g[g.ep.isin(ACUTE)].logTox.to_numpy()
            if len(a)==0: continue
            val=float(np.median(a))
        elif kind=='depth_min':
            if len(v)<k: continue
            val=float(np.mean([rng.choice(v,k,replace=False).min() for _ in range(nrep)]))
        rows.append((p,val,len(v)))
    return pd.DataFrame(rows,columns=['parent','y','n_rec'])

def ci(y,p,B=3000):
    r=np.random.default_rng(SEED); o=[]
    for _ in range(B):
        i=r.integers(0,len(y),len(y))
        if len(np.unique(y[i]))>2: o.append(spearmanr(y[i],p[i]).statistic)
    return float(np.nanpercentile(o,2.5)), float(np.nanpercentile(o,97.5))

res={}
print(f"{'endpoint':15s} {'n':>4s} {'R2':>7s} {'R2_bc':>7s} {'R2_mean':>8s} {'rho':>7s} {'95% CI':>16s} "
      f"{'bias':>7s} {'sd(y)':>6s} {'sdratio':>7s} {'rho_eff':>8s}")
for kind in ['minimum','median','acute_median','depth_min']:
    d=build(kind); X,k=feats(d.parent); y=d.y.to_numpy()[k]; n=d.n_rec.to_numpy()[k]
    pred=model.predict(X)
    ss=float(np.sum((y-y.mean())**2))
    r2=1-float(np.sum((y-pred)**2))/ss
    bias=float(pred.mean()-y.mean())
    r2bc=1-float(np.sum((y-(pred-bias))**2))/ss
    r2mean=1-float(np.sum((y-ybar)**2))/ss
    rho=float(spearmanr(y,pred).statistic); lo,hi=ci(y,pred)
    reff=float(spearmanr(y,np.log10(n)).statistic)
    res[kind]=dict(n=int(len(y)),r2=r2,r2_bias_corrected=r2bc,r2_trivial_mean=r2mean,
                   rho=rho,ci=[lo,hi],bias=bias,sd_y=float(y.std()),
                   sd_ratio=float(pred.std()/y.std()),rho_effort=reff)
    print(f'{kind:15s} {len(y):>4d} {r2:>7.3f} {r2bc:>7.3f} {r2mean:>8.3f} {rho:>+7.3f} '
          f'[{lo:+.2f},{hi:+.2f}] {bias:>+7.3f} {y.std():>6.2f} {pred.std()/y.std():>7.3f} {reff:>+8.3f}')

# matched sample
dmin=build('minimum').set_index('parent'); dacu=build('acute_median').set_index('parent')
common=dmin.index.intersection(dacu.index)
print(f'\n=== matched sample, the same {len(common)} parents ===')
match={}
for lab,d in [('minimum',dmin),('acute_median',dacu)]:
    dd=d.loc[common].reset_index(); X,k=feats(dd.parent); y=dd.y.to_numpy()[k]
    pred=model.predict(X); ss=float(np.sum((y-y.mean())**2))
    rho=float(spearmanr(y,pred).statistic); lo,hi=ci(y,pred)
    match[lab]=dict(n=int(len(y)),r2=1-float(np.sum((y-pred)**2))/ss,rho=rho,ci=[lo,hi])
    print(f'  {lab:14s} n={len(y)} R2={match[lab]["r2"]:+.3f} rho={rho:+.3f} [{lo:+.2f},{hi:+.2f}]')
# paired bootstrap on the difference
dm=dmin.loc[common].reset_index(); da=dacu.loc[common].reset_index()
Xm,km=feats(dm.parent); Xa,ka=feats(da.parent)
pm=model.predict(Xm); pa=model.predict(Xa)
ym=dm.y.to_numpy()[km]; ya=da.y.to_numpy()[ka]
r=np.random.default_rng(SEED); diff=[]
for _ in range(4000):
    i=r.integers(0,len(ym),len(ym))
    diff.append(spearmanr(ya[i],pa[i]).statistic-spearmanr(ym[i],pm[i]).statistic)
diff=np.array(diff)
print(f'  difference in rho = {diff.mean():+.3f}  95% CI [{np.percentile(diff,2.5):+.3f},'
      f'{np.percentile(diff,97.5):+.3f}]  P(diff<=0) = {np.mean(diff<=0):.4f}')
match['diff']=dict(mean=float(diff.mean()),ci=[float(np.percentile(diff,2.5)),
                   float(np.percentile(diff,97.5))],p=float(np.mean(diff<=0)))
json.dump(dict(endpoints=res,matched=match,broad_n=int(len(yb))),
          open('cr/endpoint_fix_recurated.json','w'),indent=2)
print('\nwrote cr/endpoint_fix_recurated.json')
