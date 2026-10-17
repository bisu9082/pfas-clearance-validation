"""Class-resolved external transfer on the ACUTE MEDIAN (the endpoint this paper endorses)."""
import numpy as np, pandas as pd, json
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')
SEED=20260620
def feats(sm):
    X,k=[],[]
    for i,s in enumerate(sm):
        m=Chem.MolFromSmiles(str(s))
        if m is not None:
            X.append(np.array(AllChem.GetMorganFingerprintAsBitVect(m,2,2048),dtype=np.int8)); k.append(i)
    return np.vstack(X), np.array(k)
broad=pd.read_csv('cr/broad_nonpfas_structure_locked.csv')
Xb,kb=feats(broad.canonical_smiles)
model=RandomForestRegressor(n_estimators=300,max_features='sqrt',random_state=SEED,n_jobs=-1
       ).fit(Xb,broad.logTox.to_numpy()[kb])
d=pd.read_csv('cr/oecd_103set.csv')
d=d[d.oecd_pfas & d.logTox_median.notna()].reset_index(drop=True)
X,k=feats(d.canonical_smiles); d=d.iloc[k].reset_index(drop=True)
d['pred']=model.predict(X)
def ci(y,p,B=4000):
    r=np.random.default_rng(SEED); o=[]
    for _ in range(B):
        i=r.integers(0,len(y),len(y))
        if len(np.unique(y[i]))>2: o.append(spearmanr(y[i],p[i]).statistic)
    return float(np.nanpercentile(o,2.5)), float(np.nanpercentile(o,97.5))
def power_floor(n,alpha=0.05,power=0.80):
    from scipy.stats import norm
    if n<5: return np.nan
    z=(norm.ppf(1-alpha/2)+norm.ppf(power))/np.sqrt(n-3)
    return float(np.tanh(z))
rows=[]
y=d.logTox_median.to_numpy(); p=d.pred.to_numpy()
lo,hi=ci(y,p)
rows.append(dict(cls='All OECD-2021 PFAS',n=len(d),rho=float(spearmanr(y,p).statistic),
                 lo=lo,hi=hi,floor=power_floor(len(d)),eff=float(np.log10(d.n_records).mean())))
for c,g in d.groupby('oecd_class'):
    if len(g)<3: 
        rows.append(dict(cls=c,n=len(g),rho=np.nan,lo=np.nan,hi=np.nan,
                         floor=np.nan,eff=float(np.log10(g.n_records).mean()))); continue
    yy=g.logTox_median.to_numpy(); pp=g.pred.to_numpy()
    l,h=ci(yy,pp)
    rows.append(dict(cls=c,n=len(g),rho=float(spearmanr(yy,pp).statistic),lo=l,hi=h,
                     floor=power_floor(len(g)),eff=float(np.log10(g.n_records).mean())))
f=pd.DataFrame(rows).sort_values('rho')
f.to_csv('cr/forest_acute_median.csv',index=False)
print(f.to_string(index=False,float_format=lambda x:f'{x:6.3f}'))
sub=f[f.cls.ne('All OECD-2021 PFAS') & f.rho.notna()]
print(f'\nclasses with an interval excluding zero: {int(((sub.lo>0)|(sub.hi<0)).sum())} of {len(sub)}')
r,pv=spearmanr(sub.rho,sub.eff)
print(f'class rho vs class mean log10(records): rho={r:+.3f} p={pv:.3f}  '
      f'(tests whether class transfer merely tracks testing effort)')
json.dump(dict(rows=rows,n_excluding_zero=int(((sub.lo>0)|(sub.hi<0)).sum()),
               rho_vs_effort=[float(r),float(pv)]),open('cr/forest.json','w'),indent=2)
