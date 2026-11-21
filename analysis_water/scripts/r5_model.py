"""Does a model fitted to an effort-contaminated target learn the effort?"""
import sys,functools
print=functools.partial(print,flush=True)
import numpy as np, pandas as pd, json
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from scipy import stats
RDLogger.DisableLog('rdApp.*')
B=pd.read_csv("/home/claude/pfasv32/canonical_results/ecotox_acute_smiles.csv")
B=B[B.smiles.notna()].copy()
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=1024)
X=[];keep=[]
for i,s in enumerate(B.smiles):
    m=Chem.MolFromSmiles(str(s))
    if m is None: continue
    X.append(np.frombuffer(gen.GetFingerprintAsNumPy(m).tobytes(),dtype=np.uint8))
    keep.append(i)
X=np.array(X,dtype=np.float32); B=B.iloc[keep].reset_index(drop=True)
print("chemicals with a parsable structure:",len(B),"  fingerprint dim:",X.shape[1])
B["log_n"]=np.log10(B.n_rec)

def cv(y,tag,n_splits=5):
    kf=KFold(n_splits=n_splits,shuffle=True,random_state=1121)
    pred=np.zeros(len(y))
    for tr,te in kf.split(X):
        rf=RandomForestRegressor(n_estimators=150,min_samples_leaf=3,n_jobs=-1,random_state=0)
        rf.fit(X[tr],y[tr]); pred[te]=rf.predict(X[te])
    r2=1-np.sum((y-pred)**2)/np.sum((y-y.mean())**2)
    rho=stats.spearmanr(y,pred)[0]
    print(f"  {tag:46s} R2={r2:+.3f}  rho={rho:+.3f}")
    return pred,r2,rho

print("\n"+"="*94)
print("DOES MOLECULAR STRUCTURE ENCODE TESTING EFFORT?")
print("="*94)
p_eff,r2_eff,rho_eff=cv(B.log_n.values,"structure -> log10(number of records)")

print("\n"+"="*94)
print("MODELS OF THE TWO TARGETS")
print("="*94)
p_min,r2_min,rho_min=cv(B.acute_min.values,"structure -> MINIMUM acute value")
p_med,r2_med,rho_med=cv(B.acute_median.values,"structure -> MEDIAN acute value")

print("\n"+"="*94)
print("HOW MUCH OF EACH MODEL IS EFFORT?")
print("="*94)
for nm,y,p in [("MINIMUM",B.acute_min.values,p_min),("MEDIAN",B.acute_median.values,p_med)]:
    # correlation between what the model predicts and the chemical's record count
    r=stats.spearmanr(p,B.log_n)[0]
    # performance after removing the effort component from the target
    res=y-np.polyval(np.polyfit(B.log_n,y,1),B.log_n)
    pr,r2r,rhor=cv(res,f"structure -> {nm} with the effort trend removed")
    print(f"  {nm}: model prediction vs record count rho={r:+.3f}"
          f"   |  R2 falls {globals()['r2_'+('min' if nm=='MINIMUM' else 'med')]:+.3f} -> {r2r:+.3f}"
          f"  ({100*(1-r2r/globals()['r2_'+('min' if nm=='MINIMUM' else 'med')]):.0f}% of the fit was the effort trend)")

print("\n"+"="*94)
print("PERFORMANCE STRATIFIED BY HOW MUCH THE CHEMICAL HAS BEEN TESTED")
print("="*94)
bins=[(1,2),(2,5),(5,20),(20,100),(100,10**9)]
print(f"  {'records':>12} {'n':>6} | {'MIN  R2':>9} {'MIN rho':>9} | {'MED  R2':>9} {'MED rho':>9}")
strat=[]
for lo,hi in bins:
    k=(B.n_rec>=lo)&(B.n_rec<hi)
    if k.sum()<40: continue
    row={"bin":f"{lo}-{hi if hi<10**9 else '+'}","n":int(k.sum())}
    for nm,y,p in [("MIN",B.acute_min.values,p_min),("MED",B.acute_median.values,p_med)]:
        yy,pp=y[k.values],p[k.values]
        row[nm+"_r2"]=float(1-np.sum((yy-pp)**2)/np.sum((yy-yy.mean())**2))
        row[nm+"_rho"]=float(stats.spearmanr(yy,pp)[0])
    print(f"  {row['bin']:>12} {row['n']:>6} | {row['MIN_r2']:>9.3f} {row['MIN_rho']:>9.3f}"
          f" | {row['MED_r2']:>9.3f} {row['MED_rho']:>9.3f}")
    strat.append(row)
json.dump(dict(r2_effort=r2_eff,rho_effort=rho_eff,r2_min=r2_min,r2_med=r2_med,
               rho_min=rho_min,rho_med=rho_med,strata=strat),open("r5_model.json","w"),indent=1)
B.assign(pred_min=p_min,pred_med=p_med,pred_eff=p_eff).to_csv("r5_pred.csv",index=False)
