import pandas as pd, numpy as np, json
from scipy import stats
R=pd.read_csv("/home/claude/wq/results_parent.csv.gz")
APICAL={"MOR","GRO","REP","DVP","POP"}; C7EP={"NOEC","NOEL","MATC","EC10"}
sel=R[R.ep.isin(C7EP)&R.effect.isin(APICAL)&(R.conc1_type=="A")].copy()
smcv=sel.groupby(["parent","latin_name"]).logTox.mean().reset_index()
L=pd.read_csv("ladder.csv")

print("="*92); print("S1  SSD FIT DIAGNOSTICS AND DISTRIBUTIONAL SENSITIVITY"); print("="*92)
rows=[]
for p,g in smcv.groupby("parent"):
    v=np.sort(g.logTox.values); n=len(v)
    if n<5: continue
    mu,sd=v.mean(),v.std(ddof=1)
    ad=stats.anderson(v,"norm"); crit5=ad.critical_values[2]
    sh=stats.shapiro(v) if n>=3 else (np.nan,np.nan)
    hc5_ln=mu+stats.norm.ppf(.05)*sd
    # log-logistic: fit location/scale by MLE on logistic
    lo,sc=stats.logistic.fit(v)
    hc5_ll=stats.logistic.ppf(.05,lo,sc)
    rows.append(dict(parent=p,n=n,mu=mu,sd=sd,AD=ad.statistic,AD_crit5=crit5,
        AD_reject=ad.statistic>crit5,shapiro_p=sh[1],
        HC5_lognormal=hc5_ln,HC5_loglogistic=hc5_ll,diff_log=hc5_ll-hc5_ln))
D=pd.DataFrame(rows)
pd.set_option("display.width",220)
print(D[["n","mu","sd","AD","AD_crit5","AD_reject","shapiro_p","HC5_lognormal","HC5_loglogistic","diff_log"]]
      .to_string(index=False,float_format=lambda x:f"{x:.3f}"))
print(f"\n  normality rejected at 5%: {int(D.AD_reject.sum())} / {len(D)} substances")
print(f"  |HC5 log-logistic - log-normal|: median {D.diff_log.abs().median():.3f} log10, max {D.diff_log.abs().max():.3f}")
D.to_csv("si_ssd_diag.csv",index=False)

print("\n"+"="*92); print("S2  SENSITIVITY OF THE EFFORT CORRELATION TO SSD CHOICES"); print("="*92)
rng=np.random.default_rng(1121)
def rho(x,y):
    k=(~pd.isna(x))&(~pd.isna(y)); x,y=np.asarray(x)[k],np.asarray(y)[k]
    r,_=stats.spearmanr(x,y)
    bs=[stats.spearmanr(x[i],y[i])[0] for i in (rng.integers(0,len(y),len(y)) for _ in range(4000))]
    return r,*np.nanpercentile(bs,[2.5,97.5]),len(y)
for minsp in (3,5,8,10):
    rr=[]
    for p,g in smcv.groupby("parent"):
        v=g.logTox.values; n=len(v)
        if n<minsp: continue
        mu,sd=v.mean(),v.std(ddof=1)
        rr.append(dict(parent=p,n=n,hc5=mu+stats.norm.ppf(.05)*sd))
    t=pd.DataFrame(rr).merge(L[["parent","n_records"]],on="parent")
    if len(t)<5: print(f"  min species {minsp:2d}: n={len(t)} too few"); continue
    r,lo,hi,nn=rho(np.log10(t.n_records),t.hc5)
    print(f"  min species {minsp:2d}  (log-normal)   n={nn:3d}  rho={r:+.3f}  CI [{lo:+.3f},{hi:+.3f}]")
for minsp in (5,):
    rr=[]
    for p,g in smcv.groupby("parent"):
        v=g.logTox.values
        if len(v)<minsp: continue
        lo_,sc=stats.logistic.fit(v)
        rr.append(dict(parent=p,hc5=stats.logistic.ppf(.05,lo_,sc)))
    t=pd.DataFrame(rr).merge(L[["parent","n_records"]],on="parent")
    r,lo,hi,nn=rho(np.log10(t.n_records),t.hc5)
    print(f"  min species  5  (log-logistic) n={nn:3d}  rho={r:+.3f}  CI [{lo:+.3f},{hi:+.3f}]")
