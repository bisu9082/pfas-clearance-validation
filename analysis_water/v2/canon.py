import pandas as pd, numpy as np, json
from scipy import stats
out=pd.read_csv("ladder_v2.csv"); ci=pd.read_csv("hc5_ci.csv")
R=pd.read_pickle("/home/claude/wq/broad/broad_records.pkl")[["test_cas","species_number","logTox"]].dropna()
rng=np.random.default_rng(20261122)
g=R.groupby("test_cas")
S=pd.DataFrame({"n_rec":g.size(),"mn":g.logTox.min(),"md":g.logTox.median(),"sd":g.logTox.std()}).reset_index()
S["logn"]=np.log10(S.n_rec)
SP=R.groupby(["test_cas","species_number"],as_index=False).logTox.mean()
def partial(x,y,z):
    rx,ry,rz=(stats.rankdata(v) for v in (x,y,z))
    res=lambda a,b: a-np.polyval(np.polyfit(b,a,1),b)
    return np.corrcoef(res(rx,rz),res(ry,rz))[0,1]
# species-equalised
res={}
for s in (3,5,10):
    el=SP.groupby("test_cas").size(); el=el[el>=s].index
    pool={c:v.logTox.values for c,v in SP[SP.test_cas.isin(el)].groupby("test_cas")}
    base=S[S.test_cas.isin(el)].set_index("test_cas").loc[list(pool)]
    ln=np.log10(base.n_rec.values)
    b=stats.spearmanr(ln,[pool[c].min() for c in pool])[0]
    e=np.mean([stats.spearmanr(ln,[rng.choice(pool[c],s,replace=False).min() for c in pool])[0] for _ in range(40)])
    res[f"species_equalised_s{s}"]=dict(n=len(pool),baseline=round(b,3),equalised=round(e,3),
        depth_component=round(e-b,3),pct_removed=round(abs(e-b)/abs(b)*100,1))
el=S[S.n_rec>=5].test_cas.values
pool={c:v.logTox.values for c,v in R[R.test_cas.isin(el)].groupby("test_cas")}
b2=S[S.test_cas.isin(el)].set_index("test_cas").loc[list(pool)]
eq=np.mean([[rng.choice(pool[c],5,replace=False).min() for c in pool] for _ in range(40)],axis=0)
msk=np.isfinite(b2.sd.values)
res["sd_control"]=dict(
  equalised_min_rho=round(stats.spearmanr(np.log10(b2.n_rec),eq)[0],3),
  partial_given_sd=round(partial(np.log10(b2.n_rec.values)[msk],eq[msk],b2.sd.values[msk]),3),
  rho_logn_sd=round(stats.spearmanr(S.logn,S.sd.fillna(0))[0],3))
z=stats.norm.ppf(.05)
for tag in ("nd","d"):
    s=out.dropna(subset=[f"HC5_{tag}"]); ln=np.log10(s.n_records)
    bmu=np.polyfit(ln,s[f"mu_{tag}"],1)[0]; bsd=np.polyfit(ln,s[f"sd_{tag}"],1)[0]
    res[f"hc5_decomp_{tag}"]=dict(n=len(s),rho_mu=round(stats.spearmanr(ln,s[f"mu_{tag}"])[0],3),
      rho_sd=round(stats.spearmanr(ln,s[f"sd_{tag}"])[0],3),d_mu=round(bmu,3),
      z_d_sd=round(z*bsd,3),total=round(bmu+z*bsd,3),sigma_share_pct=round(abs(z*bsd)/abs(bmu+z*bsd)*100))
res["hc5_ci_width_log"]=dict(min=round(ci.width.min(),2),max=round(ci.width.max(),2),median=round(ci.width.median(),2))
res["n_ssd_able_duration_filtered"]=int(out.HC5_d.notna().sum())
json.dump(res,open("canonical_v2.json","w"),indent=1)
print(json.dumps(res,indent=1))
