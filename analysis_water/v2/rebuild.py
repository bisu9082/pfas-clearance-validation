"""v2 rebuild: EU-TGD-faithful C6', duration-filtered C7, HC5 CIs, common subset."""
import pandas as pd, numpy as np, json
from scipy import stats
from scipy.stats import nct
rng=np.random.default_rng(20261122)

d=pd.read_csv("/home/claude/wq/results_parent.csv.gz",low_memory=False)
d=d[~d.censored.astype(bool)].copy()
d["dur"]=pd.to_numeric(d.obs_duration_mean,errors="coerce")
u=d.obs_duration_unit.astype(str).str.lower().str.strip()
d["hrs"]=d.dur*u.map({"h":1.0,"d":24.0,"wk":168.0,"mo":720.0,"min":1/60}).astype(float)
d["apical"]=d.apical.astype(bool)

# ---- taxon-specific chronic duration minima (EU TGD style) ----
def min_hours(grp):
    g=str(grp).lower()
    if "algae" in g or "plant" in g: return 72.0
    if "crustacean" in g or "invertebrate" in g or "insect" in g or "worm" in g or "mollusc" in g: return 168.0
    if "fish" in g or "amphibian" in g: return 168.0
    return 168.0
d["min_hrs"]=d.ecotox_group.map(min_hours)
d["chronic_ok"]=d.hrs.notna() & (d.hrs>=d.min_hrs)

CHRON={"NOEC","NOEL","MATC","EC10"}
AFm={3:10.0,2:50.0,1:100.0}
lad=pd.read_csv("/home/claude/wq/wqp/ladder.csv")[["parent","n_records","chronic_tl"]]
ntest=d.groupby("parent").test_id.nunique().rename("n_tests")
lad=lad.join(ntest,on="parent")

base=d[d.conc1_type=="A"]
pools={
 "C6_pub":  d[d.ep.isin(["NOEC","NOEL"])],
 "C6_ai":   base[base.ep.isin(["NOEC","NOEL"])],
 "C6_ap":   base[base.ep.isin(["NOEC","NOEL"]) & base.apical],
 "C6_prime":base[base.ep.isin(["NOEC","NOEL"]) & base.apical & base.chronic_ok],
}
out=lad.copy()
for k,p in pools.items():
    out=out.join(p.groupby("parent").logTox.min().rename(k),on="parent")
out["AF"]=out.chronic_tl.map(AFm)

print("="*78); print("C6 LADDER: progressive application of the EU TGD screen"); print("="*78)
for k in pools:
    s=out.dropna(subset=[k,"AF"])
    r,p=stats.spearmanr(np.log10(s.n_records),s[k])
    rt,_=stats.spearmanr(np.log10(s.n_tests),s[k])
    print(f"  {k:9s} n={len(s):3d}  rho(records)={r:+.3f} (p={p:.2g})   rho(tests)={rt:+.3f}")

# ---- C7 with duration filter ----
sel7=base[base.ep.isin(CHRON) & base.apical]
sel7d=sel7[sel7.chronic_ok]
def ssd(sel,tag):
    sm=sel.groupby(["parent","latin_name"]).logTox.mean().reset_index()
    g=sm.groupby("parent").logTox
    t=pd.DataFrame({"n_sp":g.size(),"mu":g.mean(),"sd":g.std()}).reset_index()
    t["HC5"]=t.mu+stats.norm.ppf(.05)*t.sd
    t.loc[(t.n_sp<5)|t.sd.isna(),"HC5"]=np.nan
    t["minSMCV"]=g.min().values
    return t.rename(columns={c:f"{c}_{tag}" for c in ["n_sp","mu","sd","HC5","minSMCV"]})
a=ssd(sel7,"nd"); b=ssd(sel7d,"d")
out=out.merge(a,on="parent",how="left").merge(b,on="parent",how="left")

print("\n"+"="*78); print("C7 SSD: without vs with the chronic-duration requirement"); print("="*78)
for tag,nm in [("nd","no duration filter (as published)"),("d","duration-filtered")]:
    s=out.dropna(subset=[f"HC5_{tag}"])
    r,p=stats.spearmanr(np.log10(s.n_records),s[f"HC5_{tag}"])
    print(f"  {nm:34s} n={len(s):3d}  rho={r:+.3f}  p={p:.4f}")
    sm=out.dropna(subset=[f"minSMCV_{tag}"])
    r2,p2=stats.spearmanr(np.log10(sm.n_records),sm[f"minSMCV_{tag}"])
    print(f"    C7a lowest apical SMCV           n={len(sm):3d}  rho={r2:+.3f}  p={p2:.2g}")

# ---- sigma decomposition ----
print("\n"+"="*78); print("MECHANISM: mu vs sigma channel"); print("="*78)
for tag in ("nd","d"):
    s=out.dropna(subset=[f"HC5_{tag}"]); ln=np.log10(s.n_records)
    z=stats.norm.ppf(.05)
    bmu=np.polyfit(ln,s[f"mu_{tag}"],1)[0]; bsd=np.polyfit(ln,s[f"sd_{tag}"],1)[0]
    print(f"  [{tag}] n={len(s)}  rho(mu)={stats.spearmanr(ln,s[f'mu_{tag}'])[0]:+.3f}  "
          f"rho(sd)={stats.spearmanr(ln,s[f'sd_{tag}'])[0]:+.3f}")
    print(f"       d(mu)={bmu:+.3f}  z*d(sd)={z*bsd:+.3f}  total={bmu+z*bsd:+.3f}  "
          f"sigma share={abs(z*bsd)/abs(bmu+z*bsd)*100:.0f}%")

# ---- HC5 confidence intervals ----
print("\n"+"="*78); print("HC5 95% CI (noncentral-t tolerance limit), duration-filtered"); print("="*78)
z=stats.norm.ppf(.05)
rows=[]
for _,r in out.dropna(subset=["HC5_d"]).sort_values("n_sp_d").iterrows():
    n=int(r.n_sp_d); s=r.sd_d; dnc=-z*np.sqrt(n)
    lo=r.mu_d-nct.ppf(.95,n-1,dnc)/np.sqrt(n)*s
    hi=r.mu_d-nct.ppf(.05,n-1,dnc)/np.sqrt(n)*s
    rows.append(dict(parent=r.parent,n_sp=n,HC5=r.HC5_d,lo=lo,hi=hi,width=hi-lo,
                     HC5_ngL=10**r.HC5_d*1e6,lo_ngL=10**lo*1e6,hi_ngL=10**hi*1e6))
ci=pd.DataFrame(rows)
print(ci[["n_sp","HC5_ngL","lo_ngL","hi_ngL","width"]].to_string(index=False,
      float_format=lambda x:f"{x:,.3g}"))
ci.to_csv("hc5_ci.csv",index=False)
out.to_csv("ladder_v2.csv",index=False)
print("\nwrote ladder_v2.csv, hc5_ci.csv")
