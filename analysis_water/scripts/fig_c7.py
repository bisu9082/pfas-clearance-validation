import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
C=pd.read_csv("/home/claude/wq/c7_ssd.csv"); L=pd.read_csv("ladder.csv")
ex=pd.read_csv("exceedance_c7.csv")
m=C.merge(L[["parent","n_records","C1_min_all","C3_min_acute","C5_spmed_acute",
             "C6_min_noec","AF_C6"]],on="parent")
rng=np.random.default_rng(20261121)
def rho(x,y):
    k=(~pd.isna(x))&(~pd.isna(y)); x,y=np.asarray(x)[k],np.asarray(y)[k]
    r,_=stats.spearmanr(x,y)
    bs=[stats.spearmanr(x[i],y[i])[0] for i in (rng.integers(0,len(y),len(y)) for _ in range(4000))]
    return r,*np.nanpercentile(bs,[2.5,97.5]),len(y)
x=np.log10(m.n_records)
ITEMS=[("lowest, all endpoints\n(as published)",rho(x,m.C1_min_all)),
       ("lowest acute",rho(x,m.C3_min_acute)),
       ("lowest chronic NOEC,\nEU-TGD AF",rho(x,np.log10((10**m.C6_min_noec)/m.AF_C6))),
       ("lowest apical chronic\nSMCV (screened)",rho(x,m.C7_min)),
       ("SSD HC5 of apical\nchronic SMCVs",rho(np.log10(m[m.C7_HC5.notna()].n_records),
                                               m[m.C7_HC5.notna()].C7_HC5)),
       ("median acute,\nspecies-resolved",rho(x,m.C5_spmed_acute))]
fig,ax=plt.subplots(2,2,figsize=(13.6,9.4))

# A exceedance
a=ax[0,0]; lb=["lowest,\nall endpoints","lowest\nacute","lowest chronic\nNOEC +AF",
               "lowest apical\nSMCV +AF","SSD HC5\n(regulatory)","median acute,\nspecies"]
o=[0,1,3,4,5,2]
y1=ex.pct_samples.values[o]; y2=ex.pct_sites.values[o]
xs=np.arange(6); w=.38
a.bar(xs-w/2,y1,w,color="#2c5f8a",label="samples"); a.bar(xs+w/2,y2,w,color="#9ec5e0",edgecolor="#2c5f8a",label="sites")
for i,(u,v) in enumerate(zip(y1,y2)):
    a.text(i-w/2,u+1.7,f"{u:.1f}",ha="center",fontsize=7.6); a.text(i+w/2,v+1.7,f"{v:.1f}",ha="center",fontsize=7.6)
a.set_xticks(xs); a.set_xticklabels(lb,fontsize=7.3); a.set_ylim(0,112)
a.set_ylabel("% exceeding the derived criterion"); a.legend(fontsize=8,frameon=False,loc="upper right")
a.set_title("A  21,346 measurements, 1,180 US sites: the verdict depends on the convention",fontsize=9.4,loc="left")

# B forest
b=ax[0,1]; ys=np.arange(len(ITEMS))[::-1]
for yy,(lbl,(r,lo,hi,n)) in zip(ys,ITEMS):
    sig = not (lo<=0<=hi)
    c="#b4442e" if sig else "#2f7d4f"
    b.plot([lo,hi],[yy,yy],c=c,lw=2.3,solid_capstyle="round")
    b.plot(r,yy,"o",c=c,ms=8,zorder=3)
    b.text(.06,yy+.26,f"n={n}",fontsize=7,color="#555")
b.axvline(0,c="#888",lw=1,ls="--")
b.set_yticks(ys); b.set_yticklabels([i[0] for i in ITEMS],fontsize=8)
b.set_xlabel(r"Spearman $\rho$ between testing effort (log records) and the derived criterion")
b.set_xlim(-1.02,.62)
b.set_title("B  Screening to regulatory standards does not remove the artefact\n"
            "red = CI excludes zero; green = effort-independent",fontsize=9.4,loc="left")

# C mechanism
c=ax[1,0]
c.scatter(m.n_records,m.n_sp_C7,s=40,c="#2c5f8a",alpha=.8,edgecolor="w",lw=.5)
r,lo,hi,n=rho(x,m.n_sp_C7)
sl,ic,_,_,_=stats.linregress(x,np.log10(m.n_sp_C7.clip(lower=.5)))
gx=np.linspace(x.min(),x.max(),50); c.plot(10**gx,10**(sl*gx+ic),c="#2c5f8a",lw=1.8)
c.set_xscale("log"); c.set_yscale("log")
c.set_xlabel("ECOTOX records for the substance"); c.set_ylabel("species surviving the regulatory screen")
c.set_title(f"C  Why: more study means a longer left tail to the SSD\n"
            rf"$\rho$ = {r:+.2f}, 95% CI [{lo:+.2f}, {hi:+.2f}]",fontsize=9.4,loc="left")

# D matched pair
d=ax[1,1]; mm=pd.read_pickle("m_ladder.pkl")
mm["subst"]=mm.cas.map({"335-67-1":"PFOA","375-85-9":"PFHpA"})
for i,(s,col) in enumerate([("PFOA","#b4442e"),("PFHpA","#2c5f8a")]):
    g=mm[mm.subst==s]
    d.scatter(np.random.default_rng(7+i).normal(i,.075,len(g)),g.ngL,s=6,alpha=.2,c=col)
    p6=(10**g.C6_min_noec.iloc[0])*1e6/g.AF_C6.iloc[0]
    d.hlines(p6,i-.3,i+.3,color=col,lw=2.6,zorder=4)
    d.text(i+.33,p6,f"PNEC {p6:,.3g} ng/L",fontsize=8,color=col,va="center")
    d.text(i,3e4,f"{s}\n{int(g.n_records.iloc[0]):,} records\n{100*(g.ngL>p6).mean():.1f}% exceed",
           ha="center",fontsize=8.6,color=col,weight="bold")
d.set_yscale("log"); d.set_xticks([0,1]); d.set_xticklabels(["PFOA","PFHpA"])
d.set_xlim(-.55,2.0); d.set_ylim(1e-2,3e5); d.set_ylabel("measured concentration (ng/L)")
d.set_title("D  Two PFAS measured at the same 819 sites, opposite verdicts",fontsize=9.4,loc="left")
plt.tight_layout(); plt.savefig("fig_c7.png",dpi=160)
for lbl,(r,lo,hi,n) in ITEMS: print(f"{lbl.replace(chr(10),' '):42s} rho={r:+.3f} [{lo:+.3f},{hi:+.3f}] n={n}")
