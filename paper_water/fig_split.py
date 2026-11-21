import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({"font.size":11,"axes.titlesize":11.5,"axes.labelsize":11,
 "xtick.labelsize":10,"ytick.labelsize":10,"legend.fontsize":10})
C=pd.read_csv("/home/claude/wq/c7_ssd.csv"); L=pd.read_csv("/home/claude/wq/wqp/ladder.csv")
ex=pd.read_csv("/home/claude/wq/wqp/exceedance_c7.csv")
m=C.merge(L[["parent","n_records","C1_min_all","C3_min_acute","C5_spmed_acute","C6_min_noec","AF_C6"]],on="parent")
rng=np.random.default_rng(20261121)
def rho(x,y):
    k=(~pd.isna(x))&(~pd.isna(y)); x,y=np.asarray(x)[k],np.asarray(y)[k]
    r,_=stats.spearmanr(x,y)
    bs=[stats.spearmanr(x[i],y[i])[0] for i in (rng.integers(0,len(y),len(y)) for _ in range(4000))]
    return r,*np.nanpercentile(bs,[2.5,97.5]),len(y)
x=np.log10(m.n_records)

# ---- Fig 2: exceedance ladder
fig,a=plt.subplots(figsize=(8.6,4.6))
lb=["lowest,\nall endpoints","lowest\nacute","lowest chronic\nNOEC + AF",
    "lowest apical\nSMCV + AF","SSD HC5\n(regulatory)","median acute,\nspecies"]
o=[0,1,3,4,5,2]; y1=ex.pct_samples.values[o]; y2=ex.pct_sites.values[o]
xs=np.arange(6); w=.38
a.bar(xs-w/2,y1,w,color="#2c5f8a",label="samples (n = 21,346)")
a.bar(xs+w/2,y2,w,color="#9ec5e0",edgecolor="#2c5f8a",label="stations (n = 1,180)")
for i,(u,v) in enumerate(zip(y1,y2)):
    a.text(i-w/2,u+1.6,f"{u:.1f}",ha="center",fontsize=9); a.text(i+w/2,v+1.6,f"{v:.1f}",ha="center",fontsize=9)
a.set_xticks(xs); a.set_xticklabels(lb,fontsize=9); a.set_ylim(0,114)
a.set_ylabel("% exceeding the derived criterion"); a.legend(frameon=False,loc="upper right",fontsize=9.5)
a.set_xlabel("derivation convention")
plt.tight_layout(); plt.savefig("fig_02_exceedance.png",dpi=300,bbox_inches="tight"); plt.close()

# ---- Fig 3: forest + mechanism
ITEMS=[("lowest, all endpoints\n(as published)",rho(x,m.C1_min_all)),
       ("lowest acute",rho(x,m.C3_min_acute)),
       ("lowest chronic NOEC,\nEU-TGD AF",rho(x,np.log10((10**m.C6_min_noec)/m.AF_C6))),
       ("lowest apical chronic\nSMCV (screened)",rho(x,m.C7_min)),
       ("SSD HC5 of apical\nchronic SMCVs",rho(np.log10(m[m.C7_HC5.notna()].n_records),m[m.C7_HC5.notna()].C7_HC5)),
       ("median acute,\nspecies-resolved",rho(x,m.C5_spmed_acute))]
fig,ax=plt.subplots(1,2,figsize=(12.6,4.9))
b=ax[0]; ys=np.arange(len(ITEMS))[::-1]
for yy,(lbl,(r,lo,hi,n)) in zip(ys,ITEMS):
    c="#b4442e" if not (lo<=0<=hi) else "#2f7d4f"
    b.plot([lo,hi],[yy,yy],c=c,lw=2.4,solid_capstyle="round")
    b.plot(r,yy,"o",c=c,ms=8.5,zorder=3); b.text(.08,yy+.24,f"n={n}",fontsize=8.4,color="#555")
b.axvline(0,c="#888",lw=1,ls="--")
b.set_yticks(ys); b.set_yticklabels([i[0] for i in ITEMS],fontsize=9.4)
b.set_xlim(-1.02,.62); b.set_ylim(-0.6,len(ITEMS)-0.25)
b.set_xlabel("Spearman $\\rho$: testing effort vs criterion\n(red = CI excludes zero)")
b.set_title("A  Screening does not remove the artefact",loc="left")
c=ax[1]
c.scatter(m.n_records,m.n_sp_C7,s=42,c="#2c5f8a",alpha=.8,edgecolor="w",lw=.5)
r,lo,hi,n=rho(x,m.n_sp_C7)
sl,ic,_,_,_=stats.linregress(x,np.log10(m.n_sp_C7.clip(lower=.5)))
gx=np.linspace(x.min(),x.max(),50); c.plot(10**gx,10**(sl*gx+ic),c="#2c5f8a",lw=1.8)
c.set_xscale("log"); c.set_yscale("log")
c.set_xlabel("ECOTOX records for the substance"); c.set_ylabel("species surviving the screen")
c.set_title(rf"B  More study, longer left tail ($\rho$ = {r:+.2f}, CI [{lo:+.2f}, {hi:+.2f}])",loc="left",fontsize=10.6)
plt.tight_layout(); plt.savefig("fig_03_effort.png",dpi=300,bbox_inches="tight"); plt.close()

# ---- Fig 5: matched pair
mm=pd.read_pickle("/home/claude/wq/wqp/m_ladder.pkl")
mm["s"]=mm.cas.map({"335-67-1":"PFOA","375-85-9":"PFHpA"})
fig,d=plt.subplots(figsize=(7.4,5.2))
for i,(s,col) in enumerate([("PFOA","#b4442e"),("PFHpA","#2c5f8a")]):
    g=mm[mm.s==s]
    d.scatter(np.random.default_rng(7+i).normal(i,.08,len(g)),g.ngL,s=7,alpha=.22,c=col)
    p6=(10**g.C6_min_noec.iloc[0])*1e6/g.AF_C6.iloc[0]
    d.hlines(p6,i-.32,i+.32,color=col,lw=2.8,zorder=4)
    d.text(i+.36,p6,f"criterion {p6:,.3g} ng L$^{{-1}}$",fontsize=9.4,color=col,va="center")
    d.text(i,1.2e5,f"{s}\n{int(g.n_records.iloc[0]):,} toxicity records\nmedian {g.ngL.median():.1f} ng L$^{{-1}}$\n"
           f"{100*(g.ngL>p6).mean():.2f}% exceed",ha="center",fontsize=9.6,color=col,weight="bold")
d.set_yscale("log"); d.set_xticks([0,1]); d.set_xticklabels(["PFOA","PFHpA"])
d.set_xlim(-.55,2.1); d.set_ylim(1e-2,4e6); d.set_ylabel("measured concentration (ng L$^{-1}$)")
plt.tight_layout(); plt.savefig("fig_05_pair.png",dpi=300,bbox_inches="tight"); plt.close()
print("3 figures written")
