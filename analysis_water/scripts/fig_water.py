import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
m=pd.read_pickle("m_ladder.pkl"); ex=pd.read_csv("exceedance.csv")
NAME={"375-22-4":"PFBA","335-67-1":"PFOA","1763-23-1":"PFOS","307-24-4":"PFHxA",
"375-73-5":"PFBS","2706-90-3":"PFPeA","375-85-9":"PFHpA","355-46-4":"PFHxS",
"375-95-1":"PFNA","754-91-6":"FOSA","335-76-2":"PFDA","27619-97-2":"6:2 FTS",
"307-55-1":"PFDoA","2058-94-8":"PFUnA","72629-94-8":"PFTrDA","30334-69-1":"FBSA",
"41997-13-1":"FHxSA","13252-13-6":"HFPO-DA","4151-50-2":"N-EtFOSA"}
m["subst"]=m.cas.map(NAME).fillna(m.cas)
m["P6"]=(10.0**m.C6_min_noec)*1000.0/m.AF_C6
fig,ax=plt.subplots(1,3,figsize=(17.2,5.3))

# A: ladder
a=ax[0]; lbl=["lowest, all\nendpoints\n(as published)","lowest,\nLOECs\nremoved","lowest\nacute","median acute,\nrecord","median acute,\nspecies","lowest chronic\nNOEC,\nEU-TGD AF"]
y1=ex.pct_samples_RQgt1.values; y2=ex.pct_sites_RQgt1.values
x=np.arange(6); w=.38
a.bar(x-w/2,y1,w,color="#2c5f8a",label="samples")
a.bar(x+w/2,y2,w,color="#9ec5e0",edgecolor="#2c5f8a",label="monitoring sites")
for i,(u,v) in enumerate(zip(y1,y2)):
    a.text(i-w/2,u+1.6,f"{u:.1f}",ha="center",fontsize=8)
    a.text(i+w/2,v+1.6,f"{v:.1f}",ha="center",fontsize=8)
a.set_xticks(x); a.set_xticklabels(lbl,fontsize=7.4)
a.set_ylabel("% exceeding the derived criterion (RQ > 1)")
a.set_ylim(0,112); a.legend(fontsize=8,frameon=False,loc="upper center")
a.set_title("A  Same data, same water, six curation conventions",fontsize=9.6,loc="left")
a.text(.50,.60,"21,346 samples · 1,180 US sites · 2016–2026",transform=a.transAxes,
       ha="center",fontsize=7.8,color="#444")

# B: PNEC_C6 vs effort
b=ax[1]
d=m.groupby("subst").agg(P6=("P6","first"),recs=("n_records","first")).dropna()
b.scatter(d.recs,d.P6*1000,s=46,c="#b4442e",zorder=3,edgecolor="w",lw=.6)
OFF={"PFPeA":(-30,7),"PFHpA":(6,-10),"PFTrDA":(-34,-9),"PFUnA":(-6,9),"PFDoA":(6,6),
     "PFOA":(-30,-11),"PFOS":(5,4),"PFBA":(-28,-10),"HFPO-DA":(5,-10),"PFNA":(5,-9),
     "6:2 FTS":(7,-2),"FOSA":(-32,4),"FBSA":(5,5),"N-EtFOSA":(6,2),"PFDA":(5,5),
     "PFHxS":(5,4),"PFHxA":(5,-9),"PFBS":(5,5)}
for s,r in d.iterrows():
    b.annotate(s,(r.recs,r.P6*1000),fontsize=7,xytext=OFF.get(s,(4,3)),textcoords="offset points")
b.set_xlim(5.5,3200)
xs=np.log10(d.recs); ys=np.log10(d.P6*1000)
sl,ic,_,pv,_=stats.linregress(xs,ys); rho,prho=stats.spearmanr(d.recs,d.P6)
gx=np.linspace(xs.min(),xs.max(),50)
b.plot(10**gx,10**(sl*gx+ic),c="#b4442e",lw=1.8,zorder=2,ls="-" if prho<.05 else (0,(5,3)))
b.axhline(.65,c="#2c5f8a",ls="--",lw=1.3); b.text(6.6,.78,"EU EQS, PFOS (0.65 ng/L)",fontsize=7.4,color="#2c5f8a")
b.set_xscale("log"); b.set_yscale("log")
b.set_xlabel("number of ECOTOX toxicity records for the substance")
b.set_ylabel("PNEC, proper EU-TGD chronic route (ng/L)")
b.set_title(f"B  Doing it correctly does not remove the artefact\n"
            f"Spearman $\\rho$ = {rho:+.2f} (p = {prho:.1e})",fontsize=9.6,loc="left")

# C: the matched pair
c=ax[2]
pair=[("PFOA","#b4442e"),("PFHpA","#2c5f8a")]
for i,(s,col) in enumerate(pair):
    d2=m[m.subst==s]
    c.scatter(np.random.default_rng(7+i).normal(i,.075,len(d2)),d2.ngL,s=6,alpha=.22,c=col,zorder=2)
    p6=d2.P6.iloc[0]*1000
    c.hlines(p6,i-.3,i+.3,color=col,lw=2.6,zorder=4)
    c.text(i+.33,p6,f"PNEC {p6:,.3g} ng/L",fontsize=8,color=col,va="center")
    c.text(i,2.3e4,f"{s}\n{int(d2.n_records.iloc[0]):,} tox records\n{100*(d2.ugL>d2.P6.iloc[0]).mean():.1f}% exceed",
           ha="center",fontsize=8.4,color=col,weight="bold")
c.set_yscale("log"); c.set_xticks([0,1]); c.set_xticklabels(["PFOA","PFHpA"])
c.set_xlim(-.55,2.05); c.set_ylim(1e-2,2e5)
c.set_ylabel("measured concentration (ng/L)")
c.set_title("C  Two PFAS, the same rivers, opposite verdicts",fontsize=9.6,loc="left")
plt.tight_layout(); plt.savefig("fig_water.png",dpi=165)
print("medians: PFOA %.2f ng/L, PFHpA %.2f ng/L"%(m[m.subst=="PFOA"].ngL.median(),m[m.subst=="PFHpA"].ngL.median()))
print("shared sites:",len(set(m[m.subst=='PFOA'].MonitoringLocationIdentifier)&set(m[m.subst=='PFHpA'].MonitoringLocationIdentifier)))
print("saved fig_water.png")
