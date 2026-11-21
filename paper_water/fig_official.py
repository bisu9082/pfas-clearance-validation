import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11,"axes.titlesize":11.5,"axes.labelsize":11,
 "xtick.labelsize":10,"ytick.labelsize":10})
m=pd.read_pickle("/home/claude/wq/wqp/m_ladder.pkl")
C=pd.read_csv("/home/claude/wq/c7_ssd.csv")
m=m.merge(C,on="parent",how="left")
NAME={"335-67-1":"PFOA","1763-23-1":"PFOS"}; m["s"]=m.cas.map(NAME)
CONV=[("C1_min_all","C1  lowest, all endpoints",1000.0),
      ("C3_min_acute","C3  lowest acute",1000.0),
      ("C6_min_noec","C6  lowest chronic NOEC + AF","AF_C6"),
      ("C7_min","C7a lowest apical chronic SMCV + AF",10.0),
      ("C7_HC5","C7b SSD HC5 (regulatory style)",1.0),
      ("C5_spmed_acute","C5  median acute, species",1000.0)]
OFF={"PFOS":[("US EPA 2024 final CCC",250),("EU AA-EQS, direct toxicity",230),
             ("Australia DGV 95% species",480),("Australia DGV 99% species",9.1),
             ("Canada FEQG, SSD 5th pct",6800)],
     "PFOA":[("US EPA 2024 final CCC",100000)]}
fig,ax=plt.subplots(1,2,figsize=(12.4,5.4),sharey=True)
for k,(s,A) in enumerate(zip(["PFOS","PFOA"],ax)):
    d=m[m.s==s].iloc[0]
    ov=[v for _,v in OFF[s]]
    lo,hi=min(ov),max(ov)
    A.axhspan(lo,hi,color="#2c5f8a",alpha=.13,zorder=1)
    for lbl,v in OFF[s]:
        A.axhline(v,color="#2c5f8a",lw=1.1,ls="--",zorder=2,label=f"{lbl}  ({v:,.3g})")
    xs,ys,lab=[],[],[]
    for i,(col,l,AF) in enumerate(CONV):
        a=d[AF] if isinstance(AF,str) else AF
        if pd.isna(d[col]) or pd.isna(a): continue
        xs.append(i); ys.append((10**d[col])*1e6/a); lab.append(l)
    A.scatter(xs,ys,s=95,c="#b4442e",zorder=5,edgecolor="w",lw=1.0)
    A.plot(xs,ys,c="#b4442e",lw=1.0,alpha=.45,zorder=4)
    for x,y in zip(xs,ys):
        A.annotate(f"{y:,.3g}",(x,y),xytext=(7,-3),textcoords="offset points",fontsize=8.2,color="#b4442e")
    A.set_yscale("log"); A.set_ylim(1e-4,1e7); A.set_xlim(-.5,5.5)
    A.set_xticks(range(6)); A.set_xticklabels(["C1","C3","C6","C7a","C7b","C5"])
    A.legend(fontsize=7.6,frameon=True,framealpha=.92,loc="lower right",
             title="official aquatic-life criteria (ng L$^{-1}$)",title_fontsize=7.8)
    A.set_title(f"{'AB'[k]}  {s}   ({int(d.n_records):,} curated records)",loc="left")
    A.set_xlabel("derivation convention")
    if k==0: A.set_ylabel("derived criterion (ng L$^{-1}$)")
plt.tight_layout(); plt.savefig("fig_04_official.png",dpi=300,bbox_inches="tight")
print("saved")
