import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11,"axes.titlesize":11.5,"axes.labelsize":11,
 "xtick.labelsize":10,"ytick.labelsize":10,"legend.fontsize":9.5})
R=pd.read_csv("/home/claude/wq/results_parent.csv.gz")
DEP={"EC50","IC50","LC50","LOEC","LOEL","MATC","NOEC","NOEL"}
G=R[R.ep.isin(DEP)].copy()
d=G.drop_duplicates(subset=["test_id","ep","conc1_mean","conc1_unit"])
AP={"MOR","GRO","REP","DVP","POP"}
fig,ax=plt.subplots(1,3,figsize=(13.2,4.3))

a=ax[0]
ec=d.ep.value_counts()
cols=["#b4442e" if k=="LOEC" else ("#2c5f8a" if k in ("NOEC","NOEL","MATC") else "#9aa7b2") for k in ec.index]
a.bar(range(len(ec)),ec.values,color=cols)
a.set_xticks(range(len(ec))); a.set_xticklabels(ec.index,rotation=40,ha="right")
for i,v in enumerate(ec.values): a.text(i,v+30,f"{v:,}",ha="center",fontsize=8.6)
a.set_ylabel("curated records"); a.set_ylim(0,ec.max()*1.16)
a.set_title("A  Endpoint composition (n = 4,744)",loc="left")
a.text(.97,.72,"red: LOEC, 38.2% of records,\nan effect level used as a\nno-effect basis under C1",
       transform=a.transAxes,ha="right",fontsize=8.4,color="#b4442e")

b=ax[1]
fc=G.effect.value_counts().head(12)
cb=["#2c5f8a" if k in AP else "#c9a227" for k in fc.index]
b.barh(range(len(fc))[::-1],fc.values,color=cb)
b.set_yticks(range(len(fc))[::-1]); b.set_yticklabels(fc.index)
b.set_xlabel("result rows")
b.set_title("B  Effect types; blue = apical",loc="left")
b.text(.96,.08,f"apical share {100*G.effect.isin(AP).mean():.1f}%",transform=b.transAxes,
       ha="right",fontsize=9)

c=ax[2]
n=d.groupby("parent").size().sort_values(ascending=False).values
c.plot(np.arange(1,len(n)+1),n,"o-",ms=3.6,lw=1.1,color="#2c5f8a")
c.set_yscale("log"); c.set_xlabel("substance, ranked"); c.set_ylabel("curated records")
c.set_title("C  Testing effort is extremely uneven",loc="left")
top=n[:2].sum(); c.text(.96,.84,f"87 substances\ntop 2 hold {100*top/n.sum():.0f}% of records\nmedian {int(np.median(n))} records",
     transform=c.transAxes,ha="right",fontsize=8.8)
plt.tight_layout(); plt.savefig("fig_01_data.png",dpi=300,bbox_inches="tight")
print("endpoint mix:",ec.to_dict())
print("top2 share %.1f%%  median %d  max %d"%(100*top/n.sum(),np.median(n),n.max()))
