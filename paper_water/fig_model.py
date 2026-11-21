import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({"font.size":11,"axes.titlesize":11.2,"axes.labelsize":11,
 "xtick.labelsize":10,"ytick.labelsize":10,"legend.fontsize":9.8})
P=pd.read_csv("/home/claude/wq/broad/r5_pred.csv")
bins=[(1,2),(2,5),(5,20),(20,100),(100,10**9)]
lab=["1","2–4","5–19","20–99","100+"]
rows=[]
for (lo,hi),l in zip(bins,lab):
    d=P[(P.n_rec>=lo)&(P.n_rec<hi)]
    r=dict(lab=l,n=len(d))
    for nm,o,p in [("min",d.acute_min,d.pred_min),("med",d.acute_median,d.pred_med)]:
        r[nm+"_bias"]=float((p-o).mean())
        r[nm+"_r2"]=float(1-((o-p)**2).sum()/((o-o.mean())**2).sum())
        r[nm+"_rho"]=float(stats.spearmanr(o,p)[0])
    rows.append(r)
R=pd.DataFrame(rows)
fig,ax=plt.subplots(1,3,figsize=(15.0,4.7))
x=np.arange(len(R))

a=ax[0]
a.plot(x,R.min_bias,"o-",c="#b4442e",ms=8,lw=2.4,label="minimum as target")
a.plot(x,R.med_bias,"s-",c="#2f7d4f",ms=7,lw=2.4,label="median as target")
a.axhline(0,c="#888",lw=1,ls="--")
for i,v in enumerate(R.min_bias): a.annotate(f"{v:+.2f}",(i,v),xytext=(0,9),textcoords="offset points",
                                             ha="center",fontsize=9,color="#b4442e")
a.set_xticks(x); a.set_xticklabels(R.lab); a.set_xlabel("acute records for the chemical")
a.set_ylabel("prediction bias (log$_{10}$ units)")
a.set_title("A  The minimum-target model over-predicts\nwell-studied chemicals by 2.3 orders",loc="left")
a.legend(frameon=False,loc="upper left"); a.set_ylim(-1.1,3.0)

b=ax[1]
b.plot(x,R.min_r2,"o-",c="#b4442e",ms=8,lw=2.4)
b.plot(x,R.med_r2,"s-",c="#2f7d4f",ms=7,lw=2.4)
b.axhline(0,c="#888",lw=1,ls="--")
b.set_xticks(x); b.set_xticklabels(R.lab); b.set_xlabel("acute records for the chemical")
b.set_ylabel("$R^2$ within the stratum")
b.set_title("B  Aggregate $R^2$ = 0.37 hides a collapse\nto $-$0.95 where study is deepest",loc="left")
b.annotate(f"{R.min_r2.iloc[-1]:.2f}",(x[-1],R.min_r2.iloc[-1]),xytext=(-6,12),
           textcoords="offset points",ha="right",fontsize=9.6,color="#b4442e")

c=ax[2]
c.plot(x,R.min_rho,"o-",c="#b4442e",ms=8,lw=2.4)
c.plot(x,R.med_rho,"s-",c="#2f7d4f",ms=7,lw=2.4)
c.set_xticks(x); c.set_xticklabels(R.lab); c.set_xlabel("acute records for the chemical")
c.set_ylabel(r"Spearman $\rho$ within the stratum")
c.set_ylim(0,0.95)
c.set_title("C  Rank correlation is flat, so a $\\rho$-based\nvalidation does not see the failure",loc="left")
plt.tight_layout(); plt.savefig("fig_07_model.png",dpi=300,bbox_inches="tight")
print(R.to_string(index=False,float_format=lambda v:f"{v:.3f}"))
print("\nbias swing  min: %.2f log   median: %.2f log"%(R.min_bias.max()-R.min_bias.min(),
                                                        R.med_bias.max()-R.med_bias.min()))
R.to_csv("/home/claude/wq/broad/r5_strata.csv",index=False)
