import pandas as pd, numpy as np, json, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({"font.size":11,"axes.titlesize":11.3,"axes.labelsize":11,
 "xtick.labelsize":10,"ytick.labelsize":10,"legend.fontsize":9.6})
d=pd.read_csv("/home/claude/wq/broad/broad_summary.csv")
rows=json.load(open("/home/claude/wq/broad/broad_subsample.json"))
fig,ax=plt.subplots(1,2,figsize=(12.8,5.6))

a=ax[0]
a.scatter(d.n_rec,d.acute_min,s=5,alpha=.13,c="#b4442e",lw=0,label="minimum")
a.scatter(d.n_rec,d.acute_median,s=5,alpha=.13,c="#2f7d4f",lw=0,label="median")
x=np.log10(d.n_rec)
for col,c in [("acute_min","#b4442e"),("acute_median","#2f7d4f")]:
    sl,ic,_,_,_=stats.linregress(x,d[col])
    gx=np.linspace(x.min(),x.max(),60); a.plot(10**gx,sl*gx+ic,c=c,lw=2.6,zorder=4)
rmin=stats.spearmanr(x,d.acute_min)[0]; rmed=stats.spearmanr(x,d.acute_median)[0]
a.set_xscale("log"); a.set_xlabel("acute records for the chemical")
a.set_ylabel("log$_{10}$ acute value (mg L$^{-1}$)")
a.set_title("A  4,757 ECOTOX chemicals, all classes",loc="left")
a.text(.025,.055,f"minimum  $\\rho$ = {rmin:+.3f}\nmedian   $\\rho$ = {rmed:+.3f}",
       transform=a.transAxes,fontsize=10.4,va="bottom",
       bbox=dict(fc="white",ec="#ccc",boxstyle="round,pad=0.4"))
h=[plt.Line2D([],[],marker="o",ls="",c="#b4442e",ms=7,label="minimum acute value"),
   plt.Line2D([],[],marker="o",ls="",c="#2f7d4f",ms=7,label="median acute value")]
a.legend(handles=h,loc="upper right",frameon=False)

b=ax[1]
k=[r["k"] for r in rows]; xs=np.arange(len(k)); w=.2
b.bar(xs-1.5*w,[r["min_baseline"] for r in rows],w,color="#b4442e",label="minimum, all records")
b.bar(xs-0.5*w,[r["min_depth_equal"] for r in rows],w,color="#e8a596",edgecolor="#b4442e",
      label="minimum, equal depth $k$")
b.bar(xs+0.5*w,[r["med_baseline"] for r in rows],w,color="#2f7d4f",label="median, all records")
b.bar(xs+1.5*w,[r["med_depth_equal"] for r in rows],w,color="#a9ccae",edgecolor="#2f7d4f",
      label="median, equal depth $k$")
for i,r in enumerate(rows):
    b.annotate("",xy=(i-0.5*w,r["min_depth_equal"]-0.012),xytext=(i-1.5*w,r["min_baseline"]-0.012),
               arrowprops=dict(arrowstyle="->",color="#444",lw=1.4,
                               connectionstyle="arc3,rad=-0.45"))
    b.text(i-w,0.022,f"{r['depth_effect']:+.2f}",ha="center",fontsize=9.6,color="#b4442e",weight="bold")
    b.text(i+w,0.022,f"{r['med_depth_equal']-r['med_baseline']:+.2f}",ha="center",fontsize=9.6,color="#2f7d4f",weight="bold")
b.text(-0.62,0.062,"change when depth is equalised:",fontsize=9.4,color="#333")
b.axhline(0,c="#888",lw=1)
b.set_xticks(xs); b.set_xticklabels([f"$k$={kk}\n$n$={r['n_eligible']:,}" for kk,r in zip(k,rows)],fontsize=9.4)
b.set_ylabel(r"Spearman $\rho$ with log$_{10}$(records)")
b.set_ylim(-0.52,0.095); b.set_xlim(-0.65,len(k)-0.35)
b.legend(loc="lower center",bbox_to_anchor=(0.5,-0.33),ncol=2,frameon=False,fontsize=9.4)
b.set_title("B  Equalising sampling depth removes half of the minimum's\n"
            "correlation and almost none of the median's",loc="left",fontsize=10.8)
plt.tight_layout(); plt.savefig("fig_06_broad.png",dpi=300,bbox_inches="tight")
print("rho min %.3f  median %.3f"%(rmin,rmed))
for r in rows: print(r)
