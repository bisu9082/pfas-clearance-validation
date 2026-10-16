"""Final figure: endpoint repair restores predictive validity, yet clearance
still fails at an equivalently protective threshold."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, json, sys
TEAL,AMBER,BLUE,INK,MUTED,GRID="#0D9488","#B45309","#1D4ED8","#1A1A1A","#5A5A5A","#DDDDDD"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7.2,"text.color":INK,
    "axes.edgecolor":"#9A9A9A","axes.linewidth":0.7,"xtick.color":MUTED,
    "ytick.color":MUTED,"axes.labelcolor":INK})
IV,REC,MET,OUT=sys.argv[1:5]
iv=pd.read_csv(IV); d=iv[iv.in_domain]; m=json.load(open(MET))
rec=pd.read_csv(REC,dtype={"test_cas":str}); ACUTE={"LC50","EC50","IC50"}
am=rec[rec.ep.isin(ACUTE)].groupby("test_cas").logTox.median()

fig,axes=plt.subplots(1,3,figsize=(7.1,2.6),dpi=400); fig.patch.set_facecolor("#FFFFFF")
def tidy(ax):
    ax.spines[["top","right"]].set_visible(False); ax.set_axisbelow(True)
    ax.grid(color=GRID,lw=0.6); ax.tick_params(length=2.2,pad=1.8)

# (a) before/after on the two headline quantities
ax=axes[0]; x=np.arange(2); w=0.34
pub=[-0.580,10.53]; new=[m["ext_r2"],m["cond_w"]]
ax2=ax.twinx() if False else None
ax.bar(x[0]-w/2,pub[0],w,color=AMBER,zorder=3)
ax.bar(x[0]+w/2,new[0],w,color=TEAL,zorder=3)
ax.axhline(0,color=INK,lw=0.8,zorder=4)
ax.set_xticks([x[0]-w/2,x[0]+w/2])
ax.set_xticklabels(["minimum\n(published)","acute median\n(rebuilt)"],fontsize=6.3)
ax.set_ylim(-0.72,0.5); tidy(ax); ax.set_ylabel("external $R^2$",fontsize=6.9)
for xi,v in ((x[0]-w/2,pub[0]),(x[0]+w/2,new[0])):
    ax.text(xi,v+(0.04 if v>0 else -0.08),f"{v:+.2f}",ha="center",fontsize=6.6,fontweight="bold")
ax.set_title("(a)  Predictive validity restored",fontsize=7.2,fontweight="bold",loc="left",pad=5)

# (b) interval width
ax=axes[1]
ax.bar([0],[10.53],0.42,color=AMBER,zorder=3)
ax.bar([1],[m["cond_w"]],0.42,color=TEAL,zorder=3)
for xi,v in ((0,10.53),(1,m["cond_w"])):
    ax.text(xi,v+0.3,f"{v:.2f}",ha="center",fontsize=6.8,fontweight="bold")
tidy(ax); ax.set_xticks([0,1])
ax.set_xticklabels(["minimum\n(published)","acute median\n(rebuilt)"],fontsize=6.3)
ax.set_ylabel("conditional 90% width (log$_{10}$ mg/L)",fontsize=6.6); ax.set_ylim(0,12.6)
ax.set_title("(b)  Intervals halve",fontsize=7.2,fontweight="bold",loc="left",pad=5)

# (c) clearance vs threshold
ax=axes[2]
th=np.linspace(-3,3.2,160)
frac=[(d.lower>t).mean()*100 for t in th]
ax.plot(th,frac,color=TEAL,lw=1.6,zorder=4)
ax.fill_between(th,0,frac,color=TEAL,alpha=0.13,zorder=2)
for t,lab,c,dx,dy,ha,va in ((-1.0,"0.1 mg/L\npublished",AMBER,-0.20,-3.2,"right","top"),
                            (0.0,"1 mg/L\nGHS acute 1",MUTED,0.16,2.4,"left","bottom"),
                            (1.69,"49 mg/L\nequal stringency",BLUE,0.16,2.4,"left","bottom")):
    ax.axvline(t,color=c,lw=1.0,ls=(0,(3,2)),zorder=3)
    v=(d.lower>t).mean()*100
    ax.text(t-0.1,43,lab,rotation=90,ha="right",va="top",fontsize=5.7,color=c)
    ax.plot([t],[v],"o",ms=4,color=c,zorder=5)
    ax.text(t+dx,v+dy,f"{v:.1f}%",fontsize=6.3,color=c,fontweight="bold",ha=ha,va=va)
tidy(ax); ax.set_xlabel("hazard threshold (log$_{10}$ mg/L)",fontsize=6.9)
ax.set_ylabel("in-domain structures cleared (%)",fontsize=6.6)
ax.set_ylim(0,44); ax.set_xlim(-3,3.2)
ax.set_title("(c)  Clearance still fails when matched",fontsize=7.2,fontweight="bold",loc="left",pad=5)

fig.tight_layout(pad=0.5,w_pad=1.7)
fig.savefig(OUT,bbox_inches="tight",facecolor="#FFFFFF",dpi=400)
print("wrote",OUT)
