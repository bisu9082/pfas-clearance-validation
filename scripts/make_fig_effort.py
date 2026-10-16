"""Key figure: the protective minimum is confounded with testing effort, and the
model transfers exactly on the endpoint where that confound is absent."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, sys
from scipy.stats import spearmanr
MEAS, FIXRES, SUBS, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
TEAL, AMBER, BLUE, INK, MUTED, GRID = "#0D9488","#B45309","#1D4ED8","#1A1A1A","#5A5A5A","#DDDDDD"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7.2,"text.color":INK,
                     "axes.edgecolor":"#9A9A9A","axes.linewidth":0.7,
                     "xtick.color":MUTED,"ytick.color":MUTED,"axes.labelcolor":INK})

d=pd.read_csv(MEAS); fx=pd.read_csv(FIXRES); sb=pd.read_csv(SUBS)
fig,axes=plt.subplots(1,3,figsize=(7.1,2.55),dpi=400); fig.patch.set_facecolor("#FFFFFF")
def tidy(ax):
    ax.spines[["top","right"]].set_visible(False); ax.set_axisbelow(True)
    ax.grid(color=GRID,lw=0.6); ax.tick_params(length=2.2,pad=1.8)

# (a) effort vs protective minimum
ax=axes[0]; x=np.log10(d.n_records); y=d.logTox
ax.scatter(x,y,s=11,color=AMBER,alpha=0.72,lw=0,zorder=3)
b=np.polyfit(x,y,1); xs=np.linspace(x.min(),x.max(),50)
ax.plot(xs,np.polyval(b,xs),color=INK,lw=1.2,zorder=4)
r=spearmanr(x,y).statistic
ax.text(0.04,0.06,f"$\\rho$ = {r:+.2f}",transform=ax.transAxes,fontsize=7.4,fontweight="bold")
tidy(ax); ax.set_xlabel("records per substance (log$_{10}$)",fontsize=6.9)
ax.set_ylabel("protective minimum (log$_{10}$ mg/L)",fontsize=6.9)
ax.set_title("(a)  The minimum tracks testing effort",fontsize=7.2,fontweight="bold",loc="left",pad=5)

# (b) subsample depth
ax=axes[1]
ax.axhline(r,color=AMBER,lw=1.1,ls=(0,(3,2)),zorder=2)
ax.text(sb.k.max(),r+0.025,"all records",ha="right",va="bottom",fontsize=6.2,color=AMBER)
ax.errorbar(sb.k,sb.rho_mean,yerr=[sb.rho_mean-sb.lo,sb.hi-sb.rho_mean],
            fmt="o-",color=TEAL,ms=4.4,lw=1.3,capsize=2.4,zorder=4)
ax.axhline(0,color=INK,lw=0.8,zorder=3)
tidy(ax); ax.set_xlabel("records retained per substance ($k$)",fontsize=6.9)
ax.set_ylabel("$\\rho$ (minimum, log records)",fontsize=6.9)
ax.set_xticks(sb.k); ax.set_ylim(-0.85,0.12)
ax.set_title("(b)  Half of it is sampling depth",fontsize=7.2,fontweight="bold",loc="left",pad=5)

# (c) transfer by endpoint
ax=axes[2]
lab={"minimum":"minimum\n(published)","median":"median","acute_median":"acute\nmedian","depth_min":"depth-equal\nminimum"}
o=["minimum","median","depth_min","acute_median"]
f=fx.set_index("endpoint").loc[o]
xp=np.arange(len(o))
cols=[TEAL if l>0 else AMBER for l in f.lo]
for i,c in enumerate(cols):
    ax.errorbar([xp[i]],[f.rho.iloc[i]],
                yerr=[[f.rho.iloc[i]-f.lo.iloc[i]],[f.hi.iloc[i]-f.rho.iloc[i]]],
                fmt="none",ecolor=c,lw=1.5,capsize=2.6,zorder=3)
ax.scatter(xp,f.rho,s=34,c=cols,zorder=4,lw=0.8,edgecolors="#FFFFFF")
ax.axhline(0,color=INK,lw=0.8,zorder=2)
for i,(_,r_) in enumerate(f.iterrows()):
    ax.text(i,f.hi.iloc[i]+0.045,f"$R^2$={r_.r2:+.2f}",ha="center",fontsize=6.0,color=MUTED)
tidy(ax); ax.set_xticks(xp); ax.set_xticklabels([lab[c] for c in o],fontsize=6.0)
ax.set_ylabel("external Spearman $\\rho$",fontsize=6.9); ax.set_ylim(-0.42,0.82)
ax.set_title("(c)  Transfer returns without the confound",fontsize=7.2,fontweight="bold",loc="left",pad=5)

fig.tight_layout(pad=0.5,w_pad=1.6)
fig.savefig(OUT,bbox_inches="tight",facecolor="#FFFFFF",dpi=400)
print("wrote",OUT)
