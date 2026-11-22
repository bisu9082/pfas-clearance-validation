import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from scipy import stats
rng=np.random.default_rng(20261122)
out=pd.read_csv("/home/claude/wq/v2/ladder_v2.csv")
L=pd.read_csv("/home/claude/wq/wqp/ladder.csv")
out=out.merge(L[["parent","C1_min_all","C3_min_acute","C5_spmed_acute"]],on="parent",how="left")

def rho(x,y,seed=1,B=4000):
    m=np.isfinite(x)&np.isfinite(y); x,y=np.asarray(x)[m],np.asarray(y)[m]
    r=stats.spearmanr(x,y)[0]; g=np.random.default_rng(seed); acc=[]
    for _ in range(B):
        j=g.integers(0,len(x),len(x))
        if len(set(x[j]))>2: acc.append(stats.spearmanr(x[j],y[j])[0])
    return r,np.percentile(acc,2.5),np.percentile(acc,97.5),len(x)

ln=np.log10(out.n_records)
ITEMS=[("C1  lowest, all endpoints",      rho(ln,out.C1_min_all,1)),
       ("C3  lowest acute",               rho(ln,out.C3_min_acute,2)),
       ("C6  lowest chronic NOEC + AF",   rho(ln,out.C6_pub,3)),
       ("C6$'$  + apical + duration",     rho(ln,out.C6_prime,4)),
       ("C7a  lowest apical SMCV",        rho(ln,out.minSMCV_d,5)),
       ("C7b  SSD fifth percentile",      rho(ln,out.HC5_d,6)),
       ("C5  median acute (central)",     rho(ln,out.C5_spmed_acute,7))]

fig,(a,b)=plt.subplots(1,2,figsize=(11.0,4.4),gridspec_kw={"width_ratios":[1.25,1]})
y=np.arange(len(ITEMS))[::-1]
for yi,(lab,(r,lo,hi,n)) in zip(y,ITEMS):
    col="#b4442e" if hi<0 else "#8a8a8a"
    a.plot([lo,hi],[yi,yi],color=col,lw=2.4,solid_capstyle="round")
    a.plot(r,yi,"o",color=col,ms=7,zorder=3)
    a.text(-0.99,yi+0.26,f"n={n}",fontsize=7.4,color="#777",va="center")
a.axvline(0,color="#444",lw=.8)
a.set_yticks(y); a.set_yticklabels([i[0] for i in ITEMS],fontsize=9)
a.set_xlim(-1.02,0.22); a.set_xlabel("Spearman $\\rho$: testing effort vs derived criterion",fontsize=9.5)
a.set_title("A  Every convention anchored on an extreme moves with effort",loc="left",fontsize=9.8)
a.grid(axis="x",alpha=.25,lw=.5)
for s in ("top","right"): a.spines[s].set_visible(False)

# B: the sigma channel, which is the operative one for the SSD
d=out.dropna(subset=["HC5_d"])
x=np.log10(d.n_records)
b.scatter(x,d.sd_d,s=46,color="#2c5f8a",edgecolor="w",lw=.8,zorder=3,label="fitted $\\hat\\sigma$")
k=np.polyfit(x,d.sd_d,1); xs=np.linspace(x.min(),x.max(),50)
b.plot(xs,np.polyval(k,xs),color="#2c5f8a",lw=1.4,ls="--")
r1=stats.spearmanr(x,d.sd_d)[0]; r2=stats.spearmanr(x,d.mu_d)[0]
b2=b.twinx()
b2.scatter(x,d.mu_d,s=42,marker="s",color="#c08a3e",edgecolor="w",lw=.8,zorder=3)
k2=np.polyfit(x,d.mu_d,1); b2.plot(xs,np.polyval(k2,xs),color="#c08a3e",lw=1.4,ls=":")
b.set_xlabel("$\\log_{10}$ toxicity records",fontsize=9.5)
b.set_ylabel("fitted $\\hat\\sigma$  (spread)",color="#2c5f8a",fontsize=9.5)
b2.set_ylabel("fitted $\\hat\\mu$  (centre)",color="#c08a3e",fontsize=9.5)
b.tick_params(axis="y",colors="#2c5f8a"); b2.tick_params(axis="y",colors="#c08a3e")
b.set_title("B  Mostly through its spread (63 % of slope), not its centre (37 %)\n"
            f"$\\hat\\sigma$: $\\rho$={r1:+.2f}      $\\hat\\mu$: $\\rho$={r2:+.2f}",
            loc="left",fontsize=9.8)
b.grid(alpha=.22,lw=.5)
plt.tight_layout(); plt.savefig("fig_03_effort.png",dpi=260,bbox_inches="tight"); plt.close()
print("fig_03_effort rebuilt")
for lab,(r,lo,hi,n) in ITEMS: print(f"  {lab:34s} rho={r:+.3f} [{lo:+.2f},{hi:+.2f}] n={n}")

# ================= Fig 5: matched pair, criterion-ratio framing =================
w=pd.read_pickle("/home/claude/wq/wqp/m_ladder.pkl")
PFOA="O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
PFHpA="O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
sh=set(w[w.parent==PFOA].MonitoringLocationIdentifier)&set(w[w.parent==PFHpA].MonitoringLocationIdentifier)
o=out.set_index("parent")
AF={p:{3:10.,2:50.,1:100.}[o.loc[p,"chronic_tl"]] for p in (PFOA,PFHpA)}
cr={p:10**o.loc[p,"C6_prime"]*1e6/AF[p] for p in (PFOA,PFHpA)}
cr_un={p:10**o.loc[p,"C6_pub"]*1e6/AF[p] for p in (PFOA,PFHpA)}
fig,ax=plt.subplots(figsize=(7.2,4.6))
cols={"PFOA":"#b4442e","PFHpA":"#2c5f8a"}
for i,(p,nm) in enumerate([(PFOA,"PFOA"),(PFHpA,"PFHpA")]):
    v=w[(w.parent==p)&w.MonitoringLocationIdentifier.isin(sh)].ngL.values
    xj=i+rng.normal(0,.055,len(v))
    ax.scatter(xj,v,s=7,alpha=.22,color=cols[nm],edgecolor="none",zorder=2)
    ax.plot([i-.3,i+.3],[np.median(v)]*2,color=cols[nm],lw=2.6,zorder=4)
    ax.text(i+.33,np.median(v),f"median\n{np.median(v):.1f}",fontsize=8,color=cols[nm],va="center")
    ax.plot([i-.34,i+.34],[cr[p]]*2,color=cols[nm],lw=2.2,ls="--",zorder=5)
    ax.text(i,cr[p]*1.45,f"C6$'$ criterion {cr[p]:,.0f} ng/L",fontsize=8.6,color=cols[nm],
            ha="center",fontweight="bold")
    ax.plot([i-.26,i+.26],[cr_un[p]]*2,color=cols[nm],lw=1.0,ls=":",zorder=5,alpha=.75)
    ax.text(i,cr_un[p]*1.35,f"unscreened {cr_un[p]:,.3g}",fontsize=7.2,color=cols[nm],
            ha="center",alpha=.85)
ax.annotate("",xy=(0.42,cr[PFOA]),xytext=(0.42,cr[PFHpA]),
            arrowprops=dict(arrowstyle="<->",color="#444",lw=1.5))
ax.text(0.50,np.sqrt(cr[PFOA]*cr[PFHpA]),f"{cr[PFHpA]/cr[PFOA]:,.0f}$\\times$\nfrom 313 vs 6\napical chronic\nrecords",
        fontsize=9,ha="left",va="center",fontweight="bold",color="#333")
ax.set_yscale("log"); ax.set_xlim(-.55,1.75); ax.set_xticks([0,1]); ax.set_xticklabels(["PFOA","PFHpA"],fontsize=11)
ax.set_ylabel("concentration / criterion (ng L$^{-1}$)",fontsize=10)
ax.set_title(f"Two adjacent homologues at the same {len(sh)} stations",fontsize=11,loc="left")
ax.grid(axis="y",alpha=.25,lw=.5)
for s_ in ("top","right"): ax.spines[s_].set_visible(False)
ax.set_ylim(0.04,2.2e5)
plt.tight_layout(); plt.savefig("fig_05_pair.png",dpi=260,bbox_inches="tight"); plt.close()
print(f"fig_05_pair rebuilt: PFOA {cr[PFOA]:.1f} vs PFHpA {cr[PFHpA]:,.0f} ng/L = {cr[PFHpA]/cr[PFOA]:,.0f}x; stations {len(sh)}")
