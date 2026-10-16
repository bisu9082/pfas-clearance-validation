"""Per-OECD-class external rank correlation with bootstrap CIs - the figure the
reviewers asked for: uncertainty shown, not point estimates alone."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')
BROAD, OUT = sys.argv[1], sys.argv[2]
SEED, MIN_N = 20260620, 10
TEAL, AMBER, INK, MUTED, GRID = "#0D9488", "#B45309", "#1A1A1A", "#5A5A5A", "#DDDDDD"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7.2,"text.color":INK,
                     "axes.edgecolor":"#9A9A9A","axes.linewidth":0.7,
                     "xtick.color":MUTED,"ytick.color":MUTED})

def feats(sm):
    X,k=[],[]
    for i,s in enumerate(sm):
        m=Chem.MolFromSmiles(str(s))
        if m is not None:
            X.append(np.array(AllChem.GetMorganFingerprintAsBitVect(m,2,2048),dtype=np.int8)); k.append(i)
    return np.vstack(X), np.array(k)

def rho_ci(y,p,B=4000,rng=None):
    rng=rng or np.random.default_rng(SEED); n=len(y); o=[]
    for _ in range(B):
        i=rng.integers(0,n,n)
        if len(np.unique(y[i]))<3: continue
        o.append(spearmanr(y[i],p[i]).statistic)
    return float(spearmanr(y,p).statistic), float(np.nanpercentile(o,2.5)), float(np.nanpercentile(o,97.5))

broad=pd.read_csv(BROAD); Xb,kb=feats(broad.canonical_smiles); yb=broad.logTox.to_numpy()[kb]
full=RandomForestRegressor(n_estimators=300,max_features="sqrt",random_state=SEED,n_jobs=-1).fit(Xb,yb)

d=pd.read_csv("/home/claude/oecd/ext_R1_oecd.csv"); d=d[d.oecd_pfas].reset_index(drop=True)
Xe,ke=feats(d.canonical_smiles); d=d.iloc[ke].reset_index(drop=True)
ye=d.logTox.to_numpy(); pr=full.predict(Xe)

rows=[]
r,l,h=rho_ci(ye,pr); rows.append(("All OECD PFAS", len(ye), r, l, h, True))
for c,g in d.groupby("oecd_class"):
    if len(g)<MIN_N: continue
    i=g.index.to_numpy(); r,l,h=rho_ci(ye[i],pr[i])
    rows.append((c, len(i), r, l, h, False))
rows = [rows[0]] + sorted(rows[1:], key=lambda x: x[2])

fig,ax=plt.subplots(figsize=(5.6,2.95),dpi=400); fig.patch.set_facecolor("#FFFFFF")
y=np.arange(len(rows))[::-1]
for (lab,n,r,l,h,is_all),yy in zip(rows,y):
    col = INK if is_all else (TEAL if l>0 or h<0 else AMBER)
    ax.plot([l,h],[yy,yy],color=col,lw=1.6,solid_capstyle="round",zorder=3)
    ax.plot([r],[yy],"o",ms=4.6,color=col,mec="#FFFFFF",mew=0.8,zorder=4)
ax.axvline(0,color=INK,lw=0.9,zorder=2)
ax.set_yticks(y); ax.set_yticklabels([f"{l}  ($n$={n})" for l,n,*_ in rows],fontsize=6.6)
ax.set_xlabel("external Spearman $\\rho$ on measured PFAS  (95% bootstrap CI)",fontsize=7.0,labelpad=9)
ax.set_xlim(-1.05,1.05); ax.grid(axis="x",color=GRID,lw=0.6); ax.set_axisbelow(True)
ax.spines[["top","right","left"]].set_visible(False); ax.tick_params(axis="y",length=0)
ax.axhline(len(rows)-1.5,color="#BBBBBB",lw=0.7,ls=(0,(3,2)))
ax.text(1.02,len(rows)-1.0,"",ha="right")
ax.set_title("Transfer to PFAS, resolved by OECD (2021) structure class",
             fontsize=7.6,fontweight="bold",loc="left",pad=7)
ax.text(0,-1.45,"teal = CI excludes zero   ·   amber = CI spans zero",
        transform=ax.get_yaxis_transform(), fontsize=6.3, color=MUTED, va="top")
fig.tight_layout(pad=0.45); fig.savefig(OUT,bbox_inches="tight",facecolor="#FFFFFF",dpi=400)
print("wrote",OUT)
for lab,n,r,l,h,_ in rows: print(f"  {lab:38s} n={n:3d}  rho={r:+.3f} [{l:+.2f},{h:+.2f}]")
