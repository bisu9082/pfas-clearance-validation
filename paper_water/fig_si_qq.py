import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
RDLogger.DisableLog('rdApp.*')
plt.rcParams.update({"font.size":9,"axes.titlesize":9,"axes.labelsize":8.6,
 "xtick.labelsize":7.6,"ytick.labelsize":7.6})
R=pd.read_csv("/home/claude/wq/results_parent.csv.gz")
AP={"MOR","GRO","REP","DVP","POP"}; EP={"NOEC","NOEL","MATC","EC10"}
sel=R[R.ep.isin(EP)&R.effect.isin(AP)&(R.conc1_type=="A")]
smcv=sel.groupby(["parent","latin_name"]).logTox.mean().reset_index()
g=[(p,d.logTox.values) for p,d in smcv.groupby("parent") if len(d)>=5]
g.sort(key=lambda t:-len(t[1]))
NAME={"O=S(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFOS",
"O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFOA",
"O=C(O)C(F)(OC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F":"HFPO-DA",
"O=S(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFBS",
"O=C(O)C(F)(F)C(F)(F)C(F)(F)F":"PFBA",
"O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFNA",
"O=S(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFHxS",
"O=S(=O)(O)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"6:2 FTS"}
def lbl(p):
    if p in NAME: return NAME[p]
    m=Chem.MolFromSmiles(p)
    return rdMolDescriptors.CalcMolFormula(m) if m else p[:16]
n=len(g); ncol=4; nrow=int(np.ceil(n/ncol))
fig,ax=plt.subplots(nrow,ncol,figsize=(10.6,2.5*nrow))
ax=np.atleast_1d(ax).ravel()
for i,(p,v) in enumerate(g):
    a=ax[i]; v=np.sort(v); k=len(v)
    (osm,osr),(sl,ic,r)=stats.probplot(v,dist="norm")
    a.plot(osm,osr,"o",ms=4,c="#2c5f8a")
    a.plot(osm,sl*osm+ic,"-",c="#b4442e",lw=1.3)
    ad=stats.anderson(v,"norm")
    rej=ad.statistic>ad.critical_values[2]
    a.set_title(f"{lbl(p)}  n={k}\nAD={ad.statistic:.2f} {'(reject)' if rej else '(ok)'}  $R^2$={r**2:.3f}",
                fontsize=8.2,color="#b4442e" if rej else "#333")
    a.set_xlabel("theoretical quantile"); a.set_ylabel("log$_{10}$ SMCV (mg L$^{-1}$)")
for j in range(n,len(ax)): ax[j].axis("off")
plt.tight_layout(); plt.savefig("fig_si_qq.png",dpi=250,bbox_inches="tight")
print("panels:",n)
