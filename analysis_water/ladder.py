import pandas as pd, numpy as np, json
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from scipy import stats
RDLogger.DisableLog('rdApp.*'); UNC=rdMolStandardize.Uncharger()
def canon(s):
    m=Chem.MolFromSmiles(str(s))
    if m is None: return None
    fr=Chem.GetMolFrags(m,asMols=True,sanitizeFrags=False)
    fl=[f for f in fr if any(a.GetSymbol()=="F" for a in f.GetAtoms())]
    big=max(fl if fl else fr,key=lambda f:f.GetNumHeavyAtoms())
    try: big=UNC.uncharge(big); Chem.SanitizeMol(big)
    except Exception: pass
    return Chem.MolToSmiles(big)
meas=pd.read_csv("/home/claude/pfasv32/archive_v2/ecotox_pfas_measured.csv")
meas["parent"]=[canon(s) for s in meas.canonical_smiles]; meas=meas[meas.parent.notna()]
rec=pd.read_csv("/home/claude/wq/records_with_species.csv.gz")
rec["test_cas"]=rec.test_cas.astype(str); meas["test_cas"]=meas.test_cas.astype(str)
rec=rec.drop_duplicates(subset=["test_id","ep","conc1_mean","conc1_unit"])
rec=rec.merge(meas[["test_cas","parent"]],on="test_cas",how="inner")
ACUTE={"LC50","EC50","IC50"}
rec["is_acute"]=rec.ep.isin(ACUTE)
print("endpoint mix (all %d records):"%len(rec)); print(rec.ep.value_counts().to_string())
print("\nLOEC share: %.1f%%  of non-acute: %.1f%%"%(
      100*(rec.ep=="LOEC").mean(), 100*(rec.loc[~rec.is_acute,"ep"]=="LOEC").mean()))

def trophic(g):
    g=str(g).split(";")[0].strip()
    if g=="Fish": return "fish"
    if g in ("Algae","Moss, Hornworts"): return "algae"
    if g in ("Crustaceans","Insects/Spiders","Molluscs","Worms","Other Invertebrates","Invertebrates"): return "invertebrate"
    return "other"
rec["troph"]=rec.ecotox_group.map(trophic)

rows=[]
for p,g in rec.groupby("parent"):
    ga=g[g.is_acute]; gc=g[~g.is_acute]
    noec=gc[gc.ep=="NOEC"]
    ctl=len(set(noec.troph)-{"other"})
    rows.append(dict(parent=p,n_records=len(g),n_acute=len(ga),n_noec=len(noec),
      C1_min_all=g.logTox.min(),
      C2_min_noLOEC=g.loc[~g.ep.isin(["LOEC","LOEL"]),"logTox"].min(),
      C3_min_acute=ga.logTox.min() if len(ga) else np.nan,
      C4_med_acute=ga.logTox.median() if len(ga) else np.nan,
      C5_spmed_acute=(ga.groupby("latin_name").logTox.median().median() if len(ga) else np.nan),
      C6_min_noec=noec.logTox.min() if len(noec) else np.nan,
      chronic_tl=ctl,n_noec_sp=noec.latin_name.nunique()))
D=pd.DataFrame(rows)
# EU TGD AF for the chronic route
def af_chronic(r):
    if r.chronic_tl>=3: return 10.0
    if r.chronic_tl==2: return 50.0
    if r.chronic_tl==1: return 100.0
    return np.nan
D["AF_C6"]=D.apply(af_chronic,axis=1)
A=D[D.n_acute>0].copy()
print("\nparents with acute: %d   with >=1 NOEC: %d"%(len(A),A.n_noec.gt(0).sum()))
print("\n=== convention ladder: gap to the species-resolved acute median (log10) ===")
for c in ["C1_min_all","C2_min_noLOEC","C3_min_acute","C4_med_acute","C6_min_noec"]:
    d=(A.C5_spmed_acute-A[c]).dropna()
    x=np.log10(A.loc[d.index,"n_records"]); y=A.loc[d.index,c]
    r,pv=stats.spearmanr(x,y)
    print(f"  {c:16s} n={len(d):3d}  median gap {d.median():+.2f}  max {d.max():+.2f}   rho(effort)={r:+.3f} p={pv:.2g}")
D.to_csv("ladder.csv",index=False); print("\nwrote ladder.csv")
