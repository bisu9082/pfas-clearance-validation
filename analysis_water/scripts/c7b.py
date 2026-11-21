import pandas as pd, numpy as np
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
meas["test_cas"]=meas.test_cas.astype(str)
rec=pd.read_csv("records_with_species.csv.gz"); rec["test_cas"]=rec.test_cas.astype(str)
R=pd.read_csv("results_enriched.csv.gz",dtype=str)

# unit -> multiplier (to mg/L), learned from the deposited records
um=rec.dropna(subset=["mult"]).groupby("conc1_unit").mult.agg(lambda s:s.iloc[0]).to_dict()
chk=rec.groupby("conc1_unit").mult.nunique()
assert (chk<=1).all(), "unit->mult not unique"
print("units:",um)

R["conc"]=pd.to_numeric(R.conc1_mean,errors="coerce")
R["unit"]=R.conc1_unit.astype(str).str.strip()
R["mult"]=R.unit.map(um)
R=R[R.conc.notna() & R["mult"].notna() & (R.conc>0)].copy()
R["logTox"]=np.log10(R.conc*R["mult"])
R["test_id"]=R.test_id.astype(str)

sp=rec.drop_duplicates("test_id")[["test_id","test_cas","latin_name","ecotox_group","organism_habitat","media_type"]]
sp["test_id"]=sp.test_id.astype(str)
R=R.merge(sp,on="test_id",how="inner").merge(meas[["test_cas","parent"]],on="test_cas",how="inner")
R["ep"]=R.endpoint.astype(str).str.strip().str.upper().str.replace("*","",regex=False)
R["censored"]=R.conc1_mean_op.astype(str).str.strip().isin([">","<",">=","<=","~"])
print("censored rows dropped:",int(R.censored.sum()))
R=R[~R.censored].copy()
R=R.drop_duplicates(subset=["test_id","ep","conc1_mean","conc1_unit","effect","measurement"])
print("\nresult-level rows with parent+species+logTox: %d  (%d parents)"%(len(R),R.parent.nunique()))
print("effect assigned: %.1f%%"%(100*R.effect.notna().mean()))

# ---- GATE: the deposited minimum must be reproducible from this result-level table
dep=pd.read_csv("/home/claude/pfasv32/canonical_results/ecotox_pfas_measured_recurated.csv")
mn=R.drop_duplicates(subset=["test_id","ep","conc1_mean","conc1_unit"]).groupby("parent").logTox.min()
cmp=dep.set_index("parent").logTox.to_dict()
agree=sum(1 for p,v in mn.items() if p in cmp and abs(v-cmp[p])<1e-9)
print("GATE: parents whose deposited logTox min is reproduced: %d / %d"%(agree,len(cmp)))

APICAL={"MOR","GRO","REP","DVP","POP"}
CHRON={"NOEC","EC10","MATC","NOEL"}
R["apical"]=R.effect.isin(APICAL); R["chron"]=R.ep.isin(CHRON)
print("\napical %.1f%%   chronic-endpoint %.1f%%"%(100*R.apical.mean(),100*R.chron.mean()))
lo=R[R.logTox<-3]; hi=R[R.logTox>=-3]
print("\nextreme tail (<1 ug/L): n=%d, apical %.1f%%"%(len(lo),100*lo.apical.mean()))
print("rest:                   n=%d, apical %.1f%%"%(len(hi),100*hi.apical.mean()))
print("\neffect mix in the extreme tail:"); print(lo.effect.value_counts().head(10).to_string())
R.to_csv("results_parent.csv.gz",index=False,compression="gzip")
