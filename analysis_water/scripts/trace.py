import pandas as pd, numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
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
meas["parent"]=[canon(s) for s in meas.canonical_smiles]
rec=pd.read_csv("/home/claude/wq/records_with_species.csv.gz")
rec["test_cas"]=rec.test_cas.astype(str); meas["test_cas"]=meas.test_cas.astype(str)
rec=rec.drop_duplicates(subset=["test_id","ep","conc1_mean","conc1_unit"])
rec=rec.merge(meas[["test_cas","parent","name"]],on="test_cas",how="inner")
ACUTE={"LC50","EC50","IC50"}; rec["is_acute"]=rec.ep.isin(ACUTE)
pd.set_option("display.width",260); pd.set_option("display.max_colwidth",34)

for cas,lbl in [("335671","PFOA"),("1763231","PFOS"),("375735","PFBS"),("13252136","HFPO-DA"),("307244","PFHxA"),("375224","PFBA")]:
    p=meas.loc[meas.test_cas==cas,"parent"]
    if not len(p): print(lbl,"no parent"); continue
    g=rec[rec.parent==p.iloc[0]]
    print(f"\n{'='*100}\n{lbl}  n_records={len(g)}  acute={g.is_acute.sum()}  "
          f"min(all)={g.logTox.min():.3f}  min(acute)={g[g.is_acute].logTox.min():.3f}  "
          f"median(acute)={g[g.is_acute].logTox.median():.3f}")
    print("  5 lowest records (all endpoints):")
    print(g.nsmallest(5,"logTox")[["ep","logTox","conc1_mean","conc1_unit","latin_name","ecotox_group"]].to_string(index=False))
