import pandas as pd, numpy as np, re, sys
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog('rdApp.*')
UNC = rdMolStandardize.Uncharger()
def canon(s):
    m = Chem.MolFromSmiles(str(s))
    if m is None: return None
    fr = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    fl = [f for f in fr if any(a.GetSymbol()=="F" for a in f.GetAtoms())]
    big = max(fl if fl else fr, key=lambda f: f.GetNumHeavyAtoms())
    try: big = UNC.uncharge(big); Chem.SanitizeMol(big)
    except Exception: pass
    return Chem.MolToSmiles(big)

M = {
"Perfluorobutanoate":"375-22-4","Perfluorobutanoic acid":"375-22-4",
"Heptafluorobutyric acid***retired***use Perfluorobutanoic acid":"375-22-4",
"PFOA ion":"335-67-1","Perfluorooctanoic acid":"335-67-1","Perfluorooctanoate (anionic form)":"335-67-1",
"Perfluorooctanesulfonate":"1763-23-1","Perfluorooctane sulfonic acid":"1763-23-1",
"Perfluorohexanoate":"307-24-4","Perfluorohexanoic acid":"307-24-4",
"Perfluorobutanesulfonate":"375-73-5","Perfluorobutanesulfonic acid":"375-73-5",
"Perfluoropentanoate":"2706-90-3","Perfluorovaleric acid":"2706-90-3","Perfluoropentanoic acid":"2706-90-3",
"Perfluoroheptanoate":"375-85-9","Perfluoroheptanoic acid":"375-85-9",
"Perfluorohexanesulfonate":"355-46-4","Perfluorohexanesulfonic acid":"355-46-4","Perfluorohexane sulfonate":"355-46-4",
"Perfluorononanoate":"375-95-1","Perfluorononanoic acid":"375-95-1",
"Perfluorooctanesulfonamide":"754-91-6",
"Perfluorodecanoic acid":"335-76-2","Perfluorodecanoate":"335-76-2",
"Perfluoropentane sulfonic acid":"2706-91-4","Perfluoropentanesulfonate":"2706-91-4",
"6:2 Fluorotelomer sulfonate acid":"27619-97-2",
"6:2 Fluorotelomer sulfonate***retired***use 6:2 Fluorotelomer sulfonate acid":"27619-97-2",
"1-Octanesulfonic acid, 3,3,4,4,5,5,6,6,7,7,8,8,8-tridecafluoro-***retired***use 6:2 Fluorotelomer sulfonate acid":"27619-97-2",
"Perfluorododecanoic acid":"307-55-1","Perfluorododecanoate":"307-55-1",
"Perfluoro-3-methoxypropanoic acid":"377-73-1",
"Perfluoroundecanoate":"2058-94-8","Perfluoroundecanoic acid":"2058-94-8",
"3:3 Fluorotelomer carboxylic acid":"356-02-5",
"Perfluoro(4-methoxybutanoic) acid":"863090-89-5",
"Fluorotelomer sulfonate 4:2":"757124-72-4",
"4,8-Dioxa-3H-perfluorononanoic acid":"919005-14-4",
"Perfluoro-3,6-dioxaheptanoic acid":"151772-58-6",
"Nonafluoro-3,6-dioxaheptanoic acid***retired***use Perfluoro-3,6-dioxaheptanoic acid":"151772-58-6",
"N-ethyl perfluorooctanesulfonamidoacetic acid":"2991-50-6",
"N-Ethylperfluorooctanesulfonamidoacetate":"2991-50-6",
"N-Ethyl-perfluorooctane sulfonamido acetic acid":"2991-50-6",
"Fluorotelomer sulfonate 8:2":"39108-34-4",
"Perfluorotridecanoic acid":"72629-94-8","Perfluorotridecanoate":"72629-94-8",
"Perfluorotetradecanoic acid":"376-06-7","Perfluorotetradecanoate":"376-06-7",
"N-methyl perfluorooctanesulfonamidoacetic acid":"2355-31-9",
"N-Methylperfluorooctanesulfonamidoacetate":"2355-31-9",
"Perfluoroheptanesulfonate":"375-92-8",
"Perfluoro(2-ethoxyethane)sulfonic acid":"113507-82-7",
"Perfluoro-1-nonanesulfonic acid":"68259-12-1","Perfluorononanesulfonate":"68259-12-1",
"Perfluorohexanesulfonamide":"41997-13-1",
"3-Perfluoroheptylpropanoic acid (7:3 FTCA)":"812-70-4",
"Perfluorobutylsulfonamide":"30334-69-1",
"3-Perfluoropentylpropanoic acid (5:3 FTCA)":"914637-49-3","3-perfluoropentyl propanoate":"914637-49-3",
"N-Ethyl-N-(2-hydroxyethyl)perfluorooctanesulfonamide":"1691-99-2",
"N-Ethyl perfluorooctane sulfonamido ethanol":"1691-99-2",
"2-(N-methylperfluoro-1-octanesulfonamido)-ethanol":"24448-09-7",
"N-methylperfluoro-1-octanesulfonamide":"31506-32-8",
"N-ethyl perfluoro-1-octanesulfonamide":"4151-50-2",
"Perfluoro(2-propoxypropanoate)":"13252-13-6",
"Perfluoropalmitic acid":"67905-19-5",
"Perfluorodecanesulfonate":"335-77-3",
}
DROP = ["N-Methylperfluorooctanesulfonamide (N-MeFOSA) (D3, 98%)",
        "N-Ethylperfluorooctanesulfonamidoacetic acid (N-EtFOSAA) (N-ethyl-D5, 98%)"]

meas = pd.read_csv("/home/claude/pfasv32/archive_v2/ecotox_pfas_measured.csv")
meas["test_cas"]=meas.test_cas.astype(str)
# ECOTOX test_cas is CAS without dashes
def nodash(c): return c.replace("-","")
smi = dict(zip(meas.test_cas, meas.canonical_smiles))

conc = pd.read_pickle("conc.pkl")
conc = conc[~conc.CharacteristicName.isin(DROP)].copy()
chars = sorted(conc.CharacteristicName.unique())
rows=[]
for c in chars:
    cas = M.get(c)
    s = smi.get(nodash(cas)) if cas else None
    p = canon(s) if s else None
    f = rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(p)) if p else ""
    rows.append(dict(char=c,cas=cas,in_ecotox=s is not None,parent=p,formula=f,
                     n=int((conc.CharacteristicName==c).sum())))
t=pd.DataFrame(rows)
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",60)
print(t[["n","char","cas","in_ecotox","formula"]].sort_values("n",ascending=False).to_string(index=False))
print("\nunmapped names:", t[t.cas.isna()].char.tolist())
print("mapped but absent from ECOTOX:", t[(~t.cas.isna())&(~t.in_ecotox)].cas.nunique(),
      sorted(set(t[(~t.cas.isna())&(~t.in_ecotox)].cas)))
t.to_csv("char_map.csv",index=False)
