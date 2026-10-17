"""Equal-depth subsampling with the eligible-subset baseline the original omitted."""
import numpy as np, pandas as pd, json
from scipy.stats import spearmanr
rng = np.random.default_rng(20260620)
rec = pd.read_csv('cr/ecotox_pfas_records.csv.gz').drop_duplicates(
        subset=['test_id','ep','conc1_mean','conc1_unit'])
meas = pd.read_csv('cr/ecotox_pfas_measured_recurated.csv')
src  = pd.read_csv('cr/ecotox_pfas_measured.csv')
from rdkit import Chem, RDLogger; from rdkit.Chem.MolStandardize import rdMolStandardize
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
src['parent']=[canon(s) for s in src.canonical_smiles]
rec=rec.merge(src[['test_cas','parent']].assign(test_cas=lambda d:d.test_cas.astype(str)),
              left_on=rec.test_cas.astype(str), right_on='test_cas', how='inner')
pools={p:g.logTox.values for p,g in rec.groupby('parent')}
full={p:v.min() for p,v in pools.items()}
nrec={p:len(v) for p,v in pools.items()}

print(f'{"k":>3} {"eligible":>9} {"baseline rho":>13} {"depth-equal rho":>16} {"depth effect":>13}')
rows=[]
for k in (3,5,10,20):
    elig=[p for p,v in pools.items() if len(v)>=k]
    base=spearmanr([full[p] for p in elig], np.log10([nrec[p] for p in elig]))[0]
    rs=[]
    for _ in range(400):
        mins=[rng.choice(pools[p],k,replace=False).min() for p in elig]
        rs.append(spearmanr(mins, np.log10([nrec[p] for p in elig]))[0])
    sub=float(np.mean(rs)); sd=float(np.std(rs))
    print(f'{k:>3} {len(elig):>9} {base:>13.3f} {sub:>12.3f}±{sd:.3f} {base-sub:>13.3f}')
    rows.append(dict(k=k,n_eligible=len(elig),rho_baseline=float(base),
                     rho_depth_equal=sub,sd=sd,depth_effect=float(base-sub)))
allrho=spearmanr([full[p] for p in pools], np.log10([nrec[p] for p in pools]))[0]
print(f'\nall {len(pools)} substances, full records: rho = {allrho:.3f}')
sh=[r['depth_effect']/(allrho-0) for r in rows]
print('depth effect as share of the full-sample rho:', ' '.join(f'{abs(s):.0%}' for s in sh))
json.dump(dict(rho_all=float(allrho), rows=rows), open('cr/subsample.json','w'), indent=2)
