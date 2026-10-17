"""Regenerate the per-substance PFAS table with the curation rule the Methods state.

Fixes: (1) parents were never collapsed in ecotox_pfas_measured.csv;
       (2) effort was counted on raw rows, 70% of which are exact duplicates.
"""
import gzip, json, sys
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')
UNC = rdMolStandardize.Uncharger()

def canon(s):
    m = Chem.MolFromSmiles(str(s))
    if m is None: return None
    fr = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    fl = [f for f in fr if any(a.GetSymbol()=="F" for a in f.GetAtoms())]
    big = max(fl if fl else fr, key=lambda f: f.GetNumHeavyAtoms())
    try:
        big = UNC.uncharge(big); Chem.SanitizeMol(big)
    except Exception: pass
    return Chem.MolToSmiles(big)

meas = pd.read_csv('cr/ecotox_pfas_measured.csv')
rec  = pd.read_csv('cr/ecotox_pfas_records.csv.gz')
print(f'input: {len(meas)} CAS rows, {len(rec)} record rows')

meas['parent'] = [canon(s) for s in meas.canonical_smiles]
meas = meas[meas.parent.notna()].copy()
print(f'distinct parents: {meas.parent.nunique()}  (collapsed {len(meas)-meas.parent.nunique()} redundant rows)')

# --- deduplicate records -------------------------------------------------
rec['test_cas'] = rec.test_cas.astype(str)
meas['test_cas'] = meas.test_cas.astype(str)
key = ['test_id','ep','conc1_mean','conc1_unit']
before = len(rec)
rec = rec.drop_duplicates(subset=key).copy()
print(f'records: {before} -> {len(rec)} after exact-duplicate removal '
      f'({100*(1-len(rec)/before):.1f}% removed); unique tests {rec.test_id.nunique()}')

rec = rec.merge(meas[['test_cas','parent']], on='test_cas', how='inner')
print(f'records mapped to a parent: {len(rec)}')

ACUTE = {'LC50','EC50','IC50'}   # matches Methods &2.2; no LD50 rows occur in the freshwater extraction anyway
rows = []
for p, g in rec.groupby('parent'):
    ga = g[g.ep.isin(ACUTE)]
    rows.append(dict(
        parent=p,
        n_records=len(g),
        n_tests=g.test_id.nunique(),
        n_acute=len(ga),
        logTox=g.logTox.min(),                       # protective minimum
        logTox_median=ga.logTox.median() if len(ga) else np.nan,
        name=meas.loc[meas.parent==p,'name'].iloc[0],
    ))
out = pd.DataFrame(rows)
out.to_csv('cr/ecotox_pfas_measured_recurated.csv', index=False)
print(f'\nparent-level table: {len(out)} substances, {out.n_acute.gt(0).sum()} with acute records')

def rho(y, x, label):
    m = np.isfinite(y) & np.isfinite(x)
    r, p = spearmanr(y[m], x[m])
    print(f'  {label:<48} rho={r:+.3f}  p={p:.3g}  n={m.sum()}')
    return r, p, int(m.sum())

print('\n=== effort associations (parent-collapsed, de-duplicated) ===')
res = {}
res['min_vs_records'] = rho(out.logTox.values, np.log10(out.n_records.values),
                            'minimum vs log10(records)')
res['min_vs_tests']   = rho(out.logTox.values, np.log10(out.n_tests.values),
                            'minimum vs log10(unique tests)')
res['med_vs_records'] = rho(out.logTox_median.values, np.log10(out.n_records.values),
                            'acute median vs log10(records)')
sub = out[out.n_acute>0]
res['med_vs_acute']   = rho(sub.logTox_median.values, np.log10(sub.n_acute.values),
                            'acute median vs log10(ACUTE records)')
res['min_vs_acute']   = rho(sub.logTox.values, np.log10(sub.n_acute.values),
                            'minimum vs log10(ACUTE records)')

print('\n=== depth strata (minimum) ===')
strata = {}
for lab, m in [('<=3 records', out.n_records<=3), ('4-30', out.n_records.between(4,30)),
               ('31-300', out.n_records.between(31,300)), ('>300', out.n_records>300)]:
    v = out.loc[m,'logTox']
    strata[lab] = dict(n=int(m.sum()), median=float(v.median()))
    print(f'  {lab:<14} n={m.sum():<4} median={v.median():+.2f}')
gap = strata['<=3 records']['median'] - strata['>300']['median']
print(f'  gap = {gap:.2f} log10 units')

print('\n=== offset: acute median minus minimum (paired) ===')
pr = out.dropna(subset=['logTox_median'])
d = pr.logTox_median - pr.logTox
print(f'  n={len(d)}  mean={d.mean():+.3f}  median={d.median():+.3f}')

json.dump(dict(n_parents=int(len(out)), n_with_acute=int(out.n_acute.gt(0).sum()),
               rho=res, strata=strata, gap=float(gap),
               offset_mean=float(d.mean()), offset_median=float(d.median()),
               records_before=int(before), records_after=int(len(rec))),
          open('cr/recuration_metrics.json','w'), indent=2)
print('\nwrote cr/ecotox_pfas_measured_recurated.csv and cr/recuration_metrics.json')
