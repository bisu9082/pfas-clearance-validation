#!/usr/bin/env python3
"""Parent-level recuration WITH species + ECOTOX trophic group, and EU-TGD AF tiering.
   Step 0 reproduces the deposited 87-parent table exactly, as a correctness gate."""
import pandas as pd, numpy as np, json
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from scipy import stats
RDLogger.DisableLog('rdApp.*')
UNC = rdMolStandardize.Uncharger()
R = "/home/claude/pfasv32"

def canon(s):
    m = Chem.MolFromSmiles(str(s))
    if m is None: return None
    fr = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    fl = [f for f in fr if any(a.GetSymbol()=="F" for a in f.GetAtoms())]
    big = max(fl if fl else fr, key=lambda f: f.GetNumHeavyAtoms())
    try: big = UNC.uncharge(big); Chem.SanitizeMol(big)
    except Exception: pass
    return Chem.MolToSmiles(big)

meas = pd.read_csv(f"{R}/archive_v2/ecotox_pfas_measured.csv")
rec  = pd.read_csv("/home/claude/wq/records_with_species.csv.gz")   # 15,991 + species
meas["parent"] = [canon(s) for s in meas.canonical_smiles]
meas = meas[meas.parent.notna()].copy()
rec["test_cas"]  = rec.test_cas.astype(str)
meas["test_cas"] = meas.test_cas.astype(str)

key = ["test_id","ep","conc1_mean","conc1_unit"]
rec = rec.drop_duplicates(subset=key).copy()
rec = rec.merge(meas[["test_cas","parent"]], on="test_cas", how="inner")
print("records after dedup + parent map: %d" % len(rec))

ACUTE = {"LC50","EC50","IC50"}
rec["is_acute"] = rec.ep.isin(ACUTE)
rec["grp"] = rec.ecotox_group.fillna("").str.split(";").str[0].str.strip()

# EU-TGD style trophic buckets
def trophic(g):
    if g in ("Fish",): return "fish"
    if g in ("Crustaceans",): return "invertebrate"
    if g in ("Algae","Moss, Hornworts"): return "algae"
    if g in ("Insects/Spiders","Molluscs","Worms","Other Invertebrates",
             "Invertebrates"): return "invertebrate"
    return "other"
rec["troph"] = rec.grp.map(trophic)

rows = []
for p, g in rec.groupby("parent"):
    ga = g[g.is_acute]
    sp_med = (ga.groupby("latin_name").logTox.median().median()
              if len(ga) else np.nan)
    tl = set(ga.troph) - {"other"}
    rows.append(dict(parent=p, n_records=len(g), n_tests=g.test_id.nunique(),
        n_acute=len(ga), logTox=g.logTox.min(),
        logTox_median=ga.logTox.median() if len(ga) else np.nan,
        logTox_sp_median=sp_med,
        n_species=ga.latin_name.nunique(), n_trophic=len(tl),
        has_chronic=int((~g.is_acute).any()),
        n_chronic_species=g[~g.is_acute].latin_name.nunique(),
        name=meas.loc[meas.parent==p,"name"].iloc[0]))
D = pd.DataFrame(rows)

# ---------- GATE: reproduce the deposited table ----------
dep = pd.read_csv(f"{R}/canonical_results/ecotox_pfas_measured_recurated.csv")
mg = D.merge(dep, on="parent", suffixes=("_new","_dep"))
ok = (len(D)==len(dep) and len(mg)==len(dep)
      and (mg.n_records_new==mg.n_records_dep).all()
      and (mg.logTox_new==mg.logTox_dep).all()
      and np.allclose(mg.logTox_median_new, mg.logTox_median_dep, equal_nan=True))
print("GATE reproduce deposited 87-parent table: %s  (%d parents, %d with acute)"
      % ("PASS" if ok else "FAIL", len(D), int(D.n_acute.gt(0).sum())))
assert ok, "recuration does not reproduce the deposited table"

A = D[D.n_acute>0].copy()
print("\n" + "="*70)
print("Species resolution at PARENT level (n = %d)" % len(A))
print("="*70)
A["delta"] = A.logTox_sp_median - A.logTox_median
print("  |species-resolved - record-resolved| median %.3f log10, 90th %.3f, max %.3f"
      % (A.delta.abs().median(), A.delta.abs().quantile(.9), A.delta.abs().max()))
print("  identical (<0.01 log): %d / %d" % ((A.delta.abs()<.01).sum(), len(A)))
x = np.log10(A.n_records); rng = np.random.default_rng(20260620)
for c,l in [("logTox","protective minimum"),
            ("logTox_median","acute median, record-resolved"),
            ("logTox_sp_median","acute median, SPECIES-resolved")]:
    y=A[c].values; r,p = stats.spearmanr(x,y)
    bs=[stats.spearmanr(x.values[i],y[i])[0] for i in (rng.integers(0,len(y),len(y)) for _ in range(4000))]
    lo,hi=np.percentile(bs,[2.5,97.5])
    print("  %-34s rho=%+.3f  95%% CI [%+.3f,%+.3f]  p=%.3g" % (l,r,lo,hi,p))

print("\n" + "="*70)
print("EU-TGD assessment factor tiering (what the data can actually support)")
print("="*70)
def af(r):
    if r.n_chronic_species >= 3: return 10
    if r.n_chronic_species == 2: return 50
    if r.n_chronic_species == 1: return 100
    return 1000
A["AF"] = A.apply(af, axis=1)
print(A.AF.value_counts().sort_index().rename("substances").to_string())
print("\n  trophic levels covered by ACUTE data:")
print(A.n_trophic.value_counts().sort_index().rename("substances").to_string())
print("  substances meeting the 3-trophic-level base requirement: %d / %d"
      % ((A.n_trophic>=3).sum(), len(A)))
r,p = stats.spearmanr(np.log10(A.n_records), A.AF)
print("\n  Spearman(log records, AF) = %+.3f  p = %.3g   <- does AF itself track effort?" % (r,p))

# PNEC with tiered AF vs fixed AF
for tag, afv in [("fixed AF=1000", np.full(len(A),1000.0)), ("tiered AF", A.AF.values.astype(float))]:
    pmin = 10**A.logTox.values*1000/afv
    pmed = 10**A.logTox_sp_median.values*1000/afv
    rmin = stats.spearmanr(np.log10(A.n_records), np.log10(pmin))
    rmed = stats.spearmanr(np.log10(A.n_records), np.log10(pmed))
    print("\n  [%s]  PNEC(protective min) vs effort rho=%+.3f (p=%.2g)"
          "\n  %s  PNEC(species median) vs effort rho=%+.3f (p=%.2g)"
          % (tag, rmin[0], rmin[1], " "*(len(tag)+4), rmed[0], rmed[1]))

A.to_csv("/home/claude/wq/parent_species_af.csv", index=False)
print("\nwrote /home/claude/wq/parent_species_af.csv")
