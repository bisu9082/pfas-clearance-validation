#!/usr/bin/env python3
"""EU TGD / ECHA R.10 assessment factors, implemented on trophic-level coverage.

  AF 10   : long-term NOECs from >=3 species representing 3 trophic levels
  AF 50   : long-term NOECs from 2 species representing 2 trophic levels
  AF 100  : one long-term NOEC
  AF 1000 : short-term L(E)C50 from each of 3 trophic levels (base set)
  base set not met -> PNEC not derivable under the standard scheme
"""
import pandas as pd, numpy as np, json
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from scipy import stats
RDLogger.DisableLog('rdApp.*'); UNC = rdMolStandardize.Uncharger()
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
rec  = pd.read_csv("/home/claude/wq/records_with_species.csv.gz")
meas["parent"] = [canon(s) for s in meas.canonical_smiles]
meas = meas[meas.parent.notna()].copy()
for df in (rec, meas): df["test_cas"] = df.test_cas.astype(str)
rec = rec.drop_duplicates(subset=["test_id","ep","conc1_mean","conc1_unit"])
rec = rec.merge(meas[["test_cas","parent"]], on="test_cas", how="inner")
assert len(rec) == 4744, len(rec)

ACUTE = {"LC50","EC50","IC50"}
rec["is_acute"] = rec.ep.isin(ACUTE)
rec["grp"] = rec.ecotox_group.fillna("").str.split(";").str[0].str.strip()

STRICT = {"Fish":"fish", "Crustaceans":"crustacean", "Algae":"algae"}
BROAD  = dict(STRICT, **{"Insects/Spiders":"invertebrate","Molluscs":"invertebrate",
              "Worms":"invertebrate","Other Invertebrates":"invertebrate",
              "Moss, Hornworts":"algae"})
def levels(g, mapping):
    return {mapping[x] for x in g.grp if x in mapping}

rows=[]
for p, g in rec.groupby("parent"):
    ga, gc = g[g.is_acute], g[~g.is_acute]
    if not len(ga): continue
    d = dict(parent=p, n_records=len(g), n_acute=len(ga),
             logTox=g.logTox.min(), logTox_median=ga.logTox.median(),
             logTox_sp_median=ga.groupby("latin_name").logTox.median().median())
    for tag, mp in [("strict", STRICT), ("broad", BROAD)]:
        d[f"acute_tl_{tag}"]   = len(levels(ga, mp))
        d[f"chronic_tl_{tag}"] = len(levels(gc, mp))
    d["n_chronic_species"] = gc.latin_name.nunique()
    rows.append(d)
A = pd.DataFrame(rows)
print("substances with acute data: %d" % len(A))

def assign_af(r, tag):
    ctl = r[f"chronic_tl_{tag}"]
    if ctl >= 3: return 10.0
    if ctl == 2: return 50.0
    if ctl == 1 or r.n_chronic_species >= 1: return 100.0
    return 1000.0

for tag in ("strict","broad"):
    A[f"AF_{tag}"]   = A.apply(lambda r: assign_af(r, tag), axis=1)
    A[f"base_{tag}"] = A[f"acute_tl_{tag}"] >= 3

for tag in ("strict","broad"):
    print("\n" + "="*70)
    print("TROPHIC MAPPING: %s   (%s)" % (tag.upper(),
          "fish / crustacean / algae" if tag=="strict" else "+ other invertebrates, moss"))
    print("="*70)
    print("  acute trophic levels covered:")
    print(A[f"acute_tl_{tag}"].value_counts().sort_index().to_string())
    print("  base set (>=3 acute trophic levels) met: %d / %d"
          % (A[f"base_{tag}"].sum(), len(A)))
    print("\n  assigned AF:")
    print(A[f"AF_{tag}"].value_counts().sort_index().astype(int).to_string())
    r,p = stats.spearmanr(np.log10(A.n_records), A[f"AF_{tag}"])
    print("\n  Spearman(log records, AF) = %+.3f  p = %.3g" % (r,p))

    for sub_lab, sub in [("ALL substances with acute data", A),
                         ("base-set-compliant only", A[A[f"base_{tag}"]])]:
        if len(sub) < 8: 
            print("  [%s] n=%d too small" % (sub_lab, len(sub))); continue
        lx = np.log10(sub.n_records)
        print("\n  [%s]  n = %d" % (sub_lab, len(sub)))
        for col, lab in [("logTox","protective minimum"),
                         ("logTox_sp_median","species-resolved acute median")]:
            for aftag, afv in [("AF=1000 fixed", np.full(len(sub),1000.0)),
                               ("AF tiered",     sub[f"AF_{tag}"].values)]:
                pn = 10**sub[col].values*1000/afv
                rr,pp = stats.spearmanr(lx, np.log10(pn))
                print("     PNEC %-30s %-14s rho=%+.3f  p=%.3g" % (lab, aftag, rr, pp))

A.to_csv("/home/claude/wq/af_proper.csv", index=False)
print("\nwrote /home/claude/wq/af_proper.csv")
