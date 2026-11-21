#!/usr/bin/env python3
"""R4-1: is the acute median driven by repeated testing in one species?
   And does a species-resolved endpoint still track testing effort?"""
import pandas as pd, numpy as np
from scipy import stats

m = pd.read_csv("/home/claude/wq/records_with_species.csv.gz")
cur = pd.read_csv("/home/claude/pfasv32/canonical_results/ecotox_pfas_measured_recurated.csv")

ACUTE = {"LC50","EC50","IC50"}
m["is_acute"] = m.ep.isin(ACUTE)
m["grp"] = m.ecotox_group.fillna("").str.split(";").str[0]

print("="*70); print("1. Which species actually dominate?"); print("="*70)
a = m[m.is_acute]
print("acute records: %d across %d species" % (len(a), a.latin_name.nunique()))
top = a.latin_name.value_counts().head(8)
for s, n in top.items():
    print("   %-34s %5d  (%4.1f%%)" % (s, n, 100*n/len(a)))
print("\n   R4 named Daphnia magna. Its share of acute records: %.1f%%"
      % (100*(a.latin_name=="Daphnia magna").sum()/len(a)))

print("\n" + "="*70)
print("2. Per-substance dominance of the single most-tested species")
print("="*70)
dom = (a.groupby(["test_cas","latin_name"]).size().rename("n").reset_index()
        .sort_values("n", ascending=False).groupby("test_cas")
        .agg(top_n=("n","first"), tot=("n","sum"), nsp=("n","size")).reset_index())
dom["share"] = dom.top_n/dom.tot
print("substances with acute data: %d" % len(dom))
print("   median share of the single most-tested species: %.0f%%" % (100*dom.share.median()))
print("   substances where one species is >=50%% of acute records: %d / %d"
      % ((dom.share>=.5).sum(), len(dom)))
print("   substances tested on only ONE species:            %d / %d"
      % ((dom.nsp==1).sum(), len(dom)))

print("\n" + "="*70)
print("3. Record-median vs species-resolved median")
print("="*70)
rec_med = a.groupby("test_cas").logTox.median().rename("rec_median")
sp_med  = (a.groupby(["test_cas","latin_name"]).logTox.median()
            .groupby("test_cas").median().rename("sp_median"))
nrec    = m.groupby("test_cas").size().rename("n_records")
nsp     = a.groupby("test_cas").latin_name.nunique().rename("n_species")
prot    = m.groupby("test_cas").logTox.min().rename("prot_min")
D = pd.concat([rec_med, sp_med, nrec, nsp, prot], axis=1).dropna(subset=["rec_median"])
D["delta"] = D.sp_median - D.rec_median
print("n = %d substances" % len(D))
print("   |species-resolved - record| median %.3f log10, 90th pct %.3f, max %.3f"
      % (D.delta.abs().median(), D.delta.abs().quantile(.9), D.delta.abs().max()))
print("   identical (<0.01 log): %d / %d" % ((D.delta.abs()<.01).sum(), len(D)))

print("\n--- effort dependence of each endpoint (log10 record count) ---")
x = np.log10(D.n_records)
rng = np.random.default_rng(20260620)
for col, lab in [("prot_min","protective minimum"),
                 ("rec_median","acute median, record-resolved"),
                 ("sp_median","acute median, SPECIES-resolved")]:
    y = D[col].values
    r, p = stats.spearmanr(x, y)
    bs = [stats.spearmanr(x.values[i], y[i])[0]
          for i in (rng.integers(0, len(y), len(y)) for _ in range(4000))]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    print("   %-34s rho = %+.3f  95%% CI [%+.3f, %+.3f]  p = %.3g"
          % (lab, r, lo, hi, p))

print("\n--- does taxonomic coverage itself track effort? (R4's actual worry) ---")
r, p = stats.spearmanr(np.log10(D.n_records), D.n_species)
print("   Spearman(log records, n species) = %+.3f  p = %.3g" % (r, p))
r2, p2 = stats.spearmanr(D.n_species, D.sp_median)
print("   Spearman(n species, species-resolved median) = %+.3f  p = %.3g" % (r2, p2))

D.to_csv("/home/claude/wq/species_resolved.csv")
print("\nwrote /home/claude/wq/species_resolved.csv")
