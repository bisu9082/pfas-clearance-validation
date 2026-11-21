#!/usr/bin/env python3
"""Water-quality criterion derivation under two endpoint conventions."""
import pandas as pd, numpy as np, json
from scipy import stats

D = "/home/claude/pfasv32/canonical_results/ecotox_pfas_measured_recurated.csv"
d = pd.read_csv(D)
d = d[d["n_acute"] > 0].copy()          # acute median only defined here
AF = 1000.0                              # EU TGD: acute-only data, 3 trophic levels

# logTox / logTox_median are log10(mg/L). PNEC in ug/L = 10^log * 1000 / AF
d["pnec_min_ugL"]    = 10**d["logTox"]        * 1000.0 / AF
d["pnec_median_ugL"] = 10**d["logTox_median"] * 1000.0 / AF
d["ratio"]           = d["pnec_median_ugL"] / d["pnec_min_ugL"]
d["log_ratio"]       = d["logTox_median"] - d["logTox"]
d["log_eff"]         = np.log10(d["n_records"])

print("=" * 68)
print("PNEC under two endpoint conventions   (AF = %d, n = %d PFAS)" % (AF, len(d)))
print("=" * 68)
q = d["ratio"].describe(percentiles=[.25, .5, .75])
print("ratio PNEC(acute median) / PNEC(protective minimum)")
for k in ["min", "25%", "50%", "75%", "max"]:
    print("   %-5s %12.4g x" % (k, q[k]))
print("   geometric mean %.4g x" % 10**d["log_ratio"].mean())

print("\nPNEC distribution (ug/L)")
for c, lab in [("pnec_min_ugL", "protective minimum"), ("pnec_median_ugL", "acute median")]:
    s = d[c].describe(percentiles=[.25, .5, .75])
    print("   %-20s median %10.4g   IQR %10.4g - %-10.4g" % (lab, s["50%"], s["25%"], s["75%"]))

# --- does the gap grow with testing effort? ---
r, p = stats.spearmanr(d["log_eff"], d["log_ratio"])
print("\nSpearman(log record count, log10 PNEC ratio) = %+.3f   p = %.3g" % (r, p))
rmin, pmin = stats.spearmanr(d["log_eff"], np.log10(d["pnec_min_ugL"]))
rmed, pmed = stats.spearmanr(d["log_eff"], np.log10(d["pnec_median_ugL"]))
print("Spearman(log record count, log10 PNEC)")
print("   protective minimum  %+.3f  (p = %.3g)" % (rmin, pmin))
print("   acute median        %+.3f  (p = %.3g)" % (rmed, pmed))

# --- decision reversal against real benchmarks ---
print("\n" + "=" * 68)
print("Classification against published water benchmarks")
print("=" * 68)
BM = [("EU EQS, PFOS (Dir. 2013/39/EU)", 0.00065),
      ("US EPA MCL, PFOA/PFOS (2024)",   0.004),
      ("1 ug/L screening level",         1.0),
      ("10 ug/L screening level",        10.0)]
rows = []
for lab, thr in BM:
    a = (d["pnec_min_ugL"]    > thr).sum()    # "clears" the benchmark
    b = (d["pnec_median_ugL"] > thr).sum()
    flip = ((d["pnec_min_ugL"] <= thr) & (d["pnec_median_ugL"] > thr)).sum()
    rows.append((lab, thr, a, b, flip))
    print("  %-32s thr %9.5g ug/L" % (lab, thr))
    print("     clears under protective minimum : %2d / %d" % (a, len(d)))
    print("     clears under acute median       : %2d / %d" % (b, len(d)))
    print("     classification REVERSED         : %2d  (%.0f%%)" % (flip, 100*flip/len(d)))

out = {"n": int(len(d)), "AF": AF,
       "ratio_median": float(d["ratio"].median()),
       "ratio_geomean": float(10**d["log_ratio"].mean()),
       "ratio_max": float(d["ratio"].max()),
       "spearman_effort_vs_logratio": [float(r), float(p)],
       "spearman_effort_vs_pnec_min": [float(rmin), float(pmin)],
       "spearman_effort_vs_pnec_median": [float(rmed), float(pmed)],
       "benchmarks": [{"label": l, "thr_ugL": t, "clear_min": int(a),
                       "clear_median": int(b), "reversed": int(f)} for l, t, a, b, f in rows]}
json.dump(out, open("/home/claude/wq/criterion.json", "w"), indent=1)
d.to_csv("/home/claude/wq/criterion_per_substance.csv", index=False)
print("\nwrote /home/claude/wq/criterion.json and criterion_per_substance.csv")
