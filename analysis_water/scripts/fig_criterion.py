#!/usr/bin/env python3
import pandas as pd, numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

d = pd.read_csv("/home/claude/wq/criterion_per_substance.csv")
J = json.load(open("/home/claude/wq/criterion.json"))
x = d["n_records"].values
pmin, pmed = d["pnec_min_ugL"].values, d["pnec_median_ugL"].values

fig, (a, b) = plt.subplots(1, 2, figsize=(11.0, 4.3))

# ---- (a) PNEC vs testing effort, both conventions
for y, c, m, lab, rho in [
    (pmin, "#C0392B", "o", "protective minimum", J["spearman_effort_vs_pnec_min"][0]),
    (pmed, "#2471A3", "s", "effort-stable acute median", J["spearman_effort_vs_pnec_median"][0])]:
    a.scatter(x, y, s=26, c=c, marker=m, alpha=.72, edgecolors="white", linewidths=.5,
              label=r"%s  ($\rho$ = %+.2f)" % (lab, rho), zorder=3)
    lx, ly = np.log10(x), np.log10(y)
    sl, ic = np.polyfit(lx, ly, 1)
    xs = np.linspace(lx.min(), lx.max(), 50)
    a.plot(10**xs, 10**(sl*xs + ic), c=c, lw=1.8, zorder=2)
a.set_xscale("log"); a.set_yscale("log")
a.set_xlabel("ECOTOX records per substance")
a.set_ylabel(r"derived PNEC ($\mu$g L$^{-1}$),  AF = 1000")
a.set_title("(a) The derived criterion inherits the artefact", fontsize=10.5, loc="left")
a.legend(fontsize=8, loc="upper right", framealpha=.95)
a.grid(alpha=.25, lw=.5, zorder=0)

# ---- (b) classification reversal
bm = J["benchmarks"]
labs = ["EU EQS\nPFOS\n0.65 ng/L", "US EPA MCL\nPFOA/PFOS\n4 ng/L",
        "screening\n1 $\\mu$g/L", "screening\n10 $\\mu$g/L"]
n = J["n"]
rev = [k["reversed"] / n * 100 for k in bm]
bars = b.bar(range(4), rev, color="#7D3C98", alpha=.85, width=.62, zorder=3)
for i, (bar, k) in enumerate(zip(bars, bm)):
    b.text(bar.get_x() + bar.get_width()/2, rev[i] + 1.4,
           "%d / %d" % (k["reversed"], n), ha="center", fontsize=8.5, fontweight="bold")
b.set_xticks(range(4)); b.set_xticklabels(labs, fontsize=8)
b.set_ylabel("substances whose classification reverses (%)")
b.set_ylim(0, max(rev) * 1.28)
b.set_title("(b) Which side of the benchmark depends on the convention",
            fontsize=10.5, loc="left")
b.grid(axis="y", alpha=.25, lw=.5, zorder=0)

plt.tight_layout()
plt.savefig("/home/claude/wq/fig_criterion.png", dpi=190)

# --- assertion: plotted rho must match the computed one
rr = stats.spearmanr(np.log10(x), np.log10(pmin))[0]
assert abs(rr - J["spearman_effort_vs_pnec_min"][0]) < 1e-9, "panel (a) rho mismatch"
assert abs(rev[0] - bm[0]["reversed"]/n*100) < 1e-9, "panel (b) height mismatch"
print("figure written; label-value assertions passed")
print("reversal %%:", [round(v,1) for v in rev])
