#!/usr/bin/env python3
import pandas as pd, numpy as np, json, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
A = pd.read_csv("/home/claude/wq/af_proper.csv")
rng = np.random.default_rng(20260620)
def boot(x,y,n=6000):
    r=stats.spearmanr(x,y)[0]
    bs=[stats.spearmanr(x[i],y[i])[0] for i in (rng.integers(0,len(x),len(x)) for _ in range(n))]
    return (r,)+tuple(np.percentile(bs,[2.5,97.5]))

lx = np.log10(A.n_records).values
af = A.AF_strict.values
pmin = 10**A.logTox.values*1000/af
pmed = 10**A.logTox_sp_median.values*1000/af
r1 = boot(lx, np.log10(pmin)); r2 = boot(lx, np.log10(pmed))

# --- E1's demand: CI on the base-set-compliant subset ---
B = A[A.base_strict]
lb = np.log10(B.n_records).values
pb = 10**B.logTox.values*1000/B.AF_strict.values
r3 = boot(lb, np.log10(pb))
print("base-set-compliant (n=%d) protective-min PNEC vs effort: rho=%+.3f CI [%+.3f, %+.3f]"
      % (len(B), *r3))

fig, ax = plt.subplots(1, 3, figsize=(13.6, 4.2))

# (a) PNEC vs effort, tiered AF
a = ax[0]
for y, c, m, lab, rr in [(pmin,"#C0392B","o","protective minimum", r1),
                         (pmed,"#2471A3","s","effort-stable acute median\n(species-resolved)", r2)]:
    a.scatter(A.n_records, y, s=26, c=c, marker=m, alpha=.72,
              edgecolors="white", linewidths=.5, zorder=3,
              label=r"%s" "\n" r"$\rho$ = %+.2f  [%+.2f, %+.2f]" % (lab, rr[0], rr[1], rr[2]))
    sl, ic = np.polyfit(lx, np.log10(y), 1)
    xs = np.linspace(lx.min(), lx.max(), 50)
    sig = not (rr[1] <= 0 <= rr[2])          # dashed when the CI spans zero
    a.plot(10**xs, 10**(sl*xs+ic), c=c, lw=1.8, zorder=2,
           ls="-" if sig else (0, (5, 3)))
a.set_xscale("log"); a.set_yscale("log")
a.set_xlabel("ECOTOX records per substance")
a.set_ylabel(r"derived PNEC ($\mu$g L$^{-1}$), tiered AF")
a.set_title("(a) The derived criterion tracks testing effort", fontsize=10.5, loc="left")
a.legend(fontsize=7.2, loc="lower left", framealpha=.95)
a.grid(alpha=.25, lw=.5, zorder=0)

# (b) AF itself vs effort
b = ax[1]
rA = boot(lx, af.astype(float))
for lv, c in zip([10, 50, 100, 1000], ["#1B4F72","#2E86C1","#A9CCE3","#D5D8DC"]):
    s = A[A.AF_strict == lv]
    b.scatter(s.n_records, np.full(len(s), lv), s=40, c=c, alpha=.85,
              edgecolors="#34495E", linewidths=.6, zorder=3, label="AF %d  (n=%d)" % (lv, len(s)))
b.set_xscale("log"); b.set_yscale("log")
b.set_yticks([10,50,100,1000]); b.set_yticklabels(["10","50","100","1000"])
b.set_xlabel("ECOTOX records per substance"); b.set_ylabel("EU TGD assessment factor")
b.set_title(r"(b) The AF rule is itself effort-dependent" "\n"
            r"    $\rho$ = %+.2f  [%+.2f, %+.2f]" % (rA[0], rA[1], rA[2]),
            fontsize=10.5, loc="left")
b.legend(fontsize=7.5, loc="upper right", framealpha=.95)
b.grid(alpha=.25, lw=.5, zorder=0)

# (c) base-set compliance
c = ax[2]
w = .38; idx = np.arange(4)
for off, tag, col, lab in [(-w/2,"strict","#884EA0","fish / crustacean / algae"),
                           ( w/2,"broad", "#48C9B0","+ other invertebrates")]:
    v = A[f"acute_tl_{tag}"].clip(upper=3).value_counts().reindex(range(4), fill_value=0)
    c.bar(idx+off, v.values, w, color=col, alpha=.88, zorder=3, label=lab)
    for i, n in enumerate(v.values):
        if n: c.text(i+off, n+.6, str(n), ha="center", fontsize=8)
c.set_xticks(idx); c.set_xticklabels(["0", "1", "2", "3 (base set met)"])
c.set_xlabel("acute trophic levels covered"); c.set_ylabel("substances")
c.set_title("(c) Most PFAS cannot support a standard PNEC\n"
            "    base set met: 10/58 strict, 17/58 broad", fontsize=10.5, loc="left")
c.legend(fontsize=7.5, framealpha=.95); c.grid(axis="y", alpha=.25, lw=.5, zorder=0)
c.set_ylim(0, 36)

plt.tight_layout(); plt.savefig("/home/claude/wq/fig_pnec.png", dpi=190)
# assertions: plotted rho == computed rho; bar heights == counts
assert abs(stats.spearmanr(lx, np.log10(pmin))[0] - r1[0]) < 1e-12
assert A[A.acute_tl_strict>=3].shape[0] == 10 and A[A.acute_tl_broad>=3].shape[0] == 17
print("figure written; assertions passed")
