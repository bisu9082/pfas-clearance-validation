# analysis_water — effort dependence of PFAS water-quality criteria

Companion analysis to the structure-based PFAS clearance work in this repository.
It asks a separate question: **how much does the derived aquatic criterion for a
PFAS depend on how much that PFAS has been tested, rather than on how toxic it is?**

## Data

| source | what | n |
|---|---|---|
| US EPA ECOTOX `ecotox_ascii_09_15_2026` | aquatic toxicity records, species-joined | 4,744 deduplicated records, 87 parent substances |
| Water Quality Portal (NWIS + WQX), `www.waterqualitydata.us` | measured surface-water concentrations, Stream + Lake/Reservoir/Impoundment, sampleMedia=Water, 2015-01-01 onward | 21,346 quantified values, 1,180 sites, 2016-09-20 .. 2026-06-25 |

Toxicity records are collapsed to parent substances with the same RDKit
`Uncharger` + largest-fluorinated-fragment rule used elsewhere in this repository,
so acid and anion forms of the same PFAS (e.g. "Perfluorooctanoate", "PFOA ion",
"Perfluorooctanoic acid") become one substance. `scripts/parent_species_af.py`
contains an assertion gate that reproduces the deposited 87-parent table exactly
before any new analysis runs.

## Convention ladder

PNEC = 10^(logTox, mg/L) x 1000 / AF, in ug/L. Six ways of choosing `logTox`:

| id | convention | AF |
|---|---|---|
| C1 | lowest value over all endpoints (acute LC50/EC50/IC50 + chronic NOEC/LOEC) | 1000 |
| C2 | as C1 with LOECs removed | 1000 |
| C3 | lowest acute value | 1000 |
| C4 | median acute, record-resolved | 1000 |
| C5 | median acute, species-resolved | 1000 |
| C6 | lowest chronic NOEC, EU-TGD assessment factor on trophic coverage | 10 / 50 / 100 |

`results/exceedance.csv` holds the fraction of measured samples and of monitoring
sites with risk quotient (measured / PNEC) above 1 under each convention.

## Scripts

Run in this order; paths are absolute to the working tree used during the analysis.

    join_species.py        ECOTOX test_id -> species_number -> latin_name, ecotox_group
    parent_species_af.py   parent collapse + reproduction gate + effort correlations
    af_proper.py           EU-TGD assessment factors from trophic coverage
    getchars.py            enumerate the WQP PFAS characteristic-name vocabulary
    pull.py                bulk WQP result download (139 characteristic names)
    map_cas.py             characteristic name -> CAS -> ECOTOX SMILES -> parent
    ladder.py              the six conventions, per parent
    exceed.py              risk quotients and exceedance against measured water
    fig_water.py           figures/fig_water.png

`pull.py` uses `www.waterqualitydata.us`; the `api.` host is not reachable from
every network.

## Note on the characteristic-name mapping

`scripts/map_cas.py` contains a hand-written characteristic-name -> CAS table.
Structures are not taken from that table: the CAS is used only as a key into
ECOTOX's own substance file, and `results/char_map.csv` records the resulting
molecular formula for every mapped name so that name-structure agreement can be
checked directly.

## C7 — reproducing the regulatory screen

Regulators do not use ECOTOX as deposited. They screen it. `scripts/c7c.py`
implements the screen explicitly:

* chronic endpoints only: NOEC, NOEL, MATC, **EC10** (EC10 is not in the deposited
  endpoint set; including it is a deliberate departure, since regulatory chronic
  derivations prefer it)
* **apical effects only**: mortality, growth, reproduction, development, population.
  Gene expression (GEN), biochemistry (BCM), enzyme activity (ENZ), hormones (HRM),
  behaviour (BEH), histology and accumulation are excluded.
* active ingredient only (`conc1_type == "A"`)
* species mean chronic value = geometric mean within species
* HC5 = 5th percentile of a log-normal SSD fitted to the species mean chronic values
  (>= 5 species)

`scripts/c7b.py` carries a second assertion gate: the deposited per-parent minimum
`logTox` must be reproducible from the result-level table for all 87 parents. This
requires the same two rules the deposited pipeline used — strip the `*` qualifier
from endpoint codes, and drop censored values (`conc1_mean_op` in `> < >= <= ~`).

### Result

The screen brings the derived value into the regulatory range: the PFOS HC5 is
4,738 ng/L against Canada's SSD-derived FEQG of 6,800 ng/L, and the PFOA HC5 is
9,516 ng/L against the US EPA 2024 final chronic criterion of 100,000 ng/L.

It does not remove the effort dependence. Spearman rho between log(records) and the
derived criterion:

| convention | rho | 95% CI | n |
|---|---|---|---|
| lowest, all endpoints (as published) | -0.783 | -0.903, -0.586 | 52 |
| lowest acute | -0.196 | -0.509, +0.159 | 41 |
| lowest chronic NOEC, EU-TGD AF | -0.706 | -0.856, -0.479 | 46 |
| lowest apical chronic SMCV (screened) | -0.787 | -0.887, -0.624 | 52 |
| SSD HC5 of apical chronic SMCVs | -0.691 | -0.940, -0.306 | 16 |
| median acute, species-resolved | +0.115 | -0.240, +0.424 | 41 |

The mechanism is in `figures/fig_c7.png` panel C: the number of species surviving
the screen rises with the number of records (rho = +0.76), so a better-studied
substance gets a longer left tail on its SSD and a lower 5th percentile. Every
criterion anchored on an extreme of the available data inherits this. Only the
central-tendency statistic is effort-stable, and it is not protective.
