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
