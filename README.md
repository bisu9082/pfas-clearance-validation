# pfas-clearance-validation

Data and code for:

**A testing-effort artefact in the most-sensitive aquatic endpoint: structure-based
PFAS screening recovers predictive validity but still cannot support clearance**
Jeongyun Kim\*, Jin Yoo\*, Ku Kang (\*equal contribution).
Submitted to *Journal of Hazardous Materials*.

> Release tags: `v3.0` is this work. `v2.0` and `v1.0` are earlier states of the
> project, built on a curated subset of the aquatic record and reaching a
> different conclusion by a different route; they are kept for provenance.

## What this paper reports

Regulatory hazard screening summarises a substance by the **lowest** effect
concentration reported for it. That makes the response variable depend on how
often the substance has been tested, because the minimum of a sample falls as
the sample grows.

Across 103 PFAS with freshwater data in ECOTOX, the protective minimum correlates
with the log number of available records at **ρ = −0.73** (p = 1×10⁻¹⁸).
Substances with >300 records have a median protective value **6.4 log₁₀ units
below** those with ≤3. Subsampling every substance to equal depth halves the
association (ρ ≈ −0.29), so roughly half is a pure extreme-order-statistic
artefact and roughly half is genuine reverse causation. The same holds across all
4,757 chemicals with freshwater acute data: the acute **minimum** carries ρ = −0.46
against record count, the acute **median** ρ = −0.16.

A structure-based model cannot predict testing history. Holding one published
model fixed and changing only how the same records are summarised:

| external target | R² | ρ | effort confound |
|---|---|---|---|
| minimum (protective convention) | −0.745 | −0.015 | −0.734 |
| median | −0.147 | +0.076 | −0.442 |
| **acute median** | **+0.216** | **+0.424** | **−0.036** |

Rebuilding training and external sets consistently on the acute median raises
broad CV to R² = 0.443 ± 0.029, restores external validity on 58 structure-locked
PFAS (R² = +0.367; ρ = +0.460, 95% CI [+0.22, +0.65]), and halves the
PFAS-conditional 90% interval from 10.53 to 4.46 log₁₀ mg/L.

**Clearance still fails.** The acute median is 2.69 log₁₀ units less protective
than the minimum on the same substances, so the 0.1 mg/L threshold is not
transferable between them. At a stringency-matched threshold, 0 of 16,070
in-domain structures clear; at the GHS acute Category 1 boundary (1 mg/L), 1.9% do.

The published conclusion that structure-based screening cannot clear untested
PFAS therefore stands — but not because such models cannot predict PFAS aquatic
toxicity. They can, on an endpoint built without an effort confound. What blocks
clearance is interval width relative to the decision threshold, and at protective
stringency the point predictions themselves.

## Repository structure

```
scripts/            extraction, diagnostics, rebuild, figures
canonical_results/  per-substance summaries, metrics, inventory predictions
manuscript/         LaTeX sources, figures, compiled PDFs
```

## Reproducing

Source data are public and not redistributed here:
ECOTOX ASCII release `ecotox_ascii_09_15_2026` (https://gaftp.epa.gov/ecotox/),
EPA PFASSTRUCTV6, and PubChem for CAS→SMILES. Seed 20260620 throughout.

```bash
python scripts/extract.py        <ecotox_raw_dir> <pfas_cas_dir> canonical_results
python scripts/effort_test.py    <ecotox_raw_dir> <pfas_cas_dir> canonical_results
python scripts/build_broad.py    <ecotox_raw_dir> canonical_results
python scripts/resolve_smiles.py canonical_results/ecotox_acute_all.csv canonical_results
python scripts/endpoint_fix.py   <broad.csv> canonical_results/ecotox_pfas_records.csv \
                                 canonical_results/ecotox_pfas_measured.csv canonical_results
python scripts/final_pipeline.py canonical_results/ecotox_acute_smiles.csv \
                                 canonical_results/ecotox_pfas_records.csv \
                                 canonical_results/ecotox_pfas_measured.csv \
                                 <inventory.csv> canonical_results
```

`resolve_smiles.py` takes ~35 min at PubChem's 5 req/s limit and caches to
`cas_smiles_cache.json`, so it is interruptible and resumable. 4,311 of 4,757
chemicals (90.6%) resolve.

Requirements: `rdkit`, `scikit-learn`, `pandas`, `numpy`, `scipy`, `matplotlib`,
`requests`.

## Choices that change the data, stated explicitly

- **Censored records excluded.** Values carrying `>`, `<`, `>=`, `<=` are not
  used as point observations. Left-censored entries have maximum leverage on a
  minimum, which is the quantity under examination.
- **Desalting keeps the fluorine-bearing fragment, not the largest.** For
  tetrabutylammonium perfluoro-1-butanesulfonate the counterion is the larger
  fragment, so the usual largest-fragment rule discards the analyte.
- **Fragments are charge-neutralised** before canonicalisation. Without it,
  perfluorooctanesulfonate persists as seven separate "substances" whose assigned
  values span 7.93 log₁₀ units, and they can fall on both sides of a calibration
  split.
- **Subclass results are reported with bootstrap intervals.** At the available
  per-class sizes (6–14 on the acute-median set) 80% power detects only
  |ρ| ≥ 0.69. The manuscript reports that limit rather than the point estimates
  alone; resolving ρ = 0.3 across seven OECD classes would need ~85 substances
  per class.

## Licence

See `LICENSE`.
