# pfas-clearance-validation

Data and code for the manuscript:

**Can Structure-Based Models Clear Untested PFAS? External Validation Reveals
Endpoint-, Representation-, and Distribution-Dependent Reliability**
Jeongyun Kim*, Jin Yoo*, Ku Kang (\*equal contribution). Submitted to
*Environmental Science & Technology*.

## Summary

A structure-based aquatic-toxicity model that is moderately accurate on general
chemicals (5×5 repeated CV R² = 0.39 ± 0.02) **fails on fully unseen PFAS**
(external R² = −0.58, Spearman ρ = −0.07). Standard (marginal) conformal
intervals under-cover PFAS (0.659 vs nominal 0.90); PFAS-conditional (Mondrian)
calibration restores held-out coverage (0.922) but the intervals become so wide
(median ≈10.5 log₁₀ mg/L) that **no untested PFAS can be cleared** as low-hazard
at 90% confidence. The failure is governed by the endpoint and the molecular
representation, not by the algorithm: a clean acute endpoint and physicochemical
descriptors recover *rank* prediction (ρ up to +0.48 / +0.42), though not
clearance-grade accuracy. The message is that structure-based scores support
**prioritization, not low-hazard clearance**.

## Repository structure

```
manuscript/            LaTeX sources, figures, and compiled PDFs
  main.tex, SI.tex, cover_letter.tex, myref.bib
  figures/fig_01..04.png
  main.pdf, SI.pdf, cover_letter.pdf
canonical_results/     All authoritative numbers behind the manuscript
  canonical_metrics.json          headline metrics
  broad_nonpfas_structure_locked.csv   3,788 unique non-PFAS training structures
  measured_pfas_locked.csv             138 measured PFAS (external test)
  untested_pfas_structure_locked.csv   18,851 unique untested inventory structures
  conformal_repeated_split.csv / conformal_summary.csv   40-repeat coverage/width
  endpoint_*_predictions.csv, endpoint_representation_results.csv,
  endpoint_paired_same_pfas.csv, endpoint_paired_summary.json
  external_pfas_predictions.csv, inventory_measured_overlap.csv
  untested_pfas_(marginal_and_)conditional90.csv
  algorithm_comparison.csv, ad_threshold_sensitivity.csv,
  deployment_loo_stability.csv         SI robustness tables
scripts/               Reproduction scripts
  reproduce_robustness.py    regenerates algorithm / AD-threshold / LOO tables
  run_gnn_chemprop.py        optional D-MPNN supplementary check (GPU)
  gnn_input_broad.csv, gnn_input_pfas.csv
```

## Data sources (not redistributed)

- **EnviroTox** v2.0 — curated freshwater aquatic-toxicity records.
- **U.S. EPA ECOTOX** Knowledgebase — measured PFAS aquatic effects.
- **U.S. EPA CompTox Chemicals Dashboard, PFASSTRUCTV6** — screening inventory.

Raw third-party databases are not redistributed here. `canonical_results/`
contains the processed, structure-locked tables that the manuscript uses; the
scripts document the exact, reproducible processing rules.

## Reproduce

```bash
pip install numpy pandas scikit-learn scipy rdkit
python scripts/reproduce_robustness.py        # algorithm / AD / LOO robustness
# optional learned-representation check (needs a GPU + chemprop):
pip install chemprop
python scripts/run_gnn_chemprop.py
```

Random seed 20260620; ECFP4 = 2048-bit Morgan radius 2; random forest 300 trees,
max_features = sqrt.

## Citation

Kim J.*, Yoo J.*, Kang K. *Can Structure-Based Models Clear Untested PFAS? …*
(2026). (\*equal contribution.) Update with DOI on acceptance.

## License

Code: MIT (see `LICENSE`). Processed result tables may be reused with attribution.
