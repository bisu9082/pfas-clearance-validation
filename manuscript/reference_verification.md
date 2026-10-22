# Reference verification report

Bibliography: `myref.bib` — 72 entries, all 72 cited in the manuscript or SI.

Every DOI was resolved against the Crossref REST API and compared on **title, author surnames, first author, author count, year, journal, volume and page range**. An entry passes only if all available fields agree.

| | count |
|---|---|
| DOI resolved and all fields verified | 62 |
| No DOI (software, reports, regulatory documents — none is a journal article) | 10 |
| Mismatched or unresolvable | 0 |

## Entries without a DOI

| key | what it is |
|---|---|
| `epacomptox2026` | EPA CompTox Dashboard (web resource) |
| `oecd2021pfas` | OECD report |
| `epa2021roadmap` | US EPA report |
| `settles2009active` | Univ. Wisconsin tech. report |
| `sklearn` | JMLR (no DOI assigned) |
| `rdkit` | open-source software |
| `shafer2008tutorial` | JMLR (no DOI assigned) |
| `echa2008r10` | ECHA guidance |
| `eu2013directive` | EU legislation |
| `epa2024mcl` | US federal rule |

## Corrections made during this audit

| key | problem found | resolution |
|---|---|---|
| `saouter2017usetox` | DOI `10.1002/etc.3921` resolved to a different paper; title and author list were also wrong | replaced with `10.1002/etc.3889`, entry regenerated from Crossref |
| `vandenheuvel2006ppar` | DOI `10.1093/toxsci/kfj128` resolved to a different paper | corrected to `10.1093/toxsci/kfl014` (title, authors, volume and pages were already right) |
| `vandijk2023dataselection` | correct title, but wrong authors, volume and pages — the paper is by Aggarwal et al. | rekeyed to `aggarwal2024cfpm`, regenerated from Crossref |
| `mceachran2018msready` | three author surnames did not belong to this paper | regenerated from Crossref |
| `mansouri2018opera` | one author surname did not belong to this paper | regenerated from Crossref |
| `olker2022ecotox` | author count did not match the record | regenerated from Crossref |
| `liang2024datagaps` | year was 2024; the article is 2025 | rekeyed to `liang2025datagaps`, regenerated from Crossref |

## Unverifiable quotation removed

The Introduction previously quoted Saouter et al. (2017) verbatim. The full text could not be retrieved, so the quotation was removed and replaced with a statement supported by the paper's abstract, which was verified word for word. The two quotations from Schur et al. (2023) were checked against the published full text and are verbatim and in context.
