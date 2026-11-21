#!/usr/bin/env python3
"""Attach species and ECOTOX taxonomic group to the deposited 15,991 PFAS records."""
import pandas as pd, numpy as np
RAW="/home/claude/ecotox_raw/ecotox_ascii_09_15_2026"
REC="/home/claude/pfasv32/canonical_results/ecotox_pfas_records.csv.gz"

rec = pd.read_csv(REC)
print("deposited records      :", len(rec))

tests = pd.read_csv(f"{RAW}/tests.txt", sep="|", usecols=["test_id","species_number"],
                    dtype=str, low_memory=False, on_bad_lines="skip")
print("tests.txt rows         :", len(tests))

sp = pd.read_csv(f"{RAW}/validation/species.txt", sep="|", dtype=str,
                 usecols=["species_number","latin_name","common_name",
                          "kingdom","phylum_division","class","tax_order","ecotox_group"],
                 low_memory=False, on_bad_lines="skip")
print("species.txt rows       :", len(sp))

rec["test_id"] = rec.test_id.astype(str).str.strip()
tests["test_id"] = tests.test_id.astype(str).str.strip()
tests["species_number"] = tests.species_number.astype(str).str.strip()
sp["species_number"] = sp.species_number.astype(str).str.strip()

m = rec.merge(tests.drop_duplicates("test_id"), on="test_id", how="left")
print("after test join        : %d   species_number missing: %d"
      % (len(m), m.species_number.isna().sum()))
m = m.merge(sp.drop_duplicates("species_number"), on="species_number", how="left")
print("after species join     : %d   latin_name missing:     %d"
      % (len(m), m.latin_name.isna().sum()))

print("\ndistinct species in the PFAS set:", m.latin_name.nunique())
print("\nECOTOX group distribution (records):")
print(m.ecotox_group.fillna("(missing)").value_counts().head(15).to_string())

m.to_csv("/home/claude/wq/records_with_species.csv.gz", index=False)
print("\nwrote /home/claude/wq/records_with_species.csv.gz")
