"""Rebuild the measured-PFAS aquatic toxicity set from the full ECOTOX dump.

Protocol mirrors the manuscript: freshwater, aquatic organisms, mass-concentration
units, acute LC50/EC50 plus chronic NOEC/LOEC/MATC, protective (minimum) target.
Censored records (>, <, >=, <=) are kept separately rather than silently used.
"""
import sys, re
import pandas as pd, numpy as np

RAW, REPO, OUT = sys.argv[1], sys.argv[2], sys.argv[3]

# ---- PFAS CAS universe (PFASSTRUCTV6) ----------------------------------
inv = pd.read_csv(f"{REPO}/untested_pfas_structure_locked.csv",
                  usecols=["CASRN", "cas_norm", "canonical_smiles", "PREFERRED NAME"])
mea = pd.read_csv(f"{REPO}/measured_pfas_locked.csv", usecols=["cas", "name", "canonical_smiles"])
mea["cas_norm"] = mea.cas.astype(str).str.replace(r"\D", "", regex=True)
pfas = pd.concat([
    inv.rename(columns={"PREFERRED NAME": "name"})[["cas_norm", "name", "canonical_smiles"]],
    mea[["cas_norm", "name", "canonical_smiles"]]]).drop_duplicates("cas_norm")
pfas["cas_norm"] = pfas.cas_norm.astype(str)
print(f"PFAS CAS universe: {len(pfas)}", flush=True)

# ---- tests: freshwater, aquatic, PFAS ----------------------------------
tcols = ["test_id", "test_cas", "species_number", "organism_habitat", "media_type"]
tests = pd.read_csv(f"{RAW}/tests.txt", sep="|", usecols=tcols, dtype=str,
                    low_memory=False, on_bad_lines="skip")
tests["test_cas"] = tests.test_cas.astype(str).str.strip()
t = tests[tests.test_cas.isin(set(pfas.cas_norm))]
print(f"tests on PFAS CAS: {len(t)}  (unique CAS {t.test_cas.nunique()})", flush=True)
t = t[t.organism_habitat.eq("Water") &
      t.media_type.fillna("").str.upper().str.startswith("FW")]
print(f"  freshwater + aquatic: {len(t)}  (unique CAS {t.test_cas.nunique()})", flush=True)

# ---- results -----------------------------------------------------------
rcols = ["result_id", "test_id", "endpoint", "effect",
         "conc1_mean_op", "conc1_mean", "conc1_unit"]
res = pd.read_csv(f"{RAW}/results.txt", sep="|", usecols=rcols, dtype=str,
                  low_memory=False, on_bad_lines="skip")
m = res.merge(t, on="test_id", how="inner")
print(f"results joined: {len(m)}", flush=True)

ACUTE = {"LC50", "EC50", "IC50"}
CHRON = {"NOEC", "LOEC", "MATC", "NOEL", "LOEL", "MATC"}
m["ep"] = m.endpoint.astype(str).str.strip().str.upper().str.replace("*", "", regex=False)
m = m[m.ep.isin(ACUTE | CHRON)]
print(f"  with a usable endpoint: {len(m)}", flush=True)

UNIT = {  # -> multiplier to mg/L
    "mg/L": 1.0, "ug/L": 1e-3, "ng/L": 1e-6, "g/L": 1e3, "pg/L": 1e-9,
    "ppm": 1.0, "ppb": 1e-3, "ppt": 1e-6, "mg/dm3": 1.0, "ug/dm3": 1e-3,
    "AI mg/L": 1.0, "AI ug/L": 1e-3, "AE mg/L": 1.0, "AE ug/L": 1e-3,
    "mg/l": 1.0, "ug/l": 1e-3, "ng/l": 1e-6,
}
m["u"] = m.conc1_unit.astype(str).str.strip()
m["mult"] = m.u.map(UNIT)
unknown = m[m.mult.isna()].u.value_counts().head(8)
m = m[m.mult.notna()]
m["val"] = pd.to_numeric(m.conc1_mean, errors="coerce")
m = m[m.val.notna() & (m.val > 0)]
m["mgL"] = m.val * m.mult
print(f"  with mass-concentration units and a numeric value: {len(m)}", flush=True)
print(f"  (most common rejected units: {dict(unknown)})", flush=True)

m["censored"] = m.conc1_mean_op.astype(str).str.strip().isin([">", "<", ">=", "<=", "~"])
print(f"  censored: {m.censored.sum()}   uncensored: {(~m.censored).sum()}", flush=True)

unc = m[~m.censored].copy()
unc["logTox"] = np.log10(unc.mgL)

agg = (unc.groupby("test_cas")
          .agg(logTox=("logTox", "min"), n_records=("logTox", "size"),
               n_acute=("ep", lambda s: s.isin(ACUTE).sum()),
               logTox_median=("logTox", "median"),
               logTox_max=("logTox", "max"))
          .reset_index())
agg["spread"] = agg.logTox_max - agg.logTox
out = agg.merge(pfas, left_on="test_cas", right_on="cas_norm", how="left")
out.to_csv(f"{OUT}/ecotox_pfas_measured.csv", index=False)

print(f"\n=== UNIQUE PFAS WITH FRESHWATER AQUATIC DATA: {len(out)} ===")
print(f"   (manuscript used 138 rows / 119 unique substances)")
print(f"   with >=3 records: {(out.n_records>=3).sum()}   with acute data: {(out.n_acute>0).sum()}")
print(f"   median records per substance: {out.n_records.median():.0f}")
print("\n--- widest within-substance spread (log10 mg/L) ---")
print(out.nlargest(8, "spread")[["name","n_records","logTox","logTox_max","spread"]]
         .to_string(index=False, float_format=lambda x: f"{x:7.2f}", max_colwidth=46))
