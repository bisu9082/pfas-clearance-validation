import pandas as pd, numpy as np
rec=pd.read_csv("records_with_species.csv.gz")
tids=set(rec.test_id.astype(str))
cols=["result_id","test_id","obs_duration_mean","obs_duration_unit","endpoint",
      "effect","measurement","conc1_type","conc1_mean_op","conc1_mean","conc1_unit"]
keep=[]
for ch in pd.read_csv("/home/claude/ecotox_raw/ecotox_ascii_09_15_2026/results.txt",
        sep="|",usecols=cols,dtype=str,chunksize=400000,on_bad_lines="skip",quoting=3):
    keep.append(ch[ch.test_id.isin(tids)])
R=pd.concat(keep,ignore_index=True)
R.to_csv("results_enriched.csv.gz",index=False,compression="gzip")
print("rows",len(R))
