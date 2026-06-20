#!/usr/bin/env python3
"""
GNN (D-MPNN, Chemprop v2.2.x) external validation for the ticpool PFAS study.
Run on a GPU machine:  pip install chemprop ; python run_gnn_chemprop.py
Trains a directed message-passing NN on non-PFAS broad EnviroTox chemicals and
evaluates on held-out measured PFAS, for the most-sensitive and acute endpoints.
Inputs: data/gnn_input_broad.csv, data/gnn_input_pfas.csv
"""
import subprocess, tempfile, glob
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import r2_score

broad = pd.read_csv("data/gnn_input_broad.csv")
pfas  = pd.read_csv("data/gnn_input_pfas.csv")

def run_endpoint(target):
    tr = broad[(~broad.is_pfas) & broad[target].notna()][["smiles", target]].dropna()
    te = pfas[pfas[target].notna()][["smiles", target]].dropna()
    wd = tempfile.mkdtemp()
    tr_p, te_p, md, pp = f"{wd}/train.csv", f"{wd}/test.csv", f"{wd}/model", f"{wd}/preds.csv"
    tr.to_csv(tr_p, index=False); te.to_csv(te_p, index=False)
    subprocess.run(["chemprop","train","-i",tr_p,"-t","regression","-s","smiles",
        "--target-columns",target,"-o",md,"-n","0","--epochs","50",
        "--accelerator","gpu","--devices","1"], check=True)
    ck = (glob.glob(f"{md}/**/best.pt", recursive=True)
          or glob.glob(f"{md}/**/*.ckpt", recursive=True)
          or glob.glob(f"{md}/**/*.pt", recursive=True))[0]
    subprocess.run(["chemprop","predict","-i",te_p,"-s","smiles",
        "--model-paths",ck,"-o",pp], check=True)
    pr = pd.read_csv(pp)
    pcol = [c for c in pr.columns if c != "smiles"][-1]
    y, yhat = te[target].values, pr[pcol].values
    return r2_score(y, yhat), spearmanr(y, yhat).correlation, len(te)

for ep in ["most_sens", "acute"]:
    r2, rho, n = run_endpoint(ep)
    print(f"GNN (D-MPNN)  {ep:10s}  external PFAS (n={n}):  R2={r2:+.2f}  rho={rho:+.2f}")
# Result (RTX 3090, 2026-06-20): most_sens R2=-2.04 rho=+0.34 ; acute R2=-1.97 rho=+0.08
