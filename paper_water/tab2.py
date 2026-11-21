import pandas as pd, numpy as np
m=pd.read_pickle("/home/claude/wq/wqp/m_ladder.pkl")
C=pd.read_csv("/home/claude/wq/c7_ssd.csv"); m=m.merge(C,on="parent",how="left")
NAME={"375-22-4":"PFBA","335-67-1":"PFOA","1763-23-1":"PFOS","307-24-4":"PFHxA",
"375-73-5":"PFBS","2706-90-3":"PFPeA","375-85-9":"PFHpA","355-46-4":"PFHxS",
"375-95-1":"PFNA","754-91-6":"FOSA","335-76-2":"PFDA","27619-97-2":"6:2 FTS",
"307-55-1":"PFDoA","2058-94-8":"PFUnA","72629-94-8":"PFTrDA","30334-69-1":"FBSA",
"812-70-4":"7:3 FTCA","914637-49-3":"5:3 FTCA","41997-13-1":"FHxSA",
"24448-09-7":"N-MeFOSE","13252-13-6":"HFPO-DA","4151-50-2":"N-EtFOSA"}
m["s"]=m.cas.map(NAME).fillna(m.cas)
def pn(c,af): return (10.0**c)*1e6/af
rows=[]
for s,g in m.groupby("s"):
    d=g.iloc[0]
    rows.append(dict(s=s,n=len(g),sites=g.MonitoringLocationIdentifier.nunique(),
      med=g.ngL.median(),recs=int(d.n_records),
      C1=pn(d.C1_min_all,1000.0),
      C6=pn(d.C6_min_noec,d.AF_C6) if pd.notna(d.C6_min_noec) and pd.notna(d.AF_C6) else np.nan,
      C7=pn(d.C7_HC5,1.0) if pd.notna(d.C7_HC5) else np.nan,
      e1=100*(g.ngL>pn(d.C1_min_all,1000.0)).mean(),
      e6=100*(g.ngL>pn(d.C6_min_noec,d.AF_C6)).mean() if pd.notna(d.C6_min_noec) and pd.notna(d.AF_C6) else np.nan,
      e7=100*(g.ngL>pn(d.C7_HC5,1.0)).mean() if pd.notna(d.C7_HC5) else np.nan))
t=pd.DataFrame(rows).sort_values("n",ascending=False)
def f(v,d=3):
    if pd.isna(v): return "---"
    if v>=1e4: return f"\\num{{{v:.2e}}}"
    if v>=100: return f"{v:,.0f}"
    if v>=1:   return f"{v:.3g}"
    return f"{v:.3g}"
L=[]
for _,r in t.iterrows():
    L.append(f"{r.s} & {int(r.n):,} & {int(r.sites)} & {r.med:.2f} & {int(r.recs):,} & "
             f"{f(r.C1)} & {f(r.C6)} & {f(r.C7)} & "
             f"{'---' if pd.isna(r.e1) else f'{r.e1:.2f}'} & "
             f"{'---' if pd.isna(r.e6) else f'{r.e6:.2f}'} & "
             f"{'---' if pd.isna(r.e7) else f'{r.e7:.2f}'} \\\\")
tex=r"""\begin{table}[p]
\centering
\caption{Measured concentrations, derived criteria and exceedance by substance.
Measurements are Water Quality Portal results for United States stream and lake
stations, 2016--2026. Criteria are in \si{\nano\gram\per\litre}: C1 the lowest
value over all endpoints with an assessment factor of 1000, C6 the lowest
chronic \textsf{NOEC} with an EU Technical Guidance factor, C7b the fifth
percentile of a log-normal species sensitivity distribution fitted to apical
chronic species mean values. A dash means the convention cannot be applied for
want of qualifying data.}
\label{tab:persub}
\scriptsize
\setlength{\tabcolsep}{3.2pt}
\begin{tabular}{@{}l r r r r r r r r r r@{}}
\toprule
 & \multicolumn{3}{c}{measured} & tox & \multicolumn{3}{c}{criterion (\si{\nano\gram\per\litre})}
 & \multicolumn{3}{c}{\% samples $>$ criterion} \\
\cmidrule(lr){2-4}\cmidrule(lr){6-8}\cmidrule(lr){9-11}
Substance & $n$ & sites & median & recs & C1 & C6 & C7b & C1 & C6 & C7b \\
\midrule
""" + "\n".join(L) + r"""
\bottomrule
\end{tabular}
\end{table}"""
open("tab_persub.tex","w").write(tex)
print(t[["s","n","recs","C1","C6","C7","e1","e6","e7"]].to_string(index=False,float_format=lambda x:f"{x:,.4g}"))
