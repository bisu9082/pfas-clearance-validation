import pandas as pd, numpy as np
def esc(t): return str(t).replace("&","\\&").replace("_","\\_").replace("%","\\%").replace("***","\\,***\\,")
# ---- S1 mapping
cm=pd.read_csv("/home/claude/wq/wqp/char_map.csv").sort_values("n",ascending=False)
L=[]
for _,r in cm.iterrows():
    L.append(f"{esc(r['char'])[:64]} & {r['cas']} & {int(r['n']):,} & "
             f"{'yes' if r['in_ecotox'] else 'no'} & {('$\\mathrm{'+str(r['formula'])+'}$') if isinstance(r['formula'],str) and r['formula']==r['formula'] else '---'} \\\\")
open("si_tab1.tex","w").write(
"\\begin{longtable}{@{}p{66mm} l r c l@{}}\n\\caption{Water Quality Portal characteristic names with quantified data, the CAS "
"registry number assigned to each, the number of detections, whether that CAS appears in the ECOTOX substance file, and the "
"molecular formula of the parent structure retrieved from ECOTOX. The formula is shown so that the name-to-structure "
"assignment can be checked independently of the mapping table.}\\label{tab:s1}\\\\\n\\toprule\n"
"Characteristic name & CAS & $n$ & in ECOTOX & parent formula \\\\\n\\midrule\n\\endfirsthead\n"
"\\toprule\nCharacteristic name & CAS & $n$ & in ECOTOX & parent formula \\\\\n\\midrule\n\\endhead\n"
+ "\n".join(L) + "\n\\bottomrule\n\\end{longtable}\n")
print("S1 rows:",len(L))

# ---- S2 SSD diagnostics
D=pd.read_csv("/home/claude/wq/wqp/si_ssd_diag.csv")
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
RDLogger.DisableLog('rdApp.*')
NAME={"O=S(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFOS",
"O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFOA",
"O=C(O)C(F)(OC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F":"HFPO-DA",
"O=S(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFBS","O=C(O)C(F)(F)C(F)(F)C(F)(F)F":"PFBA",
"O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFNA",
"O=S(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFHxS",
"O=S(=O)(O)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"6:2 FTS"}
def lb(p):
    if p in NAME: return NAME[p]
    m=Chem.MolFromSmiles(p); return "$\\mathrm{"+rdMolDescriptors.CalcMolFormula(m)+"}$" if m else "?"
D["lab"]=[lb(p) for p in D.parent]
D=D.sort_values("n",ascending=False)
rows=[f"{r.lab} & {int(r.n)} & {r.mu:.2f} & {r.sd:.2f} & {r.AD:.2f} & {r.AD_crit5:.2f} & "
      f"{'reject' if r.AD_reject else 'ok'} & {r.shapiro_p:.3f} & {r.HC5_lognormal:.2f} & "
      f"{r.HC5_loglogistic:.2f} & {r.diff_log:+.2f} \\\\" for _,r in D.iterrows()]
open("si_tab2.tex","w").write(
"\\begin{table}[htbp]\\centering\\caption{Fit diagnostics for every species sensitivity distribution with at least five "
"species. HC5 values are $\\log_{10}$ of the fifth percentile in \\si{\\milli\\gram\\per\\litre}. Anderson--Darling is "
"compared with its 5\\,\\% critical value.}\\label{tab:s2}\\footnotesize\\setlength{\\tabcolsep}{3.4pt}\n"
"\\begin{tabular}{@{}l r r r r r l r r r r@{}}\\toprule\n"
"Substance & $n$ & mean & s.d. & AD & crit & normality & Shapiro $p$ & HC5 ln & HC5 ll & $\\Delta$ \\\\\\midrule\n"
+ "\n".join(rows) + "\n\\bottomrule\\end{tabular}\\end{table}\n")
print("S2 rows:",len(rows))

# ---- S4 all-results exceedance
E=pd.read_csv("/home/claude/wq/wqp/exceedance_c7.csv")
rows=[f"{esc(r.convention)} & {int(r.n_det):,} & {int(r.n_all):,} & {r.pct_samples:.2f} & {r.pct_all:.2f} \\\\"
      for _,r in E.iterrows()]
open("si_tab4.tex","w").write(
"\\begin{table}[htbp]\\centering\\caption{Exceedance computed over detections only, as reported in the main text, and over "
"all returned Water Quality Portal results for the same substances, with non-detections counted as not exceeding. The "
"latter is a lower bound: reporting limits are not carried in the result profile used, so a non-detection cannot be placed "
"relative to a criterion that lies below the limit.}\\label{tab:s4}\\footnotesize\n"
"\\begin{tabular}{@{}l r r r r@{}}\\toprule\nConvention & detections & all results & \\% of detections & \\% of all \\\\\\midrule\n"
+ "\n".join(rows) + "\n\\bottomrule\\end{tabular}\\end{table}\n")
print("S4 rows:",len(rows))
