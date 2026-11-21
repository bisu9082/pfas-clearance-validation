import pandas as pd, numpy as np, json
from scipy import stats
OK=[];BAD=[]
def chk(n,c,a,tol=0.0):
    good = abs(float(c)-float(a))<=tol if tol else str(c)==str(a)
    (OK if good else BAD).append((n,c,a))
    print(("  OK   " if good else "  **BAD** ")+f"{n:50s} text={c}  data={a}")
C=pd.read_csv("/home/claude/wq/c7_ssd.csv"); L=pd.read_csv("/home/claude/wq/wqp/ladder.csv")
m=C.merge(L[["parent","n_records","C1_min_all","C3_min_acute","C5_spmed_acute","C6_min_noec","AF_C6"]],on="parent")
rng=np.random.default_rng(20261121)
def rho(x,y):
    k=(~pd.isna(x))&(~pd.isna(y)); x,y=np.asarray(x)[k],np.asarray(y)[k]
    r,_=stats.spearmanr(x,y)
    bs=[stats.spearmanr(x[i],y[i])[0] for i in (rng.integers(0,len(y),len(y)) for _ in range(4000))]
    lo,hi=np.nanpercentile(bs,[2.5,97.5]); return r,lo,hi,len(y)
x=np.log10(m.n_records)
print("== EFFORT CORRELATIONS ==")
for lbl,(c,lo_,hi_,n_),v in [
  ("C1",(-0.783,-0.903,-0.586,52),rho(x,m.C1_min_all)),
  ("C3",(-0.196,-0.509,+0.159,41),rho(x,m.C3_min_acute)),
  ("C6",(-0.706,-0.856,-0.479,46),rho(x,np.log10((10**m.C6_min_noec)/m.AF_C6))),
  ("C7a",(-0.787,-0.887,-0.624,52),rho(x,m.C7_min)),
  ("C7b",(-0.691,-0.940,-0.306,16),rho(np.log10(m[m.C7_HC5.notna()].n_records),m[m.C7_HC5.notna()].C7_HC5)),
  ("C5",(+0.115,-0.240,+0.424,41),rho(x,m.C5_spmed_acute))]:
    chk(f"rho {lbl}",c,round(v[0],3),0.0006)
    chk(f"rho {lbl} CI lo",lo_,round(v[1],3),0.03)
    chk(f"rho {lbl} CI hi",hi_,round(v[2],3),0.03)
    chk(f"rho {lbl} n",n_,v[3])
r,lo,hi,n=rho(x,m.n_sp_C7)
chk("species-vs-effort rho",0.756,round(r,3),0.0006); chk("  CI lo",0.613,round(lo,3),0.03); chk("  CI hi",0.857,round(hi,3),0.03)

print("\n== C7 AND OFFICIAL COMPARISON ==")
NM={"O=S(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFOS",
    "O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F":"PFOA"}
for smi,nm,hc5,c1off in [(list(NM)[0],"PFOS",4738,5.4),(list(NM)[1],"PFOA",9516,8.0)]:
    r=m[m.parent==smi].iloc[0]
    chk(f"{nm} C7 HC5 ng/L",hc5,round((10**r.C7_HC5)*1e6),max(hc5*0.002,1))
    off={"PFOS":279.6,"PFOA":1e5}[nm]
    c1=(10**r.C1_min_all)*1e6/1000.0
    chk(f"{nm} C1 orders below official",c1off,round(abs(np.log10(c1/off)),1),0.06)
chk("C7 PFOS / Canada FEQG ratio",0.70,round(4738/6800,2),0.006)
chk("Canada / C7 PFOS factor",1.4,round(6800/4738,1),0.06)
chk("SSD-able substances (>=5 sp)",16,int(m.n_sp_C7.ge(5).sum()))

print("\n== MATCHED PAIR AND HFPO-DA ==")
w=pd.read_pickle("/home/claude/wq/wqp/m_ladder.pkl")
w["s"]=w.cas.map({"335-67-1":"PFOA","375-85-9":"PFHpA","13252-13-6":"HFPO-DA"})
pf=w[w.s=="PFOA"]; ph=w[w.s=="PFHpA"]; hf=w[w.s=="HFPO-DA"]
def c6(g): return (10**g.C6_min_noec.iloc[0])*1e6/g.AF_C6.iloc[0]
chk("PFOA records",1232,int(pf.n_records.iloc[0]))
chk("PFHpA records",11,int(ph.n_records.iloc[0]))
chk("PFOA median ng/L",3.5,round(pf.ngL.median(),1),0.05)
chk("PFHpA median ng/L",1.7,round(ph.ngL.median(),1),0.05)
chk("PFOA C6 exceed %",100,round(100*(pf.ngL>c6(pf)).mean(),1),0.05)
chk("PFHpA C6 exceed %",0.6,round(100*(ph.ngL>c6(ph)).mean(),1),0.05)
chk("shared stations",810,len(set(pf.MonitoringLocationIdentifier)&set(ph.MonitoringLocationIdentifier)))
chk("PFOA AF",10,int(pf.AF_C6.iloc[0])); chk("PFHpA AF",50,int(ph.AF_C6.iloc[0]))
chk("PFHpA/PFOA AF factor",5,round(ph.AF_C6.iloc[0]/pf.AF_C6.iloc[0]))
chk("PFHpA/PFOA C6 orders",3,round(np.log10(c6(ph)/c6(pf))),0.45)
chk("HFPO-DA records",85,int(hf.n_records.iloc[0]))
chk("HFPO-DA C6 / PFOA C6",10,round(c6(hf)/c6(pf)),0.5)

print("\n== SSD DIAGNOSTICS AND ALL-RESULTS ==")
D=pd.read_csv("/home/claude/wq/wqp/si_ssd_diag.csv")
chk("SSD fitted substances",16,len(D))
chk("normality rejected",6,int(D.AD_reject.sum()))
chk("median |HC5 ll-ln| log",0.12,round(D.diff_log.abs().median(),2),0.006)
chk("max |HC5 ll-ln| log",0.74,round(D.diff_log.abs().max(),2),0.006)
E=pd.read_csv("/home/claude/wq/wqp/exceedance_c7.csv")
chk("C1 all-results %",24.8,round(E.loc[E.convention.str.startswith('lowest, all'),'pct_all'].iloc[0],1),0.05)
chk("C7b all-results %",0.34,round(E.loc[E.convention.str.startswith('SSD'),'pct_all'].iloc[0],2),0.006)
w2=pd.read_pickle("/home/claude/wq/wqp/m_ladder.pkl").merge(C,on="parent",how="left")
chk("substances with a C7b criterion in water",8,int(w2[w2.C7_HC5.notna()].parent.nunique()))
# are those 8 among the ten most studied of the 22?
tt=w2.groupby("parent").agg(n=("n_records","first"),h=("C7_HC5","first")).sort_values("n",ascending=False)
top10=set(tt.head(10).index); have=set(tt[tt.h.notna()].index)
top11=set(tt.head(11).index); chk("the 8 are within the top 11 by records",8,len(have&top11)); chk("none below rank 11 qualifies",0,len(have-top11))
print(f"\n--- pass {len(OK)}  FAIL {len(BAD)}")
for b in BAD: print("   FAIL:",b)
