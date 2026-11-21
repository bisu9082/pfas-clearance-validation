import pandas as pd, numpy as np, re, json
from scipy import stats
OK=[];BAD=[]
def chk(name,claim,actual,tol=0.0):
    try:
        good = abs(float(claim)-float(actual))<=tol if tol else (str(claim)==str(actual))
    except Exception:
        good = str(claim)==str(actual)
    (OK if good else BAD).append((name,claim,actual))
    print(("  OK   " if good else "  **BAD** ")+f"{name:52s} text={claim}  data={actual}")

R=pd.read_csv("/home/claude/wq/results_parent.csv.gz")
DEP={"EC50","IC50","LC50","LOEC","LOEL","MATC","NOEC","NOEL"}
G=R[R.ep.isin(DEP)].copy()
d=G.drop_duplicates(subset=["test_id","ep","conc1_mean","conc1_unit"])
AP={"MOR","GRO","REP","DVP","POP"}
print("\n== TOXICITY SIDE ==")
chk("curated records",4744,len(d))
chk("parent substances",87,d.parent.nunique())
chk("distinct species",151,d.latin_name.nunique())
chk("NOEC records",1918,int((d.ep=="NOEC").sum()))
chk("LOEC records",1811,int((d.ep=="LOEC").sum()))
chk("LOEC share %",38.2,round(100*(d.ep=="LOEC").mean(),1),0.05)
chk("acute median-effect records",953,int(d.ep.isin({"EC50","LC50","IC50"}).sum()))
for e,v in [("GEN",3501),("BCM",1613),("ENZ",867),("MOR",1011),("GRO",1002),
            ("POP",564),("REP",520),("DVP",311)]:
    chk(f"effect {e} rows",v,int((G.effect==e).sum()))
chk("apical share of rows %",29.8,round(100*G.effect.isin(AP).mean(),1),0.05)
g=G.groupby("parent").agg(n_all=("logTox","size"),ap=("effect",lambda s:s.isin(AP).sum()))
g["f"]=g.ap/g.n_all; g=g[g.n_all>=10]
chk("apical frac, <50 recs (median)",0.607,round(g[g.n_all<50].f.median(),3),0.0006)
chk("apical frac, >=200 recs (median)",0.216,round(g[g.n_all>=200].f.median(),3),0.0006)
n=d.groupby("parent").size().sort_values(ascending=False).values
chk("top-2 share of records %",59,round(100*n[:2].sum()/n.sum()),0)
chk("median records per substance",8,int(np.median(n)))

print("\n== WATER SIDE ==")
conc=pd.read_pickle("/home/claude/wq/wqp/conc.pkl")
raw=pd.read_pickle("/home/claude/wq/wqp/all_raw.pkl")
m=pd.read_pickle("/home/claude/wq/wqp/m_ladder.pkl")
cm=pd.read_csv("/home/claude/wq/wqp/char_map.csv")
chk("WQP results returned",117034,len(raw))
chk("detections",22255,len(conc))
chk("detection rate %",19.0,round(100*len(conc)/len(raw),1),0.05)
chk("non-detect flagged",89857,int(raw.ResultDetectionConditionText.notna().sum()))
chk("characteristic names with detections",66,len(cm))
chk("distinct CAS",40,cm.cas.nunique())
chk("CAS absent from ECOTOX",18,cm[~cm.in_ecotox].cas.nunique())
chk("measurements joined",20682,len(m))
chk("stations",1180,m.MonitoringLocationIdentifier.nunique())
chk("substances",22,m.parent.nunique())
chk("first sample date","2016-09-20",str(m.ActivityStartDate.min()))
chk("last sample date","2026-06-25",str(m.ActivityStartDate.max()))

print("\n== TABLE 1 ==")
E=pd.read_csv("/home/claude/wq/wqp/exceedance_c7.csv")
want={"lowest, all endpoints (as published)":(0.05,76.5,99.3),
      "lowest, LOECs removed":(0.05,76.0,99.3),
      "lowest acute":(2740,22.7,87.3),
      "lowest chronic NOEC, EU-TGD AF":(1,55.5,97.7),
      "lowest apical chronic SMCV, AF 10":(430,8.5,55.4),
      "SSD HC5 of apical chronic SMCVs":(4738,0.78,3.6),
      "median acute, record-resolved":(1.46e5,0.11,1.2),
      "median acute, species-resolved":(2.97e5,0.11,1.2)}
for _,r in E.iterrows():
    w=want[r.convention]
    chk(f"T1 {r.convention[:30]} PNEC",w[0],round(r.med_PNEC_ngL,2 if w[0]<10 else 0),max(abs(w[0])*0.005,0.005))
    chk(f"T1 {r.convention[:30]} samples",w[1],round(r.pct_samples,2),0.05)
    chk(f"T1 {r.convention[:30]} sites",w[2],round(r.pct_sites,2),0.05)
json.dump({"ok":len(OK),"bad":BAD},open("audit1.json","w"),indent=1,default=str)
print(f"\n--- pass {len(OK)}  FAIL {len(BAD)}")
