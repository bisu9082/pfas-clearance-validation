import pandas as pd, numpy as np
C=pd.read_csv("/home/claude/wq/c7_ssd.csv")
m=pd.read_pickle("m_ladder.pkl").merge(C,on="parent",how="left")
m["ugL"]=m.ngL/1000.0
CONV=[("C1_min_all","lowest, all endpoints (as published)",1000.0),
      ("C3_min_acute","lowest acute",1000.0),
      ("C5_spmed_acute","median acute, species-resolved",1000.0),
      ("C6_min_noec","lowest chronic NOEC, EU-TGD AF","AF_C6"),
      ("C7_min","lowest apical chronic SMCV, AF 10",10.0),
      ("C7_HC5","SSD HC5 of apical chronic SMCVs (regulatory style)",1.0)]
rows=[]
for col,lbl,AF in CONV:
    d=m[m[col].notna()].copy()
    a=d[AF] if isinstance(AF,str) else AF
    if isinstance(AF,str): d=d[d[AF].notna()]; a=d[AF]
    d["PNEC"]=(10.0**d[col])*1000.0/a
    d["RQ"]=d.ugL/d.PNEC
    st=d.groupby("MonitoringLocationIdentifier").RQ.max()
    rows.append(dict(convention=lbl,n_samples=len(d),n_sub=d.parent.nunique(),
        med_PNEC_ngL=1000*d.PNEC.median(),
        pct_samples=100*(d.RQ>1).mean(),pct_sites=100*(st>1).mean()))
r=pd.DataFrame(rows)
pd.set_option("display.width",200); pd.set_option("display.max_colwidth",52)
print(r.to_string(index=False,float_format=lambda x:f"{x:,.4g}"))
r.to_csv("exceedance_c7.csv",index=False)
