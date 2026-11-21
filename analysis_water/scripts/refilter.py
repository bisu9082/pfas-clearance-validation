import pandas as pd, numpy as np
df=pd.read_pickle("all_raw.pkl")
df["unit"]=df["ResultMeasure/MeasureUnitCode"].astype(str).str.strip().str.lower()
df["val"]=pd.to_numeric(df.ResultMeasureValue,errors="coerce")
dc=df.ResultDetectionConditionText.fillna("")
NONDET=dc.str.contains("Not Detected|Below|Not Reported|Systematic Contamination",case=False,regex=True)
conc=df[df.unit.isin(["ng/l","ug/l"]) & df.val.notna() & (df.val>0) & (~NONDET)].copy()
conc["ngL"]=np.where(conc.unit=="ug/l",conc.val*1000.0,conc.val)
iso=conc.CharacteristicName.str.contains(r"13C|\(D3|\(N-ethyl-D5|\[",regex=True,na=False)
conc=conc[~iso]
DROP=["N-Methylperfluorooctanesulfonamide (N-MeFOSA) (D3, 98%)",
      "N-Ethylperfluorooctanesulfonamidoacetic acid (N-EtFOSAA) (N-ethyl-D5, 98%)"]
conc=conc[~conc.CharacteristicName.isin(DROP)]
print("detections after the corrected filter: %d  (was 23,237 incl. %d non-detect zeros)"%(len(conc),615))
print("sites:",conc.MonitoringLocationIdentifier.nunique())
conc.to_pickle("conc.pkl")

cm=pd.read_csv("char_map.csv"); L=pd.read_csv("ladder.csv")
c=conc.merge(cm[["char","cas","parent"]],left_on="CharacteristicName",right_on="char",how="left")
m=c.merge(L,on="parent",how="inner"); m["ugL"]=m.ngL/1000.0
print("joined to a criterion-bearing parent: %d measurements, %d stations, %d substances"%(
 len(m),m.MonitoringLocationIdentifier.nunique(),m.parent.nunique()))
m.to_pickle("m_ladder.pkl")

C=pd.read_csv("/home/claude/wq/c7_ssd.csv"); mm=m.merge(C,on="parent",how="left")
CONV=[("C1_min_all","lowest, all endpoints (as published)",1000.0),
      ("C2_min_noLOEC","lowest, LOECs removed",1000.0),
      ("C3_min_acute","lowest acute",1000.0),
      ("C4_med_acute","median acute, record-resolved",1000.0),
      ("C5_spmed_acute","median acute, species-resolved",1000.0),
      ("C6_min_noec","lowest chronic NOEC, EU-TGD AF","AF_C6"),
      ("C7_min","lowest apical chronic SMCV, AF 10",10.0),
      ("C7_HC5","SSD HC5 of apical chronic SMCVs",1.0)]
# denominator for the all-results view
raw=df.merge(cm[["char","parent"]],left_on="CharacteristicName",right_on="char",how="left")
rows=[]
for col,lbl,AF in CONV:
    d=mm[mm[col].notna()].copy()
    if isinstance(AF,str): d=d[d[AF].notna()]; a=d[AF]
    else: a=AF
    d["PNEC"]=(10.0**d[col])*1000.0/a; d["RQ"]=d.ugL/d.PNEC
    st=d.groupby("MonitoringLocationIdentifier").RQ.max()
    nall=len(raw[raw.parent.isin(d.parent.unique())])
    rows.append(dict(convention=lbl,n_det=len(d),n_sub=d.parent.nunique(),
      med_PNEC_ngL=1000*d.PNEC.median(),
      pct_samples=100*(d.RQ>1).mean(), pct_sites=100*(st>1).mean(),
      n_all=nall, pct_all=100*(d.RQ>1).sum()/nall))
r=pd.DataFrame(rows)
pd.set_option("display.width",230); pd.set_option("display.max_colwidth",44)
print("\n"+r.to_string(index=False,float_format=lambda x:f"{x:,.4g}"))
r.to_csv("exceedance_c7.csv",index=False)
