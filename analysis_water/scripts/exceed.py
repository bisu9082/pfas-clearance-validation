import pandas as pd, numpy as np, json
L=pd.read_csv("ladder.csv"); cm=pd.read_csv("char_map.csv")
conc=pd.read_pickle("conc.pkl")
DROP=["N-Methylperfluorooctanesulfonamide (N-MeFOSA) (D3, 98%)",
      "N-Ethylperfluorooctanesulfonamidoacetic acid (N-EtFOSAA) (N-ethyl-D5, 98%)"]
conc=conc[~conc.CharacteristicName.isin(DROP)]
conc=conc.merge(cm[["char","cas","parent"]],left_on="CharacteristicName",right_on="char",how="left")
m=conc.merge(L,on="parent",how="inner"); m["ugL"]=m.ngL/1000.0
print("samples %d   sites %d   substances %d   %s..%s"%(
 len(m),m.MonitoringLocationIdentifier.nunique(),m.parent.nunique(),
 m.ActivityStartDate.min(),m.ActivityStartDate.max()))

CONV=[("C1_min_all"   ,"lowest value, all endpoints (as published)",1000.0),
      ("C2_min_noLOEC","lowest value, LOECs removed"              ,1000.0),
      ("C3_min_acute" ,"lowest acute value"                       ,1000.0),
      ("C4_med_acute" ,"median acute, record-resolved"            ,1000.0),
      ("C5_spmed_acute","median acute, species-resolved"          ,1000.0),
      ("C6_min_noec"  ,"lowest chronic NOEC, EU-TGD AF"           ,"AF_C6")]
rows=[]
for col,lbl,AF in CONV:
    d=m[m[col].notna()].copy()
    if isinstance(AF,str): d=d[d[AF].notna()]; a=d[AF]
    else: a=AF
    d["PNEC"]=(10.0**d[col])*1000.0/a
    d["RQ"]=d.ugL/d.PNEC
    st=d.groupby("MonitoringLocationIdentifier").RQ.max()
    rows.append(dict(convention=lbl,n_samples=len(d),n_sub=d.parent.nunique(),
      med_PNEC_ngL=1000*d.PNEC.median(),
      pct_samples_RQgt1=100*(d.RQ>1).mean(),
      pct_sites_RQgt1=100*(st>1).mean(),
      median_RQ=d.RQ.median()))
r=pd.DataFrame(rows)
pd.set_option("display.width",230); pd.set_option("display.max_colwidth",46)
print("\n"+"="*118)
print("EXCEEDANCE OF THE DERIVED CRITERION IN US SURFACE WATER, BY CURATION CONVENTION")
print("="*118)
print(r.to_string(index=False,float_format=lambda x:f"{x:,.4g}"))
print("\nspread, published vs species-resolved median (same AF):")
a=r.iloc[0]; b=r.iloc[4]
print("  samples exceeding: %.2f%% -> %.2f%%  (%.0fx)"%(a.pct_samples_RQgt1,b.pct_samples_RQgt1,
      a.pct_samples_RQgt1/max(b.pct_samples_RQgt1,1e-9)))
print("  sites exceeding:   %.2f%% -> %.2f%%"%(a.pct_sites_RQgt1,b.pct_sites_RQgt1))
print("\nbenchmarks for orientation: EU EQS PFOS inland SW 0.65 ng/L; US EPA MCL PFOA/PFOS 4 ng/L")
r.to_csv("exceedance.csv",index=False); m.to_pickle("m_ladder.pkl")
