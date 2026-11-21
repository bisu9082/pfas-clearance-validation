import pandas as pd, numpy as np, json
pd.set_option("display.width",250)

af  = pd.read_csv("/home/claude/wq/af_proper.csv")          # 58 parents w/ acute
cm  = pd.read_csv("char_map.csv")
conc= pd.read_pickle("conc.pkl")
DROP=["N-Methylperfluorooctanesulfonamide (N-MeFOSA) (D3, 98%)",
      "N-Ethylperfluorooctanesulfonamidoacetic acid (N-EtFOSAA) (N-ethyl-D5, 98%)"]
conc=conc[~conc.CharacteristicName.isin(DROP)].copy()
conc=conc.merge(cm[["char","cas","parent"]],left_on="CharacteristicName",right_on="char",how="left")
print("measured rows:",len(conc),"  with parent SMILES:",conc.parent.notna().sum())

m = conc.merge(af,on="parent",how="inner")
print("measured rows joined to a PNEC-bearing parent:",len(m))
print("substances (parents) matched:",m.parent.nunique(),"of",af.parent.nunique())
print("sites:",m.MonitoringLocationIdentifier.nunique())

# PNEC (ug/L): 10^logTox (mg/L) * 1000 / AF
def pnec(lt,AF): return (10.0**lt)*1000.0/AF
m["ugL"]=m.ngL/1000.0

sch = {
 "A_min_fixed"   : ("logTox",           1000.0),   # published convention, fixed AF
 "B_spmed_fixed" : ("logTox_sp_median", 1000.0),   # species-resolved median, fixed AF
 "A_min_tiered"  : ("logTox",           "AF_strict"),
 "B_spmed_tiered": ("logTox_sp_median", "AF_strict"),
}
out={}
for k,(col,AF) in sch.items():
    A = m[AF] if isinstance(AF,str) else AF
    m["PNEC_"+k]=pnec(m[col],A)
    m["RQ_"+k]=m.ugL/m["PNEC_"+k]
    ex=(m["RQ_"+k]>1)
    sites=m.groupby("MonitoringLocationIdentifier")["RQ_"+k].max()
    out[k]=dict(n_samples=int(len(m)),n_exceed=int(ex.sum()),pct_samples=100*ex.mean(),
                n_sites=int(sites.size),n_sites_exceed=int((sites>1).sum()),
                pct_sites=100*(sites>1).mean(),
                median_PNEC_ugL=float(m["PNEC_"+k].median()),
                median_RQ=float(m["RQ_"+k].median()))
r=pd.DataFrame(out).T
print("\n=== exceedance (RQ>1) by endpoint convention x AF scheme ===")
print(r.to_string())
m.to_pickle("rq.pkl"); json.dump(out,open("rq.json","w"),indent=1)
