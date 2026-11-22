import pandas as pd, numpy as np
m=pd.read_pickle("m_ladder.pkl"); C=pd.read_csv("/home/claude/wq/c7_ssd.csv")
m=m.merge(C,on="parent",how="left")
raw=pd.read_pickle("all_raw.pkl"); cm=pd.read_csv("char_map.csv")
raw=raw.merge(cm[["char","cas","parent"]],left_on="CharacteristicName",right_on="char",how="left")
raw["unit"]=raw["ResultMeasure/MeasureUnitCode"].astype(str).str.strip().str.lower()
raw["val"]=pd.to_numeric(raw.ResultMeasureValue,errors="coerce")
raw["quant"]=raw.unit.isin(["ng/l","ug/l"])&raw.val.notna()
NAME={"375-22-4":"PFBA","335-67-1":"PFOA","1763-23-1":"PFOS","307-24-4":"PFHxA",
"375-73-5":"PFBS","2706-90-3":"PFPeA","375-85-9":"PFHpA","355-46-4":"PFHxS",
"375-95-1":"PFNA","754-91-6":"FOSA","335-76-2":"PFDA","27619-97-2":"6:2 FTS",
"307-55-1":"PFDoA","2058-94-8":"PFUnA","72629-94-8":"PFTrDA","30334-69-1":"FBSA",
"812-70-4":"7:3 FTCA","914637-49-3":"5:3 FTCA","41997-13-1":"FHxSA",
"24448-09-7":"N-MeFOSE","13252-13-6":"HFPO-DA","4151-50-2":"N-EtFOSA"}
m["s"]=m.cas.map(NAME).fillna(m.cas)
def pn(c,af): return (10.0**c)*1e6/af   # ng/L
print("="*100)
print("S3  CRITERIA VERSUS THE ANALYTICAL FLOOR, AND THE EFFECT OF NON-DETECTS")
print("="*100)
rows=[]
for s,g in m.groupby("s"):
    d=g.iloc[0]; par=d.parent
    rr=raw[raw.parent==par]
    det=rr.quant.mean() if len(rr) else np.nan
    pos=g.ngL[g.ngL>0]
    rl=pos.min() if len(pos) else np.nan   # smallest quantified value; a weak proxy for the reporting limit
    c1=pn(d.C1_min_all,1000.0)
    c6=pn(d.C6_min_noec,d.AF_C6) if pd.notna(d.C6_min_noec) and pd.notna(d.AF_C6) else np.nan
    c7=pn(d.C7_HC5,1.0) if pd.notna(d.C7_HC5) else np.nan
    rows.append(dict(s=s,n_all=len(rr),det_pct=100*det,rl_proxy=rl,
        C1=c1,C6=c6,C7=c7,C1_below_rl=c1<rl,
        C6_below_rl=(c6<rl) if pd.notna(c6) else None,
        exc1_quant=100*(g.ngL>c1).mean(),
        exc1_all=100*(g.ngL>c1).sum()/len(rr) if len(rr) else np.nan))
t=pd.DataFrame(rows).sort_values("n_all",ascending=False)
pd.set_option("display.width",230)
print(t.to_string(index=False,float_format=lambda x:f"{x:,.4g}"))
nb1=t.C1_below_rl.sum(); nb6=t.C6_below_rl.dropna().sum()
print(f"\n  substances whose C1 criterion lies BELOW the smallest value ever quantified: {nb1} / {len(t)}")
print(f"  same for C6: {int(nb6)} / {t.C6_below_rl.notna().sum()}")
print(f"\n  overall detection frequency across all returned results: {100*raw.quant.mean():.1f}%")
tot=len(raw[raw.parent.notna()])
print(f"  results for PNEC-bearing substances: {tot:,}  of which quantified {int(raw[raw.parent.notna()].quant.sum()):,}")
print("\n  exceedance of C1 as a fraction of ALL results (non-detects counted as not exceeding):")
q=m.copy(); q["c1"]=pn(q.C1_min_all,1000.0)
ex=(q.ngL>q.c1).sum()
print(f"    {ex:,} / {tot:,} = {100*ex/tot:.1f}%   (vs {100*(q.ngL>q.c1).mean():.1f}% of quantified results)")
t.to_csv("si_nd.csv",index=False)
