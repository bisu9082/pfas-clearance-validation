import pandas as pd, numpy as np, re
m=pd.read_pickle("/home/claude/wq/wqp/m_ladder.pkl")
C=pd.read_csv("/home/claude/wq/c7_ssd.csv"); m=m.merge(C,on="parent",how="left")
NAME={"375-22-4":"PFBA","335-67-1":"PFOA","1763-23-1":"PFOS","307-24-4":"PFHxA",
"375-73-5":"PFBS","2706-90-3":"PFPeA","375-85-9":"PFHpA","355-46-4":"PFHxS",
"375-95-1":"PFNA","754-91-6":"FOSA","335-76-2":"PFDA","27619-97-2":"6:2 FTS",
"307-55-1":"PFDoA","2058-94-8":"PFUnA","72629-94-8":"PFTrDA","30334-69-1":"FBSA",
"812-70-4":"7:3 FTCA","914637-49-3":"5:3 FTCA","41997-13-1":"FHxSA",
"24448-09-7":"N-MeFOSE","13252-13-6":"HFPO-DA","4151-50-2":"N-EtFOSA"}
m["s"]=m.cas.map(NAME)
tex=open("tab_persub.tex").read()
rows=[l for l in tex.split("\n") if "&" in l and "\\\\" in l and not l.startswith("\\")]
bad=0;n=0
for l in rows:
    c=[x.strip() for x in l.replace("\\\\","").split("&")]
    if len(c)<11 or c[0] in ("Substance",""): continue
    s=c[0]; g=m[m.s==s]
    if not len(g): print("  **BAD** table row with no data:",s); bad+=1; continue
    n+=1
    def num(x): 
        x=x.replace(",","").replace("\\num{","").replace("}","")
        try: return float(x)
        except: return np.nan
    if int(num(c[1]))!=len(g): print(f"  **BAD** {s} n: table={c[1]} data={len(g)}"); bad+=1
    if int(num(c[2]))!=g.MonitoringLocationIdentifier.nunique():
        print(f"  **BAD** {s} sites: table={c[2]} data={g.MonitoringLocationIdentifier.nunique()}"); bad+=1
    if abs(num(c[3])-g.ngL.median())>0.006: print(f"  **BAD** {s} median: table={c[3]} data={g.ngL.median():.3f}"); bad+=1
    if int(num(c[4]))!=int(g.n_records.iloc[0]): print(f"  **BAD** {s} recs: table={c[4]} data={int(g.n_records.iloc[0])}"); bad+=1
    c1=(10**g.C1_min_all.iloc[0])*1e6/1000
    e1=100*(g.ngL>c1).mean()
    if abs(num(c[8])-e1)>0.06: print(f"  **BAD** {s} exc C1: table={c[8]} data={e1:.2f}"); bad+=1
print(f"\nTable 2: {n} substance rows checked, {bad} mismatches")

print("\n== ABSTRACT vs BODY ==")
t=open("main.tex").read()
ab=t[t.index("\\begin{abstract}"):t.index("\\end{abstract}")]
for pat,lbl in [(r"4744","curated records"),(r"87 PFAS","parents"),(r"20\\,682","measurements"),
                (r"1180","stations"),(r"76.5","C1 samples"),(r"0.11","C5 samples"),
                (r"4738","PFOS HC5"),(r"6800","Canada FEQG"),(r"-0.69","rho C7b"),
                (r"\+0.76","rho species"),(r"\+0.12","rho C5"),(r"810","shared stations"),
                (r"1232","PFOA records"),(r"11 toxicity records|1232 versus 11","PFHpA records"),
                (r"3.5","PFOA median"),(r"1.7","PFHpA median"),(r"98.2|100","pair exceed")]:
    print(("  OK  " if re.search(pat,ab) else "  --  ")+f"{lbl:22s} /{pat}/ in abstract")
