import json,subprocess,urllib.parse
CA="/root/.ccr/ca-bundle.crt"
def q(query,n=6):
    u=("https://api.crossref.org/works?query.bibliographic="+urllib.parse.quote(query)
       +f"&rows={n}&filter=type:journal-article&select=DOI,title,container-title,issued,is-referenced-by-count")
    o=subprocess.run(["curl","-sS","--cacert",CA,u],capture_output=True,text=True,timeout=120).stdout
    try: return json.loads(o)["message"]["items"]
    except Exception: return []
T=["short-chain perfluoroalkyl acids environmental occurrence persistence surface water",
   "fluorinated alternatives GenX HFPO-DA aquatic toxicity",
   "predicted no effect concentration derivation assessment factor uncertainty chemicals",
   "PFAS surface water monitoring national rivers streams United States Geological Survey",
   "ecotoxicity database ECOTOX data curation quality aquatic",
   "6:2 fluorotelomer sulfonate toxicity aquatic organisms",
   "water framework directive environmental quality standards derivation chemicals",
   "PFOS species sensitivity distribution aquatic life criteria derivation",
   "testing effort bias chemical hazard data richness regulatory attention",
   "QSAR applicability domain aquatic toxicity PFAS prediction",
   "zebrafish Danio rerio PFAS sublethal gene expression low concentration",
   "chemical mixture PFAS cumulative risk surface water exceedance"]
for t in T:
    print("\n### "+t[:70])
    for it in q(t):
        print(f"  {it['DOI']}  {(it.get('title') or [''])[0][:80]}")
        print(f"      {(it.get('container-title') or [''])[0][:50]} | {it.get('issued',{}).get('date-parts',[[None]])[0][0]} | cit {it.get('is-referenced-by-count')}")
