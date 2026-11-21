import json,subprocess,urllib.parse,sys
CA="/root/.ccr/ca-bundle.crt"
def q(query,n=6,filt="type:journal-article"):
    u=("https://api.crossref.org/works?query.bibliographic="+urllib.parse.quote(query)
       +f"&rows={n}&filter={filt}&select=DOI,title,container-title,issued,volume,page,author,is-referenced-by-count")
    o=subprocess.run(["curl","-sS","--cacert",CA,u],capture_output=True,text=True,timeout=120).stdout
    try: return json.loads(o)["message"]["items"]
    except Exception as e: return []
TOPICS=[
 "per- and polyfluoroalkyl substances occurrence United States surface water rivers national survey",
 "PFAS aquatic toxicity review ecological risk assessment freshwater organisms",
 "species sensitivity distribution sample size uncertainty HC5 estimation number of species",
 "no observed effect concentration NOEC criticism EC10 regression ecotoxicology statistics",
 "publication bias data availability ecotoxicological databases chemical risk assessment",
 "short chain PFAS replacement alternatives environmental occurrence toxicity",
 "risk quotient screening level ecological risk assessment surface water monitoring",
 "PFAS grouping read-across regulatory assessment thousands of substances",
 "extreme value order statistics environmental quality standard derivation",
 "machine learning prediction aquatic toxicity data bias training set composition",
]
out={}
for t in TOPICS:
    items=q(t)
    print("\n### "+t[:72])
    for it in items:
        au=it.get("author",[])
        print(f"  {it['DOI']}  {(it.get('title') or [''])[0][:84]}")
        print(f"      {(it.get('container-title') or [''])[0][:56]} | {it.get('issued',{}).get('date-parts',[[None]])[0][0]} | cited {it.get('is-referenced-by-count')}")
    out[t]=items
json.dump(out,open("cand.json","w"))
