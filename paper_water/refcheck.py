import json,subprocess
CA="/root/.ccr/ca-bundle.crt"
def cr(doi):
    o=subprocess.run(["curl","-sS","--cacert",CA,f"https://api.crossref.org/works/{doi}"],
                     capture_output=True,text=True,timeout=90).stdout
    try: return json.loads(o)["message"]
    except Exception as e: return {"_err":str(e),"_raw":o[:200]}
DOIS={"ruffle2024criteria":"10.1002/ieam.4776",
      "ding2026pfaswater":"10.1016/j.watres.2025.125054"}
for k,d in DOIS.items():
    m=cr(d)
    if "_err" in m: print(k,"FAIL",m); continue
    au=m.get("author",[])
    print(f"\n--- {k}  {d}")
    print("  title  :",(m.get('title') or [''])[0])
    print("  authors:",len(au),"| surnames:",", ".join(a.get('family','?') for a in au[:10]))
    print("  journal:",(m.get('container-title') or [''])[0])
    print("  year   :",m.get("issued",{}).get("date-parts",[[None]])[0][0],
          "vol",m.get("volume"),"pages",m.get("page"),"artno",m.get("article-number"))
