import json,subprocess,sys
pats=["perfluoro","fluorotelomer","PFOS","PFOA","sulfonamido","HFPO","polyfluoro","fluorobutane","fluorohexane","fluorooctane","fluoroheptane","fluorodecane","fluorononane"]
names=set()
for p in pats:
    u=f"https://www.waterqualitydata.us/Codes/characteristicname?text={p}&mimeType=json"
    o=subprocess.run(["curl","-sS",u],capture_output=True,text=True,timeout=90).stdout
    try: d=json.loads(o)["codes"]
    except Exception as e: print("fail",p,e); continue
    for x in d: names.add(x["value"])
    print(f"{p:18s} {len(d):4d}  cum={len(names)}")
names=sorted(names)
open("charnames.txt","w").write("\n".join(names))
print("TOTAL unique characteristic names:",len(names))
