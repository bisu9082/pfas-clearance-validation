import json,subprocess,re,sys
CA="/root/.ccr/ca-bundle.crt"
def get(doi):
    o=subprocess.run(["curl","-sS","--cacert",CA,f"https://api.crossref.org/works/{doi}"],
                     capture_output=True,text=True,timeout=120).stdout
    try: return json.loads(o)["message"]
    except Exception: return None
def esc(t):
    t=re.sub(r"<[^>]+>","",t or "")
    t=t.replace("&amp;","&")
    for a,b in [("&","\\&"),("%","\\%"),("_","\\_"),("#","\\#"),
                ("–","--"),("—","---"),("’","'"),("“","``"),("”","''")]:
        t=t.replace(a,b)
    return re.sub(r"\s+"," ",t).strip()
def prot(t):
    # brace acronyms so BibTeX keeps their case
    return re.sub(r"\b([A-Z]{2,}(?:-[A-Z0-9]+)?|PFOA|PFOS|PFAS|HFPO-DA|NOEC|EC10|HC5|SSD|QSAR|EU|US)\b",
                  lambda m:"{"+m.group(0)+"}", t)
KEYS={
# --- SSD / criteria derivation
"belanger2019ssdsize":"10.1002/etc.4445",
"iwasaki2022afssd":"10.1016/j.ecoenv.2022.114170",
"iwasaki2015noecec10":"10.1021/acs.est.5b02069",
"hiki2018mindata":"10.1016/j.envpol.2018.03.003",
"hoondert2014pnecvar":"10.1002/ieam.1473",
"giesy2011pfospnec":"10.1002/etc.460",
"qi2014pfoscriteria":"10.1016/j.scitotenv.2013.09.094",
"zhao2016ssdmodel":"10.1016/j.ecoenv.2015.11.039",
"maltby2007eqswfd":"10.1897/ieam_2006-045.1",
"brack2023wfdcomment":"10.1186/s12302-023-00726-3",
# --- PFAS occurrence and alternatives
"salvatore2023pfasoccur":"10.1016/j.watres.2023.120528",
"zhao2013shortchain":"10.1021/es402120y",
"gaballah2020genxtox":"10.1016/j.aquatox.2023.106655",
"brandsma2022genxsub":"10.1016/j.envres.2021.111995",
"hagenaars2017pfoalow":"10.1016/j.reprotox.2017.01.009",
# --- data and databases
"laan2024ecotoxr":"10.1016/j.chemosphere.2024.143078",
"hull2000dataquality":"10.1080/10807030091124176",
"yang2015samplingbias":"10.3897/bdj.3.e5361",
# --- author's own work (new)
"lee2027uavforensics":"10.1016/j.trac.2026.119163",
"bae2025depdrivers":"10.1016/j.jenvrad.2025.107725",
"kim2026multitaskecotox":"10.1021/acs.est.6c12447",
}
out=[]
print(f"{'key':24s} {'status':9s} title / journal / year / vol / pages")
for k,d in KEYS.items():
    m=get(d)
    if m is None: print(f"{k:24s} FAIL      {d}"); continue
    au=m.get("author",[])
    names=" and ".join(
        (f"{a.get('family','')}, {a.get('given','')}".strip(", ") if a.get("family") else a.get("name",""))
        for a in au)
    ti=esc((m.get("title") or [""])[0]); jr=esc((m.get("container-title") or [""])[0])
    yr=m.get("issued",{}).get("date-parts",[[None]])[0][0]
    vol=m.get("volume"); pg=m.get("page") or m.get("article-number")
    e=f"@article{{{k},\n  author  = {{{esc(names)}}},\n  title   = {{{prot(ti)}}},\n  journal = {{{jr}}},\n"
    if vol: e+=f"  volume  = {{{vol}}},\n"
    if pg:  e+=f"  pages   = {{{pg}}},\n"
    e+=f"  year    = {{{yr}}},\n  doi     = {{{d}}}\n}}\n"
    out.append(e)
    print(f"{k:24s} VERIFIED  {ti[:52]} | {jr[:26]} | {yr} | {vol} | {pg} | {len(au)} au")
open("newrefs.bib","w").write("\n".join(out))
print("\nwrote newrefs.bib with %d entries"%len(out))
