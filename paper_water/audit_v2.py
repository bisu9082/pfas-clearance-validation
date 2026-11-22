import re,json,sys
import pandas as pd, numpy as np
from scipy import stats
s=open("main.tex").read()+open("si.tex").read()
out=pd.read_csv("/home/claude/wq/v2/ladder_v2.csv")
canon=json.load(open("/home/claude/wq/v2/canonical_v2.json"))
ci=pd.read_csv("/home/claude/wq/v2/hc5_ci.csv")
fails=[];checks=0
def chk(name,claim_in_text,computed,tol=0.02):
    global checks; checks+=1
    ok = abs(claim_in_text-computed)<=tol*max(1e-9,abs(computed))
    if not ok: fails.append(f"{name}: text {claim_in_text} vs computed {computed}")
def present(name,pat):
    global checks; checks+=1
    if not re.search(pat,s): fails.append(f"{name}: pattern not found -> {pat}")

ln=np.log10(out.n_records)
r_c6p=stats.spearmanr(ln[out.C6_prime.notna()],out.C6_prime.dropna())[0]
r_c7b=stats.spearmanr(np.log10(out[out.HC5_d.notna()].n_records),out.HC5_d.dropna())[0]
r_c7a=stats.spearmanr(np.log10(out[out.minSMCV_d.notna()].n_records),out.minSMCV_d.dropna())[0]
chk("C6' rho",-0.590,r_c6p); chk("C7b rho",-0.681,r_c7b); chk("C7a rho",-0.609,r_c7a)
chk("sigma share d",63,canon["hc5_decomp_d"]["sigma_share_pct"],0.01)
chk("sigma share nd",69,canon["hc5_decomp_nd"]["sigma_share_pct"],0.01)
chk("partial|SD",-0.004,canon["sd_control"]["partial_given_sd"],0.5)
chk("rho logn sd",0.783,canon["sd_control"]["rho_logn_sd"])
chk("n SSD-able",13,canon["n_ssd_able_duration_filtered"],0)
chk("HC5 CI max width",7.21,ci.width.max(),0.01)
chk("HC5 CI min width",0.78,ci.width.min(),0.02)
# species-equalised range quoted as 50-61%
lo=min(canon[f"species_equalised_s{k}"]["pct_removed"] for k in (3,5,10))
hi=max(canon[f"species_equalised_s{k}"]["pct_removed"] for k in (3,5,10))
chk("species-eq low",50,lo,0.02); chk("species-eq high",61,hi,0.02)
# text-presence checks for the corrected claims
present("title says support-dependent",r"support-dependent artefact that regulatory screening narrows")
present("no 'same data and the same water'",r"^(?!.*same data and the same water)")
for bad,label in [(r"extreme-order-statistic artefact that regulatory screening does not remove","old title"),
                  (r"applied to the same data and the same water","false sameness claim"),
                  (r"slightly more than half","withdrawn half-and-half claim"),
                  (r"a real association between hazard and attention, which a central","withdrawn hazard claim"),
                  (r"PNEC is the median across substances\.\}","old Table 1 caption"),
                  (r"do report developmental responses","wrong endpoint description")]:
    checks+=1
    if re.search(bad,s): fails.append(f"STALE CLAIM present: {label}")
# cross-document consistency
m=open("main.tex").read(); si=open("si.tex").read()
checks+=1
if re.search(r'ead\{([^}]*)\}',m).group(1)!=re.search(r'ead\{([^}]*)\}',si).group(1):
    fails.append("corresponding e-mail differs between main and SI")
print(f"checks run: {checks}   failures: {len(fails)}")
for f in fails: print("  FAIL:",f)
sys.exit(1 if fails else 0)
