import pandas as pd, numpy as np, json
from scipy import stats
OK=[];BAD=[]
def chk(n,c,a,tol=0.0):
    good=abs(float(c)-float(a))<=tol if tol else str(c)==str(a)
    (OK if good else BAD).append((n,c,a)); print(("  OK   " if good else "  **BAD** ")+f"{n:44s} text={c} data={a}")
d=pd.read_csv("/home/claude/wq/broad/broad_summary.csv")
rows=json.load(open("/home/claude/wq/broad/broad_subsample.json"))
chk("broad chemicals",4757,len(d))
x=np.log10(d.n_rec)
chk("broad min rho",-0.462,round(stats.spearmanr(x,d.acute_min)[0],3),0.0006)
chk("broad median rho",-0.160,round(stats.spearmanr(x,d.acute_median)[0],3),0.0006)
r5=[r for r in rows if r["k"]==5][0]
chk("k=5 eligible",1977,r5["n_eligible"])
chk("k=5 min baseline",-0.409,round(r5["min_baseline"],3),0.0006)
chk("k=5 min equal depth",-0.192,round(r5["min_depth_equal"],3),0.004)
chk("k=5 depth effect",-0.217,round(r5["depth_effect"],3),0.004)
chk("k=5 median baseline",-0.140,round(r5["med_baseline"],3),0.0006)
chk("k=5 median equal depth",-0.134,round(r5["med_depth_equal"],3),0.004)
de=[abs(r["depth_effect"]) for r in rows]
chk("depth effect on min, min across k",0.16,round(min(de),2),0.006)
chk("depth effect on min, max across k",0.24,round(max(de),2),0.006)
dm=[abs(r["med_depth_equal"]-r["med_baseline"]) for r in rows]
chk("depth effect on median, max across k",0.01,round(max(dm),2),0.006)
import subprocess
nrec=subprocess.run(["grep","-c","acute freshwater records","/home/claude/wq/broad/broad.log"],capture_output=True,text=True)
log=open("/home/claude/wq/broad/broad.log").read()
chk("broad acute records",106095,int([l for l in log.split("\n") if "acute freshwater records" in l][0].split(":")[1]))
print(f"\n--- pass {len(OK)}  FAIL {len(BAD)}")
