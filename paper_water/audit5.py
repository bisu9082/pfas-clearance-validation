import pandas as pd, numpy as np, json
from scipy import stats
OK=[];BAD=[]
def chk(n,c,a,tol=0.0):
    good=abs(float(c)-float(a))<=tol if tol else str(c)==str(a)
    (OK if good else BAD).append((n,c,a)); print(("  OK   " if good else "  **BAD** ")+f"{n:40s} text={c} data={a}")
# --- R2 assessment-factor decomposition
A=json.load(open("/home/claude/wq/wqp/af_compensation.json"))
chk("AF decomposition n",62,A["n"])
chk("slope NOEC",-2.123,round(A["slope_noec"],3),0.0006)
chk("slope AF",-0.302,round(A["slope_af"],3),0.0006)
chk("slope criterion",-1.821,round(A["slope_criterion"],3),0.0006)
chk("compensation delivered %",14,round(A["pct_compensation_delivered"]),0.5)
chk("passthrough %",86,round(A["pct_passthrough"]),0.5)
chk("AF range log",1.0,round(A["af_range_log"],1),0.05)
chk("NOEC range log",9.2,round(A["noec_range_log"],1),0.05)
# --- R5 model experiment
M=json.load(open("/home/claude/wq/broad/r5_model.json"))
S=pd.read_csv("/home/claude/wq/broad/r5_strata.csv")
chk("structure->effort R2",0.085,round(M["r2_effort"],3),0.0006)
chk("min target R2",0.371,round(M["r2_min"],3),0.0006)
chk("median target R2",0.422,round(M["r2_med"],3),0.0006)
chk("bias, 1 record",-0.66,round(S.min_bias.iloc[0],2),0.006)
chk("bias, 100+ records",2.33,round(S.min_bias.iloc[-1],2),0.006)
chk("bias swing, min",3.0,round(S.min_bias.max()-S.min_bias.min(),1),0.06)
chk("bias swing, median",0.47,round(S.med_bias.max()-S.med_bias.min(),2),0.006)
chk("min R2 in the deepest stratum",-0.951,round(S.min_r2.iloc[-1],3),0.0006)
chk("min rho range lo",0.57,round(S.min_rho.min(),2),0.006)
chk("min rho range hi",0.63,round(S.min_rho.max(),2),0.006)
chk("chemicals with a parsable structure",4311,len(pd.read_csv("/home/claude/wq/broad/r5_pred.csv")))
print(f"\n--- pass {len(OK)}  FAIL {len(BAD)}")
for b in BAD: print("   FAIL:",b)
