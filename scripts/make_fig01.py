import json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
plt.rcParams.update({'font.size':13,'font.family':'DejaVu Sans','axes.grid':True,
                     'grid.alpha':.3,'axes.spines.top':False,'axes.spines.right':False})
OR,TE,GY='#C1562B','#128A7D','#4A4A4A'
d=pd.read_csv('cr/ecotox_pfas_measured_recurated.csv')
ss=json.load(open('cr/subsample.json')); ep=json.load(open('cr/endpoint_fix_recurated.json'))
fig,ax=plt.subplots(1,3,figsize=(16.2,5.3))

# (a)
x=np.log10(d.n_records); y=d.logTox
ax[0].scatter(x,y,s=34,color=OR,alpha=.72,edgecolor='none')
b=np.polyfit(x,y,1); xs=np.linspace(x.min(),x.max(),50)
ax[0].plot(xs,np.polyval(b,xs),color='k',lw=2)
r=spearmanr(y,x).statistic
ax[0].text(.03,.05,f'$\\rho$ = {r:.2f}',transform=ax[0].transAxes,fontsize=15,fontweight='bold')
ax[0].set_xlabel('records per substance (log$_{10}$)')
ax[0].set_ylabel('protective minimum (log$_{10}$ mg/L)')
ax[0].set_title('(a)  The minimum tracks testing effort',fontweight='bold',loc='left')

# (b) two lines: eligible-subset baseline and depth-equalised
k=[r['k'] for r in ss['rows']]
base=[r['rho_baseline'] for r in ss['rows']]
sub=[r['rho_depth_equal'] for r in ss['rows']]
sd=[r['sd'] for r in ss['rows']]
pos=np.arange(len(k))
ax[1].plot(pos,base,'s--',color=GY,lw=2,ms=8,label='same substances, all their records')
ax[1].errorbar(pos,sub,yerr=sd,fmt='o-',color=TE,lw=2.4,ms=9,capsize=4,
               label='same substances, $k$ records each')
ax[1].axhline(ss['rho_all'],color=OR,ls=':',lw=2)
ax[1].annotate('all 87 substances',xy=(pos[0],ss['rho_all']),xytext=(pos[0]-0.08,ss['rho_all']+0.016),
               color=OR,fontsize=11,va='bottom',ha='left')
ax[1].set_xticks(pos); ax[1].set_xticklabels([str(v) for v in k])
ax[1].set_xlim(-0.45,len(k)-0.3)
ax[1].set_ylim(-0.82,-0.17)
ax[1].set_xlabel('records retained per substance ($k$)')
ax[1].set_ylabel(r'$\rho$ (minimum, log records)')
ax[1].set_title('(b)  Depth explains 40–48% at every $k$',fontweight='bold',loc='left')
ax[1].legend(loc='upper center',bbox_to_anchor=(0.5,1.0),fontsize=10.5,frameon=False)

# (c) rank, matched sample
lab=['minimum\n(endpoint swap)','acute median\n(endpoint swap)']
m=ep['matched']
vals=[m['minimum']['rho'],m['acute_median']['rho']]
los=[m['minimum']['ci'][0],m['acute_median']['ci'][0]]
his=[m['minimum']['ci'][1],m['acute_median']['ci'][1]]
cols=[OR,TE]
for i,(v,l,h,c) in enumerate(zip(vals,los,his,cols)):
    ax[2].errorbar([i],[v],yerr=[[v-l],[h-v]],fmt='o',color=c,ms=12,lw=2.6,capsize=6)
ax[2].axhline(0,color='k',lw=1.2)
ax[2].set_xticks(range(2)); ax[2].set_xticklabels(lab)
ax[2].set_xlim(-.55,1.55); ax[2].set_ylim(-.42,.82)
ax[2].set_ylabel(r'external Spearman $\rho$  (same 58 parents)')
ax[2].set_title('(c)  Rank transfer returns on the matched set',fontweight='bold',loc='left')
dd=m['diff']
ax[2].text(.5,.74,f"$\\Delta\\rho$ = {dd['mean']:+.2f}  [{dd['ci'][0]:+.2f}, {dd['ci'][1]:+.2f}]",
           ha='center',fontsize=12.5,fontweight='bold',color=GY)
plt.tight_layout()
fig.canvas.draw()
# ---- layout verification: measured extents, not guesses ----
boxes=[]
for n,a in enumerate(ax):
    items=[t for t in [a.title,a.xaxis.label,a.yaxis.label]+list(a.texts) if t.get_text()]
    lg=a.get_legend()
    if lg: items+=[t for t in lg.get_texts() if t.get_text()]
    # tick labels, restricted to ticks actually inside the current view
    for axis,lim in ((a.xaxis,a.get_xlim()),(a.yaxis,a.get_ylim())):
        locs=axis.get_ticklocs(); labs=axis.get_ticklabels()
        for loc,lab in zip(locs,labs):
            if min(lim)-1e-9 <= loc <= max(lim)+1e-9 and lab.get_text():
                items.append(lab)
    for t in items:
        boxes.append((f'ax{n}:{t.get_text()[:24]}', t.get_window_extent()))
bad=[(p[0],q[0]) for i_,p in enumerate(boxes) for q in boxes[i_+1:] if p[1].overlaps(q[1])]
W,H=fig.canvas.get_width_height()
off=[b[0] for b in boxes if b[1].x0<-1 or b[1].y0<-1 or b[1].x1>W+1 or b[1].y1>H+1]
print('text items checked   :',len(boxes))
print('overlapping text pairs:',bad if bad else 'none')
print('off-canvas text      :',off if off else 'none')
assert not bad and not off, 'LAYOUT FAILURE'
fig.savefig('fig_01_effort.png',dpi=200,bbox_inches='tight')
print('wrote fig_01_effort.png')
