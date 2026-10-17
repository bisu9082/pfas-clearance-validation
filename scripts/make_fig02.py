import json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':13,'font.family':'DejaVu Sans','axes.grid':True,
                     'grid.alpha':.3,'axes.spines.top':False,'axes.spines.right':False})
OR,TE,GY,BL='#C1562B','#128A7D','#4A4A4A','#2B5CA8'
inv=pd.read_csv('cr/final_inventory_predictions.csv'); d=inv[inv.in_domain]
rep=pd.read_csv('cr/rebuild_v2_conformal_repeats.csv')
fig,ax=plt.subplots(1,3,figsize=(16.6,5.4))

# (a) rebuilt external validation
ext=pd.read_csv('cr/external_pfas_predictions.csv') if False else None
vals=[('external $R^2$',0.366,None,None),(r'external $\rho$',0.460,0.23,0.65)]
for i,(lab,v,lo,hi) in enumerate(vals):
    if lo is None: ax[0].bar([i],[v],color=TE,width=.55)
    else: ax[0].errorbar([i],[v],yerr=[[v-lo],[hi-v]],fmt='o',color=TE,ms=12,lw=2.6,capsize=6)
    ax[0].text(i,v+.045,f'{v:+.2f}',ha='center',va='bottom',fontweight='bold',fontsize=13)
ax[0].axhline(0,color='k',lw=1.2)
ax[0].set_xticks([0,1]); ax[0].set_xticklabels(['$R^2$','Spearman $\\rho$\n(95% CI)'])
ax[0].set_xlim(-.6,1.6); ax[0].set_ylim(-.12,.82)
ax[0].set_ylabel('rebuilt model, 58 external PFAS')
ax[0].set_title('(a)  External validity of the rebuild',fontweight='bold',loc='left')

# (b) conditional width with its spread, against what the decision needs
w=rep.cond_width.values
ax[1].boxplot([w],positions=[0],widths=.45,patch_artist=True,
              boxprops=dict(facecolor=TE,alpha=.6,edgecolor=TE),
              medianprops=dict(color='k',lw=2),whiskerprops=dict(color=TE),
              capprops=dict(color=TE),showfliers=False)
# required half-width at a threshold T is (median in-domain prediction) - log10(T)
med=float(d.pred.median())
need1 = 2*(med-0.0)        # 1 mg/L  (GHS acute Category 1)
ax[1].axhline(need1,color=OR,ls='--',lw=2.2)
ax[1].text(.42,need1,'  width needed to clear\n  at 1 mg/L',color=OR,va='center',fontsize=11)
ax[1].set_xticks([0]); ax[1].set_xticklabels(['PFAS-conditional\n90% interval'])
ax[1].set_xlim(-.5,1.25); ax[1].set_ylim(1.9,6.4)
ax[1].set_ylabel('interval width (log$_{10}$ mg/L)')
ax[1].set_title('(b)  Achieved width is still too wide',fontweight='bold',loc='left')
ax[1].text(0,w.max()+0.22,f'{w.mean():.2f}',ha='center',va='bottom',fontweight='bold',fontsize=13)

# (c) clearance curve, with the band in which PFAS actually occur
ts=np.linspace(-6.4,3.2,400)
frac=[100*(d.lower>t).mean() for t in ts]
ax[2].axvspan(-6.4,-3,color='0.88',zorder=0)
ax[2].text(-4.7,41,'PFAS occurrence\nin surface water',ha='center',fontsize=10.5,color=GY)
ax[2].plot(ts,frac,color=TE,lw=2.6)
for t,lab,c,dy in [(-1,'0.1 mg/L',OR,0),(0,'1 mg/L (GHS acute 1)',GY,0),
                   (np.log10(46),'46 mg/L (matched)',BL,0)]:
    ax[2].axvline(t,color=c,ls='--',lw=1.8)
    f=100*(d.lower>t).mean()
    ax[2].plot([t],[f],'o',color=c,ms=9)
    # percentage to the RIGHT of its own line, at the marker; rotated name pinned to the top
    side = -1 if f > 20 else 1          # the 0.1 mg/L marker is high: label it on the left
    ax[2].text(t+0.18*side, f+1.6, f'{f:.1f}%', color=c, fontweight='bold', fontsize=12,
               ha='left' if side>0 else 'right', va='bottom')
    ax[2].text(t-0.16,47.0,lab,rotation=90,va='top',ha='right',fontsize=10,color=c)
ax[2].set_xlim(-6.4,3.2); ax[2].set_ylim(-2,48)
ax[2].set_xlabel('hazard threshold (log$_{10}$ mg/L)')
ax[2].set_ylabel('in-domain structures cleared (%)')
ax[2].set_title('(c)  Clearance fails when the threshold is matched',fontweight='bold',loc='left')
plt.tight_layout(); fig.canvas.draw()
boxes=[]
for n,a in enumerate(ax):
    items=[t for t in [a.title,a.xaxis.label,a.yaxis.label]+list(a.texts) if t.get_text()]
    for axis,lim in ((a.xaxis,a.get_xlim()),(a.yaxis,a.get_ylim())):
        for loc,lab in zip(axis.get_ticklocs(),axis.get_ticklabels()):
            if min(lim)-1e-9<=loc<=max(lim)+1e-9 and lab.get_text(): items.append(lab)
    for t in items: boxes.append((f'ax{n}:{t.get_text()[:22]}',t.get_window_extent()))
bad=[(p[0],q[0]) for i,p in enumerate(boxes) for q in boxes[i+1:] if p[1].overlaps(q[1])]
W,H=fig.canvas.get_width_height()
off=[b[0] for b in boxes if b[1].x0<-1 or b[1].y0<-1 or b[1].x1>W+1 or b[1].y1>H+1]
print('checked',len(boxes)); print('overlaps:',bad if bad else 'none'); print('off-canvas:',off if off else 'none')
assert not bad and not off,'LAYOUT FAILURE'
fig.savefig('fig_02_rebuild.png',dpi=200,bbox_inches='tight'); print('wrote fig_02_rebuild.png')
