import json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':13,'font.family':'DejaVu Sans','axes.grid':True,
                     'grid.alpha':.3,'axes.spines.top':False,'axes.spines.right':False})
OR,TE,GY='#C1562B','#128A7D','#2B2B2B'
f=pd.read_csv('cr/forest_final.csv')
pooled=f[f.cls=='All OECD-2021 PFAS'].iloc[0]
cl=f[f.cls!='All OECD-2021 PFAS'].sort_values('rho').reset_index(drop=True)
fig,ax=plt.subplots(figsize=(12.4,6.4))
ys=list(range(len(cl)))[::-1]
for y,(_,r) in zip(ys,cl.iterrows()):
    ax.plot([r.lo,r.hi],[y,y],color=OR,lw=3.2,solid_capstyle='round')
    ax.plot([r.rho],[y],'o',color=OR,ms=10)
yp=len(cl)+0.8
ax.plot([pooled.lo,pooled.hi],[yp,yp],color=TE,lw=3.6,solid_capstyle='round')
ax.plot([pooled.rho],[yp],'o',color=TE,ms=11)
ax.axhline(len(cl)+0.1,color='0.75',ls='--',lw=1.2)
ax.axvline(0,color='k',lw=1.4)
# ys[i] is the row that cl.iloc[i] was drawn on: label in the SAME order, not reversed
labels=[f'{r.cls}  ($n$={int(r.n)})' for _,r in cl.iterrows()]
ax.set_yticks(ys+[yp])
ax.set_yticklabels(labels+[f'All OECD-2021 PFAS  ($n$={int(pooled.n)})'])
ax.set_ylim(-0.8,yp+0.8)
ax.set_xlim(-1.08,1.08)
ax.set_xlabel(r'external Spearman $\rho$ on the acute median   (95% bootstrap CI)')
ax.set_title('Transfer to PFAS, resolved by OECD (2021) structure class',
             fontweight='bold',loc='left',pad=14)
ax.text(0.5,-0.155,'No single class resolves: every class interval covers zero, and at these '
        'sample sizes\n80% power detects only '
        f'$|\\rho|\\geq${cl.floor.min():.2f}–{cl.floor.max():.2f}.',
        transform=ax.transAxes,ha='center',va='top',fontsize=11.5,color=GY)
plt.tight_layout()
fig.canvas.draw()
boxes=[]
items=[t for t in [ax.title,ax.xaxis.label,ax.yaxis.label]+list(ax.texts) if t.get_text()]
for axis,lim in ((ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())):
    for loc,lab in zip(axis.get_ticklocs(),axis.get_ticklabels()):
        if min(lim)-1e-9<=loc<=max(lim)+1e-9 and lab.get_text(): items.append(lab)
for t in items: boxes.append((t.get_text()[:30],t.get_window_extent()))
bad=[(p[0],q[0]) for i,p in enumerate(boxes) for q in boxes[i+1:] if p[1].overlaps(q[1])]
W,H=fig.canvas.get_width_height()
off=[b[0] for b in boxes if b[1].x0<-1 or b[1].y0<-1 or b[1].x1>W+1 or b[1].y1>H+1]
print('checked',len(boxes),'items'); print('overlaps:',bad if bad else 'none')
print('off-canvas:',off if off else 'none')
assert not bad and not off,'LAYOUT FAILURE'
# verify each tick label sits on the row of the class it names
tickmap={round(t,3):l.get_text() for t,l in zip(ax.get_yticks(),ax.get_yticklabels())}
for y,(_,r) in zip(ys,cl.iterrows()):
    assert r.cls in tickmap[round(float(y),3)], (
        f'label/value mismatch: row {y} drawn for {r.cls} but labelled {tickmap[round(float(y),3)]}')
assert 'All OECD' in tickmap[round(float(yp),3)]
print('label-to-value mapping verified for all', len(cl)+1, 'rows')
fig.savefig('fig_03_subclass.png',dpi=200,bbox_inches='tight')
print('wrote fig_03_subclass.png')
