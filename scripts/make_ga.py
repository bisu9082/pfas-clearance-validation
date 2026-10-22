"""JHM/Elsevier graphical abstract: >=1328 px wide, readable at 13 cm, longest side horizontal."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
plt.rcParams.update({'font.family':'DejaVu Sans'})
OR,TE,GY='#C1562B','#128A7D','#3A3A3A'

CM=1/2.54
W,H=13*CM,5.5*CM          # readable at 13 cm wide
DPI=320                    # 13 cm -> 1638 px, comfortably over the 1328 minimum
fig=plt.figure(figsize=(W,H),dpi=DPI)
ax=fig.add_axes([0.105,0.185,0.495,0.605])

d=pd.read_csv('/home/claude/fix/cr/ecotox_pfas_measured_recurated.csv')
x=np.log10(d.n_records); y=d.logTox
ax.scatter(x,y,s=7.5,color=OR,alpha=.62,edgecolor='none',zorder=3)
b=np.polyfit(x,y,1); xs=np.linspace(x.min(),x.max(),50)
ax.plot(xs,np.polyval(b,xs),color=GY,lw=1.9,zorder=4)

ax.set_xlim(-0.18,3.45); ax.set_ylim(-7.4,5.2)
ax.set_xticks([0,1,2,3]); ax.set_xticklabels(['1','10','100','1000'],fontsize=6.4)
ax.set_yticks([4,0,-4]); ax.set_yticklabels(['4','0','−4'],fontsize=6.4)
ax.tick_params(length=2.2,width=.6,pad=1.6)
for s in ('top','right'): ax.spines[s].set_visible(False)
for s in ('left','bottom'): ax.spines[s].set_linewidth(.7)
ax.set_xlabel('toxicity records per substance',fontsize=7.0,labelpad=2.0)
ax.set_ylabel('"most sensitive" value\n(log$_{10}$ mg/L)',fontsize=7.0,labelpad=2.0)
ax.grid(axis='y',alpha=.18,lw=.45)
ax.text(0.035,0.045,r'$\rho$ = −0.75',transform=ax.transAxes,fontsize=8.2,
        fontweight='bold',color=OR,ha='left',va='bottom')

fig.text(0.5,0.955,'A protective endpoint that measures testing history',
         ha='center',va='top',fontsize=9.0,fontweight='bold',color=GY)

# right-hand panel: short lines, then auto-shrink until every one fits the margin
bx, RIGHT = 0.655, 0.985
# three evenly spaced blocks of two lines, then a footer: the grouping must read
# at a glance, which an overlap assertion cannot check for us.
lines = [
    (0.745, 'Test a substance more often', GY, 'normal'),
    (0.655, 'and its "protective" value falls.', GY, 'normal'),
    (0.505, '40–48% is sampling depth,', TE, 'bold'),
    (0.415, 'not chemistry.', TE, 'bold'),
    (0.265, 'Fix the endpoint: the model transfers,', GY, 'normal'),
    (0.175, 'yet no untested PFAS clears.', GY, 'normal'),
]
objs = [fig.text(bx, y, t, fontsize=7.1, color=c, va='top', ha='left', fontweight=w)
        for y, t, c, w in lines]
foot = fig.text(bx, 0.065, '87 PFAS · 4,744 ECOTOX records', fontsize=6.3,
                color='#8A8A8A', va='top', ha='left')
fig.canvas.draw()
CWtmp = fig.canvas.get_width_height()[0]
for o in objs + [foot]:              # shrink any line that crosses the right margin
    while o.get_window_extent().x1 > RIGHT * CWtmp and o.get_fontsize() > 4.5:
        o.set_fontsize(o.get_fontsize() - 0.15)
        fig.canvas.draw()

fig.canvas.draw()
items=[t for t in [ax.xaxis.label,ax.yaxis.label]+list(ax.texts)+list(fig.texts) if t.get_text()]
for axis,lim in ((ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())):
    for loc,lab in zip(axis.get_ticklocs(),axis.get_ticklabels()):
        if min(lim)-1e-9<=loc<=max(lim)+1e-9 and lab.get_text(): items.append(lab)
boxes=[(t.get_text()[:30],t.get_window_extent()) for t in items]
bad=[(p[0],q[0]) for i,p in enumerate(boxes) for q in boxes[i+1:] if p[1].overlaps(q[1])]
CW,CH=fig.canvas.get_width_height()
off=[b[0] for b in boxes if b[1].x0<-1 or b[1].y0<-1 or b[1].x1>CW+1 or b[1].y1>CH+1]
print('text items :',len(boxes)); print('overlaps   :',bad if bad else 'none')
print('off-canvas :',off if off else 'none')
assert not bad and not off,'LAYOUT FAILURE'
for ext in ('png','tif','pdf'):
    fig.savefig(f'graphical_abstract.{ext}',dpi=DPI)
print(f'canvas {CW} x {CH} px   (minimum width 1328)')
assert CW>=1328, 'TOO NARROW FOR ELSEVIER'
print('wrote graphical_abstract.png / .tif / .pdf')
