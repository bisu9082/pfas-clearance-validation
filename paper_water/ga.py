import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import pandas as pd
plt.rcParams.update({"font.family":"DejaVu Sans"})
BLUE="#2c5f8a"; LBLUE="#9ec5e0"; RED="#b4442e"; GREY="#6b7280"; INK="#1f2937"
fig=plt.figure(figsize=(16.0,6.4),dpi=100)
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,160); ax.set_ylim(0,64); ax.axis("off")
ax.add_patch(Rectangle((0,0),160,64,fc="white",ec="none"))

def box(x,y,w,h,fc,ec,lw=1.4,r=1.4,a=1.0):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0,rounding_size={r}",
                 fc=fc,ec=ec,lw=lw,alpha=a,zorder=2))
def arrow(x1,y1,x2,y2,c=GREY,lw=2.4):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=22,
                 color=c,lw=lw,zorder=3,shrinkA=0,shrinkB=0))

ax.text(80,60.0,"One toxicity database, one monitoring record, two defensible conventions",
        ha="center",va="center",fontsize=17.5,color=INK,weight="bold")

# ---------- LEFT: the shared inputs
box(3,12,34,40,"#f4f6f8",GREY,1.2)
ax.text(20,48.4,"SAME INPUTS",ha="center",fontsize=11.5,color=GREY,weight="bold")
ax.text(20,43.2,"4,744",ha="center",fontsize=25,color=BLUE,weight="bold")
ax.text(20,39.6,"curated ECOTOX records",ha="center",fontsize=10.6,color=INK)
ax.text(20,36.6,"87 PFAS · 151 species",ha="center",fontsize=9.6,color=GREY)
ax.plot([8,32],[33.6,33.6],c="#d6dbe1",lw=1.2)
ax.text(20,29.0,"20,682",ha="center",fontsize=25,color=BLUE,weight="bold")
ax.text(20,25.4,"measured concentrations",ha="center",fontsize=10.6,color=INK)
ax.text(20,22.4,"1,180 US stream and lake",ha="center",fontsize=9.6,color=GREY)
ax.text(20,19.8,"stations, 2016–2026",ha="center",fontsize=9.6,color=GREY)
ax.text(20,15.4,"22 PFAS measured",ha="center",fontsize=9.6,color=GREY)

arrow(38,32,45,32)

# ---------- MIDDLE: the fork
ax.text(66,51.5,"THE ONLY DIFFERENCE",ha="center",fontsize=11.5,color=GREY,weight="bold")
box(46,34.5,40,13.5,"#fdf1ee",RED,1.6)
ax.text(66,44.2,"anchor on an extreme",ha="center",fontsize=12.4,color=RED,weight="bold")
ax.text(66,40.6,"lowest value · lowest NOEC · SSD 5th pct",ha="center",fontsize=9.5,color=INK)
ax.text(66,37.2,"what every regulator does",ha="center",fontsize=9.2,color=GREY,style="italic")

box(46,16.5,40,13.5,"#eef5ee","#2f7d4f",1.6)
ax.text(66,26.2,"anchor on the centre",ha="center",fontsize=12.4,color="#2f7d4f",weight="bold")
ax.text(66,22.6,"median acute, species-resolved",ha="center",fontsize=9.5,color=INK)
ax.text(66,19.2,"no detectable effort bias, but not protective",ha="center",fontsize=9.2,color=GREY,style="italic")

arrow(87,41,95,45,RED); arrow(87,23,95,19,"#2f7d4f")

# ---------- RIGHT: the verdicts
box(95,34.5,62,17.0,"#fdf1ee",RED,1.6)
ax.text(100,47.2,"76.5%",ha="left",fontsize=30,color=RED,weight="bold")
ax.text(126,47.6,"of measurements exceed",ha="left",fontsize=11.6,color=INK,va="center")
ax.text(126,44.2,"99.3% of the 1,180 stations",ha="left",fontsize=10.2,color=GREY,va="center")
ax.text(100,38.2,"median criterion 0.05 ng L$^{-1}$ — below the 1st percentile of detected concentrations",
        ha="left",fontsize=9.3,color=GREY)

box(95,13.0,62,17.0,"#eef5ee","#2f7d4f",1.6)
ax.text(100,25.7,"0.11%",ha="left",fontsize=30,color="#2f7d4f",weight="bold")
ax.text(128,26.1,"of measurements exceed",ha="left",fontsize=11.6,color=INK,va="center")
ax.text(128,22.7,"1.2% of the 1,180 stations",ha="left",fontsize=10.2,color=GREY,va="center")
ax.text(100,16.7,"median criterion 296,556 ng L$^{-1}$ — above every official aquatic-life value",
        ha="left",fontsize=9.3,color=GREY)

# ---------- bottom strip: the finding
box(3,2.0,154,8.2,"#f4f6f8",INK,1.3)
ax.text(80,7.6,"Reproducing the full regulatory screen reaches the official criteria "
        "(PFOS 4,738 vs Canada 6,800 ng L$^{-1}$) and leaves the bias intact: "
        r"$\rho$ = $-$0.69",
        ha="center",fontsize=11.6,color=INK)
ax.text(80,4.1,"Across 4,757 ECOTOX chemicals of all classes, equalising sampling depth removes half of the "
        "minimum's bias and none of the median's.",
        ha="center",fontsize=11.0,color=RED)
plt.savefig("graphical_abstract.png",dpi=100)
plt.savefig("graphical_abstract.tif",dpi=300,pil_kwargs={"compression":"tiff_lzw"})
print("saved")
