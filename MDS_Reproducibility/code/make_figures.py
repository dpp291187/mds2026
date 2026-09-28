"""Draw both main-paper figures from the verified data for the Springer version (PyCharm friendly).

Run this file after installing requirements-figures.txt. No LaTeX or experiment
recomputation is required. PDF, PNG (600 dpi), and SVG are written to figures/.
Supplementary traversal plots are in make_supplementary_figures.py.
"""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch,FancyBboxPatch
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'figures';OUTPUT.mkdir(exist_ok=True)
DPI=600
FORMATS=('pdf','png','svg','eps')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,
                    'legend.fontsize':8,'pdf.fonttype':42,'ps.fonttype':42,
                    'svg.fonttype':'path','text.usetex':False,
                    'axes.spines.top':False,'axes.spines.right':False})


def save(fig,name):
    for ext in FORMATS:
        fig.savefig(OUTPUT/f'{name}.{ext}',dpi=DPI,format=ext)
        if ext in ('pdf','eps'):
            number={'vertex_cover_gadget':1,'reduction_counts':2}[name]
            (ROOT/f'Fig{number}.{ext}').write_bytes((OUTPUT/f'{name}.{ext}').read_bytes())
    plt.close(fig)


def gadget():
    fig,ax=plt.subplots(figsize=(3.55,2.55))
    fig.subplots_adjust(left=.025,right=.975,bottom=.025,top=.975)
    ax.set(xlim=(-.15,6.65),ylim=(-.15,3.9));ax.set_axis_off()
    def box(cx,cy,w,h,lines):
        ax.add_patch(FancyBboxPatch((cx-w/2,cy-h/2),w,h,
                     boxstyle='round,pad=.015,rounding_size=.055',
                     linewidth=.9,edgecolor='#565AC0',facecolor='white',zorder=2))
        ax.text(cx,cy,'\n'.join(lines),ha='center',va='center',fontsize=9,linespacing=1.5,zorder=3)
    def arrow(a,b):
        ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=9,
                     linewidth=.9,color='#333333',shrinkA=1,shrinkB=1,zorder=1))
    box(1.45,3.20,2.85,.92,[r'$v_i=x_i+y_i$',r'$h_i\in\{2,1\}$'])
    box(5.05,3.20,2.85,.92,[r'$v_j=x_j+y_j$',r'$h_j\in\{2,1\}$'])
    box(3.25,1.70,3.45,.98,[r'$e_{ij}=v_i+v_j$',r'$r_{ij}=h_i+h_j\leq3$'])
    arrow((1.45,2.72),(2.72,2.21));arrow((5.05,2.72),(3.78,2.21));arrow((3.25,1.19),(3.25,.68))
    ax.text(3.25,.32,'Canonical output\none edge normalization',ha='center',va='center',fontsize=8.5,linespacing=1.35)
    return fig


def sweep():
    original=json.loads((ROOT/'data/extended_results.json').read_text())['grid']
    new=json.loads((ROOT/'data/class_witnesses.json').read_text())['grid']
    rows=original+new
    names=[('HorizenLabs_M4','HL (14 ops)','#0072B2','o','-'),
           ('Plonky3_M4','P3 (11 ops)','#D55E00','s','-'),
           ('Signed_A49',r'$A_{49}$ (12 ops)','#7A3E9D','^','--'),
           ('A7_eleven',r'$A_7,\widetilde A_7$ (11 ops)','#009E73','D',':'),
           ('A17_eleven',r'$A_{17}$ (11 ops)','#333333','x','-.')]
    fig,axes=plt.subplots(2,1,figsize=(4.95,5.0),layout='constrained')
    for name,label,col,marker,style in names:
        values=sorted((v for v in rows if v['circuit']==name),key=lambda v:v['H'])
        assert [v['H'] for v in values]==list(range(2,17))
        for ax in axes:
            ax.step([v['H'] for v in values],[v['cost'] for v in values],
                    where='post',color=col,lw=1.45,ls=style,label=label)
            ax.plot([v['H'] for v in values],[v['cost'] for v in values],marker,
                    color=col,ms=4,markerfacecolor='none' if name=='A7_eleven' else col)
    axes[0].set(xlabel='Certificate capacity H',ylabel='Minimum normalizations',
                xlim=(1.7,16.3),ylim=(3.4,14.8),xticks=[2,3,4,5,7,10,12,16],yticks=[4,6,8,10,12,14])
    axes[0].legend(frameon=False,loc='upper right',fontsize=8.0)
    axes[1].set(xlabel='Certificate capacity H',ylabel='Minimum normalizations',
                xlim=(2.85,6.15),ylim=(3.7,10.4),xticks=[3,4,5,6],yticks=list(range(4,11)))
    for ax in axes:ax.grid(axis='y',color='#dedede',alpha=1.0)
    return fig


def main():
    save(gadget(),'vertex_cover_gadget')
    save(sweep(),'reduction_counts')
    print('Generated both main-paper figures in PDF, PNG (600 dpi), SVG, and EPS:',OUTPUT)

if __name__=='__main__':main()
