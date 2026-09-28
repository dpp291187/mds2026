"""Regenerate publication figures from exact JSON experiment data."""
from pathlib import Path
import json, statistics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'data/extended_results.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'legend.fontsize':8,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
names=[('HorizenLabs_M4','HL','#0072B2','o'),('Plonky3_M4','P3','#D55E00','s'),('Signed_A49','Signed (12 ops)','#7A3E9D','^')]
fig,axs=plt.subplots(1,3,figsize=(7.3,2.35),layout='constrained',sharey=True)
for ax,(n,label,color,m) in zip(axs,names):
 z=[t for t in r['order_experiments'] if t['circuit']==n]
 hs=[t['H'] for t in z]
 for field,style,tag,alpha in [('peak_states','-','Exact merging',.12),('dominance_peak','--','+ dominance',.06)]:
  med=[statistics.median(a[field] for a in t['records']) for t in z]
  low=[min(a[field] for a in t['records']) for t in z]
  high=[max(a[field] for a in t['records']) for t in z]
  ax.plot(hs,med,color=color,marker=m,ms=3,lw=1.5,ls=style,label=tag)
  ax.fill_between(hs,low,high,color=color,alpha=alpha)
 ax.set(title=label,xlabel='Capacity H',xticks=[3,5,8],yscale='log');ax.grid(axis='y',alpha=.2)
axs[0].set_ylabel('Peak live states');axs[1].legend(frameon=False,loc='upper left',fontsize=7)
fig.savefig(ROOT/'figures/state_counts.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(3.55,2.7),layout='constrained')
data=[[a['peak_states'] for a in next(t for t in r['order_experiments'] if t['circuit']==n and t['H']==5)['records']] for n,_,_,_ in names]
b=ax.boxplot(data,tick_labels=[n[1] for n in names],patch_artist=True,widths=.45,showfliers=True,medianprops={'color':'black','lw':1.2})
for p,(_,_,color,_) in zip(b['boxes'],names):p.set_facecolor(color);p.set_alpha(.42)
for i,((n,label,color,m),vals) in enumerate(zip(names,data),1):
 base=next(t['plain_peak'] for t in r['grid'] if t['circuit']==n and t['H']==5)
 ax.plot(i,base,'D',color='black',ms=4,label='Original order' if i==1 else None)
ax.set(ylabel='Peak live states',xlabel='100 topological orders per circuit (H = 5)',ylim=(0,480));ax.grid(axis='y',alpha=.2);ax.legend(frameon=False,loc='upper left')
fig.savefig(ROOT/'figures/order_states.pdf');plt.close(fig)
print('Generated two supplementary order-statistics figures.')
