#!/usr/bin/env python3
"""Reproduce the attached reference plots from saved data only."""
from pathlib import Path
from types import SimpleNamespace
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/opalx-original-model-plots')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import analytic_zero_sc as a
ROOT=Path(__file__).resolve().parent
h=SimpleNamespace(FIELDS=['emit_n_u_m','emit_n_y_m','emit_sE_m_eV'])
data=pd.read_csv(ROOT/'series.csv',float_precision='round_trip')
analytic=a.evaluate(np.linspace(0,10,1001))
station=a.evaluate(analytic.s_m,convention='station')
assert np.max(abs(analytic.emit_n_u_m-a.evaluate(analytic.s_m,128).emit_n_u_m))<1e-20
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.5,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
for kind,fields,factors,units,titles in [
 ('rms',[f'sigma_{c}_m' for c in ['u','y','s']],[1e3,1e3,1],['mm','mm','m'],[r'$\sigma_u$',r'$\sigma_y$',r'$\sigma_s$']),
 ('emittance',h.FIELDS,[1e6,1e6,1e3],['mm mrad','mm mrad','mm eV'],[r'$\epsilon_{n,u}$',r'$\epsilon_{n,y}$',r'$\epsilon_{s,E}$'])]:
 fig,axes=plt.subplots(2,3,figsize=(7.3,4.3));fig.subplots_adjust(left=.11,right=.99,top=.88,bottom=.25,wspace=.43,hspace=.35)
 for i,mode in enumerate(['ZeroSC','10pC']):
  for j,(field,factor,unit,title) in enumerate(zip(fields,factors,units,titles)):
   ax=axes[i,j]
   for code,color,marker in [('Committed 658','#0072B2','o'),('PyORBIT','#D55E00','s'),('Local original','#009E73','x')]:
    q=data[(data.code==code)&(data['mode']==mode)];values=q[field]*factor
    if kind=='emittance' and i==0 and j==2:values=values.where(values>0)
    ax.plot(q.s_m,values,color=color,marker=marker,markevery=10,ms=3,mfc='white',lw=1.1,ls='--' if code=='Local original' else '-',label=code)
   if kind=='emittance' and i==0 and j==0:ax.plot(analytic.s_m,analytic.emit_n_u_m*1e6,'k:',lw=1.5,label='Analytic: common time'
   )
   if kind=='emittance' and i==0 and j==0:ax.plot(station.s_m,station.emit_n_u_m*1e6,color='#CC79A7',ls=':',lw=1.5,label='Analytic: common station')
   maximum=float(data.loc[data['mode']==mode,field].max())*factor
   if maximum==0:maximum=float(data[field].max())*factor
   if kind=='emittance' and i==0 and j==0:maximum=max(maximum,float(analytic.emit_n_u_m.max())*1e6)
   ax.set_xlim(0,10);ax.grid(alpha=.16);ax.set_ylim(0,1.05*maximum if maximum>0 else 1)
   ax.set_ylabel((mode.replace('10pC','10 pC')+'\n' if j==0 else '')+('RMS size' if kind=='rms' else 'RMS emittance')+f' ({unit})')
   if i==0:ax.set_title(title)
   if i==1:ax.set_xlabel(r'$s$ (m)')
   if kind=='emittance' and i==0 and j==2:
    positive=data.loc[(data['mode']=='ZeroSC')&(data[field]>0),field]*factor
    ax.set_yscale('log');ax.set_ylim(positive.min()*.5,positive.max()*2);ax.text(.5,.07,'Zero values omitted',transform=ax.transAxes,ha='center',fontsize=6)
 handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.015),ncol=2,frameon=False,fontsize=7)
 fig.suptitle('$R=2$ m bend, $s=10$ m, 100 steps: local original / committed 658 / PyORBIT',fontsize=10)
 folder=ROOT
 for ext in ['png','svg']:fig.savefig(folder/f'{kind}.{ext}',dpi=300,bbox_inches='tight')
 plt.close(fig)
print(data[data.step==100].to_string(index=False))
