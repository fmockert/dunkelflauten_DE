"""
Plot the Dunkelflauten throughout the 40year period (1979-2018) as a calendar 
plot with the Dunkelflauten coloured with the dominant weather regime colour
"""

__author__ = "Fabian Mockert"
__copyright__ = "Copyright 2023, Fabian Mockert"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import matplotlib.pylab as pylab

#Size configuration
size = 23
params = {'legend.fontsize': size, 'figure.figsize': (15, 5), 'axes.labelsize': size, 'axes.titlesize': size,
         'font.size': size, 'xtick.labelsize': size, 'ytick.labelsize': size}
pylab.rcParams.update(params)

if __name__ == "__main__":

  regime_name = ['No', 'AT', 'GL', 'AR', 'ScTr', 'EuBL', 'ScBL', 'ZO']
  regions = {'regime': [-80, 40, 30 ,90]}
  color_dict = {'No':{'old':'#6B6B6B', 'rgba':'666666ff', 'hex':'#666666', 'opacity':1, 'rgb':[102, 102, 102]}, #ff means opacity 1 (255), rest of rgba=hex
                'AT':{'old':'#6100B3', 'rgba':'551b8aff', 'hex':'#551b8a', 'opacity':1, 'rgb':[85, 27, 138]}, 
                'GL':{'old':'#0000FE', 'rgba':'0000ffff', 'hex':'#0000ff', 'opacity':1, 'rgb':[0, 0, 255]}, 
                'AR':{'old':'#FECF0A', 'rgba':'ffd700ff', 'hex':'#ffd700', 'opacity':1, 'rgb':[255, 215, 0]}, 
                'ScTr':{'old':'#FB6207', 'rgba':'cd8500ff', 'hex':'#cd8500', 'opacity':1, 'rgb':[205, 133, 0]}, 
                'EuBL':{'old':'#117B00', 'rgba':'9acd32ff', 'hex':'#9acd32', 'opacity':1, 'rgb':[154, 205, 50]}, 
                'ScBL':{'old':'#0B5300', 'rgba':'006400ff', 'hex':'#006400', 'opacity':1, 'rgb':[0, 100, 0]}, 
                'ZO':{'old':'#FB0005', 'rgba':'cd2626ff', 'hex':'#cd2626', 'opacity':1, 'rgb':[205, 38, 38]}, 
                'all':{'old':'#6B6B6B', 'rgba':'666666ff', 'hex':'#666666', 'opacity':1, 'rgb':[102, 102, 102]}} #all needs to be adapted to not be same as No!
  color_dict_pale = {'No':{'old':'#7f7f7f', 'hex':'#cfcfcf', 'opacity':1}, #in rgba: last two ff means opacity 1 (255), rest of rgba=hex
                'AT':{'old':'#ab82ff', 'hex':'#c482f9', 'opacity':1}, 
                'GL':{'old':'#b0e2ff', 'hex':'#6666ff', 'opacity':1}, 
                'AR':{'old':'#f5deb3', 'hex':'#ffff94', 'opacity':1}, 
                'ScTr':{'old':'#ffa07a', 'hex':'#ffd662', 'opacity':1}, 
                'EuBL':{'old':'#90ee90', 'hex':'#d4ff6c', 'opacity':1}, 
                'ScBL':{'old':'#9bcd9b', 'hex':'#7fcf6b', 'opacity':1}, 
                'ZO':{'old':'#fa8072', 'hex':'#ff9e87', 'opacity':1}, 
                'all':{'old':'#cccccc', 'hex':'#cfcfcf', 'opacity':1}} #all needs to be adapted to not be same as No!
  
  # Load dunkelflauten
  dunkelflauten = pd.read_csv(snakemake.input.intervals, index_col=[0], parse_dates=[1,2,4,5])
  dunkelflauten_times = pd.read_csv(snakemake.input.timestamps, parse_dates=[0], header=None)
  # Load WR index
  iwr = pd.read_csv(snakemake.input.wr, index_col=[0], parse_dates=[0]).squeeze().resample('1h').pad()
  iwr_name = iwr
  
  for i in dunkelflauten.index:
    if dunkelflauten.loc[i,'start'].month in [9,10,11,12]: cal_start = pd.Timestamp('%s-09-01'%(dunkelflauten.loc[i,'start'].year))
    elif dunkelflauten.loc[i,'start'].month in [1,2,3,4,5,5]: cal_start = pd.Timestamp('%s-09-01'%(dunkelflauten.loc[i,'start'].year-1))
    else: print('Something is wrong!')
    dunkelflauten.loc[i, 'wdoy_start'] = pd.Timedelta(dunkelflauten.loc[i, 'start']-cal_start)/pd.to_timedelta(1, unit='D')
    dunkelflauten.loc[i, 'wdoy_end'] = pd.Timedelta(dunkelflauten.loc[i, 'end']-cal_start)/pd.to_timedelta(1, unit='D')
    dunkelflauten.loc[i, 'wyear'] = int(cal_start.year)
  
  year_start, year_end = snakemake.config["years"]
  dunkelflauten = dunkelflauten[dunkelflauten.start.dt.year.isin(np.arange(year_start, year_end+1))]
  fig,ax = plt.subplots(figsize=(15,10))
  for i in dunkelflauten.index:
    ax.fill_between(np.arange(dunkelflauten.loc[i, 'wdoy_start'], dunkelflauten.loc[i, 'wdoy_end']+1), dunkelflauten.loc[i, 'wyear']-0.4, dunkelflauten.loc[i, 'wyear']+0.4, color=color_dict[dunkelflauten.loc[i, 'regime']]['hex'])
  ax.vlines((pd.Timestamp('2000-01-01 00:00')-pd.Timestamp('1999-09-01 00:00')) / pd.to_timedelta(1, unit='D'), year_start-1, year_end+1, color='black', linestyle='dashed')
  
  ax.set_xticks([(pd.Timestamp(year=y, month=m, day=1, hour=0)-pd.Timestamp(year=2000, month=9, day=1, hour=0)) / pd.to_timedelta(1, unit='D') for y,m in zip([2000]*4+[2001]*5, [9,10,11,12,1,2,3,4,5])])
  ax.set_xticklabels(['Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May'])
  ax.set_xlabel('Months')
  ax.set_yticks(np.arange(1980,2025,5))
  ax.set_yticks(np.arange(1978,2025,1), minor=True)
  ax.set_yticklabels(['80/81', '85/86', '90/91', '95/96', '00/01', '05/06', '10/11', '15/16', '20/21'])
  ax.set_ylabel('Extended winter periods')
  ax.set_xlim(0,260)
  ax.set_ylim(year_start-1, year_end+1)
  
  
  ax.grid(which='major', linestyle='-')
  ax.grid(which='minor', linestyle='--')
  
  custom_lines = [
                  Line2D([0], [0], color=color_dict['AT']['hex'], label='AT', lw=8),
                  Line2D([0], [0], color=color_dict['ZO']['hex'], label='ZO', lw=8),
                  Line2D([0], [0], color=color_dict['ScTr']['hex'], label='ScTr', lw=8),
                  Line2D([0], [0], color=color_dict['AR']['hex'], label='AR', lw=8),
                  Line2D([0], [0], color=color_dict['EuBL']['hex'], label='EuBL', lw=8),
                  Line2D([0], [0], color=color_dict['ScBL']['hex'], label='ScBL', lw=8),
                  Line2D([0], [0], color=color_dict['GL']['hex'], label='GL', lw=8),
                  Line2D([0], [0], color=color_dict['No']['hex'], label='No', lw=8)]
  ax.legend(handles=custom_lines, bbox_to_anchor=(1,1))
  plt.savefig(snakemake.output.plot_calendar, bbox_inches='tight')
  