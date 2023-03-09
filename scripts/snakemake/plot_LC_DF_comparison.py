"""
Plot the time period between the onset of Dunkelflauten and 
weather regimes against the decay of Dunkelflauten and weather regimes 
"""

__author__ = "Fabian Mockert"
__copyright__ = "Copyright 2023, Fabian Mockert"

import numpy as np
import pandas as pd
from datetime import timedelta
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import matplotlib.pylab as pylab
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import os

# Size configuration
size = 30
params = {'legend.fontsize': size,
          'figure.figsize': (10, 10),
         'axes.labelsize': size,
         'axes.titlesize': size,
         'xtick.labelsize': size,
         'ytick.labelsize': size}
pylab.rcParams.update(params)

def create_colormaps(rgb1, rgb2, num, name):
  r_values = np.linspace(rgb1[0]/255, rgb2[0]/255, num)
  g_values = np.linspace(rgb1[1]/255, rgb2[1]/255, num)
  b_values = np.linspace(rgb1[2]/255, rgb2[2]/255, num)
  rgb_list = np.dstack([r_values, g_values, b_values])[0,:,:].tolist()

  return(LinearSegmentedColormap.from_list('%s_colors'%(name), rgb_list, N=num))

def plot_DF_LC_comparison(onset, decay, mode):
  fig, ax = plt.subplots()
  linestyle = {"EuBL": "dashed", "ScBL": "dashed", "GL": "dashed"}
  markerstyle = {"EuBL": "o", "ScBL": "o", "GL": "o"}
  
  # Plot quadrants
  ax.text(2.5, 27.5, 'I', size=40, va='center', ha='center', fontweight='bold', color='grey')
  ax.text(-25, 27.5, 'II', size=40, va='center', ha='center', fontweight='bold', color='grey')
  ax.text(-25, -2.5, 'III', size=40, va='center', ha='center', fontweight='bold', color='grey')
  ax.text(2.5, -2.5, 'IV', size=40, va='center', ha='center', fontweight='bold', color='grey')
  
  for regime in ["EuBL", "ScBL", "GL"]:
      scatters = ax.scatter(deltat_onset[regime][regime], deltat_decay[regime][regime], marker=markerstyle[regime], color=color_dict[regime]['hex'], label=regime, zorder=2) 
      if mode=='median':
        modesv = ax.axvline(np.median(onset[regime][regime]), linestyle=linestyle[regime], color=color_dict[regime]['hex'], zorder=3)#, label=regime)
        modesh = ax.axhline(np.median(decay[regime][regime]), linestyle=linestyle[regime], color=color_dict[regime]['hex'], zorder=3)
        print(regime, ', median: Onset: ', np.median(onset[regime][regime]), ', decay: ', np.median(decay[regime][regime]))
      elif mode=='mean':
        modesv = ax.axvline(np.mean(onset[regime][regime]), linestyle=linestyle[regime], color=color_dict[regime]['hex'], zorder=3)#, label=regime)
        modesh = ax.axhline(np.mean(decay[regime][regime]), linestyle=linestyle[regime], color=color_dict[regime]['hex'], zorder=3)
        print(regime, ', mean: Onset: ', np.mean(onset[regime][regime]), ', decay: ', np.mean(decay[regime][regime]))
  
  x_limit, y_limit = ax.get_xlim(), ax.get_ylim()
  ax.set_xlim(x_limit); ax.set_ylim(y_limit)
  
  ax.legend(loc='upper left', markerscale=3, borderaxespad=0.4, handletextpad=0.1, borderpad=0.1)
  x_ticks, y_ticks = np.arange(-40,10,10), np.arange(0,50,10)
  ax.set_xticks(x_ticks)
  ax.set_yticks(y_ticks)
  ax.grid(True)
  x_tick_labels = x_ticks.tolist()
  ind_x = x_tick_labels.index(0)
  gridlines_x = ax.xaxis.get_gridlines()
  gridlines_x[ind_x].set_linewidth(5)
  y_tick_labels = y_ticks.tolist()
  ind_y = y_tick_labels.index(0)
  gridlines_y = ax.yaxis.get_gridlines()
  gridlines_y[ind_y].set_linewidth(5)

  ax.set_xlabel("DF to LC onset [days]")
  ax.set_ylabel("DF to LC decay [days]")
  ax2 = ax.twinx()
  ax2.set_ylim(ax.get_ylim()[0]*24, ax.get_ylim()[1]*24)
  ax2.set_ylabel("DF to LC decay [hours]")
  ax3 = ax.twiny()
  ax3.set_xlim(ax.get_xlim()[0]*24, ax.get_xlim()[1]*24)
  ax3.set_xlabel("DF to LC onset [hours]")

  plt.tight_layout()

  name_split = snakemake.output.plot_mean.split('/')
  n=name_split[0]
  for i in range(1, len(name_split)):
    if os.path.isdir(n)==False: os.system('mkdir %s'%(n))
    n+='/%s'%(name_split[i])

  if mode=='median': 
    plt.savefig(dir_plot_median, bbox_inches='tight')
    # Graphic as text
    text = open('%s.txt'%(dir_plot_median[:-4]), mode='w')
    for regime in ["EuBL", "ScBL", "GL"]:
      text.writelines([regime, ', median: Onset: ', str(np.median(onset[regime][regime])), ', decay: ', str(np.median(decay[regime][regime])), '\n'])
    text.close()
  elif mode=='mean': 
    plt.savefig(dir_plot_mean, bbox_inches='tight')
    # Graphic as text
    text = open('%s.txt'%(dir_plot_mean[:-4]), mode='w')
    for regime in ["EuBL", "ScBL", "GL"]:
      text.writelines([regime, ', mean: Onset: ', str(np.mean(onset[regime][regime])), ', decay: ', str(np.mean(decay[regime][regime])), '\n'])
    text.close()


if __name__ == "__main__":
   
  dir_df_lc_intervals=snakemake.input.df_lc_intervals
  dir_plot_mean=snakemake.output.plot_mean
  dir_plot_median=snakemake.output.plot_median
  
  regime_names_sorted = ["AT", "ZO", "ScTr", "AR", "EuBL", "ScBL", "GL", "No"]
  
  # Regime colours
  color_dict = {'No':{'old':'#6B6B6B', 'rgba':'666666ff', 'hex':'#666666', 'opacity':1, 'rgb':[102, 102, 102]},
                'AT':{'old':'#6100B3', 'rgba':'551b8aff', 'hex':'#551b8a', 'opacity':1, 'rgb':[85, 27, 138]}, 
                'GL':{'old':'#0000FE', 'rgba':'0000ffff', 'hex':'#0000ff', 'opacity':1, 'rgb':[0, 0, 255]}, 
                'AR':{'old':'#FECF0A', 'rgba':'ffd700ff', 'hex':'#ffd700', 'opacity':1, 'rgb':[255, 215, 0]}, 
                'ScTr':{'old':'#FB6207', 'rgba':'cd8500ff', 'hex':'#cd8500', 'opacity':1, 'rgb':[205, 133, 0]}, 
                'EuBL':{'old':'#117B00', 'rgba':'9acd32ff', 'hex':'#9acd32', 'opacity':1, 'rgb':[154, 205, 50]}, 
                'ScBL':{'old':'#0B5300', 'rgba':'006400ff', 'hex':'#006400', 'opacity':1, 'rgb':[0, 100, 0]}, 
                'ZO':{'old':'#FB0005', 'rgba':'cd2626ff', 'hex':'#cd2626', 'opacity':1, 'rgb':[205, 38, 38]}, 
                'all':{'old':'#6B6B6B', 'rgba':'666666ff', 'hex':'#666666', 'opacity':1, 'rgb':[102, 102, 102]}} 
  color_dict_pale = {'No':{'old':'#7f7f7f', 'hex':'#cfcfcf', 'opacity':1}, 
                'AT':{'old':'#ab82ff', 'hex':'#c482f9', 'opacity':1}, 
                'GL':{'old':'#b0e2ff', 'hex':'#6666ff', 'opacity':1}, 
                'AR':{'old':'#f5deb3', 'hex':'#ffff94', 'opacity':1}, 
                'ScTr':{'old':'#ffa07a', 'hex':'#ffd662', 'opacity':1}, 
                'EuBL':{'old':'#90ee90', 'hex':'#d4ff6c', 'opacity':1}, 
                'ScBL':{'old':'#9bcd9b', 'hex':'#7fcf6b', 'opacity':1}, 
                'ZO':{'old':'#fa8072', 'hex':'#ff9e87', 'opacity':1}, 
                'all':{'old':'#cccccc', 'hex':'#cfcfcf', 'opacity':1}} 
  
  # Loading new df_lc_intervals file
  df_lc_intervals = pd.read_csv(dir_df_lc_intervals, index_col=0, parse_dates=True)
  for column in ["DF_start", "DF_end", "AT_LC_onset", "AT_LC_decay", "ZO_LC_onset", "ZO_LC_decay", 
                 "ScTr_LC_onset", "ScTr_LC_decay", "AR_LC_onset", "AR_LC_decay", "EuBL_LC_onset", 
                 "EuBL_LC_decay", "ScBL_LC_onset", "ScBL_LC_decay", "GL_LC_onset", "GL_LC_decay"]:
      df_lc_intervals[column] = df_lc_intervals[column].apply(pd.to_datetime)
      
  # Configure dictionary
  deltat_onset = {"AT": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}, 
                  "ZO": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}, 
                  "ScTr": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}, 
                  "AR": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}, 
                  "EuBL": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}, 
                  "ScBL": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}, 
                  "GL": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}}
  
  deltat_decay = {"AT": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}, 
                  "ZO": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}, 
                  "ScTr": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}, 
                  "AR": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}, 
                  "EuBL": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}, 
                  "ScBL": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}, 
                  "GL": {"AT":[], "ZO":[], "ScTr":[], "AR":[], "EuBL":[], "ScBL":[], "GL":[]}}
               
  for df_regime in regime_names_sorted[0:7]:
      df_lc_intervals_regime = df_lc_intervals[df_lc_intervals["DF_regime"]==df_regime]
      df_lc_intervals_regime = df_lc_intervals_regime.reset_index(drop=True)
      for i in np.arange(0, len(df_lc_intervals_regime)):
          for lc_regime in regime_names_sorted[0:7]:
              if np.isnan(df_lc_intervals_regime.at[i, lc_regime+"_LC_num"])==False:
                  deltat_onset_s = (df_lc_intervals_regime.at[i, lc_regime+"_LC_onset"]-df_lc_intervals_regime.at[i,"DF_start"])/pd.Timedelta(hours=1)/24
                  deltat_decay_s = (df_lc_intervals_regime.at[i, lc_regime+"_LC_decay"]-df_lc_intervals_regime.at[i,"DF_end"])/pd.Timedelta(hours=1)/24
                  deltat_onset[df_regime][lc_regime].append(deltat_onset_s)
                  deltat_decay[df_regime][lc_regime].append(deltat_decay_s)
  
  # Plot comparison DF to LC with mean/median
  plot_DF_LC_comparison(onset=deltat_onset, decay=deltat_decay, mode='mean')
  plot_DF_LC_comparison(onset=deltat_onset, decay=deltat_decay, mode='median')
  
  