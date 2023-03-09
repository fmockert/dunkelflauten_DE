"""Plot histograms with the duration of Dunkelflauten, associated weather regimes, Dunkelflauten per winter season"""

__author__ = "Fabian Mockert"
__copyright__ = "Copyright 2023, Fabian Mockert"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import matplotlib
import os

#Size configuration
size = 23
params = {'legend.fontsize': size, 'figure.figsize': (15, 5), 'axes.labelsize': size, 'axes.titlesize': size,
         'font.size': size, 'xtick.labelsize': size, 'ytick.labelsize': size}
pylab.rcParams.update(params)

# Separate legend
def export_legend(legend, filename):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig('%s_legend.png'%(filename), dpi="figure", bbox_inches='tight')

# Create a climatology
def create_climatology(input_data):
  clim = input_data.groupby([input_data.index.month, input_data.index.day, input_data.index.hour]).mean()
  clim_series = pd.Series(np.nan, index=input_data.index)

  for time in input_data.index:
    clim_series[time] = clim.loc[time.month, time.day, time.hour]['actual']
  return(clim_series)

# Plot Dunkelflauten per month with length indicated
def df_len_plot(dict):
  fig,ax = plt.subplots(figsize=(15,5))
  bottom = [0,0,0,0,0,0,0,0,0,0,0,0]
  label_ = False
  for month in months:
    for day in days:
      if label_==False:
        ax.bar(months.index(month), dict[month][day], bottom=bottom[months.index(month)], label='<%s days'%(day), color=day_colour[days.index(day)]) 
      else: ax.bar(months.index(month), dict[month][day], bottom=bottom[months.index(month)], color=day_colour[days.index(day)])
      bottom[months.index(month)] += dict[month][day]
    label_=True
  ax.set_xticks(np.arange(0,12,1))
  ax.set_xlim(-0.5,11.5)
  ax.set_xticklabels(['Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov'])
  ax.set_xlabel('Months')
  ax.set_yticks(np.arange(0,50,10))
  ax.set_ylabel('Number of Dunkelflauten')
  ax.grid(axis='y')
  ax.set_axisbelow(True)
  handles, labels = ax.get_legend_handles_labels()
  name_split = snakemake.output.df_amount_len.split('/')
  n=name_split[0]
  for i in range(1, len(name_split)):
    if os.path.isdir(n)==False: os.system('mkdir %s'%(n))
    n+='/%s'%(name_split[i])

  ax.legend(handles[::-1], labels[::-1], ncol=4, bbox_to_anchor=(0.97, -0.2)) #title='Length')
  plt.savefig(snakemake.output.df_amount_len, bbox_inches='tight')

  #Graphic as text
  text = open('%s.txt'%(snakemake.output.df_amount_len[:-4]), mode='w')
  for month in df_len.keys():
    text.writelines(['Month: ', str(month), ' [max days:amount]: ', str(df_len[month]), '\n'])
  text.close()

# Plot Dunkelflauten per month with regime indicated
def df_regime_plot(dict):
  fig,ax = plt.subplots(figsize=(15,5))
  bottom = [0,0,0,0,0,0,0,0,0,0,0,0]
  label_ = False
  r= [regime_name[i] for i in regime_order]
  r.insert(0, r.pop())
  for month in months:
    for regime in r:
      if label_==False:
        ax.bar(months.index(month), dict[month][regime], bottom=bottom[months.index(month)], label=regime, color=color_dict[regime]['hex'])  
      else: ax.bar(months.index(month), dict[month][regime], bottom=bottom[months.index(month)], color=color_dict[regime]['hex']) 
      bottom[months.index(month)] += dict[month][regime]
    label_=True
  ax.legend(bbox_to_anchor=(1,1))
  ax.set_xticks(np.arange(0,12,1))
  ax.set_xlim(-0.5,11.5)
  ax.set_xticklabels(['Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov'])
  ax.set_xlabel('Months')
  ax.set_yticks(np.arange(0,50,10))
  ax.set_ylabel('Number of Dunkelflauten')
  ax.grid(axis='y')
  ax.set_axisbelow(True)
  handles, labels = ax.get_legend_handles_labels()
  labels = [labels[i] for i in [0,4,1,5,2,6,3,7]]
  handles = [handles[i] for i in [0,4,1,5,2,6,3,7]]
  ax.legend(handles, labels, ncol=4, bbox_to_anchor=(0.9, -0.2)) #title='Length')
  plt.savefig(snakemake.output.df_amount_regime, bbox_inches='tight')

  #Graphic as text
  text = open('%s.txt'%(snakemake.output.df_amount_regime[:-4]), mode='w')
  for month in df_len.keys():
    text.writelines(['Month: ', str(month), ' [regime:amount]: ', str(df_regime[month]), '\n'])
  text.close()

# Plot Dunkelflauten per year/season with length indicated
def df_len_over_season_plot(dict):
  fig,ax = plt.subplots(figsize=(15,5))
  bottom = [0]*len(dict.keys())
  label_ = False
  for year in np.arange(year_start, year_end+1,1):
    for day in [3,4,5,6,7,8,9]:
      if label_==False:
        ax.bar(years.index(year), dict[year][day], bottom=bottom[years.index(year)], label=r'$\leq$%s days'%(day), color=day_colour[days.index(day)])
      else: ax.bar(years.index(year), dict[year][day], bottom=bottom[years.index(year)], color=day_colour[days.index(day)])
      bottom[years.index(year)] += dict[year][day]
    label_=True
  ax.set_xticks(np.arange(0, len(years) ,5))
  ax.set_xticklabels(['%s/%s'%(str(x)[2:], str(x+1)[2:]) for x in np.arange(year_start, year_end,5)], rotation=45, ha='right')
  ax.set_xlabel('Extended winter periods')
  ax.set_yticks(np.arange(0,11,2))
  ax.set_ylabel('Number of Dunkelflauten')
  ax.set_xlim(-0.5, len(years)-0.5)
  ax.grid(axis='y')
  ax.set_axisbelow(True)
  handles, labels = ax.get_legend_handles_labels()
  ax.legend([handles[i] for i in [0,4,1,5,2,6,3]], [labels[i] for i in [0,4,1,5,2,6,3]], ncol=4, bbox_to_anchor=(0.97, -0.4)) #title='Length')

  plt.savefig(snakemake.output.df_amount_len_season, bbox_inches='tight')

  # Graphic as text
  text = open('%s.txt'%(snakemake.output.df_amount_len_season[:-4]), mode='w')
  for season_start in df_len_season.keys():
    text.writelines(['Season: ', '%s/%s'%(season_start, season_start+1), ' [max days:amount]: ', str(df_len_season[season_start]), '\n'])
  text.close()

# Plot histograms of the Dunkelflauten duration including Box and Whisker plots
def plot_distribution_dunkelflauten(dunkelflauten):
  sns.set_style("whitegrid")
  sns.set_context(rc=params)
  
  g = sns.JointGrid(xlim=(48, 204), ylim=(0,11))				
  g.ax_marg_y.remove()
  g.fig.set_figwidth(20)
  g.fig.set_figheight(10)
  g.ax_joint.set_xticks(np.arange(48,205,12))
  
  # Main KDE plot
  main = sns.histplot(data=dunkelflauten[dunkelflauten.regime.isin(['EuBL', 'ScBL', 'GL'])].rename(columns={'regime':'Regime'}), x='hours', hue='Regime', 
                 binwidth=12, binrange=[48,204], multiple='dodge', palette={x:color_dict[x]['hex'] for x in ['EuBL', 'ScBL', 'GL']},
                 kde=True, ax=g.ax_joint, legend=False)
  
  # X marginal
  d_all = dunkelflauten.copy()
  d_all['regime']='all'
  d_marg = pd.concat([dunkelflauten, d_all])
  
  medianprops = dict(linestyle='-', linewidth=2.5, color='firebrick')
  meanpointprops = dict(marker='D', markeredgecolor='black',markerfacecolor='#FF6EB4')
  boxprops = dict(linestyle='-', linewidth=3, color='black')
  flierprops = dict(marker='x', markerfacecolor='none', markersize=12,
                    linestyle='none')
  custom_lines = [Line2D([0], [0], color=color_dict['EuBL']['hex'], label='EuBL', lw=8),
                  Line2D([0], [0], color=color_dict['ScBL']['hex'], label='ScBL', lw=8),
                  Line2D([0], [0], color=color_dict['GL']['hex'], label='GL', lw=8),
                  Line2D([0], [0], color=color_dict['all']['hex'], label='All', lw=8)]
  main.legend(handles=custom_lines, bbox_to_anchor=(1,1))
  marginal_x = sns.boxplot(data=d_marg[d_marg.regime.isin(['all', 'EuBL', 'ScBL', 'GL'])], x='hours', y='regime', 
                           ax=g.ax_marg_x, orient='h', palette={x:color_dict[x]['hex'] for x in ['all', 'EuBL', 'ScBL', 'GL']},
                           whis=[5,95], showmeans=True, meanprops=meanpointprops, medianprops=medianprops, 
                           flierprops=flierprops)
  
  g.set_axis_labels(xlabel='Dunkelflauten length [hours]', ylabel='Number of Dunkelflauten')
  g.savefig(snakemake.output.df_distribution, bbox_inches='tight')

if __name__ == "__main__":

  regime_name = ['No', 'AT', 'GL', 'AR', 'ScTr', 'EuBL', 'ScBL', 'ZO']
  regime_order=[1,7,4,3,5,6,2,0]
  year_start, year_end = snakemake.config["years"]
  years = np.arange(year_start, year_end+1,1).tolist()
  # Load dunkelflauten
  dunkelflauten = pd.read_csv(snakemake.input.dunkelflauten, index_col=[0], parse_dates=[1,2,4,5])
  
  color_dict = {'No':{'old':'#6B6B6B', 'rgba':'666666ff', 'hex':'#666666', 'opacity':1, 'rgb':[102, 102, 102]}, #ff means opacity 1 (255), rest of rgba=hex
                'AT':{'old':'#6100B3', 'rgba':'551b8aff', 'hex':'#551b8a', 'opacity':1, 'rgb':[85, 27, 138]}, 
                'GL':{'old':'#0000FE', 'rgba':'0000ffff', 'hex':'#0000ff', 'opacity':1, 'rgb':[0, 0, 255]}, 
                'AR':{'old':'#FECF0A', 'rgba':'ffd700ff', 'hex':'#ffd700', 'opacity':1, 'rgb':[255, 215, 0]}, 
                'ScTr':{'old':'#FB6207', 'rgba':'cd8500ff', 'hex':'#cd8500', 'opacity':1, 'rgb':[205, 133, 0]}, 
                'EuBL':{'old':'#117B00', 'rgba':'9acd32ff', 'hex':'#9acd32', 'opacity':1, 'rgb':[154, 205, 50]}, 
                'ScBL':{'old':'#0B5300', 'rgba':'006400ff', 'hex':'#006400', 'opacity':1, 'rgb':[0, 100, 0]}, 
                'ZO':{'old':'#FB0005', 'rgba':'cd2626ff', 'hex':'#cd2626', 'opacity':1, 'rgb':[205, 38, 38]}, 
                'all':{'old':'#6B6B6B', 'rgba':'666666ff', 'hex':'#666666', 'opacity':1, 'rgb':[102, 102, 102]}} #all needs to be adapted to not be same as No!
  
  #Save the length of Dunkelflauten per month
  df_len = {}
  months = [12,1,2,3,4,5,6,7,8,9,10,11]
  days = [3,4,5,6,7,8,9]
  for month in months:
    df_len[month] = {}
    df_month = dunkelflauten[dunkelflauten.start.dt.month==month]	
    for day in days:
      df_len[month][day] = len(df_month[(df_month.hours>=(day-1)*24) & (df_month.hours<(day*24))])
  
  # Plot Dunkelflauten per month with length as colour
  day_colour = ['#377eb8', '#ff7f00', '#4daf4a', '#a65628', '#984ea3', '#e41a1c', '#dede00'] 
  df_len_plot(df_len)
  
  # Save the regime of Dunkelflauten per month
  df_regime = {}
  months = [12,1,2,3,4,5,6,7,8,9,10,11]
  for month in months:
    df_regime[month] = {}
    df_month = dunkelflauten[dunkelflauten.start.dt.month==month]	
    for regime in regime_name:
      df_regime[month][regime] = len(df_month[df_month.regime==regime])
  
  # Plot Dunkelflauten per month with regime as colour
  r= [regime_name[i] for i in regime_order]
  r.insert(0, r.pop())
  df_regime_plot(df_regime)
  
  # Save the length of Dunkelflauten per season
  df_len_season = {}
  season_starts = np.arange(year_start,year_end+1,1)
  df_length = [3,4,5,6,7,8,9]
  dunkelflauten['days'] = (np.ceil(dunkelflauten['hours']/24)).apply(lambda x: int(x))
  dunkelflauten.loc[dunkelflauten[dunkelflauten.hours==48].index, 'days']=3
  for season_start in season_starts:
    df_len_season[season_start] = {key: 0 for key in df_length}
    df_tmp = dunkelflauten[(dunkelflauten.start>=pd.Timestamp(year=season_start, month=7, day=1)) & (dunkelflauten.start<=pd.Timestamp(year=season_start+1, month=6, day=30))].groupby('days').count()['start']
    for index in df_tmp.index:
      df_len_season[season_start][index] = df_tmp[index]
  
  # Plot Dunkelflauten per year/season with regime as colour
  df_len_over_season_plot(df_len_season)
  plot_distribution_dunkelflauten(dunkelflauten=dunkelflauten)