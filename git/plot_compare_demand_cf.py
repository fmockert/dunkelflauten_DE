import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import os
#Size configuration
size = 23
params = {'legend.fontsize': size, 'figure.figsize': (15, 5), 'axes.labelsize': size, 'axes.titlesize': size,
         'font.size': size, 'xtick.labelsize': size, 'ytick.labelsize': size}
pylab.rcParams.update(params)

def export_legend(legend, filename):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig('%s_legend.pdf'%(filename), dpi="figure", bbox_inches=bbox)

def create_climatology(input_data):
  clim = input_data.groupby([input_data.index.month, input_data.index.day, input_data.index.hour]).mean()
  clim_series = pd.Series(np.nan, index=input_data.index)

  for time in input_data.index:
    clim_series[time] = clim.loc[time.month, time.day, time.hour]['actual']
  return(clim_series)

#main_dir = '../../' #'~/Documents/HiWi_IAI_20/dig/' #'/home/fabian/Documents/UniKIT/Travel_Grant/Research/Reading_attempt/dunkelflauten_analysis/'

#load demand
demand_era = pd.read_csv(snakemake.input.demand, index_col=[0], parse_dates=[0])['Germany_wd_demand_no_pop_weights_no_time_trend_1979_2018.dat']
demand_era = demand_era.resample("H").ffill()
demand = demand_era.to_frame().rename(columns={'Germany_wd_demand_no_pop_weights_no_time_trend_1979_2018.dat':'actual'})
demand['climatology'] = create_climatology(input_data=demand)
demand['difference'] = demand['actual'] - demand['climatology']

#load capacity factor
capacityfactor = pd.read_csv(snakemake.input.capacityfactor, parse_dates=[0], index_col=[0]).rename(columns={'DE':'actual'})
capacityfactor['climatology']=create_climatology(input_data=capacityfactor)
capacityfactor['difference'] = capacityfactor['actual'] - capacityfactor['climatology']

#load dunkelflauten
dunkelflauten = pd.read_csv(snakemake.input.dunkelflauten, index_col=[0], parse_dates=[1,2,4,5])

#Calculate demand/capacityfactor difference to climatology for separate dunkelflauten phases
for df_num in dunkelflauten.index:
  df_start = dunkelflauten.loc[df_num, 'start']
  df_end = dunkelflauten.loc[df_num, 'end']
  dunkelflauten.loc[df_num, 'demand_diff']= demand['difference'][demand.index.isin(pd.date_range(df_start, df_end, freq='h'))].mean()
  dunkelflauten.loc[df_num, 'capacityfactor_diff']= capacityfactor['difference'][capacityfactor.index.isin(pd.date_range(df_start, df_end, freq='h'))].mean()

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

sns.set_style("whitegrid")
sns.set_context(rc=params)
#Select the dunkelflauten which shall be used for the distribution plot by regime (EuBL, ScBL, GL) and season
density_season = 'all' #'DJF' # 'SONDJF'
density_seasons = {'all':[1,2,3,4,5,6,7,8,9,10,11,12], 'DJF':[12,1,2], 'SONDJF':[9,10,11,12,1]}
dunkelflauten_plot = dunkelflauten[dunkelflauten.regime.isin(['EuBL', 'ScBL', 'GL'])]
#dunkelflauten_plot.loc[:,'capacityfactor_diff'] *= 100
dunkelflauten_plot = dunkelflauten_plot[dunkelflauten_plot.start.dt.month.isin(density_seasons[density_season])]
g = sns.JointGrid(xlim=(dunkelflauten_plot['demand_diff'].min()-0.1, dunkelflauten_plot['demand_diff'].max()+0.1),
									ylim=(dunkelflauten_plot['capacityfactor_diff'].min()-0.1, dunkelflauten_plot['capacityfactor_diff'].max()+0.1))
g.fig.set_figwidth(20)
g.fig.set_figheight(10)

#Main KDE plot
sns.kdeplot(data=dunkelflauten_plot, x='demand_diff', y='capacityfactor_diff', hue='regime',
            palette={x:color_dict[x]['hex'] for x in ['EuBL', 'ScBL', 'GL']},
						n_levels=4,
						bw_adjust=0.8,
						legend=False, #shade=True, alpha=0.3, #add shading to the 2D distribution plot
						ax=g.ax_joint
						)

#X marginal = demand
marginal_x = marginal_x = sns.kdeplot(x=dunkelflauten_plot.demand_diff, hue=dunkelflauten_plot.regime,
            palette={x:color_dict[x]['hex'] for x in ['EuBL', 'ScBL', 'GL']},
						bw_adjust=0.4, #smoothing parameter
            ax=g.ax_marg_x, legend=False, fill=True
						)

x_ticks = np.arange(-2,4,1)
x_tick_labels = x_ticks.tolist()
marginal_x.set_xticks(x_ticks)

#Y marginal = capacityfactor
marginal_y = sns.kdeplot(y=dunkelflauten_plot.capacityfactor_diff, hue=dunkelflauten_plot.regime,
            palette={x:color_dict[x]['hex'] for x in ['EuBL', 'ScBL', 'GL']},
						bw_adjust=0.4,
            ax=g.ax_marg_y, legend=False, fill=True
						)

#Main Scatterplot
#add season selection and therefore different markers
season_marker={'SON':'^', 'DJF':'o', 'MAM':'s'}
season_months={'SON':[9,10,11], 'DJF':[12,1,2], 'MAM':[3,4,5]}
scattered_regimes = ['EuBL', 'ScBL', 'GL'] #regime_names
for season in season_marker.keys():
  df_season_selection = dunkelflauten_plot[(dunkelflauten_plot['start'].dt.month.isin(season_months[season])) & (dunkelflauten_plot['regime'].isin(scattered_regimes))]
  scatt = sns.scatterplot(data=df_season_selection, x='demand_diff', y='capacityfactor_diff', hue='regime',
            palette={x:color_dict[x]['hex'] for x in scattered_regimes},
						legend=False,
						marker=season_marker[season], s=70,
						ax=g.ax_joint
						)

#Set X and Y axis labels
g.set_axis_labels(xlabel='Average electricity demand anomaly [GW]', ylabel='Average capacity factor anomaly [ratio]')
#add 0 demand line

scatt.set(xlim=(dunkelflauten['demand_diff'].min()-0.1, dunkelflauten['demand_diff'].max()+0.1),
  ylim=(dunkelflauten['capacityfactor_diff'].min()-0.01, dunkelflauten['capacityfactor_diff'].max()+0.01))

#Plot the legends
if len(scattered_regimes)==3:
	custom_lines = [Line2D([0], [0], color=color_dict['EuBL']['hex'], label='EuBL', lw=8),
                Line2D([0], [0], color=color_dict['ScBL']['hex'], label='ScBL', lw=8),
                Line2D([0], [0], color=color_dict['GL']['hex'], label='GL', lw=8)]
else:
	custom_lines = [Line2D([0], [0], color=color_dict['ZO']['hex'], label='ZO', lw=8),
                Line2D([0], [0], color=color_dict['ScBL']['hex'], label='ScBL', lw=8),
                Line2D([0], [0], color=color_dict['EuBL']['hex'], label='EuBL', lw=8),
                Line2D([0], [0], color=color_dict['ScTr']['hex'], label='ScTr', lw=8),
                Line2D([0], [0], color=color_dict['AR']['hex'], label='AR', lw=8),
                Line2D([0], [0], color=color_dict['GL']['hex'], label='GL', lw=8),
                Line2D([0], [0], color=color_dict['AT']['hex'], label='AT', lw=8),
                Line2D([0], [0], color=color_dict['No']['hex'], label='No', lw=8)]
custom_lines2 = [Line2D([0], [0], marker='^', color='w', label='SON', markerfacecolor='k', markersize=12),
                Line2D([0], [0], marker='o', color='w', label='DJF', markerfacecolor='k', markersize=12),
                Line2D([0], [0], marker='s', color='w', label='MAM', markerfacecolor='k', markersize=12)]

legend1 = scatt.legend(handles=custom_lines, bbox_to_anchor=(1,1), fontsize=20, borderaxespad=0.8, handletextpad=0.8, borderpad=0.8)
legend2 = scatt.legend(handles=custom_lines2, bbox_to_anchor=(0.8,1), fontsize=20, borderaxespad=0.8, handletextpad=0.8, borderpad=0.8)
scatt.add_artist(legend1)
scatt.add_artist(legend2)

scatt.set_xticks(x_ticks)
scatt.grid(True)

ind_x = x_tick_labels.index(0)
gridlines_x = scatt.xaxis.get_gridlines()
gridlines_x[ind_x].set_linewidth(3)
#Save the figure
if os.path.isdir('../../outputs/plots/demand')==False: os.system('mkdir ../../outputs/plots/demand')
g.savefig(snakemake.output.plot_demand, bbox_inches='tight')
#g.savefig('%s/outputs/plots/demand/demand_capacityfactor_comparison_%s_extended.png'%(main_dir, density_season), bbox_inches='tight')
#g.savefig('%s/outputs/plots/presentation/demand_capacityfactor_comparison_%s_extended.png'%(main_dir, density_season), bbox_inches='tight')

#Graphic as text
text = open('%s.txt'%(snakemake.output.plot_demand[:-4]), mode='w')

#Dunkelflauten absolute values
text.write('Absolute values\n')
for regime in ['EuBL', 'ScBL', 'GL']:
  regime_df = dunkelflauten[dunkelflauten.regime==regime]
  dates = pd.Series()
  for index in regime_df.index:
    date_range_tmp = pd.date_range(regime_df.loc[index, 'start'], regime_df.loc[index, 'end'], freq='h')
    dates = dates.append(pd.Series(date_range_tmp))
    dates = dates.reset_index(drop=True)
    mean_cf = capacityfactor[capacityfactor.index.isin(dates.values)]['actual'].mean()
    mean_demand = demand[demand.index.isin(dates.values)]['actual'].mean()
  text.writelines(['Mean CF for %s DFs: %s'%(regime, mean_cf), '. Mean demand for %s DFs: %s'%(regime, mean_demand), '\n'])

#Daily climatology absolute values
text.write('Daily Climatology absolute values\n')
for regime in ['EuBL', 'ScBL', 'GL']:
  regime_df = dunkelflauten[dunkelflauten.regime==regime]
  dates = pd.Series()
  for index in regime_df.index:
    date_range_tmp = pd.date_range(regime_df.loc[index, 'start'], regime_df.loc[index, 'end'], freq='h')
    dates = dates.append(pd.Series(date_range_tmp))
    dates = dates.reset_index(drop=True)
    mean_cf = capacityfactor[capacityfactor.index.isin(dates.values)]['climatology'].mean()
    mean_demand = demand[demand.index.isin(dates.values)]['climatology'].mean()
  text.writelines(['Mean CF for %s DFs: %s'%(regime, mean_cf), '. Mean demand for %s DFs: %s'%(regime, mean_demand), '\n'])

#Anomaly values
text.write('Anomaly values\n')
for regime in ['EuBL', 'ScBL', 'GL']:
  regime_df = dunkelflauten[dunkelflauten.regime==regime]
  dates = pd.Series()
  for index in regime_df.index:
    date_range_tmp = pd.date_range(regime_df.loc[index, 'start'], regime_df.loc[index, 'end'], freq='h')
    dates = dates.append(pd.Series(date_range_tmp))
    dates = dates.reset_index(drop=True)
    mean_cf = capacityfactor[capacityfactor.index.isin(dates.values)]['difference'].mean()
    mean_demand = demand[demand.index.isin(dates.values)]['difference'].mean()
  text.writelines(['Mean CF for %s DFs: %s'%(regime, mean_cf), '. Mean demand for %s DFs: %s'%(regime, mean_demand), '\n'])
text.writelines(['Mean NDJFM demand: %s'%(demand[demand.index.month.isin([11,12,1,2,3])]['actual'].mean())]) 
winter_demand = demand[demand.index.month.isin([11,12,1,2,3])]
text.writelines(['Min NDJFM demand: %s, Max NDJFM demand: %s'%(winter_demand['actual'].min(), winter_demand['actual'].max())])
text.writelines(['Mean NDJFM capacity factor: %s'%(capacityfactor[capacityfactor.index.month.isin([11,12,1,2,3])]['actual'].mean())])
text.close()