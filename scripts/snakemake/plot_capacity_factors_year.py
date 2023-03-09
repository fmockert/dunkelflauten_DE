"""Creates plots regarding capacity factors of wind and solar."""

__author__      = "Fabian Mockert, Fabian Neumann"
__copyright__   = "Copyright 2023, Fabian Mockert, Fabian Neumann"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pylab
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

def export_legend(legend):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(snakemake.input.legend, dpi="figure", bbox_inches=bbox)

def plot_cf(clim_dict, tech_colors, tech_style, include):
  ys, ye = snakemake.config["years"]
  
  fig, ax = plt.subplots(figsize=(20,7))
  width=2
  if 'offwind' in include: ax.plot(clim_dict['offwind']['mean'].values, c=tech_colors['offwind'], linestyle=tech_style['offwind'], label=f"Offshore wind", linewidth=width, zorder=3)
  if 'onwind' in include: ax.plot(clim_dict['onwind']['mean'].values, c=tech_colors['onwind'], linestyle=tech_style['onwind'], label=f"Onshore wind", linewidth=width, zorder=3)
  if 'solar' in include: ax.plot(clim_dict['solar']['mean'].values, c=tech_colors['solar'], linestyle=tech_style['solar'], label=f"Solar", linewidth=width, zorder=3)
  if 'combined' in include: ax.plot(clim_dict['combined']['mean'].values, c=tech_colors['combined'], linestyle=tech_style['combined'], label=f"Combined", linewidth=width, zorder=3)
  ax.plot(cf_to_annual_mean(cf_tech['offwind'].rolling(24*7, min_periods=0, center=True).mean()).values, c=tech_colors['offwind'], zorder=1, linestyle='dotted')
  ax.plot(cf_to_annual_mean(cf_tech['onwind'].rolling(24*7, min_periods=0, center=True).mean()).values, c=tech_colors['onwind'], zorder=1, linestyle='dotted')
  ax.plot(cf_to_annual_mean(cf_tech['solar'].rolling(24*7, min_periods=0, center=True).mean()).values, c=tech_colors['solar'], zorder=1, linestyle='dotted')
  ax.plot(cf_to_annual_mean(cf_tech['combined'].rolling(24*7, min_periods=0, center=True).mean()).values, c=tech_colors['combined'], zorder=1, linestyle='dotted')
  
  if 'minmax' in include: ax.fill_between(np.arange(0,366*24,1), clim_dict['combined']['min'].values.reshape(366*24), clim_dict['combined']['max'].values.reshape(366*24), color=tech_colors['range'], alpha=0.4, label=f"Combined range")
  if 'threshold' in include: ax.hlines(0.05, 0, 366*24, color='red', linestyle='dashed', label=f"Threshold")
  
  ax.set_xlabel('Months')
  ax.set_ylabel('Capacity factor [ratio]')
  
  x_labels = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
  x_ticks = [i*24 for i in [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]]
  ax.set_xticks(x_ticks)
  ax.set_xticklabels(x_labels)
  ax.set_ylim(0,0.60)
  ax.set_xlim(0,366*24)
  ax.grid(True)
  handles, labels = ax.get_legend_handles_labels()
  ax.legend(handles, labels, ncol=3, bbox_to_anchor=(0.95,-0.2))
  plt.savefig(snakemake.output.cf_means, bbox_inches='tight')
  #Graphic as text
  text = open('%s.txt'%(snakemake.output.cf_means[:-4]), mode='w')
  for tech in tech_style.keys():
    text.write("%s\n"%(tech))
    text.write("year: min: %s, mean: %s, max: %s\n"%(clim_dict[tech]['min']['DE'].min(), clim_dict[tech]['mean']['DE'].mean(), clim_dict[tech]['max']['DE'].max()))
    text.write("DJF: min: %s, mean: %s, max: %s\n"%(clim_dict[tech]['min']['DE'][[12,1,2]].min(), clim_dict[tech]['mean']['DE'][[12,1,2]].mean(), clim_dict[tech]['max']['DE'][[12,1,2]].max()))
    text.write("JJA: min: %s, mean: %s, max: %s\n\n"%(clim_dict[tech]['min']['DE'][[6,7,8]].min(), clim_dict[tech]['mean']['DE'][[6,7,8]].mean(), clim_dict[tech]['max']['DE'][[6,7,8]].max()))
  text.close()

def cf_to_annual_min(df):
    return df.groupby([df.index.month,df.index.day]).min()
def cf_to_annual_mean(df):
    return df.groupby([df.index.month,df.index.day]).mean()
def cf_to_annual_max(df):
    return df.groupby([df.index.month,df.index.day]).max()    

if __name__ == "__main__":

  countries = ['DE']
  # load data from atlite
  cf_tech = {'onwind': pd.read_csv(snakemake.input.onwind, index_col=0, parse_dates=True),
              'offwind': pd.read_csv(snakemake.input.offwind, index_col=0, parse_dates=True),
              'solar': pd.read_csv(snakemake.input.solar, index_col=0, parse_dates=True),
              'combined': pd.read_csv(snakemake.input.combined, index_col=0, parse_dates=True)}
  
  # 48 hour rolling mean capacity factors
  clim_dict = {key:{subkey: None for subkey in ['min', 'mean', 'max']} for key in ['onwind', 'offwind', 'solar', 'combined']}
  for key in clim_dict.keys():
      clim_dict[key]['min'] = cf_to_annual_min(cf_tech[key].rolling(48, min_periods=0, center=True).mean())
      clim_dict[key]['mean'] = cf_to_annual_mean(cf_tech[key].rolling(48, min_periods=0, center=True).mean())
      clim_dict[key]['max'] = cf_to_annual_max(cf_tech[key].rolling(48, min_periods=0, center=True).mean())
  
  onwind_w = snakemake.config['technology']['onwind']['weight']
  offwind_w = snakemake.config['technology']['offwind']['weight']
  solar_w = snakemake.config['technology']['solar']['weight']
  tech_colors = {"solar": "gold", "onwind": "royalblue", "offwind": "#6823bc",
                  "wind": "#2393bc", "combined": "black", "range": "palegreen"}
  tech_style = {"solar": 'solid', "onwind": 'solid', "offwind": 'solid', "combined": (0, (7, 1))}

  name_split = snakemake.output.cf_means.split('/')
  n=name_split[0]
  for i in range(1, len(name_split)):
    if os.path.isdir(n)==False: os.system('mkdir %s'%(n))
    n+='/%s'%(name_split[i])

  plot_cf(clim_dict, tech_colors=tech_colors, tech_style=tech_style, include=['threshold', 'onwind', 'offwind', 'solar', 'combined', 'minmax'])