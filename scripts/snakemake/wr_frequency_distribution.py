import pandas as pd
import numpy as np
import matplotlib.pylab as pylab
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib as mpl
import os
from matplotlib.colors import ListedColormap
# Size configuration
size = 20
params = {'legend.fontsize': size,
          'figure.figsize': (10, 10),
         'axes.labelsize': size,
         'axes.titlesize': size,
         'xtick.labelsize': size,
         'ytick.labelsize': size}
pylab.rcParams.update(params)


def plot_distribution_heatmap(distribution, name):
  if os.path.isdir('../../outputs/plots/WR')==False: os.system('mkdir ../../outputs/plots/WR')
  fig,ax=plt.subplots()
  vmin=5; vcenter=10; vmax=15; cnum=10
  heat = sns.heatmap(distribution, cmap=sns.color_palette("seismic",cnum), norm=colors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax), annot=True, annot_kws={'size':size+10}, fmt='.1f', linewidths=1, 
              square=True, cbar=False, cbar_kws={'extend':'both', 'label':'Frequency [%]'})
  ax.set_xlabel('Weather regimes')
  ax.set_ylabel('Season')
  #ax.set_title('%s weather regimes'%(name))
  ax.set_yticklabels(ax.get_yticklabels(), va = "center")
  plt.savefig('%s_%s.pdf'%(snakemake.output.wr_dis_plot[:-13], name), bbox_inches='tight')
  #Colorbar
  fig2, ax2 = plt.subplots(figsize=(6, 0.5))
  fig2.subplots_adjust(bottom=0.5)
  cb1 = mpl.colorbar.ColorbarBase(ax2, extend='both', cmap=ListedColormap(sns.color_palette("seismic",cnum).as_hex()), norm=colors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax), orientation='horizontal')
  cb1.set_label('Frequency [%]')
  plt.savefig(snakemake.output.wr_dis_plot, bbox_inches='tight')


#Load DF data
dunkelflauten = pd.read_csv(snakemake.input.dunkelflauten, index_col=[0], parse_dates=[1,2,4,5])
dunkelflauten_times = pd.read_csv(snakemake.input.dunkelflauten_times, parse_dates=[0], header=None)
dunkelflauten_times.index = pd.to_datetime(dunkelflauten_times[0], format='%Y%m%d_%H').values

#load WR index
#ERA5
iwr = pd.read_csv(snakemake.input.iwr, sep=',', index_col=[0], parse_dates=[0]).squeeze()

#WR all year
wr_dis_year = iwr.value_counts()/len(iwr)*100

#WR SON
wr_dis_son = iwr[iwr.index.month.isin([9,10,11])].value_counts()/len(iwr[iwr.index.month.isin([9,10,11])])*100

#WR DJF
wr_dis_djf = iwr[iwr.index.month.isin([12,1,2])].value_counts()/len(iwr[iwr.index.month.isin([12,1,2])])*100

#WR NDJFM
wr_dis_ndjfm = iwr[iwr.index.month.isin([11,12,1,2,3])].value_counts()/len(iwr[iwr.index.month.isin([11,12,1,2,3])])*100

#WR in DF
wr_dis_df = iwr[iwr.index.isin(dunkelflauten_times.index)].value_counts()/len(iwr[iwr.index.isin(dunkelflauten_times.index)])*100

#WR in DF in SON
wr_dis_df_son = iwr[(iwr.index.isin(dunkelflauten_times.index)) & (iwr.index.month.isin([9,10,11]))].value_counts()/len(iwr[(iwr.index.isin(dunkelflauten_times.index)) & (iwr.index.month.isin([9,10,11]))])*100

#WR in DF in DJF
wr_dis_df_djf = iwr[(iwr.index.isin(dunkelflauten_times.index)) & (iwr.index.month.isin([12,1,2]))].value_counts()/len(iwr[(iwr.index.isin(dunkelflauten_times.index)) & (iwr.index.month.isin([12,1,2]))])*100

wr_dis = pd.DataFrame([wr_dis_year, wr_dis_df, wr_dis_son, wr_dis_df_son, wr_dis_djf, wr_dis_df_djf], index=[['All days', 'Dunkelflauten days', 'All days', 'Dunkelflauten days', 'All days', 'Dunkelflauten days'],['Year',  'Year', 'SON', 'SON', 'DJF', 'DJF']])

regime_names = ['AT', 'ZO', 'ScTr', 'AR', 'EuBL', 'ScBL', 'GL', 'No']
wr_dis = wr_dis[regime_names]
wr_dis_table = open(snakemake.output.wr_dis_table, 'w')
wr_dis_table.write('\multicolumn{1}{r}{Day selection} & {Season} ')
for wr_i in regime_names: wr_dis_table.write('& {%s}'%(wr_i))
wr_dis_table.write('\\\\\n\midrule\n')
for category in wr_dis.index:
  wr_dis_table.write('\multicolumn{1}{r}{%s} & {%s}'%(category[0], category[1]))
  for wr_i in regime_names: wr_dis_table.write(' & %s'%(np.round(wr_dis.loc[category, wr_i], 1)))
  wr_dis_table.write('\\\\\n')
wr_dis_table.close()

all_wr = wr_dis.xs('All days')
df_wr = wr_dis.xs('Dunkelflauten days')


plot_distribution_heatmap(distribution=all_wr, name='All',)
plot_distribution_heatmap(distribution=df_wr, name='Dunkelflauten')