import numpy as np
import pandas as pd
from datetime import timedelta
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

# Size configuration
size = 30
params = {'legend.fontsize': size,
          'figure.figsize': (10, 10),
         'axes.labelsize': size,
         'axes.titlesize': size,
         'xtick.labelsize': size,
         'ytick.labelsize': size}
pylab.rcParams.update(params)

dir_df_intervals=snakemake.input.df_intervals
dir_bootstrap=snakemake.input.bootstrap
dir_subsampling=snakemake.input.subsampling
dir_plot=snakemake.output.plot_bootstrap


# Color setup
# Regime colours
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

# Regimes
regime_names_sorted = ["AT", "ZO", "ScTr", "AR", "EuBL", "ScBL", "GL", "No"]

# Data input
df_intervals = pd.read_csv(dir_df_intervals, index_col=0, parse_dates=True)
df_intervals[["start", "end", "start_wr", "end_wr"]] = df_intervals[["start", "end", "start_wr", "end_wr"]].apply(pd.to_datetime)

bootstrapLC = pd.read_csv(dir_bootstrap, index_col=0, parse_dates=True)
subsamplingLC = pd.read_csv(dir_subsampling, index_col=0, parse_dates=True)
print("reference:\n"+str(bootstrapLC.mean()))
print("Dunkelflaute:\n"+str(subsamplingLC.mean()))
#purpose = 'presentation_reference'
purpose = 'presentation_full'

# Plot
#left: subsampling, right: bootstrap
fig, ax = plt.subplots()
regimes_shown = ["EuBL", "ScBL", "GL"]

medianprops = dict(linestyle='-', linewidth=2.5, color='firebrick')
meanpointprops = dict(marker='D', markeredgecolor='black',markerfacecolor='#FF6EB4')
boxprops = dict(linestyle='-', linewidth=3, color='black')
flierprops = dict(marker='x', markerfacecolor='none', markersize=12,
                  linestyle='none')
#Reference
bp_bootstrap = ax.boxplot(bootstrapLC[regimes_shown].to_dict('list').values(), positions = [1,3,5], widths=0.9, whis=[5,95], showmeans=True, 
                   meanprops=meanpointprops, medianprops=medianprops, boxprops=boxprops, 
                   flierprops=flierprops, patch_artist=True)
#only show min and max fliers
fliers = bp_bootstrap['fliers'] 
for i in range(len(fliers)):
  box = fliers[i]
  box.set_data([[box.get_xdata()[0],box.get_xdata()[0]],[np.min(box.get_ydata()), np.max(box.get_ydata())]]) 

    #set colors
for patch, color in zip(bp_bootstrap['boxes'], [color_dict_pale['EuBL']['hex'], color_dict_pale['ScBL']['hex'], color_dict_pale['GL']['hex']]):
  patch.set_facecolor(color)

#Reference+Dunkelflauten
if purpose!='presentation_reference':
  bp_subsampling = ax.boxplot(subsamplingLC[regimes_shown].to_dict('list').values(), positions=[2,4,6], widths=0.9, whis=[5,95], showmeans=True, 
                   meanprops=meanpointprops, medianprops=medianprops, boxprops=boxprops, 
                   flierprops=flierprops, patch_artist=True)
  #only show min and max fliers
  fliers = bp_subsampling['fliers'] 
  for i in range(len(fliers)):
    box = fliers[i]
    box.set_data([[box.get_xdata()[0],box.get_xdata()[0]],[np.min(box.get_ydata()), np.max(box.get_ydata())]]) 

    #set colors
  for patch, color in zip(bp_subsampling['boxes'], [color_dict['EuBL']['hex'], color_dict['ScBL']['hex'], color_dict['GL']['hex']]):
    patch.set_facecolor(color)

  #Plot Overlap Percentile
  x_pos = 1
  for regime in ['EuBL', 'ScBL', 'GL']:
    boot_ordered = bootstrapLC[regime].sort_values(ascending=True).reset_index(drop=True) #rechts
    subsamp_ordered = subsamplingLC[regime].sort_values(ascending=False).reset_index(drop=True) #links
    difference = abs(subsamp_ordered-boot_ordered)
    #(difference.idxmin()+1)/1000
    ax.hlines(boot_ordered[difference.idxmin()], xmin=x_pos-0.5, xmax=x_pos+1.5, color='black', linestyle='dashed', zorder=1)
    ax.text(x_pos-0.25, boot_ordered[difference.idxmin()]+0.25, '%s%s'%(int(np.round(100-(difference.idxmin()+1)/10, 0)), '%'), va='center', ha='center', fontsize=15)
    x_pos+=2
#  dir_plot = '../../outputs/plots/bootstrapLC_DE_0.06_48.png'
#else: dir_plot = '../../outputs/plots/bootstrapLC_DE_0.06_48_reference.png'
ax.set_xticks([1.5, 3.5, 5.5])
ax.set_xticklabels(regimes_shown)
ax.grid(axis='y')
#ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax.set_xlabel("Life cycle")
ax.set_ylabel("Duration [days]")
ax.set_ylim(7, 30)
ax.set_xlim(0.5, 6.5)

for x,text in zip(np.arange(1,7,1), ['Reference', 'Dunkelflaute']*3):
  plt.text(x, 7.5, text, ha='center', va='center', size=13)

ax2 = ax.twinx()
ax2.set_ylim(ax.get_ylim()[0]*24, ax.get_ylim()[1]*24)
ax2.set_ylabel("Duration [hours]")
#ax.set_yticks(np.arange(50,450, 50))
plt.tight_layout()
plt.savefig(dir_plot, bbox_inches='tight')

#Graphic as text
text = open('%s.txt'%(dir_plot[:-4]), mode='w')
text.write("Reference sample:\n")
text.writelines(["min:          ", str(bootstrapLC[['EuBL', 'ScBL', 'GL']].min()).replace('\n', '; ')[:-16], "\n"])
text.writelines(["0.05quantile: ", str(bootstrapLC[['EuBL', 'ScBL', 'GL']].quantile(0.05)).replace('\n', '; ')[:-28], "\n"])
text.writelines(["0.25quantile: ", str(bootstrapLC[['EuBL', 'ScBL', 'GL']].quantile(0.25)).replace('\n', '; ')[:-28], "\n"])
text.writelines(["mean:         ", str(bootstrapLC[['EuBL', 'ScBL', 'GL']].mean()).replace('\n', '; ')[:-16], "\n"])
text.writelines(["median:       ", str(bootstrapLC[['EuBL', 'ScBL', 'GL']].median()).replace('\n', '; ')[:-16], "\n"])
text.writelines(["0.75quantile: ", str(bootstrapLC[['EuBL', 'ScBL', 'GL']].quantile(0.75)).replace('\n', '; ')[:-28], "\n"])
text.writelines(["0.95quantile: ", str(bootstrapLC[['EuBL', 'ScBL', 'GL']].quantile(0.95)).replace('\n', '; ')[:-28], "\n"])
text.writelines(["max:          ", str(bootstrapLC[['EuBL', 'ScBL', 'GL']].max()).replace('\n', '; ')[:-16], "\n"])
text.writelines("Dunkelflauten sample:\n")
text.writelines(["min:          ", str(subsamplingLC[['EuBL', 'ScBL', 'GL']].min()).replace('\n', '; ')[:-16], "\n"])
text.writelines(["0.05quantile: ", str(subsamplingLC[['EuBL', 'ScBL', 'GL']].quantile(0.05)).replace('\n', '; ')[:-28], "\n"])
text.writelines(["0.25quantile: ", str(subsamplingLC[['EuBL', 'ScBL', 'GL']].quantile(0.25)).replace('\n', '; ')[:-28], "\n"])
text.writelines(["mean:         ", str(subsamplingLC[['EuBL', 'ScBL', 'GL']].mean()).replace('\n', '; ')[:-16], "\n"])
text.writelines(["median:       ", str(subsamplingLC[['EuBL', 'ScBL', 'GL']].median()).replace('\n', '; ')[:-16], "\n"])
text.writelines(["0.75quantile: ", str(subsamplingLC[['EuBL', 'ScBL', 'GL']].quantile(0.75)).replace('\n', '; ')[:-28], "\n"])
text.writelines(["0.95quantile: ", str(subsamplingLC[['EuBL', 'ScBL', 'GL']].quantile(0.95)).replace('\n', '; ')[:-28], "\n"])
text.writelines(["max:          ", str(subsamplingLC[['EuBL', 'ScBL', 'GL']].max()).replace('\n', '; ')[:-16], "\n"])
text.close()
# Confidence
for regime in regimes_shown:
    a=bootstrapLC[regime].sort_values(ascending=True).reset_index(drop=True)
    b=subsamplingLC[regime].sort_values(ascending=True).reset_index(drop=True)
    meet = 1000
    for i in np.arange(0,10):#00):
        dist = a[9-i]-b[i]
        #dist = a[999-i]-b[i]
        if np.abs(dist)<=meet:
            meet_point = i
            meet = np.abs(dist)
#    print("Regime: "+regime+", confidence: "+str(100-(meet_point+1)/10)+"%")

