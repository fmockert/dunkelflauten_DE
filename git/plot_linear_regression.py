"""Plot the linear regression of the capacity factors"""

__author__ = "Fabian Mockert"
__copyright__ = "Copyright 2022, Fabian Mockert"

import pandas as pd
import numpy as np
from datetime import timedelta
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import logging
# Size configuration
size = 30; params = {'legend.fontsize': size,'figure.figsize': (10, 10), 'axes.labelsize': size, 'axes.titlesize': size, 'xtick.labelsize': size, 'ytick.labelsize': size}
import matplotlib.pylab as pylab; pylab.rcParams.update(params)

country = "DE"#snakemake.wildcards.country
source = "atlite"#snakemake.wildcards.source
carrier = snakemake.wildcards.carrier

cf = {}

opsd_raw = pd.read_csv(snakemake.input.opsd)
opsd_date_range = pd.to_datetime(
    opsd_raw["utc_timestamp"], format="%Y-%m-%d"
).values
opsd_names = {
    "solar": "solar",
    "onwind": "wind_onshore",
    "offwind": "wind_offshore",
}

#if carrier=='onwind': 
#    uncor=snakemake.input.onwind_uncor
#    cor=snakemake.input.onwind_cor
#elif carrier=='offwind':
#    uncor=snakemake.input.offwind_uncor
#    cor=snakemake.input.offwind_cor
#elif carrier=='solar':
#    uncor=snakemake.input.solar_uncor
#    cor=snakemake.input.solar_cor
uncor=snakemake.input.cf_uncor
cor=snakemake.input.cf_cor
cf_uncor = pd.read_csv(uncor, index_col=0, parse_dates=True)[country].drop_duplicates()
cf_cor = pd.read_csv(cor, index_col=0, parse_dates=True)[country].drop_duplicates()
carrier_opsd = opsd_names[carrier]
cf_opsd = pd.Series((opsd_raw[f"{country}_{carrier_opsd}_generation_actual"] / opsd_raw[f"{country}_{carrier_opsd}_capacity"]).values,index=opsd_date_range).dropna()
cf[carrier] = pd.DataFrame([cf_uncor, cf_cor, cf_opsd], index=['uncor', 'cor', 'opsd']).T.dropna()
cf[carrier] = cf[carrier][cf[carrier].index<='2018-12-31 23:00']

tech_colors = {"solar": "gold", "onwind": "royalblue", "offwind": "#6823bc",
                "wind": "#2393bc", "combined": "black", "range": "palegreen"}


fig,ax = plt.subplots(figsize=(10,10))
ax.scatter(cf[carrier].opsd, cf[carrier].uncor, color='grey', marker='x', label='uncorrected', zorder=3)
ax.scatter(cf[carrier].opsd, cf[carrier].cor, color=tech_colors[carrier], marker='+', label='corrected', alpha=0.8, zorder=3)
#ax.plot([0,1], [0,1], color='black', label='diagonal', zorder=3, linewidth=3)
slope_u, intercept_u, r_value, p_value, std_err = stats.linregress(cf[carrier].uncor, cf[carrier].opsd)
ax.plot([intercept_u, intercept_u+slope_u], [0,1], color='red', label='uncorrected', zorder=3, linewidth=3)  
slope_c, intercept_c, r_value, p_value, std_err = stats.linregress(cf[carrier].cor, cf[carrier].opsd)
ax.plot([intercept_c, intercept_c+slope_c], [0,1], color='black', label='corrected', zorder=3, linewidth=3)  
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_xlabel('Capacity factor [OPSD]')
#ax.set_title(opsd_names[carrier])
ax.set_yticks(np.arange(0.2,1.2,0.2))
#if carrier != 'onwind':
#  ax.set_yticklabels(['']*5)
#else:
ax.set_ylabel('Capacity factor [Atlite]')

legend_elements = [Line2D([0], [0], color='red', lw=4, label='uncorrected', marker='x', markeredgecolor='grey', markersize=20, markeredgewidth=4),
                 Line2D([0], [0], color='black', lw=4, label='corrected', marker='+', markeredgecolor=tech_colors[carrier], markersize=20, markeredgewidth=4),
                 #Line2D([0], [0], color='black', lw=4, label='diagonal')
                ]
ax.legend(handles=legend_elements, loc='lower right', fontsize=25)
ax.grid(zorder=1)
plt.savefig(snakemake.output.regression, bbox_inches='tight')
