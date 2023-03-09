"""Plot the capacity factor weighting across Germany"""

__author__ = "Fabian Mockert"
__copyright__ = "Copyright 2023, Fabian Mockert"

import atlite
import geopandas as gpd
from atlite.gis import ExclusionContainer
import matplotlib.pyplot as plt
from shapely.geometry import mapping


# Size configuration
size = 30; params = {'legend.fontsize': size,'figure.figsize': (10, 10), 'axes.labelsize': size, 'axes.titlesize': size, 'xtick.labelsize': size, 'ytick.labelsize': size}
import matplotlib.pylab as pylab; pylab.rcParams.update(params)

import xarray as xr
import numpy as np
import pandas as pd
import os

if __name__ == "__main__":
  country = "DE"#snakemake.wildcards.country
  source = "atlite"#snakemake.wildcards.source
  carrier = snakemake.wildcards.carrier
  start_year, end_year = snakemake.config["years"]
  
  technology = snakemake.config["technology"]
  weights_year = {t:{} for t in technology.keys()}
  weights_clipped = {t:{} for t in technology.keys()}
  weighted = {t:pd.Series([], name='DE', dtype="float") for t in technology.keys()}
  
  # reshaping on and off shore with excluder
  excluder = ExclusionContainer(crs=3035)
  weights_clipped_concat = {}
  weights_clipped_m = {}
  for carrier, weights_data in zip(['onwind', 'offwind', 'solar'],[snakemake.input.weights_on, snakemake.input.weights_off, snakemake.input.weights_pv]):
    weights_clipped_concat[carrier] = xr.open_dataset(weights_data)
    weights_clipped_m[carrier] = weights_clipped_concat[carrier].loc[dict(year=np.arange(start_year,end_year+1))].mean(dim='year')
  
  # Plot weights clipped merged to mean
  countries_land=gpd.read_file(snakemake.input.countries_land).set_index('name')
  countries_ocean=gpd.read_file(snakemake.input.countries_ocean).set_index('name')
  excluder = ExclusionContainer(crs=3035)
  shape_land = countries_land.to_crs(excluder.crs).loc[["DE"]].geometry
  shape_ocean = countries_ocean.to_crs(excluder.crs).loc[["DE"]].geometry
  
  # Wind
  fig, ax = plt.subplots() 
  wmin = min(weights_clipped_m['onwind'].min(), weights_clipped_m['offwind'].min())
  wmax = max(weights_clipped_m['onwind'].max(), weights_clipped_m['offwind'].max())
  
  on = ax.pcolormesh(weights_clipped_m['onwind'].x,weights_clipped_m['onwind'].y,weights_clipped_m['onwind'], vmin=wmin, vmax=wmax, cmap='jet')
  off = ax.pcolormesh(weights_clipped_m['offwind'].x,weights_clipped_m['offwind'].y,weights_clipped_m['offwind'], vmin=wmin, vmax=wmax, cmap='jet')
  shape_land.to_crs(4326).plot(ax=ax, edgecolor='k', color='none')
  shape_ocean.to_crs(4326).plot(ax=ax, edgecolor='k', color='none')
  cbar = fig.colorbar(on).set_label(label='Capacity factor weighting', size=size)
  ax.set_xlabel('Longitude [°E]', fontsize=size)
  ax.set_ylabel('Latitude [°N]', fontsize=size)
  ax.set_xlim(3,15.2)
  ax.set_ylim(47.2,56)
  ax.grid(color='lightgrey', linestyle='dotted')
  plt.savefig(snakemake.output.wind_weighting, bbox_inches='tight')
  
  # Solar
  fig, ax = plt.subplots() 
  smin = weights_clipped_m['solar'].min()
  smax = weights_clipped_m['solar'].max()
  pv = ax.pcolormesh(weights_clipped_m['solar'].x,weights_clipped_m['solar'].y,weights_clipped_m['solar'], vmin=smin, vmax=smax, cmap='jet')
  
  shape_land.to_crs(4326).plot(ax=ax, edgecolor='k', color='none')
  shape_ocean.to_crs(4326).plot(ax=ax, edgecolor='k', color='none')
  
  cbar = fig.colorbar(pv).set_label(label='Capacity factor weighting', size=size)
  
  ax.set_xlabel('Longitude [°E]', fontsize=size)
  ax.set_ylabel('Latitude [°N]', fontsize=size)
  ax.set_xlim(3,15.2)
  ax.set_ylim(47.2,56)
  ax.grid(color='lightgrey', linestyle='dotted')
  ax.set_yticklabels(['']*9)
  plt.savefig(snakemake.output.solar_weighting, bbox_inches='tight')
