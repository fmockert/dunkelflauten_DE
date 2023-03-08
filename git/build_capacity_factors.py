"""Calculates wind and solar capacity factors per country."""

__author__      = "Fabian Mockert, Fabian Neumann"
__copyright__   = "Copyright 2019-2020, Fabian Mockert, Fabian Neumann"

import atlite
import geopandas as gpd
from atlite.gis import ExclusionContainer
import matplotlib.pyplot as plt
from shapely.geometry import mapping
import xarray as xr
import numpy as np
import pandas as pd

import logging
logging.basicConfig(level="INFO")

if __name__ == "__main__":

  config = snakemake.config["technology"]
  weights_year = {t:{} for t in config.keys()}
  weights_clipped = {t:{} for t in config.keys()}
  weighted = {t:pd.Series([], name='DE', dtype="float") for t in config.keys()}

  year_start, year_end = snakemake.config["years"]
  
  #reshaping on and off shore with excluder
  excluder = ExclusionContainer(crs=3035)

  for year in np.arange(year_start, year_end+1): #1979,2018+1):
    # Load cutout
    cutout = atlite.Cutout(snakemake.input.cutouts+"/cutouts_germany_%s.nc"%(year))
    # Load region shape
    if carrier=='offwind': countries=gpd.read_file(snakemake.input.offshore_shapes).set_index('name')
    else: countries=gpd.read_file(snakemake.input.country_shapes).set_index('name')
    shape = countries.to_crs(excluder.crs).loc[["DE"]].geometry
    
    matrix = cutout.indicatormatrix([countries.squeeze()["DE"]])
    cutout.prepare()
    # Create weights
    if 'wind' in carrier:
      weights = cutout.wind(turbine=config[carrier]['turbine'], per_unit=False, capacity_factor=True)
      weights_year[carrier][year] = cutout.wind(matrix=matrix, turbine=config[carrier]['turbine'], layout=weights, per_unit=True, capacity_factor=False)
    elif 'solar' in carrier:
      weights = cutout.pv(panel=config[carrier]['panel'], orientation=config[carrier]['orientation'], per_unit=False, capacity_factor=True)
      weights_year[carrier][year] = cutout.pv(matrix=matrix, panel=config[carrier]['panel'], orientation=config[carrier]['orientation'], layout=weights, per_unit=True, capacity_factor=False)
    else: raise NotImplementedError("Can not calculate capacity factors for carrier %s. Only 'solar', 'onwind', 'offwind' are implemented."%(carrier))
    # Masking weights
    weights_clipped[carrier][year] = weights.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False).rio.write_crs("epsg:4326", inplace=False).rio.clip(shape.geometry.apply(mapping), shape.crs, drop=False)
    
    weighted[carrier] = weighted[carrier].append(pd.Series(weights_year[carrier][year].values.flatten(), index=weights_year[carrier][year].time.values))
    weighted[carrier].name='DE'
  weighted[carrier].to_csv(snakemake.ouput.cf, header=True)

  # Save clipped weights
  weights_clipped_concat = {}
  weights_clipped_concat[carrier] = xr.concat([weights_clipped[carrier][i] for i in np.arange(year_start,year_end+1)], dim='year') #change year indexing from 0-44 to actual years
  weights_clipped_concat[carrier] = weights_clipped_concat[carrier].assign_coords(year=('year', np.arange(year_start,year_end+1)))
  weights_clipped_concat[carrier].to_netcdf(snakemake.output.cf_weights)