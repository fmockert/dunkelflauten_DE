"""Compare capacity factors with OPSD data and make linear regression
OPSD data from https://data.open-power-system-data.org/time_series/2019-06-05

Includes a fallback in case OPSD does not include the chosen country.
"""

__author__ = "Fabian Mockert, Fabian Neumann"
__copyright__ = "Copyright 2023, Fabian Mockert, Fabian Neumann"

import pandas as pd
import numpy as np
from datetime import timedelta
from scipy import stats

import logging

logging.basicConfig(level="INFO")

if __name__ == "__main__":

  country = "DE"#snakemake.wildcards.country
  source = "atlite"#snakemake.wildcards.source
  carrier = snakemake.wildcards.carrier
  
  cf = pd.read_csv(snakemake.input.cf, index_col=0, parse_dates=True)[country]
  
  opsd_raw = pd.read_csv(snakemake.input.opsd)
  opsd_date_range = pd.to_datetime(
      opsd_raw["utc_timestamp"], format="%Y-%m-%d"
  ).values
  
  opsd_names = {
      "solar": "solar",
      "onwind": "wind_onshore",
      "offwind": "wind_offshore",
  }
  carrier_opsd = opsd_names[carrier]
  
  cf_mod = pd.DataFrame()
  
  try:
      cf_mod[carrier] = (
          opsd_raw[f"{country}_{carrier_opsd}_generation_actual"]
          / opsd_raw[f"{country}_{carrier_opsd}_capacity"]
      )
      cf_mod.index = opsd_date_range
      cf_mod = cf_mod[cf_mod.index.isin(cf.index)]
      cf_mod = cf_mod.dropna()
      cf_mod[f"{carrier}_{source}"] = cf[cf.index.isin(cf_mod.index)]
  
      slope, intercept, r_value, p_value, std_err = stats.linregress(
          cf_mod[f"{carrier}_{source}"], cf_mod[f"{carrier}"]
      )
  except:
      slope = snakemake.config['technology'][carrier]["fallback"]["slope"]
      intercept = snakemake.config['technology'][carrier]["fallback"]["intercept"]
  
  cf = cf * slope + intercept
  
  cf.to_csv(snakemake.output.cfcorr)
