"""Calculates Dunkelflauten periods per country for different thresholds/durations."""

__author__ = "Fabian Mockert, Fabian Neumann"
__copyright__ = "Copyright 2023, Fabian Mockert, Fabian Neumann"

import pandas as pd
import numpy as np
import math
from datetime import timedelta

import logging

logging.basicConfig(level="INFO")


def map_to_wr_timestamps(ts):
    if ts.hour in range(0, 3):
        ts = ts.replace(hour=0)
    elif ts.hour in range(3, 9):
        ts = ts.replace(hour=6)
    elif ts.hour in range(9, 15):
        ts = ts.replace(hour=12)
    elif ts.hour in range(15, 21):
        ts = ts.replace(hour=18)
    elif ts.hour in range(21, 24):
        ts = ts.replace(hour=0) + timedelta(days=1)

    return ts

def dunkelflaute_moving_average(combined_cf, level, window):
    """
    Smoothing the raw capacity factors `combined_cf` with the
    rolling function over a time period of `window`. `window` 
    symbolises the minimum period a DF needs to be. As soon as
    `combined` (smoothened) has a value below `level`, there is
    a DF.
    """
    combined = combined_cf.rolling(window, min_periods=0, center=True).mean()
    
    dunkelflaute = pd.DataFrame(columns=["start", "end", "hours"])
    t = 0 
    while t < len(combined)-int((window/2-1)):
        if combined[t] < level:
            T_count = 1
            while (t + T_count < len(combined)) & (
                combined[t + T_count - 1] < level
            ):
                T_count += 1
            
            ts_start = combined.index[t-int((window/2))]
            ts_end = combined.index[t + T_count - 1 + int((window/2-1))]
            hours = [T_count - 1 + int((window-1))]
            dunkelflaute = dunkelflaute.append(
                pd.DataFrame(
                    dict(zip(dunkelflaute.columns, (ts_start, ts_end, hours)))
                ),
                ignore_index=True,
            )
            
            t += T_count+48
        else:
            t += 1

    return dunkelflaute

def dunkelflaute_join_no_distance(dunkelflaute):
    """
    When applying the map_to_wr_timestamps function, there is the 
    possibility that 2 consecutive dunkelflauten have no distance 
    between each other, those dunkelflauten will then be joined.
    """
    dunkelflaute["start_wr"] = dunkelflaute.start.apply(map_to_wr_timestamps)
    dunkelflaute["end_wr"] = dunkelflaute.end.apply(map_to_wr_timestamps)
    dunkelflaute["hours_wr"] = (dunkelflaute.end_wr - dunkelflaute.start_wr).astype(
        "timedelta64[h]"
    )

    i = 0
    while i < len(dunkelflaute) - 1:
        if (
            dunkelflaute.at[i + 1, "start_wr"] - dunkelflaute.at[i, "end_wr"]
        ) <= pd.Timedelta(0, unit="h"):
            dunkelflaute.at[i, "end_wr"] = dunkelflaute.at[i + 1, "end_wr"]
            dunkelflaute.at[i, "end"] = dunkelflaute.at[i + 1, "end"]
            dunkelflaute.at[i, "hours_wr"] = (
                dunkelflaute.at[i + 1, "end_wr"] - dunkelflaute.at[i, "start_wr"]
            ) / np.timedelta64(1, "h")
            dunkelflaute.at[i, "hours"] = (
                dunkelflaute.at[i + 1, "end"] - dunkelflaute.at[i, "start"]
            ) / np.timedelta64(1, "h")
            dunkelflaute = dunkelflaute.drop(i + 1).reset_index(drop=True)
        else:
            i += 1
    return dunkelflaute


if __name__ == "__main__":

  config = snakemake.config
  level = 0.06#float(snakemake.wildcards.level)
  duration = 48#int(snakemake.wildcards.duration)
  country = "DE" #snakemake.wildcards.country
  method = "moving_average"#snakemake.wildcards.method
  
  logging.info(
      f"Calculate Dunkelflauten for country {country} with definition: level={level}, duration={duration}"
  )
  
  combined_cf = pd.read_csv(
      snakemake.input.capacity_factors, index_col=0, parse_dates=True
  ).squeeze()
  
  dunkelflaute = dunkelflaute_moving_average(combined_cf, level, duration)
  dunkelflaute = dunkelflaute_join_no_distance(dunkelflaute)
  
  timestamps = pd.Series(dtype="timedelta64[ns]")
  for i, flaute in dunkelflaute.iterrows():
      to_append = pd.Series(
          pd.date_range(flaute["start_wr"], flaute["end_wr"], freq="6H")
      )
      timestamps = timestamps.append(to_append, ignore_index=True)
  
  dateformat = "%Y%m%d_%H"
  timestamps = timestamps.apply(lambda t: t.strftime(dateformat))
  timestamps.to_csv(snakemake.output.timestamps, index=False, header=False)
  dunkelflaute.to_csv(snakemake.output.intervals)
