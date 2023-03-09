"""Linking Dunkelflauten with weather regimes."""

__author__ = "Fabian Mockert"
__copyright__ = "Copyright 2023, Fabian Mockert"

import numpy as np
import pandas as pd
from datetime import timedelta

# Weather regime
def load_startdates(filename):
    dataframe = pd.read_csv(filename, header=None)
    dataframe = dataframe.apply(lambda x: pd.to_datetime(x, format='%Y%m%d_%H'))
    return dataframe

if __name__ == "__main__":    
  wr = pd.read_csv(snakemake.input.weatherregime, index_col =[0], parse_dates=[0]).squeeze()
  
  regimes = {"No": {"number":0, "color":"#6B6B6B", "startdates": load_startdates(snakemake.input.No)},
             "AT": {"number":1, "color":"#6100B3", "startdates": load_startdates(snakemake.input.AT)},
             "GL": {"number":3, "color":"#0000FE", "startdates": load_startdates(snakemake.input.GL)},
             "AR": {"number":2, "color":"#FECF0A", "startdates": load_startdates(snakemake.input.AR)},
             "ScTr": {"number":7, "color":"#FB6207", "startdates": load_startdates(snakemake.input.ScTr)},
             "EuBL": {"number":4, "color":"#117B00", "startdates": load_startdates(snakemake.input.EuBL)},
             "ScBL": {"number":5, "color":"#0B5300", "startdates": load_startdates(snakemake.input.ScBL)},
             "ZO": {"number":6, "color":"#FB0005", "startdates": load_startdates(snakemake.input.ZO)}
            }
            
  regime_transform = {0: "No", 1: "AT", 3: "GL", 2: "AR", 7: "ScTr", 4: "EuBL", 5: "ScBL", 6: "ZO"}
  
  # Creating the dominant weather regime intervals
  wr_intervals = pd.DataFrame(columns=["start_date", "end_date", "regime", "regime_number"])
  i=0
  while i in range(0,len(wr)-1):
      start_date = wr.index[i]
      n=i
      while ((wr[n]==wr[i]) & (n<len(wr)-1)):
          end_date = wr.index[n]
          n += 1
      regime_name = wr[i]
      wr_intervals = wr_intervals.append(pd.DataFrame([[start_date, end_date, regime_name, regimes[regime_name]['number']]], columns= ["start_date", "end_date", "regime", "regime_number"]), ignore_index=True)
      i = n
  wr_intervals.to_csv(snakemake.output.wr_intervals)
  
  # Creating df_intervals including regimes
  df_timestamps = pd.read_csv(snakemake.input.df_timestamps, header=None, index_col=0)
  df_intervals = pd.read_csv(snakemake.input.df_intervals, index_col=0, parse_dates=True)
  df_intervals[["start", "end", "start_wr", "end_wr"]] = df_intervals[["start", "end", "start_wr", "end_wr"]].apply(pd.to_datetime)
  
  df_intervals["regime"] = ""
  for regime in ["AT", "ZO", "ScTr", "AR", "EuBL", "ScBL", "GL", "No"]:
      for i in range(0,len(regimes[regime]["startdates"].squeeze())):
          df_intervals.loc[df_intervals.start_wr==regimes[regime]["startdates"].squeeze()[i], "regime"] = regime
          
  df_intervals.to_csv(snakemake.output.df_intervals)
  
  
  
  