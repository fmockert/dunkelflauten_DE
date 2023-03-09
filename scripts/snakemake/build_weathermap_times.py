"""Build the specific times to plot composites of atmospheric field variables for Dunkelflauten. 
Those scripts are not involved in the snakemake process.
"""

__author__ = "Fabian Mockert, Fabian Neumann"
__copyright__ = "Copyright 2023, Fabian Mockert, Fabian Neumann"

import datetime

import pandas as pd
import numpy as np


if __name__ == "__main__":
  wr = pd.read_csv(snakemake.input.wr, index_col=[0], parse_dates=[0]).squeeze()
  start_year, end_year = snakemake.config['years']
  try:
      df_intervals = pd.read_csv(snakemake.input.df_intervals, index_col=0, parse_dates=[1,2,4,5])
      initial_dataframe = pd.DataFrame(columns=['start_wr', 'end_wr'])
      df_sorted = {0: {'df_start_end': initial_dataframe, 'wr': 'No'},
  		    1: {'df_start_end': initial_dataframe, 'wr': 'AT'},
  		    3: {'df_start_end': initial_dataframe, 'wr': 'GL'},
  		    2: {'df_start_end': initial_dataframe, 'wr': 'AR'},
  	    7: {'df_start_end': initial_dataframe, 'wr': 'ScTr'},
  	    4: {'df_start_end': initial_dataframe, 'wr': 'EuBL'},
  	    5: {'df_start_end': initial_dataframe, 'wr': 'ScBL'},
  	    6: {'df_start_end': initial_dataframe, 'wr':'ZO'}
  	    }
      df_sorted = {'No': initial_dataframe.copy(), 'AT': initial_dataframe.copy(), 'GL': initial_dataframe.copy(), 'AR': initial_dataframe.copy(), 'ScTr': initial_dataframe.copy(),
  	    'EuBL': initial_dataframe.copy(), 'ScBL': initial_dataframe.copy(), 'ZO': initial_dataframe.copy()}  
      for i in df_intervals.index:
          start = df_intervals.start_wr[i]
          end = df_intervals.end_wr[i]
          df_actual = wr[(wr.index>=start) & (wr.index<=end)]
          wr_max = df_actual.value_counts().idxmax()
          df_sorted[wr_max] = df_sorted[wr_max].append(dict(zip(['start_wr', 'end_wr'], [start, end])), ignore_index=True)
      timestamps_all = pd.Series()
      for regime in df_sorted.keys():
          timestamps = pd.Series()
          for j in df_sorted[regime].index:
              timestamps = timestamps.append((pd.Series(pd.date_range(df_sorted[regime].start_wr[j], df_sorted[regime].end_wr[j], freq='3h'))).apply(lambda x: x.strftime('%Y%m%d_%H')))
          timestamps_all = timestamps_all.append(timestamps)
          timestamps.to_csv(snakemake.output["times"+str(regime)], index=False, header=False)  
      timestamps_all.to_csv(snakemake.output["timesall"], index=False, header=False)  
      for regime in df_sorted.keys():
          start_wr = df_sorted[regime].start_wr.apply(lambda x: x.strftime('%Y%m%d_%H'))
          start_wr.to_csv(snakemake.output[regime], index=False, header=False)
      start_all = df_intervals.start_wr.apply(lambda x: x.strftime('%Y%m%d_%H'))
      start_all.to_csv(snakemake.output['regime_all'], index=False, header=False)  
  except pd.errors.EmptyDataError:  
      # create dummy outputs
      open(snakemake.output.No, 'a').close()
      open(snakemake.output.AT, 'a').close()
      open(snakemake.output.GL, 'a').close()
      open(snakemake.output.AR, 'a').close()
      open(snakemake.output.ScTr, 'a').close()
      open(snakemake.output.EuBL, 'a').close()
      open(snakemake.output.ScBL, 'a').close()
      open(snakemake.output.ZO, 'a').close()  
