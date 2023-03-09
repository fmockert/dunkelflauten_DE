"""Build life cycle intervals that include a Dunkelflaute"""

__author__ = "Fabian Mockert"
__copyright__ = "Copyright 2023, Fabian Mockert"

import numpy as np
import pandas as pd
from datetime import timedelta


# Weather regime lifecycles separated
def lifecycle_import(number):
    #filename = "Z0500_N81_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local_cluster"
    #lifecycles = pd.read_csv(snakemake.input.lifecycles+str(number)+".txt", skiprows=10, header=None, delim_whitespace=True)[[1,5]]#.squeeze()
    lifecycles = pd.read_csv(number, skiprows=10, header=None, delim_whitespace=True)[[1,5]]#.squeeze()
    #lifecycles = pd.read_csv("../../data/bundle/regimes/Z0500_N81_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local_cluster"+str(number)+".txt", skiprows=10, header=None, delim_whitespace=True)[[1,5]]#.squeeze()
    lifecycles.columns= ["onset", "decay"]
    lifecycles[["onset", "decay"]] = lifecycles[["onset", "decay"]].apply(pd.to_datetime, format='%Y%m%d_%H')
    return lifecycles

if __name__ == "__main__":  

  converter = {0:"No", 1:"AT", 2:"AR", 3:"GL", 4:"EuBL", 5:"ScBL", 6:"ZO", 7:"ScTr"}
  era5_reconvert = {"No":0, "AT":1, "AR":2, "GL":3, "EuBL":4, "ScBL":5, "ZO":6, "ScTr":7}
  regime_names_sorted = ["AT", "ZO", "ScTr", "AR", "EuBL", "ScBL", "GL", "No"]
  
  # Import dunkelflauten
  df_intervals = pd.read_csv(snakemake.input.df_intervals, index_col=0, parse_dates=[1,2,4,5])
  
  regimes = {"AT": {"LC": lifecycle_import(snakemake.input.LC1), "color": '#6100B3'},
        "GL": {"LC": lifecycle_import(snakemake.input.LC3), "color": '#0000FE'},
        "AR": {"LC": lifecycle_import(snakemake.input.LC2), "color": '#FECF0A'},
        "ScTr": {"LC": lifecycle_import(snakemake.input.LC7), "color": '#FB6207'},
        "EuBL": {"LC": lifecycle_import(snakemake.input.LC4), "color": '#117B00'},
        "ScBL": {"LC": lifecycle_import(snakemake.input.LC5), "color": '#0B5300'},
        "ZO": {"LC": lifecycle_import(snakemake.input.LC6), "color": '#FB0005'}
            }
  
  # Create df_lc_intervals
  columns = ["DF_start", "DF_end", "DF_regime", 
             "AT_LC_num", "AT_LC_onset", "AT_LC_decay", 
             "ZO_LC_num", "ZO_LC_onset", "ZO_LC_decay", 
             "ScTr_LC_num", "ScTr_LC_onset", "ScTr_LC_decay", 
             "AR_LC_num", "AR_LC_onset", "AR_LC_decay", 
             "EuBL_LC_num", "EuBL_LC_onset", "EuBL_LC_decay", 
             "ScBL_LC_num", "ScBL_LC_onset", "ScBL_LC_decay", 
             "GL_LC_num", "GL_LC_onset", "GL_LC_decay"]
  df_lc_intervals = pd.DataFrame(columns=columns)
  df_lc_intervals[["DF_start", "DF_end", "DF_regime"]]= df_intervals[["start","end","regime"]]
  
  for df_num in np.arange(0, len(df_lc_intervals)):
      df_start = df_lc_intervals["DF_start"][df_num]
      df_end = df_lc_intervals["DF_end"][df_num]
      
      for regime_name in regime_names_sorted[:7]:
          lifecycles = regimes[regime_name]["LC"]
          lc_num_current=0
          for lc_num in np.arange(0, len(lifecycles)):
              if lc_num>=lc_num_current-1:
                  # increasing the computation speed by not going through obviously past lifecycles
                  lc_onset = lifecycles["onset"][lc_num]
                  lc_decay = lifecycles["decay"][lc_num]
                  lc_num_current = lc_num
                  if (((df_start>=lc_onset) & (df_start<=lc_decay)) | ((df_end>=lc_onset) & (df_end<=lc_decay))):
                      df_lc_intervals.at[df_num, regime_name+"_LC_num"] = lc_num
                      df_lc_intervals.at[df_num, regime_name+"_LC_onset"] = lc_onset
                      df_lc_intervals.at[df_num, regime_name+"_LC_decay"] = lc_decay
                      break
  df_lc_intervals.to_csv(snakemake.output.df_lc_intervals)
  