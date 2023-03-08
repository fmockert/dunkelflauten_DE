import numpy as np
import pandas as pd
from datetime import timedelta

# Random pick function
#create function with input data and select_variable, double=True allows LC to appear multiple times
def random_pick(input_data, select_variable, double):
    LC_random = dict.fromkeys(lifecycles, {})
    for regime in regime_names_sorted[0:7]:
        LC_random[regime] = pd.DataFrame(columns=["onset", "decay"])
        LC = input_data[regime]#["LC"]["onset"]
        DF = df_intervals[df_intervals["regime"]==regime]
        for i in DF.index:
            LC_choosefrom = pd.DataFrame(columns=["onset", "decay"])
            for year in np.arange(year_start, year_end+1):
                date = pd.Timestamp(str(year)+"-"+str(DF["start"][i].month)+"-"+str(DF["start"][i].day)+" "+str(DF["start"][i].hour)+":00:00")
                start = date-pd.Timedelta(days=45)
                end = date+pd.Timedelta(days=45)
                LC_choosefrom = LC_choosefrom.append(LC[(LC[select_variable]>=start) & (LC[select_variable]<=end)])
            sample = LC_choosefrom.sample()
            if double==False:
                while sample.index[0] in LC_random[regime].index:
                    sample = LC_choosefrom.sample()
            LC_random[regime] = LC_random[regime].append(sample)
        LC_random[regime] = LC_random[regime].reset_index(drop=True)
    output_data = LC_random
    return output_data
# Lifecycle import
def lifecycle_import(number):
    filename = snakemake.input.lifecycle_name
    lifecycle = pd.read_csv(snakemake.input.lifecycle+"/"+filename+str(number)+".txt", skiprows=10, header=None, delim_whitespace=True)[[1,5]]#.squeeze()
    lifecycle.columns= ["onset", "decay"]
    lifecycle[["onset", "decay"]] = lifecycle[["onset", "decay"]].apply(pd.to_datetime, format='%Y%m%d_%H')
    return lifecycle
#{"No":0, "AT":1, "AR":2, "GL":3, "EuBL":4, "ScBL":5, "ZO":6, "ScTr":7}
lifecycles = {"AT": lifecycle_import(1),
      "GL": lifecycle_import(3),
      "AR": lifecycle_import(2),
      "ScTr": lifecycle_import(7),
      "EuBL": lifecycle_import(4),
      "ScBL": lifecycle_import(5),
      "ZO": lifecycle_import(6)}

year_start, year_end = snakemake.config["years"]      
# Regime names
regime_names_sorted = ["AT", "ZO", "ScTr", "AR", "EuBL", "ScBL", "GL", "No"]

# Data input
df_intervals = pd.read_csv(snakemake.input.df_intervals, index_col=0, parse_dates=True)
df_intervals[["start", "end", "start_wr", "end_wr"]] = df_intervals[["start", "end", "start_wr", "end_wr"]].apply(pd.to_datetime)

df_lc_intervals = pd.read_csv(snakemake.input.df_lc_intervals, index_col=0, parse_dates=True)
for column in ["DF_start", "DF_end", "AT_LC_onset", "AT_LC_decay", "ZO_LC_onset", "ZO_LC_decay", 
               "ScTr_LC_onset", "ScTr_LC_decay", "AR_LC_onset", "AR_LC_decay", "EuBL_LC_onset", 
               "EuBL_LC_decay", "ScBL_LC_onset", "ScBL_LC_decay", "GL_LC_onset", "GL_LC_decay"]:
    df_lc_intervals[column] = df_lc_intervals[column].apply(pd.to_datetime)
    
df_regimes = dict.fromkeys(lifecycles, {})
for regime in regime_names_sorted[0:7]:
    df_regimes[regime] = df_lc_intervals[df_lc_intervals["DF_regime"]==regime][[regime+"_LC_onset", regime+"_LC_decay"]]
    df_regimes[regime] = df_regimes[regime].rename(columns={regime+"_LC_onset": "onset", regime+"_LC_decay": "decay"}).reset_index(drop=True)

# Bootstraping    
#True = Doubles allowed, random pick with putting back
#Attention: the for loops are taking long to compute with 1000cycles.
bootstrapLC = dict.fromkeys(lifecycles,pd.Series())
for i in np.arange(0,1000):
        random_regimes = random_pick(lifecycles, "onset", True)
        for regime in regime_names_sorted[0:7]:
            randomLC = (random_regimes[regime].decay-random_regimes[regime].onset)/pd.Timedelta(hours=1)/24
            bootstrapLC[regime] = bootstrapLC[regime].append(pd.Series(np.mean(randomLC)))
        bootstrapLC_frame = pd.DataFrame.from_dict(bootstrapLC)
        bootstrapLC_frame = bootstrapLC_frame.reset_index(drop=True)
        bootstrapLC_frame.to_csv(snakemake.output.bootstrap) #saving file after every step and overwriting previous
        
# Subsampling
#True = Doubles allowed, random pick with putting back
subsamplingLC = dict.fromkeys(df_regimes,pd.Series())
for i in np.arange(0,1000):
        random_df_regimes = random_pick(df_regimes, "onset", True)
        for regime in regime_names_sorted[0:7]:
            random_df_LC = (random_df_regimes[regime].decay-random_df_regimes[regime].onset)/pd.Timedelta(hours=1)/24
            subsamplingLC[regime] = subsamplingLC[regime].append(pd.Series(np.mean(random_df_LC)))
        subsamplingLC_frame = pd.DataFrame.from_dict(subsamplingLC)
        subsamplingLC_frame = subsamplingLC_frame.reset_index(drop=True)
        subsamplingLC_frame.to_csv(snakemake.output.subsampling) #saving file after every step and overwriting previous
