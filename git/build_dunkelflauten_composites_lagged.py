import xarray as xr
import numpy as np
import pandas as pd
import os
import socket
from sys import exit

dir_out_variables=os.getenv('dir_out_variables')

dir_composites="%s/df_composites/"%(dir_out_variables)
if os.path.isdir('%s'%(dir_composites))==False: os.system('mkdir %s'%(dir_composites))

wr = os.getenv('regime').lower()
variable=os.getenv('variable')
running_days=int(os.getenv('running_days'))
clim_type=os.getenv('clim_type')

lag_inputs = [int(x) for x in os.getenv('lag_inputs').split(';')] #'-4;-2;0;2;4'
interval_input = int(os.getenv('interval_input')) #2
date_file = open("%s/wr/%s_startdates_wr_DE_0.06_48.txt"%(dir_out_variables, wr), "r")

#Read in start dates of Dunkelflauten
start_dates = date_file.read().splitlines()
start_dates = pd.to_datetime(start_dates, format='%Y%m%d_%H')

for lag_input in lag_inputs:
  lag = pd.Timedelta(lag_input, 'days')
  interval = pd.Timedelta(interval_input, 'days')
  lag_dates = [pd.date_range(x+lag, x+lag+interval, freq='3h') for x in start_dates]
  lag_dates = [subitem for item in [x.strftime('%Y%m%d_%H') for x in lag_dates] for subitem in item]
 
  #Absolute
  abs_files = ['%s/%s/%s/%s/%s_%s'%(dir_out_variables, date[0:4], date[4:6], variable, variable, date) for date in list(lag_dates) if os.path.isfile('%s/%s/%s/%s/%s_%s'%(dir_out_variables, date[0:4], date[4:6], variable, variable, date))]
  print('Absolute composite of %s files for %s Dunkelflauten and variable: %s'%(len(abs_files), wr, variable))
  abs_fields = xr.open_mfdataset(abs_files)
  abs_fields = abs_fields.where(abs_fields.apply(np.isfinite)).fillna(np.nan)
  abs_composite = abs_fields.mean(dim='time', skipna=True)
  abs_composite.to_netcdf("%s/%s_%sabs_lag%s_%s"%(dir_composites, wr, variable, lag_input, lag_input+interval_input))
  
#  #Anomaly
#  if variable=='SLP': break #do not compute anomalies for SLP as not necessary and also files do not exist!
#  ano_files = ['%s/anomaly/%s/%s/%sano_%s_%s'%(dir_out_variables, clim_type, variable, variable, date, clim_type) for date in list(lag_dates) if os.path.isfile('%s/anomaly/%s/%s/%sano_%s_%s'%(dir_out_variables, clim_type, variable, variable, date, clim_type))]
#  print('Anomaly composite of %s files for %s Dunkelflauten and variable: %s'%(len(ano_files), wr, variable))
#  ano_fields = xr.open_mfdataset(ano_files)
#  ano_fields = ano_fields.where(ano_fields.apply(np.isfinite)).fillna(np.nan)
#  ano_composite = ano_fields.mean(dim='time', skipna=True)
#  ano_composite.to_netcdf("%s/%s_%sano_%s_clim_lag%s_%s"%(dir_composites, wr, variable, clim_type, lag_input, lag_input+interval_input))
#
