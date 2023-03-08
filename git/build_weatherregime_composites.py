import xarray as xr
import numpy as np
import os
import socket
import pandas as pd
from sys import exit

#if socket.gethostname() == "imk-tnn-famo": !Get from main_composite_workflow.sh
#    #from worklaptop
#    dir_out_variables="/home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/"
#    dir_dunkelflauten="/home/qb3829/Documents/HiWi_IAI_20/dig/outputs/wr/" 
#else:
#    #from uc2
#    dir_out_variables="/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/dunkelflauten_variables/"
#
dir_out_variables=os.getenv('dir_out_variables')

dir_composites="%s/wr_composites/"%(dir_out_variables)
if os.path.isdir('%s'%(dir_composites))==False: os.system('mkdir %s'%(dir_composites))

#Load dunkelflauten timestamps
#   SELECTION OF WEATHER REGIME BY BASH OR SBATCH TO RUN IN PARALLEL
wr = os.getenv('regime').lower()
season = 'NDJFM'
#file path local to local file, not to file on UC2

date_file = open("%s/wr/%s_regime_times_79_18_%s.txt"%(dir_out_variables, wr, season), "r")

times = date_file.read().splitlines()

#Generate file names of dunkelflauten times and variable, for anomaly and climatology
variable=os.getenv('variable')
running_days=int(os.getenv('running_days'))
clim_type=os.getenv('clim_type')

#Absolute
abs_files = ['%s/%s/%s/%s/%s_%s'%(dir_out_variables, date[0:4], date[4:6], variable, variable, date) for date in list(times) if os.path.isfile('%s/%s/%s/%s/%s_%s'%(dir_out_variables, date[0:4], date[4:6], variable, variable, date))]
print('Absolute composite of %s files for %s Dunkelflauten and variable: %s'%(len(abs_files), wr, variable))
abs_fields = xr.open_mfdataset(abs_files)
abs_fields = abs_fields.where(abs_fields.apply(np.isfinite)).fillna(np.nan)
abs_composite = abs_fields.mean(dim='time', skipna=True)
abs_composite.to_netcdf("%s/%s_%sabs_%s"%(dir_composites, wr, variable, season))

#Anomaly
if variable=='SLP': exit() #do not compute anomalies for SLP as not necessary and also files do not exist!
ano_files = ['%s/anomaly/%s/%s/%sano_%s_%s'%(dir_out_variables, clim_type, variable, variable, date, clim_type) for date in list(times) if os.path.isfile('%s/anomaly/%s/%s/%sano_%s_%s'%(dir_out_variables, clim_type, variable, variable, date, clim_type))]
print('Anomaly composite of %s files for %s Dunkelflauten and variable: %s'%(len(ano_files), wr, variable))
ano_fields = xr.open_mfdataset(ano_files)
ano_fields = ano_fields.where(ano_fields.apply(np.isfinite)).fillna(np.nan)
ano_composite = ano_fields.mean(dim='time', skipna=True)
ano_composite.to_netcdf("%s/%s_%sano_%s_clim_%s"%(dir_composites, wr, variable, clim_type, season))
