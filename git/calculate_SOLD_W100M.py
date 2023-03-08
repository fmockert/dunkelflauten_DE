import xarray as xr
import pandas as pd
import os
import socket
import numpy as np
import xarray.ufuncs as xu #deprecated, use np.ufuncs directly
#rsync /home/qb3829/Documents/HiWi_IAI_20/dig/scripts/data_and_analysis/calculate_SOLD_W100M.py /home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/
dir_era5=os.getenv('dir_era5')
dir_out_variables=os.getenv('dir_out_variables')

idx = int(os.environ["SLURM_ARRAY_TASK_ID"])

dates = pd.date_range('1979-01-01', '1979-01-02', freq='D').strftime('%Y%m%d')
#dates = np.array_split(dates, 40)[idx]
for date in dates:
  year=date[0:4]; month=date[4:6]; day=date[6:8]
  if os.path.isdir('%s/%s/%s/SOLD'%(dir_out_variables, year, month))==False: os.system(' mkdir %s/%s/%s/SOLD'%(dir_out_variables, year, month))
  if os.path.isdir('%s/%s/%s/W100M'%(dir_out_variables, year, month))==False: os.system(' mkdir %s/%s/%s/W100M'%(dir_out_variables, year, month))
  #SOLD
  SSR_tmp = xr.open_mfdataset("%s/%s/%s/SSR/SSR_%s%s%s_*"%(dir_out_variables, year, month, year, month, day))['SSR'].sum(dim='time')
  SSRC_tmp = xr.open_mfdataset("%s/%s/%s/SSRC/SSRC_%s%s%s_*"%(dir_out_variables, year, month, year, month, day))['SSRC'].sum(dim='time')
  SOLD_tmp = (SSR_tmp/SSRC_tmp).rename('SOLD') #.rename(name_dict={'__xarray_dataarray_variable__':'SOLD'})
  for hour in ["00", "03", "06", "09", "12", "15", "18", "21"]:
    SOLD = SOLD_tmp.expand_dims(time=[pd.Timestamp('%s-%s-%sT%s:00:00'%(year, month, day, hour))])
    SOLD.to_netcdf("%s/%s/%s/SOLD/SOLD_%s%s%s_%s"%(dir_out_variables, year, month, year, month, day, hour))
    #W100M
    U_tmp = xr.open_mfdataset("%s/%s/%s/U100/U100_%s%s%s_%s"%(dir_out_variables, year, month, year, month, day, hour))['U100']
    V_tmp = xr.open_mfdataset("%s/%s/%s/V100/V100_%s%s%s_%s"%(dir_out_variables, year, month, year, month, day, hour))['V100']
    W_tmp = xu.sqrt(U_tmp**2+V_tmp**2).rename('W100M')
    W_tmp.to_netcdf("%s/%s/%s/W100M/W100M_%s%s%s_%s"%(dir_out_variables, year, month, year, month, day, hour))
