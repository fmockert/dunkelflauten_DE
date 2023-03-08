import xarray as xr
import pandas as pd
import os
import socket
import xarray.ufuncs as xu #deprecated, use np.ufuncs directly
import numpy as np

if socket.gethostname() == "imk-tnn-famo":
    #from worklaptop
    dir_era5="/home/qb3829/Documents/PhD/mnt/ec.era5/"
    dir_out_variables="/home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/"
    dir_dunkelflauten="/home/qb3829/Documents/HiWi_IAI_20/dig/outputs/wr/" 
else:
    #from uc2
    dir_era5="/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/ec.era5/"
    dir_out_variables="/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/dunkelflauten_variables/"

idx = int(os.environ["SLURM_ARRAY_TASK_ID"])
variable = ["T2M", "W100M", "SOLD", "SLP"][idx]
#variable=os.getenv('variable')
season=os.getenv('season')
running_days=int(os.getenv('running_days'))
clim_type=os.getenv('clim_type')
fields_to_compute=os.getenv('fields_to_compute').split(',')

dates = pd.date_range('1979-01-01T00:00:00', '1979-12-31T21:00:00', freq='3h')

#DAILY CLIMATOLOGY
if 'hourly_climatology' in fields_to_compute:
  if os.path.isdir('%s/climatology'%(dir_out_variables))==False: os.system('mkdir %s/climatology'%(dir_out_variables))
  if os.path.isdir('%s/climatology/hourly'%(dir_out_variables))==False: os.system('mkdir %s/climatology/hourly'%(dir_out_variables))
  if os.path.isdir('%s/climatology/hourly/%s'%(dir_out_variables, variable))==False: os.system('mkdir %s/climatology/hourly/%s'%(dir_out_variables, variable))
  
  #Go through 3hourly date list of one representative year
  for date in dates.strftime('%m%d_%H'):
    #iterate over years to get the 40 dates to join
    #generate the file names and select existing
    file_names = ['%s/%s/%s/%s/%s_%s%s%s_%s'%(dir_out_variables, year, date[0:2], variable, variable, year, date[0:2], date[2:4], date[5:7]) for year in np.arange(1979,2018+1) if os.path.isfile('%s/%s/%s/%s/%s_%s%s%s_%s'%(dir_out_variables, year, date[0:2], variable, variable, year, date[0:2], date[2:4], date[5:7]))]
    #load files and replace inf values
    fields = xr.open_mfdataset(file_names)
    fields = fields.where(fields.apply(np.isfinite)).fillna(np.nan)
    #take mean over time and save
    fields_mean = fields.mean(dim='time', skipna=True) 
    # ADD DATE INFORMATION TO NC FILE (YEAR ARBITRARY EG 1979)
    fields_mean = fields_mean.expand_dims(time=[pd.Timestamp('%s-%s-%sT%s:00:00'%(1979, date[0:2], date[2:4], date[5:7]))])
    fields_mean.to_netcdf('%s/climatology/hourly/%s/%s_%s_clim'%(dir_out_variables, variable, variable, date))
  #Generate the 29th of February Climatology, daily climatology
  for hour in ['00', '03', '06', '09', '12', '15', '18', '21']:
    clim_0228 = xr.open_dataset('%s/climatology/hourly/%s/%s_%s_%s_clim'%(dir_out_variables, variable, variable, '0228', hour))
    clim_0228['time'] = [pd.Timestamp('1980-02-29T%s:00:00'%(hour))]
    clim_0228.to_netcdf('%s/climatology/hourly/%s/%s_%s_%s_clim'%(dir_out_variables, variable, variable, '0229', hour))
  

#RUNNING CLIMATOLOGY
if 'running_climatology' in fields_to_compute:
  if os.path.isdir('%s/climatology/%sd_rmean'%(dir_out_variables, running_days))==False: os.system('mkdir %s/climatology/%sd_rmean'%(dir_out_variables, running_days))
  if os.path.isdir('%s/climatology/%sd_rmean/%s'%(dir_out_variables, running_days, variable))==False: os.system('mkdir %s/climatology/%sd_rmean/%s'%(dir_out_variables, running_days, variable))
  #read in start_dates/files
  dates_rolling = pd.date_range(dates[0]-pd.Timedelta(running_days/2*24, 'h'), dates[0]+pd.Timedelta(running_days/2*24, 'h'), freq='3h').strftime('%m%d_%H')
  files_rolling = ['%s/climatology/hourly/%s/%s_%s_clim'%(dir_out_variables, variable, variable, date) for date in dates_rolling if os.path.isfile('%s/climatology/hourly/%s/%s_%s_clim'%(dir_out_variables, variable, variable, date))]
  #files_rolling = files_rolling[:120]+files_rolling[124:] #CHANGE IN THE END AGAIN ONCE PLEVEL CORRECT
  running_field = xr.open_mfdataset(files_rolling)
  
  for date in dates:
    #Generate mean of running field and add time dimension
    running_field_mean = running_field.mean(dim='time', skipna=True).expand_dims(time=[pd.Timestamp('%s-%s-%sT%s:00:00'%(1979, date.strftime('%m'), date.strftime('%d'), date.strftime('%H')))])
    #Save running climatology field
    running_field_mean.to_netcdf('%s/climatology/%sd_rmean/%s/%s_%s%s_%s_%sd_rmean_clim'%(dir_out_variables, running_days, variable, variable, date.strftime('%m'), date.strftime('%d'), date.strftime('%H'), running_days))
  
    #Drop and add dates
    date_drop=date-pd.Timedelta(running_days/2, 'd')
    if pd.Timestamp('1979-%s-%sT%s:00:00'%(date_drop.strftime('%m'), date_drop.strftime('%d'), date_drop.strftime('%H'))) in running_field.time.values: 
      running_field = running_field.drop_sel(time=pd.Timestamp('1979-%s-%sT%s:00:00'%(date_drop.strftime('%m'), date_drop.strftime('%d'), date_drop.strftime('%H'))))
    date_add = date+pd.Timedelta(running_days/2, 'd')+pd.Timedelta(3, 'h')
    #Controlling if file is available (every day should be available for hourly clim)?
    running_field = xr.concat([running_field, xr.open_dataset('%s/climatology/hourly/%s/%s_%s_clim'%(dir_out_variables, variable, variable, date_add.strftime('%m%d_%H')))], dim='time')
  
  #Generate the 29th of February Climatology, running_climatology
  for hour in ['00', '03', '06', '09', '12', '15', '18', '21']:
    clim_0228 = xr.open_dataset('%s/climatology/%sd_rmean/%s/%s_%s_%s_%sd_rmean_clim'%(dir_out_variables, running_days, variable, variable, '0228', hour, running_days))
    clim_0228['time'] = [pd.Timestamp('1980-02-29T%s:00:00'%(hour))]
    clim_0228.to_netcdf('%s/climatology/%sd_rmean/%s/%s_%s_%s_%sd_rmean_clim'%(dir_out_variables, running_days, variable, variable, '0229', hour, running_days))

#SEASON CLIMATOLOGY
print('The clim_type is %s, and the season is: %s'%(clim_type, season))
if 'season_climatology' in fields_to_compute:
  if os.path.isdir('%s/climatology/%s'%(dir_out_variables, clim_type))==False: os.system('mkdir %s/climatology/%s'%(dir_out_variables, clim_type))
  if os.path.isdir('%s/climatology/%s/%s/%s'%(dir_out_variables, clim_type, variable))==False: os.system('mkdir %s/climatology/%s/%s/%s'%(dir_out_variables, clim_type, variable))
  season_dates = {'NDJFM': pd.date_range('1979-11-01T00:00:00', '1980-03-31T21:00:00', freq='3h'), #['11', '12', '01', '02', '03'],
                   'DJF': pd.date_range('1979-12-01T00:00:00', '1980-02-29T21:00:00', freq='3h'), #['12', '01', '02'],
                   'SON': pd.date_range('1979-09-01T00:00:00', '1979-11-30T21:00:00', freq='3h')} #['09', '10', '11']}

  #read in dates/files
  files_season = ['%s/climatology/hourly/%s/%s_%s_clim'%(dir_out_variables, variable, variable, date) for date in season_dates[season].strftime('%m%d_%H') if os.path.isfile('%s/climatology/hourly/%s/%s_%s_clim'%(dir_out_variables, variable, variable, date))]
  #files_rolling = files_rolling[:120]+files_rolling[124:] #CHANGE IN THE END AGAIN ONCE PLEVEL CORRECT
  season_field = xr.open_mfdataset(files_season)
  
  #Generate mean of season field and add arbitrary time (first time of season)
  season_field_mean = season_field.mean(dim='time', skipna=True).expand_dims(time=[season_dates[season].strftime('%m%d_%H')[0]])
  #Save seasonal climatology field
  season_field_mean.to_netcdf('%s/climatology/%s/%s/%s/%s_%s_clim'%(dir_out_variables, clim_type, variable, variable, season))

#CREATE ANOMALIES
if 'anomalies' in fields_to_compute:
  if os.path.isdir('%s/anomaly/%s'%(dir_out_variables, clim_type))==False: os.system('mkdir %s/anomaly/%s'%(dir_out_variables, clim_type))
  if os.path.isdir('%s/anomaly/%s/%s'%(dir_out_variables, clim_type, variable))==False: os.system('mkdir %s/anomaly/%s/%s'%(dir_out_variables, clim_type, variable))

  if 'rmean' in clim_type:
    for date in dates.append(pd.date_range('1980-02-29T00:00:00', '1980-02-29T21:00:00', freq='3h')): #adding leap date
      #Load the climatology files one after the other through the representative year  
      clim_field = xr.open_dataset('%s/climatology/%sd_rmean/%s/%s_%s_%s_clim'%(dir_out_variables, running_days, variable, variable, date.strftime('%m%d_%H'), clim_type))
      #Load the absolute fields for all years for that climatology date
      for year in np.arange(1979, 2018+1):
        if ((date.month==2) & (date.day==29) & (year%4!=0)): continue #excluding nonexisting leap year/date combinations
        absolute_field = xr.open_dataset('%s/%s/%s/%s/%s_%s%s'%(dir_out_variables, year, date.strftime('%m'), variable, variable, year, date.strftime('%m%d_%H')), engine="netcdf4")
        #Change the date of the climatology file according to the year so no problem with time occurs
        clim_field_tmp = clim_field.copy()
        clim_field_tmp['time'] = [absolute_field.time.values[0]]
        #Subtract the climatology from the absolute field
        anomaly_field = absolute_field-clim_field_tmp
        #Save absolute field
        anomaly_field.to_netcdf('%s/anomaly/%s/%s/%s/%sano_%s%s_%s'%(dir_out_variables, clim_type, season, variable, variable, year, date.strftime('%m%d_%H'), clim_type))

  elif 'season' in clim_type:
    season_dates = {'NDJFM': pd.date_range('1979-11-01T00:00:00', '1980-03-31T21:00:00', freq='3h'), #['11', '12', '01', '02', '03'],
                    'DJF': pd.date_range('1979-12-01T00:00:00', '1980-02-29T21:00:00', freq='3h'), #['12', '01', '02'],
                    'SON': pd.date_range('1979-09-01T00:00:00', '1979-11-30T21:00:00', freq='3h')} #['09', '10', '11']}
    clim_field = xr.open_dataset('%s/climatology/%s/%s/%s_%s_clim'%(dir_out_variables, clim_type, variable, variable, season))
    for date in season_dates[season]:
      for year in np.arange(1979, 2018+1):
        if ((date.month==2) & (date.day==29) & (year%4!=0)): continue #excluding nonexisting leap year/date combinations
        absolute_field = xr.open_dataset('%s/%s/%s/%s/%s_%s%s'%(dir_out_variables, year, date.strftime('%m'), variable, variable, year, date.strftime('%m%d_%H')), engine="netcdf4")
        #Change the date of the climatology file according to the year so no problem with time occurs
        clim_field_tmp = clim_field.copy()
        clim_field_tmp['time'] = [absolute_field.time.values[0]]
        #Subtract the climatology from the absolute field
        anomaly_field = absolute_field-clim_field_tmp
        #Save absolute field
        anomaly_field.to_netcdf('%s/anomaly/%s/%s/%sano_%s%s_%s'%(dir_out_variables, clim_type, variable, variable, year, date.strftime('%m%d_%H'), season))
