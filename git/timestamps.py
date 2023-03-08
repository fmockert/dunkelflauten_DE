#Generate a time series with YYYYMMDD_HH of all times (40years)
import pandas as pd
import os
#rsync /home/qb3829/Documents/HiWi_IAI_20/dig/scripts/data_and_analysis/timestamps.py /home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/

dir_era5=os.getenv('dir_era5')
dir_out_variables=os.getenv('dir_out_variables')

times = pd.date_range('1979-01-01 00:00:00', '2022-12-31 21:00:00', freq='3h')
times = pd.Series(times.strftime('%Y%m%d_%H'))
print(dir_out_variables)
times.to_csv('%s/timestamps.txt'%(dir_out_variables), header=None, index=None)

