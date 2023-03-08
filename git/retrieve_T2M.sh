#!/bin/bash
#SBATCH --partition=single
#SBATCH --time=04:00:00
#SBATCH --mem=10000mb
#SBATCH --job-name=T2M
#SBATCH --constraint=LSDF
#SBATCH --output=/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/uc2_logfiles/T2M.out

#Compute T2M variable for European region and 1979-2018
#Main usage for composites of Mockert et al. Dunkelflauten research

#rsync /home/qb3829/Documents/HiWi_IAI_20/dig/scripts/data_and_analysis/retrieve_T2M.sh /home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/

#directories
if [ "$HOSTNAME" = "imk-tnn-famo" ]; then
    #from worklaptop
    dir_era5="/home/qb3829/Documents/PhD/mnt/ec.era5/" 
    dir_out_variables="/home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/"
    dir_dunkelflauten="/home/qb3829/Documents/HiWi_IAI_20/dig/outputs/wr/" 
else
    #from uc2
    dir_era5="/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/ec.era5/"
    dir_out_variables="/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/dunkelflauten_variables/"
fi

##Select time period
for year in {1979..2018}
do
  [ ! -d "${dir_out_variables}/${year}" ] && mkdir "${dir_out_variables}/${year}"
  for month in {01..12}
  do
    echo ${year}${month}
    [ ! -d "${dir_out_variables}/${year}/${month}" ] && mkdir "${dir_out_variables}/${year}/${month}"
      #Select variable, region
      cdo -sellonlatbox,-20,40,30,70 -select,name=T2M "${dir_era5}/${year}/${month}/N${year}${month}*_*" "${dir_out_variables}/${year}/${month}/T2M_${year}${month}"
  done #month
done #year
echo 'loops done'
#
##Join T2M of all dates into one file
[ -f "${dir_out_variables}/T2M_daily" ] && rm "${dir_out_variables}/T2M_daily"
cdo -mergetime "${dir_out_variables}/*/*/T2M_**" "${dir_out_variables}/T2M_3hourly"
echo 'merge done'

#Create running_days running climatology
running_days=30
running_steps=$(( ${running_days} * 8 +1 )) #8 timesteps per day!
cdo ydaymean -runmean,${running_steps} "${dir_out_variables}/T2M_3hourly" "${dir_out_variables}/T2M_3hourly_${running_days}d_climatology"

#Create Anomaly, subtract from each month file individual by looping over year and month
for year in {1979..2018}
do
  for month in {01..12}
  do
    cdo -sub "${dir_out_variables}/${year}/${month}/T2M_${year}${month}" "${dir_out_variables}/T2M_3hourly_${running_days}d_climatology" "${dir_out_variables}/${year}/${month}/T2M${running_days}runmean_clim_ano_${year}${month}"
  done #month
done #year

#Merge anomalies
cdo -mergetime "${dir_out_variables}/*/*/T2M_${running_days}runmean_clim_ano_**" "${dir_out_variables}/T2M_${running_days}runmean_clim_ano_3hourly" 

