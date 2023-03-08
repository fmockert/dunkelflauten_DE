#!/bin/bash
#SBATCH --partition=single
#SBATCH --time=04:00:00
#SBATCH --mem=10000mb
#SBATCH --job-name=SLP
#SBATCH --constraint=LSDF
#SBATCH --output=/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/uc2_logfiles/SLP.out

#Compute SLP variable for European region and 1979-2018
#Main usage for composites of Mockert et al. Dunkelflauten research

#rsync /home/qb3829/Documents/HiWi_IAI_20/dig/scripts/data_and_analysis/retrieve_SLP.sh /home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/

#directories
#from worklaptop
dir_era5="/home/qb3829/Documents/PhD/mnt/ec.era5/" 
dir_out_variables="/home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/" 

#from uc2
dir_era5="/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/ec.era5/"
dir_out_variables="/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/dunkelflauten_variables/"

##Select time period
for year in {1979..2018}
do
  [ ! -d "${dir_out_variables}/${year}" ] && mkdir "${dir_out_variables}/${year}"
  for month in {01..12}
  do
    echo ${year}${month}
    [ ! -d "${dir_out_variables}/${year}/${month}" ] && mkdir "${dir_out_variables}/${year}/${month}"
      #Select variable, region
      cdo -sellonlatbox,-20,40,30,70 -select,name=MSL "${dir_era5}/${year}/${month}/N${year}${month}*_*" "${dir_out_variables}/${year}/${month}/SLP_${year}${month}"
  done #month
done #year
echo 'loops done'
#
##Join SLP of all dates into one file
cdo -mergetime "${dir_out_variables}/*/*/SLP_**" "${dir_out_variables}/SLP_3hourly"
echo 'merge done'

#Anomaly not needed for SLP
running_days=30
running_steps=$(( ${running_days} * 8 )) #8 timesteps per day!
##Create running_days running climatology
#cdo ydaymean -runmean,${running_steps} "${dir_out_variables}/SLP_3hourly" "${dir_out_variables}/SLP_3hourly_${running_days}d_climatology"
##Subtract from each month file individual by looping over year and month
#for year in {1979..2018}
#do
#  for month in {01..12}
#  do
#    cdo -sub "${dir_out_variables}/${year}/${month}/SLP_${year}${month}" "${dir_out_variables}/SLP_3hourly_${running_days}d_climatology" "${dir_out_variables}/${year}/${month}/SLP${running_days}runmean_clim_ano_${year}${month}"
#  done #month
#done #year
#
##Merge anomalies
#cdo -mergetime "${dir_out_variables}/*/*/SLP_${running_days}runmean_clim_ano_**" "${dir_out_variables}/SLP_${running_days}runmean_clim_ano_3hourly" 
#
#