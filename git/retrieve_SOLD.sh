#!/bin/bash
#SBATCH --partition=single
#SBATCH --time=04:00:00
#SBATCH --mem=10000mb
#SBATCH --job-name=SOLD_small
#SBATCH --constraint=LSDF
#SBATCH --output=/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/uc2_logfiles/SOLD_small.out

#Compute SOLD variable from SSR and SSRC for European region and 1979-2018
#Main usage for composites of Mockert et al. Dunkelflauten research

#rsync /home/qb3829/Documents/HiWi_IAI_20/dig/scripts/data_and_analysis/retrieve_SOLD.sh /home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/

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
#Select time period
for year in {1979..2018}
do
  [ ! -d "${dir_out_variables}/${year}" ] && mkdir "${dir_out_variables}/${year}"
  for month in {01..12}
  do
    echo ${year}${month}
    [ ! -d "${dir_out_variables}/${year}/${month}" ] && mkdir "${dir_out_variables}/${year}/${month}"
      #Select variable, region, make daily sum of SSR and SSRC
#      cdo -daysum -sellonlatbox,-20,40,30,70 -select,name=SSR "${dir_era5}/${year}/${month}/N${year}${month}*_*" "${dir_out_variables}/${year}/${month}/SSRdsum_${year}${month}"
#      cdo -daysum -sellonlatbox,-20,40,30,70 -select,name=SSRC "${dir_era5}/${year}/${month}/N${year}${month}*_*" "${dir_out_variables}/${year}/${month}/SSRCdsum_${year}${month}"
      #Compute daily SOLD
      cdo -chname,SSR,SOLD -div "${dir_out_variables}/${year}/${month}/SSRdsum_${year}${month}" "${dir_out_variables}/${year}/${month}/SSRCdsum_${year}${month}" "${dir_out_variables}/${year}/${month}/SOLD_${year}${month}daily"
       #Convert SOLD from daily to 3hourly
       for time in 00 03 06 09 12 15 18 21
       do
         cdo settime,${time}:00:00 "${dir_out_variables}/${year}/${month}/SOLD_${year}${month}_daily" "${dir_out_variables}/${year}/${month}/SOLD_${year}${month}_${time}"
       done
       cdo mergetime "${dir_out_variables}/SOLD_${year}${month}_*" "${dir_out_variables}/SOLD_${year}${month}"
       rm ${dir_out_variables}/SOLD_${year}${month}_*
  done #month
done #year
echo 'create variable done'

#Join SOLD of all dates into one file
[ ! -d "${dir_out_variables}/absolute" ] && mkdir "${dir_out_variables}/absolute"
[ -f "${dir_out_variables}/SOLD" ] && rm "${dir_out_variables}/SOLD"
cdo -mergetime "${dir_out_variables}/*/*/SOLD_*" "${dir_out_variables}/absolute/SOLD"
echo 'merge done'

#Create running_days running climatology
[ ! -d "${dir_out_variables}/climatology" ] && mkdir "${dir_out_variables}/climatology"
running_days=30
running_steps=$(( ${running_days} * 8 + 1))
CDO_TIMESTAT_DATE=middle
cdo ydaymean -runmean,${running_steps} "${dir_out_variables}/absolute/SOLD" "${dir_out_variables}/climatology/SOLD_${running_days}d_climatology"
echo 'climatology done'

#Create Anomaly, subtract from each month file individual by looping over year and month
for year in {1979..2018}
do
  for month in {01..12}
  do
    cdo -sub "${dir_out_variables}/${year}/${month}/SOLD_${year}${month}" "${dir_out_variables}/climatology/SOLD_${running_days}d_climatology" "${dir_out_variables}/${year}/${month}/SOLD${running_days}rmean_clim_ano_${year}${month}"
  done #month
done #year

#Merge anomalies
[ ! -d "${dir_out_variables}/anomaly" ] && mkdir "${dir_out_variables}/anomaly"
rm ${dir_out_variables}/SOLD_${running_days}runmean_clim_ano_daily
cdo -mergetime "${dir_out_variables}/*/*/SOLD_${running_days}rmean_clim_ano_*" "${dir_out_variables}/anomaly/SOLD_${running_days}d_rmean_clim_anomaly" 
echo 'anomaly done'