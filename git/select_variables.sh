#!/bin/bash
##SBATCH --partition=single
##SBATCH --time=20:00:00
##SBATCH --mem=10000mb
##SBATCH --job-name=select
##SBATCH --constraint=LSDF
##SBATCH --output=/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/uc2_logfiles/select.out

#Compute variables for European region and 1979-2018
#Main usage for composites of Mockert et al. Dunkelflauten research

#rsync /home/qb3829/Documents/HiWi_IAI_20/dig/scripts/data_and_analysis/select_variables.sh /home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/

#directories !get it from the main_composite_workflow.sh
if [ "$HOSTNAME" = "imk-tnn-famo" ]; then
    #from worklaptop
    export dir_era5="/home/qb3829/Documents/PhD/mnt/ec.era5/" 
    export dir_out_variables="/home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/"
    export dir_dunkelflauten="/home/qb3829/Documents/HiWi_IAI_20/dig/outputs/wr/" 
else
    #from uc2
    export dir_era5="/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/ec.era5/"
    export dir_out_variables="/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/dunkelflauten_variables/"
fi

##do this on LSDF!
#for variable in "U100" "V100" "SSR" "SSRC" "T2M" "SLP" #"W100M" "SOLD" 
#do
#  for year in {2019..2019} # {1979..2018}
#  do
#    for month in {01..12} # {01..12}
#    do
#      mkdir ${dir_out_variables}/${year}/${month}/${variable}
#      mv ${dir_out_variables}/${year}/${month}/${variable}_* ${dir_out_variables}/${year}/${month}/${variable}
#    done
#  done
#done

#Generate timestamps of 40years
python "timestamps.py" #${dir_out_variables}/timestamps.py"

##Generate variable files
#while read -r date; 
#do
#  year=${date:0:4}; month=${date:4:2}; day=${date:6:2}; hour=${date:9:2}
#  echo "$year/$month/$day/$hour"
#  [ ! -d "${dir_out_variables}/${year}" ] && mkdir "${dir_out_variables}/${year}"
#  [ ! -d "${dir_out_variables}/${year}/${month}" ] && mkdir "${dir_out_variables}/${year}/${month}"
#  if [[ -f "${dir_era5}/${year}/${month}/N${year}${month}${day}_${hour}" ]]; then
#    for variable in "W100M" "SOLD" "T2M" "SLP"
#    do
#      echo $variable
#       [ ! -d "${dir_out_variables}/${year}/${month}/${variable}" ] && mkdir "${dir_out_variables}/${year}/${month}/${variable}"
#      #SOLD
#      if [ "$variable" = "SOLD" ]; then
#        if [[ ! -f "${dir_out_variables}/${year}/${month}/SSR_${year}${month}${day}_${hour}" ]]; then
#          cdo -select,name=SSR "${dir_era5}/${year}/${month}/N${year}${month}${day}_${hour}" "${dir_out_variables}/${year}/${month}/SSR/SSR_${year}${month}${day}_${hour}"
#        fi
#        if [[ ! -f "${dir_out_variables}/${year}/${month}/SSRC_${year}${month}${day}_${hour}" ]]; then
#          cdo -select,name=SSRC "${dir_era5}/${year}/${month}/N${year}${month}${day}_${hour}" "${dir_out_variables}/${year}/${month}/SSRC/SSRC_${year}${month}${day}_${hour}"
#        fi
#      #T2M
#      elif [ "$variable" = "T2M" ]; then
#        if [[ ! -f "${dir_out_variables}/${year}/${month}/T2M_${year}${month}${day}_${hour}" ]]; then
#          cdo -select,name=T2M "${dir_era5}/${year}/${month}/N${year}${month}${day}_${hour}" "${dir_out_variables}/${year}/${month}/T2M/T2M_${year}${month}${day}_${hour}"
#        fi
#      #W100M
#      elif [ "$variable" = "W100M" ]; then
#        if [[ ! -f "${dir_out_variables}/${year}/${month}/U100_${year}${month}${day}_${hour}" ]]; then
#          cdo -select,name=U100 "${dir_era5}/${year}/${month}/N${year}${month}${day}_${hour}" "${dir_out_variables}/${year}/${month}/U100/U100_${year}${month}${day}_${hour}"
#        fi
#        if [[ ! -f "${dir_out_variables}/${year}/${month}/V100_${year}${month}${day}_${hour}" ]]; then
#          cdo -select,name=V100 "${dir_era5}/${year}/${month}/N${year}${month}${day}_${hour}" "${dir_out_variables}/${year}/${month}/V100/V100_${year}${month}${day}_${hour}"
#        fi
#      #SLP
#      elif [ "$variable" = "SLP" ]; then
#        if [[ ! -f "${dir_out_variables}/${year}/${month}/SLP_${year}${month}${day}_${hour}" ]]; then
#          cdo -select,name=MSL "${dir_era5}/${year}/${month}/N${year}${month}${day}_${hour}" "${dir_out_variables}/${year}/${month}/SLP/SLP_${year}${month}${day}_${hour}"    
#        fi
#      fi #variable selection  
#    done #variable 
#  fi # N file exists
#done #<"${dir_out_variables}/timestamps.txt" #read in dates
#
#