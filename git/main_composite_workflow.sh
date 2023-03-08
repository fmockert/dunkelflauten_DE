#!/bin/bash
#SBATCH --partition=single
#SBATCH --time=10:00:00
#SBATCH --mem=10000mb
#SBATCH --job-name=main
#SBATCH --constraint=LSDF
#SBATCH --output=/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/uc2_logfiles/main.out

#need to program in a way that next computation waits till earlier computation is done!

#directories
if [ "$HOSTNAME" = "imk-tnn-famo" ]; then
    #from worklaptop
    export dir_era5="/home/qb3829/Documents/PhD/mnt/ec.era5/" 
    export dir_out_variables="/home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/"
    export dir_dunkelflauten="/home/qb3829/Documents/HiWi_IAI_20/dig/outputs/wr/" 
else
    #from uc2
    export dir_era5="/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/ec.era5/"
    export dir_out_variables="/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/dunkelflauten_variables/"
    export dir_dunkelflauten="${dir_out_variables}/wr/"
fi

#Generate/Get meteorological field variables from ERA5 for 1979-2018
bash select_variables.sh #(also calling timestamps.py to get a list of all timestamps)
#sbatch select_variables.sh #(also calling timestamps.py to get a list of all timestamps)
#
##Calculate the variables that are not directly available by ERA5 (SOLD=SSR/SSRC and W100M=sqrt(U100^2+V100^2))
#sbatch calculate_SOLD_W100M.sh #(calling calculate_SOLD_W100M.py)
#
##Build climatology and anomaly
#sbatch build_climatology_anomaly.sh #(need to select what to calculate hourly_clim, running_clim and anomaly, and select running_days. Calls build_climatology_anomaly.py)
##hourly clim takes 2h10min
##30d rmean clim takes 6h
##anomaly don't know, ran over night
#
##Build composites of WR-Dunkelflauten and WR, needs build_weathermap_times.py from snakefile to run to get the timestamps of the WR-DF and WR
#export season='NDJFM'
#export running_days=30
#export clim_type="season" #"${running_days}d_rmean" #"season"
#for regime in 'eubl' #at zo sctr ar eubl scbl gl all
#do
#  export regime
#  for variable in 'T2M' #T2M W100M SOLD SLP
#  do
#    export variable
#    python build_dunkelflauten_composites.py #quick also locally
#    #python build_weatherregime_composites.py #takes few minutes locally for 1regime and 1field
#  done #variable
#done #regime
#
##Plot the composites
#export season='NDJFM'
#export category='dunkelflaute' #'regime' 'dunkelflaute'
#export running_days=30
#export clim_type="${running_days}d_rmean"
#for regime in all eubl scbl gl #at zo sctr ar eubl scbl gl 
#do
#  export regime
#  for variable in T2M W100M SOLD
#  do
#    export variable
#    python plot_composites.py    
#  done #variable
#done #regime
#
##Build LAGGED composites of WR-Dunkelflauten, needs build_weathermap_times.py from snakefile to run to get the timestamps of the WR-DF
#export season='NDJFM'
#export running_days=30
#export clim_type="${running_days}d_rmean" #"${running_days}d_rmean" #"season"
#export lag_inputs='-6;-4;-2;0;2;4' #;-2;0;2;4'
#export interval_input=2
#
#for regime in 'gl' #at zo sctr ar eubl scbl gl all
#do
#  export regime
#  for variable in 'SLP' #T2M W100M SOLD SLP
#  do
#    export variable
#    python build_dunkelflauten_composites_lagged.py #quick also locally
#  done #variable
#done #regime
##Plot the LAGGED composites
#for regime in 'gl' #all eubl scbl gl #at zo sctr ar eubl scbl gl 
#do
#  export regime
#  for variable in 'T2M' #T2M W100M SOLD
#  do
#    export variable
#    python plot_composites_lagged.py    
#  done #variable
#done #regime