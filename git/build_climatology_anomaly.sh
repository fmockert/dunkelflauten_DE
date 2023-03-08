#!/bin/bash
#SBATCH --partition=single
#SBATCH --time=10:00:00
#SBATCH --mem=10000mb
#SBATCH --job-name=seasonano
#SBATCH --array=0,1,2,3 #["T2M", "W100M", "SOLD", "SLP"]
#SBATCH --constraint=LSDF
#SBATCH --output=/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/uc2_logfiles/seasonano.out

#Compute variables for European region and 1979-2018
#Main usage for composites of Mockert et al. Dunkelflauten research

#rsync /home/qb3829/Documents/HiWi_IAI_20/dig/scripts/data_and_analysis/select_variables.sh /home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/

##directories !Get from main_composite_workflow.sh
#if [ "$HOSTNAME" = "imk-tnn-famo" ]; then
#    #from worklaptop
#    export dir_era5="/home/qb3829/Documents/PhD/mnt/ec.era5/" 
#    export dir_out_variables="/home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/"
#    export dir_dunkelflauten="/home/qb3829/Documents/HiWi_IAI_20/dig/outputs/wr/" 
#else
#    #from uc2
#    export dir_era5="/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/ec.era5/"
#    export dir_out_variables="/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/dunkelflauten_variables/"
#fi

export season='NDJFM'
export running_days=30
export clim_type='season' #'30d_rmean' 'NDJFM'
export fields_to_compute='anomalies'  #'hourly_climatology;season_climatology;running_climatology;anomalies'
python build_climatology_anomaly.py