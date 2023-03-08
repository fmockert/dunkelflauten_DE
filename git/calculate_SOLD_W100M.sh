#!/bin/bash
#SBATCH --partition=single
#SBATCH --time=20:00:00
#SBATCH --mem=10000mb
#SBATCH --job-name=soldw100m
#SBATCH --constraint=LSDF
#SBATCH --array=0-39
#SBATCH --output=/lsdf/kit/imk-tro/projects/MOD/Gruppe_Grams/qb3829/uc2_logfiles/soldw100m.out

#Compute variables for European region and 1979-2018
#Main usage for composites of Mockert et al. Dunkelflauten research

#rsync /home/qb3829/Documents/HiWi_IAI_20/dig/scripts/data_and_analysis/calculate_SOLD_W100M.sh /home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/
#rsync /home/qb3829/Documents/HiWi_IAI_20/dig/scripts/data_and_analysis/calculate_SOLD_W100M.py /home/qb3829/Documents/PhD/mnt/qb3829/dunkelflauten_variables/

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

#SOLD calculation
#W100M calculation
python calculate_SOLD_W100M.py