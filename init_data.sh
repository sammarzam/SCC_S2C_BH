#!/bin/bash

##----------------------- Start job description -----------------------
#SBATCH --partition=standard
#SBATCH --job-name=copy_input_data
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:20:00
#SBATCH --mail-type=END
#SBATCH --mail-user=samuel.martinez@upm.es
##------------------------ End job description ------------------------

module purge

module load rclone

rclone copy onedrive:SCC_S2C_BH_Data /home/w384/w384256/SCC_S2C_BH/Input_Data

rclone copy onedrive:SCC_S2C_BH_Normalization_Data /home/w384/w384256/SCC_S2C_BH/Normalization_Data
