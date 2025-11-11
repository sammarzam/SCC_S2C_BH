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

module load rclone

rclone copy onedrive:TN_NTN_Data_A /home/w384/w384256/TN_NTN_A/Input_Data
