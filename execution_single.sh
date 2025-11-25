#!/bin/bash

##----------------------- Start job description -----------------------
#SBATCH --partition=standard
#SBATCH --job-name=bh_array_execution
#SBATCH --array=1-1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=60:15:00
#SBATCH --time-min=60:15:00
#SBATCH --mail-type=END
#SBATCH --mail-user=samuel.martinez@upm.es
##------------------------ End job description ------------------------

module purge
module load Gurobi 
module load MATLAB

# ----------------------- Simulation parameters -----------------------
use_case="iridium"
Hcap=10 
m_continuous=1    # 1: continuous, 0: discrete m 
h3_resolution=3
r0=3
rmax=3
d_threshold=5000

# MIPGap=0.04 RES2
MIPGap=0.2

# ----------------------- Define simulation ranges --------------------
# POWER VALUES (12) 
P_T_values=750
#P_T_values=(1500 2500 2750)
 
# SCENARIOS (6)
SCENARIOS=('D')
#SCENARIOS=('D')
num_s=${#SCENARIOS[@]}
num_pt=${#P_T_values[@]}

# ----------------------- Compute indices -----------------------------
task_id=$((SLURM_ARRAY_TASK_ID - 1))

s_index=$(( task_id / num_pt ))
pt_index=$(( task_id % num_pt ))

scenario=${SCENARIOS[$s_index]}
pt=${P_T_values[$pt_index]}

# ----------------------- Run MATLAB function -------------------------
matlab -nosplash -nojvm -nodisplay -r "BH_main_fixed_normalization('$scenario','$use_case',$h3_resolution,$r0,$rmax,$d_threshold,$Hcap,$pt,$m_continuous,$MIPGap); exit"

# ----------------------- Sync results --------------------------------
module load rclone

#rclone copy /home/w384/w384256/SCC_S2C_BH/Output_Data onedrive:SCC_S2C_BH_Results
