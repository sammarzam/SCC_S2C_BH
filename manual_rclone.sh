#!/bin/bash

export LC_ALL=C  # <<--- IMPORTANT: ensures beta = 0.70 (not 0,70)

OUTPUT_DIR="/home/w384/w384256/SCC_S2C_BH/Output_Data"
REMOTE="onedrive:SCC_S2C_BH_Results"

# ----------------------- Evaluation range -----------------------
betta=0.7
t0_file=1191
t1_file=1200

# ----------------------- Simulation parameters -----------------------
use_case="iridium"
m_continuous=1
h3_resolution=3

# ----------------------- Define simulation ranges --------------------
P_T_values=(250 500 750 1000 1250 1500 1750 2000 2250 2500 2750 3000)
SCENARIOS=(A B C D E F)

for scenario in "${SCENARIOS[@]}"; do
    for P_T in "${P_T_values[@]}"; do
        
        fname=$(printf "%s_BH_[%s_res%d]_P_%d_mC_%d_beta_%0.2f_win_%d_%d.mat" \
                "$scenario" "$use_case" "$h3_resolution" "$P_T" "$m_continuous" \
                "$betta" "$t0_file" "$t1_file")

        fullpath="${OUTPUT_DIR}/${fname}"

        if [[ -f "$fullpath" ]]; then
            echo "Copying $fname ..."
            rclone copy "$fullpath" "$REMOTE"
        else
            echo "Missing file: $fname"
        fi
    done
done
