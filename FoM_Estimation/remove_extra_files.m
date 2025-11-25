clear all
close all


global PWD;
PWD=pwd;

outputDir='C:\Users\DOC06\OneDrive - Universidad Politécnica de Madrid\SCC_S2C_BH_Results';

% ----------------------- Evaluation range -----------------------
betta=0.7;
% t0_file = 1191; 
% t1_file = 1200;

% ----------------------- Simulation parameters -----------------------
use_case="iridium"
Hcap=10;
m_continuous=1;   % 1: continuous, 0: discrete m 
h3_resolution=3;
r0=2;
rmax=2;
d_threshold=5000;

% ----------------------- Define simulation ranges --------------------

% POWER VALUES (12)
P_T_values=[250 500 750 1000 1250 1500 1750 2000 2250 2500 2750 3000] %
% SCENARIOS (6)
SCENARIOS={'A'; 'B'; 'C'; 'D'; 'E'; 'F'}%{'A'; 'B'; 'C'}%

num_s=length(SCENARIOS);
num_pt=length(P_T_values);

W=10;
for t_idx=1:W:1190
    t0_file=t_idx;
    t1_file=t_idx+W-1;

    for s_idx=1:num_s
        scenario=SCENARIOS{s_idx};
        for p_idx=1:num_pt
            P_T=P_T_values(p_idx);
    
            fname = sprintf('%s_BH_[%s_res%d]_P_%d_mC_%d_beta_%0.2f_win_%d_%d.mat',scenario,use_case, h3_resolution, P_T, m_continuous, betta, t0_file, t1_file);
            if isfile(fullfile(outputDir, fname))
                delete((fullfile(outputDir, fname)));
            end
        end
    end
end

