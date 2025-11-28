clc
clear all
close all


global PWD;
PWD=pwd;

inputDir='C:\Users\DOC06\OneDrive - Universidad Politécnica de Madrid\SCC_S2C_BH_Data';
outputDir='C:\Users\DOC06\OneDrive - Universidad Politécnica de Madrid\SCC_S2C_BH_Results_Full';

% ----------------------- Evaluation range -----------------------
betta=0.7;
t0_file = 1191; 
t1_file = 1200;

t0 = 1; 
t1 = 1200;
win = t0:t1;

% ----------------------- Simulation parameters -----------------------
use_case="iridium"
Hcap=10;
m_continuous=1;   % 1: continuous, 0: discrete m 
h3_resolution=2;
r0=2;
rmax=2;
d_threshold=5000;

beams=16;

frame_dur=0.1;

% ----------------------- Define simulation ranges --------------------

% POWER VALUES (12)
%P_T_values=[250 500 750 1000 1250 1500 1750 2000 2250 2500 2750] %3000
P_T_values=[50 150 250 500 750 1000 1500 2000 2500]

% SCENARIOS (6)
SCENARIOS={'A'; 'B'; 'C'; 'D'; 'E'; 'F'}%{'A'; 'B'; 'C'}%

num_s=length(SCENARIOS);
num_pt=length(P_T_values);

% KPI MATRIX INITIALIZATION 
UC_w_ALL=zeros(num_s,num_pt);
UC_g_ALL=zeros(num_s,num_pt);
EC_all=zeros(num_s,num_pt);
TTS_all=zeros(num_s,num_pt);
TTS_act_all=zeros(num_s,num_pt);

Handover_all=zeros(num_s,num_pt);
Lisl_mean_all=zeros(num_s,num_pt);
Lisl_mean_util=zeros(num_s,num_pt);
Lisl_max_all=zeros(num_s,num_pt);
Lisl_max_util=zeros(num_s,num_pt);

Pused_all=zeros(num_s,num_pt);
Bused_all=zeros(num_s,num_pt);
SE_eff_all=zeros(num_s,num_pt);
Lights_all=zeros(num_s,num_pt);
Served_all=zeros(num_s,num_pt);

addpath(fullfile(inputDir,'SAT_ROUTING')) % R for ISL load
load(strcat('EUR_R_K_res_',num2str(h3_resolution),'_',use_case,'.mat'));


for s_idx=1:num_s
    % General SCENARIO parameters:
    scenario=SCENARIOS{s_idx};

    % SCC: B
    addpath(fullfile(inputDir,'CELL_COLOURING_MAX_BAND_ALLOCATION_REPAIR_FORWARD'))
    tag = sprintf('r0%d_rmax%d_thr%g', r0, rmax, d_threshold); % y, Btot, K
    load(['EUR_B_h3_adaptive_' tag '_V5.mat']);
    B=y;
    B_T=Btot;
    N=K;
    
    if scenario=='A' || scenario=='B' || scenario=='C' % UNIFORM Band distribution
        B=B.*0;
        reuse_factor=3;
    
        band_bins_uniform=floor((B_T/reuse_factor)/(B_T/N));
        % band_bins_uniform=12
        for idx=1:size(B,1)
            B(idx,1:band_bins_uniform)=1;
        end
    end
    clear K
    clear Btot

    % S2C: Demand d, sol, beams
    addpath(fullfile(inputDir,'SAT_CELL_ASSOCIATION'))
    if scenario=='A' || scenario=='D' % CLOSEST
        load(strcat('EUR_X_res_',num2str(h3_resolution),'_',use_case,'_beams_',num2str(beams),'_closest.mat')); % sol.x, associate satellite s to cell i at time t
        X=sol.x;
    end
    if scenario=='B' || scenario=='E' % DEMAND ONLY
         load(strcat('EUR_X_res_',num2str(h3_resolution),'_',use_case,'_beams_',num2str(beams),'_d_only.mat')); % sol.x, associate satellite s to cell i at time t
        X=sol.x;
    end

    if scenario=='C' || scenario=='F' % DEMAND + ISL+ HANDOVERS
         load(strcat('EUR_X_res_',num2str(h3_resolution),'_',use_case,'_beams_',num2str(beams),'_d.mat')); % sol.x, associate satellite s to cell i at time t
        X=sol.x;
    end

    X=(X==1); % logical!
    
    % Load pre-BH
    fname = sprintf('%s_pre_BH_[%s_res%d_beams%d]_data.mat', scenario, use_case, h3_resolution,beams);
    load(fullfile(outputDir, fname));

    for p_idx=1:num_pt
        % Specific POWER file:
        P_T=P_T_values(p_idx);

        fname = sprintf('%s_BH_[%s_res%d_beams%d]_P_%d_mC_%d_beta_%0.2f_win_%d_%d.mat',scenario,use_case, h3_resolution,beams, P_T, m_continuous, betta, t0_file, t1_file);
        load(fullfile(outputDir, fname)); % i,m,g,f,c,d(?)
     
        % KPI calculation:
       
        %% ===============================================================
        %  Assumes in workspace:
        %    - R_D  [cells x Hf]  : requested demand (inner time)
        %    - d    [cells x Hf]  : served demand (inner time)
        %    - g    [cells x Hf]  : backlog ledger (inner time)
        %    - f    [cells x Hf]  : oversupply ledger (inner time)
        %    - c    [cells x Hf]  : time-to-serve counter (inner time)
        %    - X    [cells x nSats x H] : responsibility (OUTER time)
        %    - frame (scalar)     : inner slots per outer slot
        %    - t0, t1 (scalars)   : inner-time window (inclusive)
        % ===============================================================
        
        [C, Hf] = size(R_D);
        W = numel(win);
        
        %% -------------------------------
        %  UC (global, demand-weighted)
        Rcum = cumsum(R_D, 2);        % [C x Hf]
        den  = Rcum(:, win);          % denominator for r_{c,t}
        num  = g(:,    win);          % numerator for r_{c,t}
        mask = den > 0;
        
        r_inst = zeros(C, W);
        r_inst(mask) = num(mask) ./ den(mask);
        
        w = R_D(:, win);              % demand weights per slot (zero when no demand)
        w(~mask) = 0;
        
        UC_inst_weighted      = sum(mean(r_inst,2).*sum(R_D,2))./sum(R_D,'all');     % MAIN UC KPI WEIGHTED
        UC_inst_global        = mean(mean(r_inst,2));     % MAIN UC KPI non-WEIGHTED
        
        %% -------------------------------
        %  EC 
        EC_global     = sum(f(:, win), 'all') / max(sum(R_D(:, win), 'all'), 1);
        
        %% -------------------------------
        %  TTS 
        TTS_mean_c         = mean(c(:, win), 2);
        TTS_active_rate_c  = sum(c(:, win) > 0, 2) / W;
        
        %% -------------------------------
        %  Handover metrics from X (responsibility)
        %  X is OUTER time -> expand to inner time via frame
        % -------------------------------
        if exist('X','var')==1 && exist('frame','var')==1
            X_hf  = repelem(X, 1, 1, frame);      % [cells x nSats x Hf]
            X_win = X_hf(:, :, win);
            [~, s_win] = max(X_win, [], 2);       % [cells x 1 x W]
            s_win = squeeze(s_win);               % [cells x W]
        
            handover_c = (W>1) .* sum(s_win(:,2:end) ~= s_win(:,1:end-1), 2);
            handover_rate_c = zeros(C,1);
            if W > 1, handover_rate_c = handover_c / (W-1); end
        
            avg_dwell_c = zeros(C,1);
            for cc = 1:C
                seq = s_win(cc,:);
                if all(seq==seq(1))
                    avg_dwell_c(cc) = W;
                else
                    b  = [true, diff(seq)~=0, true];
                    re = find(b);
                    rl = diff(re);
                    avg_dwell_c(cc) = mean(rl);
                end
            end
        
            HO_mean       = mean(handover_c);
            HO_rate_mean  = mean(handover_rate_c);
            HO_dwell_mean = mean(avg_dwell_c);
        else
            warning('X/frame not found — skipping handover metrics.');
            handover_c = []; handover_rate_c = []; avg_dwell_c = [];
            HO_mean = NaN; HO_rate_mean = NaN; HO_dwell_mean = NaN;
        end
        
        %% ISL LOAD CALCULATION FOR Hcap <- R, X, and d (given by BH)
        d_bh=squeeze(sum(reshape(d, size(d,1), frame, []), 2))./(frame*frame_dur);
        [nV, Kmax, H] = size(X);
        [~, Hcap_check] = size(d_bh);
        if Hcap_check ~= Hcap
            error('d_bh has %d columns but Hcap=%d', Hcap_check, Hcap);
        end
        
        Cisl=4*2.0e3.*ones(Kmax,H); % per-(k,t) ISL capacity <- fixed (ISL terminal)
    
        nISL_per_sat  = 4;
        Cisl_total_sat = Cisl;                    % total per satellite
        cap_perISL    = Cisl_total_sat / nISL_per_sat;   % per-ISL capacity (scalar)
        
        % ISL load per satellite and time, restricted to BH slots 1..Hcap
        Lisl_bh = zeros(Kmax, Hcap);              % Kmax x Hcap
        
        for tt = 1:Hcap
            Rt = R(:,:,tt);                       % Kmax x Kmax, same t index
            [kNZ, sNZ, rcoeff] = find(Rt);        % nonzeros of R(:, :, tt)
            if isempty(kNZ), continue; end
        
            for idx = 1:numel(kNZ)
                k   = kNZ(idx);                   % ISL node (satellite) index
                s   = sNZ(idx);                   % serving satellite index
                rks = rcoeff(idx);                % R(k,s,tt)
        
                % cells i served by satellite s at time tt
                cells_s = find(X(:,s,tt) > 0.5);
                if isempty(cells_s), continue; end
        
                % total demand served by satellite s at time tt (using BH demand)
                dem_s = sum(d_bh(cells_s, tt));   % same units as Cisl
        
                % contribution to ISL load on satellite k
                Lisl_bh(k, tt) = Lisl_bh(k, tt) + rks * dem_s;
            end
        end
        
        % Mean load per ISL (divide satellite total load by 4)
        meanLoad_perISL_bh = Lisl_bh / nISL_per_sat;    % Kmax x Hcap
        
        % Per-ISL utilization for each satellite/time in BH window
        util_perISL_bh = meanLoad_perISL_bh ./ cap_perISL(:,1:Hcap);   % Kmax x Hcap
        
        %% ===== NON-ZERO LOAD FILTERING =====
        % Satellites that ever carry *any* ISL load over BH
        sat_has_load = any(meanLoad_perISL_bh > 0, 2);   % Kmax x 1 logical
        
        % Time–ISL points that actually carry load
        nz_mask = meanLoad_perISL_bh > 0;                % Kmax x Hcap logical
        
        % --- Per-satellite stats BUT only over sats that have load ---
        meanLoad_sat_all = mean(meanLoad_perISL_bh, 2);       % as before (includes zeros)
        maxLoad_sat_all  = max(meanLoad_perISL_bh, [], 2);
        
        meanUtil_sat_all = mean(util_perISL_bh, 2);
        maxUtil_sat_all  = max(util_perISL_bh, [], 2);
        
        % Restrict to active satellites
        meanLoad_sat_nz = meanLoad_sat_all(sat_has_load);
        maxLoad_sat_nz  = maxLoad_sat_all(sat_has_load);
        
        meanUtil_sat_nz = meanUtil_sat_all(sat_has_load);
        maxUtil_sat_nz  = maxUtil_sat_all(sat_has_load);
        
        % --- Global mean/max over non-zero ISL samples only ---
        if any(nz_mask(:))
            nz_load = meanLoad_perISL_bh(nz_mask);   % vector of >0 loads
            nz_util = util_perISL_bh(nz_mask);       % vector of >0 utils
        
            meanLoad_global_nz = mean(nz_load);
            maxLoad_global_nz  = max(nz_load);
        
            meanUtil_global_nz = mean(nz_util);
            maxUtil_global_nz  = max(nz_util);
        else
            meanLoad_global_nz = 0;
            maxLoad_global_nz  = 0;
            meanUtil_global_nz = 0;
            maxUtil_global_nz  = 0;
        end
        
        % %% ===== PRINT NON-ZERO SUMMARY =====
        % fprintf('\n===== ISL BUSY-HOUR SUMMARY (only non-zero loads) =====\n');
        % fprintf('Mean load   (global, >0): %.2f Mbps\n', meanLoad_global_nz);
        % fprintf('Max  load   (global, >0): %.2f Mbps\n', maxLoad_global_nz);
        % fprintf('Mean util   (global, >0): %.2f\n',     meanUtil_global_nz);
        % fprintf('Max  util   (global, >0): %.2f\n\n',   maxUtil_global_nz);
        %% ===== SAVE NON-ZERO SUMMARY =====
        Lisl_mean_all(s_idx,p_idx)=meanLoad_global_nz;
        Lisl_mean_util(s_idx,p_idx)=meanUtil_global_nz;
        Lisl_max_all(s_idx,p_idx)=maxLoad_global_nz;
        Lisl_max_util(s_idx,p_idx)=maxUtil_global_nz;

        
        %% POWER PROFILEs
        X_hf = repelem(X,1,1,frame);                  % [C x S x Hf]
        [~,s_resp] = max(X_hf,[],2); s_resp=squeeze(s_resp);
        
        % Compute per-(s,t) power usage with incumbent m and compare to P_T
        % s_star: [C x Hf] responsible sat index
        pow_used = zeros(nSats, Hf);
        
        P_inner = cell(1,nSats);
        D_inner = cell(1,nSats);
        for s = 1:nSats
           Psz = size(P{s});  % [C x M x H]
           Dsz = size(D{s});  % [C x M x H]
           P_inner{s} = repelem(P{s}, 1, 1, frame);  % -> [C x M x Hf]
           D_inner{s} = repelem(D{s}, 1, 1, frame);  % -> [C x M x Hf]
        end
        
        for t = 1:Hf
           for c = 1:number_cells
               ss = s_resp(c,t);
               % served power at (c,t) from incumbent m
               pow_used(ss,t) = pow_used(ss,t) + squeeze(P_inner{ss}(c,:,t)) * squeeze(m(c,:,t)).';
           end
        end
        tight = pow_used ./ P_T(:);
        tight(~isfinite(tight)) = 0;
        % fprintf('Power util: mean=%.3f, p95=%.3f, max=%.3f\n', mean(tight(:)), prctile(tight(:),95), max(tight(:)));
        
        % figure
        % plot(pow_used(51,:))
        % hold on
        % % legend('D','A')
        
        %% ILLUMINATION
        % Count beams used per (s,t)
        beams_used = zeros(nSats, Hf);
        for t=1:Hf
           for c=1:number_cells
               ss = s_resp(c,t);
               beams_used(ss,t) = beams_used(ss,t) + (i(c,t) > 0.5);
           end
        end
        beam_util = beams_used / beams;
        % fprintf('Beams util: mean=%.3f, p95=%.3f, max=%.3f\n', mean(beam_util(:)), prctile(beam_util(:),95), max(beam_util(:)));
        
        % figure
        % plot(beams_used(51,:))
        % hold on
        % legend('D','A')
        
        %% SERVED ILL AND BW -> EFFECTIVE SE
        
        % Served s(c,t) = sum_m D_s*(c,m,t) * m(c,m,t)
        served = zeros(C,t1);
        for t=1:t1
           ss = s_resp(:,t);
           for c=1:C
               served(c,t) = squeeze(D_inner{ss(c)}(c,:,t)) * squeeze(m(c,:,t)).';
           end
        end
        % Effective SE and served-per-light
        lights = sum(i(:,win),'all');
        served_total = sum(served(:));
        B_used = (sum(B,2) .* i(:,win));              % count BW only when illuminated
        B_used = B_used.*(B_T/N);              % * 5MHz

        SE_eff = served_total / max(sum(B_used,'all'),1e-12);    % bits/s/Hz (units follow D)
        served_per_lit = served_total / max(lights,1);
        
        %% -------------------------------
        %  RESULT SAVING 
        % -------------------------------
        UC_w_ALL(s_idx,p_idx)=100*UC_inst_weighted; %
        UC_g_ALL(s_idx,p_idx)=100*UC_inst_global; %
        EC_all(s_idx,p_idx)=100*EC_global; %
        TTS_all(s_idx,p_idx)=mean(TTS_mean_c)*0.1;
        TTS_act_all(s_idx,p_idx)= 100*mean(TTS_active_rate_c); %
        
        Handover_all(s_idx,p_idx)=HO_mean;
        
        Pused_all(s_idx,p_idx)=sum(pow_used,'all')./1e6; % W -> MW (10^6)
        Bused_all(s_idx,p_idx)=sum(B_used,'all')./1e6; % MHz (10*6) -> THz (10^12)

        SE_eff_all(s_idx,p_idx)=SE_eff;
        Served_all(s_idx,p_idx)=served_total./1e6; % Mb (10^6) -> Tb (10^12)


    end

end
% 
% % TRAFFIC (UC/EC) + TIME TOGETHER:
% figure
% hold on
% yyaxis left
% set(gca, 'YColor', '#ea633e','FontSize',12,'YLim',[0,100],'FontName','Computer Modern');
% plot(P_T_values,(UC_w_ALL(1,:)),'Color', '#d7191c', 'LineStyle', '--','LineWidth',1.5,'Marker', 'none');
% plot(P_T_values,(UC_w_ALL(2,:)),'Color','#d7191c','LineStyle',':','LineWidth',1.5,'Marker', 'none')
% plot(P_T_values,(UC_w_ALL(3,:)),'Color','#d7191c','LineStyle','-.','LineWidth',1.5,'Marker', 'none')
% %plot(users,100*(UC_A_array./RC_array),'Color','#d7191c','LineStyle','-','LineWidth',1.5,'Marker', 'none')
% 
% plot(P_T_values,(EC_all(1,:)),'Color','#fdae61','LineStyle','--','LineWidth',1.5,'Marker', 'none')
% plot(P_T_values,(EC_all(2,:)),'Color','#fdae61','LineStyle',':','LineWidth',1.5,'Marker', 'none')
% plot(P_T_values,(EC_all(3,:)),'Color','#fdae61','LineStyle','-.','LineWidth',1.5,'Marker', 'none')
% %plot(users,100*(EC_A_array./RC_array),'Color','#fdae61','LineStyle','-','LineWidth',1.5,'Marker', 'none')
% 
% ylabel ('Traffic [%]','FontName','Computer Modern','FontSize',12)
% yyaxis right
% set(gca, 'YColor', '#2b83ba','FontSize',12,'FontName','Computer Modern'); % Set left y-axis color to blue
% plot(P_T_values,(TTS_all(1,:)),'Color','#2b83ba','LineStyle','--','LineWidth',1.5,'Marker', 'none')
% plot(P_T_values,(TTS_all(2,:)),'Color','#2b83ba','LineStyle',':','LineWidth',1.5,'Marker', 'none')
% plot(P_T_values,(TTS_all(3,:)),'Color','#2b83ba','LineStyle','-.','LineWidth',1.5,'Marker', 'none')
% %plot(users,TTS_A_array,'Color','#2b83ba','LineStyle','-','LineWidth',1.5,'Marker', 'none')
% % title('DB vs GA vs MILP vs MILP_{split} Results')
% xlabel('Power [W]','FontSize',12,'FontName','Computer Modern')
% ylabel ('Time [s]','FontName','Computer Modern','FontSize',12)
% % legend({'UC', 'EC', 'TTS'})
% legend({'UC A', 'UC B','UC C','EC A', 'EC B','EC C','TTS A', 'TTS B','TTS C'})
% grid on
% box on
% hold off
% 
% xlim([250 2500])


% TRAFFIC (UC/UC weighted) + TIME TOGETHER:
figure
hold on
yyaxis left
set(gca, 'YColor', '#ea633e','FontSize',12,'YLim',[0,100],'FontName','Computer Modern');
plot(P_T_values,(UC_w_ALL(1,:)),'Color', '#d7191c', 'LineStyle', '--','LineWidth',1.5,'Marker', 'none');
plot(P_T_values,(UC_w_ALL(2,:)),'Color','#d7191c','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(UC_w_ALL(3,:)),'Color','#d7191c','LineStyle','-.','LineWidth',1.5,'Marker', 'none')

plot(P_T_values,(UC_w_ALL(4,:)),'Color','#fdae61','LineStyle','--','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(UC_w_ALL(5,:)),'Color','#fdae61','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(UC_w_ALL(6,:)),'Color','#fdae61','LineStyle','-.','LineWidth',1.5,'Marker', 'none')

ylabel ('Traffic [%]','FontName','Computer Modern','FontSize',12)
yyaxis right
set(gca, 'YColor', '#2b83ba','FontSize',12,'FontName','Computer Modern'); % Set left y-axis color to blue
plot(P_T_values,(TTS_all(1,:)),'Color','#2b83ba','LineStyle','--','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(TTS_all(2,:)),'Color','#2b83ba','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(TTS_all(3,:)),'Color','#2b83ba','LineStyle','-.','LineWidth',1.5,'Marker', 'none')

plot(P_T_values,(TTS_all(4,:)),'Color','#abd9e9','LineStyle','--','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(TTS_all(5,:)),'Color','#abd9e9','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(TTS_all(6,:)),'Color','#abd9e9','LineStyle','-.','LineWidth',1.5,'Marker', 'none')

xlabel('Power [W]','FontSize',12,'FontName','Computer Modern')
ylabel ('Time [s]','FontName','Computer Modern','FontSize',12)
% legend({'UC', 'EC', 'TTS'})
legend({'UC A', 'UC B','UC C','UC D', 'UC E','UC F','TTS A', 'TTS B','TTS C','TTS D', 'TTS E','TTS F'})
grid on
box on
hold off

xlim([250 2500])

% TRAFFIC (UC/UC global) + TIME TOGETHER:
figure
hold on
yyaxis left
set(gca, 'YColor', '#ea633e','FontSize',12,'YLim',[0,100],'FontName','Computer Modern');
plot(P_T_values,(UC_g_ALL(1,:)),'Color', '#d7191c', 'LineStyle', '--','LineWidth',1.5,'Marker', 'none');
plot(P_T_values,(UC_g_ALL(2,:)),'Color','#d7191c','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(UC_g_ALL(3,:)),'Color','#d7191c','LineStyle','-.','LineWidth',1.5,'Marker', 'none')

plot(P_T_values,(UC_g_ALL(4,:)),'Color','#fdae61','LineStyle','--','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(UC_g_ALL(5,:)),'Color','#fdae61','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(UC_g_ALL(6,:)),'Color','#fdae61','LineStyle','-.','LineWidth',1.5,'Marker', 'none')

ylabel ('Traffic [%]','FontName','Computer Modern','FontSize',12)
yyaxis right
set(gca, 'YColor', '#2b83ba','FontSize',12,'FontName','Computer Modern'); % Set left y-axis color to blue
plot(P_T_values,(TTS_all(1,:)),'Color','#2b83ba','LineStyle','--','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(TTS_all(2,:)),'Color','#2b83ba','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(TTS_all(3,:)),'Color','#2b83ba','LineStyle','-.','LineWidth',1.5,'Marker', 'none')

plot(P_T_values,(TTS_all(4,:)),'Color','#abd9e9','LineStyle','--','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(TTS_all(5,:)),'Color','#abd9e9','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(TTS_all(6,:)),'Color','#abd9e9','LineStyle','-.','LineWidth',1.5,'Marker', 'none')

xlabel('Power [W]','FontSize',12,'FontName','Computer Modern')
ylabel ('Time [s]','FontName','Computer Modern','FontSize',12)
% legend({'UC', 'EC', 'TTS'})
legend({'UC A', 'UC B','UC C','UC D', 'UC E','UC F','TTS A', 'TTS B','TTS C','TTS D', 'TTS E','TTS F'})
grid on
box on
hold off

xlim([250 2500])



% SERVED TRAFFIC:
figure
hold on
plot(P_T_values,(Served_all(1,:)),'Color','#d7191c', 'LineStyle','--','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(Served_all(2,:)),'Color','#d7191c','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(Served_all(3,:)),'Color','#d7191c','LineStyle','-.','LineWidth',1.5,'Marker', 'none')

plot(P_T_values,(Served_all(4,:)),'Color','#fdae61', 'LineStyle','--','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(Served_all(5,:)),'Color','#fdae61','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(Served_all(6,:)),'Color','#fdae61','LineStyle','-.','LineWidth',1.5,'Marker', 'none')

hold off
%title('DB vs GA vs MILP vs MILP_{split} Exectution Time')
xlabel('Power [W]','FontSize',12,'FontName','Computer Modern')
ylabel ('Total Served Traffic [Tb]','FontName','Computer Modern','FontSize',12)
legend({'A', 'B', 'C','D', 'E', 'F'})
set(gca, 'FontSize',12,'FontName','Computer Modern');
grid on
box on
hold off

xlim([250 2500])

% P USED:
figure
hold on
plot(P_T_values,(Pused_all(1,:)),'Color','#d7191c','LineStyle','--','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(Pused_all(2,:)),'Color','#d7191c','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(Pused_all(3,:)),'Color','#d7191c','LineStyle','-.','LineWidth',1.5,'Marker', 'none')

plot(P_T_values,(Pused_all(4,:)),'Color','#fdae61','LineStyle','--','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(Pused_all(5,:)),'Color','#fdae61','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(Pused_all(6,:)),'Color','#fdae61','LineStyle','-.','LineWidth',1.5,'Marker', 'none')

hold off
%title('DB vs GA vs MILP vs MILP_{split} Exectution Time')
xlabel('Power [W]','FontSize',12,'FontName','Computer Modern')
ylabel ('Total Power Used [MW]','FontName','Computer Modern','FontSize',12)
legend({'A', 'B', 'C','D', 'E', 'F'})
set(gca, 'FontSize',12,'FontName','Computer Modern');
grid on
box on
hold off

xlim([250 2500])

% BW USED:
figure
hold on
plot(P_T_values,(Bused_all(1,:)),'Color','#d7191c', 'LineStyle','--','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(Bused_all(2,:)),'Color','#d7191c','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(Bused_all(3,:)),'Color','#d7191c','LineStyle','-.','LineWidth',1.5,'Marker', 'none')

plot(P_T_values,(Bused_all(4,:)),'Color','#fdae61', 'LineStyle','--','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(Bused_all(5,:)),'Color','#fdae61','LineStyle',':','LineWidth',1.5,'Marker', 'none')
plot(P_T_values,(Bused_all(6,:)),'Color','#fdae61','LineStyle','-.','LineWidth',1.5,'Marker', 'none')
hold off
%title('DB vs GA vs MILP vs MILP_{split} Exectution Time')
xlabel('Power [W]','FontSize',12,'FontName','Computer Modern')
ylabel ('Total Bandwidth Used [THz]','FontName','Computer Modern','FontSize',12)
legend({'A', 'B', 'C','D', 'E', 'F'})
set(gca, 'FontSize',12,'FontName','Computer Modern');
grid on
box on
hold off

xlim([250 2500])



