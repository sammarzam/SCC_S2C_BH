function BH_main_fixed_normalization(scenario,use_case,beams,h3_resolution,r0,rmax,d_threshold,Hcap,P_T,m_continuous,MIPGap)

    global PWD;
    PWD=pwd;

    inputDir = fullfile(PWD, 'Input_Data');
    outputDir = fullfile(PWD, 'Output_Data_Full');

    normalizationDir = fullfile(PWD, 'Normalization_Data'); %Simulation_Results_A_framex10

    if strcmp(scenario,'A')
        %% ========== IMPORT A)==========
        % --- Scenario ---
        % use_case='iridium'; % 'starlink'; %'iris2';
    
        % h3_resolution=2;
        % r0 = 2;                 % base res (start here)
        % rmax = 2;               % max res
        % d_threshold = 5e3;      % per-cell demand cap (tune to your units)
    
        % --- Inputs ---
        % Cell set:
        addpath(fullfile(inputDir,'COMMS_JUAN','scripts'))
        load(strcat('EUR_Dcelda_h3_',num2str(h3_resolution),'.mat')) % Dcell
    
        % Cell coloring: B
        addpath(fullfile(inputDir,'CELL_COLOURING_MAX_BAND_ALLOCATION_REPAIR_FORWARD'))
        tag = sprintf('r0%d_rmax%d_thr%g', r0, rmax, d_threshold); % y, Btot, K
        load(['EUR_B_h3_adaptive_' tag '_V5.mat'])
        B=y;
        B_T=Btot;
        N=K;
        
        % UNIFORM Band distribution
        B=B.*0;
        reuse_factor=3;
    
        band_bins_uniform=floor((B_T/reuse_factor)/(B_T/N));
        % band_bins_uniform=12
        for i=1:size(B,1)
            B(i,1:band_bins_uniform)=1;
        end
        clear K
        clear Btot
    
        % Visibility (K_i,t) and routing R(s,k,t), h_sat, position,
        addpath(fullfile(inputDir,'SAT_ROUTING'))
        load(strcat('EUR_R_K_res_',num2str(h3_resolution),'_',use_case,'.mat')) % cell(nV, H), each K{i,t} is vector of SAT IDs visible to cell i at time t
        H=size(K,2);
        nSats=size(R,1);
    
        % Demand d, sol, beams
        addpath(fullfile(inputDir,'SAT_CELL_ASSOCIATION'))
        load(strcat('EUR_X_res_',num2str(h3_resolution),'_',use_case,'_beams_',num2str(beams),'_closest.mat')) % sol.x, associate satellite s to cell i at time t
        X=sol.x;



    elseif strcmp(scenario,'B')
    
        %% ========== IMPORT B)==========
        % --- Scenario ---
        % use_case='iridium'; % 'starlink'; %'iris2';
    
        % h3_resolution=2;
        % r0 = 2;                 % base res (start here)
        % rmax = 2;               % max res
        % d_threshold = 5e3;      % per-cell demand cap (tune to your units)
    
        % --- Inputs ---
        % Cell set:
        addpath(fullfile(inputDir,'COMMS_JUAN','scripts'))
        load(strcat('EUR_Dcelda_h3_',num2str(h3_resolution),'.mat')) % Dcell
    
        % Cell coloring: B
        addpath(fullfile(inputDir,'CELL_COLOURING_MAX_BAND_ALLOCATION_REPAIR_FORWARD'))
        tag = sprintf('r0%d_rmax%d_thr%g', r0, rmax, d_threshold); % y, Btot, K
        load(['EUR_B_h3_adaptive_' tag '_V5.mat'])
        B=y;
        B_T=Btot;
        N=K;
        % UNIFORM Band distribution
        B=B.*0;
        reuse_factor=3;
    
        band_bins_uniform=floor((B_T/reuse_factor)/(B_T/N));
        for i=1:size(B,1)
            B(i,1:band_bins_uniform)=1;
        end
        clear K
        clear Btot
    
        % Visibility (K_i,t) and routing R(s,k,t), h_sat, position,
        addpath(fullfile(inputDir,'SAT_ROUTING'))
        load(strcat('EUR_R_K_res_',num2str(h3_resolution),'_',use_case,'.mat')) % cell(nV, H), each K{i,t} is vector of SAT IDs visible to cell i at time t
        H=size(K,2);
        nSats=size(R,1);
    
        % Demand d, sol, beams
        addpath(fullfile(inputDir,'SAT_CELL_ASSOCIATION'))
        load(strcat('EUR_X_res_',num2str(h3_resolution),'_',use_case,'_beams_',num2str(beams),'_d_only.mat')) % sol.x, associate satellite s to cell i at time t
        X=sol.x;

    elseif strcmp(scenario,'C')    
        %% ========== IMPORT C)==========
        % --- Scenario ---
        % use_case='iridium'; % 'starlink'; %'iris2';
    
        % h3_resolution=2;
        % r0 = 2;                 % base res (start here)
        % rmax = 2;               % max res
        % d_threshold = 5e3;      % per-cell demand cap (tune to your units)
    
        % --- Inputs ---
        % Cell set:
        addpath(fullfile(inputDir,'COMMS_JUAN','scripts'))
        load(strcat('EUR_Dcelda_h3_',num2str(h3_resolution),'.mat')) % Dcell
    
        % Cell coloring: B
        addpath(fullfile(inputDir,'CELL_COLOURING_MAX_BAND_ALLOCATION_REPAIR_FORWARD'))
        tag = sprintf('r0%d_rmax%d_thr%g', r0, rmax, d_threshold); % y, Btot, K
        load(['EUR_B_h3_adaptive_' tag '_V5.mat'])
        B=y;
        B_T=Btot;
        N=K;
        % UNIFORM Band distribution
        B=B.*0;
        reuse_factor=3;
    
        band_bins_uniform=floor((B_T/reuse_factor)/(B_T/N));
        for i=1:size(B,1)
            B(i,1:band_bins_uniform)=1;
        end
        clear K
        clear Btot
    
        % Visibility (K_i,t) and routing R(s,k,t), h_sat, position,
        addpath(fullfile(inputDir,'SAT_ROUTING'))
        load(strcat('EUR_R_K_res_',num2str(h3_resolution),'_',use_case,'.mat')) % cell(nV, H), each K{i,t} is vector of SAT IDs visible to cell i at time t
        H=size(K,2);
        nSats=size(R,1);
    
        % Demand d, sol, beams
        addpath(fullfile(inputDir,'SAT_CELL_ASSOCIATION'))
        load(strcat('EUR_X_res_',num2str(h3_resolution),'_',use_case,'_beams_',num2str(beams),'_d.mat')) % sol.x, associate satellite s to cell i at time t
        X=sol.x;
   
    elseif strcmp(scenario,'D')
        %% ========== IMPORT D)==========
        %--- Scenario ---
        % use_case='iridium'; % 'starlink'; %'iris2';
    
        % h3_resolution=2;
        % r0 = 2;                 % base res (start here)
        % rmax = 2;               % max res
        % d_threshold = 5e3;      % per-cell demand cap (tune to your units)
    
        %--- Inputs ---
        %Cell set:
        addpath(fullfile(inputDir,'COMMS_JUAN','scripts'))
        load(strcat('EUR_Dcelda_h3_',num2str(h3_resolution),'.mat')) % Dcell
    
        %Cell coloring: B
        addpath(fullfile(inputDir,'CELL_COLOURING_MAX_BAND_ALLOCATION_REPAIR_FORWARD'))
        tag = sprintf('r0%d_rmax%d_thr%g', r0, rmax, d_threshold); % y, Btot, K
        load(['EUR_B_h3_adaptive_' tag '_V5.mat'])
        B=y;
        B_T=Btot;
        N=K;
        clear K
        clear Btot
    
        %Visibility (K_i,t) and routing R(s,k,t), h_sat, position,
        addpath(fullfile(inputDir,'SAT_ROUTING'))
        load(strcat('EUR_R_K_res_',num2str(h3_resolution),'_',use_case,'.mat')) % cell(nV, H), each K{i,t} is vector of SAT IDs visible to cell i at time t
        H=size(K,2);
        nSats=size(R,1);
    
        %Demand d, sol, beams
        addpath(fullfile(inputDir,'SAT_CELL_ASSOCIATION'))
        load(strcat('EUR_X_res_',num2str(h3_resolution),'_',use_case,'_beams_',num2str(beams),'_closest.mat')) % sol.x, associate satellite s to cell i at time t
        X=sol.x;
    
        % ASSIGN FULL BAND TO NON ILLUMINATED CELLS TO FORCE THE DEACTIVATION
        no_demand_c=find(sum(d,2)==0);
        assigned_with_demand_0=find(sum(B(no_demand_c,:),2)>0);
        assigned_with_demand_0_cells=no_demand_c(assigned_with_demand_0);
    
        B(no_demand_c,:)=1; % CHECK! FORCE NOT TO ILLUMINATE NON DEMANDED CELLS!!
    
        % d=d./2;

    elseif strcmp(scenario,'E')
        %% ========== IMPORT E)==========
        % --- Scenario ---
        % use_case='iridium'; % 'starlink'; %'iris2';
        
        % h3_resolution=2;
        % r0 = 2;                 % base res (start here)
        % rmax = 2;               % max res
        % d_threshold = 5e3;      % per-cell demand cap (tune to your units)
        
        % --- Inputs ---
        % Cell set:
        addpath(fullfile(inputDir,'COMMS_JUAN','scripts'))
        load(strcat('EUR_Dcelda_h3_',num2str(h3_resolution),'.mat')) % Dcell
        
        % Cell coloring: B
        addpath(fullfile(inputDir,'CELL_COLOURING_MAX_BAND_ALLOCATION_REPAIR_FORWARD'))
        tag = sprintf('r0%d_rmax%d_thr%g', r0, rmax, d_threshold); % y, Btot, K
        load(['EUR_B_h3_adaptive_' tag '_V5.mat'])
        B=y;
        B_T=Btot;
        N=K;
        clear K
        clear Btot
        
        % Visibility (K_i,t) and routing R(s,k,t), h_sat, position,
        addpath(fullfile(inputDir,'SAT_ROUTING'))
        load(strcat('EUR_R_K_res_',num2str(h3_resolution),'_',use_case,'.mat')) % cell(nV, H), each K{i,t} is vector of SAT IDs visible to cell i at time t
        H=size(K,2);
        nSats=size(R,1);
        
        % Demand d, sol, beams
        addpath(fullfile(inputDir,'SAT_CELL_ASSOCIATION'))
        load(strcat('EUR_X_res_',num2str(h3_resolution),'_',use_case,'_beams_',num2str(beams),'_d_only.mat')) % sol.x, associate satellite s to cell i at time t
        X=sol.x;
        
        % ASSIGN FULL BAND TO NON ILLUMINATED CELLS TO FORCE THE DEACTIVATION
        no_demand_c=find(sum(d,2)==0);
        assigned_with_demand_0=find(sum(B(no_demand_c,:),2)>0);
        assigned_with_demand_0_cells=no_demand_c(assigned_with_demand_0);

        B(no_demand_c,:)=1; % CHECK! FORCE NOT TO ILLUMINATE NON DEMANDED CELLS!!

    elseif strcmp(scenario,'F')
        %% ========== IMPORT F)==========
        % --- Scenario ---
        % use_case='iridium'; % 'starlink'; %'iris2';
        
        % h3_resolution=2;
        % r0 = 2;                 % base res (start here)
        % rmax = 2;               % max res
        % d_threshold = 5e3;      % per-cell demand cap (tune to your units)
        
        % --- Inputs ---
        % Cell set:
        addpath(fullfile(inputDir,'COMMS_JUAN','scripts'))
        load(strcat('EUR_Dcelda_h3_',num2str(h3_resolution),'.mat')) % Dcell
        
        % Cell coloring: B
        addpath(fullfile(inputDir,'CELL_COLOURING_MAX_BAND_ALLOCATION_REPAIR_FORWARD'))
        tag = sprintf('r0%d_rmax%d_thr%g', r0, rmax, d_threshold); % y, Btot, K
        load(['EUR_B_h3_adaptive_' tag '_V5.mat'])
        B=y;
        B_T=Btot;
        N=K;
        clear K
        clear Btot
        
        % Visibility (K_i,t) and routing R(s,k,t), h_sat, position,
        addpath(fullfile(inputDir,'SAT_ROUTING'))
        load(strcat('EUR_R_K_res_',num2str(h3_resolution),'_',use_case,'.mat')) % cell(nV, H), each K{i,t} is vector of SAT IDs visible to cell i at time t
        H=size(K,2);
        nSats=size(R,1);
        
        % Demand d, sol, beams
        addpath(fullfile(inputDir,'SAT_CELL_ASSOCIATION'))
        load(strcat('EUR_X_res_',num2str(h3_resolution),'_',use_case,'_beams_',num2str(beams),'_d.mat')) % sol.x, associate satellite s to cell i at time t
        X=sol.x;
        
        % ASSIGN FULL BAND TO NON ILLUMINATED CELLS TO FORCE THE DEACTIVATION
        no_demand_c=find(sum(d,2)==0);
        assigned_with_demand_0=find(sum(B(no_demand_c,:),2)>0);
        assigned_with_demand_0_cells=no_demand_c(assigned_with_demand_0);
    
        B(no_demand_c,:)=1; % CHECK! FORCE NOT TO ILLUMINATE NON DEMANDED CELLS!!

    end

    % %% ========== IMPORT F_BW)==========
    % % --- Scenario ---
    % use_case='iridium'; % 'starlink'; %'iris2';
    % 
    % h3_resolution=2;
    % r0 = 2;                 % base res (start here)
    % rmax = 2;               % max res
    % d_threshold = 5e3;      % per-cell demand cap (tune to your units)
    % 
    % % --- Inputs ---
    % % Cell set:
    % addpath([PWD,'\..\COMMS_JUAN\scripts'])
    % load(strcat('EUR_Dcelda_h3_',num2str(h3_resolution),'.mat')) % Dcell
    % 
    % % Cell coloring: B
    % addpath([PWD,'\..\CELL_COLOURING_MAX_BAND_ALLOCATION_REPAIR_FORWARD'])
    % tag = sprintf('r0%d_rmax%d_thr%g', r0, rmax, d_threshold); % y, Btot, K
    % load(['EUR_B_h3_adaptive_' tag '_V4.mat'])
    % B=y;
    % B_T=Btot;
    % N=K;
    % clear K
    % clear Btot
    % 
    % % Visibility (K_i,t) and routing R(s,k,t), h_sat, position,
    % addpath([PWD,'\..\SAT_ROUTING'])
    % load(strcat('EUR_R_K_res_',num2str(h3_resolution),'_',use_case,'.mat')) % cell(nV, H), each K{i,t} is vector of SAT IDs visible to cell i at time t
    % H=size(K,2);
    % nSats=size(R,1);
    % 
    % % Demand d, sol, beams
    % addpath([PWD,'\..\SAT_CELL_ASSOCIATION'])
    % load(strcat('EUR_X_res_',num2str(h3_resolution),'_',use_case,'_beams_',num2str(beams),'_bw.mat')) % sol.x, associate satellite s to cell i at time t
    % X=sol.x;
    % 
    % Result_Folder='Simulation_Results_F_V4_framex10_d0_mC_normA';
    % Normalization_Ref_Folder='Simulation_Results_A_framex10';
    
    % 
    % %-----------
    % eta=[0.1523, 0.3770, 0.8770, 1.4766, 1.9141, 2.4063, 2.7305, 3.3223, 3.9023, 4.5234, 5.1152, 5.5547, 6.2266, 6.9141, 7.4063];
    % figure
    % plot(sum(max(d(:,1)-eta.*sum(B,2)*(B_T/N),0)),'LineWidth',2)
    % hold on
    % plot(sum(max(d(:,1)-eta*(B_T/3),0)),'LineWidth',2)
    % legend('BW^*','BW^r')
    % ylabel('UC')
    % xlabel('MCS')
    % 
    % figure
    % plot((-1)*sum(min(d(:,1)-eta.*sum(B,2)*(B_T/N),0)),'LineWidth',2)
    % hold on
    % plot((-1)*sum(min(d(:,1)-eta.*(B_T/3),0)),'LineWidth',2)
    % legend('BW^*','BW^r')
    % ylabel('EC')
    % xlabel('MCS')
    
    % disp('Generic overall UC - BW Method:')
    % sum(max(d(:,1)-eta*sum(B,2)*(B_T/N),0))
    % disp('Generic overall EC - BW Method:')
    % (-1)*sum(min(d(:,1)-eta*sum(B,2)*(B_T/N),0))
    % 
    % disp('Generic overall UC - BW Re-use3:')
    % sum(max(d(:,1)-eta*(B_T/3),0))
    % disp('Generic overall EC - BW Re-use3:')
    
    
    
    %% COMPUTE BH WITHIN 𐤃t FOR EACH SATELLITE VISIBLE IN REGION (K{i,t}) FOR A TOTAL DURATION OF H INSTANTS! 
    freq=20e9; % [Hz] Ka Band (DL)
    % h_sat=XXXX; % {km] <- LOADED FROM constellation routing.
    % el_min=25; % Minimum Elevation Angle [º]  <- NOT NEEDED 
    Re=6378; % Earth Radius [km]
    T_T=floor(samplin);
    frame_dur=0.1;%0.1; % [s] -> 100ms/frame * 100 frames = T_T=10s
    frame=T_T/frame_dur; % 
    % colours=1; <- NOT NEEDED
    % beams=16;
    % P_T=2000; % [W]
    % n_users=50; <- NOT NEEDED
    % traffic_model='uniform'; %'linear' %'hotspot' <- NOT NEEDED
    % B_T=250/2; % [MHz/colour] <- DEFINED ABOVE
    % N=B_T/25; % Frequency bins <- DEFINED ABOVE
    % rings=10; % Fixed grid: EFC <- IMPORTED
    % number_cells=0;
    % for i=1:rings
    %     number_cells=number_cells+6*(i-1);
    % end
    % number_cells=number_cells+1;
    
    if Hcap<H
        H=Hcap
        X=X(:,:,1:Hcap);
        R=R(:,:,1:Hcap);
        K=K(:,1:Hcap);
        d=d(:,1:Hcap);
        position=position(:,1:Hcap,:);
    
    end
        
    % SELECTRESOLUTION METHODOLOGY: TIME-SPLIT MILP
    beta=0.7;
    
    disp('Loading MILP...');
    cd (fullfile(PWD,'Analytical'))
    
    P=cell(1,size(R,1));
    D=cell(1,size(R,1));
        
    % Pre-BH:
    for sat_idx=1:size(R,1)
        [D{sat_idx},P{sat_idx},R_D,c_scenario,number_cells,theta,M]=pre_BH_computation(freq,h_sat,Re,frame,frame_dur,B,B_T,N,Dcell,h3_resolution,position,sat_idx,d,H);  
    end
    
    
    % BH: Build global Ω-compact model ONCE.
    % [A,total_range,total_constraints,b, vtype,i_index,m_index,d_index,g_index,f_index,c_index,i_range,m_range,d_range,g_range,f_range,c_range,i_aux,m_aux,d_aux,g_aux,f_aux,c_aux, t_constraint_label]=BH(H,frame,beams,P_T,D,P,R_D,M, number_cells,nSats,X);
    % [A,total_range,total_constraints,b,vtype,i_index,m_index,d_index,g_index,f_index,c_index, i_range,m_range,d_range,g_range,f_range,c_range,i_aux,m_aux,d_aux,g_aux,f_aux,c_aux,t_constraint_label,sense] = BH_vectored_OLD(H,frame,beams,P_T,D,P,R_D,M, number_cells, nSats, X)
    [A,total_range,total_constraints,b,sense,vtype,i_index,m_index,g_index,f_index,c_index, i_range,m_range,g_range,f_range,c_range,i_aux,m_aux,g_aux,f_aux,c_aux,t_row_label,tlabel_i, tlabel_m, tlabel_g, tlabel_f, tlabel_c] = BH_vectored(H,frame,beams,P_T,D,P,R_D,M, number_cells, nSats, X, m_continuous);
    
    % Initialize lb and ub + warm start:
    lb = zeros(total_range,1);       % will be overwritten per window
    ub = inf(total_range,1);         % nice to have a default, will be overwritten per window
    ub(vtype=='B') = 1;              % harmless default
    mip_start = zeros(total_range,1);
    
    % Choose window size W (in inner slots)
    W = 10*1;  
        
    fname = sprintf('%s_pre_BH_[%s_res%d_beams%d]_data.mat', scenario, use_case, h3_resolution,beams);
    
    save(fullfile(outputDir, fname), 'D','P','R_D','W','H','frame','P_T','number_cells','nSats','M','beams');
    
    use_case_normalization="iridium";

    for t_current = 1:W:(H*frame)
            fprintf(['Evaluation of t_current=', num2str(t_current),'\n'])
            for betta=[0,1,beta]
                if betta==0 % NORMALIZATION: max_UC_betta_0
                    % fname = sprintf('%s_BH_[%s_res%d]_P_%d_mC_%d_beta_%0.2f_win_%d_%d.mat', scenario, use_case, h3_resolution, P_T, m_continuous, betta, t_current, t_current+W-1);
                    fname = sprintf('BH_[%s_res%d]_beta_%0.2f_win_%d_%d.mat',use_case_normalization, h3_resolution, betta, t_current, t_current+W-1); % A

                    load(fullfile(normalizationDir, fname),'objectives');
                    max_UC_betta_0=objectives.UC;
                    clear objectives
     
                elseif betta==1 % NORMALIZATION: max_EC_betta_1, max_TTS_betta_1
                    % fname = sprintf('%s_BH_[%s_res%d]_P_%d_mC_%d_beta_%0.2f_win_%d_%d.mat', scenario, use_case, h3_resolution, P_T, m_continuous, betta, t_current, t_current+W-1);
                    fname = sprintf('BH_[%s_res%d]_beta_%0.2f_win_%d_%d.mat',use_case_normalization, h3_resolution, betta, t_current, t_current+W-1); % A

                    load(fullfile(normalizationDir, fname),'objectives');
                    max_EC_betta_1=objectives.EC;
                    max_TTS_betta_1=objectives.TTS;
                    clear objectives
    
                else % FULL MILP EXECUTION
                    %MIPGap=0.01;
                    normUC=max_UC_betta_0;
                    if max_EC_betta_1==0
                        normEC=1;
                    else
                        normEC=max_EC_betta_1;
                    end
                    if max_TTS_betta_1==0
                        normTTS=1;
                    else
                        normTTS=max_TTS_betta_1;
                    end
                    cd (PWD)
                    [res, data, solutions] = gurobi_execution_BH_windowed(betta, normUC, normEC, normTTS, total_range, A, b, sense, vtype, MIPGap, mip_start, lb, ub, i_index, i_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, t_row_label, tlabel_i, tlabel_m, tlabel_g, tlabel_f, tlabel_c, t_current, W, R_D);
                    cd (PWD)
                    xrow = res.x(:)';  % incumbent row (1 × total_range)
                    [~,~,~]= solution_plot_saving(outputDir, scenario, use_case, h3_resolution, beams, P_T, m_continuous, betta, normUC, normEC, normTTS, data, xrow, i_index,i_range,m_index,m_range,g_index,g_range,f_index,f_range,c_index,c_range, tlabel_g, tlabel_f, tlabel_c,t_current, min(t_current + W - 1, H*frame), number_cells, H, frame, M, X, D, nSats);
    
    
                    % Update incumbent (carry to next window)
                    if ~isempty(xrow)
                        mip_start = xrow;   % this fixes all outside-window vars next call
                    else
                        error('No incumbent returned in window starting at t=%d', t_current);
                    end
    
                end  
            end
    end
end
