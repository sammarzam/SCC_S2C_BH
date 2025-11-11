function [D,P,R_D,c_scenario,number_cells,theta,M]=pre_BH_computation(freq,h_sat,Re,frame,frame_dur,B,B_T,N,Dcell,h3_resolution,position,sat_idx,d,H) % Precompute scenery based matrixes.

global PWD;

SE=[0.1523, 0.3770, 0.8770, 1.4766, 1.9141, 2.4063, 2.7305, 3.3223, 3.9023, 4.5234, 5.1152, 5.5547, 6.2266, 6.9141, 7.4063]; % MCS index table 2 for PDSCH
M=length(SE);

for CQI=1:length(SE)
    C_N(CQI)=(CQI-2.9142)/0.2789;
    if C_N(CQI)>=4
        C_N(CQI)=(CQI-1.8749)/0.5246;
        if C_N(CQI)>=14
            C_N(CQI)=(CQI-1.4946)/0.5369;
        end
    end
end

%  PAINT
% C/N vs SE plot:
% figure
% plot(C_N, SE,'bx')
% hold on
% stairs(C_N, SE)
% xlabel('C/N [dB]')
% ylabel('SE [bit/s/Hz]')
% title('MCS PDSCH 5G')
% hold off

% Edge lengths [km] H3 resolutions from 0 -> 15
Edge_lengths_H3 = [ ...
    1281.256011, 483.0568391, 182.5129565, 68.97922179, ...
    26.07175968, 9.854090990, 3.724532667, 1.406475763, ...
    0.531414010, 0.200786148, 0.075863783, 0.028663897, ...
    0.010830188, 0.004092010, 0.001546100, 0.000584169 ];

%Cell theta calculation as a function of resolution: en funcion de resolucion
%1)
gamma=(2*Edge_lengths_H3(h3_resolution + 1)*360)/(2*pi*Re*2);
%2)
h_prime=Re*(1-cos(gamma*pi/180));
z=sin(gamma*pi/180)*Re;
%3)
betta=(atan(z/(h_prime+h_sat)))*180/pi;
%4)
theta=2*betta;

% UE assumption:
antenna_diameter=0.49;  %Starlink
g_rx=10*log10(0.65*((pi*antenna_diameter/((3*10^8)/freq)))^2); %dB
lnb=2; %dB
T_noise_rx=((10^(lnb/10)-1)*290); %K

% Cell Footprint 
c_scenario=[];
number_cells=size(Dcell,1);
for i=1:number_cells
    c_scenario=[c_scenario c(i,[Dcell.Latitude(i), Dcell.Longitude(i)],Edge_lengths_H3(h3_resolution + 1),[])]; %constructor
    lla = ecef2lla(position(:,:,sat_idx)'); 
    lat_sat = lla(:,1);
    lon_sat = lla(:,2);
    c_scenario(i).compute_betta_to_sat(h_sat, lat_sat, lon_sat); % Compute betta to the cell center so that the gain loss due to beam scanning (moving ftom boresight) can be then accounted: non-ideal isotropic behavior of the embedded element gain.
    c_scenario(i).adduser(u([Dcell.Latitude(i), Dcell.Longitude(i)],d(i,:),g_rx,T_noise_rx));  %constructor
    c_scenario(i).users(length(c_scenario(i).users)).compute_distance_elevation_betta_to_sat(h_sat, lat_sat, lon_sat); % Compute distance, elevation and betta angle to the satellite of the the last added user
    c_scenario(i).users(length(c_scenario(i).users)).compute_betta_to_cell_center(h_sat,c_scenario(i)); % Compute the betta angle with respect to the cell center so that the gain loss can be accounted by the fact of not being at the cell center    
    % % PAINT
   % c_scenario(i).draw(1);
   % hold on
end
% % PAINT
% title(strcat('Cell Scenario:  ',num2str(rings),' rings'))


%3)
% Scenery based MATRIX generation initialization:
P=zeros(number_cells,M,H);
D=zeros(number_cells,M,H);
R_D=zeros(number_cells,H*frame);

% Link Budget Caculation for each M (MODCOD and therefore C/N->SE), N (numer of assigned BAND) and u (user, depending on the scenario):
% The Demand matrix (D) does not change, depending on the user. It is the
% power the one varies from user to user depending on the particular link
% conditions. 

% Input Kink Budget data:
% CANCELAR EXTRA LOSSES:
Extra_losses.At=0;
G_t=10*log10(0.65*48360/theta^2); % It needs to be corrected by user's position.
k=-228.6012; %dBW/(HzK^-1)
T_ant=30; %K
alpha=0.1;

% Antenna pattern
q=1.3;
np=1001;
theta_scanning=linspace(-pi/2,pi/2,np);

for t=1:H
    for m=1:M
        for i=1:number_cells
            D(i,m,t)=SE(m)*((sum(B(i,:))*(B_T/N))/(1+alpha))*frame_dur; % [Mbits/s]*frame_dur
                for k_u=1:length(c_scenario(i).users) % 1! 
                    %R(c_scenario(i).users(k_u).id,:)=c_scenario(i).users(k_u).traffic_demand*frame_dur; % [Mbits/s]*frame_dur
                    % C_N=(G_tx+P_tx-Lfsl-Extra_losses.At+G_rx)-(k+T_eq+BW_rx)->Determine:P_tx
                     % Gain TX:
                     G_tx=G_t-12*(c_scenario(i).users(k_u).betta_to_center/theta)^2; % Gain loss approximation when moving away from the cell center. <- ALWAYS 0
                     [minScanning,closestScanningIndex]=min(abs(theta_scanning(1:((np-1)/2+1))+deg2rad(c_scenario(i).betta_to_sat(t)))); % Gain loss due to scanning
                     G_tx=G_tx+20*log10(cos(abs(theta_scanning(closestScanningIndex)))^q);
                     % Gain RX:
                     G_rx=c_scenario(i).users(k_u).gain; %dB           
                     % T_eq RX:
                     T_eq=10*log10(T_ant+c_scenario(i).users(k_u).T_noise+280*(1-10^(-Extra_losses.At/10))); %[dBK] -> Contributors: captured by antenna + rx equivalent + attenuation due to extra losses
                     Lfsl=20*log10(freq)+20*log10(c_scenario(i).users(k_u).distace_to_sat(t)*10^3)-147.55;
                     BW_rx=10*log10((sum(B(i,:)))*(B_T/N)*10^6); 
                     P(i,m,t)=10^((C_N(m)-G_tx+Lfsl+Extra_losses.At-G_rx+(k+T_eq+BW_rx))/10); % [W]
                end
        end
    end
end

% 
for t_sat=1:H
    for t_bh=1:frame %1frame=10 subframe by 10subframe (beamswithching is performed each frame, for 10 slots constant illumination constant)
        % Traffic Pending Adjustment based on frame length: from second to second add the total demand/s to the pending:
        if mod(t_bh-1,1/frame_dur)==0
            for i=1:number_cells
                for j=1:length(c_scenario(i).users)
                    R_D(i,t_bh+(t_sat-1)*frame)=c_scenario(i).users(j).traffic_demand(t_sat);
                end
            end
        end
    end
end


end
