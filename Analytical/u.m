classdef u < handle
   properties (Access=public)
      % Static
      location
      type_station
      traffic_demand
      gain
      T_noise
      % Computed
      betta_to_sat
      distace_to_sat
      elevation_to_sat
      betta_to_center
      % Dynamic (t=1:frame)
      C_N=[]
      P=[]
      BW=[]
      I
      traffic_served=[]
      traffic_pending=0
   end
   methods
       function obj = u(location,traffic_demand,gain,T_noise) %constructor
         obj.location = location;
         obj.traffic_demand = traffic_demand;
         obj.gain=gain;
         obj.T_noise=T_noise;
       end
       % function obj = compute_distance_elevation_betta_to_sat(obj, h)  % Distance and Elevation to Satellite (0,0):
       %     Re=6378; % Earth Radius [km]
       %     %1)
       %     gamma=(sqrt(obj.location(1)^2+obj.location(2)^2)*360)/(2*pi*Re);
       %     %2)
       %     h_prime=Re*(1-cos(gamma*pi/180));
       %     z=sin(gamma*pi/180)*Re;
       %     %3)
       %     betta=(atan(z/(h_prime+h)))*180/pi;
       %     obj.betta_to_sat=betta;
       %     %4)
       %     obj.distace_to_sat=(h+h_prime)/cos(betta*pi/180);
       %     obj.elevation_to_sat=90-betta-gamma;
       % end
       function obj = compute_distance_elevation_betta_to_sat(obj, h, lat_sp, lon_sp)
            Re = 6378;  % Radio de la Tierra [km]
        
            % Punto de observación (ej: receptor o celda)
            lat_user = obj.location(1);
            lon_user = obj.location(2);
        
            % Convertir ambos puntos a coordenadas ECEF
            ecef_user = lla2ecef([lat_user, lon_user, 0]);
            
            for t=1:length(lat_sp)

                ecef_sat  = lla2ecef([lat_sp(t), lon_sp(t), 0]);  % subsatélite sobre la Tierra
            
                % Calcular ángulo central γ (rad) entre los dos vectores desde el centro de la Tierra
                u = ecef_user / norm(ecef_user);
                v = ecef_sat  / norm(ecef_sat);
                gamma_rad = acos(dot(u, v));  % en radianes
            
                % Calcular componentes geométricos
                h_prime = Re * (1 - cos(gamma_rad));
                z = Re * sin(gamma_rad);
                betta = atan(z / (h_prime + h)) * 180 / pi;  % en grados
            
                % Calcular distancia satelital y ángulo de elevación
                distance = (h + h_prime) / cosd(betta);  % distancia slant
                elevation = 90 - betta;  % en grados
            
                % Guardar en el objeto
                obj.betta_to_sat = [obj.betta_to_sat, betta];
                obj.distace_to_sat = [obj.distace_to_sat, distance];
                obj.elevation_to_sat =[obj.elevation_to_sat, elevation];
            end
       end


       function obj = compute_betta_to_cell_center(obj, h, c)  % Betta computation as if center of the cell would have been the user, then compare the difference with respect to the real user and evaluate the real apperture from cell center.
           Re=6378; % Earth Radius [km]
           %1)
           gamma=(sqrt((obj.location(1)-c.center(1))^2+(obj.location(2)-c.center(2))^2)*360)/(2*pi*Re);
           %2)
           h_prime=Re*(1-cos(gamma*pi/180));
           z=sin(gamma*pi/180)*Re;
           %3)
           betta=(atan(z/(h_prime+h)))*180/pi;
           %4)
           obj.betta_to_center=betta; %abs(obj.betta_to_sat-betta)
       end
       function draw(obj)
       plot(obj.location(1),obj.location(2),'-s','MarkerSize',7,'MarkerEdgeColor','black','MarkerFaceColor','black') 
       end
   end
end

