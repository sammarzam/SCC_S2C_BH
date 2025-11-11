<<<<<<< HEAD
<<<<<<< HEAD
classdef c < handle
   properties (Access=public)
      number
      center
      radius
      users=[]
      % Computed
      betta_to_sat
      interfering
      % Dynamic (t=1:frame)
      aggregated_traffic=[]
      BW=[0]
      P=[0]
      active=[0]
      colour=[0]
   end
   properties (Constant)
   colour_colours=[{'blue'}, {'red'}, {'green'}, {'magenta'},{'cyan'}, {'black'}, {'yellow'}]
   end
   methods 
       function obj = c(number,center,radius,interfering) %constructor 
         obj.number = number;
         obj.center = center;
         obj.radius = radius;   
         obj.interfering = interfering;   
       end
       function pgon=draw(obj, t)
           pgon = nsidedpoly(1000, 'Center',  [0, 0], 'Radius', obj.radius);
           Latitude=obj.center(1)+(pgon.Vertices(:,2))/110.574; %approximated conversion from km diff. in x and y to degrees in lat. and long.
           Longitude=obj.center(2)+(pgon.Vertices(:,1))./(111.32.*cosd(Latitude));
           %[xc,yc] = centroid(pgon);
           if obj.colour(t)==0
            %plot(pgon, 'FaceColor', char("white"))
            plot(polyshape(Longitude',Latitude'),'FaceColor',char("white"))
           else
            %plot(pgon, 'FaceColor', char(obj.colour_colours(obj.colour(t)))) 
            plot(polyshape(Longitude',Latitude'),'FaceColor',char(obj.colour_colours(obj.colour(t))))
           end
           axis equal
           %text(yc,xc, sprintf(num2str(obj.number), xc, yc), 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontSize',14)
           text(obj.center(2), obj.center(1), sprintf(num2str(obj.number)), 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontSize',14)

       end
       function pgon=draw_latlong(obj, t)
           pgon = nsidedpoly(1000, 'Center',  [0, 0], 'Radius', obj.radius);
           Latitude=obj.center(1)+(pgon.Vertices(:,2))/110.574; %approximated conversion from km diff. in x and y to degrees in lat. and long.
           Longitude=obj.center(2)+(pgon.Vertices(:,1))./(111.32.*cosd(Latitude));
           plot(polyshape(Longitude',Latitude'),'FaceColor','none','EdgeColor', [mod(obj.active(t)+1,2),obj.active(t),0])
       end
       function obj=adduser(obj,u)
           obj.users=[obj.users u];
       end
       function obj=aggregatetraffic(obj)  % FALTA RESTA DE LO YA SERVIDO!
           sum=0;
           for u=1:length(obj.users)
               sum=sum+obj.users(u).traffic_pending;
           end
           obj.aggregated_traffic=[obj.aggregated_traffic sum];
       end
       % function obj = compute_betta_to_sat(obj, h, lat_sp, lon_sp)  % Betta to Satellite from cell center (0,0): meter punto subsatelite como paramtero
       %     Re=6378; % Earth Radius [km]
       %     %traducir el center y el ss de lat-lon  a xyz en ecef y hacer
       %     %sqrt de la resta de cada elemento
       %     %1)
       %     gamma=(sqrt(obj.center(1)^2+obj.center(2)^2)*360)/(2*pi*Re);
       %     %2)
       %     h_prime=Re*(1-cos(gamma*pi/180));
       %     z=sin(gamma*pi/180)*Re;
       %     %3)
       %     betta=(atan(z/(h_prime+h)))*180/pi;
       %     obj.betta_to_sat=betta;
       % end
       function obj = compute_betta_to_sat(obj, h, lat_sp, lon_sp)
        Re = 6378;  % Radio terrestre en km

        % --- Paso 1: convertir centro de celda a ECEF ---
        lat_c = obj.center(1);
        lon_c = obj.center(2);
        alt = 0;  % Asumimos celda sobre la superficie

        % Coordenadas ECEF del centro de celda
        ecef_c = lla2ecef([lat_c, lon_c, alt]);

            for t=1:length(lat_sp)
                % Coordenadas ECEF del subsatélite
                ecef_sp = lla2ecef([lat_sp(t), lon_sp(t), 0]);
        
                % --- Paso 2: calcular ángulo central γ entre los dos puntos ---
                % Normalizar ambos vectores (solo dirección)
                u = ecef_c / norm(ecef_c);
                v = ecef_sp / norm(ecef_sp);
        
                % Ángulo entre ambos vectores (en radianes) con producto escalar
                gamma_rad = acos(dot(u, v));
                gamma = gamma_rad*180/pi;
        
                % --- Paso 3: calcular h' y z ---
                h_prime = Re * (1 - cos(gamma_rad));
                z = sin(gamma_rad) * Re;
        
                % --- Paso 4: calcular β ---
                betta = atan(z / (h_prime + h)) * 180 / pi;
    
                % Guardar en el objeto
                obj.betta_to_sat = [obj.betta_to_sat betta];
            end

        end
   end
end

=======
classdef c < handle
   properties (Access=public)
      number
      center
      radius
      users=[]
      % Computed
      betta_to_sat
      interfering
      % Dynamic (t=1:frame)
      aggregated_traffic=[]
      BW=[0]
      P=[0]
      active=[0]
      colour=[0]
   end
   properties (Constant)
   colour_colours=[{'blue'}, {'red'}, {'green'}, {'magenta'},{'cyan'}, {'black'}, {'yellow'}]
   end
   methods 
       function obj = c(number,center,radius,interfering) %constructor 
         obj.number = number;
         obj.center = center;
         obj.radius = radius;   
         obj.interfering = interfering;   
       end
       function pgon=draw(obj, t)
           pgon = nsidedpoly(1000, 'Center',  [0, 0], 'Radius', obj.radius);
           Latitude=obj.center(1)+(pgon.Vertices(:,2))/110.574; %approximated conversion from km diff. in x and y to degrees in lat. and long.
           Longitude=obj.center(2)+(pgon.Vertices(:,1))./(111.32.*cosd(Latitude));
           %[xc,yc] = centroid(pgon);
           if obj.colour(t)==0
            %plot(pgon, 'FaceColor', char("white"))
            plot(polyshape(Longitude',Latitude'),'FaceColor',char("white"))
           else
            %plot(pgon, 'FaceColor', char(obj.colour_colours(obj.colour(t)))) 
            plot(polyshape(Longitude',Latitude'),'FaceColor',char(obj.colour_colours(obj.colour(t))))
           end
           axis equal
           %text(yc,xc, sprintf(num2str(obj.number), xc, yc), 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontSize',14)
           text(obj.center(2), obj.center(1), sprintf(num2str(obj.number)), 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontSize',14)

       end
       function pgon=draw_latlong(obj, t)
           pgon = nsidedpoly(1000, 'Center',  [0, 0], 'Radius', obj.radius);
           Latitude=obj.center(1)+(pgon.Vertices(:,2))/110.574; %approximated conversion from km diff. in x and y to degrees in lat. and long.
           Longitude=obj.center(2)+(pgon.Vertices(:,1))./(111.32.*cosd(Latitude));
           plot(polyshape(Longitude',Latitude'),'FaceColor','none','EdgeColor', [mod(obj.active(t)+1,2),obj.active(t),0])
       end
       function obj=adduser(obj,u)
           obj.users=[obj.users u];
       end
       function obj=aggregatetraffic(obj)  % FALTA RESTA DE LO YA SERVIDO!
           sum=0;
           for u=1:length(obj.users)
               sum=sum+obj.users(u).traffic_pending;
           end
           obj.aggregated_traffic=[obj.aggregated_traffic sum];
       end
       % function obj = compute_betta_to_sat(obj, h, lat_sp, lon_sp)  % Betta to Satellite from cell center (0,0): meter punto subsatelite como paramtero
       %     Re=6378; % Earth Radius [km]
       %     %traducir el center y el ss de lat-lon  a xyz en ecef y hacer
       %     %sqrt de la resta de cada elemento
       %     %1)
       %     gamma=(sqrt(obj.center(1)^2+obj.center(2)^2)*360)/(2*pi*Re);
       %     %2)
       %     h_prime=Re*(1-cos(gamma*pi/180));
       %     z=sin(gamma*pi/180)*Re;
       %     %3)
       %     betta=(atan(z/(h_prime+h)))*180/pi;
       %     obj.betta_to_sat=betta;
       % end
       function obj = compute_betta_to_sat(obj, h, lat_sp, lon_sp)
        Re = 6378;  % Radio terrestre en km

        % --- Paso 1: convertir centro de celda a ECEF ---
        lat_c = obj.center(1);
        lon_c = obj.center(2);
        alt = 0;  % Asumimos celda sobre la superficie

        % Coordenadas ECEF del centro de celda
        ecef_c = lla2ecef([lat_c, lon_c, alt]);

            for t=1:length(lat_sp)
                % Coordenadas ECEF del subsatélite
                ecef_sp = lla2ecef([lat_sp(t), lon_sp(t), 0]);
        
                % --- Paso 2: calcular ángulo central γ entre los dos puntos ---
                % Normalizar ambos vectores (solo dirección)
                u = ecef_c / norm(ecef_c);
                v = ecef_sp / norm(ecef_sp);
        
                % Ángulo entre ambos vectores (en radianes) con producto escalar
                gamma_rad = acos(dot(u, v));
                gamma = gamma_rad*180/pi;
        
                % --- Paso 3: calcular h' y z ---
                h_prime = Re * (1 - cos(gamma_rad));
                z = sin(gamma_rad) * Re;
        
                % --- Paso 4: calcular β ---
                betta = atan(z / (h_prime + h)) * 180 / pi;
    
                % Guardar en el objeto
                obj.betta_to_sat = [obj.betta_to_sat betta];
            end

        end
   end
end

>>>>>>> e2b063cc1b089e5fb1b6d4fe46a42b06d939b922
=======
classdef c < handle
   properties (Access=public)
      number
      center
      radius
      users=[]
      % Computed
      betta_to_sat
      interfering
      % Dynamic (t=1:frame)
      aggregated_traffic=[]
      BW=[0]
      P=[0]
      active=[0]
      colour=[0]
   end
   properties (Constant)
   colour_colours=[{'blue'}, {'red'}, {'green'}, {'magenta'},{'cyan'}, {'black'}, {'yellow'}]
   end
   methods 
       function obj = c(number,center,radius,interfering) %constructor 
         obj.number = number;
         obj.center = center;
         obj.radius = radius;   
         obj.interfering = interfering;   
       end
       function pgon=draw(obj, t)
           pgon = nsidedpoly(1000, 'Center',  [0, 0], 'Radius', obj.radius);
           Latitude=obj.center(1)+(pgon.Vertices(:,2))/110.574; %approximated conversion from km diff. in x and y to degrees in lat. and long.
           Longitude=obj.center(2)+(pgon.Vertices(:,1))./(111.32.*cosd(Latitude));
           %[xc,yc] = centroid(pgon);
           if obj.colour(t)==0
            %plot(pgon, 'FaceColor', char("white"))
            plot(polyshape(Longitude',Latitude'),'FaceColor',char("white"))
           else
            %plot(pgon, 'FaceColor', char(obj.colour_colours(obj.colour(t)))) 
            plot(polyshape(Longitude',Latitude'),'FaceColor',char(obj.colour_colours(obj.colour(t))))
           end
           axis equal
           %text(yc,xc, sprintf(num2str(obj.number), xc, yc), 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontSize',14)
           text(obj.center(2), obj.center(1), sprintf(num2str(obj.number)), 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontSize',14)

       end
       function pgon=draw_latlong(obj, t)
           pgon = nsidedpoly(1000, 'Center',  [0, 0], 'Radius', obj.radius);
           Latitude=obj.center(1)+(pgon.Vertices(:,2))/110.574; %approximated conversion from km diff. in x and y to degrees in lat. and long.
           Longitude=obj.center(2)+(pgon.Vertices(:,1))./(111.32.*cosd(Latitude));
           plot(polyshape(Longitude',Latitude'),'FaceColor','none','EdgeColor', [mod(obj.active(t)+1,2),obj.active(t),0])
       end
       function obj=adduser(obj,u)
           obj.users=[obj.users u];
       end
       function obj=aggregatetraffic(obj)  % FALTA RESTA DE LO YA SERVIDO!
           sum=0;
           for u=1:length(obj.users)
               sum=sum+obj.users(u).traffic_pending;
           end
           obj.aggregated_traffic=[obj.aggregated_traffic sum];
       end
       % function obj = compute_betta_to_sat(obj, h, lat_sp, lon_sp)  % Betta to Satellite from cell center (0,0): meter punto subsatelite como paramtero
       %     Re=6378; % Earth Radius [km]
       %     %traducir el center y el ss de lat-lon  a xyz en ecef y hacer
       %     %sqrt de la resta de cada elemento
       %     %1)
       %     gamma=(sqrt(obj.center(1)^2+obj.center(2)^2)*360)/(2*pi*Re);
       %     %2)
       %     h_prime=Re*(1-cos(gamma*pi/180));
       %     z=sin(gamma*pi/180)*Re;
       %     %3)
       %     betta=(atan(z/(h_prime+h)))*180/pi;
       %     obj.betta_to_sat=betta;
       % end
       function obj = compute_betta_to_sat(obj, h, lat_sp, lon_sp)
        Re = 6378;  % Radio terrestre en km

        % --- Paso 1: convertir centro de celda a ECEF ---
        lat_c = obj.center(1);
        lon_c = obj.center(2);
        alt = 0;  % Asumimos celda sobre la superficie

        % Coordenadas ECEF del centro de celda
        ecef_c = lla2ecef([lat_c, lon_c, alt]);

            for t=1:length(lat_sp)
                % Coordenadas ECEF del subsatélite
                ecef_sp = lla2ecef([lat_sp(t), lon_sp(t), 0]);
        
                % --- Paso 2: calcular ángulo central γ entre los dos puntos ---
                % Normalizar ambos vectores (solo dirección)
                u = ecef_c / norm(ecef_c);
                v = ecef_sp / norm(ecef_sp);
        
                % Ángulo entre ambos vectores (en radianes) con producto escalar
                gamma_rad = acos(dot(u, v));
                gamma = gamma_rad*180/pi;
        
                % --- Paso 3: calcular h' y z ---
                h_prime = Re * (1 - cos(gamma_rad));
                z = sin(gamma_rad) * Re;
        
                % --- Paso 4: calcular β ---
                betta = atan(z / (h_prime + h)) * 180 / pi;
    
                % Guardar en el objeto
                obj.betta_to_sat = [obj.betta_to_sat betta];
            end

        end
   end
end

>>>>>>> e2b063cc1b089e5fb1b6d4fe46a42b06d939b922
