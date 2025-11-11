function [A,total_range,total_constraints,b, vtype,i_index,m_index,d_index,g_index,f_index,c_index,i_range,m_range,d_range,g_range,f_range,c_range,i_aux,m_aux,d_aux,g_aux,f_aux,c_aux, t_constraint_label]=BH(H,frame,beams,P_T,D,P,R,M, number_cells, nSats, X)

global PWD;
% 
% 
% % MIP Start:
% cd (strcat(PWD,'/FoM_Estimation/'))
% TTL=5;
% 
% [~,~,~,~,~, Ill_db,B_db,P_db]=FOM_calculation_demand_based_band_slots_fixed_MODCODs(B_T,N, P_T, zeros(number_cells,frame), B, zeros(number_cells,frame),c_scenario, beams, theta, frame, frame_dur, TTL, freq,P,sat_idx,K);
% 
% i_start=Ill_db(1:number_cells*frame);
% b_start=B_db(1:n_users*N*frame);
% 
% X_db=zeros(n_users,N,frame);
% M_db=zeros(n_users,M,frame);
% Z_db=zeros(n_users,M,N,frame);
% for t_idx=1:frame
%     for u_idx=1:n_users
%        if P_db(u_idx,t_idx)>0 && sum(B_db(u_idx,:,t_idx))>0
%             X_db(u_idx,sum(B_db(u_idx,:,t_idx)),t_idx)=1;
% 
%             p_idx=find(P_db(u_idx,t_idx)==P(:,sum(B_db(u_idx,:,t_idx)),u_idx));
%             if p_idx==0
%             disp()
%             end
%             M_db(u_idx,p_idx,t_idx)=1;
%             Z_db(u_idx,p_idx,sum(B_db(u_idx,:,t_idx)),t_idx)=1;
%         end
%     end
% end
% 
% x_start=X_db(1:n_users*N*frame);
% m_start=M_db(1:n_users*M*frame);
% z_start=Z_db(1:n_users*M*N*frame);

cd (strcat(PWD,'/Analytical/'))

gamma0=zeros(1,number_cells);
ctts0=zeros(1,number_cells);


% ===== Variables in a single vector x =====
% x = [ i , m , d , gamma , f , ctts ]

% i(c,s,t)  binaries
i_index = 1;
i_range = number_cells * nSats * H*frame;
i_aux   = reshape(i_index:(i_index+i_range-1), [number_cells, nSats, H*frame]);
vtype( i_index:(i_index+i_range-1) ) = 'B';

% m(c,s,m,t)  binaries
m_index = i_index + i_range;
m_range = number_cells * nSats * M * H*frame;
m_aux   = reshape(m_index:(m_index+m_range-1), [number_cells, nSats, M, H*frame]);
vtype( m_index:(m_index+m_range-1) ) = 'B';

% d(c,t)  continuous
d_index = m_index + m_range;
d_range = number_cells * H*frame;
d_aux   = reshape(d_index:(d_index+d_range-1), [number_cells, H*frame]);
vtype( d_index:(d_index+d_range-1) ) = 'C';

% gamma(c,t)  continuous
g_index = d_index + d_range;
g_range = number_cells * H*frame;
g_aux   = reshape(g_index:(g_index+g_range-1), [number_cells, H*frame]);
vtype( g_index:(g_index+g_range-1) ) = 'C';

% f(c,t)  continuous
f_index = g_index + g_range;
f_range = number_cells * H*frame;
f_aux   = reshape(f_index:(f_index+f_range-1), [number_cells, H*frame]);
vtype( f_index:(f_index+f_range-1) ) = 'C';

% ctts(c,t)  continuous
c_index = f_index + f_range;
c_range = number_cells * H*frame;
c_aux   = reshape(c_index:(c_index+c_range-1), [number_cells, H*frame]);
vtype( c_index:(c_index+c_range-1) ) = 'C';

vtype = vtype';
total_range = i_range + m_range + d_range + g_range + f_range + c_range;
lb = zeros(total_range,1);
ub = inf(total_range,1);
ub( vtype=='B' ) = 1;   % binaries in [0,1]


% ===== Constraints matrix Ax, b =====

% Rough upper bound of constraints for spalloc; 
rough_rows = ...
  (nnz(X) * (1 + 1)) ...    % gating: i<=X, sum_m m <= i
+ (nSats*H*frame) ...             % illum budget
+ (nSats*H*frame) ...             % power budget
+ (number_cells*H*frame) ...      % d definition (=)
+ 3*(number_cells*H*frame);       % f,gamma,ctts updates

A = spalloc(rough_rows, total_range, rough_rows*5);
b = zeros(rough_rows,1);
sense = repmat('<', rough_rows, 1);   % will modify per-row
t_constraint_label = zeros(rough_rows,1);

nf = 1;

% (1) Gating + single-MCS (only if X(c,s,t)==1)

for t = 1:H*frame
  for s = 1:nSats
    for c = 1:number_cells
      if X(c,s,t)==1
        % i(c,s,t) <= X(c,s,t)
        A(nf, i_aux(c,s,t)) = 1;
        b(nf) = X(c,s,t);
        sense(nf) = '<';
        t_constraint_label(nf) = t; nf = nf+1;

        % sum_m m(c,s,m,t) - i(c,s,t) <= 0
        for mm = 1:M
          A(nf, m_aux(c,s,mm,t)) = 1;
        end
        A(nf, i_aux(c,s,t)) = -1;
        b(nf) = 0; sense(nf) = '<'; t_constraint_label(nf) = t; nf = nf+1;
      end
    end
  end
end

% (2) Simultaneous illumination limit: ∑c​i(c,s,t)≤Nill(s)​

for t = 1:H*frame
  for s = 1:nSats
    for c = 1:number_cells
      if X(c,s,t)==1
        A(nf, i_aux(c,s,t)) = 1;
      end
    end
    b(nf) = beams;   
    sense(nf) = '<'; t_constraint_label(nf) = t; nf = nf+1;
  end
end

% (3) Power budget ∑c,mP(c,m,t)m(c,s,m,t)≤P_T
for t = 1:H*frame
  for s = 1:nSats
    for c = 1:number_cells
      if X(c,s,t)==1
        for mm = 1:M
          A(nf, m_aux(c,s,mm,t)) = P{s}(c,mm,t);
        end
      end
    end
    b(nf) = P_T(s);    % scalar or vector per sat
    sense(nf) = '<'; t_constraint_label(nf) = t; nf = nf+1;
  end
end

% (4) Served capacity definition d(c,t)=∑s,mD(c,m,t)m(c,s,m,t)
for t = 1:H*frame
  for c = 1:number_cells
    % left: d(c,t)
    A(nf, d_aux(c,t)) = 1;
    % right: - sum_{s: X=1} sum_m D * m
    for s = 1:nSats
      if X(c,s,t)==1
        for mm = 1:M
          A(nf, m_aux(c,s,mm,t)) = -D{s}(c,mm,t);
        end
      end
    end
    b(nf) = 0;
    sense(nf) = '='; t_constraint_label(nf) = t; nf = nf+1;
  end
end
% (5) Global state updates: γ(c,t), f(c,t), c(c,t)
for t = 1:H*frame
  for c = 1:number_cells
        % ---- f rows ----
        if t==1
          % f(c,1) - d(c,1) >= -(R + gamma0(c))
          A(nf, f_aux(c,1)) =  1;
          A(nf, d_aux(c,1)) = -1;
          b(nf) = -( R(c,1) + gamma0(c) );
          sense(nf) = '>'; nf = nf+1;
        else
          % f(c,t) + gamma(c,t-1) - d(c,t) >= R(c,t)
          A(nf, f_aux(c,t))     =  1;
          A(nf, g_aux(c,t-1))   =  1;
          A(nf, d_aux(c,t))     = -1;
          b(nf) = R(c,t);
          sense(nf) = '>'; nf = nf+1;
        end
        
        % ---- gamma rows ----
        if t==1
          % gamma(c,1) + d(c,1) >= R + gamma0(c)
          A(nf, g_aux(c,1)) = 1;
          A(nf, d_aux(c,1)) = 1;
          b(nf) = R(c,1) + gamma0(c);
          sense(nf) = '>'; nf = nf+1;
        else
          % gamma(c,t) + d(c,t) - gamma(c,t-1) >= R(c,t)
          A(nf, g_aux(c,t))   =  1;
          A(nf, d_aux(c,t))   =  1;
          A(nf, g_aux(c,t-1)) = -1;
          b(nf) = R(c,t);
          sense(nf) = '>'; nf = nf+1;
        end
        
        % ---- TTS rows ----
        Rt = max(R(c,t), epsR);
        if t==1
          % ctts(c,1) - (1/Rt)*gamma(c,1) + N_T * sum_s i(c,s,1) >= ctts0(c)
          A(nf, c_aux(c,1)) = 1;
          A(nf, g_aux(c,1)) = -1.0/Rt;
          for s=1:nSats
            if X(c,s,1)==1, A(nf, i_aux(c,s,1)) = H*frame; end
          end
          b(nf) = ctts0(c);
          sense(nf) = '>'; nf = nf+1;
        else
          % ctts(c,t) - (1/Rt)*gamma(c,t) - ctts(c,t-1) + N_T * sum_s i(c,s,t) >= 0
          A(nf, c_aux(c,t))   =  1;
          A(nf, g_aux(c,t))   = -1.0/Rt;
          A(nf, c_aux(c,t-1)) = -1;
          for s=1:nSats
            if X(c,s,t)==1, A(nf, i_aux(c,s,t)) = H*frame; end
          end
          b(nf) = 0;
          sense(nf) = '>'; nf = nf+1;
        end
  end
end
total_constraints=nf-1;
end