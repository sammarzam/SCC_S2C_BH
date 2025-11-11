function [A,total_range,total_constraints,b,sense,vtype, ...
          i_index,m_index,g_index,f_index,c_index, ...
          i_range,m_range,g_range,f_range,c_range, ...
          i_aux,m_aux,g_aux,f_aux,c_aux, ...
          t_constraint_label, ...
          tlabel_i, tlabel_m, tlabel_g, tlabel_f, tlabel_c] = ...
          BH_vectored(H,frame,beams,P_T,D,P,R,M, number_cells, nSats, X,m_continuous)
% BH_vectored
% Builds the global MILP (Ω-compact) for the BH planning:
%   Variables (single long x):
%     i(c,t)           ∈ {0,1}                  illumination decision for (cell,time)
%     m(c,mcs,t)       ∈ {0,1}                  modulation choice for illuminated (cell,time)
%     gamma(c,t) ≥ 0   (continuous)             backlog-like ledger
%     f(c,t)     ≥ 0   (continuous)             oversupply ledger
%     ctts(c,t)  ≥ 0   (continuous by default)  timing ledger (can be set integer if desired)
%
%   We use “Ω-compact” indexing: only one i per (c,t) and only M m’s per (c,t)
%   for the “responsible” satellite s* (taken from argmax(X)).
%
% Inputs must be in “outer” time H and are expanded to inner time Hf = H*frame.

% ===== Sizes =====
Hf   = H * frame;              % inner BH slots
CT   = number_cells * Hf;      % total (c,t) pairs

% ===== Normalize inputs to Hf =====
% X: [cells x nSats x H] -> repeat per inner slot
if ~isequal(size(X), [number_cells, nSats, H])
  error('X must be [number_cells x nSats x H].');
end
Xh  = X;
X   = repelem(Xh, 1, 1, frame);     % -> [cells x nSats x Hf]

% R: [cells x Hf] (already inner-time)
if ~isequal(size(R), [number_cells, Hf])
  error('R must be [number_cells x H*frame].');
end

% D, P are cell arrays per satellite: each [cells x M x H], expand to Hf
if ~iscell(D) || numel(D)~=nSats, error('D must be cell array (nSats).'); end
if ~iscell(P) || numel(P)~=nSats, error('P must be cell array (nSats).'); end
for s=1:nSats
  if ~isequal(size(D{s}), [number_cells, M, H]), error('D{%d} size mismatch', s); end
  if ~isequal(size(P{s}), [number_cells, M, H]), error('P{%d} size mismatch', s); end
  D{s} = repelem(D{s}, 1,1, frame); % -> [cells x M x Hf]
  P{s} = repelem(P{s}, 1,1, frame); % -> [cells x M x Hf]
end

% Power budget per satellite
if numel(P_T)==1
  P_T_vec = repmat(P_T, nSats,1);
elseif numel(P_T)==nSats
  P_T_vec = P_T(:);
else
  error('P_T must be scalar or length nSats.');
end

% ===== Responsibility s*(c,t) and (c,t) helpers =====
% Responsible satellite at each (c,t): s* = argmax_s X(c,s,t)
[~, s_resp] = max(X, [], 2);                % [cells x 1 x Hf]
s_resp = squeeze(s_resp);                    % [cells x Hf]
s_star = s_resp(:);                          % [CT x 1]

% Build (c,t) linearization consistent with reshape(...,[number_cells,Hf])
[Cmat, Tmat] = ndgrid(1:number_cells, 1:Hf); % [cells,Hf]
c_rep = Cmat(:);                              % [CT x 1] cells 1..C repeated across t
t_rep = Tmat(:);                              % [CT x 1] times 1..Hf, each repeated for all cells

% Flatten demand to match CT linearization
Rvec = reshape(R, [CT,1]);

% ===== Variables (Ω-compact) =====
% NOTE: by default, gamma, f, ctts are continuous ('C').
% If you want ctts to be an integer “counter”, flip vtype for the c-block to 'I'.

% -- i(c,t) in {0,1}
i_index = 1;
i_range = CT;
i_aux   = reshape(i_index:(i_index+i_range-1), [number_cells, Hf]);
tlabel_i = t_rep;

% -- m(c,mcs,t) in {0,1}, only for responsible sat s*
m_index = i_index + i_range;
m_range = CT * M;
m_aux   = reshape(m_index:(m_index+m_range-1), [number_cells, M, Hf]);
tlabel_m = repelem(t_rep, M, 1);

% -- gamma(c,t) >= 0 (continuous)
g_index = m_index + m_range;
g_range = CT;
g_aux   = reshape(g_index:(g_index+g_range-1), [number_cells, Hf]);
tlabel_g = t_rep;

% -- f(c,t) >= 0 (continuous)
f_index = g_index + g_range;
f_range = CT;
f_aux   = reshape(f_index:(f_index+f_range-1), [number_cells, Hf]);
tlabel_f = t_rep;

% -- ctts(c,t) >= 0 (continuous by default)
c_index = f_index + f_range;
c_range = CT;
c_aux   = reshape(c_index:(c_index+c_range-1), [number_cells, Hf]);
tlabel_c = t_rep;

% total variables
total_range = i_range + m_range + g_range + f_range + c_range;

% vtypes
vtype = repmat('C', total_range, 1);
vtype(i_index:(i_index+i_range-1)) = 'B';
if m_continuous==1
    vtype(m_index:(m_index+m_range-1)) = 'C';
else
    vtype(m_index:(m_index+m_range-1)) = 'B';
end
%  ctts integer:
vtype(c_index:(c_index+c_range-1)) = 'I';

% ===== Assemble constraints A x (<=,=,>=) b =====
I = []; J = []; V = [];   % sparse triplets
B = []; S = []; TL = [];  % rhs, sense, time label per row
row = 0;

% ----------------------------
% (1) MCS gating (exactly one MCS if illuminated)
% ----------------------------
% (1a)  sum_m m(c,·,t) - i(c,t) <= 0
r1 = (row+1 : row+CT)'; row = row + CT;

% + sum_m m
I = [I; repelem(r1, M)];
J = [J; m_aux(sub2ind(size(m_aux), ...
         repelem(c_rep, M,1), ...
         repmat((1:M)', CT,1), ...
         repelem(t_rep, M,1)))];
V = [V; ones(CT*M,1)];

% - i(c,t)
I = [I; r1];
J = [J; i_aux(sub2ind(size(i_aux), c_rep, t_rep))];
V = [V; -ones(CT,1)];

B = [B; zeros(CT,1)];
S = [S; repmat('<', CT,1)];
TL= [TL; t_rep];

% (1b)  i(c,t) - sum_m m(c,·,t) <= 0    ⇒ together with (1a):  sum_m m = i
r1b = (row+1 : row+CT)'; row = row + CT;

% + i(c,t)
I = [I; r1b];
J = [J; i_aux(sub2ind(size(i_aux), c_rep, t_rep))];
V = [V; ones(CT,1)];

% - sum_m m
I = [I; repelem(r1b, M)];
J = [J; m_aux(sub2ind(size(m_aux), ...
         repelem(c_rep, M,1), ...
         repmat((1:M)', CT,1), ...
         repelem(t_rep, M,1)))];
V = [V; -ones(CT*M,1)];

B = [B; zeros(CT,1)];
S = [S; repmat('<', CT,1)];
TL= [TL; t_rep];

% ----------------------------
% (2) Illumination budget per (s,t): sum_{c: s*(c,t)=s} i(c,t) <= beams
% ----------------------------
key_st = (t_rep-1)*nSats + s_star;          % encode (t,s) group
[uniq_st, ~, grp_st] = unique(key_st);
nGroups = numel(uniq_st);

r2 = (row+1 : row+nGroups)'; row = row + nGroups;

I = [I; r2(grp_st)];
J = [J; i_aux(sub2ind(size(i_aux), c_rep, t_rep))];
V = [V; ones(CT,1)];

B = [B; beams * ones(nGroups,1)];
S = [S; repmat('<', nGroups,1)];
TL= [TL; floor((uniq_st-1)/nSats)+1];   % label by time

% ----------------------------
% (3) Power budget per (s,t): sum_{c,m} P_s(c,m,t) m(c,m,t) <= P_T(s)
% ----------------------------
r3 = (row+1 : row+nGroups)'; row = row + nGroups;

% repeat row IDs per each (c,m,t) in that (s,t) group
I = [I; repelem(r3(grp_st), M)];

% same m indices as above, grouped by (s,t)
Jm_all = m_aux(sub2ind(size(m_aux), ...
         repelem(c_rep, M,1), ...
         repmat((1:M)', CT,1), ...
         repelem(t_rep, M,1)));
J = [J; Jm_all];

% coefficients from P for the responsible satellite s*
Pcoef = zeros(CT*M,1);
ptr = 1;
for idx=1:CT
  cc = c_rep(idx); ss = s_star(idx); tt = t_rep(idx);
  Pcoef(ptr:ptr+M-1) = reshape(P{ss}(cc, :, tt), [M,1]);
  ptr = ptr + M;
end
V = [V; Pcoef];

% RHS per group uses P_T(s)
s_of_group = mod(uniq_st-1, nSats) + 1;
B = [B; P_T_vec(s_of_group)];
S = [S; repmat('<', nGroups,1)];
TL= [TL; floor((uniq_st-1)/nSats)+1];

% ----------------------------
% Helper: served capacity Σ_m D * m for a set of (c,t)
% ----------------------------
getDm = @(ct_vec) deal( ...
  m_aux(sub2ind(size(m_aux), ...
        repelem(c_rep(ct_vec), M,1), ...
        repmat((1:M)', numel(ct_vec),1), ...
        repelem(t_rep(ct_vec), M,1))), ...
  localDcoeff(ct_vec, c_rep, s_star, t_rep, D, M) );

% Basic (c,t) indexing blocks
ct_idx = reshape(1:CT, number_cells, Hf);  % matrix of linear ct indices
t1_idx = ct_idx(:,1);                      % all (c,1)
curr   = ct_idx(:,2:end);                  % all (c,t>=2)
prev   = ct_idx(:,1:end-1);                % aligned (c,t-1)
curr_v = curr(:);
prev_v = prev(:);

% Flattened variable indices for quick addressing
idx_g = g_aux(:);
idx_f = f_aux(:);
idx_c = c_aux(:);

% Optional initial conditions
gamma0 = zeros(number_cells,1);
ctts0  = zeros(number_cells,1);

% ----------------------------
% (4) Ledgers using Σ D m (no explicit d variable)
% ----------------------------
% 4a) f rows  (oversupply)
% Desired balance:  f_t  ≥  s_t - (γ_{t-1} + R_t)
%                  ⇔ f_t + γ_{t-1} - s_t ≥ R_t
% Below we implement the equivalent:  f + γ_prev - ΣDm - R ≥ 0  ⇒ RHS = -R

% -- t = 1: f(1) + gamma0 - ΣDm - R(1) ≥ 0
r4a = (row+1 : row+numel(t1_idx))'; row = row + numel(t1_idx);
I = [I; r4a];               J = [J; idx_f(t1_idx)];      V = [V;  ones(numel(t1_idx),1)];
[JDm, VDm] = getDm(t1_idx(:));
I = [I; r4a(repelem(1:numel(t1_idx), M)')];  J = [J; JDm];  V = [V; -VDm];
% put +gamma0 on LHS by leaving it implicit (since gamma0 is data), move R to RHS as -R
B = [B; -( Rvec(t1_idx) + gamma0((1:number_cells)') )];
S = [S; repmat('>', numel(t1_idx),1)];
TL= [TL; ones(numel(t1_idx),1)];

% -- t ≥ 2: f(t) + γ(t-1) - ΣDm ≥ R(t)  ⇔ f + γ_prev - ΣDm ≥ -RHS with RHS = -R
r4b = (row+1 : row+numel(curr_v))'; row = row + numel(curr_v);
I = [I; r4b; r4b]; 
J = [J; idx_f(curr_v); idx_g(prev_v)];
V = [V;  ones(numel(curr_v),1);  ones(numel(curr_v),1)];
[JDm, VDm] = getDm(curr_v);
I = [I; r4b(repelem(1:numel(curr_v), M)')];  J = [J; JDm];  V = [V; -VDm];
B = [B; -Rvec(curr_v)];                      % <— your requested form
S = [S; repmat('>', numel(curr_v),1)];
TL= [TL; t_rep(curr_v)];

% 4b) gamma rows  (backlog)
% Desired balance: γ_t ≥ γ_{t-1} + R_t - s_t
% Implemented as:  γ_t + s_t - γ_{t-1} ≥ R_t

% -- t = 1: γ(1) + ΣDm ≥ R(1) + gamma0
r5a = (row+1 : row+numel(t1_idx))'; row = row + numel(t1_idx);
I = [I; r5a];               J = [J; idx_g(t1_idx)];      V = [V;  ones(numel(t1_idx),1)];
[JDm, VDm] = getDm(t1_idx(:));
I = [I; r5a(repelem(1:numel(t1_idx), M)')];  J = [J; JDm];  V = [V;  VDm];
B = [B; Rvec(t1_idx) + gamma0((1:number_cells)')];
S = [S; repmat('>', numel(t1_idx),1)];
TL= [TL; ones(numel(t1_idx),1)];

% -- t ≥ 2: γ(t) + ΣDm - γ(t-1) ≥ R(t)
r5b = (row+1 : row+numel(curr_v))'; row = row + numel(curr_v);
I = [I; r5b; r5b]; 
J = [J; idx_g(curr_v); idx_g(prev_v)];
V = [V;  ones(numel(curr_v),1); -ones(numel(curr_v),1)];
[JDm, VDm] = getDm(curr_v);
I = [I; r5b(repelem(1:numel(curr_v), M)')];  J = [J; JDm];  V = [V;  VDm];
B = [B; Rvec(curr_v)];
S = [S; repmat('>', numel(curr_v),1)];
TL= [TL; t_rep(curr_v)];

% NEW! TESTING:
% --- (4b.1) gamma "no-overshoot": γ_t - γ_{t-1} ≤ R_t  (caps growth)
% t = 1:  γ(c,1) - gamma0(c) ≤ R(c,1)
r5c1 = (row+1 : row+numel(t1_idx))'; row = row + numel(t1_idx);
I = [I; r5c1];                      J = [J; idx_g(t1_idx)];   V = [V;  ones(numel(t1_idx),1)];
% -gamma0 goes to RHS:  RHS = R + gamma0
B = [B; Rvec(t1_idx) + gamma0((1:number_cells)')];
S = [S; repmat('<', numel(t1_idx),1)];
TL= [TL; ones(numel(t1_idx),1)];

% t ≥ 2:  γ(c,t) - γ(c,t-1) ≤ R(c,t)
r5c2 = (row+1 : row+numel(curr_v))'; row = row + numel(curr_v);
I = [I; r5c2; r5c2];
J = [J; idx_g(curr_v);  idx_g(prev_v)];
V = [V;  ones(numel(curr_v),1);  -ones(numel(curr_v),1)];
B = [B; Rvec(curr_v)];
S = [S; repmat('<', numel(curr_v),1)];
TL= [TL; t_rep(curr_v)];

% 4c) ctts rows  (R-aware, avoids 1/R when R=0)
% We keep your split:
%   t==1 and t>=2, each with cases R>0 and R==0

% R: [number_cells x Hf]
Rcum      = cumsum(R, 2);                 % cumulative per cell over time
Rcum_vec  = reshape(Rcum, [], 1);         % [CT x 1], aligned with (c,t) linearization
pos_mask  = (Rcum_vec > 0);               % guard where cumulative is positive

% -- t == 1 subsets
t1_pos  = t1_idx( pos_mask(t1_idx) );
t1_zero = t1_idx(~pos_mask(t1_idx));

% i(c,1) indices only for those subsets (not all cells)
i_t1_pos  = i_aux(sub2ind(size(i_aux), c_rep(t1_pos),  t_rep(t1_pos)));
i_t1_zero = i_aux(sub2ind(size(i_aux), c_rep(t1_zero), t_rep(t1_zero)));

% t==1, R>0:  c - (1/R)·γ + Hf·i >= ctts0
r6a_pos = (row+1 : row+numel(t1_pos))'; row = row + numel(t1_pos);
I = [I; r6a_pos; r6a_pos];
J = [J; idx_c(t1_pos); idx_g(t1_pos)];
V = [V; ones(numel(t1_pos),1); -1 ./ Rcum_vec(t1_pos)];
I = [I; r6a_pos];
J = [J; i_t1_pos];
V = [V; Hf * ones(numel(t1_pos),1)];
B = [B; ctts0( c_rep(t1_pos) )];
S = [S; repmat('>', numel(t1_pos),1)];
TL= [TL; t_rep(t1_pos)];

% t==1, R==0:  c + Hf·i >= ctts0 (no gamma term)
r6a_zero = (row+1 : row+numel(t1_zero))'; row = row + numel(t1_zero);
I = [I; r6a_zero];
J = [J; idx_c(t1_zero)];
V = [V; ones(numel(t1_zero),1)];
I = [I; r6a_zero];
J = [J; i_t1_zero];
V = [V; Hf * ones(numel(t1_zero),1)];
B = [B; ctts0( c_rep(t1_zero) )];
S = [S; repmat('>', numel(t1_zero),1)];
TL= [TL; t_rep(t1_zero)];

% -- t >= 2 subsets
curr_pos  = curr_v( pos_mask(curr_v) );
prev_pos  = prev_v( pos_mask(curr_v) );
curr_zero = curr_v(~pos_mask(curr_v));
prev_zero = prev_v(~pos_mask(curr_v));

% t>=2, R>0:  c - (1/R)·γ - c_prev + Hf·i >= 0
r6b_pos = (row+1 : row+numel(curr_pos))'; row = row + numel(curr_pos);
I = [I; r6b_pos; r6b_pos; r6b_pos];
J = [J; idx_c(curr_pos); idx_g(curr_pos); idx_c(prev_pos)];
V = [V; ones(numel(curr_pos),1); -1 ./ Rcum_vec(curr_pos); -ones(numel(curr_pos),1)];
i_curr = i_aux(sub2ind(size(i_aux), c_rep(curr_pos), t_rep(curr_pos)));
I = [I; r6b_pos];
J = [J; i_curr];
V = [V; Hf * ones(numel(curr_pos),1)];
B = [B; zeros(numel(curr_pos),1)];
S = [S; repmat('>', numel(curr_pos),1)];
TL= [TL; t_rep(curr_pos)];

% t>=2, R==0:  c - c_prev + Hf·i >= 0
r6b_zero = (row+1 : row+numel(curr_zero))'; row = row + numel(curr_zero);
I = [I; r6b_zero; r6b_zero];
J = [J; idx_c(curr_zero); idx_c(prev_zero)];
V = [V; ones(numel(curr_zero),1); -ones(numel(curr_zero),1)];
i_curr0 = i_aux(sub2ind(size(i_aux), c_rep(curr_zero), t_rep(curr_zero)));
I = [I; r6b_zero];
J = [J; i_curr0];
V = [V; Hf * ones(numel(curr_zero),1)];
B = [B; zeros(numel(curr_zero),1)];
S = [S; repmat('>', numel(curr_zero),1)];
TL= [TL; t_rep(curr_zero)];

% ===== Finalize =====
total_constraints = row;

% Safety checks
assert(numel(I)==numel(J) && numel(J)==numel(V), 'Triplet length mismatch');
assert(length(B)==total_constraints && length(S)==total_constraints && ...
       length(TL)==total_constraints, 'Row arrays mismatch');

A = sparse(I, J, V, total_constraints, total_range);
b = B;
sense = S.';                 % Gurobi expects char row
t_constraint_label = TL;     % per-row time label

end  % BH_vectored


% ---- helper: pull D coefficients for responsible sat on a list of ct ----
function VDm = localDcoeff(ct_vec, c_rep, s_star, t_rep, D, M)
  VDm = zeros(numel(ct_vec)*M,1);
  ptr = 1;
  for k=1:numel(ct_vec)
    idx = ct_vec(k);
    cc  = c_rep(idx); ss = s_star(idx); tt = t_rep(idx);
    VDm(ptr:ptr+M-1) = reshape(D{ss}(cc,:,tt), [M,1]);
    ptr = ptr + M;
  end
end
