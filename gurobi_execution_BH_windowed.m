<<<<<<< HEAD
<<<<<<< HEAD
function [result, data, solutions] = gurobi_execution_BH_windowed( ...
    betta, normalization_UC, normalization_EC, normalization_time, ...
    total_range, A, b, sense, vtype, MIPGap, ...
    mip_start, lb, ub, ...
    i_index, i_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, ...
    t_row_label, tlabel_i, tlabel_m, tlabel_g, tlabel_f, tlabel_c, ...
    t_current, W, R)

% -------- weights from beta (same logic as before) ----------
weight_traffic_EC = (1 - betta)/2;
weight_traffic_UC = betta;
weight_timing     = (1 - betta)/2;

% -------- window bounds ----------
Hf = max(t_row_label);                    % total inner slots in the built model
t0 = t_current;
t1 = min(t_current + W - 1, Hf);

% -------- select window constraints ----------
rows_keep = (t_row_label >= t0) & (t_row_label <= t1);
Awin  = A(rows_keep, :);
bwin  = b(rows_keep, :);
sense_win = sense(rows_keep);

% -------- decision variable selection per window ----------
%+ past (< t0): fix to incumbent (ub==lb==mip_start)
%+ in window [t0..t1]: free (binaries [0,1], continuous [0,∞))
%+ future (> t1): hard-off (ub==lb==0) — not strictly required since their rows are dropped, but it tightens the model and speeds the solver
past_i   = (tlabel_i <  t0);   in_i = (tlabel_i >= t0 & tlabel_i <= t1);   fut_i = (tlabel_i >  t1);
past_m   = (tlabel_m <  t0);   in_m = (tlabel_m >= t0 & tlabel_m <= t1);   fut_m = (tlabel_m >  t1);
past_g   = (tlabel_g <  t0);   in_g = (tlabel_g >= t0 & tlabel_g <= t1);   fut_g = (tlabel_g >  t1);
past_f   = (tlabel_f <  t0);   in_f = (tlabel_f >= t0 & tlabel_f <= t1);   fut_f = (tlabel_f >  t1);
past_c   = (tlabel_c <  t0);   in_c = (tlabel_c >= t0 & tlabel_c <= t1);   fut_c = (tlabel_c >  t1);

% build full-length logical masks aligned with x for window
mask_in_i = false(total_range,1);  mask_in_i(i_index:(i_index+i_range-1)) = in_i;
mask_in_m = false(total_range,1);  mask_in_m(m_index:(m_index+m_range-1)) = in_m;
mask_in_g = false(total_range,1);  mask_in_g(g_index:(g_index+g_range-1)) = in_g;
mask_in_f = false(total_range,1);  mask_in_f(f_index:(f_index+f_range-1)) = in_f;
mask_in_c = false(total_range,1);  mask_in_c(c_index:(c_index+c_range-1)) = in_c;

mask_past = false(total_range,1);
mask_past(i_index:(i_index+i_range-1)) = past_i;
mask_past(m_index:(m_index+m_range-1)) = past_m;
mask_past(g_index:(g_index+g_range-1)) = past_g;
mask_past(f_index:(f_index+f_range-1)) = past_f;
mask_past(c_index:(c_index+c_range-1)) = past_c;

mask_fut = false(total_range,1);
mask_fut(i_index:(i_index+i_range-1)) = fut_i;
mask_fut(m_index:(m_index+m_range-1)) = fut_m;
mask_fut(g_index:(g_index+g_range-1)) = fut_g;
mask_fut(f_index:(f_index+f_range-1)) = fut_f;
mask_fut(c_index:(c_index+c_range-1)) = fut_c;

% initialize bounds from incumbent (carry past decisions by default)
if isempty(mip_start), mip_start = zeros(total_range,1); end
lb(:) = mip_start;     % default: fixed to incumbent
ub(:) = mip_start;

% free in-window blocks
lb(mask_in_i) = 0.0;           ub(mask_in_i) = 1.0;          % binaries
lb(mask_in_m) = 0.0;           ub(mask_in_m) = 1.0;          % binaries

% choose safe caps for γ and f (scalar is OK)
U_cap = max(sum(R(:,1:t1),2)) + 1;   % scalar upper bound
lb(mask_in_g) = 0.0;           ub(mask_in_g) = U_cap; %inf;       % g
lb(mask_in_f) = 0.0;           ub(mask_in_f) = inf;       % f

% for ctts, allow growth up to current time
lb(mask_in_c) = 0.0;           ub(mask_in_c) = t1 + 1; %inf      % c 

% hard-off future vars
lb(mask_fut) = 0.0;
ub(mask_fut) = 0.0;

% % (optional) explicitly fix the carry-in ledgers at t0-1 (redundant but clear)
% if t0 > 1
%     fix_g_prev = (tlabel_g == (t0-1));
%     fix_c_prev = (tlabel_c == (t0-1));
%     idx_g_prev = g_index - 1 + find(fix_g_prev);
%     idx_c_prev = c_index - 1 + find(fix_c_prev);
%     lb(idx_g_prev) = mip_start(idx_g_prev);
%     ub(idx_g_prev) = mip_start(idx_g_prev);
%     lb(idx_c_prev) = mip_start(idx_c_prev);
%     ub(idx_c_prev) = mip_start(idx_c_prev);
% end


% -------- objective (only on window slots) ----------
obj = zeros(total_range,1);

obj(g_index:(g_index+g_range-1)) = (weight_traffic_UC / normalization_UC) * double(tlabel_g >= t0 & tlabel_g <= t1);
obj(f_index:(f_index+f_range-1)) = (weight_traffic_EC / normalization_EC) * double(tlabel_f >= t0 & tlabel_f <= t1);
obj(c_index:(c_index+c_range-1)) = (weight_timing     / normalization_time) * double(tlabel_c >= t0 & tlabel_c <= t1);

% -------- Build Gurobi model ----------
model.A          = Awin;          % sparse
model.rhs        = bwin;          % double
model.sense      = sense_win;     % char row
model.obj        = obj;
model.modelsense = 'min';
model.vtype      = vtype;         % char column
model.lb         = lb;
model.ub         = ub;
% model.start      = mip_start;

% Parameters
params.OutputFlag = 1;
params.MIPGap     = MIPGap;
% (you can also set: params.MIPFocus=1; params.Heuristics=0.3; params.Cuts=1;)

% ---- Solve (MATLAB API) ----
% Set path if needed:
if ispc
    addpath('C:\gurobi1100\win64\matlab');
elseif isunix
    addpath('/media/apps/avx512-2021/software/Gurobi/10.0.1-GCCcore-12.2.0/matlab');
end

result = gurobi(model, params);

% Cleanup occasionally: unloads all compiled MEX binaries from memory (Matlab MexFunction crash) 
if (mod(t_current-1, W*10) == 0) && (betta~=0) && (betta~=1) 
    clear mex; clear gurobi;
    fprintf('Gurobi reloaded at iteration %d\n', t_current);
end

% Make outputs consistent with your previous signature
data      = [];     % (no callback trace here)
solutions = [];     % (if you want incumbent history, we can add a callback later)

% Convert to doubles (defensive)
if isfield(result,'x'), result.x = double(result.x); end
end
=======
function [result, data, solutions] = gurobi_execution_BH_windowed( ...
    betta, normalization_UC, normalization_EC, normalization_time, ...
    total_range, A, b, sense, vtype, MIPGap, ...
    mip_start, lb, ub, ...
    i_index, i_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, ...
    t_row_label, tlabel_i, tlabel_m, tlabel_g, tlabel_f, tlabel_c, ...
    t_current, W, R)

% -------- weights from beta (same logic as before) ----------
weight_traffic_EC = (1 - betta)/2;
weight_traffic_UC = betta;
weight_timing     = (1 - betta)/2;

% -------- window bounds ----------
Hf = max(t_row_label);                    % total inner slots in the built model
t0 = t_current;
t1 = min(t_current + W - 1, Hf);

% -------- select window constraints ----------
rows_keep = (t_row_label >= t0) & (t_row_label <= t1);
Awin  = A(rows_keep, :);
bwin  = b(rows_keep, :);
sense_win = sense(rows_keep);

% -------- decision variable selection per window ----------
%+ past (< t0): fix to incumbent (ub==lb==mip_start)
%+ in window [t0..t1]: free (binaries [0,1], continuous [0,∞))
%+ future (> t1): hard-off (ub==lb==0) — not strictly required since their rows are dropped, but it tightens the model and speeds the solver
past_i   = (tlabel_i <  t0);   in_i = (tlabel_i >= t0 & tlabel_i <= t1);   fut_i = (tlabel_i >  t1);
past_m   = (tlabel_m <  t0);   in_m = (tlabel_m >= t0 & tlabel_m <= t1);   fut_m = (tlabel_m >  t1);
past_g   = (tlabel_g <  t0);   in_g = (tlabel_g >= t0 & tlabel_g <= t1);   fut_g = (tlabel_g >  t1);
past_f   = (tlabel_f <  t0);   in_f = (tlabel_f >= t0 & tlabel_f <= t1);   fut_f = (tlabel_f >  t1);
past_c   = (tlabel_c <  t0);   in_c = (tlabel_c >= t0 & tlabel_c <= t1);   fut_c = (tlabel_c >  t1);

% build full-length logical masks aligned with x for window
mask_in_i = false(total_range,1);  mask_in_i(i_index:(i_index+i_range-1)) = in_i;
mask_in_m = false(total_range,1);  mask_in_m(m_index:(m_index+m_range-1)) = in_m;
mask_in_g = false(total_range,1);  mask_in_g(g_index:(g_index+g_range-1)) = in_g;
mask_in_f = false(total_range,1);  mask_in_f(f_index:(f_index+f_range-1)) = in_f;
mask_in_c = false(total_range,1);  mask_in_c(c_index:(c_index+c_range-1)) = in_c;

mask_past = false(total_range,1);
mask_past(i_index:(i_index+i_range-1)) = past_i;
mask_past(m_index:(m_index+m_range-1)) = past_m;
mask_past(g_index:(g_index+g_range-1)) = past_g;
mask_past(f_index:(f_index+f_range-1)) = past_f;
mask_past(c_index:(c_index+c_range-1)) = past_c;

mask_fut = false(total_range,1);
mask_fut(i_index:(i_index+i_range-1)) = fut_i;
mask_fut(m_index:(m_index+m_range-1)) = fut_m;
mask_fut(g_index:(g_index+g_range-1)) = fut_g;
mask_fut(f_index:(f_index+f_range-1)) = fut_f;
mask_fut(c_index:(c_index+c_range-1)) = fut_c;

% initialize bounds from incumbent (carry past decisions by default)
if isempty(mip_start), mip_start = zeros(total_range,1); end
lb(:) = mip_start;     % default: fixed to incumbent
ub(:) = mip_start;

% free in-window blocks
lb(mask_in_i) = 0.0;           ub(mask_in_i) = 1.0;          % binaries
lb(mask_in_m) = 0.0;           ub(mask_in_m) = 1.0;          % binaries

% choose safe caps for γ and f (scalar is OK)
U_cap = max(sum(R(:,1:t1),2)) + 1;   % scalar upper bound
lb(mask_in_g) = 0.0;           ub(mask_in_g) = U_cap; %inf;       % g
lb(mask_in_f) = 0.0;           ub(mask_in_f) = inf;       % f

% for ctts, allow growth up to current time
lb(mask_in_c) = 0.0;           ub(mask_in_c) = t1 + 1; %inf      % c 

% hard-off future vars
lb(mask_fut) = 0.0;
ub(mask_fut) = 0.0;

% % (optional) explicitly fix the carry-in ledgers at t0-1 (redundant but clear)
% if t0 > 1
%     fix_g_prev = (tlabel_g == (t0-1));
%     fix_c_prev = (tlabel_c == (t0-1));
%     idx_g_prev = g_index - 1 + find(fix_g_prev);
%     idx_c_prev = c_index - 1 + find(fix_c_prev);
%     lb(idx_g_prev) = mip_start(idx_g_prev);
%     ub(idx_g_prev) = mip_start(idx_g_prev);
%     lb(idx_c_prev) = mip_start(idx_c_prev);
%     ub(idx_c_prev) = mip_start(idx_c_prev);
% end


% -------- objective (only on window slots) ----------
obj = zeros(total_range,1);

obj(g_index:(g_index+g_range-1)) = (weight_traffic_UC / normalization_UC) * double(tlabel_g >= t0 & tlabel_g <= t1);
obj(f_index:(f_index+f_range-1)) = (weight_traffic_EC / normalization_EC) * double(tlabel_f >= t0 & tlabel_f <= t1);
obj(c_index:(c_index+c_range-1)) = (weight_timing     / normalization_time) * double(tlabel_c >= t0 & tlabel_c <= t1);

% -------- Build Gurobi model ----------
model.A          = Awin;          % sparse
model.rhs        = bwin;          % double
model.sense      = sense_win;     % char row
model.obj        = obj;
model.modelsense = 'min';
model.vtype      = vtype;         % char column
model.lb         = lb;
model.ub         = ub;
% model.start      = mip_start;

% Parameters
params.OutputFlag = 1;
params.MIPGap     = MIPGap;
% (you can also set: params.MIPFocus=1; params.Heuristics=0.3; params.Cuts=1;)

% ---- Solve (MATLAB API) ----
% Set path if needed:
if ispc
    addpath('C:\gurobi1100\win64\matlab');
elseif isunix
    addpath('/media/apps/avx512-2021/software/Gurobi/10.0.1-GCCcore-12.2.0/matlab');
end

result = gurobi(model, params);

% Cleanup occasionally: unloads all compiled MEX binaries from memory (Matlab MexFunction crash) 
if (mod(t_current-1, W*10) == 0) && (betta~=0) && (betta~=1) 
    clear mex; clear gurobi;
    fprintf('Gurobi reloaded at iteration %d\n', t_current);
end

% Make outputs consistent with your previous signature
data      = [];     % (no callback trace here)
solutions = [];     % (if you want incumbent history, we can add a callback later)

% Convert to doubles (defensive)
if isfield(result,'x'), result.x = double(result.x); end
end
>>>>>>> e2b063cc1b089e5fb1b6d4fe46a42b06d939b922
=======
function [result, data, solutions] = gurobi_execution_BH_windowed( ...
    betta, normalization_UC, normalization_EC, normalization_time, ...
    total_range, A, b, sense, vtype, MIPGap, ...
    mip_start, lb, ub, ...
    i_index, i_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, ...
    t_row_label, tlabel_i, tlabel_m, tlabel_g, tlabel_f, tlabel_c, ...
    t_current, W, R)

% -------- weights from beta (same logic as before) ----------
weight_traffic_EC = (1 - betta)/2;
weight_traffic_UC = betta;
weight_timing     = (1 - betta)/2;

% -------- window bounds ----------
Hf = max(t_row_label);                    % total inner slots in the built model
t0 = t_current;
t1 = min(t_current + W - 1, Hf);

% -------- select window constraints ----------
rows_keep = (t_row_label >= t0) & (t_row_label <= t1);
Awin  = A(rows_keep, :);
bwin  = b(rows_keep, :);
sense_win = sense(rows_keep);

% -------- decision variable selection per window ----------
%+ past (< t0): fix to incumbent (ub==lb==mip_start)
%+ in window [t0..t1]: free (binaries [0,1], continuous [0,∞))
%+ future (> t1): hard-off (ub==lb==0) — not strictly required since their rows are dropped, but it tightens the model and speeds the solver
past_i   = (tlabel_i <  t0);   in_i = (tlabel_i >= t0 & tlabel_i <= t1);   fut_i = (tlabel_i >  t1);
past_m   = (tlabel_m <  t0);   in_m = (tlabel_m >= t0 & tlabel_m <= t1);   fut_m = (tlabel_m >  t1);
past_g   = (tlabel_g <  t0);   in_g = (tlabel_g >= t0 & tlabel_g <= t1);   fut_g = (tlabel_g >  t1);
past_f   = (tlabel_f <  t0);   in_f = (tlabel_f >= t0 & tlabel_f <= t1);   fut_f = (tlabel_f >  t1);
past_c   = (tlabel_c <  t0);   in_c = (tlabel_c >= t0 & tlabel_c <= t1);   fut_c = (tlabel_c >  t1);

% build full-length logical masks aligned with x for window
mask_in_i = false(total_range,1);  mask_in_i(i_index:(i_index+i_range-1)) = in_i;
mask_in_m = false(total_range,1);  mask_in_m(m_index:(m_index+m_range-1)) = in_m;
mask_in_g = false(total_range,1);  mask_in_g(g_index:(g_index+g_range-1)) = in_g;
mask_in_f = false(total_range,1);  mask_in_f(f_index:(f_index+f_range-1)) = in_f;
mask_in_c = false(total_range,1);  mask_in_c(c_index:(c_index+c_range-1)) = in_c;

mask_past = false(total_range,1);
mask_past(i_index:(i_index+i_range-1)) = past_i;
mask_past(m_index:(m_index+m_range-1)) = past_m;
mask_past(g_index:(g_index+g_range-1)) = past_g;
mask_past(f_index:(f_index+f_range-1)) = past_f;
mask_past(c_index:(c_index+c_range-1)) = past_c;

mask_fut = false(total_range,1);
mask_fut(i_index:(i_index+i_range-1)) = fut_i;
mask_fut(m_index:(m_index+m_range-1)) = fut_m;
mask_fut(g_index:(g_index+g_range-1)) = fut_g;
mask_fut(f_index:(f_index+f_range-1)) = fut_f;
mask_fut(c_index:(c_index+c_range-1)) = fut_c;

% initialize bounds from incumbent (carry past decisions by default)
if isempty(mip_start), mip_start = zeros(total_range,1); end
lb(:) = mip_start;     % default: fixed to incumbent
ub(:) = mip_start;

% free in-window blocks
lb(mask_in_i) = 0.0;           ub(mask_in_i) = 1.0;          % binaries
lb(mask_in_m) = 0.0;           ub(mask_in_m) = 1.0;          % binaries

% choose safe caps for γ and f (scalar is OK)
U_cap = max(sum(R(:,1:t1),2)) + 1;   % scalar upper bound
lb(mask_in_g) = 0.0;           ub(mask_in_g) = U_cap; %inf;       % g
lb(mask_in_f) = 0.0;           ub(mask_in_f) = inf;       % f

% for ctts, allow growth up to current time
lb(mask_in_c) = 0.0;           ub(mask_in_c) = t1 + 1; %inf      % c 

% hard-off future vars
lb(mask_fut) = 0.0;
ub(mask_fut) = 0.0;

% % (optional) explicitly fix the carry-in ledgers at t0-1 (redundant but clear)
% if t0 > 1
%     fix_g_prev = (tlabel_g == (t0-1));
%     fix_c_prev = (tlabel_c == (t0-1));
%     idx_g_prev = g_index - 1 + find(fix_g_prev);
%     idx_c_prev = c_index - 1 + find(fix_c_prev);
%     lb(idx_g_prev) = mip_start(idx_g_prev);
%     ub(idx_g_prev) = mip_start(idx_g_prev);
%     lb(idx_c_prev) = mip_start(idx_c_prev);
%     ub(idx_c_prev) = mip_start(idx_c_prev);
% end


% -------- objective (only on window slots) ----------
obj = zeros(total_range,1);

obj(g_index:(g_index+g_range-1)) = (weight_traffic_UC / normalization_UC) * double(tlabel_g >= t0 & tlabel_g <= t1);
obj(f_index:(f_index+f_range-1)) = (weight_traffic_EC / normalization_EC) * double(tlabel_f >= t0 & tlabel_f <= t1);
obj(c_index:(c_index+c_range-1)) = (weight_timing     / normalization_time) * double(tlabel_c >= t0 & tlabel_c <= t1);

% -------- Build Gurobi model ----------
model.A          = Awin;          % sparse
model.rhs        = bwin;          % double
model.sense      = sense_win;     % char row
model.obj        = obj;
model.modelsense = 'min';
model.vtype      = vtype;         % char column
model.lb         = lb;
model.ub         = ub;
% model.start      = mip_start;

% Parameters
params.OutputFlag = 1;
params.MIPGap     = MIPGap;
% (you can also set: params.MIPFocus=1; params.Heuristics=0.3; params.Cuts=1;)

% ---- Solve (MATLAB API) ----
% Set path if needed:
if ispc
    addpath('C:\gurobi1100\win64\matlab');
elseif isunix
    addpath('/media/apps/avx512-2021/software/Gurobi/10.0.1-GCCcore-12.2.0/matlab');
end

result = gurobi(model, params);

% Cleanup occasionally: unloads all compiled MEX binaries from memory (Matlab MexFunction crash) 
if (mod(t_current-1, W*10) == 0) && (betta~=0) && (betta~=1) 
    clear mex; clear gurobi;
    fprintf('Gurobi reloaded at iteration %d\n', t_current);
end

% Make outputs consistent with your previous signature
data      = [];     % (no callback trace here)
solutions = [];     % (if you want incumbent history, we can add a callback later)

% Convert to doubles (defensive)
if isfield(result,'x'), result.x = double(result.x); end
end
>>>>>>> e2b063cc1b089e5fb1b6d4fe46a42b06d939b922
