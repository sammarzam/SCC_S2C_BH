%% SOLUTION PLOTTING:
function [max_UC_betta_0, max_EC_betta_1, max_TTS_betta_1] = solution_plot_saving(save_dir, scenario, use_case, h3_resolution, beams, P_T, m_continuous, betta, normUC, normEC, normTTS, data_trace, x_incumbent, i_index,i_range,m_index,m_range,g_index,g_range,f_index,f_range,c_index,c_range, tlabel_g, tlabel_f, tlabel_c, t0, t1, number_cells, H, frame, M, X, D, nSats)

% save_dir: folder to write .mat (string)
% data_trace: optional solver trace ([], or your callback data)
% x_incumbent: result.x (vector, doubles)

if nargin<1 || isempty(save_dir), save_dir = pwd; end
if ~isfolder(save_dir), mkdir(save_dir); end

% --- compute raw objectives on the current window ---
% ----- variable slices -----
i_all = x_incumbent(i_index:(i_index+i_range-1));
m_all = x_incumbent(m_index:(m_index+m_range-1));
g_all = x_incumbent(g_index:(g_index+g_range-1));
f_all = x_incumbent(f_index:(f_index+f_range-1));
c_all = x_incumbent(c_index:(c_index+c_range-1));

% ----- window masks -----
mask_g = (tlabel_g >= t0) & (tlabel_g <= t1);
mask_f = (tlabel_f >= t0) & (tlabel_f <= t1);
mask_c = (tlabel_c >= t0) & (tlabel_c <= t1);

% ----- objective components (sums over window) -----
sum_UC  = sum( g_all(mask_g) );   % O1
O1=sum_UC;
sum_EC  = sum( f_all(mask_f) );   % O2
O2=sum_EC;
sum_TTS = sum( c_all(mask_c) );   % O3
O3=sum_TTS;

% --- normalized combined objective (for logging only) ---
w_EC  = (1-betta)/2;
w_UC  = betta;
w_TTS = (1-betta)/2;
if normUC<=0,  normUC=1;  end
if normEC<=0,  normEC=1;  end
if normTTS<=0, normTTS=1; end
Obj_norm = (w_UC*O1)/normUC + (w_EC*O2)/normEC + (w_TTS*O3)/normTTS;

% ---- Reshape back to their natural tensor shapes ----
Hf=H*frame;
i_val = reshape(i_all, [number_cells, Hf]);          % i(c,t)
m_val = reshape(m_all, [number_cells, M, Hf]);       % m(c,mcs,t)
g_val = reshape(g_all, [number_cells, Hf]);      % γ(c,t)
f_val = reshape(f_all, [number_cells, Hf]);          % f(c,t)
c_val = reshape(c_all, [number_cells, Hf]);       % c(c,t)

% ---- Rebuild d_{c,t} from m and D under Ω-compact indexing ----
% Expand responsibility to inner slots and pick responsible sat per (c,t)
X_exp = repelem(X, 1, 1, frame);                      % [cells x nSats x Hf]
[~, s_star_mat] = max(X_exp, [], 2);                  % [cells x 1 x Hf]
s_star_mat = squeeze(s_star_mat);                     % [cells x Hf]
% Repeat D from H to Hf for each satellite
D_exp = D;                                           % 1 x nSats cell
for s = 1:nSats
   % D{s}: [cells x M x H] -> [cells x M x Hf]
   D_exp{s} = repelem(D{s}, 1, 1, frame);
end
% Compute d(c,t) = sum_m D_{s*(c,t)}(c,m,t) * m(c,m,t)
d_val = zeros(number_cells, Hf);
for s = 1:nSats
   mask_s = (s_star_mat == s);                      % logical [cells x Hf]
   if any(mask_s(:))
       % contribution for this satellite: sum over MCS
       contrib_s = squeeze(sum(m_val .* D_exp{s}, 2));  % [cells x Hf]
       d_val(mask_s) = contrib_s(mask_s);
   end
end

% --- file name & payload ---
fname = sprintf('%s_BH_[%s_res%d_beams%d]_P_%d_mC_%d_beta_%0.2f_win_%d_%d.mat',scenario,use_case, h3_resolution,beams, P_T, m_continuous, betta, t0, t1);


payload = struct();
payload.beta          = betta;
payload.weights       = struct('w_UC',w_UC,'w_EC',w_EC,'w_TTS',w_TTS);
payload.norm          = struct('UC',normUC,'EC',normEC,'TTS',normTTS);
payload.window        = struct('t0',t0,'t1',t1);
payload.objectives    = struct('UC',O1,'EC',O2,'TTS',O3,'Obj_norm',Obj_norm);
payload.trace         = data_trace;
payload.x             = x_incumbent;
payload.i             = i_val;
payload.m             = m_val;
payload.g             = g_val;
payload.f             = f_val;
payload.c             = c_val;

% derived served demand
payload.d = d_val;     % d(c,t) reconstructed

% ONLY SAVE WHEN FINAL FILE!
if t1==Hf
    save(fullfile(save_dir, fname), '-struct', 'payload');
end

% --- return values to the caller (for normalization passes) ---
% Caller decides how to use them depending on beta:
%   beta==0: use max_UC_betta_0 = O1
%   beta==1: use max_EC_betta_1 = O2 and max_TTS_betta_1 = O3
max_UC_betta_0  = O1;
max_EC_betta_1  = O2;
max_TTS_betta_1 = O3;

end

