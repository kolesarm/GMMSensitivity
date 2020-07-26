% Run this by typing get_ags into matlab's REPL

% 0. Unzip https://doi.org/10.7910/DVN/LLARSN/2KFPRA
% into ags/ directory (the files not needed for code below to execute were deleted)

% 1. Run sections of
%
% "ags/analysis/Transparent Identification (BLP Visualization)/code/sensitivity_to_moments.m"
%
% to retrieve H, G, W, Omega, h_init

addpath(genpath('ags/analysis/Transparent Identification (BLP Visualization)/external'))
load 'ags/analysis/Transparent Identification (BLP Visualization)/external/data/blp_estimation'

data = BlpData('blp_1999_data.csv', 'meanincome.csv', 'sdincome.csv');
data = data.LoadUnobservablesFromEstimate(est);

m_fun_param       = @(param) get_mean_markup(est, data, param);
mjacobian_param   = NumJacob(m_fun_param, est.param, 10^-4);
H                 = [mjacobian_param, zeros(1, length(est.beta))];
h_init            = est.GetMeanMarkup(data); % estimate of markup
G=est.gjacobian;
W=est.wmatrix;

Omega=est.Omega;
names_iv=est.model.iv_varnames;
names_th=[est.model.paramlist, est.model.beta_paramlist];
g_init = est.g;

% from function get_IV_Sensitivity in sensitivity_to_moments.m
demand_iv_var = data.GetArray(est.model.demand_iv_varlist);
supply_iv_var = data.GetArray(est.model.supply_iv_varlist);
z = blkdiag(demand_iv_var, supply_iv_var);
Om_ZZ = (z' * z) ./ est.nmodels;
sd_Z = sqrt(var([demand_iv_var, supply_iv_var]));

% 2. Run section of
%
% "ags/analysis/Transparent Identification (BLP Visualization)/code/misspecification.m"
%
% to calculate perturbation scaling from
[avg_price, avg_mc] = est.GetAvgPriceMC(data);
k_xi                = est.GetWTPDerivative(data);
demand_perturb =  0.01 .* avg_price ./ k_xi;
supply_perturb = -0.01 .* avg_price ./ avg_mc;

n = est.nmodels;

save('agm_data', 'Om_ZZ', 'sd_Z', 'W', 'H', 'G', 'Omega', 'demand_perturb', ...
     'supply_perturb', 'g_init', 'h_init', 'names_iv', 'names_th', 'n');

% Print thattheta (se)
[[est.param; est.beta], sqrt(diag(get_vcov(est.gjacobian, est.Omega, est.wmatrix))/n)]
% Print estimate normalized in terms of % willingness to pay for 1sd increase
(est.beta'./[demand_perturb./sd_Z(1:5), supply_perturb./sd_Z(14:19)])'

% Markup
[h_init, sqrt(H * get_vcov(est.gjacobian, est.Omega, est.wmatrix) * H' ./n)]

% Copied from "sensitivity_to_moments.m"
function mean_markup = get_mean_markup(est, data, param)
    price = data.GetArray(data.varlist.price);
    [~, mc] = est.model.ComputeModelOutputs(data, param);
    markup = (price - mc) ./ price;
    mean_markup = mean(markup);
end

% Relative to Soonwoo's data directory:
% Z.csv -> sdZ = sd_Z
% ZZmean.csv = Om_ZZ
% g.csv -> g_init
% gjacob.csv = G
% ivnames.csv = ivnames
% mjacob.csv = H
% weights.csv = W
% Omega is new
