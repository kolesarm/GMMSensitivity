% Run this by typing get_ags into matlab's REPL

% 0. Unzip https://doi.org/10.7910/DVN/LLARSN/2KFPRA
% into ags/ directory (the files not needed for code below to execute were deleted)

% 1. Run sections of
%
% "ags/analysis/Transparent Identification (BLP Visualization)/code/sensitivity_to_moments.m"
%
% to retrieve H, G, W, Omega, h_init

% Move the data needed to data/ and the functions to functions/
addpath(genpath('ags/functions'))
addpath('ags/data')
load 'ags/data/blp_estimation'

data = BlpData('blp_1999_data.csv', 'meanincome.csv', 'sdincome.csv');
data = data.LoadUnobservablesFromEstimate(est);

G = est.gjacobian;
W = est.wmatrix;
Omega=est.Omega; % max(max(abs(cov(est.g_model)-est.Omega))) = 1e-11
names_iv=est.model.iv_varnames;
% paramlist: sigma's and alpha; beta: betabar and gamma;
names_th=[est.model.paramlist, est.model.beta_paramlist];
% This is average of moment condition, mean(est.g_model)-est.g' is 0
g_init = est.g;

% Print thattheta (se)
[[est.param; est.beta], sqrt(diag(get_vcov(est.gjacobian, est.Omega, est.wmatrix))/n)]
%   2.52188717037472          3.77917192869453
%   3.52454664272272          4.23608134396864
%   4.16663786353709          2.10594354664105
%   0.39290543029902         0.419166781167063
%   1.93661212235348          0.88853292704118
%    42.870219350619          8.27984710282358
%  -7.72836850083045          1.72230726521818
%   4.62042079204929          1.68180934532446
%  -1.22659416183146          2.05927031944406
%  0.293176102249026         0.233429754679872
%   3.99186552658677         0.526861969488043
%   2.75052070150035         0.125121812242542
%  0.812222486417721        0.0894122902461129
%  0.430134306782569        0.0790501959720123
% -0.610079037870289        0.0725284356317645
% -0.352031684256758         0.163635519972971
% 0.0268821385305752       0.00231182747003052

% Markup
m_fun_param       = @(param) get_mean_markup(est, data, param);
H                 = [NumJacob(m_fun_param, est.param, 10^-4), zeros(1, length(est.beta))];
h_init            = est.GetMeanMarkup(data); % estimate of markup

[h_init, sqrt(H * get_vcov(est.gjacobian, est.Omega, est.wmatrix) * H' ./n)]
% 0.327178898099534        0.0181566468728511

% from function get_IV_Sensitivity in sensitivity_to_moments.m
demand_iv = data.GetArray(est.model.demand_iv_varlist); % 2217 by 13
supply_iv = data.GetArray(est.model.supply_iv_varlist); % 2217 by 18
z = blkdiag(demand_iv, supply_iv);
Om_ZZ = (z' * z) ./ est.nmodels;
sd_Z = sqrt(var([demand_iv, supply_iv]));

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

% Print estimate normalized in terms of % willingness to pay for 1sd increase
(est.beta'./[demand_perturb./sd_Z(1:5), supply_perturb./sd_Z(14:19)])'

% log of output -> Table 1 in blp.pdf

% Relative to Soonwoo's data directory:
% Z.csv -> sdZ = sd_Z
% ZZmean.csv = Om_ZZ
% g.csv -> g_init
% gjacob.csv = G
% ivnames.csv = ivnames
% mjacob.csv = H
% weights.csv = W
% Omega is new
