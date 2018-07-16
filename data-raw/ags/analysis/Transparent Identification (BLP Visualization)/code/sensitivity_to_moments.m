function sensitivity_to_moments
    addpath(genpath('../external'))
    load ../external/data/blp_estimation

    blp1995.momlist = est.model.iv_varnames;
    blp1995.paramlist = [est.model.paramlist, est.model.beta_paramlist];

    % Get objects for computing senstivity from the estimation output
    blp1995.weight      = est.wmatrix;
    blp1995.jacobian    = est.gjacobian;
    blp1995.ahat        = est.ahat;
    blp1995.moment_vcov = est.Omega;

    blp1995.se_moments = sqrt(diag(blp1995.moment_vcov));
    blp1995.vcov_param = get_vcov(blp1995.jacobian, blp1995.moment_vcov, blp1995.weight);
    blp1995.se = sqrt(diag(blp1995.vcov_param));

    [blp1995.sensitivity, blp1995.standardized_sensitivity] = ...
        get_sensitivity(blp1995.jacobian, blp1995.weight, blp1995.se, blp1995.se_moments);

    [blp1995.sample_sensitivity, blp1995.standardized_sample_sensitivity] = ...
        get_sample_sensitivity(blp1995.jacobian, blp1995.ahat, ...
                               blp1995.weight, blp1995.se, blp1995.se_moments);

    write_text_table(blp1995.sensitivity', '../temp/sensitivity_matrix.tsv', ...
                     blp1995.momlist, blp1995.paramlist, '', 12);
    write_text_table(abs(blp1995.standardized_sensitivity'), ...
                     '../output/standardized_sensitivity_matrix.tsv', ...
                     blp1995.momlist, blp1995.paramlist, '', 12);

    write_text_table(blp1995.sample_sensitivity', '../temp/sample_sensitivity_matrix.tsv', ...
                     blp1995.momlist, blp1995.paramlist, '', 12);
    write_text_table(abs(blp1995.standardized_sample_sensitivity'), ...
                     '../output/standardized_sample_sensitivity_matrix.tsv', ...
                     blp1995.momlist, blp1995.paramlist, '', 12);

    % Supply and demand sensitivities
    supply = 1 - abs(cellfun(@isempty,strfind(blp1995.momlist, 'supply')));
    supply = logical(supply);

    demand = 1 - abs(cellfun(@isempty,strfind(blp1995.momlist, 'demand')));
    demand = logical(demand);

    supplymom = blp1995.standardized_sensitivity(:, supply);
    demandmom = blp1995.standardized_sensitivity(:, demand);

    abs_demandmom = abs(demandmom);
    abs_supplymom = abs(supplymom);

    absmean_supplydemand_sens = zeros(length(blp1995.paramlist), 2);
    absmean_supplydemand_sens(:, 1) = mean(abs_demandmom, 2);
    absmean_supplydemand_sens(:, 2) = mean(abs_supplymom, 2);

    save('../temp/supplydemand_sensitivities.mat', 'absmean_supplydemand_sens');

    % Load data
    data = BlpData('blp_1999_data.csv', 'meanincome.csv', 'sdincome.csv');
    data = data.LoadUnobservablesFromEstimate(est);

    % Compute standard deviation of moments
    std_momlist = blp1995.momlist;
    std_momlist(1:5)   = {'const', 'hpwt',    'air', 'mpd',    'space'};
    std_momlist(14:19) = {'const', 'loghpwt', 'air', 'logmpg', 'logspace', 'trend'};
    std_momlist(31)    = {'mpd'};
    std_momvec = zeros(length(std_momlist), 1);
    std_mom_var = {};
    for i = 1:length(std_momlist)
        varname        = sprintf('data.var.demeaned_%s', std_momlist{i});
        std_mom_var{i} = sprintf('sd_%s', blp1995.momlist{i});
        std_momvec(i)  = std(eval(varname));
    end
    write_text_table(std_momvec', '../temp/moments_standard_deviation_matrix.tsv', ...
                     'variable', std_mom_var', '', 12);
    rownames = repelem({''}, length(std_momlist));

    % Print excluded instrument standard deviation and list
    excluded_instruments = [6:13 20:31];
    std_excluded = std_momvec(excluded_instruments);
    write_text_table(std_excluded, ...
                    '../output/excluded_instrument_standard_deviation_matrix.tsv', ...
                     rownames(excluded_instruments), '<tab:instrument_std_dev>', '', 12);
    write_text_table(std_excluded, '../output/excluded_instrument_list.tsv', ...
                     std_momlist(excluded_instruments), '<tab:instrument_list>', '', 12);
    save('../temp/excluded_instrument_standard_deviations.mat', 'std_excluded');

    % Sensitivities of mean markup to moments
    m_fun_param       = @(param)get_mean_markup(est, data, param);
    mjacobian_param   = NumJacob(m_fun_param, est.param, 10^-4);
    mjacobian_beta    = zeros(1, length(est.beta));
    blp1995.mjacobian = [mjacobian_param, mjacobian_beta];

    [blp1995.m_sensitivity, blp1995.m_std_sensitivity] = ...
        get_counterfactual_sensitivity(blp1995.sensitivity, ...
            blp1995.mjacobian, blp1995.vcov_param, blp1995.moment_vcov);
    write_text_table(blp1995.m_sensitivity', '../temp/markup_sensitivity_matrix.tsv', ...
                     blp1995.momlist, 'markup', '', 12);
    write_text_table(abs(blp1995.m_std_sensitivity'), ...
                     '../output/standardized_markup_sensitivity_matrix.tsv', ...
                     blp1995.momlist, 'markup', '', 12);

    % Sample sensitivities of mean markup to moments
    [blp1995.m_sample_sensitivity, blp1995.m_std_sample_sensitivity] = ...
        get_counterfactual_sensitivity(blp1995.sample_sensitivity, ...
            blp1995.mjacobian, blp1995.vcov_param, blp1995.moment_vcov);
    write_text_table(blp1995.m_sample_sensitivity', '../temp/markup_sample_sensitivity_matrix.tsv', ...
                     blp1995.momlist, 'markup', '', 12);
    write_text_table(abs(blp1995.m_std_sample_sensitivity'), ...
                     '../output/standardized_markup_sample_sensitivity_matrix.tsv', ...
                     blp1995.momlist, 'markup', '', 12);

    % IV Sensitivity
    iv = struct();
    iv.sensitivity            = get_IV_Sensitivity(est, data, blp1995.m_sensitivity);
    iv.std_sensitivity        = get_IV_Sensitivity(est, data, blp1995.m_std_sensitivity);
    iv.sample_sensitivity     = get_IV_Sensitivity(est, data, blp1995.m_sample_sensitivity);
    iv.std_sample_sensitivity = get_IV_Sensitivity(est, data, blp1995.m_std_sample_sensitivity);
    write_text_table(iv.sensitivity', '../temp/markup_iv_sensitivity_matrix.tsv', ...
                     blp1995.momlist, 'markup', '', 12);
    write_text_table(iv.std_sensitivity', '../output/standardized_markup_iv_sensitivity_matrix.tsv', ...
                     blp1995.momlist, 'markup', '', 12);
    write_text_table(iv.sample_sensitivity', '../temp/markup_iv_sample_sensitivity_matrix.tsv', ...
                     blp1995.momlist, 'markup', '', 12);
    write_text_table(iv.std_sample_sensitivity', '../output/standardized_markup_iv_sample_sensitivity_matrix.tsv', ...
                     blp1995.momlist, 'markup', '', 12);
    save('../temp/markup_iv_sensitivity.mat', 'iv');

    % Compute estimation moments for constructing influence components
    demand_iv_var = data.GetArray(est.model.demand_iv_varlist);
    supply_iv_var = data.GetArray(est.model.supply_iv_varlist);
    z = blkdiag(demand_iv_var, supply_iv_var);
    proj_z = z * est.wmatrix * z';
    [~, g_model] = est.model.ComputeDistanceVector(data, est.param, proj_z);
    blp1995.g = g_model';
    save('../output/transp_identification_blp.mat', 'blp1995')

    exit
end

function mean_markup = get_mean_markup(est, data, param)
    price = data.GetArray(data.varlist.price);
    [~, mc] = est.model.ComputeModelOutputs(data, param);
    markup = (price - mc) ./ price;
    mean_markup = mean(markup);
end

function [sensitivity_iv] = get_IV_Sensitivity(est, data, sensitivity)
    n             = est.nmodels;
    demand_iv_var = data.GetArray(est.model.demand_iv_varlist);
    supply_iv_var = data.GetArray(est.model.supply_iv_varlist);
    z             = blkdiag(demand_iv_var, supply_iv_var);

    sensitivity_iv  = sensitivity * (z' * z) ./ n;
end
