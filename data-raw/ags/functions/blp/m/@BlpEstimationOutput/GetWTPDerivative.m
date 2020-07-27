function [k_xi] = GetWTPDerivative(est, data)
%  Calculate the derivative of willingness to pay with 
%  respect to xi for a 1980 household with mean income

    % Calculating k_xi 
    %   get scaling parameters for interpreting magnitudes (base year = tstar = 1980)
    m_tstar = est.log_income_mean(1971:1990 == 1980);
    sigma_y = est.log_income_sd;
    y_bar_tstar = exp(m_tstar + sigma_y^2 / 2);
    %   gets parameter value for 'alpha_price'
    alpha = est.param(find(~cellfun('isempty', strfind(est.model.paramlist, 'alpha_price')))); 
    k_xi = y_bar_tstar / alpha;

end