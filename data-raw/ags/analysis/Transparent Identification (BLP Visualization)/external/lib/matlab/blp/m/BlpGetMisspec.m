function [misspec_perturb] = BlpGetMisspec(demand_supply, same_rival, delta,...
                                           avg_price, avg_mc, k_xi, ...
                                           M_same, M_rival)

% This helper function calculate the degree of misspecification for the BLP Model,
% needed to calculate either first order asymptotic bias or re-estimate alternatives.

    % Check symmetry
    if ~(isvector(M_same) && isvector(M_rival) && numel(M_same) == numel(M_rival))
        error('The vectors M_same and M_rival have different dimensions!');
    end
    n = length(M_same);
    zero_vector = zeros([n, 1]);

    % Get number of cars for same/rival firm 
    if strcmp(same_rival, 'same')
        numb_cars = M_same;
    end
    if strcmp(same_rival, 'rival')
        numb_cars = M_rival;
    end

    % Calculate misspecification
    if strcmp(demand_supply, 'supply')
        supply_perturb = (-delta .* (avg_price ./ avg_mc)) .* numb_cars;
        misspec_perturb = [zero_vector; supply_perturb];
    end
    if strcmp(demand_supply, 'demand')
        demand_perturb = ( delta .* (avg_price ./ k_xi))   .* numb_cars;
        misspec_perturb = [demand_perturb; zero_vector];
    end
end