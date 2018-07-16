function numb_cars = get_cars_rival_firm(obj)
%
% Calculates number of rival cars produced each firm-year
%

    numb_cars = zeros([length(obj.var(:, 1)), 1]);
    firms = unique(obj.var.firm_id);
    years = unique(obj.var.year);
    for i = 1:length(firms)    % over firms
        for j = 1:length(years)    % over years
            numb_cars(obj.var.firm_id == firms(i) & obj.var.year == years(j), 1) = ...
                 length(obj.var(obj.var.firm_id ~= firms(i) & obj.var.year == years(j), 1));
        end
    end
end
