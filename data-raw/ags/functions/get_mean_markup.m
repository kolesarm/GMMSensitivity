% Copied from "sensitivity_to_moments.m"
function mean_markup = get_mean_markup(est, data, param)
    price = data.GetArray(data.varlist.price);
    [~, mc] = est.model.ComputeModelOutputs(data, param);
    markup = (price - mc) ./ price;
    mean_markup = mean(markup);
end
