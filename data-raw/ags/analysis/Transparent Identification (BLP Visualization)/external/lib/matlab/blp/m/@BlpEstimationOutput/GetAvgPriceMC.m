function [avg_price, avg_mc] = GetAvgPriceMC(obj, data)
% INPUTS
%       - BLPData object

        param = obj.param;
        model = obj.model;
    [~, mc] = model.ComputeModelOutputs(data, param);

    sales_weights  = data.var.quantity(data.var.year == 80) ./ sum(data.var.quantity(data.var.year == 80));
    avg_price      = sum(data.var.price(data.var.year == 80) .* sales_weights);
    avg_mc         = sum(mc(data.var.year == 80) .* sales_weights);

end
