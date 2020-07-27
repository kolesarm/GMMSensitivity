function mean_markup = GetMeanMarkup(obj, data, param)
	if nargin <= 2
		param = obj.param;
	end
    
    [~, mc] = obj.model.ComputeModelOutputs(data, param);
    price = data.GetArray(data.varlist.price);
    markup = (price - mc) ./ price;
    mean_markup = mean(markup);
end