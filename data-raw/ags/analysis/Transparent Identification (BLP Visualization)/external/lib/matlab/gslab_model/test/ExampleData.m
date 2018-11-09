classdef ExampleData < ModelData

    methods
        function obj = ExampleData(varargin)
            % Extract the input variable names
            variable_names = cell(size(varargin));
            for j = 1:nargin
                variable_names{j} = inputname(j);
            end
        
            inputlist = ModelData.ParseInputList(varargin, variable_names);
            obj.var = dataset(inputlist{:});
        end    
    end

end