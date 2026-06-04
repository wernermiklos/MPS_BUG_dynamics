function Model = Set_Model(varargin)
%SET_MODEL Summary of this function goes here
%   Detailed explanation goes here
  try
    ModelParams = struct(varargin{:});
  catch
    error('Input must be given in Name-Value pairs, or as a struct')
  end
  if ~isfield(ModelParams,'ModelInput')
    error('aaa')
  end
  ModelFunction = str2func(ModelParams.ModelInput);
  Model = ModelFunction(ModelParams);
end

