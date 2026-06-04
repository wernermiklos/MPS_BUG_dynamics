function Model = model_TFIM_2D_OBC(varargin)
%MODEL_TFIM_2D Summary of this function goes here
%   Detailed explanation goes here
  if nargin == 1
    ModelParams = varargin{1};
  else
    try
      ModelParams = struct(varargin{:});
    catch
      error('name value pairs')
    end
  end
  Model.Lx = ModelParams.Lx;
  Model.Ly = ModelParams.Ly;
  Model.h = ModelParams.h;
  Model.ModelInput = ModelParams.ModelInput;
  Model = internal_set_MPO_TFIM_2D_OBC(Model);
end

function Model = internal_set_MPO_TFIM_2D_OBC(Model)


end

