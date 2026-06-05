function Model = model_TFIM_1D_OBC(varargin)
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
  Model.ModelParams = ModelParams;
  Model = internal_set_MPO_TFIM_1D_OBC(Model);
end

function Model = internal_set_MPO_TFIM_1D_OBC(Model)
  [Sx,~,Sz] = SpinOperators(Model.ModelParams.S);
  d = 2*Model.ModelParams.S + 1;
  Model.CL = Model.ModelParams.L;
  h = Model.ModelParams.h;
  Model.MPO = cell(1,Model.CL);
  Model.MPOlogical = cell(1,Model.CL);
  Model.MPOtasktable = cell(1,Model.CL);
  Model.MPO{1} = cell(1,3);
  Model.MPO{1}{1} = eye(d);
  Model.MPO{1}{2} = Sz;
  Model.MPO{1}{3} = h*Sx;
  Model.MPOlogical{1} = cellfun(@(x) ~isempty(x),Model.MPO{1});
  [Model.MPOtasktable{1}(:,1), Model.MPOtasktable{1}(:,2)] = find(Model.MPOlogical{1});
  Model.MPO{end} = cell(3,1);
  Model.MPO{end}{1} = h*Sx;
  Model.MPO{end}{2} = Sz;
  Model.MPO{end}{3} = eye(d);
  Model.MPOlogical{end} = cellfun(@(x) ~isempty(x),Model.MPO{end});
  [Model.MPOtasktable{end}(:,1), Model.MPOtasktable{end}(:,2)] = find(Model.MPOlogical{end});
  for pos = 2:Model.CL-1
    Model.MPO{pos} = cell(3,3);
    Model.MPO{pos}{1,1} = eye(d);
    Model.MPO{pos}{1,2} = Sz;
    Model.MPO{pos}{1,3} = h*Sx;
    Model.MPO{pos}{2,3} = Sz;
    Model.MPO{pos}{3,3} = eye(d);
    Model.MPOlogical{pos} = cellfun(@(x) ~isempty(x),Model.MPO{pos});
    [Model.MPOtasktable{pos}(:,1), Model.MPOtasktable{pos}(:,2)] = find(Model.MPOlogical{pos});
  end
end

