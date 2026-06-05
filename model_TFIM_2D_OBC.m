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
  Model.ModelParams = ModelParams;
  Model = internal_set_MPO_TFIM_2D_OBC(Model);
end

function Model = internal_set_MPO_TFIM_2D_OBC(Model)
  S = Model.ModelParams.S;
  [Sx,~,Sz] = SpinOperators(S);
  d = 2*S+1;
  Model.CL = Model.ModelParams.Lx*Model.ModelParams.Ly;
  Lx = Model.ModelParams.Lx;
  Ly = Model.ModelParams.Ly;
  CL = Model.CL;
  h = Model.ModelParams.h;
  J = 1;
  % first column of sites
  x = 1; y = 1;
  pos = 1;  %position index within MPS chain
  Model.MPO{pos} = cell(1,y+2);
  Model.MPO{pos}{1} = h*Sx;
  Model.MPO{pos}{2} = eye(d);
  Model.MPO{pos}{3} = Sz;
  for y = 2:Ly
    pos = pos + 1;
    Model.MPO{pos} = cell(y+1,y+2);
    Model.MPO{pos}{1,1} = eye(d);
    Model.MPO{pos}{2,1} = h*Sx;
    Model.MPO{pos}{end,1} = J*Sz;
    Model.MPO{pos}{2,2} = eye(d);
    for i=3:(y+1)
      Model.MPO{pos}{i,i} = eye(d);
    end
    Model.MPO{pos}{2,end} = Sz;
  end
  % from second column till the end  (the last column will be modified
  % afterwards by erasing unnecessary rows and columns)
  for x = 2:Lx
    y = 1;
    pos = pos + 1;
    Model.MPO{pos} = cell(Ly+2,Ly+2);
    Model.MPO{pos}{1,1} = eye(d);
    Model.MPO{pos}{2,1} = h*Sx;
    Model.MPO{pos}{3,1} = J*Sz;
    Model.MPO{pos}{2,2} = eye(d);
    for i=3:(Ly+1)
      Model.MPO{pos}{i+1,i} = eye(d);
    end
    Model.MPO{pos}{2,Ly+2} = Sz;
    for y = 2:Ly
      pos = pos + 1;
      Model.MPO{pos} = cell(Ly+2,Ly+2);
      Model.MPO{pos}{1,1} = eye(d);
      Model.MPO{pos}{2,1} = h*Sx;
      Model.MPO{pos}{3,1} = J*Sz;
      Model.MPO{pos}{Ly+2,1} = J*Sz;
      Model.MPO{pos}{2,2} = eye(d);
      for i=3:(Ly+1)
        Model.MPO{pos}{i+1,i} = eye(d);
      end
      Model.MPO{pos}{2,Ly+2} = Sz;
    end
  end
  %Erase unnecessary rows/columns from MPO at the last column of sites
  pos = pos - Ly + 1;   %x = Lx, y = 1 here we do not need to erase anything
  for y = 2:Ly
    pos = pos + 1;
    Model.MPO{pos}(:,end-(y-1):end-1) = [];
    Model.MPO{pos}(end-(y-2):end-1,:) = [];
  end
  Model.MPO{pos}(:,2:end) = [];  % at the last site, only the Block Hamiltonian is needed.
  for pos = 1:CL
    Model.MPOlogical{pos} = cellfun(@(x) ~isempty(x),Model.MPO{pos});
    [Model.MPOtasktable{pos}(:,1), Model.MPOtasktable{pos}(:,2)] = find(Model.MPOlogical{pos});
  end







end

