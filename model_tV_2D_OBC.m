function Model = model_tV_2D_OBC(varargin)
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
  Model = internal_set_MPO_tV_2D_OBC(Model);
end

function Model = internal_set_MPO_tV_2D_OBC(Model)
  myOP = spinless_fermi_operators();
  d = 2;
  Model.CL = Model.ModelParams.Lx*Model.ModelParams.Ly;
  Lx = Model.ModelParams.Lx;
  Ly = Model.ModelParams.Ly;
  V = Model.ModelParams.V;
  J = -1;
  % first column of sites
  x = 1; y = 1;
  pos = 1;  %position index within MPS chain
  Model.MPO{pos} = cell(1,3*y+2);
  Model.MPO{pos}{2} = myOP.id;
  Model.MPO{pos}{3} = myOP.cdag_ph;
  Model.MPO{pos}{4} = myOP.c_ph;
  Model.MPO{pos}{5} = myOP.n;
  for y = 2:Ly
    pos = pos + 1;
    Model.MPO{pos} = cell(2+3*(y-1),2+3*y);
    Model.MPO{pos}{1,1} = eye(d);
    Model.MPO{pos}{2+(y-1),1} = J*myOP.c;
    Model.MPO{pos}{2+2*(y-1),1} = -J*myOP.cdag;
    Model.MPO{pos}{2+3*(y-1),1} = V*myOP.n;
    Model.MPO{pos}{2,2} = eye(d);
    for i=1:(y-1);
      Model.MPO{pos}{2+i,2+i} = myOP.ph;
      Model.MPO{pos}{2+(y-1)+i,2+y+i} = myOP.ph;
      Model.MPO{pos}{2+2*(y-1)+i,2+2*y+i} = myOP.id;
    end
    Model.MPO{pos}{2,2+y} = myOP.cdag_ph;
    Model.MPO{pos}{2,2+2*y} = myOP.c_ph;
    Model.MPO{pos}{2,2+3*y} = myOP.n;
  end
  % from second column till the end  (the last column will be modified
  % afterwards by erasing unnecessary rows and columns)
  for x = 2:Lx
    y = 1;
    pos = pos + 1;
    Model.MPO{pos} = cell(2+3*Ly,2+3*Ly);
    Model.MPO{pos}{1,1} = eye(d);
    Model.MPO{pos}{3,1} = J*myOP.c;
    Model.MPO{pos}{3+Ly,1} = -J*myOP.cdag;
    Model.MPO{pos}{3+2*Ly,1} = V*myOP.n;
    Model.MPO{pos}{2,2} = eye(d);
    for i=1:(Ly-1)
      Model.MPO{pos}{2+i+1,2+i} = myOP.ph;
      Model.MPO{pos}{2+Ly+i+1,2+Ly+i} = myOP.ph;
      Model.MPO{pos}{2+2*Ly+i+1,2+2*Ly+i} = myOP.id;
    end
    Model.MPO{pos}{2,2+Ly} = myOP.cdag_ph;
    Model.MPO{pos}{2,2+2*Ly} = myOP.c_ph;
    Model.MPO{pos}{2,2+3*Ly} = myOP.n;
    for y = 2:Ly
      pos = pos + 1;
      Model.MPO{pos} = cell(2+3*Ly,2+3*Ly);
      Model.MPO{pos}{1,1} = eye(d);
      Model.MPO{pos}{3,1} = J*myOP.c;
      Model.MPO{pos}{2+Ly,1} = J*myOP.c;
      Model.MPO{pos}{3+Ly,1} = -J*myOP.cdag;
      Model.MPO{pos}{2+2*Ly,1} = -J*myOP.cdag;
      Model.MPO{pos}{3+2*Ly,1} = V*myOP.n;
      Model.MPO{pos}{2+3*Ly,1} = V*myOP.n;
      Model.MPO{pos}{2,2} = eye(d);
      for i=1:(Ly-1)
        Model.MPO{pos}{2+i+1,2+i} = myOP.ph;
        Model.MPO{pos}{2+Ly+i+1,2+Ly+i} = myOP.ph;
        Model.MPO{pos}{2+2*Ly+i+1,2+2*Ly+i} = myOP.id;
      end
      Model.MPO{pos}{2,2+Ly} = myOP.cdag_ph;
      Model.MPO{pos}{2,2+2*Ly} = myOP.c_ph;
      Model.MPO{pos}{2,2+3*Ly} = myOP.n;
    end
  end
  %Erase unnecessary rows/columns from MPO at the last column of sites
  Model.MPO{pos}(:,2:end) = [];  % at the last site, only the Block Hamiltonian is needed.
  Model.MPO{end}(:,2:end) = [];  % at the last site, only the Block Hamiltonian is needed.
  Model = MPO_setupLogical(Model);
end

