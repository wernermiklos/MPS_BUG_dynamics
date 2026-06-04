function [BUGrun,Model,MPS,Blocks] = Setup_BUG_withMPO(BUGparams)
  arguments
    BUGparams.dt (1,1) double = 0.01;
    BUGparams.tmax (1,1) double = 0.01;
    BUGparams.ModelParams (1,1) struct;
    BUGparams.InitMPS struct;
  end

  BUGrun.t = 0;
  BUGrun.dt = BUGparams.dt;
  BUGrun.tmax = BUGparams.tmax;
  
  Model = Set_Model(BUGparams.ModelParams);
  
  nvpairs = namedargs2cell(BUGparams.InitMPS);
  MPS = Set_MPS(nvpairs{:});

  if MPS.CL ~= Model.CL
    error('ajjaj');
  end
  Blocks = struct();
  Blocks.Left = cell(1,length(MPS.LeftMatrices)+1);
  Blocks.Right = cell(1,length(MPS.RightMatrices)+1);
  Blocks.Left{1} = struct();
  Blocks.Left{1}.OP = {[1]};
  Blocks.Left{1}.OP_logical = true;
  Blocks.Left{1}.side = 'left';
  for pos = 1:length(MPS.LeftMatrices)
    Blocks.Left{pos+1} = renorm(Blocks.Left{pos},Model.MPO{pos},Model.MPOlogical{pos},MPS.LeftMatrices{pos},'SIDE','left');
  end
  Blocks.Right{1} = struct();
  Blocks.Right{1}.OP = {[1]};
  Blocks.Right{1}.OP_logical = true;
  Blocks.Right{1}.side = 'right';
  for rightpos = 1:length(MPS.RightMatrices)
    MPOpos = Model.CL-rightpos+1;
    Blocks.Right{rightpos+1} = renorm(Blocks.Right{rightpos},Model.MPO{MPOpos},...
                                      Model.MPOlogical{MPOpos},MPS.RightMatrices{rightpos},'SIDE','right');
  end
end