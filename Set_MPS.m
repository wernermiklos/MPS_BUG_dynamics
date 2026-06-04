function MPS = Set_MPS(MPSparams)
%SET_MPS Summary of this function goes here
%   Detailed explanation goes here
  arguments
    MPSparams.LeftMatrices (1,:) cell;
    MPSparams.RightMatrices (1,:) cell;
    MPSparams.CoreMatrix (:,:) double;
  end
  MPS = MPSparams;
  MPS.left = size(MPSparams.LeftMatrices,2);
  MPS.right = size(MPSparams.RightMatrices,2);
  MPS.CL = MPS.left + MPS.right;
end

