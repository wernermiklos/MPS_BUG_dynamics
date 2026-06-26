function [MPS, Blocks, N_iter] = BUG_fullstep(MPS, Blocks, Model, dt, M)
%BUG_FULLSTEP Summary of this function goes here
%   Detailed explanation goes here
  [MPS_left,Blocks_left,ROTmatrix_left] = BUG_augmentsweep_left(MPS,Blocks,dt,Model,M);
  [MPS,Blocks,~] = BUG_augmentsweep_right(MPS,Blocks,dt,Model,M);
  MPS.CoreMatrix = ROTmatrix_left * MPS.CoreMatrix;
  MPS.LeftMatrices = MPS_left.LeftMatrices;
  Blocks.Left = Blocks_left.Left;
  Blocks_left = [];
  MPS_left = [];
  myRightBlock = Blocks.Right{end};
  myLeftBlock = Blocks.Left{end};
  myHdotPsi = @(x) HdotPsi_LR(myLeftBlock,myRightBlock,x);
  [newcorematrix,N_iter] = expevolv_matrix_Krylov(MPS.CoreMatrix,dt,myHdotPsi,NaN);
  %truncation
  [U,S,V] = svd(newcorematrix,'econ');
  if size(S,1) > M
    S = S(1:M,1:M);
    U = U(:,1:M);
    V = V(:,1:M);
    newcorematrix = S;
    for opID = 1:length(Blocks.Left{end}.OP)
      if Blocks.Left{end}.OP_logical(opID)
        Blocks.Left{end}.OP{opID} = U'*Blocks.Left{end}.OP{opID}*U;
      end
    end
    for opID = 1:length(Blocks.Right{end}.OP)
      if Blocks.Right{end}.OP_logical(opID)
        Blocks.Right{end}.OP{opID} = conj(V)'*Blocks.Right{end}.OP{opID}*conj(V);
      end
    end
    MPS.LeftMatrices{end} = tensorprod(MPS.LeftMatrices{end},U,3,1);
    MPS.RightMatrices{end} = tensorprod(MPS.RightMatrices{end},conj(V),3,1);
  end
  MPS.CoreMatrix = newcorematrix;
end

