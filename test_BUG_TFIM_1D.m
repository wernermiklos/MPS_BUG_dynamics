% ----- run parameters ----
dt = 0.01;
tmax = 1.0;
L = 40;
h = 2;
M=64;



% ----- Neél state: ------
LeftMatrices = cell(1,L/2);
for pos = 1:(L/2)
  if mod(pos,2)
    LeftMatrices{pos} = [1,0];
  else
    LeftMatrices{pos} = [0,1];
  end
end

RightMatrices = cell(1,L/2);
for pos = 1:(L/2)
  if mod(pos,2)
    RightMatrices{pos} = [0,1];
  else
    RightMatrices{pos} = [1,0];
  end
end

CoreMatrix = 1;

% ----- Initialization ------



[BUGrun,Model,MPS,Blocks] = Setup_BUG_withMPO('dt',0.01,...
                                               'tmax',1,...
                                               'ModelParams',struct('L',L,...
                                                                    'h',h,...
                                                                     'ModelInput','model_TFIM_1D_OBC'),...
                                               'InitMPS',struct('LeftMatrices',{LeftMatrices},...
                                                                'RightMatrices',{RightMatrices},...
                                                                'CoreMatrix',CoreMatrix));
[MPS,Blocks] = move_orthogonality_center(MPS,Blocks,Model,-(length(MPS.LeftMatrices)),M);
[MPS,Blocks] = move_orthogonality_center(MPS,Blocks,Model,Model.CL/2,M);
disp('----------')
[MPS,Blocks] = BUG_augmentsweep_left(MPS,Blocks,dt,Model,M);
[MPS,Blocks] = BUG_augmentsweep_right(MPS,Blocks,dt,Model,M);
[MPS,Blocks] = BUG_augmentsweep_left(MPS,Blocks,dt,Model,M);
[MPS,Blocks] = BUG_augmentsweep_right(MPS,Blocks,dt,Model,M);
[MPS,Blocks] = BUG_augmentsweep_left(MPS,Blocks,dt,Model,M);
[MPS,Blocks] = BUG_augmentsweep_right(MPS,Blocks,dt,Model,M);
[MPS,Blocks] = BUG_augmentsweep_left(MPS,Blocks,dt,Model,M);
[MPS,Blocks] = BUG_augmentsweep_right(MPS,Blocks,dt,Model,M);

for t = dt:dt:tmax
  disp('----------')
  [MPS_left,Blocks_left,ROTmatrix_left] = BUG_augmentsweep_left(MPS,Blocks,dt,Model,M);
  [MPS,Blocks,ROTmatrix_right] = BUG_augmentsweep_right(MPS,Blocks,dt,Model,M);
  MPS.CoreMatrix = ROTmatrix_left * MPS.CoreMatrix;
  MPS.LeftMatrices = MPS_left.LeftMatrices;
  Blocks.Left = Blocks_left.Left;
  Blocks_left = [];
  MPS_left = [];
  myRightBlock = Blocks.Right{end};
  myLeftBlock = Blocks.Left{end};
  myHdotPsi = @(x) HdotPsi_LR(myLeftBlock,myRightBlock,x);
  [newcorematrix,N_iter] = expevolv_matrix_shifttrick(MPS.CoreMatrix,dt,myHdotPsi,NaN);
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
  disp([' Sx_BUG = ', num2str(tensorprod(Blocks.Left{end}.OP{2}*MPS.CoreMatrix, conj(MPS.CoreMatrix),[1,2],[1,2]))])
  disp(['t = ', num2str(t)]);
end



