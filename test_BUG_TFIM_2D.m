% ----- run parameters ----
dt = 0.01;
tmax = 1.0;
Lx = 4;
Ly = 4;
CL = Lx*Ly;
h = 2;
S = 1/2;
M=64;



% ----- Checkerboard state: ------
cut = floor(CL/2);
LeftMatrices = cell(1,cut);
RightMatrices = cell(1,CL-cut);
pos = 0;
for x = 1:Lx
  for y = 1:Ly
    pos = pos + 1;
    if pos <= cut
      LeftMatrices{pos} = zeros(1,2*S+1);
      if mod(x+y,2)
        LeftMatrices{pos}(1) = 1;
      else
        LeftMatrices{pos}(end) = 1;
      end
    else
      rightpos = CL-pos+1;
      RightMatrices{rightpos} = zeros(1,2*S+1);
      if mod(x+y,2)
        RightMatrices{rightpos}(1) = 1;
      else
        RightMatrices{rightpos}(end) = 1;
      end
    end

  end
end

CoreMatrix = 1;

% ----- Initialization ------



[BUGrun,Model,MPS,Blocks] = Setup_BUG_withMPO('dt',0.01,...
                                               'tmax',1,...
                                               'ModelParams',struct('Lx',Lx,...
                                                                    'Ly',Ly,...
                                                                    'h',h,...
                                                                    'S',S,...
                                                                     'ModelInput','model_TFIM_2D_OBC'),...
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
  disp([' Sx_BUG = ', num2str(tensorprod(Blocks.Left{end}.OP{end}*MPS.CoreMatrix, conj(MPS.CoreMatrix),[1,2],[1,2]))])
  disp(['t = ', num2str(t)]);
end



