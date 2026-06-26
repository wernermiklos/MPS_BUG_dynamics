% ----- run parameters ----
dt = 0.01;
tmax = 1.0;
L = 40;
h = 2;
S = 5/2;
M=64;
[~,~,Sz_loc] = SpinOperators(S);



% ----- Neél state: ------
LeftMatrices = cell(1,L/2);
for pos = 1:(L/2)
  LeftMatrices{pos} = zeros(1,2*S+1);
  if mod(pos,2)
    LeftMatrices{pos}(1) = 1;
  else
    LeftMatrices{pos}(end) = 1;
  end
end

RightMatrices = cell(1,L/2);
for pos = 1:(L/2)
  RightMatrices{pos} = zeros(1,2*S+1);
  if mod(pos,2)
    RightMatrices{pos}(end) = 1;
  else
    RightMatrices{pos}(1) = 1;
  end
end

CoreMatrix = 1;

% ----- Initialization ------



[BUGrun,Model,MPS,Blocks] = Setup_BUG_withMPO('dt',0.01,...
                                               'tmax',1,...
                                               'ModelParams',struct('L',L,...
                                                                    'h',h,...
                                                                    'S',S,...
                                                                     'ModelInput','model_TFIM_1D_OBC'),...
                                               'InitMPS',struct('LeftMatrices',{LeftMatrices},...
                                                                'RightMatrices',{RightMatrices},...
                                                                'CoreMatrix',CoreMatrix));
disp('----------')
while size(MPS.CoreMatrix,1) < M
  [MPS,Blocks] = BUG_augmentsweep_left(MPS,Blocks,dt,Model,M);
  if size(MPS.CoreMatrix,2) < M
    [MPS,Blocks] = BUG_augmentsweep_right(MPS,Blocks,dt,Model,M);
  end
  disp(['Initial augmentation M = ', num2str(size(MPS.CoreMatrix))]);
end

for t = dt:dt:tmax
  disp('----------')
  tic; [MPS, Blocks, N_iter] = BUG_fullstep(MPS,Blocks,Model,dt,M); toc;
  SzData = measure_local_expval(MPS,Sz_loc);
  disp(['t = ', num2str(t)]);
  fprintf('<Sz> = ')
  fprintf('%.3f ',SzData);
  fprintf('\n')
end


% [MPS,Blocks] = move_orthogonality_center(MPS,Blocks,Model,-(length(MPS.LeftMatrices)),M);
% [MPS,Blocks] = move_orthogonality_center(MPS,Blocks,Model,Model.CL/2,M);
% disp('----------')
% [MPS,Blocks] = BUG_augmentsweep_left(MPS,Blocks,dt,Model,M);
% [MPS,Blocks] = BUG_augmentsweep_right(MPS,Blocks,dt,Model,M);
% [MPS,Blocks] = BUG_augmentsweep_left(MPS,Blocks,dt,Model,M);
% [MPS,Blocks] = BUG_augmentsweep_right(MPS,Blocks,dt,Model,M);
% [MPS,Blocks] = BUG_augmentsweep_left(MPS,Blocks,dt,Model,M);
% [MPS,Blocks] = BUG_augmentsweep_right(MPS,Blocks,dt,Model,M);
% [MPS,Blocks] = BUG_augmentsweep_left(MPS,Blocks,dt,Model,M);
% [MPS,Blocks] = BUG_augmentsweep_right(MPS,Blocks,dt,Model,M);
% 
% for t = dt:dt:tmax
%   disp('----------')
%   [MPS_left,Blocks_left,ROTmatrix_left] = BUG_augmentsweep_left(MPS,Blocks,dt,Model,M);
%   [MPS,Blocks,ROTmatrix_right] = BUG_augmentsweep_right(MPS,Blocks,dt,Model,M);
%   MPS.CoreMatrix = ROTmatrix_left * MPS.CoreMatrix;
%   MPS.LeftMatrices = MPS_left.LeftMatrices;
%   Blocks.Left = Blocks_left.Left;
%   Blocks_left = [];
%   MPS_left = [];
%   myRightBlock = Blocks.Right{end};
%   myLeftBlock = Blocks.Left{end};
%   myHdotPsi = @(x) HdotPsi_LR(myLeftBlock,myRightBlock,x);
%   [newcorematrix,N_iter] = expevolv_matrix_shifttrick(MPS.CoreMatrix,dt,myHdotPsi,NaN);
%   %truncation
%   [U,S,V] = svd(newcorematrix,'econ');
%   if size(S,1) > M
%     S = S(1:M,1:M);
%     U = U(:,1:M);
%     V = V(:,1:M);
%     newcorematrix = S;
%     for opID = 1:length(Blocks.Left{end}.OP)
%       if Blocks.Left{end}.OP_logical(opID)
%         Blocks.Left{end}.OP{opID} = U'*Blocks.Left{end}.OP{opID}*U;
%       end
%     end
%     for opID = 1:length(Blocks.Right{end}.OP)
%       if Blocks.Right{end}.OP_logical(opID)
%         Blocks.Right{end}.OP{opID} = conj(V)'*Blocks.Right{end}.OP{opID}*conj(V);
%       end
%     end
%     MPS.LeftMatrices{end} = tensorprod(MPS.LeftMatrices{end},U,3,1);
%     MPS.RightMatrices{end} = tensorprod(MPS.RightMatrices{end},conj(V),3,1);
%   end
% 
%   MPS.CoreMatrix = newcorematrix;
%   disp([' Sx_BUG = ', num2str(tensorprod(Blocks.Left{end}.OP{2}*MPS.CoreMatrix, conj(MPS.CoreMatrix),[1,2],[1,2]))])
%   disp(['t = ', num2str(t)]);
% end



