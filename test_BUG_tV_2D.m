% ----- run parameters ----
dt = 0.01;
tmax = 1.0;
Lx = 4;
Ly = 4;
CL = Lx*Ly;
V = 1;
M=256;



% ----- Checkerboard (CDW) state: ------
cut = floor(CL/2);
LeftMatrices = cell(1,cut);
RightMatrices = cell(1,CL-cut);
pos = 0;
for x = 1:Lx
  for y = 1:Ly
    pos = pos + 1;
    if pos <= cut
      LeftMatrices{pos} = zeros(1,2);
      if mod(x+y,2)
        LeftMatrices{pos}(1) = 1;
      else
        LeftMatrices{pos}(2) = 1;
      end
    else
      rightpos = CL-pos+1;
      RightMatrices{rightpos} = zeros(1,2);
      if mod(x+y,2)
        RightMatrices{rightpos}(1) = 1;
      else
        RightMatrices{rightpos}(2) = 1;
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
                                                                    'V',V,...
                                                                     'ModelInput','model_tV_2D_OBC'),...
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
  disp([' Sx_BUG = ', num2str(tensorprod(Blocks.Left{end}.OP{end}*MPS.CoreMatrix, conj(MPS.CoreMatrix),[1,2],[1,2]))])
  disp(['t = ', num2str(t)]);
end



