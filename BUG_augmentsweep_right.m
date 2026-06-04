function [MPS,Blocks, ROTmatrix_R] = BUG_augmentsweep_right(MPS,Blocks,dt,Model, bond_dim, options)
%BUG_AUGMENTSWEEP_LEFT Summary of this function goes here
%   Detailed explanation goes here
  arguments
    MPS struct
    Blocks struct
    dt (1,1) double
    Model struct
    bond_dim (1,1) double = NaN;
    options.ExpOrder (1,1) double = NaN;   %if NaN, then adaptive
    options.KeepPSI0 (1,1) logical = false;   %if true, then keeps the 'old' vecs first, and adds the news next to it.
  end
  right_orig = length(MPS.RightMatrices);
  RightMatrices_orig = MPS.RightMatrices;
  [MPS, Blocks] = move_orthogonality_center(MPS,Blocks,Model,+right_orig,bond_dim); 

  ROTmatrix_R = 1;

  for rightpos = 1:right_orig
    pos = Model.CL-rightpos+1;
    PSI0_lqr = tensorprod(MPS.LeftMatrices{end},MPS.CoreMatrix,3,1);
    dim_l = size(PSI0_lqr,1);
    dim_q = size(PSI0_lqr,2);
    dim_r = size(PSI0_lqr,3);
    PSI0_LR = reshape(PSI0_lqr,[dim_l, dim_q*dim_r]); 
    MPS.LeftMatrices(end) = [];
    Blocks.Left(end) = [];
    RightMPSMatrix_full = permute(reshape(eye(dim_r*dim_q),[dim_q,dim_r,dim_r*dim_q]),[2,1,3]);
    myRightBlock = renorm(Blocks.Right{rightpos},Model.MPO{pos},Model.MPOlogical{pos},RightMPSMatrix_full,"SIDE",'right');
    myLeftBlock = Blocks.Left{end};

    myHdotPsi = @(x) HdotPsi_LR(myLeftBlock,myRightBlock,x);
    [PSInew_LR, N_iter] = expevolv_matrix_shifttrick(PSI0_LR,dt,myHdotPsi,options.ExpOrder);
    if options.KeepPSI0
      RightMPSMatrix_new = orth_keepPSI0(PSI0_LR.',PSInew_LR.');
    else
      [RightMPSMatrix_new,~] = qr([PSI0_LR.',PSInew_LR.'],'econ');
    end
    RightMPSMatrix_new = permute(reshape(RightMPSMatrix_new,[dim_q,dim_r,size(RightMPSMatrix_new,2)]),[2,1,3]);
    CoreMatrix_new = tensorprod(PSI0_lqr,conj(RightMPSMatrix_new),[3,2],[1,2],NumDimensionsA=3);
    ROTmatrix_R_tmp = tensorprod(ROTmatrix_R,RightMatrices_orig{rightpos},2,1,NumDimensionsA=2);
    ROTmatrix_R = tensorprod(conj(RightMPSMatrix_new),ROTmatrix_R_tmp,[1,2],[1,2]);
    MPS.RightMatrices{end+1} = RightMPSMatrix_new;
    MPS.CoreMatrix = CoreMatrix_new;
    Blocks.Right{end+1} = renorm(Blocks.Right{end},Model.MPO{pos}, Model.MPOlogical{pos},RightMPSMatrix_new, "SIDE", 'right');
  end
  MPS.left = length(MPS.LeftMatrices);
  MPS.right = length(MPS.RightMatrices);
end

