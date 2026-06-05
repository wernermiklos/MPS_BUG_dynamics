function [MPS,Blocks,ROTmatrix_L] = BUG_augmentsweep_left(MPS,Blocks,dt,Model, bond_dim, options)
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
  left_orig = length(MPS.LeftMatrices);
  LeftMatrices_orig = MPS.LeftMatrices;
  [MPS, Blocks] = move_orthogonality_center(MPS,Blocks,Model,-left_orig,bond_dim); 
  ROTmatrix_L = 1;
  for pos = 1:left_orig
    PSI0_lqr = permute(tensorprod(MPS.CoreMatrix,MPS.RightMatrices{end},2,3),[1,3,2]);
    MPS.RightMatrices(end) = [];
    Blocks.Right(end) = [];
    myLeftBlock = Blocks.Left{pos};
    myRightBlock = Blocks.Right{end};
    myHdotPsi = @(x) HdotPsi_lqr(myLeftBlock,myRightBlock,Model.MPO{pos}, Model.MPOtasktable{pos},x);
    [PSInew_lqr, N_iter] = expevolv_matrix_shifttrick(PSI0_lqr,dt,myHdotPsi,options.ExpOrder);
    if options.KeepPSI0
      LeftMPSMatrix_new = orth_keepPSI0(PSI0_lqr,PSInew_lqr);
    else
      [LeftMPSMatrix_new,~] = tensor_qr(cat(3,PSI0_lqr,PSInew_lqr),[1,2],[3],'econ');
    end
    CoreMatrix_new = tensorprod(conj(LeftMPSMatrix_new),PSI0_lqr,[1,2],[1,2],NumDimensionsA=3);
    ROTmatrix_L_tmp = tensorprod(ROTmatrix_L,LeftMatrices_orig{pos},2,1,NumDimensionsA=2);
    ROTmatrix_L = tensorprod(conj(LeftMPSMatrix_new),ROTmatrix_L_tmp,[1,2],[1,2]);
    MPS.LeftMatrices{end+1} = LeftMPSMatrix_new;
    MPS.CoreMatrix = CoreMatrix_new;
    Blocks.Left{end+1} = renorm(Blocks.Left{end},Model.MPO{pos}, Model.MPOlogical{pos},LeftMPSMatrix_new, "SIDE", 'left');
  end
  MPS.left = length(MPS.LeftMatrices);
  MPS.right = length(MPS.RightMatrices);
end

