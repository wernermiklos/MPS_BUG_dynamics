function Psi_out = HdotPsi_lqr(LeftBlock, RightBlock, MPO_matrix, MPO_tasklist, Psi_in)
%HDOTPSI_LR Summary of this function goes here
%   Detailed explanation goes here
  %OpIDs = find(LeftBlock.OP_logical & RightBlock.OP_logical)';
  Psi_out = 0;
  something_done = false;
  for taskrow = 1:size(MPO_tasklist,1)
    OP_loc = MPO_matrix{MPO_tasklist(taskrow,1),MPO_tasklist(taskrow,2)};
    OP_left = LeftBlock.OP{MPO_tasklist(taskrow,1)};
    OP_right = RightBlock.OP{MPO_tasklist(taskrow,2)};
    if ~isempty(OP_loc) && ~isempty(OP_left) && ~isempty(OP_right)
      %tmp = tensorprod(OP_right,Psi_in,2,3);
      %tmp = tensorprod(OP_loc,tmp,2,3);
      %tmp = tensorprod(OP_left,tmp,2,3);
      
      tmp = tensorprod(Psi_in,OP_left,1,2,"NumDimensionsA",3);
      tmp = tensorprod(tmp,OP_loc,1,2,"NumDimensionsA",3);
      tmp = tensorprod(tmp,OP_right,1,2,"NumDimensionsA",3);
      
      %tmp = permute(tensorprod(Psi_in,OP_loc,2,2,"NumDimensionsA",3),[1,3,2]);
      %tmp = tensorprod(OP_left,tmp,2,1);
      %tmp = tensorprod(tmp,OP_right,3,1,"NumDimensionsA",3);
      
      Psi_out = Psi_out + tmp;
      something_done = true;
    end
  end
  if ~something_done
    Psi_out = zeros(size(Psi_in));
  end
end

