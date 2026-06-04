function Psi_out = HdotPsi_LR(LeftBlock, RightBlock,Psi_in)
%HDOTPSI_LR Summary of this function goes here
%   Detailed explanation goes here
  OpIDs = find(LeftBlock.OP_logical & RightBlock.OP_logical)';
  Psi_out = [];
  for opID = OpIDs
    if isempty(Psi_out)
      Psi_out = (LeftBlock.OP{opID}*Psi_in)*RightBlock.OP{opID}.';
    else
      Psi_out = Psi_out + (LeftBlock.OP{opID}*Psi_in)*RightBlock.OP{opID}.';
    end
  end
  if isempty(Psi_out)
    Psi_out = zeros(size(Psi_in));
  end
end

