function [Q,R] = tensor_qr(myTEN,Qlegs,Rlegs,mode)
%TENSOR_QR Summary of this function goes here
%   Detailed explanation goes here
  arguments
    myTEN double
    Qlegs (1,:) double
    Rlegs (1,:) double
    mode {mustBeMember(mode, {'econ', 'standard'})} = 'standard'

  end
  tmp = permute(myTEN,[Qlegs,Rlegs]);
  sizeQlegs_orig = size(myTEN,Qlegs);
  dimQtot = prod(sizeQlegs_orig);
  sizeRlegs_orig = size(myTEN,Rlegs);
  dimRtot = prod(sizeRlegs_orig);

  tmp = reshape(tmp,[dimQtot,dimRtot]);
  switch mode
    case 'standard'
      [Q,R] = qr(tmp);   %no econ!!!!
    case 'econ'
      [Q,R] = qr(tmp,'econ');
  end
  Q = reshape(Q,[sizeQlegs_orig,size(Q,2)]);
  R = reshape(R,[size(R,1),sizeRlegs_orig]);
end

