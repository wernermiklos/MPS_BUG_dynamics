function [U,S,V] = tensor_svd(myTEN,Ulegs,Vlegs,mode)
%TENSOR_SVD performs tensor SVD (no truncation here)
%   myTEN = sum_ij U(:,:,i) S(i,j) conj(V(:,:,j))
  arguments
    myTEN double
    Ulegs (1,:) double
    Vlegs (1,:) double
    mode {mustBeMember(mode, {'econ', 'standard'})} = 'standard'
  end

  tmp = permute(myTEN,[Ulegs,Vlegs]);
  sizeUlegs_orig = size(myTEN,Ulegs);
  dimUtot = prod(sizeUlegs_orig);
  sizeVlegs_orig = size(myTEN,Vlegs);
  dimVtot = prod(sizeVlegs_orig);

  tmp = reshape(tmp,[dimUtot,dimVtot]);
  switch mode
    case 'standard'
      [U,S,V] = svd(tmp);   %no econ!!!!
    case 'econ'
      [U,S,V] = svd(tmp,'econ');
  end
  U = reshape(U,[sizeUlegs_orig,size(U,2)]);
  V = reshape(V,[sizeVlegs_orig,size(V,2),]);
end

