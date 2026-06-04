function [MPS, Blocks] = move_orthogonality_center(MPS,Blocks,Model,shift, bond_dim, mode)
%MOVE_ORTHOGONALITY_CENTER Summary of this function goes here
%   Detailed explanation goes here
  arguments
    MPS struct
    Blocks struct
    Model struct
    shift (1,1) double
    bond_dim (1,1) double = NaN;
    mode {mustBeMember(mode, {'qr', 'svd'})} = 'qr'
  end
  
  if ~isnan(bond_dim)
    mode='svd';
  end

  if shift == 0
    return;
  elseif shift > 0
    if shift > length(MPS.RightMatrices)
      error('shift cannot be larger then length of right matrices')
    end
    for iter = 1:shift
      [MPS,Blocks] = internal_move_center_right(MPS,Blocks,Model,bond_dim,mode);
    end
  elseif shift < 0
    if abs(shift) > length(MPS.LeftMatrices)
      error('-shift cannot be larger then length of left matrices')
    end
    for iter = 1:abs(shift)
      [MPS,Blocks] = internal_move_center_left(MPS,Blocks,Model,bond_dim,mode);
    end
  end
end


function [MPS,Blocks] = internal_move_center_right(MPS,Blocks,Model,M,mode)
  pos = length(MPS.LeftMatrices)+1;
  tmp = tensorprod(MPS.RightMatrices{end},MPS.CoreMatrix,3,2);
  switch mode
    case 'qr'
      [newleftmatrix,newcorematrix] = tensor_qr(tmp,[3,2],[1],'econ');
      MPS.LeftMatrices{end+1} = newleftmatrix;
      MPS.CoreMatrix = newcorematrix;
      MPS.RightMatrices(end) = [];
    case 'svd'
      [newleftmatrix,newcorematrix,V] = tensor_svd(tmp,[3,2],[1],'standard');
      sizeCore = size(newcorematrix,[1,2]);
      if ~isnan(M) && (sizeCore(1) > M || sizeCore(2) > M)
        norm_orig = norm(diag(newcorematrix));
        newcorematrix = newcorematrix(1:min(sizeCore(1),M),1:(min(sizeCore(2),M)));
        newcorematrix = norm_orig / norm(diag(newcorematrix))*newcorematrix;
        newleftmatrix = newleftmatrix(:,:,1:min(sizeCore(1),M));
        V = V(:,1:min(sizeCore(2),M));
      end
      MPS.LeftMatrices{end+1} = newleftmatrix;
      MPS.CoreMatrix = newcorematrix*V';
      MPS.RightMatrices(end) = [];
  end
  MPS.left = length(MPS.LeftMatrices);
  MPS.right = length(MPS.RightMatrices);
  if ~isempty(Blocks)
    Blocks.Right(end) = [];
    Blocks.Left{end+1} = renorm(Blocks.Left{end},Model.MPO{pos},Model.MPOlogical{pos},newleftmatrix,'SIDE','left');
  end
end


function [MPS,Blocks] = internal_move_center_left(MPS,Blocks,Model,M,mode)
  pos = length(MPS.LeftMatrices);
  tmp = tensorprod(MPS.LeftMatrices{end},MPS.CoreMatrix,3,1);
  switch mode
    case 'qr'
      [newrightmatrix,newcorematrix] = tensor_qr(tmp,[3,2],[1],'econ');
      MPS.RightMatrices{end+1} = newrightmatrix;
      MPS.CoreMatrix = newcorematrix.';
      MPS.LeftMatrices(end) = [];
    case 'svd'
      [newrightmatrix,newcorematrix,V] = tensor_svd(tmp,[3,2],[1],'standard');
      sizeCore = size(newcorematrix,[1,2]);
      if ~isnan(M) && (sizeCore(1) > M || sizeCore(2) > M)
        norm_orig = norm(diag(newcorematrix));
        newcorematrix = newcorematrix(1:min(sizeCore(1),M),1:(min(sizeCore(2),M)));
        newcorematrix = norm_orig / norm(diag(newcorematrix))*newcorematrix;
        newrightmatrix = newrightmatrix(:,:,1:min(sizeCore(1),M));
        V = V(:,1:min(sizeCore(2),M));
      end
      MPS.RightMatrices{end+1} = newrightmatrix;
      MPS.CoreMatrix = (newcorematrix*V').';
      MPS.LeftMatrices(end) = [];
  end
  MPS.left = length(MPS.LeftMatrices);
  MPS.right = length(MPS.RightMatrices);
  if ~isempty(Blocks)
    Blocks.Left(end) = [];
    Blocks.Right{end+1} = renorm(Blocks.Right{end},Model.MPO{pos},Model.MPOlogical{pos},newrightmatrix,'SIDE','right');
  end
end
