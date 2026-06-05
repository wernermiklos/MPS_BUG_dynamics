function block_out = renorm_diffbraket(block_in, MPOmatrix, MPOmatrix_logical, MPSmatrix_ket,MPSmatrix_bra, options)
%RENORM Summary of this function goes here
%   Detailed explanation goes here
arguments
  block_in struct
  MPOmatrix (:,:) cell
  MPOmatrix_logical (:,:) logical;
  MPSmatrix_ket (:,:,:) double;   %block_in local block_out
  MPSmatrix_bra (:,:,:) double;   %block_in local block_out
  options.SIDE (1,1) string = 'left';
end
  if ~strcmp(block_in.side, options.SIDE)
    error('side mismatch in renorm')
  end
  block_out = struct();
  block_out.side = options.SIDE;
  switch options.SIDE
    case 'left'
      block_out.OP_logical = (logical(block_in.OP_logical.'*MPOmatrix_logical)).';
      block_out.OP = cell(size(MPOmatrix,2),1);
      [MPOrow,MPOcol] = find(MPOmatrix_logical);
      for taskID = 1:length(MPOrow)
        OP_block_in = block_in.OP{MPOrow(taskID)};
        OP_loc = MPOmatrix{MPOrow(taskID),MPOcol(taskID)};
        tmp = tensorprod(OP_block_in,MPSmatrix_ket,2,1,'NumDimensionsA',2);
        tmp = tensorprod(tmp,OP_loc,2,2,'NumDimensionsA',3);
        tmp = tensorprod(MPSmatrix_bra,tmp,[1,2],[1,3],'NumDimensionsA',3);
        if isempty(block_out.OP{MPOcol(taskID)})
          block_out.OP{MPOcol(taskID)} = tmp;
        else
          block_out.OP{MPOcol(taskID)} = block_out.OP{MPOcol(taskID)} + tmp;
        end
      end
    case 'right'
      block_out.OP_logical = logical(MPOmatrix_logical*block_in.OP_logical);
      block_out.OP = cell(size(MPOmatrix,1),1);
      [MPOrow,MPOcol] = find(MPOmatrix_logical);
      for taskID = 1:length(MPOrow);
        OP_block_in = block_in.OP{MPOcol(taskID)};
        OP_loc = MPOmatrix{MPOrow(taskID),MPOcol(taskID)};
        tmp = tensorprod(OP_block_in,MPSmatrix_ket,2,1,'NumDimensionsA',2);
        tmp = tensorprod(tmp,OP_loc,2,2,'NumDimensionsA',3);
        tmp = tensorprod(MPSmatrix_bra,tmp,[1,2],[1,3],'NumDimensionsA',3);
        if isempty(block_out.OP{MPOrow(taskID)})
          block_out.OP{MPOrow(taskID)} = tmp;
        else
          block_out.OP{MPOrow(taskID)} = block_out.OP{MPOrow(taskID)} + tmp;
        end
      end
    otherwise
      error('left or right renormalization...')
  end
end

