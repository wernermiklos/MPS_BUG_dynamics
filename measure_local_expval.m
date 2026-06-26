function OutData = measure_local_expval(MPS,localop)
%MEASURE_LOCAL_EXPVAL Measures expectation value of localop for every site
%   Detailed explanation goes here
CL = MPS.CL;
left = MPS.left;
right = MPS.right;
OutData = zeros(1,CL);
tmp_left = {};
tmp_right = {};
for ii = 1:left
  %op1 (x) id
  for kk = 1:(ii-1)
    tmp_tensor = tensorprod(tmp_left{kk},MPS.LeftMatrices{ii},2,1);
    tmp_left{kk} = tensorprod(conj(MPS.LeftMatrices{ii}),tmp_tensor,[1,2],[1,2],"NumDimensionsA",3);
  end
  % id (x) op1
  tmp_tensor = tensorprod(MPS.LeftMatrices{ii},localop,2,2,"NumDimensionsA",3);
  tmp_left{end+1} = tensorprod(conj(MPS.LeftMatrices{ii}),tmp_tensor,[1,2],[1,3],"NumDimensionsA",3);
end

for ii = 1:left
  tmp_mx = tmp_left{ii}*MPS.CoreMatrix;
  OutData(ii) = MPS.CoreMatrix(:)' * tmp_mx(:);
end


for ii = 1:right
  %id (x) op1
  for kk = 1:(ii-1)
    tmp_tensor = tensorprod(tmp_right{kk},MPS.RightMatrices{ii},2,1);
    tmp_right{kk} = tensorprod(conj(MPS.RightMatrices{ii}),tmp_tensor,[1,2],[1,2],"NumDimensionsA",3);
  end
  % op1 (x) id
  tmp_tensor = tensorprod(MPS.RightMatrices{ii},localop,2,2,"NumDimensionsA",3);
  tmp_right{end+1} = tensorprod(conj(MPS.RightMatrices{ii}),tmp_tensor,[1,2],[1,3],"NumDimensionsA",3);
end

for ii = 1:right 
  tmp_mx = MPS.CoreMatrix*(tmp_right{ii}.');
  OutData(CL-ii+1) = MPS.CoreMatrix(:)' * tmp_mx(:);
end


   
end

