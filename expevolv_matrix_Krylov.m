function [XM_out,i] = expevolv_matrix_Krylov(XM_in,dt,HdotX, exporder)
%EXPEVOLV_MATRIX works only for hermitian HdotX
%   E = <PSI|H|PSI>
%   exp(-i H dt) |PSI> = exp(-i E dt) exp(-i(H-E)dt) |PSI>
%   exp(-i(H-E)dt) |PSI> is first expanded in series, and the trivial phase
%   rotation is only performed at the end. 
%
  ITER_MAX = 500;
  EPSILON = 1e-13;
  i = 0;
  keepgoing = true;
  %Y = XM_in;
  norm0 = norm(XM_in(:));
  v = {XM_in/norm0};
  w = {};

  while i < exporder || (isnan(exporder) && keepgoing)
    i = i + 1;
    w{end+1} = HdotX(v{end});
    vNew = w{end};
    for j = 1:i
      T(j,i) = v{j}(:)'*vNew(:);
      vNew = vNew - T(j,i)*v{j};
      T(i,j) = conj(T(j,i));
    end
    norm_vNew = norm(vNew(:));
    %fprintf('%d -- %.3e \n',i, norm_vNew)
    v{end+1} = vNew / norm_vNew;
    expT = expm(-1i*dt*T);
    if abs(dt)*norm_vNew*abs(expT(i,1)) < EPSILON || norm_vNew < EPSILON || i > ITER_MAX
      keepgoing = false;
    end
  end
  XM_out = 0;
  for j = 1:i
    XM_out = XM_out + expT(j,1)*v{j};
  end
  XM_out = norm0*XM_out;
  fprintf('%d \n',i);
end





  
%     i = i + 1;
%     if i ~= 1
%       Y = -1i*dt*HdotX_act(Y)/i;
%     else
%       Y = -1i*dt*Ytmp + 1i*dt*E*XM_in;
%     end
%     XM_out = XM_out + Y;
%     if isnan(exporder)
%       mynorm = norm(Y(:));
%       if mynorm/norm0 < 1e-15
%         keepgoing = false;
%       end
%     end
%   end
% 
% end


%  E = XM_in(:)'*Ytmp(:);
%HdotX_act = @(x) HdotX(x) - E*x;
%  XM_out = exp(-1i*E*dt)*XM_out;