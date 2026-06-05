function [XM_out,iter] = expevolv_matrix_shifttrick(XM_in,dt,HdotX, exporder)
%EXPEVOLV_MATRIX Summary of this function goes here
%   E = <PSI|H|PSI>
%   exp(-i H dt) |PSI> = exp(-i E dt) exp(-i(H-E)dt) |PSI>
%   exp(-i(H-E)dt) |PSI> is first expanded in series, and the trivial phase
%   rotation is only performed at the end. 
%
  iter = 0;
  keepgoing = true;
  %Y = XM_in;
  XM_out = XM_in;
  norm0 = norm(XM_in(:));
  Ytmp = HdotX(XM_in);
  E = XM_in(:)'*Ytmp(:);
  HdotX_act = @(x) HdotX(x) - E*x;
  while iter < exporder || (isnan(exporder) && keepgoing)
    iter = iter + 1;
    if iter ~= 1
      Y = -1i*dt*HdotX_act(Y)/iter;
    else
      Y = -1i*dt*Ytmp + 1i*dt*E*XM_in;
    end
    XM_out = XM_out + Y;
    if isnan(exporder)
      mynorm = norm(Y(:));
      if mynorm/norm0 < 1e-15
        keepgoing = false;
      end
    end
  end
  XM_out = exp(-1i*E*dt)*XM_out;
end

