
function [out] = coeff2der(alpha,rho,x)
% [out] = coeff2func(alpha,rho,x) 
% Returns the values at x of the function 
% \sum_{p} alpha_p e_p(r) 
% Where e_p(r) = C_p J_0(rho_p r) (see function 'Cp' for a definition of 
% C_p). 

y = x(:)';
C = Cp(rho);

J1vals = -besselj(1,rho(:)*y);
out = J1vals'*(C(:).*rho(:).*alpha(:));

reshape(out,size(x,1),size(x,2));

end
