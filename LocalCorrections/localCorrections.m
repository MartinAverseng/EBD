function[C] = localCorrections(x,y,a,G,rq,rMax,tol)
% C = localCorrections(x,y,a,G,rq,rMax,tol)
% Local correction matrix for the EBD.
% inputs : 
% - x and y : arrays of size N1x2 and N2x2 (clouds of points) inside the
% unit disk. 
% - a : the cutoff parameter of the EBD method
% - G : the kernel
% - rq : the radial quadrature that has been computed for G
% - rMax : such that initially, the data was X = rMax*x, Y = rMax*y (before
% rescaling). 
% - tol : required tolerance. 
% output : the sparse correction matrix C. 

N1 = size(x,1);
N2 = size(y,1);
[I,rxyTemp] = rangesearch(y,x,a*1.05);
jdx = cell2mat(I')';
rxy = cell2mat(rxyTemp')';
idx = zeros(size(jdx));
sp_ind = 1;
for x_ind=1:length(I);
    idx(sp_ind:(sp_ind+length(I{x_ind})-1)) = x_ind;
    sp_ind = sp_ind + length(I{x_ind});
end
% Save memory
clear I;
NCI = length(rxy);

if NCI ~= 0
    rxyApply = rxy;
    rxyApply(rxyApply < 1e-8) = 1e-8/rMax;
    exact_Interactions = G.func(rxyApply); % Exact local interactions
    C1 = sum(abs(rq.alpha0.*Cp(rq.rho)).*rq.rho.^2); % Bound for the second derivative. 
    Ninterp = fix(sqrt(C1)*a/sqrt(8*tol))+10; % Guarantees interpolation error < tol.
    xinterp = linspace(0,a*1.1,Ninterp);
    yinterp = rq.eval(xinterp);
    radial_quadratureNear0 = interp1(xinterp,yinterp,rxy); % we remove the radial 
    % quadrature contribution (using interpolation to avoid evaluating at
    % all points). 
    
    C_val = exact_Interactions - radial_quadratureNear0;
    C = sparse(idx,jdx,C_val,N1,N2);
else
    % No close interactions
    C = sparse(N1,N2); % all zeros
end