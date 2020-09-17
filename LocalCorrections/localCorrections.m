function[C] = localCorrections(x,y,a,G,rq,tol,gradOpt)
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
clear rxyTemp;
idx = zeros(size(jdx));
sp_ind = 1;
for x_ind=1:length(I)
    idx(sp_ind:(sp_ind+length(I{x_ind})-1)) = x_ind;
    sp_ind = sp_ind + length(I{x_ind});
end
% Save memory
clear I;
if gradOpt
    idx = idx(rxy> 1e-12);
    jdx = jdx(rxy> 1e-12);
    rxy = rxy(rxy> 1e-12);
end
NCI = length(rxy);

if NCI ~= 0
    
    if gradOpt
        C1 = sum(abs(rq.alpha0.*Cp(rq.rho)).*rq.rho.^3); % Bound for the second derivative.
        Ninterp = fix(sqrt(C1)*a/sqrt(8*tol))+10; % Guarantees interpolation error < tol.
        xinterp = linspace(0,a*1.1,Ninterp);
        yinterp = rq.evalDer(xinterp);
        radial_quadratureNear0 = interp1(xinterp,yinterp,rxy);
        R =  sparse(idx,jdx,1./rxy,N1,N2);
        rxy = G.evalDer(rxy) - radial_quadratureNear0;
        clear radial_quadratureNear0;
        clear xinterp
        clear yinterp
        Tx = sparse(idx,jdx,y(jdx,1) - x(idx,1),N1,N2);
        Ty = sparse(idx,jdx,y(jdx,2) - x(idx,2),N1,N2);
        C = sparse(idx,jdx,rxy,N1,N2);
        Cx = C.*R.*Tx;
        Cy = C.*R.*Ty;
        C = {Cx,Cy};
    else
        C1 = sum(abs(rq.alpha0.*Cp(rq.rho)).*rq.rho.^2); % Bound for the second derivative.
        Ninterp = fix(sqrt(C1)*a/sqrt(8*tol))+10; % Guarantees interpolation error < tol.
        xinterp = linspace(0,a*1.1,Ninterp);
        yinterp = rq.eval(xinterp);
        val = G.eval(rxy)- interp1(xinterp,yinterp,rxy); % Exact local interactions
        clear rxy
        % we remove the radial
        clear radial_quadratureNear0;
        clear xinterp
        clear yinterp
        C = sparse(idx,jdx,val,N1,N2);
        clear val
    end
    
    
    % quadrature contribution (using interpolation to avoid evaluating at
    % all points).
    
else
    % No close interactions
    if gradOpt
        C = {sparse(N1,N2),sparse(N1,N2)};
    else
        C = sparse(N1,N2); % all zeros
    end
end