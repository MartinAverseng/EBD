function[C] = localCorrections(x,y,a,kernel,rq,rMax,tol)
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
    B1_inds = kernel.func(rxyApply);
    C1 = sum(abs(rq.alpha0.*Cp(rq.rho)).*rq.rho.^2);
    Ninterp = fix(sqrt(C1)*a/sqrt(8*tol))+10;
    xinterp = linspace(0,a*1.1,Ninterp);
    yinterp = rq.eval(xinterp);
    B2_inds = interp1(xinterp,yinterp,rxy);% - op.q2d.offset;% + op.q2d.offset;
    %
    B_inds = B1_inds - B2_inds;
    C = sparse(idx,jdx,B_inds,N1,N2);
else
    % No close interactions
    C = sparse(N1,N2); % all zeros
end