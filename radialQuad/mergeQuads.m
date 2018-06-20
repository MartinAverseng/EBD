function [ rho,alpha ] = mergeQuads(rho1,alpha1,rho2,alpha2)

rho1_2 = rho1(~ismember(rho1,rho2));
alpha1_2 = alpha1(~ismember(rho1,rho2));
rho2_1 = rho2(~ismember(rho2,rho1));
alpha2_1 = alpha2(~ismember(rho2,rho1));
[rho12,I1,I2] = intersect(rho1,rho2);
alpha12 = alpha1(I1) + alpha2(I2);

rho_merge = [rho1_2; rho2_1; rho12];
alpha_merge = [alpha1_2;alpha2_1;alpha12];

[rho,I] = sort(rho_merge);
alpha = alpha_merge(I);

end

