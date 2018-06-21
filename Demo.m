
% Fast computation of the vector q defined by 
% q_k = \sum_{l = 1}^{N_y} G(X_k - Y_l) f_l, k = 1 .. Nx
% Using the "Efficient Bessel Decomposition"
% Code developped by Martin Averseng
% See also SCSD method by François Alouges and Matthieu Aussal. 


addpath(genpath(pwd)); % Add folders of the toolbox to the path. 
clear all;
close all;
clc;


Nx = 10^5;
Ny = 10^4;
% Data points
X = uniformDisk([0,0],1,Nx);
Y = uniformDisk([0.2,0],1,Ny);
f = rand(size(Y,1),1); % Vector f

% Parameter for rescaling
rMax = rMaxCalc(X,Y);


% Kernel choice:

% G = LogKernel; % G(x) = log(x)

G = Y0Kernel(0.1); % G(x) = Y0(0.1*x) => Bessel decomposition with Robin 
%conditions. 

% G = Y0Kernel(2); % G(x) = Y0(2*x) => Method of rescaling to a root of Y0

% G = Y0Kernel(1000); % G(x) = Y0(1000*x) => Selects frequencies near 0 and 1000

% G = ThinPlate(10,25); % G(x) = 10*x^2*log(25*x)

% G = Kernel(@(r)(exp(-r.^2)),@(r)(-2*r.*exp(-r.^2))); % Arbitrary (smooth)
% kernel

% G = Kernel(@(r)(1./r.^2 ),@(r)(-2./r.^3)); % Arbitrary (singular) kernel


% Choice of the cutoff parameter. 
lambda = 1;
a = lambda/sqrt(sqrt(Nx*Ny)); %this value is of the order of the optimal 
% value for data uniformly distributed in a disk. Choose lambda by trial
% and error to minimize the online time.

tol = 1e-2; % input tolerance

% Offline computations.

[onlineEBD,rq,loc] = offlineEBD(G,X,Y,a,tol); 
% show the radial quadrature : 
rq.show;


% Online procedure.
tic;
q = onlineEBD(f);
time = toc;
fprintf('Online product computed in \n %s seconds \n',num2str(time))

% Error on first entry of q 
dist = sqrt((X(1,1) - Y(:,1)).^2 + (X(1,2) - Y(:,2)).^2);
qval = sum(G.func(dist).*f);
disp('error on first entry');
disp(abs(qval - q(1))/(norm(q,1)));

% Error when f = [1 0 0 ... 0]
dist = sqrt((X(:,1) - Y(1,1)).^2 + (X(:,2) - Y(1,2)).^2);
f = [1; zeros(size(Y,1)-1,1)];
q = onlineEBD(f);
qval = G.func(dist);
fprintf('Linf error for f = [1 0 0 ... 0] \n (effective error / target accuracy) \n');
fprintf('%s / %s \n\n',num2str(max(abs(qval - q))),num2str(tol))
disp('Done');

