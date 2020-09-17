
addpath(genpath(pwd)); % Add folders of the toolbox to the path. 
clear all;%#ok
close all;
clc;

%% Small data, validate error

N = 10^3; 
X = uniformDisk([0,0],1,N); 
G = LogKernel; 
a = 1/sqrt(N); 
tol = 1e-3; 

% Offline part
tic
onlineEBD = offlineEBD(G,X,X,a,tol); 
toc

% Online part
f = rand(size(X,1),1);

tic
q = onlineEBD(f);
toc

rXX = sqrt((X(:,1) - X(:,1).').^2 + (X(:,2) - X(:,2).').^2);
Gfull = log(rXX);
Gfull(rXX == 0) = 0;
q_ref = Gfull*f;
err = max(abs(q - q_ref))./norm(f,1);
disp(err);

%%



% Parameter def
N = 10^5; 
X = uniformDisk([0,0],1,N); 
G = LogKernel; 
a = 1/sqrt(N); 
tol = 1e-3; 

% Offline part
tic
onlineEBD = offlineEBD(G,X,X,a,tol); 
toc

% Online part
f = rand(size(X,1),1);

tic
q = onlineEBD(f);
toc



