
addpath(genpath(pwd)); % Add folders of the toolbox to the path. 
clear all;%#ok
close all;
clc;

% Parameter def
N = 10^5; 
X = uniformDisk(N); 
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



