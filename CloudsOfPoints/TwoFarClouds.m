function[X,Y,V] = twoFarClouds(M,N)
    X = randn(M,2);
    Y = randn(N,2) + 10;
    V = randn(N,1);
    V = V/norm(V,1);
end