function [ X ] = uniformSine( a,b,N)

X1 = a + (b-a)*rand(N,1);
X2 = sin(X1);

X = [X1,X2];

end

