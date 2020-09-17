function [ X ] = uniformSegment( a,b,N)


X1 = a + (b-a)*rand(N,1);
X2 = 0*X1;
X = [X1,X2];

end

