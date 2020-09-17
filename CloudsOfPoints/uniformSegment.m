function [ X ] = uniformSegment( a,b,N)
% X  = uniformDisk( center,R,N)
% inputs : c = [c1,c2] center, R radius, N number of points 
% output : X of size Nx2 cloud of points uniformly distributed on the disk.


X1 = a + (b-a)*rand(N,1);
X2 = sin(X1);

X = [X1,X2];

end

