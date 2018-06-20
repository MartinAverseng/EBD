function [ X ] = uniformCircle(c,R,N)

theta = 2*pi*rand(N,1);
X = [c(1)+R*cos(theta) c(2)+R*sin(theta)];

end

