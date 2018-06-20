function [ X ] = uniformDisk( center,R,N)


x = center(1);
y = center(2);
L = 2*R;
X1 = L*(rand(fix(4*N/pi*1.2),1)-1/2);
X2 = L*(rand(fix(4*N/pi*1.2),1)-1/2);

r = sqrt(X1.^2 + X2.^2);
X1 = X1(r<R);
X2 = X2(r<R);

while length(X1) < N
    x1 = 2*L*(rand(fix(4*N/pi*1.2),1)-1/2);
    x2 = 2*L*(rand(fix(4*N/pi*1.2),1)-1/2);
    r = sqrt(x1.^2 + x2.^2);
    x1 = x1(r<R);
    x2 = x2(r<R);
    X1 = [X1;x1];
    X2 = [X2;x2]; 
end

X = [X1(1:N)+x X2(1:N)+y];

end

