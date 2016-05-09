function [ x ] = vdp( x )
    a = 0.2;
    A = [0,1;-1,-a*(x(1)^2-1)];
    x = A*x;
end

