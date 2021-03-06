function [ x ] = euler( dx, x, t)
    k1 = dx(x);
    k2 = dx(x+k1*t/2);
    k3 = dx(x+k2*t/2);
    k4 = dx(x+k3*t);
    x = x + (k1+2*k2+2*k3+k4)*t/6;
end

