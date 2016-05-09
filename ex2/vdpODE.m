function [ x ] = vdpODE(t, x, a)
    x  = [x(2,:); a*x(2,:).*(1 - x(1,:).^2) - x(1,:)];
end