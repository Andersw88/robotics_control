function [ x ] = euler( dx, x, t)
    x = x + dx(x)*t;
end

