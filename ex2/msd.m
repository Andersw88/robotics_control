function [ x ] = msd( x, k, m, c )
%     k = 0.5;
%     m = 1;
%     c = 0;
    A = [0,1;-k/m,-c/m];
    x = A*x;
end

