
clear all; close all; clc;

format short

load('tripod.mat')
t = tripod.t;
dt = t(2)-t(1);
n = 10;
d = struct('a',25^2/4,'b',25/4,'g',25/3,'tau',2)

s = exp(-d.g/d.tau*tripod.t);

g = 1:10;
f = zeros(length(t),n);

for i = 1:n
    f(:,i) = gaussmf(s, [0.5/n (i-1)/(n-1)]);
end
subplot(2,1,1)
plot(t,f)
subplot(2,1,2)
plot(s,f)