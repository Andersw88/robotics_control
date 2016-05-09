
clear all; close all; clc;

format short

load('tripod.mat');

d = tripod.data;
t = tripod.t;


p2 = 9;
p = struct('a',p2^2/4,'b',p2,'g',p2/3,'tau',1,'n',10); %totalError = 0.0085
% p = struct('a',p2^2/4,'b',p2,'g',p2/3,'tau',1,'n',20); %totalError = 0.0090
% p = struct('a',p2^2/4,'b',p2,'g',p2/3,'tau',1,'n',15); %totalError = 0.0092
% p = struct('a',p2^2/4,'b',p2,'g',p2/3,'tau',1,'n',5); %totalError = 0.0378

p.dt = t(2)-t(1);
p.end = d(end,1);
p.start = d(1,1);
p.s = exp(-p.g/p.tau*(0:p.dt:1)');

p.w = learnDMP(d,p);

% p.w

p.tau = 1;
t = 0:p.dt:p.tau;
y = zeros(size(t));
error = zeros(size(t));
x = d(1,:)';
% p.end = 0.2;

for i = 1:length(t)
    s = exp(-p.g/p.tau*t(i)');
    f = @(x)transformationSystem(x,s,p);
%     x = x + transformationSystem(x,s,p)*p.dt;
    x = rk4(f,x,p.dt);
    y(i,:) = x(1);
    if p.tau == 1
        error(i) = norm(x(1)-d(i,1))^2;
    end
end

subplot(2,1,1)
plot(tripod.t,d(:,1)); hold on;
plot(t,y(:,1)); hold on;
legend('real','x');
subplot(2,1,2)
plot(t,error);
legend('error');

totalError = sum(error)



