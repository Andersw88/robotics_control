
clear all; close all; clc;

segway=Segway; %create a 2D segway instance
axis([-3,3,-1.5,1.5]); pbaspect([3 1.5 1]); grid on; %adjust the visualization settings

% k = 1;
k = (0.80*pi/2)/0.3; %Easier to control again
x2 = [0.01; 0; 0.3; 0]*k;
segway.x_= x2; %set the initial state (x=[x; dx; theta; dtheta])
segway.dt_= 2*1e-3;   %set the sampling rate

%Simulation duration
tf=8;
t=linspace(0,tf,tf/segway.dt_);
x = zeros(size(t));
theta = zeros(size(t));

l = segway.l_;
g = segway.g_;
m = segway.m_;
M = segway.M_;

A = [0,1,0,0;
     0,0,m*g/M,0;
     0,0,0,1;
     0,0,g/l + m*g/(l*M),0];
B = [0;l/M;0;1/l*M];
P = [-10.0,-0.2,-12.0,-0.3];
K = place(A,B,P);

C = [1,0,0,0;0,0,1,0];
if length(A) == rank([C; C*A; C*A^2; C*A^3])
    disp('System is observable')
else
    disp('System is not observable')
end
o = obsv(A,C);
PL = [-500.1,-17.2,-800.3,-15.0];
L = place(A',C',PL)';


for i=1:length(t)
    x(i) = segway.x_(1);
    theta(i) = segway.x_(3);
    
    y = C*segway.x_;
    
    u = -K*x2;
    x2 = x2 + segway.dt_ * (A*x2 + B*u + L*(y - C*x2));
    segway.u_ = u; %set the control input at the current time step
	tic;
	segway.x_= segway.step; %integrate forward according to x_new=f(x,u,dt) and update the state vector
	t2=toc;
	
	pause(segway.dt_-t2); %a crude way of making the visualization appear in real-time
end

figure
subplot(2,1,1);
plot(t,x)
subplot(2,1,2);
plot(t,theta)
