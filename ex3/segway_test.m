clear all; close all; clc;

segway=Segway; %create a 2D segway instance
axis([-3,3,-1.5,1.5]); pbaspect([3 1.5 1]); grid on; %adjust the visualization settings

k = (0.99*pi/2)/0.3; % 99% of 90 a degree angle.
% k = 1
% segway.x_=[0; 0; 0.2; 0]; %set the initial state (x=[x; dx; theta; dtheta])
segway.x_=[0.01; 0; 0.3; 0]*k; %set the initial state (x=[x; dx; theta; dtheta])
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

P = [-35.0,-25.0,-0.03,-0.04];
% rank(ctrb(A,B));
if length(A) == rank([B, A*B, A^2*B, A^3*B])
    disp('System is controllable')
end
K = place(A,B,P);



for i=1:length(t)
    x(i) = segway.x_(1);
    theta(i) = segway.x_(3);
    
	segway.u_ = -K*segway.x_; %set the control input at the current time step
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
