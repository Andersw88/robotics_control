clear all; close all; clc;

inverted_pendulum=InvertedPendulum; %create an inverted pendulum instance
axis([-1.2,1.2,-1.2,1.2]); grid on; %adjust the visualization settings

inverted_pendulum.x_=[0.3; 0]; %set the initial state (x=[theta; dtheta])
inverted_pendulum.dt_= 5*1e-3;   %set the sampling rate

%Simulation duration
tf=8;
t=linspace(0,tf,tf/inverted_pendulum.dt_);

angle = zeros(size(t));

g = 9.81;
l = 1;


A = [0,1;g/l,0];
B = [0;1];
P = [-10.0,-5.5];
K = place(A,B,P);

for i=1:length(t)
	tic;
    
    angle(i) = inverted_pendulum.x_(1);
    
%     inverted_pendulum.u_ = pid_(-inverted_pendulum.x_(1),inverted_pendulum.dt_,[30,0.0,10]);
    inverted_pendulum.u_ = -K * inverted_pendulum.x_;
		
	inverted_pendulum.x_=inverted_pendulum.step; %integrate forward according to x_new=f(x,u,dt) and update the state vector
	duration=toc;
	
	pause(inverted_pendulum.dt_-duration); %a crude way of making the visualization appear in real-time
end

figure(2)
plot(t,angle);