clear all; close all; clc;
% example showing the usage of the alternative hgtransform based visualizer
% for in-the-loop visualization (much faster alternative then the
% Draw_System function)

%SP = model_SDH; %test system (SCHUNK SDH hand)
SP = model_ABB;

tf=10; %final simulation time
dt=5*1e-2; %simulation time step
gravity=[0;0;-9.81]; %gravity vector

% generate random initial configuration
SV = System_Variables(SP);
SV.q = rand(size(SV.q));
SV = calc_pos(SP,SV);

%create an instance of the visualizer
visualizer=MBSVisualizer(SP,SV); 
axis([-1.2 1.2 -1.2 1.2 -0.4 1.2]);
pbaspect([1 1 1]);

%the only external system influence is gravity
t=0:dt:tf;
for i=1:length(t)
	tic;
	SV = int_rk4(SP,SV,dt,gravity); % numerical integration
	visualizer.update(SP,SV);% update the visualization (needs calc_pos to be called before, which happens in int_rk4 in this example)
	dt_=toc;
	
	pause(dt-dt_); %crude way of making the simulation appear in real-time
end