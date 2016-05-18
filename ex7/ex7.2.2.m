
clear all; close all; clc;
format short


SP = model_ABB;
SV = System_Variables(SP);
SV.q = ones(6,1)*-pi/5;
[q, iter] = ik_e(SP,SV,SP.bN,SP.bP,inv(SP.bR),[0,0.5,0.5]');
SV.q = q;
SV = calc_pos(SP,SV);

n1 = [0;1/sqrt(2);1/sqrt(2)];
n2 = [0;0;1];
a = [2+0.3,0.25];
tf = 10;
dt=5*1e-3; 
t=0:dt:tf;

visualizer=MBSVisualizer(SP,SV);
plotPlane(n1',a(1),1); hold on;
% plotPlane(n2',a(2),1); hold on;
r1 = [0,1/sqrt(2),1/sqrt(2);-1/sqrt(2),1/sqrt(2),0;-1/sqrt(2),0,1/sqrt(2)];

visualizer.update(SP,SV);

p = 0.05;
view([90,0])
axis equal

opts = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
H = [eye(6)*p,zeros(6,1);zeros(1,6),1];
f = zeros(length(H),1);
for i=1:length(t)
    SV = calc_pos(SP,SV);  
    SV = calc_vel(SP,SV);
    Je = calc_Je(SP,SV,SP.bN,SP.bP);
    [pE,RE] = fk_e(SP,SV,SP.bN,SP.bP);	

    Jet = n1'*Je(1:3,:);
    e1 = n1'*pE - a(1);
    
    A = [Jet,1];
    b = -e1;

    x = quadprog(H,f,[],[],A,b,[],[],[],opts);
    SV.dq = x(1:6);
    
    SV = int_rk4(SP,SV,dt,[0,0,0]');
    visualizer.update(SP,SV); 
end

%You get numberical instabilities if p is to low.

