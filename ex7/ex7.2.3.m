
clear all; close all; clc;
format short


SP = model_ABB;
SV = System_Variables(SP);
[q, iter] = ik_e(SP,SV,SP.bN,SP.bP,SP.bR,[0.0,0.0,0.7]',eye(3));
SV.q = q;
SV = calc_pos(SP,SV);

n1 = [0;1/sqrt(2);1/sqrt(2)];
n2 = [0;0;1];
a = [0.3,0.25];
tf = 10;
dt=5*1e-3; 
t=0:dt:tf;

r1 = vrrotvec2mat(vrrotvec([1,0,0]',n1));
r2 = vrrotvec2mat(vrrotvec([1,0,0]',n2));

visualizer=MBSVisualizer(SP,SV);
% plotPlane(n1',a(1),1); hold on;
% plotPlane(n2',a(2),1); hold on;

environment=Environment([r1, r1*[a(1),0,0]'; zeros(1,3) 1]);
environment=Environment([r2, r2*[a(2),0,0]'; zeros(1,3) 1]);
visualizer.update(SP,SV);

p = 0.001;
view([90,0])
axis equal

opts = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
% opts = optimoptions('quadprog','Algorithm','active-set','Display','off');

dq = zeros(6,1);
H = [eye(6)*p,zeros(6,1);zeros(1,6),1];
f = zeros(length(H),1);
for i=1:length(t)
    SV = calc_pos(SP,SV);  
    SV = calc_vel(SP,SV);
    Je = calc_Je(SP,SV,SP.bN,SP.bP);
    [pE,RE] = fk_e(SP,SV,SP.bN,SP.bP);	

    Jet = n1'*Je(1:3,:);
    Jet2 = n2'*Je(1:3,:);

    e1 = n1'*pE - a(1);
    e2 = n2'*pE - a(2);
    v = Je*SV.dq;
    
    A1 = [Jet,-1];
    b1 = -e1;
    
    x = quadprog(H,f,[],[],A1,b1,[],[],[],opts);
    w = x(7);

    A1 = [-Jet,0];
    b1 = e1 - w;
    A2 = [Jet2,-1];
    b2 = -e2 ;

    x = quadprog(H,f,A1,b1,A2,b2,[],[],[],opts);
    dq = x(1:6);
    SV.dq = dq;
    
    SV = int_rk4(SP,SV,dt,[0,0,0]');
    if mod(i,10) == 0
        visualizer.update(SP,SV);
    end
    
end

