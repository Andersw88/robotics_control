clear all; close all; clc;

format short
SP = model_planar_system_n();

SV = System_Variables(SP);

SV.q = rand(SP.n,1)*10; 

% SV.q = [1.5,0.5,0.0,0.5]'; 
SV = calc_pos(SP,SV); 

Je = zeros(6,SP.n);

pE = fk_e(SP,SV,SP.bN,SP.bP);


Gravity = [0;0;-9.81];
d_time = 0.01;
t = 0:d_time:1;


% x = poly3([-1,0,0],[1,0,0],[0,0,0],[0,0,0],t);
x = poly3(pE',[4,0,0],[0,0,0],[0,0,0],t);

prevVel = zeros(SP.n,1);
for i = 1:length(t)

    [q, iter] = ik_e(SP,SV,SP.bN,SP.bP,eye(3),x(i,:)');
    SV.q = q;
    SV = calc_pos(SP,SV); 
    
    cla
    Draw_System(SP, SV, SP.bN, SP.bP,[1:SP.n]);
    axis equal
    grid on
end


