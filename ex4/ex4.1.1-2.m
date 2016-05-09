clear all; close all; clc;

format short
% SP = model_ETS7();
SP = model_planar_system_n();

SV = System_Variables(SP);

SV.q = rand(SP.n,1)*10; 
SV = calc_pos(SP,SV); 

% Je = calc_Je(SP,SV,SP.bN,SP.bP);
% n = null(Je);
% SV.dq = n;

Je = zeros(6,SP.n);
Gravity = [0;0;-9.81];
d_time = 0.01;

prevVel = zeros(SP.n,1);
for time=0:d_time:5
  
    
    pE = fk_e(SP,SV,SP.bN,SP.bP);
    pJ = fk_j(SP,SV,[1:SP.n]);
    for i = 1:SP.n
        k = SV.L(i).R(:,3);
        Je(:,i) = [cross(k,pE - pJ(:,i)); k];
    end
    n = null(Je(1:3,:));
    n = n(:,1);
    if dot(prevVel,n) < 0
        n = -n;
    end
    prevVel = n;
    SV.dq = n*5;
    
    SV = int_rk4(SP,SV,d_time,Gravity); % numerical integration

    cla
    Draw_System(SP, SV, SP.bN, SP.bP,[1:SP.n]);
    axis equal
    grid on
end


