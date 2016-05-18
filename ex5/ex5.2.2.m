clear all; close all; clc;

format short
SP = model_ABB();
load('traj.mat')

SV = System_Variables(SP);

SV.q = traj.q_i;
SV = calc_pos(SP,SV); 
% Je = zeros(6,SP.n);

pE = fk_e(SP,SV,SP.bN,SP.bP);

Gravity = [0;0;-9.81];
d_time = traj.dt;
t = (0:size(traj.p_d,2)-1)*d_time;

A = [0,1;0,0];
B = [0;1];
P = [-10,-6];
pd = place(A,B,P)

visualizer=MBSVisualizer(SP,SV);

e = zeros(6,1);
ev = zeros(6,1);
sqErrors = zeros(length(t),3);

Je = calc_Je(SP,SV,SP.bN,SP.bP);
JeInv = inv(Je);

for i = 1:size(traj.p_d,2)
    
    SV = calc_pos(SP,SV);  
    SV = calc_vel(SP,SV);

    Je = calc_Je(SP,SV,SP.bN,SP.bP);
    
    [HH1,R,M] = calc_hh_a(SP,SV);
    G = calc_G(SP,R,Gravity);
    [pE,RE] = fk_e(SP,SV,SP.bN,SP.bP);
    

    v = Je*SV.dq;
    e(1:3) = traj.p_d(:,i) - pE;
    e(4:6) = R_err(RE, traj.R_d(:,(1+(i-1)*3):((1+(i-1)*3)+2)));
    ev = [traj.dp_d(:,i) - v(1:3);  traj.w_d(:,i) - SV.L(SP.bN).w];
    
    SV.tau = JeInv*(pd(:,1).*e + pd(:,2).*ev) + G(7:end);
%     SV.tau = pd(:,1).*Je*e  + pd(:,2).*Je*ev + G(7:end);
    sqErrors(i,:) = sum(e(1:3).^2);

    
    SV = int_rk4(SP,SV,d_time,Gravity);

    if mod(i,10) == 0
        visualizer.update(SP,SV);
    end

end

figure(2);
plot(t,sqErrors); 
