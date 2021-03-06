clear all; close all; clc;

format short
SP = model_ABB();
load('traj.mat')

SV = System_Variables(SP);

SV.q = traj.q_i;

SV = calc_pos(SP,SV); 
Je = zeros(6,SP.n);

pE = fk_e(SP,SV,SP.bN,SP.bP);

Gravity = [0;0;-9.81];
d_time = traj.dt;
t = (0:size(traj.p_d,2)-1)*d_time;

A = [0,1;0,0];
B = [0;1];
P = [-10,-8];
pd = place(A,B,P);

% pd= [repmat([5,5],3,1);repmat([1,1],3,1)];
% pd = K(1)


visualizer=MBSVisualizer(SP,SV);

% pe = zeros(6,1);
pPos = zeros(3,1);
pRE = zeros(3,1);
e = zeros(6,1);
ev = zeros(6,1);
sqErrors = zeros(length(t),3);

for i = 1:size(traj.p_d,2)
    
    SV = calc_pos(SP,SV);  
    SV = calc_vel(SP,SV);
    HH = calc_hh(SP,SV); H = HH(7:end,7:end);
    N = r_ne(SP,SV,Gravity); N=N(7:end);
    Je = calc_Je(SP,SV,SP.bN,SP.bP);
    dJe = calc_dJe(SP,SV,SP.bN,SP.bP);
    [pE,RE] = fk_e(SP,SV,SP.bN,SP.bP);
    
    e(1:3) = traj.p_d(:,i) - pE;
    e(4:6) = R_err(RE, traj.R_d(:,(1+(i-1)*3):((1+(i-1)*3)+2)));
    
    v = Je*SV.dq;
%     ev = [traj.dp_d(:,i) - SV.L(SP.bN).v;  traj.w_d(:,i) - SV.L(SP.bN).w];
%     ev = [traj.dp_d(:,i) - (pE - pPos)/d_time ;  traj.w_d(:,i) - (e(4:6) - pRE)/d_time];
    ev = [traj.dp_d(:,i) - v(1:3);  traj.w_d(:,i) - SV.L(SP.bN).w];

    

    SV.tau = H*inv(Je)*([traj.ddp_d(:,i);traj.dw_d(:,i)] + pd(:,1).*e + pd(:,2).*ev - dJe*SV.dq) + N;
%     [e,SV.tau]
%     pe = e;
%     pPos = pE;
%     pRE = e(4:6);
    sqErrors(i,:) = sum(e(1:3).^2);

    
    SV = int_rk4(SP,SV,d_time,Gravity);

    if mod(i,10) == 0
        visualizer.update(SP,SV);
    end

end

figure(2);
plot(t,sqErrors); 
