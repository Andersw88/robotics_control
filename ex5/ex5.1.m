clear all; close all; clc;

format short
SP = model_ABB();

SV = System_Variables(SP);

SV.q = mod((rand(SP.n,1)*10),pi); 
% SV.q = repmat([pi/3],SP.n,1); 
SV = calc_pos(SP,SV); 
Je = zeros(6,SP.n);

pE = fk_e(SP,SV,SP.bN,SP.bP);

Gravity = [0;0;-9.81];
d_time = 0.01;
t = 0:d_time:10;

pd= [200;200];

[q,dq,ddq] = poly3(SV.q',-SV.q',zeros(SP.n,1)',zeros(SP.n,1)',t);

q = q';
dq = dq'
ddq = ddq';
visualizer=MBSVisualizer(SP,SV);

pe = zeros(6,1);
G = zeros(6,1);
e = zeros(6,3);
sqErrors = zeros(size(t));
for i = 1:length(t)
    

%     [HH1,R,M] = calc_hh_a(SP,SV);
    HH = calc_hh(SP,SV);
    H = HH(7:end,7:end);
    N = r_ne(SP,SV,Gravity); N=N(7:end);

    e = q(:,i) - SV.q
    de = dq(:,i) - SV.dq;
    SV.tau = H*(ddq(:,i) + pd(1)*e + pd(2)*de) + N;
    pe = e;
    
    sqErrors(i) = sum(e.^2);
    
    SV = int_rk4(SP,SV,d_time,Gravity);

    if mod(i,10) == 0
        visualizer.update(SP,SV);
    end

end

figure(2);
plot(t,sqErrors); 
