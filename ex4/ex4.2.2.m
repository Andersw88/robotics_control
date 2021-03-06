clear all; close all; clc;

format short
SP = model_ABB();

SV = System_Variables(SP);

SV.q = rand(SP.n,1)*10; 
% SV.q = repmat([pi/3],SP.n,1); 
SV = calc_pos(SP,SV); 
Je = zeros(6,SP.n);

pE = fk_e(SP,SV,SP.bN,SP.bP);
pJ = fk_j(SP,SV,[1:SP.n]);

Gravity = [0;0;-9.81];
% Gravity = [0;0;0];
d_time = 0.01;
t = 0:d_time:10;

pid= [50;1;20];

% x = poly3([-1,0,0],[1,0,0],[0,0,0],[0,0,0],t);
x = poly3(pE',[0.3,-0.3,0.3],[0,0,0],[0,0,0],t);


prevVel = zeros(SP.n,1);
pe = zeros(6,1);
G = zeros(6,1);
e = zeros(6,3);
sqErrors = zeros(size(t));
for i = 1:length(t)

    [q, iter] = ik_e(SP,SV,SP.bN,SP.bP,eye(3),x(i,:)');
  
    SV = calc_pos(SP,SV); 
    [HH1,R,M] = calc_hh_a(SP,SV);
    G = calc_G(SP,R,Gravity);

    e2 = q - SV.q;
    e = [e2,G(7:end), (e2- pe(:))/d_time]
    pe = e(:,1);
    
    sqErrors(i) = e(1).^2;
    
    SV.tau = e*pid;
    SV = f_dyn(SP,SV,Gravity);
    SV = int_rk4(SP,SV,d_time,Gravity);

    if mod(i,10) == 0
        cla
        Draw_System(SP, SV, SP.bN, SP.bP,[1:SP.n]);
        axis equal
        grid on
    end
end

figure(2);
plot(t,sqErrors); 
%The maximum errors are most of the time more than an order of magniture less than in 2.1


