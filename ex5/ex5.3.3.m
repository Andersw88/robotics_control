clear all; close all; clc;

SP = model_ABB;

tf=10; %final simulation time
dt=5*1e-3; %simulation time step - needs to be fairly small since large contact forces due to stiff environments might be created 
gravity=[0;0;-9.81]; %gravity vector

% generate initial configuration
SV = System_Variables(SP);


R_v=ry(pi/4);
rotd = R_v*[0,0,1;0,1,0;1,0,0];
startPos = rotd*[-0.02,0,0]';
[q, iter] = ik_e(SP,SV,SP.bN,SP.bP,inv(SP.bR),startPos,rotd);
SV.q = q;
SV = calc_pos(SP,SV); %need to call calc_pos for the visualizer


%create an instance of the multi body system visualizer
visualizer=MBSVisualizer(SP,SV);

%create an instance of a (planar) environment
v=[0.6;0;0.6]; %envrionment position
 %environment orientation (x-axis points along environment plane normal)
R_v=ry(pi/4);
environment=Environment([R_v v; zeros(1,3) 1]); %environment is instantiated with a homogeneous transform


t=0:dt:tf;
tx = 1;
r = 0.1;
% x = R_v*[zeros(size(t));cos(t*tx);sin(t*tx)]*r;
x = bsxfun(@plus,R_v*[zeros(size(t));cos(t*tx)*r;sin(t*tx)*r],startPos - r*[0,1,0]');
dx = R_v*[zeros(size(t));-tx*sin(t*tx);tx*cos(t*tx)]*r;
ddx = R_v*[zeros(size(t));-tx^2*cos(t*tx);-tx^2*sin(t*tx)]*r;
rf = - R_v * [10-exp(-t);zeros(size(t));zeros(size(t))];
drf = - R_v * [exp(-t);zeros(size(t));zeros(size(t))];
% rf = [10-exp(-t);zeros(size(t));zeros(size(t))];

A = [0,1;0,0];
B = [0;1];
P = [-30,-20];
pd = place(A,B,P)'
% pd = [0.5,0.1];

% pid_f = [2,1,1]';
pid_f = [2,2,0.004]'*5000;

%visualization tweaks
axis([-0.4 1.6 -1.0 1.0 -0.6 1.4]);
pbaspect([1 1 1]);

e = zeros(6,1);
sqErrors = zeros(length(t),1);
forceError = zeros(length(t),1);
orientationError = zeros(length(t),1);
dfError = zeros(length(t),1);
I = zeros(6,1);
Pf_v = 0;
for i=1:length(t)
    %calculate the end-effector Jacobian
    SV = calc_pos(SP,SV);  
    SV = calc_vel(SP,SV);
    Je = calc_Je(SP,SV,SP.bN,SP.bP);
	
	%calculate end-effector position 
    [pE,RE] = fk_e(SP,SV,SP.bN,SP.bP);
		
	dJe = calc_dJe(SP,SV,SP.bN,SP.bP);
    v=Je*SV.dq;
	
    H = calc_hh(SP,SV); H=H(7:end,7:end);
    N = r_ne(SP,SV,gravity); N=N(7:end);
    
    [f_v df_v]=getContactForce(environment,pE,v(1:3),dt); 
    
%     if t(i) > 5
%        environment.k_e_ = 0;
%        environment.c_e_ = 0;
%     end
    
    fe = ([rf(:,i);0;0;0] - [f_v;0;0;0]);
    fde = ([drf(:,i);0;0;0] - [df_v;0;0;0]);
    force_control = H*inv(Je)*(-[fe,I,fde]*pid_f/5000);
    
    I = I + fe*dt;
    
    e(1:3) = x(:,i) - pE;
    e(4:6) = R_err(RE, rotd);
    ev = [dx(:,i) - v(1:3); - v(4:6)];
    
    e(1:3) = inv(R_v)*e(1:3);
    e(1) = 0;
    e(1:3) = R_v(1)*e(1:3);
    ev(1:3) = inv(R_v)*ev(1:3);
    ev(1) = 0;
    ev(1:3) = R_v(1)*ev(1:3);
    
    
    SV.tau = H*inv(Je)*([ddx(:,i);0;0;0] + [e,ev]*pd - dJe*SV.dq) + force_control + N;
%     SV.tau = zeros(6,1);
  
    sqErrors(i) = sum(e(1:3).^2);
    orientationError(i) = sum(e(4:6).^2);
    forceError(i) = sum(fe.^2);
%     dfError(i) = sum(fde.^2);
    
    
    SV = int_rk4(SP,SV,dt,gravity); %call forward dynamics and integrate
    
    if mod(i,10) == 0
        visualizer.update(SP,SV);% update the visualization (needs calc_pos to be called before, which happens in int_rk4 in this example)
    end
end

%Augment the controller with an integrator – do you notice a sizable performance improvement?
%The force takes longer to settle after the initial errors when it hit the plane.

skip = 100;
figure(2);
subplot(3,1,1);
plot(t(skip:end),sqErrors(skip:end,:)); 
ylabel('position');
subplot(3,1,2);
plot(t(skip:end),forceError(skip:end)); 
ylabel('force');
subplot(3,1,3);
plot(t(skip:end),orientationError(skip:end)); 
ylabel('orientation');
