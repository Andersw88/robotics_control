clear all; close all; clc;

SP = model_ABB;

tf=10; %final simulation time
dt=5*1e-3; 
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
environment=Environment([R_v v; zeros(1,3) 1]); %environment is instantiated with a homogeneous transform

t=0:dt:tf;
tx = 1;
r = 0.1;
fK = 5000; %fC^2/4 
fC = sqrt(fK*4); %150 If i use params from enviromnet it ends up exactly at 10N
force_d = 10%N
x = bsxfun(@plus,R_v*[ones(size(t))*force_d/fK;cos(t*tx)*r;sin(t*tx)*r],startPos - r*[0,1,0]');
% x = R_v*[ones(size(t))*force_d/fK;cos(t*tx)*r;sin(t*tx)*r];
dx = R_v*[zeros(size(t));-tx*sin(t*tx);tx*cos(t*tx)]*r;
ddx = R_v*[zeros(size(t));-tx^2*cos(t*tx);-tx^2*sin(t*tx)]*r;
w = zeros(size(x));
dw = zeros(size(x));

% rf = - R_v * [10-exp(-t);zeros(size(t));zeros(size(t))];
% drf = - R_v * [exp(-t);zeros(size(t));zeros(size(t))];

A = [0,1;0,0];
B = [0;1];
P = [-30,-20];
pd = place(A,B,P)'

axis([-0.4 1.6 -1.0 1.0 -0.6 1.4]);
pbaspect([1 1 1]);

e = zeros(6,1);
sqErrors = zeros(length(t),3);
forceError = zeros(length(t),1);
orientationError = zeros(length(t),1);
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
    tau_ext=Je(1:3,:)'*f_v;

    if t(i) > 5
       environment.k_e_ = 0;
       environment.c_e_ = 0;
    end
    
    e(1:3) = x(:,i) - pE;
    e(4:6) = R_err(RE, rotd);
    ev = [dx(:,i) - v(1:3); w(:,i) - SV.L(SP.bN).w];
    
    if norm(f_v) > 0
        e(1:3) = inv(R_v)*e(1:3);
        ev(1:3) = inv(R_v)*ev(1:3);
        
        Pf_v = f_v;
%         Je2 = [R_v, zeros(size(R_v)) ; zeros(size(R_v)), R_v]*Je;

        M = 1;
        force_control = H*inv(Je)*[((f_v + R_v*[fC*ev(1) + fK*(e(1) + force_d/fK);0;0])/M);0;0;0];
        
        e(1) = 0;
        e(1:3) = R_v(1)*e(1:3);
        ev(1) = 0;
        ev(1:3) = R_v(1)*ev(1:3);
    else
        force_control = 0;
    end

    SV.tau = H*inv(Je)*([ddx(:,i);dw(:,i)] + [e,ev]*pd - dJe*SV.dq) + N + force_control + tau_ext;
    
    sqErrors(i,:) = sum(e(1:3).^2);
    orientationError(i) = sum(e(4:6).^2);
    forceError(i) = norm(f_v);
    SV = int_rk4(SP,SV,dt,gravity); %call forward dynamics and integrate
    
    if mod(i,10) == 0
        visualizer.update(SP,SV);% update the visualization (needs calc_pos to be called before, which happens in int_rk4 in this example)
    end
end


%Again unexpectedly remove the environment in a separate experiment – how does the
%resulting behavior compare to the hybrid force/position controller? 

%It no longer explodes since it will stop at the trajectory if it doesn't
%feel any force.

skip = 10;
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
