
clear all; close all; clc;
format short


SP = model_ABB;
SV = System_Variables(SP);
SV.q = ones(6,1)*-pi/5;
[q, iter] = ik_e(SP,SV,SP.bN,SP.bP,inv(SP.bR),[0,0,0]');
% SV.q = ones(6,1).*rand(6,1)*10;
SV.q = q;
SV = calc_pos(SP,SV);

n1 = [0;1/sqrt(2);1/sqrt(2)];
n2 = [0;0;1];
a = [0.3,0.25];
tf = 100;
dt=5*1e-3; 
t=0:dt:tf;

rt = vrrotvec2mat(vrrotvec([0,0,1]',n1))

visualizer=MBSVisualizer(SP,SV);
plotPlane(n1',ones(1,3)*a(1)/2,1); hold on;
% plotPlane(n1',zeros(1,3)*a(1),1); hold on;
% plotPlane(n2',ones(1,3)*a(2),1); hold on;
r1 = [0,1/sqrt(2),1/sqrt(2);-1/sqrt(2),1/sqrt(2),0;-1/sqrt(2),0,1/sqrt(2)];

visualizer.update(SP,SV);

p = 0.1;

opts = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
% opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
 
for i=1:length(t)
    SV = calc_pos(SP,SV);  
    SV = calc_vel(SP,SV);
    Je = calc_Je(SP,SV,SP.bN,SP.bP);
    [pE,RE] = fk_e(SP,SV,SP.bN,SP.bP);	
% 	dJe = calc_dJe(SP,SV,SP.bN,SP.bP);
    
%     v = Je*SV.dq;
%     Jet = diag([n1;n1])*Je;
    
    Jet = diag(n1)*Je(1:3,1:6);
    Jet2 = diag([n1;ones(3,1)])*Je;
%     Jet2 = diag([n1])*Je(1:3,1:3);
    
%     Jet = [n1*Je(1:3,1:3),zeros(3);zeros(3),n1*Je(1:3)];
%     Jet = [n1.*Je(1:3,1:3),zeros(3);zeros(3),n1.*Je(1:3)];
    
%     e1 = -Je*SV.q.*SV.dq.*[n1;n1];

    e1 = n1*(n1'*pE - a(1));
    
%     v = Jet2*SV.dq;
    v = Je*SV.dq;

%     b = Je*SV.dq + [zeros(3,1);v(4:6)];
    H = [eye(3)*p,zeros(3);zeros(3),eye(3)]/2;
    f = zeros(length(H),1);

    
    A = Jet2;
    w = R_err(RE, rt);
    b = [-e1 ; w];
%     b = [-e1 ; 0;0;0];
%     b = -e1 + v(4:6);
%     b = Jet*SV.dq 
    x = quadprog(H,f,[],[],A,b,[],[],[],opts);
    
%     SV.dq = pinv(Jet)*e1;
    SV.dq = x;
    
%     [x,SV.dq]
%     SV.dq = pinv(Jet)*[de;0;0;0];
%     SV.dq = inv(Je)*[de;0;0;0];
    SV = int_rk4(SP,SV,dt,[0,0,0]');
    visualizer.update(SP,SV);
    
end


% 
% % options = optimset('Algorithm','active-set');
% [x,FVAL] = quadprog(H,f,A,b);
% % [x,FVAL] =  quadprog(H,f,A,b,[],[],[],[],[],options);
% 
% area = [-2,10];
% X1 = linspace(area(1),area(2));
% X2 = linspace(area(1),area(2));
% [X,Y] = meshgrid(X1,X2);
% Z = sqrt(X.^2 + 4*(Y.^2 - 4).^2);
% % Z = (X.^2 + 4*(Y.^2 - 4).^2);
% contour(X,Y,Z,20); hold on;
% 
% grad = [2*x(1), 16*x(2) - 32;1,1;-1,2];
% % for i = 1:length(grad)
% %     grad(i,:) = grad(i,:)/norm(grad(i,:));
% % end
% for i = 1:length(grad)
%     grad(i,:) = grad(i,:)/sqrt(norm(grad(i,:)));
% end
% % x11 = -X1 + 7;
% % x12 = (X1 + 4)/2;
% plot(X1,-X1 + 7)
% plot(X1,(X1 + 4)/2)
% 
% for i = 1:length(grad)
%     plot([x(1),x(1) + grad(i,1)],[x(2),x(2) + grad(i,2)])
% end
% 
% legend('contour','c1: -x1 - x2 = -7','c2: x1 - 2x2 = -4','grad','grad c1','grad c2')
% % quiver(x(1),x(2),grad(1,1), grad(1,2))
% % quiver(x(1),x(2),grad(2,1), grad(2,2))
% % quiver(x(1),x(2),grad(3,1), grad(3,2))
% 
% axis equal;
% xlim(area);
% ylim(area);



