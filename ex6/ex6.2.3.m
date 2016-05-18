
clear all; close all; clc;

format short
SP = model_planar_system_n;
load('letter.mat');

% d = letter;
x1 = letter.x1.data/200;
x2 = letter.x2.data/200;
t = letter.t;
p2 = 8;
Gravity = [0;0;-9.81];

p1 = struct('a',p2^2/4,'b',p2,'g',p2/3,'tau',1,'n',20,'dt',t(2)-t(1),'start',x1(1,1),'end',x1(end,1));
p1.s = exp(-p1.g/p1.tau*(0:p1.dt:1)');
p2 = p1;
p2.start = x2(1,1);
p2.end = x2(end,1);

SV = System_Variables(SP);
[q, iter] = ik_e(SP,SV,SP.bN,SP.bP,eye(3),[p1.start,p2.start,0]');
SV.q = q;
SV = calc_pos(SP,SV); 

p1.w = learnDMP(x1,p1);
p2.w = learnDMP(x2,p2);
 
offset = [-0.4,0.4]
p1.tau = 1;
p2.tau = p1.tau;
t = 0:p1.dt:p1.tau;
p1.end = x1(end,1)+offset(1);
p2.end = x2(end,1)+offset(2);
p1.scale = (p1.start - x1(end,1))/(p1.start - x1(end,1)+offset(1)); %Improvised to keep the size of the motion similar.
p2.scale = (p2.start - x2(end,1))/(p2.start - x2(end,1)+offset(2));

y1 = x1(1,1:3)';
y2 = x2(1,1:3)';
yt = zeros(length(t),3);
error = zeros(size(t));
visualizer=MBSVisualizer(SP,SV);
e = zeros(6,1);
ev = zeros(6,1);
a = zeros(6,1);
pd = [100,50];

for i = 1:length(t)
    
    s = exp(-p1.g/p1.tau*t(i)');
    f1 = @(x)transformationSystem(x,s,p1);
    f2 = @(x)transformationSystem(x,s,p2);
    y1 = rk4(f1,y1,p1.dt);
    y2 = rk4(f2,y2,p2.dt);
%     
    
    HH = calc_hh(SP,SV); H = HH(7:end,7:end);
    N = r_ne(SP,SV,Gravity); N=N(7:end);
    Je = calc_Je(SP,SV,SP.bN,SP.bP);
    dJe = calc_dJe(SP,SV,SP.bN,SP.bP);
    [pE,RE] = fk_e(SP,SV,SP.bN,SP.bP);

    e(1:2) = [y1(1);y2(1)] - pE(1:2);
    v = Je*SV.dq;
    ev(1:2) = [y1(2);y2(2)] - v(1:2);
    a(1:2) = [y1(3);y2(3)];
    
%     [e(1:3),ev(1:3)]
    SV.tau = H*pinv(Je)*(a + pd(:,1).*e + pd(:,2).*ev - dJe*SV.dq) + N;
    
    if p1.tau == 1
        error(i) = norm([pE(1)-x1(i,1),pE(2)-x2(i,1)])^2;
    end
    yt(i,:) = pE;
    SV = int_rk4(SP,SV,p1.dt,Gravity);
    visualizer.update(SP,SV);
end

figure(2);
subplot(3,2,1)
plot(letter.t,x1(:,1)); hold on;
plot(t,yt(:,1)); hold on;
legend('real','x');
subplot(3,2,3)
plot(letter.t,x2(:,1)); hold on;
plot(t,yt(:,2)); hold on;
legend('real','x');
subplot(3,2,[2,4,6])
plot(x1(:,1),x2(:,1)); hold on;
plot(yt(:,1),yt(:,2)); hold on;
axis square

subplot(3,2,5)
plot(t,error);
legend('error');

totalError = sum(error)


