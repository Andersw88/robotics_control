
clear all; close all; clc;

format short
SP = model_planar_system_n;
load('letter.mat');

% d = letter;
x1 = letter.x1.data/200;
x2 = letter.x2.data/200;
t = letter.t;
p2 = 8;

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
 
visualizer=MBSVisualizer(SP,SV);
offset = [0.0,0.0]
for j = 1:5
    offset = offset + 0.02;
    p1.tau = 1;
    p2.tau = p1.tau;
    t = 0:p1.dt:p1.tau;
    p1.end = x1(end,1)+offset(1);
    p2.end = x2(end,1)+offset(2);
    p1.scale = (p1.start - x1(end,1))/(p1.start - x1(end,1)+offset(1)); %Improvised to keep the size of the motion similar.
    p2.scale = (p2.start - x2(end,1))/(p2.start - x2(end,1)+offset(2));
    y1 = x1(1,1:3)';
    y2 = x2(1,1:3)';
    yt = zeros(length(t),2);
    error = zeros(size(t));
    

    for i = 1:length(t)

        s = exp(-p1.g/p1.tau*t(i)');
        f1 = @(x)transformationSystem(x,s,p1);
        f2 = @(x)transformationSystem(x,s,p2);
        y1 = rk4(f1,y1,p1.dt);
        y2 = rk4(f2,y2,p2.dt);

        [q, iter] = ik_e(SP,SV,SP.bN,SP.bP,eye(3),[y1(1),y2(1),0]');
        SV.q = q;

        if p1.tau == 1
            error(i) = norm([y1(1)-x1(i,1),y2(2)-x2(i,1)])^2;
        end

        yt(i,:) = [y1(1),y2(1)];
        if mod(i,10) == 0
            SV = int_rk4(SP,SV,0,[0,0,0]'); %to make the visualizer work
            visualizer.update(SP,SV);
        end
    end

    figure(2);
    subplot(3,2,1)
    plot(t,yt(:,1)); hold on;
    subplot(3,2,3)
    plot(t,yt(:,2)); hold on;
    subplot(3,2,[2,4,6])
    plot(yt(:,1),yt(:,2)); hold on;
    axis equal
    subplot(3,2,5)
    plot(t,error); hold on;
    legend('error');
end

figure(2);
subplot(3,2,1)
plot(letter.t,x1(:,1)); hold on;
subplot(3,2,3)
plot(letter.t,x2(:,1)); hold on;
subplot(3,2,[2,4,6])
plot(x1(:,1),x2(:,1)); hold on;



