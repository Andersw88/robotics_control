
%% Mass-Spring-Damper 1
clf;

t = 0.2;
x = zeros(2,20/t);
x(:,1) = [10;0];

k = 1;
m = 1;
c = 4;
f1 = @(x)msd(x,k,m,c);

for n = 1:length(x)
    x(:,n+1) =euler(f1,x(:,n),t);
end
plot(x(1,:)), hold on;

for n = 1:length(x)
    x(:,n+1) =rk4(f1,x(:,n),t);
end
plot(x(1,:)), hold on;

k = 1;
m = 1;
c = 0.02;
f1 = @(x)msd(x,k,m,c);

for n = 1:length(x)
    x(:,n+1) =euler(f1,x(:,n),t);
end
plot(x(1,:)), hold on;

for n = 1:length(x)
    x(:,n+1) =rk4(f1,x(:,n),t);
end
plot(x(1,:)), hold on;

legend('large c/k euler','large c/k rk4','small c/k euler','small c/k rk4','Location','northwest')

%For a small c/k you get occilations around 0, with euler being unstable.
%For a large c/k you get a slope towards 0

%% Analytic vs. Numeric Solution 2
clf;
t = 0.25;
t2 = 0:t:10-t;
x = zeros(10/t,1);

for n = 1:length(x)
    x(n+1) = euler(@f2,x(n),t);
end
plot(x), hold on;

for n = 1:length(x)
    x(n+1) = rk4(@f2,x(n),t);
end

plot(x), hold on;

a = 1;
r = 1;

x = -exp(-t2) + 1;
plot(x,'o'), hold on;

legend('euler','rk4','closed form','Location','southeast')

%Rk4 is very good while euler overestimates the value due to not taking into account the decreasing derivative during each timestep.


%% Van der Pol Oscillator 3 - 2
clf;
clear all;
a = 1000;
f = @(t,x)vdpODE(t,x,a);

[t,y] = ode45(f,[0,1000],[2,0]);
plot(t,y(:,1)), hold on;

% [t,y] = ode15s(f,[0,100],[2,0]);
% plot(t,y(:,1)), hold on;

a = 0.2;
f = @(t,x)vdpODE(t,x,a);

[t,y] = ode45(f,[0,100],[2,0]);
plot(t,y(:,1)), hold on;

% [t,y] = ode15s(f,[0,100],[2,0]);
% plot(t,y(:,1)), hold on;

legend('a = 1000','a = 0.2','Location','southeast')
xlabel('t','Interpreter','latex')
ylabel('x','Interpreter','latex')

%For a = 1000, velocity is almost 0 but the acceleration slowly increase as
%the velocity increase. Then the acceleration suddenly change direction 
%as x becomes negative.
%For a = 0.2 the oscillator almost behaves like a sine function with amplitude 2.

%% Van der Pol Oscillator 3 - 3
clf;
clear all;

[X,Y] = meshgrid(-10:1:10);
X = [X(:), Y(:)]';
dx = zeros(size(X));


a2 = [0.2,1,5];

for k = 1:length(a2)
    a = a2(k);
    f = @(t,x)vdpODE(t,x,a);
    for i = 1:length(X)
        dx(:,i) = f(0,X(:,i));
    end
    figure(k)
    scale = 0.007/a;
    
    quiver(X(1,:),X(2,:),dx(1,:)*scale,dx(2,:)*scale, 'AutoScale','off'), hold on;
    xlabel('x','Interpreter','latex')
    ylabel('$\dot{x}$','Interpreter','latex')
end

[X,Y] = meshgrid(-10:0.5:10);
for k = 1:length(a2)
    a = a2(k);
    figure(k+3)
    for i = 1:numel(X)
        f = @(t,x)vdpODE(t,x,a);
        [t,y] = ode45(f,[0,0.05],[X(i),Y(i)]);
        plot(y(:,1),y(:,2)), hold on;
    end
end

[X,Y] = meshgrid(-10:0.2:10);
for k = 1:length(a2)
    a = a2(k);
    figure(k+6)
    for i = 1:numel(X)
        if  X(i) == -10 || X(i) == 10 || Y(i) == -10 || Y(i) == 10
            f = @(t,x)vdpODE(t,x,a);
            [t,y] = ode45(f,[0,10],[X(i),Y(i)]);
            plot(y(:,1),y(:,2)), hold on;
        end
    end
end