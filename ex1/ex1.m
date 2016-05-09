%% Robot Control Exercise 1: Math Basics I
% by Anders Wikström (Vikström) 
% 8806060474
% 2016-04-05

%% Matrix Manipulation
A = [0,1,2,3;4,5,6,7;8,9,10,11;12,13,14,15];

%1
M = [1,0,0,0;0,2,0,0;0,0,1,0;0,0,0,1];
A*M

%2
M = [1,0,0,0;0,0.5,0,0;0,0,1,0;0,0,0,1];
M*A

%3
M = [1,0,0,0;0,0,0,1;0,0,1,0;0,1,0,0];
A*M

%4
M = [1,0,0,0;0,0,1,0;0,1,0,0;0,0,0,1];
M*A

%5
M = [0,0,0,0;0,1,0,0;0,0,1,0;1,0,0,1];
A*M

%6
M = [1,0,0,0;0,0,1,0;0,0,0,1];
M*A

%% Projection 1 - 3

% I - P is the vector minus the projection of the vector. So you get the
%remainder of the vector after projection. (I - P)v = v - Pv

A = [1,2,3];
plotPlane(A,0,5);
x = [9,9,6];

% A2 = A/norm(A)
% proj = A2'*A2;
proj = A'*A / (A*A');
p = proj*x'
p2 = A2*x'*A2 %same as p
p3 = (eye(3) - proj)*x' 
plot3([0,x(1)],[0,x(2)],[0,x(3)])
plot3([0,p(1)],[0,p(2)],[0,p(3)])
plot3([0,p2(1)],[0,p2(2)],[0,p2(3)])
plot3([0,p3(1)],[0,p3(2)],[0,p3(3)])

%% Eigenvalues and Eigenvectors 1

A = [2,1;2,3];

% |2 - l, 1; 2, 3 - l| = 0
% (2 - l)*(3 - l) - 1*2 = 0
% l^2 - 5l + 4 = 0
% l = 1/2* (-5 +- sqrt(5^2 - 16))
l1 = 1/2* (5 + sqrt(5^2 - 16))
l2 = 1/2* (5 - sqrt(5^2 - 16))
% l1 = 4
% l2 = 1
% eig(A)

%(A - l1*I)x = 0 
%([-2,1;2,-1])*x = 0
%([-2,1;2,-1])*[1;2] %is [0,0]
%A*[1;2] %is l1*[1;2]

v1 = [1;2] %First eigenvector

%(A - l1*I)x = 0 
%([1,1;2,2])*x = 0
%([1,1;2,2])*[1;-1] %is [0,0]
%A*[1;-1] %is l2*[1;-1]

v2 = [1;-1] %Second eigenvector

l1 + l2 == trace(A) 
l1 * l2 == det(A) %these are the same

% [V,D] = eig(A)

%% Eigenvalues and Eigenvectors 2
clf
t = 0:pi/2:pi*2;
x = [cos(t);sin(t)];
% M = [3,1;4,-1];
M = [3,1;1,3];
[v,d] = eig(M);
p = M*x;
v2 = [v,-v]
p2 = M*v2

scatter(v2(1,:),v2(2,:)); hold on;
scatter(p2(1,:),p2(2,:))
for i = 1:length(v2)
    plot([v2(1,i),p2(1,i)],[v2(2,i),p2(2,i)]); hold on;
end
for i = 1:length(v2)
    plot([0,v2(1,i)],[0,v2(2,i)]); hold on;
end

scatter(x(1,:),x(2,:)); hold on;
scatter(p(1,:),p(2,:))
% for i = 1:length(x)
%     plot([x(1,i),p(1,i)],[x(2,i),p(2,i)]); hold on;
% end
for i = 1:length(x)
    plot([0,p(1,i)],[0,p(2,i)]); hold on;
end

for i = 1:length(x)
    plot([0,x(1,i)],[0,x(2,i)]); hold on;
end

%The lengths of the vectors are scaled by the eigenvalues and vectors in
%the same direction as the eigenvectors are not rotated.


%% Linear Equations 1
A = [1,2;3,4;5,6];
B = [-1.3;0.9;3.1];

% x = A\B
x2= inv(A'*A)*A'*B
%x = [3.5; -2.4] is the only solution

%% Linear Equations 2
A = [1,2,3,4;5,6,7,8];
B = [1;2]

a = [1,5];
a = a/norm(a);
b = [2,6];
b = b - ((b*a')/(a*a'))*a;
b = b/norm(b);

Q = [a;b]
% Q*Q' %is identity as is should
% Q*inv(Q)
[Q,R] = qr(A)
%Orthonormal basis for R(A)
%     0.1961    0.9806
%     0.9806   -0.1961

R = [1,2,3,4;5,6,7,8];
R(2,:) = R(2,:) - R(1,:)*5;
R(2,:) = R(2,:)/-4;
R(1,:) = R(1,:) - R(2,:)*2;
%x1 = x3 + 2x4
%x2 = -2x3 - 3x4
null = [1,-2,1,0;2,-3,0,1];

a = null(1,:);
a = a/norm(a);
b = null(2,:);
b = b - ((b*a')/(a*a'))*a;
b = b/norm(b);
null = [a;b]
% A*a'
% A*b'
% null2 = null(A)

%Orthonormal basis for N(A)
%     0.4082   -0.8165    0.4082         0
%     0.3651   -0.1826   -0.7303    0.5477

% x = A\B %just a solution
leastNorm = A'*inv(A*A')*B
%least norm solution x = [-0.0500; 0.0250; 0.1000; 0.1750];
% norm(x)

%% Taylor Polynomial Approximation 1
clf
t = 0:0.01:7;
x = zeros(6,length(t));
x(1,:) = 1 + t;
x(2,:) = 1 + t + t.^2/factorial(2);
x(3,:) = 1 + t + t.^2/factorial(2) + t.^3/factorial(3);
x(4,:) = 1 + t + t.^2/factorial(2) + t.^3/factorial(3) + t.^4/factorial(4);
x(5,:) = 1 + t + t.^2/factorial(2) + t.^3/factorial(3) + t.^4/factorial(4) + t.^5/factorial(5);
x(6,:) = exp(t);
plot(t,x)

% syms q;
% for n = 1:5
%     f(n) = taylor(exp(q),'Order', n);
% end
% f(6) = exp(q);
% fplot(f)

%% Taylor Polynomial Approximation 2
[t1,t2] = meshgrid(-2:0.1:2);

x1 = 2*t1.^2 + 3*t2.^2;
u = [-1/4, -1/4];
x2 = 2*u(1).^2 + 3*u(2).^2 + 4*u(1)*(t1 - u(1)) + 6*u(2)*(t2 - u(2));

% test1 = 2*u(1).^2 + 3*u(2).^2 
% test2 = 2*u(1).^2 + 3*u(2).^2 + 4*u(1)*(-1/4 - u(1)) + 6*u(2)*(-1/4 -u(2))

surf(x1); hold on;
surf(x2)
