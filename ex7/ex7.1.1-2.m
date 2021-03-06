
clear all; close all; clc;


% min x1, x2
% x1^2 + 4(x2 ? 4)^2 = x1^2 + 4(x2^2 - 8x2 + 16) = x1^2 + 4x2^2 - 32x2 + 64
% subject to 
% - x1 - x2 <= -7 
% x1 - 2x2 <= -4

format short
H = [1,0;0,4]*2;
f = [0,-32];
A = [1,1;-1,2];
b = [7,4];
% b = [5,4];

% options = optimset('Algorithm','active-set');
[x,FVAL] = quadprog(H,f,A,b);
% [x,FVAL] =  quadprog(H,f,A,b,[],[],[],[],[],options);

area = [-0.5,1.1]*2;
X1 = linspace(area(1) + x(1),area(2) + x(2));
X2 = linspace(area(1) + x(1),area(2) + x(2));
[X,Y] = meshgrid(X1,X2);
Z = (X.^2 + 4*(Y - 4).^2);
contour(X,Y,Z,20); hold on;
grad = [2*x(1), 8*x(2) - 32;A(1,:);A(2,:)];
for i = 1:length(grad)
    grad(i,:) = grad(i,:)/norm(grad(i,:));
end
% for i = 1:length(grad)
%     grad(i,:) = grad(i,:)/sqrt(norm(grad(i,:)));
% end
plot(X1,-X1 + b(1))
plot(X1,(X1 + b(2))/2)

for i = 1:length(grad)
    plot([x(1),x(1) + grad(i,1)],[x(2),x(2) + grad(i,2)])
end

legend('contour','c1: -x1 - x2 = -7','c2: x1 - 2x2 = -4','grad','grad c1','grad c2')
% quiver(x(1),x(2),grad(1,1), grad(1,2))
% quiver(x(1),x(2),grad(2,1), grad(2,2))
% quiver(x(1),x(2),grad(3,1), grad(3,2))

axis equal;
xlim(area + x(1));
ylim(area + x(2));

%The objective gradient is a linear combination of the active contraints.


