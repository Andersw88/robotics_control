function [ dx ] = transformationSystem(x,s,p)
   dx = [x(2); -p.a/p.tau^2*x(1) - p.b/p.tau*x(2) + forcingFunction(s,p)/p.tau^2 + p.end*p.a/p.tau; 0];
end

