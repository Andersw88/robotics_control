function [T, Y] = rk4(file_name, tspan, y, h, varargin)
%
%[T, Y] = rk4(file_name, tspan, y, h, varargin)
%
% Fixed step Runge-Kutta 4th order method
%
% Input:
% ------
% file_name     - name of a file containing the equations of motion (example: @msd)
% tspan         - a vector containing the initial and final time
% y             - initial state
%                 state = [generalized positions ; generalized velocities]
% h             - step size
%
% Output:
% ------
% T             - sequence of points in time where the state of the
%                 system was evaluated
% Y             - state history
%
% MC, Exercise 1 [2009/09/07]

n_it = (tspan(2) - tspan(1))/h + 1;  % determine number of iterations

% Memory allocation for the output
T = zeros(n_it,1);
Y = zeros(n_it,length(y));

t = tspan(1);

y = y(:);  % make sure that y is a column vector

Y(1,:) = y';
T(1)   = t;

for i = 2:1:n_it

  k1 = feval(file_name , t       , y          , varargin{:});
  k2 = feval(file_name , t + h/2 , y + h*k1/2 , varargin{:});
  k3 = feval(file_name , t + h/2 , y + h*k2/2 , varargin{:});
  k4 = feval(file_name , t + h   , y + h*k3   , varargin{:});
  
  y = y+h*(k1+2*k2+2*k3+k4)/6;
  
  t = t+h;
  
  % Assign the output
  Y(i,:) = y';
  T(i)   = t;
end

%%% EOF