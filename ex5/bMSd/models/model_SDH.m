function SP = model_SDH()
% 
% model of 
% the SCHUNK SDH hand
%
% the mass and inertia parameters are arbitrary set
% the center of mass is assumed to co-inside with the input joint
%

% definition of system structure
% ----------------------------------------
SP.C = [0 1 2 3 1 5 6 1 8];
SP.n = length(SP.C)-1;
SP.mode = 1;

% definition of joints
% ----------------------------------------
SP.J(1).t    = [ -19.05     33.0      121 ]';
SP.J(1).f    = [  0.0     0.0      0.0 ]';
SP.J(1).rpy  = [  pi     0.0      0.0 ]';
SP.J(1).type = 'R';

SP.J(2).t    = [  0.0    0.0    -17.5 ]';
SP.J(2).f    = [  0.0     0.0      0.0 ]';
SP.J(2).rpy  = [  pi/2    0.0      0.0 ]';
SP.J(2).type = 'R';

SP.J(3).t    = [  0.0     -86.5      0.0 ]';
SP.J(3).f    = [  0.0     0.0      0.0 ]';
SP.J(3).rpy  = [  0.0     0.0      0.0 ]';
SP.J(3).type = 'R';

SP.J(4).t    = [ -19.05    -33.0      121   ]';
SP.J(4).f    = [  0.0     0.0      0.0 ]';
SP.J(4).rpy  = [  0.0     0      0.0 ]';
SP.J(4).type = 'R';

SP.J(5).t    = [  0.0    0.0    17.5 ]';
SP.J(5).f    = [  0.0     0.0      0.0 ]';
SP.J(5).rpy  = [  pi/2    pi      0.0 ]';
SP.J(5).type = 'R';

SP.J(6).t    = [  0.0     86.5      0.0 ]';
SP.J(6).f    = [  0.0     0.0      0.0 ]';
SP.J(6).rpy  = [  0.0     0.0      0.0 ]';
SP.J(6).type = 'R';

SP.J(7).t    = [ 38.105     0.0      138.5 ]';
SP.J(7).f    = [  0.0     0.0      0.0 ]';
SP.J(7).rpy  = [ pi/2     0      0.0 ]';
SP.J(7).type = 'R';

SP.J(8).t    = [  0.0     86.5      0.0 ]';
SP.J(8).f    = [  0.0     0.0      0.0 ]';
SP.J(8).rpy  = [  0.0     0.0      0.0 ]';
SP.J(8).type = 'R';

% SP.J(8).t    = [  0.0     0.0      0.0 ]';
% SP.J(8).f    = [  0.0     0.0      0.0 ]';
% SP.J(8).rpy  = [  0.0    pi/2      0.0 ]';
% SP.J(8).type = 'R';

% definition of links
% ----------------------------------------
for iL = 1:SP.n+1
  SP.L(iL).m = 1;
  SP.L(iL).I = eye(3);
end

% definition of end-effectors
% ----------------------------------------
%first distal, then proximal finger pads

SP.bN = [4 7 9 3 6 8];
SP.bP = [  0.0  0.0  0.0   0.0   0.0   0.0;
         -68.5 68.5 68.5 -43.25 43.25 43.25;
           0.0  0.0  0.0   0.0   0.0   0.0];

%%%EOF