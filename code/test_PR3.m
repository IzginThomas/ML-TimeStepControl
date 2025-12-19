%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)               %%%
%%% Date: 09/29/2022                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Nonlinear problem 4D                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Right hand side
PD = @PD_PR3;

% Time interval
tstart = 0;
tend = 20*pi;
%tend = 3e-1;
% Initial values
y0 = [2; 2; 2]; % =g(0), see PD_nonlin3


% Time step size
dt0 = 1;

% Visualization amplification factors
%visfac = [1 1 1 1];

% Plot range
yrange =  [-0.1 6.1];

% Compound composition matrix
%E = [1 0 0 1 1 1; 0 1 1 0 0 0; 1 1 1 1 1 1];
E = [];

%Exact solution CHECK WITH PD_nonlin3.m !!
oscifac =@(t) 0.1*cos(0.5.*t); 
k1 = 0.3;
k2 = 1-k1; 
exsol=@(t) [2 + k1.*sin(oscifac(t).*t) 2 + k2.*sin(oscifac(t).*t) 2 - sin(oscifac(t).*t)];