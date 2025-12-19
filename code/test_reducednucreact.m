% Right hand side
PD = @PD_reducednucreact;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time interval
 tstart = 0;
 tend = 1500;


 % (Initial) Time step size
 dt0 = 1.0e-2; %FE
 %dt0 = 0.02070; %RK2
 %dt0 = 0.0287; % RK4
 %dt0 = 3;


 % Initial values

lambda = 0.07741;
beta = 0.0065;

LAMBDA = 1e-4;
kc = 0.05;
alphatilde = 5*1e-5;

Phi0 = 10;
c0 = Phi0/LAMBDA*(beta/lambda);
rho0 = 0.5*beta;
c0 = (beta-rho0)*Phi0/(LAMBDA*lambda);
T0 = 300;

y0 = [Phi0;c0;rho0]; 

% Visualization amplification factors
%visfac = [1 5 5 1 5 1 1e2 1e2 1];


% Compound composition matrix
 E=[];
% dtfunc = @(dt) 1.1*dt;
 
