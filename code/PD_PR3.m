%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)               %%%
%%% Date: 09/29/2022                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prothero & Robinson Problem with 3 Equations                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, D, f, S, r, Rrem, Rminus,varargout] = PD_PR3(t,y,varargin)

if nargin >2 && ~isempty(varargin{1})
    para = varargin{1};
else
    error('para is not specified in PD_nonlin3')
end
%para = 0.8; % between 0 and 1 CHECK WITH TEST_PR3
oscifac =@(t) 0.1*cos(0.5*t);
oscidiff =@(t) -0.05*sin(.5*t); % derivative of oscifac(t)
proddiff =@(t) oscidiff(t)*t + oscifac(t); % derivative of oscifac(t)*t
k1 = 0.3;
k2 = 1-k1; 

A = [-1 1-para para; para -1 1-para; 1-para para -1];


exsol=@(t) [2 + k1*sin(oscifac(t).*t);2 + k2*sin(oscifac(t).*t);2 - sin(oscifac(t).*t)]; %Exact solution CHECK WITH IC in test_nonlin3.m !!
g = exsol(t);
gprime = [k1*proddiff(t)*cos(oscifac(t)*t);k2*proddiff(t)*cos(oscifac(t)*t);-proddiff(t)*cos(oscifac(t)*t)];

P(1,2) = y(2) + para*y(3); 
P(1,3) = g(1) + para*g(2) + min(0,gprime(1));
P(2,1) = g(2) + para*g(3);
P(2,3) = y(3) + para*y(1) + min(0,gprime(2));   
P(3,1) = y(1) + para*y(2) + k1*min(0,gprime(3));
P(3,2) = g(3) + para*g(1) + k2*min(0,gprime(3));

D = P';


if nargout > 2
f = @(t,y) A*(y - [2 + k1*sin(oscifac(t)*t);2 + k2*sin(oscifac(t)*t);2 - sin(oscifac(t)*t)]) + [k1*proddiff(t)*cos(oscifac(t)*t);k2*proddiff(t)*cos(oscifac(t)*t);-proddiff(t)*cos(oscifac(t)*t)];
%f = @(t,y) [y(2) + para*y(3) + 2 + k1*sin(oscifac(t)*t) + para*(2 + k2*sin(oscifac(t)*t)) + min(0,k1*cos(oscifac(t)*t)) - (2 + k2*sin(oscifac(t)*t) + para*(2 - sin(oscifac(t)*t)) + y(1) + para*y(2) + min(0,k1*cos(oscifac(t)*t)));
   % 2 + k2*sin(oscifac(t)*t) + para*(2 - sin(oscifac(t)*t)) + y(3) + para*y(1) + min(0,k2*cos(oscifac(t)*t)) - (y(2) + para*y(3) + 2 - sin(oscifac(t)*t) + para*(2 + k1*sin(oscifac(t)*t)) + min(0,k2*cos(oscifac(t)*t)));
    %        y(1) + para*y(2) + min(0,k1*cos(oscifac(t)*t)) + 2 - sin(oscifac(t)*t) + para*(2 + k1*sin(oscifac(t)*t)) + min(0,k2*cos(oscifac(t)*t)) - (2 + k1*sin(oscifac(t)*t) + para*(2 + k2*sin(oscifac(t)*t)) + min(0,k1*cos(oscifac(t)*t)) + y(3) + para*y(1) + min(0,k2*cos(oscifac(t)*t)))];
end

if nargout > 3 %% Need to be done
 S=zeros(length(y));
  
  
   r = @(y) zeros(length(y),1);

end

if nargout > 5 

Rrem = zeros(length(y),1); % remainder
Rminus = Rrem;

end
if nargout > 7 

varargout{1} = exsol;

end
