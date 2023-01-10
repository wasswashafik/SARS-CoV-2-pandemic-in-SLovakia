% SIR homogeneus infection dynamics without vital dynamics
% This is just the state-space model of the dynamics
% Gergely Takacs, March-April 2020
% No guarantees given whatsoever.
% Wash your hands, cover your face, stay at home.
% See covid19.gergelytakacs.com for more

function [dx y] = SIR_ODE(t, x, u, beta , gamma, gammaTau, varargin)



N = x(1) + x(2) + x(3);                  % Computing total population from original datay=[x(1); x(2); x(3)];                     % All is measured
y=[x(1) x(2) x(3)];

dx(1) = - beta*x(2)*x(1)/(1+gammaTau*(x(2)))*N;              % (S) Differential equation for Susceptible cases
dx(2) =   beta*x(2)*x(1)/(1+gammaTau*(x(2)))*N - gamma*x(2); % (I) Differential equation for Infected cases
dx(3) =   gamma*x(2);                    % (R) Differential equation for Removed cases

end