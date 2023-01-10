% SIR homogeneus infection dynamics with vital dynamics
% This is just the state-space model of the dynamics
% Gergely Takacs, March 2020
% No guarantees given whatsoever.
% See covid19.gergelytakacs.com for more

% x(1) - S(t) - Susceptible population
% x(2) - I(t) - Active infections
% x(3) - R(t) - Removed cases (recoveries and deaths)

function [dx y] = SIR_VD_ODE(t, x, u, beta, gamma, lambda, mu, varargin)

N = x(1)+x(2)+x(3);                  % Computing total population from original data
y=[x(1); x(2); x(3)];                     % All is measured
 
dx(1) =   lambda*N - beta*x(2)*x(1)/N              - mu*x(1); % (S) Differential equation for Susceptible cases. People are born, die and get infected.
dx(2) =              beta*x(2)*x(1)/N - gamma*x(2) - mu*x(2);   % (I) Differential equation for Infected cases. People recover but also die according to normal demographics.
dx(3) =                                 gamma*x(2) - mu*x(3);   % (R) Differential equation for Removed cases. People recover, die from infection but also from natural cases.
end