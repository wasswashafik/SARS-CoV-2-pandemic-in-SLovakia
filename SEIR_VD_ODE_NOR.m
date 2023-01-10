% SEIR homogeneus infection dynamics with vital dynamics
% This is just the state-space model of the dynamics
% Gergely Takacs, March 2020
% No guarantees given whatsoever. Stay at home, wash your hands, it will
% pass sooner or later.
%
% See covid19.gergelytakacs.com for more

% x(1) - S(t) - Susceptible population
% x(2) - E(t) - Exposed cases
% x(3) - I(t) - Active infections
% x(4) - R(t) - Removed cases (recoveries and deaths)

function [dx y] = SEIR_VD_ODE(t, x, u, beta, sigma, gamma, lambda, mu, varargin)

N = x(1) + x(2) + x(3) + x(4);                  % Computing total population from original data

y=[x(1); x(3); x(4)];                     % All is measured

gammaT=0.005*(1-exp(-gamma*t));


dx(1) =   lambda*N - beta*x(3)*x(1)/N          - mu*x(1)/N;   % (S) Differential equation for Susceptible cases. People are born, die and get infected.
dx(2) =              beta*x(3)*x(1)/N -(mu/N+sigma)*x(2);     % (E) Differential equation for Exposed cases. People get exposed to the disease. This group then transfers to the infectious group but some die of natural cases.
dx(3) =             sigma*x(2)        -(gammaT+mu/N)*x(3);     % (I) Differential equation for Infected cases. People recover but also die according to normal demographics.
dx(4) =             gammaT*x(3)                 - mu*x(4)/N;   % (R) Differential equation for Removed cases. People recover, die from infection but also from natural cases.
end