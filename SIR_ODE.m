% SIR homogeneus infection dynamics without vital dynamics
% This is just the state-space model of the dynamics
% Gergely Takacs, March 2020
% No guarantees given whatsoever.
% See covid19.gergelytakacs.com for more

function [dx y] = SIR_ODE(t, x, u, R0, dR, varargin)

gamma = 1/dR;                            % Case removal rate 1/g = e.g. 14 [days]
beta=   R0*gamma;                        % Infection rate 1/bate = 5 [days]
N = x(1) + x(2) + x(3);                  % Computing total population from original data

y=[x(1); x(2); x(3)];                     % All is measured

dx(1) = - beta*x(2)*x(1)/N;              % (S) Differential equation for Susceptible cases
dx(2) =   beta*x(2)*x(1)/N - gamma*x(2); % (I) Differential equation for Infected cases
dx(3) =   gamma*x(2);                    % (R) Differential equation for Removed cases
end