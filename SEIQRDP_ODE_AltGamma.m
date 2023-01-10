% SEIQRDP homogeneus infection dynamics without vital dynamics
% This is just the state-space model of the dynamics

% Gergely Takacs, April 2020
% No guarantees given whatsoever. Stay at home, wash your hands, it will
% pass sooner or later.
%
% See covid19.gergelytakacs.com for more

% x(1) - S(t) - Susceptible population
% x(2) - E(t) - Exposed cases
% x(3) - I(t) - Infectious cases
% x(4) - Q(t) - Quarantined cases
% x(5) - R(t) - Recovered cases
% x(6) - D(t) - Deaths
% x(7) - P(t) - Insusceptible cases

% alpha       - protection rate
% beta        - infection rate
% sigma       - inverse of mean latent time
% delta       - inverse of mean quarantine time
% gamma       - cure rate
% mu          - mortality rate
% N           - population



function [dx y] = SEIQRDP_ODE_Mobility(t, x, u, alpha,beta,sigma,delta,gamma0,gamma1,gammaTau, mu, varargin)


N = x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7);                                                             % Computing total population
gamma=gamma0/(1+exp(-gamma1*(t-gammaTau)));


dx(1) =   -alpha*x(1) -beta*x(1)*x(3)/N;                                                            % (S(t)) Differential equation for Susceptible cases. People are either exposed to the virus or protected by social distancing measures, gov't. action, PPE, etc.
dx(2) =               +beta*x(1)*x(3)/N  -sigma*x(2);                                               % (E(t)) Differential equation for Exposed cases. People get exposed to the disease. This group then transfers to the actively infectious group.
dx(3) =                                   sigma*x(2)  -delta*x(3);                                  % (I(t)) Differential equation for Infected cases. Exposed people will eventually start to infect and become actively infectious.
dx(4) =                                                delta*x(3)   -gamma*x(4)  -mu*x(4);          % (Q(t)) Differential equation for Quarantined cases. These are the infections that are caught by testing, then quarantined in hospitals, facilities or at home.
dx(5) =                                                              gamma*x(4) ;                   % (R(t)) Differential equation for Recovered cases. People recover from the infection.
dx(6) =                                                                           mu*x(4);          % (D(t)) Differential equation for Deaths. Certain infections are fatal.
dx(7) =   alpha*x(1);                                                                               % (P(t)) Differential equation for Protected cases. People protected by social distancing measures, gov't. action, PPE, etc.

y=[x(4); x(5); x(6)];                                                                               % (y(t))   Output, where I(t), R(t), D(t) is measured

end