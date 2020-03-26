function [dydt, algvars] = SIRode(t,y,params)
% ODEfunc defines the system of ODEs describing the model
% 3/24/2020 JS: initial implementation

% Assign names for parameter values and state variables
[N,~,kappa,tau,gamma] = params{:};
beta = kappa*tau/N; % [1/days/person]
S = y(1); % Susceptible [# of people]
I = y(2); % Infected [# of people]
R = y(3); % Recovered [# of people]

% Differential equations;
dS = -beta*I*S;           % [people/day]             
dI = beta*I*S - gamma*I;  % [people/day]
dR = gamma*I;               % [people/day]

dydt = [dS;dI;dR];  % Reassemble differential equations
algvars = [];       % optional for seeing algebraic variables