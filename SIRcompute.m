function SIRcompute
% runs the SIR ODE model
% 3/24/2020 Jeff Saucerman: initial implementation

%% Define parameters
N = 100;    % [number of people]
initialFractionInfected = 0.05;  
kappa = 5;  % [contacts/day/person]
tau = 0.5;    % [] transmissibility fraction            
gamma = 1/5;  % [1/days] rate of recovery, gamma=1/timeRecovery
params = {N,initialFractionInfected,kappa,tau,gamma};

%% Run single simulation
I0 = initialFractionInfected*N;
S0 = N-I0;
beta = kappa*tau/N; % [1/days/person]
Re = S0*beta/gamma; % effective reproductive number

y0 = [S0;I0;0];
tspan = [0 25];
options = [];
[t,y] = ode23(@SIRode,tspan,y0,options,params);

% plot timecourse
S = y(:,1);
I = y(:,2);
R = y(:,3);
figure(1);
subplot(1,2,1);
plot(t,S,'y',t,I,'r',t,R,'b','LineWidth',2);
xlabel('Time (days)'); ylabel('Number of people'); 
legend('S','I','R');

%% Sensitivity analysis (single parameter, single output)
paramRange = [0:10];
for i=1:length(paramRange)
    kappa = paramRange(i);
    params = {N,initialFractionInfected,kappa,tau,gamma};
    I0 = initialFractionInfected*N;
    S0 = N-I0;
    y0 = [S0;I0;0]; 
    [t,y] = ode23(@SIRode,tspan,y0,options,params);
    S = y(:,1);
    I = y(:,2);
    R = y(:,3);
    Ipeak(i) = max(I);
end
subplot(1,2,2);
plot(paramRange,Ipeak,'ro-','LineWidth',2);
xlabel('\kappa (contacts/day/person)'); ylabel('Peak Infected');