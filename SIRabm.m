function SIRabm
% runs the SIR agent-based model
% 3/25/2020 Jeff Saucerman: initial implementation

%% Define parameters
N = 100;    % [number of people]
initialFractionInfected = 0.05;  
kappa = 5;  % [contacts/day/person]
tau = 0.5;    % [] transmissibility fraction            
gamma = 1/5;  % [1/days] rate of recovery, gamma=1/timeRecovery
beta = kappa*tau/N; % [1/days/person]

tspan = [0 25];
dt = 0.01;
%% Run single simulation

% initialize agents
agents = {};
for agentNum = 1:N
    if rand() < initialFractionInfected
        agents(agentNum) = {'I'};
    else
        agents(agentNum) = {'S'};
    end
end
numInfected = sum(strcmp(agents,'I'));

% run simulation
agentsTimecourse(1,:) = agents;
t = [tspan(1):dt:tspan(2)];
for tstep = 2:numel(t)
%     agentOrder = randperm(numel(agents)); % randomize order of agent updates   
    numInfected(tstep) = sum(strcmp(agents,'I')); % sums the total number of infected 
    probS2I(tstep) = beta*numInfected(tstep)*dt; % DEBUG: CHECK COMPARISON WITH ODE AND SPATIAL ABM SIMULATIONS 
    probI2R = gamma*dt;
    for agentNum = 1:numel(agents)
        if strcmp(agents{agentNum},'S')
            if rand()<probS2I(tstep)
                agents{agentNum} = 'I';
            end
        elseif strcmp(agents{agentNum},'I')
            if rand()<probI2R
                agents{agentNum} = 'R';
            end
        end
    end
    agentsTimecourse(tstep,:) = agents;
end

% plot timecourse
Sarray = strcmp(agentsTimecourse,'S');
Iarray = strcmp(agentsTimecourse,'I');
Rarray = strcmp(agentsTimecourse,'R');
S = sum(Sarray,2);
I = sum(Iarray,2);
R = sum(Rarray,2);

figure(1); clf
plot(t,S,'y',t,I,'r',t,R,'b','LineWidth',2);
xlabel('Time (days)'); ylabel('Number of people'); 
legend('S','I','R');

figure(2); clf
statesTimecourse = Sarray+2*Iarray+3*Rarray; % S: 1, I: 2, R: 3
imagesc(statesTimecourse','XData',t,'YData',N:1);
colormap([1,1,0;1,0,0;0,0,1])
lcolorbar({'S','I','R'});
xlabel('Time (days)'); ylabel('Individuals');