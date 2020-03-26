function SIRsensivity
% various sensitivity analyses using the SIR model
% 03/24/2020 Jeff Saucerman: initial implementation
% adapted from toyModelSensitivity

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
Ydefault = max(I);
figure(1);
subplot(1,2,1);
plot(t,S,'y',t,I,'r',t,R,'b','LineWidth',2);
xlabel('Time (days)'); ylabel('Number of people'); 
legend('S','I','R');

%% Sensitivity analysis coefficients
deltaParams = 1.1; % 10% increase
for paramNum = 1:length(params)    
    paramsNew = params;
    paramsNew{paramNum} = deltaParams*params{paramNum};
    N = paramsNew{1};
    initialFractionInfected = paramsNew{2};
    I0 = initialFractionInfected*N;
    S0 = N-I0;
    y0 = [S0;I0;0];
    try
        [t,y] = ode23(@SIRode,tspan,y0,options,paramsNew);
        S = y(:,1);
        I = y(:,2);
        R = y(:,3);
        Y(paramNum) = max(I);
        Sens(paramNum) = (Y(paramNum)-Ydefault)./(paramsNew{paramNum}-params{paramNum})*params{paramNum}/Ydefault;
    catch       % in case the tested parameters generate numerical problems
        Y(paramNum) = NaN;
        Sens(paramNum) = NaN;
        disp(['Error at paramNum = ',num2str(paramNum),', delta = ',num2str(deltaParams)]);
    end
end
% hold off;

% Plotting sensitivity coefficients
subplot(1,2,2); 
bar(Sens); 
paramNames = {'N','init%I','\kappa','\tau','\gamma'};
set(gca,'xticklabel',paramNames)
xlabel('Parameter'); ylabel('Sens=dY/dp*p/Y'); title('Sensitivity coefficients');

%% Response curves
deltaParams = 10.^[-1:.1:1];    % more refined range for response curves
for deltaNum=1:length(deltaParams)
    for paramNum = 1:length(params)    
        paramsNew = params;
        paramsNew{paramNum} = deltaParams(deltaNum)*params{paramNum};
        N = paramsNew{1};
        initialFractionInfected = paramsNew{2};
        I0 = initialFractionInfected*N;
        S0 = N-I0;
        y0 = [S0;I0;0];
        try
            [t,y] = ode23(@SIRode,tspan,y0,options,paramsNew);
            S = y(:,1);
            I = y(:,2);
            R = y(:,3);
            Y(paramNum,deltaNum) = max(I);
            Sens(paramNum,deltaNum) = (Y(paramNum,deltaNum)-Ydefault)./(paramsNew{paramNum}-params{paramNum})*params{paramNum}/Ydefault;            
        catch       % in case the tested parameters generate numerical problems
            Y(paramNum,deltaNum) = NaN;
            Sens(paramNum,deltaNum) = NaN;
            disp(['Error at paramNum = ',num2str(paramNum),', delta = ',num2str(deltaParams(deltaNum))]);
        end
    end
end

figure(2); clf
size = ceil(sqrt(length(params)));
for i=1:length(params)
    subplot(size,size,i);
    semilogx(deltaParams,Y(i,:),'o-','LineWidth',2);
    axis([min(deltaParams) max(deltaParams) 0 N]);
    xlabel(paramNames{i}); ylabel('Peak Infections');
end

%% 2D sensitivity matrix
I0 = initialFractionInfected*N;
S0 = N-I0;
y0 = [S0;I0;0];
tspan = [0 5];
options = [];
[t,y] = ode23(@SIRode,tspan,y0,options,params); % run default simulation
y_fin_orig = y(end,:)';
for i=1:length(params)
    paramsNew = params; % Initialize new parameters to old ones
    paramsNew{i} = params{i}*1.1; % Increase parameter by 10%
    N = paramsNew{1};
    initialFractionInfected = paramsNew{2};
    I0 = initialFractionInfected*N;
    S0 = N-I0;
    y_initial = [S0;I0;0];
    [~,y_new] = ode23(@SIRode,tspan,y_initial,options,paramsNew);
    y_final(:,i) = y_new(end,:)'; % Make column of y_final equal to steady state values of y_new
    
    % Generate normalized sensitivity coefficients
    S2d(:,i) = (y_final(:,i) - y_fin_orig(:,1))./(paramsNew{i}-params{i}).*params{i}./y_fin_orig(:,1);
end

figure(3);
imagesc(S2d)
colormap(redbluecmap);
colorbar; caxis([-2 2])
xlabel('Perturbed Parameter'); ylabel(['Sensitivity of State Variable, t = ',num2str(tspan(2))]);
yNames = {'S','I','R'};
paramNames = {'N','init%I','\kappa','\tau','\gamma'};
set(gca,'YTick',1:length(yNames));
set(gca,'YTickLabel',yNames);
set(gca,'XTick',1:length(paramNames));
set(gca,'XTickLabel',paramNames);

%% Uncertainty analysis
paramNames = {'N','init%I','\kappa','\tau','\gamma'};
cv = 0.3*[1,1,1,1,1]; % 30% coefficient of variation for all params
% cv = 0.3*[1,0,0,0,0]; % 30% CV in N
% cv = 0.3*[0,0,1,0,0]; % 30% CV in kappa
% cv = 0.3*[0,0,0,0,1]; % 30% CV in gamma
tspan = [0 25];
tinterp = [tspan(1):1:tspan(2)];
options = [];
figure(4); 
subplot(2,1,1); title('Uncertainty analysis'); hold on
for run=1:10
    paramsNew = num2cell([params{:}] + cv.*[params{:}].*randn(1,length(params)));               % randomize new params with cv
    N = paramsNew{1};
    initialFractionInfected = paramsNew{2};
    I0 = initialFractionInfected*N;
    S0 = N-I0;
    y0 = [S0;I0;0];
    try
        [t,y] = ode23(@SIRode,tspan,y0,options,paramsNew);
        I = y(:,2);                             % entire timecourse of I
        Iinterp(:,run) = interp1(t,I,tinterp);  % matrix of all I timecourses, interpolated of I
        plot(tinterp,Iinterp(:,run)); 
    catch       % in case the tested parameters generate numerical problems
        disp(['Error with params = ']);%,num2str(params)]);
    end
end

xlabel('time'); ylabel('I'); hold off;

% Plot mean and 95% confidence intervals (assuming normal distribution)
subplot(2,1,2);
Imean = mean(Iinterp,2);
Istd = std(Iinterp,0,2);
plot(tinterp,Imean,'k-',tinterp,Imean+1.96.*Istd,'k--',tinterp,Imean-1.96.*Istd,'k--'); 
xlabel('time'); ylabel('I');