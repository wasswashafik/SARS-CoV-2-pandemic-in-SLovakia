% SIR homogeneus infection dynamics without vital dynamics
% Fit from data, simulation and daily statistics
% Gergely Takacs, March 2020
% No guarantees given whatsoever.
% See covid19.gergelytakacs.com for more
%https://www.cia.gov/library/publications/the-world-factbook/geos/lo.html
popSize=5.440602;            % Population size in millions
nPop=popSize*1E6;        % Population size

importData

fitBegin=12;

I=cumsum(Confirmed);      % Cumulative sum of daily cases, transpose to make it compatible w/ E. Cheynet's code
R=cumsum(Recovered);      % Cumulative sum of daily cases, transpose to make it compatible w/ E. Cheynet's code
D=cumsum(Deaths);         % Cumulative sum of daily cases, transpose to make it compatible w/ E. Cheynet's code
R=R+D;
I=I-R;                                 % Active infections
S=nPop-I;
S=S(fitBegin:end)
I=I(fitBegin:end)
R=R(fitBegin:end)

Ts = 1;                                 % Sampling [1 day]
data = iddata([S I R],[],Ts);           % Create identification data object
data.TimeUnit='days';


data.OutputName = [{'Susceptible'};{'Infected'};{'Removed'}];              % Output name
data.OutputUnit = [{'Cases'};{'Cases'};{'Cases'}];                          % Output unit



%% Initial guess of model parameters
beta  = 1/10;             % [cases] Average base infection factor
gamma = 1/20;              % [days]  Removal rate
gammaTau=100;

%% Model structure
FileName      = 'SIR_ODE_Test2';              % File describing the SIR model structure
Order         = [3 0 3];                % Model orders [ny nu nx]
Parameters    = [beta,gamma,gammaTau];                % Initial values of parameters
InitialStates = [S(1);I(1);R(1)];           % Initial values of  [S I R] states
Ts            = 0;                      % Time-continuous system

% Set identification options
SIRinit = idnlgrey(FileName,Order,Parameters,InitialStates,Ts,'TimeUnit','days','Name','SIR Model');
SIRinit = setpar(SIRinit,'Name',{'beta','gamma','gammaTau'});
SIRinit = setinit(SIRinit,'Name',{'Susceptible' 'Infected' 'Removed'});

% beta
SIRinit.Parameters(1).Minimum = 1/31;
SIRinit.Parameters(1).Maximum = 1/1;
SIRinit.Parameters(1).Fixed = false;
    
% gamma
SIRinit.Parameters(2).Minimum = 1/100
SIRinit.Parameters(2).Maximum = 1/1; % Mean deaths 17 days, mean recoveries
SIRinit.Parameters(2).Fixed = false; 

% gammaTau
SIRinit.Parameters(3).Minimum = 1/100
SIRinit.Parameters(3).Maximum = 300; % Mean deaths 17 days, mean recoveries
SIRinit.Parameters(3).Fixed = false; 

% Susceptibles
SIRinit.InitialStates(1).Fixed = false;

SIRinit.InitialStates(2).Fixed = false;   % Let this parameter free, overall results will be better. True number unknown anyways
SIRinit.InitialStates(2).Minimum = -1000;     % Cannot be negative
SIRinit.InitialStates(2).Maximum = 1000;   % Unlikely to be more

SIRinit.InitialStates(3).Fixed = false;   % Yet again, we can let this parameter free.
SIRinit.InitialStates(3).Minimum = -1000;
SIRinit.InitialStates(3).Maximum = 10000;

% Identify model
opt = nlgreyestOptions('Display','on'); %gna
%opt.Regularization.Lambda = 0.8;
opt.OutputWeight=diag([0.5,1,0.5]);
opt.EstimateCovariance=1;
opt.SearchMethod='auto';                      %'auto' (default) | 'gn' | 'gna' | 'lm' | 'grad' | 'lsqnonlin' | 'fmincon'
opt.SearchOption.MaxIter= 100;                % Maximal number of iterations
SIR = nlgreyest(data,SIRinit,opt);              % Run identification procedure


    
%% Just internal comparison, not to output
%compare(data,SIR);                           % Compare data to model
SIR                                          % List model parameters
grid on                                        % Turn on grid

%% Simulate for generic data
udata = iddata([],zeros(365,0),1);
opt = simOptions('InitialCondition',[]);
SIRsim = sim(SIR,udata,opt);
Isim=round(SIRsim.OutputData(:,2));
%figure(1)
%SIR.Parameters(1).Value=1/3;
%SIR.Parameters(2).Value=1/10;
%sim(SIR,udata,opt)



figure(2)
compare(SIR,data)
disp(['1/beta: ',num2str(1/SIR.Parameters(1).Value),' [days]'])
disp(['1/gamma: ',num2str(1/SIR.Parameters(2).Value),' [days]'])
disp(['1/gammaTau: ',num2str(SIR.Parameters(3).Value),' [days]'])
disp(['R0: ',num2str((SIR.Parameters(1).Value)/(SIR.Parameters(2).Value)),' [cases]'])


return

%% Growth factor
growthFactorSIR= ([Isim; 0]./[0; Isim]);
growthFactorSIR= (growthFactorSIR(2:end-1)-1)*100;
gFSIR=mean(growthFactorSIR);



%% Report
R0est=SIR.Parameters(1).Value;
dRest=SIR.Parameters(2).Value;
N0est=SIR.InitialStates(2).Value;
%model.Report.Fit
MSE=SIR.Report.Fit.MSE;

%% simulate to find zero day data
opt2 = simOptions('InitialCondition',[nPop-1; 1; 0]);
SIRsim2 = sim(SIR,udata,opt2);
Isim2=round(SIRsim2.OutputData(:,2));
d0est= max(find(Isim2<N0est));