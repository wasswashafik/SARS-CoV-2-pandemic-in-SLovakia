% SIR homogeneus infection dynamics with vital dynamics
% Fit from data, simulation and daily statistics
% Gergely Takacs, March 2020
% No guarantees given whatsoever.
% See covid19.gergelytakacs.com for more
% Stay at home, wash your hands.

warning('off','Ident:general:modelDataTU'); % Stop 'sim' whining about time units, they are just fine


%% 2019 demographic data for Slovakia

% https://www.cia.gov/library/publications/the-world-factbook/geos/lo.html
demBirth= 57054; % [people] Births in 2019 (total)
demDeath= 53234; % [people] Deaths in 2019 (total)
dayYear=365;     % [days]   2019 was not a leap year
demBirthDay=demBirth/365; % [people] Average birth per day in 2019
demDeathDay=demDeath/365; % [people] Average death per day in 2019
lambda = demBirthDay/nPop; % This is adjusted for the population
mu     = demDeathDay/nPop; % This is adjusted for the population
% CIA 2020 est is 9.3 births / 1000 pop and 10.1 / 1000 pop.
% Using last year data though.
%% Prepping data

R=cumsum(Recovered)+cumsum(Deaths);  % Number of total removals, cumulative sum of new recoveries and deaths
I=Nd-R;                              % Number of actively infected is total cases - removed
S=nPop*(1+(lambda-mu)*Day)-I-R;      % Susceptible population left, adjusted for vital dynamics

SIR_fitBegin      = fitbegin;

Ts = 1;                              % Sampling [1 day]
data = iddata([S I R],[],Ts);        % Create identification data object
data.TimeUnit='days';                % Time units
data.OutputName = [{'Susceptible'};{'Infected'};{'Removed'}];              % Output name
data.OutputUnit = [{'Cases'};{'Cases'};{'Cases'}];                          % Output unit
dataToFit=data(SIR_fitBegin:end);    % Create dataset itslef.

%% Initial guess of model parameters

beta=1/10; % 1/day infection rate
gamma=1/15.5; % 1/day removal rate


%% Model structure

FileName             = 'SIR_VD_ODE';                % File describing the SIR model structure
Order                = [3 0 3];                     % Model orders [ny nu nx]
Parameters           = [beta,gamma,lambda,mu];      % Initial values of parameters
InitialStates        = [S(SIR_fitBegin);I(SIR_fitBegin);R(SIR_fitBegin )];  % Initial values of  [S I R] states
Ts                   = 0;                           % Time-continuous system

% Set identification options
SIR_VDinit = idnlgrey(FileName,Order,Parameters,InitialStates,Ts,'TimeUnit','days','Name','SIR Model');
SIR_VDinit = setpar(SIR_VDinit,'Name',{'beta (infection rate)','gamma (removal rate)','lambda (birth rate)','mu (death rate)'});
SIR_VDinit = setinit(SIR_VDinit,'Name',{'Susceptible' 'Infected' 'Removed'});

% Contact rate beta (1/beta in days)
SIR_VDinit.Parameters(1).Minimum = 1/10;
SIR_VDinit.Parameters(1).Maximum = 1/3;
SIR_VDinit.Parameters(1).Fixed = false;
    
% Removal rate gamma (1/gamma in days)
SIR_VDinit.Parameters(2).Minimum = 1/40;
SIR_VDinit.Parameters(2).Maximum = 1/10; % Mean deaths 17 days, mean recoveries
SIR_VDinit.Parameters(2).Fixed = true; 

% Birth rate lambda (1/lambda)           % Birth rate will not be changed
                                         % by the epidemic. (Yet. :)
SIR_VDinit.Parameters(3).Minimum = lambda-0.2*lambda; % -20%
SIR_VDinit.Parameters(3).Maximum = lambda+0.2*lambda; % +20%
SIR_VDinit.Parameters(3).Fixed = false; 


% Death rate mu (1/mu)                   % Consider changing this, maybe
SIR_VDinit.Parameters(4).Minimum = mu-0.2*mu; % -20%
SIR_VDinit.Parameters(4).Maximum = mu+0.2*mu; % +20%
SIR_VDinit.Parameters(4).Fixed = false; 

% --------------Initial conditions--------------------
% Susceptibles
SIR_VDinit.InitialStates(1).Fixed = false;

SIR_VDinit.InitialStates(2).Fixed = false;   % Let this parameter free, overall results will be better. True number unknown anyways
SIR_VDinit.InitialStates(2).Minimum = 0;     % Cannot be negative
SIR_VDinit.InitialStates(2).Maximum = 100;   % Unlikely to be more

SIR_VDinit.InitialStates(3).Fixed = false;   % Yet again, we can let this parameter free.
SIR_VDinit.InitialStates(3).Minimum = -1000;
SIR_VDinit.InitialStates(3).Maximum = 1000;

%% Simulate initial guess
% udata = iddata([],zeros(365,0),1);
% opt = simOptions('InitialCondition',[nPop; 0; 0]);
% sim(SIR_VDinit,udata,opt);
% RES=sim(SIR_VDinit,udata,opt);
% RES.OutputData(end,1)-RES.OutputData(1,1)
%return

% Identify model
opt = nlgreyestOptions('Display','off','EstCovar',true,'SearchMethod','Auto'); %gna
opt.SearchOption.MaxIter = 50;                % Maximal number of iterations
SIR_VD = nlgreyest(dataToFit,SIR_VDinit,opt);              % Run identification procedure


round(SIR_VD.Report.Fit.FitPercent(2)*10)/10
%round(SIR_VD.Report.Fit.FitPercent*10)/10



%% Just internal comparison, not to output
%compare(dataToFit,SIR_VD);                           % Compare data to model
%SIR_VD                                          % List model parameters
%grid on                                        % Turn on grid


%% Simulate for fit
udataFit = iddata([],zeros(max(Day)-SIR_fitBegin+1,0),1);
optFit = simOptions('InitialCondition',[SIR_VD.InitialStates(1).Value;SIR_VD.InitialStates(2).Value;SIR_VD.InitialStates(3).Value]);
SIRsimFit = sim(SIR_VD,udataFit,optFit);
IsimFit=round(SIRsimFit.OutputData(:,2));

%% Simulate for predictions
udataPred = iddata([],zeros((max(DayPred)-max(Day)),0),1);
optPred = simOptions('InitialCondition',[S(end);I(end);R(end)]);
SIRsimPred = sim(SIR_VD,udataPred,optPred);
IsimPred=round(SIRsimPred.OutputData(:,2));

%NextDay
NdSIRnext=IsimPred(2); %+R if we 



%% Growth factor
growthFactorSIR= ([IsimFit; 0]./[0; IsimFit]);
growthFactorSIR= (growthFactorSIR(2:end-1)-1)*100;
gFSIR=mean(growthFactorSIR);



%% Report
R0est=(SIR_VD.Parameters(1).Value*SIR_VD.Parameters(3).Value)/(SIR_VD.Parameters(4).Value*(SIR_VD.Parameters(4).Value+(SIR_VD.Parameters(2).Value)));
betaInvEst=1/SIR_VD.Parameters(1).Value;
gammaInvEst=1/SIR_VD.Parameters(2).Value;
N0Fitest=round(SIR_VD.InitialStates(2).Value);
%model.Report.Fit
MSE=round(SIR_VD.Report.Fit.MSE);

%% simulate to find zero cases (day)
udataZD = iddata([],zeros(max(DayPred)+10,0),1);
optZD = simOptions('InitialCondition',[nPop-1; 1; 0]);
SIRsimZD = sim(SIR_VD,udataZD,optZD);
IsimZD=round(SIRsimZD.OutputData(:,2));
d0est=abs(SIR_fitBegin-max(find(IsimZD<=Nd(SIR_fitBegin)))); % Find the same number of cases here

%% Simulate for true cases
udataSymptoms = iddata([],zeros(max(DayPred),0),1);
optSymptoms = simOptions('InitialCondition',[SIR_VD.InitialStates(1).Value;SIR_VD.InitialStates(2).Value;SIR_VD.InitialStates(3).Value]);
SIRsimSymptoms = sim(SIR_VD,udataSymptoms,optSymptoms);
IsimSymptoms=round(SIRsimSymptoms.OutputData(:,2));




%d0est= max(find(Isim2<N0est));