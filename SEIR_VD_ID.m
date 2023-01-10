% SEIR homogeneus infection dynamics with vital dynamics
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
%S=nPop*(1+(lambda-mu)*Day)-I-R;      % Susceptible population left, adjusted for vital dynamics
S=nPop-I-R;

SIR_fitBegin      = fitbegin;               % Manual override because so far it has no meaning to do this

Ts = 1;                              % Sampling [1 day]
data = iddata([S I R],[],Ts);        % Create identification data object
data.TimeUnit='days';                % Time units
data.OutputName = [{'Susceptible'};{'Infected'};{'Removed'}];              % Output name
data.OutputUnit = [{'Cases'};{'Cases'};{'Cases'}];                          % Output unit
dataToFit=data(SIR_fitBegin:end);    % Create dataset itslef.


%% Initial guess of model parameters
% Incubation period 2-14 days, 5.2 mean https://www.worldometers.info/coronavirus/coronavirus-incubation-period/
% is S-E-I the incubation period?
beta  =1/ 6;       % [1/day] exposure/contact rate
sigma =1/ 0.5;       % [1/day] infection rate (latent period)
gamma =1/ 15.5;    % [1/day] removal rate   (infectious period)
E0=3;            % [cases] Number of exposed at initial state


%% Model structure

FileName             = 'SEIR_VD_ODE';              % File describing the SIR model structure
Order                = [3 0 4];                    % Model orders [ny nu nx]
Parameters           = [beta,sigma,gamma,lambda,mu];     % Initial values of parameters
InitialStates        = [S(SIR_fitBegin); E0; I(SIR_fitBegin);R(SIR_fitBegin )];  % Initial values of  [S I R] states
Ts                   = 0;                      % Time-continuous system

% Set identification options
SEIR_VDinit = idnlgrey(FileName,Order,Parameters,InitialStates,Ts,'TimeUnit','days','Name','SIR Model');
SEIR_VDinit = setpar(SEIR_VDinit,'Name',{'beta (exposure rate)','sigma (infection rate)','gamma (removal rate)','lambda (birth rate)','mu (death rate)'});
SEIR_VDinit = setinit(SEIR_VDinit,'Name',{'Susceptible' 'Exposed' 'Infected' 'Removed'});


% --------------Parameters--------------------

% S->E Exposure (contact) rate beta (1/beta in days)
SEIR_VDinit.Parameters(1).Minimum = 1/30;      % [1/days] Cannot be reasonably more than 20 days      
SEIR_VDinit.Parameters(1).Maximum = 1/1;  % [1/days] Cannot be reasonably less than an hour days      
SEIR_VDinit.Parameters(1).Fixed = false;

% E->I Infection rate sigma (1/sigma in days)
SEIR_VDinit.Parameters(2).Minimum = 1/7;      % [1/days] Cannot be reasonably more than 20 days      
SEIR_VDinit.Parameters(2).Maximum = 1/(1/24);  % [1/days] Cannot be reasonably less than an hour days      
SEIR_VDinit.Parameters(2).Fixed = false;
    
% I->R Removal rate gamma (1/gamma in days)
%SEIR_VDinit.Parameters(3).Minimum = 1/31;
%SEIR_VDinit.Parameters(3).Maximum = 1/10; % Mean deaths 17 days, mean recoveries
SEIR_VDinit.Parameters(3).Fixed = false; 

% Birth rate lambda (1/lambda)           % Birth rate will not be changed
                                         % by the epidemic. (Yet. :)

SEIR_VDinit.Parameters(4).Value = 0;                                          
%SEIR_VDinit.Parameters(4).Minimum = lambda-0.25*lambda; % -25%
%SEIR_VDinit.Parameters(4).Maximum = lambda+0.25*lambda; % +25%
SEIR_VDinit.Parameters(4).Fixed = true; 


% Death rate mu (1/mu)                   % Consider changing this, maybe
SEIR_VDinit.Parameters(5).Value = 0;  
%SEIR_VDinit.Parameters(5).Minimum = mu-0.25*mu;   % after first death begin. Since the pandemic deaths are not counted for
%SEIR_VDinit.Parameters(5).Maximum = mu+0.25*mu;   % the normal death rate. Or as a slack variable?
SEIR_VDinit.Parameters(5).Fixed = true;

% --------------Initial conditions--------------------
% Susceptibles
SEIR_VDinit.InitialStates(1).Fixed = true;

SEIR_VDinit.InitialStates(2).Fixed = false;   % Let this parameter free, overall results will be better. True number unknown anyways
SEIR_VDinit.InitialStates(2).Minimum = 0;     % Cannot be negative
SEIR_VDinit.InitialStates(2).Maximum = 1000;   % Unlikely to be more

SEIR_VDinit.InitialStates(3).Fixed = true;   % Yet again, we can let this parameter free.


SEIR_VDinit.InitialStates(4).Fixed = true;   % Yet again, we can let this parameter free.


%% Simulate initial guess
% udata = iddata([],zeros(365,0),1);
% opt = simOptions('InitialCondition',[nPop-3; 2 ; 1; 0]);
% sim(SEIR_VDinit,udata,opt);
% RES=sim(SEIR_VDinit,udata,opt);
% RES.OutputData(end,1)-RES.OutputData(1,1);
% return

% Identify model
opt = nlgreyestOptions('Display','on','EstCovar',true); %gna
opt.SearchMethod='lsqnonlin'


opt.SearchOption.MaxIter = 50;                % Maximal number of iterations
SEIR_VD = nlgreyest(dataToFit,SEIR_VDinit,opt);              % Run identification procedure


%% Just internal comparison, not to output
disp(['SEIR model fit for I (x(3)) is: ',num2str(round(SEIR_VD.Report.Fit.FitPercent(2)*10)/10),'%'])
compare(dataToFit,SEIR_VD);                           % Compare data to model
%hold off
%grid on



%% Simulate for fit
udataFit = iddata([],zeros(max(Day)-SIR_fitBegin+1,0),1);
optFit = simOptions('InitialCondition',[SEIR_VD.InitialStates(1).Value;SEIR_VD.InitialStates(2).Value;SEIR_VD.InitialStates(3).Value;SEIR_VD.InitialStates(4).Value]);
[SIRsimFit Y_Fit X_Fit] = sim(SEIR_VD,udataFit,optFit);
IsimFit=round(SIRsimFit.OutputData(:,2));

%% Simulate for long time predictions
udataFit = iddata([],zeros(200,0),1);
optFit = simOptions('InitialCondition',[SEIR_VD.InitialStates(1).Value;SEIR_VD.InitialStates(2).Value;SEIR_VD.InitialStates(3).Value;SEIR_VD.InitialStates(4).Value]);
sim(SEIR_VD,udataFit,optFit);


%% Simulate for short time predictions
udataPred = iddata([],zeros((max(DayPred)-max(Day)),0),1);
optPred = simOptions('InitialCondition',[S(end); X_Fit(end,2);I(end);R(end)]);
SIRsimPred = sim(SEIR_VD,udataPred,optPred);
IsimPred=round(SIRsimPred.OutputData(:,2));

%NextDay
NdSIRnext=IsimPred(2); %+R if we 

%% Growth factor
growthFactorSIR= ([IsimFit; 0]./[0; IsimFit]);
growthFactorSIR= (growthFactorSIR(2:end-1)-1)*100;
gFSIR=mean(growthFactorSIR);


%% Report
%R0est=(SEIR_VD.Parameters(1).Value*SEIR_VD.Parameters(3).Value)/(SEIR_VD.Parameters(4).Value*(SEIR_VD.Parameters(4).Value+SEIR_VD.Parameters(2).Value));



betaEst=SEIR_VD.Parameters(1).Value;
sigmaEst=SEIR_VD.Parameters(2).Value;
gammaEst=SEIR_VD.Parameters(3).Value;
lambdaEst=SEIR_VD.Parameters(4).Value;
muEst=SEIR_VD.Parameters(5).Value;
N0Fitest=round(SEIR_VD.InitialStates(2).Value);
%model.Report.Fit
MSE=round(SEIR_VD.Report.Fit.MSE);
fitPerc=SEIR_VD.Report.Fit.FitPercent(2);

% Diagnostic

disp(['beta:  ',num2str(betaEst)])
disp(['sigma: ',num2str(sigmaEst)])
disp(['gamma: ',num2str(gammaEst)])




%% simulate to find zero cases (day)
udataZD = iddata([],zeros(max(DayPred)+10,0),1);
optZD = simOptions('InitialCondition',[nPop-1; 1; 0; 0]);
SIRsimZD = sim(SEIR_VD,udataZD,optZD);
IsimZD=round(SIRsimZD.OutputData(:,2));
d0est=abs(SIR_fitBegin-max(find(IsimZD<=Nd(SIR_fitBegin)))); % Find the same number of cases here

%% Simulate for true cases
udataSymptoms = iddata([],zeros(max(DayPred),0),1);
optSymptoms = simOptions('InitialCondition',[SEIR_VD.InitialStates(1).Value;SEIR_VD.InitialStates(2).Value;SEIR_VD.InitialStates(3).Value;SEIR_VD.InitialStates(4).Value]);
SIRsimSymptoms = sim(SEIR_VD,udataSymptoms,optSymptoms);
IsimSymptoms=round(SIRsimSymptoms.OutputData(:,2));

%d0est= max(find(Isim2<N0est));