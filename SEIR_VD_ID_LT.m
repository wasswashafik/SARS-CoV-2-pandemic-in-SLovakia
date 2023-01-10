% SEIR homogeneus infection dynamics with vital dynamics
% Fit from data, simulation and daily statistics
% This is an attempt to do long time predictions
% Gergely Takacs, April 2020
% No guarantees given whatsoever.
% See covid19.gergelytakacs.com for more
% Stay at home, wash your hands.

% LT means long term. This is used for long term predictions

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
S=nPop-I-R;

SIR_fitBeginLT      = fitbegin;               % Manual override because so far it has no meaning to do this

Ts = 1;                              % Sampling [1 day]
%dataLT = iddata([S I],[],Ts);        % Create identification data object
dataLT = iddata([S I R],[],Ts);        % Create identification data object
dataLT.TimeUnit='days';                % Time units
dataLT.OutputName = [{'Susceptible'};{'Infected'};{'Removed'}];
%dataLT.OutputName = [{'Susceptible'};{'Infected'}]; % Output name
dataLT.OutputUnit = [{'Cases'};{'Cases'};{'Cases'}];   
%dataLT.OutputUnit = [{'Cases'};{'Cases'}];   % Output unit
dataToFitLT=dataLT(SIR_fitBeginLT:end);    % Create dataset itslef.


%% Initial guess of model parameters
% Incubation period 2-14 days, 5.2 mean https://www.worldometers.info/coronavirus/coronavirus-incubation-period/
% is S-E-I the incubation period?
beta  =1/ 6;       % [1/day] exposure/contact rate
sigma =1/ 0.2;       % [1/day] infection rate (latent period)
gamma =1/ 20;    % [1/day] removal rate   (infectious period)
E0=3;            % [cases] Number of exposed at initial state


%% Model structure

FileName             = 'SEIR_VD_ODE_NOR';              % File describing the SIR model structure
Order                = [3 0 4];                    % Model orders [ny nu nx]
%Order                = [2 0 4];                    % Model orders [ny nu nx]
Parameters           = [beta,sigma,gamma,lambda,mu];     % Initial values of parameters
InitialStates        = [S(SIR_fitBegin); E0; I(SIR_fitBegin);R(SIR_fitBegin )];  % Initial values of  [S I R] states
Ts                   = 0;                      % Time-continuous system

% Set identification options
SEIR_VDinitLT = idnlgrey(FileName,Order,Parameters,InitialStates,Ts,'TimeUnit','days','Name','SIR Model');
SEIR_VDinitLT = setpar(SEIR_VDinitLT,'Name',{'beta (exposure rate)','sigma (infection rate)','gamma (removal rate)','lambda (birth rate)','mu (death rate)'});

% --------------Parameters--------------------

% S->E Exposure (contact) rate beta (1/beta in days)
InitialStates        = [S(SIR_fitBegin); E0; I(SIR_fitBegin);R(SIR_fitBegin )];  % Initial values of  [S I R] states

SEIR_VDinitLT.Parameters(1).Value = 1/6;
SEIR_VDinitLT.Parameters(1).Minimum = 1/30;      % [1/days] Cannot be reasonably more than 20 days      
SEIR_VDinitLT.Parameters(1).Maximum = 1/1;  % [1/days] Cannot be reasonably less than an hour days      
SEIR_VDinitLT.Parameters(1).Fixed = false;

% E->I Infection rate sigma (1/sigma in days)
SEIR_VDinitLT.Parameters(2).Value = 1/1;      % [1/days] Cannot be reasonably more than 20 days      
SEIR_VDinitLT.Parameters(2).Minimum = 1/5;      % [1/days] Cannot be reasonably more than 20 days      
SEIR_VDinitLT.Parameters(2).Maximum = 1/(1/24);  % [1/days] Cannot be reasonably less than an hour days      
SEIR_VDinitLT.Parameters(2).Fixed = false;

% Based on data the removal rate (unadjusted is) March 6-27, e.g. 20 days

% I->R Removal rate gamma (1/gamma in days)
SEIR_VDinitLT.Parameters(3).Value = 1/30; 
SEIR_VDinitLT.Parameters(3).Minimum = -inf;
SEIR_VDinitLT.Parameters(3).Maximum = inf; % Mean deaths 17 days, mean recoveries
SEIR_VDinitLT.Parameters(3).Fixed = false; 


% Turn off vital dynamics for estimation, since 
% our data in (S) does not reflect this anyways.
% Birth rate lambda (1/lambda)           
SEIR_VDinitLT.Parameters(4).Value = 0;                                          
SEIR_VDinitLT.Parameters(4).Fixed = true; 
% Death rate mu (1/mu)                   
SEIR_VDinitLT.Parameters(5).Value = 0;  
SEIR_VDinitLT.Parameters(5).Fixed = true;

% --------------Initial conditions--------------------
% Susceptibles
SEIR_VDinitLT.InitialStates(1).Name = 'Susceptible';
SEIR_VDinitLT.InitialStates(1).Fixed = true;

SEIR_VDinitLT.InitialStates(2).Name = 'Exposed';
SEIR_VDinitLT.InitialStates(2).Fixed = false;   % Let this parameter free, overall results will be better. True number unknown anyways
SEIR_VDinitLT.InitialStates(2).Minimum = 0;     % Cannot be negative
SEIR_VDinitLT.InitialStates(2).Maximum = 1000;   % Unlikely to be more

SEIR_VDinitLT.InitialStates(3).Name = 'Infected';
SEIR_VDinitLT.InitialStates(3).Fixed = false;   % Yet again, we can let this parameter free.

SEIR_VDinitLT.InitialStates(4).Name = 'Removed';
SEIR_VDinitLT.InitialStates(4).Fixed = false;   % Yet again, we can let this parameter free.


%% Simulate initial guess
% udata = iddata([],zeros(365,0),1);
% opt = simOptions('InitialCondition',[nPop-3; 2 ; 1; 0]);
% sim(SEIR_VDinit,udata,opt);
% RES=sim(SEIR_VDinit,udata,opt);
% RES.OutputData(end,1)-RES.OutputData(1,1);
% return

% Identify model
optLT = nlgreyestOptions('Display','on','EstCovar',true,'SearchMethod','Auto'); %gna
optLT.SearchOption.MaxIter = 50;                % Maximal number of iterations
SEIR_VD_LT = nlgreyest(dataToFitLT,SEIR_VDinitLT,optLT);              % Run identification procedure

%% Just internal comparison, not final public output
disp(['SEIR model fit for S (x(1)) is: ',num2str(round(SEIR_VD_LT.Report.Fit.FitPercent(1)*10)/10),'%'])
disp(['SEIR model fit for I (x(3)) is: ',num2str(round(SEIR_VD_LT.Report.Fit.FitPercent(2)*10)/10),'%'])
disp(['SEIR model fit for R (x(4)) is: ',num2str(round(SEIR_VD_LT.Report.Fit.FitPercent(3)*10)/10),'%'])
%compare(dataToFit,SEIR_VD_LT);                           % Compare data to model
hold off
grid on

% Diagnostic
%R0est=(SEIR_VD.Parameters(1).Value*SEIR_VD.Parameters(3).Value)/(SEIR_VD.Parameters(4).Value*(SEIR_VD.Parameters(4).Value+SEIR_VD.Parameters(2).Value));
%R0est=be
betaEst=SEIR_VD_LT.Parameters(1).Value;
sigmaEst=SEIR_VD_LT.Parameters(2).Value;
gammaEst=SEIR_VD_LT.Parameters(3).Value;
lambdaEst=SEIR_VD_LT.Parameters(4).Value;
muEst=SEIR_VD_LT.Parameters(5).Value;
N0Fitest=round(SEIR_VD_LT.InitialStates(2).Value);
%model.Report.Fit
MSE=round(SEIR_VD_LT.Report.Fit.MSE);
fitPerc=SEIR_VD_LT.Report.Fit.FitPercent(2);


disp(['beta:  ',num2str(betaEst),'[1/day], ',num2str(1/betaEst),'[day]'])
disp(['sigma: ',num2str(sigmaEst),'[1/day], ',num2str(1/sigmaEst*24),'[hrs]'])
disp(['gamma: ',num2str(gammaEst),'[1/day], ',num2str(1/gammaEst),'[day]'])
%disp(['R0: ',num2str(R0est),'[day]'])



%% Simulate for long time fit
udataFitLT = iddata([],zeros(200,0),1);
optFitLT = simOptions('InitialCondition',[SEIR_VD_LT.InitialStates(1).Value;SEIR_VD_LT.InitialStates(2).Value;SEIR_VD_LT.InitialStates(3).Value;SEIR_VD_LT.InitialStates(4).Value]);
%[SIRsimFit Y_Fit X_Fit] = sim(SEIR_VD,udataFit,optFit);
%IsimFit=round(SIRsimFit.OutputData(:,2));

%SEIR_VD_LT.Parameters(3).Value = 1/0.00001
sim(SEIR_VD_LT,udataFitLT,optFitLT);


return

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