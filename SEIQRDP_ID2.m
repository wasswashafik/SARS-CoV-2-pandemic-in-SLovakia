clc; clear;                   % Cleanup 
hold off;                   % Useful for tuning

% Use the wrong fits to comute a statistical uncertainity

warning('off','Ident:general:modelDataTU'); % Stop 'sim' whining about time units, they are just fine
%'MATLAB:nearlySingularMatrix'

%https://www.cia.gov/library/publications/the-world-factbook/geos/lo.html
popSize=5.440602;            % Population size in millions
NPop=popSize*1E6;        % Population size

%% Colors
blue   = [0      0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
yellow = [0.9290, 0.6940, 0.1250];
gray = [255 255 255]/255;

first=1;

fitPrev=inf;
fitBeginBest=1;


fitTestBegin=14;
fitTestSpan=31;

%fitTestSpan=31;

mod=1;
iterations=100;
method = ["gna","lsqnonlin"];
method = ["lsqnonlin"]
%method = ["fmincon"];
figure(101)


disp(['Day',' ','MSE','   ','FPE', '   AIC','    AICc','    nAIC'])
for i=1:1:length(method)
%gray=gray-1/fitTestSpan; % This is only for colors.
gray=gray*0.8;

disp(['-----------------'])
disp(['Fit method: ',method(i)])

for fitBegin=fitTestBegin:1:fitTestSpan;     % Day to begin the fit


[SEIQRDPm,InitialStates] = fitSEIRDQP2(fitBegin,iterations,mod,method(i));
disp([num2str(fitBegin),'   ',num2str(round(SEIQRDPm.Report.Fit.MSE)),'   ',num2str(round(SEIQRDPm.Report.Fit.FPE)),'   ',num2str(round(SEIQRDPm.Report.Fit.AIC)),'   ',num2str(round(SEIQRDPm.Report.Fit.AICc))])

udata = iddata([],zeros(2*365,0),1);
opt = simOptions('InitialCondition',InitialStates);
%sim(SEIQRDPinit ,udata,opt);
% RES.OutputData(end,1)-RES.OutputData(1,1);
% return
[Y A X]=sim(SEIQRDPm,udata,opt);



    
%% Simulate 

t=fitBegin:length(X)+fitBegin-1;
h1=semilogy(t,X(:,4),':','Color',gray,'LineWidth',0.5);
hold on;
h2=semilogy(t,X(:,5),':','Color',gray,'LineWidth',0.5);
h3=semilogy(t,X(:,6),':','Color',gray,'LineWidth',0.5);

set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


% If prediction does not fly off the charts
% Probably these are some local minimums, should check and compar 
% Parameters

%MSE 1E6
%FPE 1E5

% If protectin rate is not 0
if SEIQRDPm.Parameters(1).Value>eps
    
    %   9.0054e-04 alpha was pretty bad.
    
%No end in sight
if (X(end,4)<1E2) % This criterion should be improved
% If the error is smaller than the previous best
 if SEIQRDPm.Report.Fit.FPE < fitPrev
    fitBeginBest=fitBegin;
    fitPrev=SEIQRDPm.Report.Fit.FPE;
    methodBest=i;
 end
end

end %Check if not 0


end % Fit days
end % Methods



%% Best fit
runLength=2*365;

fitBegin=fitBeginBest;

[SEIQRDPm,InitialStates] = fitSEIRDQP2(fitBeginBest,iterations,mod,method(methodBest));
udata = iddata([],zeros(runLength,0),1);
opt = simOptions('InitialCondition',InitialStates);
%opt.AbsTol=1E-6;
%sim(SEIQRDPinit ,udata,opt);
% RES.OutputData(end,1)-RES.OutputData(1,1);
% return
[Y A X]=sim(SEIQRDPm,udata,opt);
figure(101)
t=fitBegin:length(X)+fitBegin-1;
semilogy(t,X(:,4),'-','Color',blue,'LineWidth',1);
hold on;
semilogy(t,X(:,5),'-','Color',orange,'LineWidth',1);
semilogy(t,X(:,6),'k-','LineWidth',1);

%pause
%% Short term, based on fit and current data
cd ..
importData;                              % Script to import the data from CSV
cd experimental

I=cumsum(Confirmed);      % Cumulative sum of daily cases, transpose to make it compatible w/ E. Cheynet's code
R=cumsum(Recovered);      % Cumulative sum of daily cases, transpose to make it compatible w/ E. Cheynet's code
D=cumsum(Deaths);         % Cumulative sum of daily cases, transpose to make it compatible w/ E. Cheynet's code
I=I-R-D;                   % Active infections


Snow = X(max(Day)-fitBegin,1);
Enow = X(max(Day)-fitBegin,2);
Inow = X(max(Day)-fitBegin,3);
Qnow = I(end);
Rnow = R(end);
Dnow = D(end);
Pnow = X(max(Day)-fitBegin,7);



udata = iddata([],zeros(14,0),1);
opt = simOptions('InitialCondition',[Snow,Enow,Inow,Qnow,Rnow,Dnow,Pnow]');
%sim(SEIQRDPinit ,udata,opt);
% RES.OutputData(end,1)-RES.OutputData(1,1);
% return
[Yst Ast Xst]=sim(SEIQRDPm,udata,opt);
InfectedTomorrow=round(Xst(2,4)+Xst(2,5)+Xst(2,6));

InfectingTomorrow=round(Xst(2,3));
ExposedTomorrow=round(Xst(2,2));
%% Report

%pause
clc

d1=datetime(2020,3,6,'Format','d.M'); % First confirmed case
DayLT=[1; 26; 56; 87; 117; 148; 179; 209; 240; 270; 301; 332; 361];
DateLT = {datestr(d1); datestr(d1+26); datestr(d1+56); datestr(d1+87); datestr(d1+117); datestr(d1+148); datestr(d1+179); datestr(d1+209); datestr(d1+240); datestr(d1+270); datestr(d1+301); datestr(d1+332); datestr(d1+361)};% Date array for predictions

disp(['SEIQRDP 7-stavový homogénný infektologický model bez vitálnej dynamiky'])
disp(['(Predikcia s parametrami na základe dostupných údajov.)'])
disp(['----------------------------'])
disp(['Overené prípady:             ',num2str(InfectedTomorrow),' (do konca dna)'])
disp(['Nové overené prípady:         ',num2str(InfectedTomorrow-(I(end)+D(end)+R(end))),' (do konca dna)'])
disp(['Celk. pocet infekcných:     ',num2str(InfectedTomorrow+InfectingTomorrow),' (celkové aktívne infekcie, odhad, ',num2str(InfectingTomorrow),' mimo testov)'])
disp(['Pocet nakazených:            ',num2str(ExposedTomorrow),' (v inkubacnej dobe, odhad)'])
[valMax indMax]=max(X(:,4));
disp(['Maximum aktívnych infekcií: ',num2str(round(valMax)),' (pre celu vlnu ochoreni, overených)'])
disp(['Vrchol overených infekcií:    ',datestr(d1+fitBegin+indMax)])
disp(['Infikovaní:                 ',num2str(round(max(X(:,5)))+round(max(X(:,6)))),' (pre celu vlnu ochoreni)'])
disp(['Vyliecení:                  ',num2str(round(max(X(:,5)))),' (pre celu vlnu ochoreni)'])
disp(['Úmrtia:                       ',num2str(round(max(X(:,6)))),' (pre celu vlnu ochoreni)'])
disp(['Koniec infekcií:            ',datestr(d1+fitBegin+min(find(X(max(Day):end,4)<10))),' (<10 aktívnych prípadov)'])


disp(['Miera ochrany alpha:           ',num2str(round(SEIQRDPm.Parameters(1).Value*1000)/1000),' '])

disp(['Infekcna doba 1/beta:         ',num2str(round(1/SEIQRDPm.Parameters(2).Value*100)/100),' [dni]'])
disp(['Inkubacna doba 1/sigma:       ',num2str(round(1/SEIQRDPm.Parameters(3).Value*100)/100),' [dni]'])
disp(['Izolacna doba 1/delta:     ',num2str(round(1/SEIQRDPm.Parameters(4).Value*100)/100),' [dni]'])
disp(['Doba vyliecenia 1/gamma0:     ',num2str(round(1/SEIQRDPm.Parameters(5).Value*100)/100),' [dni]'])
disp(['Miera umrtnosti                ',num2str(round(SEIQRDPm.Parameters(6).Value*100000)/100000),''])


%disp(['Miera odstránenia prípadov gamma: ',num2str(gammaEst),', krátkodobý fit'])
disp(['Zhoda modelu (infekcie):      ',num2str(round(SEIQRDPm.Report.Fit.FitPercent(1)*10)/10),'%'])
disp(['Zhoda modelu (vylieceni):     ',num2str(round(SEIQRDPm.Report.Fit.FitPercent(2)*10)/10),'%'])
disp(['Zhoda modelu (umrtia):        ',num2str(round(SEIQRDPm.Report.Fit.FitPercent(3)*10)/10),'%'])
disp(' ')
disp(['Najlepší fit od ',num2str(fitBeginBest),' dna nakazy v rozmedzi ',num2str(fitTestBegin),'-',num2str(fitTestSpan),' dna az dodnes, kde alpha>0.'])
first=0;


%% Finish plot

% Data reading and preparation
cd ..
importData;                              % Script to import the data from CSV
cd experimental 

I=cumsum(Confirmed);      % Cumulative sum of daily cases, transpose to make it compatible w/ E. Cheynet's code
R=cumsum(Recovered);      % Cumulative sum of daily cases, transpose to make it compatible w/ E. Cheynet's code
D=cumsum(Deaths);         % Cumulative sum of daily cases, transpose to make it compatible w/ E. Cheynet's code
I=I-R-D;                   % Active infections

semilogy(Day,I,'o','Color',blue,'MarkerSize',4,'LineWidth',0.1)
hold on
semilogy(Day,R,'o','Color',orange,'MarkerSize',4,'LineWidth',0.1)
semilogy(Day,D,'ko','MarkerSize',4,'LineWidth',0.1)

ylabel('Pocet pripadov')
xlabel('Cas (dni od prveho pripadu)')

set(gcf,'color','w')
grid on
axis tight
set(gca,'yscale','lin')
title('COVID-19 na Slovensku, SEIQRDP model')

xticks(DayLT)
xticklabels(DateLT)
xtickangle(90)
pandemicEnd=text(fitBegin+min(find(X(max(Day):end,4)<10)),50,datestr(d1+fitBegin+min(find(X(max(Day):end,4)<10))));
set(pandemicEnd,'Rotation',90)
set(pandemicEnd,'FontSize',8)

pandemicPeak=text(fitBegin+indMax,valMax*1.1,datestr(d1+fitBegin+indMax),'HorizontalAlignment','center');
set(pandemicPeak,'Rotation',0)
set(pandemicPeak,'FontSize',8)


leg = {'Infekcni (model)','Vylieceni (model)','Mrtvi (model)','Infekcni (data)','Vylieceni (data)','Mrtvi (data)'};
legend(leg{:},'location','northeastoutside')


axis([0,fitBegin+min(find(X(max(Day):end,4)<10))+20,0,max(X(:,5))*1.1])
%axis([0,270,0,max(X(:,5))*1.1])
cd out
print(['skCOVID19_SEIQRD_LongTerm'],'-dpng','-r0')
cd ..

axis([0,max(Day)*1.1,0,max(I)*1.1])

cd out
print(['skCOVID19_SEIQRD_LongTermFit'],'-dpng','-r0')
cd ..

axis([0,fitBegin+min(find(X(max(Day):end,4)<10))+20,0,max(X(:,5))*1.1])


set(gca,'yscale','log')

axis([0,fitBegin+min(find(X(max(Day):end,4)<10))+20,1,1E4])
%axis([0,270,0,max(X(:,5))*1.1])
cd out
print(['skCOVID19_SEIQRD_LongTerm_Log'],'-dpng','-r0')
cd ..

