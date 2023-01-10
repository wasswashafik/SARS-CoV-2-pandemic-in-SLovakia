close all;


%% Cases (Exponential)

h =  findobj('type','figure');
n = length(h);
figure(n+1)

patch([DayPred(max(Day):end), DayPred(end:-1:max(Day)), DayPred(max(Day))],[NdPredLow, NdPredHigh(end:-1:1),NdPredLow(1)],'r','EdgeAlpha',0,'FaceAlpha',0.2) % Confidence intervals
hold on
grid on
patch([DayPred, DayPred(end:-1:1), DayPred(1)],[NdSymptomsLow, NdSymptomsHigh(end:-1:1),NdSymptomsLow(1)],'y','EdgeAlpha',0,'FaceAlpha',0.2) % Confidence intervals

plot(Day,Nd,'o-','LineWidth',2,'Color',blue,'MarkerSize',6) % Confirmed cumulative cases
bar(Day,Confirmed) % Confirmed new cases
plot(DayPred(fitbegin:end),NdPredicted(fitbegin:end),'Color',orange,'LineWidth',1.5) % Predicted cases

plot(DayPred,NdSymptoms,'Color',yellow,'LineWidth',1) % Predicted Shifted Cases

% Previous predictions
plot(dataSKpred(:,1),dataSKpred(:,2),'k.','MarkerSize',10)
errorbar(dataSKpred(:,1),dataSKpred(:,2),dataSKpred(:,2)-dataSKpred(:,3),dataSKpred(:,2)-dataSKpred(:,4),'k')

xticks(DayPred)
xticklabels(DatePred)
xtickangle(90)
xlabel('Date')
ylabel('Cases')
legend('95% Confidence (Confirmed prediction)','95% Confidence (Total, 5 d shift)','Cumulative confirmed','New confirmed',['Exp. Approximation (Confirmed) R^2=',num2str(R2)],['Exp. Approximation (Total), 5 d shift'],'Location','northwest')
title(['SARS-CoV-2 Cases in Slovakia: Exponential model, ',datestr(dt)])

text(0.5,0.9,'covid19.gergelytakacs.com','FontSize',10,'rotation',90,'Color',[0.7 0.7 0.7])
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 20 10];

axis([0,length(Day)+7,0,NdSymptoms(length(Day)+7)])

cd out
print(['skCOVID19_Exp_Cases_',datestr(dt)],'-dpng','-r0')
%print(['skCOVID19_Exp_Cases_',datestr(dt)],'-dpdf','-r0')

axis([0,length(Day)+1,0,dataSKpred(end,4)])

print(['skCOVID19_Exp_Cases_Detail',datestr(dt)],'-dpng','-r0')
%print(['skCOVID19_Exp_Cases_Detail',datestr(dt)],'-dpdf','-r0')
cd ..


%% Growth factor

h =  findobj('type','figure');
n = length(h);
figure(n+1)

% Previous predictions
plot(dataSKpred(:,1)-1,(dataSKpred(:,5)-1)*100,'k.','MarkerSize',10)
hold on
errorbar(dataSKpred(:,1)-1,(dataSKpred(:,5)-1)*100,(dataSKpred(:,6)-1)*100-(dataSKpred(:,5)-1)*100,(dataSKpred(:,7)-1)*100-(dataSKpred(:,5)-1)*100,'k')
plot(dataSKpred(:,1)-1,(dataSKpred(:,8)-1)*100,'kx')
plot(dataSKpred(:,1)-1,(dataSKpred(:,11)),'kv')



hold on
grid on

plot(Day(2:end),growthFactor,'.-','LineWidth',2) % Predicted Shifted Cases

xticks(DayPred)
xticklabels(DatePred)
xtickangle(90)
xlabel('Date')
ylabel('Growth factor [%]')
legend('Growth factor (exponential fit)','Growth factor (exponential fit, conf.)','Growth factor (SEIR)','Growth factor (polynomial fit)','Growth factor (data)','Location','northeast')
title(['SARS-CoV-2 Growth factor in Slovakia, ',datestr(dt)])
axis([12,length(Day)+1,0,50])
text(2.5,0.9,'covid19.gergelytakacs.com','FontSize',10,'rotation',90,'Color',[0.7 0.7 0.7])
cd out
print(['skCOVID19_GrowthFactor_',datestr(dt)],'-dpng','-r0')
%print(['skCOVID19_GrowthFactor_',datestr(dt)],'-dpdf','-r0')
cd ..

%% Testing

h =  findobj('type','figure');
n = length(h);
figure(n+1)

hold on
plot(Day,totTest,'.-','LineWidth',2) % Predicted Shifted Cases
bar(Day(2:end),newTest) % Confirmed new cases
plot(Day(end-7+1:end),testA.*Day(end-7+1:end)+newTest(end-7+1),'k--','LineWidth',1.5)
grid on

xticks(DayPred)
xticklabels(DatePred)
xtickangle(90)
xlabel('Date')
ylabel('Tests')
legend('Total tests','New Tests','Trend','Location','northwest')
title(['SARS-CoV-2 Testing in Slovakia, ',datestr(dt)])
axis([1,max(Day)+1,0,max(totTest)])
text(1.5,100,'covid19.gergelytakacs.com','FontSize',10,'rotation',90,'Color',[0.7 0.7 0.7])
cd out
print(['skCOVID19_Tests_',datestr(dt)],'-dpng','-r0')
%print(['skCOVID19_Tests_',datestr(dt)],'-dpdf','-r0')
cd ..

%% Residuals
% 
% figure(4)
% 
% patch([0 max(Day), max(Day) 0 ],[mean(NdResiduals)+std(NdResiduals) mean(NdResiduals)+std(NdResiduals),mean(NdResiduals)-std(NdResiduals) mean(NdResiduals)-std(NdResiduals)],'b','EdgeAlpha',0,'FaceAlpha',0.2) % Confidence intervals
% patch([0 max(Day), max(Day) 0 ],[mean(NdResiduals)+2*std(NdResiduals) mean(NdResiduals)+2*std(NdResiduals),mean(NdResiduals)-2*std(NdResiduals) mean(NdResiduals)-2*std(NdResiduals)],'b','EdgeAlpha',0,'FaceAlpha',0.1) % Confidence intervals
% 
% patch([0 max(Day), max(Day) 0 ],[mean(Isim(1:max(Day))-Nd)+std(Isim(1:max(Day))-Nd) mean(Isim(1:max(Day))-Nd)+std(Isim(1:max(Day))-Nd),mean(Isim(1:max(Day))-Nd)-std(Isim(1:max(Day))-Nd) mean(Isim(1:max(Day))-Nd)-std(Isim(1:max(Day))-Nd)],'r','EdgeAlpha',0,'FaceAlpha',0.2) % Confidence intervals
% patch([0 max(Day), max(Day) 0 ],[mean(Isim(1:max(Day))-Nd)+2*std(Isim(1:max(Day))-Nd) mean(Isim(1:max(Day))-Nd)+2*std(Isim(1:max(Day))-Nd),mean(Isim(1:max(Day))-Nd)-2*std(Isim(1:max(Day))-Nd) mean(Isim(1:max(Day))-Nd)-2*std(Isim(1:max(Day))-Nd)],'r','EdgeAlpha',0,'FaceAlpha',0.1) % Confidence intervals
% 
% hold on
% plot(Day,NdResiduals,'.-','LineWidth',2) % Predicted Shifted Cases
% grid on
% plot(Day,Isim(1:max(Day))-Nd,'.-','LineWidth',2) % Predicted Shifted Cases
% 
% xticks(DayPred)
% xticklabels(DatePred)
% xtickangle(90)
% xlabel('Date')
% ylabel('Cases (difference)')
% legend('1 sigma (Exponential)','2 sigma (Exponential)','1 sigma (SEIR)','2 sigma (SEIR)','Case residuals (Exponential)','Case residuals (SEIR)','Location','northwest')
% title(['SARS-CoV-2 Case prediction residuals for Slovakia, ',datestr(dt)])
% %axis([1,max(Day)+1,0,max(totTest)])
% text(0.5,-40,'covid19.gergelytakacs.com','FontSize',10,'rotation',90,'Color',[0.7 0.7 0.7])
% cd out
% print(['skCOVID19_Residuals_',datestr(dt)],'-dpng','-r0')
% print(['skCOVID19_Residuals_',datestr(dt)],'-dpdf','-r0')
% cd ..

%% SEIR Output (Infected)

h =  findobj('type','figure');
n = length(h);
figure(n+1)

plot(Day,I,'-o','MarkerSize',6,'LineWidth',2)
hold on
grid on
bar(Day,Confirmed) % Confirmed new cases
plot(DayPred(max(Day):end-1),IsimPred,'--','Color',blue,'LineWidth',1.5)      % Prediction only
plot(Day(SIR_fitBegin:end),IsimFit,'LineWidth',1.5,'Color',orange)      % Prediction only
plot(-symptoms+fitbegin:1:(length(IsimSymptoms)-symptoms+fitbegin-1),IsimSymptoms,'Color',yellow,'LineWidth',1)
plot(dataSKpred(:,1),dataSKpred(:,9),'kx','MarkerSize',10)

xticks(DayPred)
xticklabels(DatePred)
xtickangle(90)
xlabel('Date')
ylabel('Cases')
legend('Infections (data)','New infections (data)','Predicted infections (SEIR model)','SEIR model fit','Total Infected (SEIR model, 5d shift)','Predictions','Location','northwest')
title(['SARS-CoV-2 Infections in Slovakia: SEIR model w/ vital dynamics, ',datestr(dt)])
axis([0,length(Day)+7,0,NdSymptoms(length(Day)+7)])
text(0.5,0.9,'covid19.gergelytakacs.com','FontSize',10,'rotation',90,'Color',[0.7 0.7 0.7])
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 20 10];

cd out
print(['skCOVID19_SEIR_Cases_',datestr(dt)],'-dpng','-r0')
%print(['skCOVID19_SEIR_Cases_',datestr(dt)],'-dpdf','-r0')

axis([0,length(Day)+1,0,dataSKpred(end,4)])

print(['skCOVID19_SEIR_Cases_Detail',datestr(dt)],'-dpng','-r0')
%print(['skCOVID19_SEIR_Cases_Detail',datestr(dt)],'-dpdf','-r0')
cd ..


%% Semilogy

h =  findobj('type','figure');
n = length(h);
figure(n+1)

semilogy(Day,Nd,'o-','LineWidth',2,'Color',blue,'MarkerSize',6) % Confirmed cumulative cases
hold on
grid on

plot(DayPred,NdPredictedPoly,'Color',orange,'LineWidth',1.5) % Predicted cases

plot(DayPred,NdPredicted,'Color',orange,'LineWidth',1,'LineStyle',':') % Predicted cases

%plot(DayPred,NdPredictedNoN0,'Color',orange,'LineWidth',1,'LineStyle','--') % Predicted cases


plot(Day(fitbegin:end),IsimFit,'LineWidth',0.5,'LineStyle','--','Color',orange) % Predicted cases

%plot(DayPred,NdSymptomsPoly,'Color',yellow,'LineWidth',1) % Predicted Shifted Cases

%plot(dataSKpred(:,1),dataSKpred(:,10),'k','Marker','v','MarkerSize',6) % Predicted Shifted Cases

% Previous predictions
%plot(dataSKpred(:,1),dataSKpred(:,2),'k.')
%errorbar(dataSKpred(:,1),dataSKpred(:,2),dataSKpred(:,2)-dataSKpred(:,3),dataSKpred(:,2)-dataSKpred(:,4),'k')

xticks(DayPred)
xticklabels(DatePred)
xtickangle(90)
xlabel('Date')
ylabel('Cases')
%legend('95% Confidence (Confirmed prediction)','95% Confidence (Total, 5 d shift)','Cumulative confirmed','New confirmed',['Poly. Approximation (Confirmed) R^2=',num2str(R2)],['Exp. Approximation (Total), 5 d shift'],'Location','northwest')
legend('Cumulative confirmed','Polynomial fit','Exponential fit','SEIR Model','Location','southeast')
title(['SARS-CoV-2 Cases in Slovakia: Logarithmic view and comparison, ',datestr(dt)])
axis([0,length(Day)+1,0,NdPredictedPoly(max(Day)+1)])
text(0.5,0.9,'covid19.gergelytakacs.com','FontSize',10,'rotation',90,'Color',[0.7 0.7 0.7])
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 20 10];

cd out
print(['skCOVID19_Logarithmic_Cases_',datestr(dt)],'-dpng','-r0')
axis([12,length(Day)+1,100,dataSKpred(end,4)])
print(['skCOVID19_Logarithmic_Cases_Detail',datestr(dt)],'-dpng','-r0')

% print(['skCOVID19_Poly_Cases_Detail',datestr(dt)],'-dpng','-r0')
% %print(['skCOVID19_Poly_Cases_Detail',datestr(dt)],'-dpdf','-r0')
cd ..


%% Tests per positive

%plot(Day(2:end),newTest/Confirmed(2:end))



