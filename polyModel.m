
[fitresultPoly,gofPoly]=polyFit(Day(fitbegin:end), Nd(fitbegin:end));

aPoly=fitresultPoly.a; % Growth factor 
bPoly=fitresultPoly.b; % Correction for zero day cases
ciPoly=confint(fitresultPoly); % Confidence intervals at 95% confidentce
R2Poly=gofPoly.rsquare;

NdPredictedPoly=round(aPoly.*DayPred.^bPoly);
NdPredPolyHigh=round(ciPoly(2,1).*DayPred(max(Day):end).^ciPoly(2,2));
NdPredPolyLow=round(ciPoly(1,1).*DayPred(max(Day):end).^ciPoly(1,2));

%NdPredHigh=round((ci(2,1).^DayPredRest*NdPredicted(max(Day)-1)));
%NdPredLow=round((ci(1,1).^DayPredRest*NdPredicted(max(Day)-1)));


%% Shift data to account for onset of symptoms

NdSymptomsPoly=round(aPoly.*(DayPred+symptoms).^bPoly);
NdSymptomsPolyHigh=round(ciPoly(2,1).*(DayPred+symptoms).^ciPoly(2,2));
NdSymptomsPolyLow=round(ciPoly(1,1).*(DayPred+symptoms).^ciPoly(1,2));

% plot(NdSymptomsPoly)
% hold on
% plot(NdSymptomsPolyHigh)
% plot(NdSymptomsPolyLow)
% 
% return

%% First case

polyFun = @(x)ceil(aPoly*(x)^bPoly);
firstCasePoly = floor(fminbnd(polyFun, -100, 0));


%% Growth factor (not really, but)....

growthFactorPoly= ([NdPredictedPoly(1:length(Day)) 0]./[0 NdPredictedPoly(1:length(Day))]);
growthFactorPoly= (growthFactorPoly(2:end-1)-1)*100;
gFPoly=mean(growthFactorPoly(end-7:end));

%% Plot cases


h =  findobj('type','figure');
n = length(h);
figure(n+1)

%patch([DayPred(max(Day):end), DayPred(end:-1:max(Day)), DayPred(max(Day))],[NdPredPolyLow, NdPredPolyHigh(end:-1:1),NdPredPolyLow(1)],'r','EdgeAlpha',0,'FaceAlpha',0.2) % Confidence intervals
hold on
grid on
%patch([DayPred, DayPred(end:-1:1), DayPred(1)],[NdSymptomsPolyLow, NdSymptomsPolyHigh(end:-1:1),NdSymptomsPolyLow(1)],'y','EdgeAlpha',0,'FaceAlpha',0.2) % Confidence intervals

plot(Day,Nd,'o-','LineWidth',2,'Color',blue,'MarkerSize',6) % Confirmed cumulative cases
bar(Day,Confirmed) % Confirmed new cases
plot(DayPred,NdPredictedPoly,'Color',orange,'LineWidth',1.5) % Predicted cases

plot(DayPred,NdSymptomsPoly,'Color',yellow,'LineWidth',1) % Predicted Shifted Cases

plot(dataSKpred(:,1),dataSKpred(:,10),'k','Marker','v','MarkerSize',6) % Predicted Shifted Cases

% Previous predictions
%plot(dataSKpred(:,1),dataSKpred(:,2),'k.')
%errorbar(dataSKpred(:,1),dataSKpred(:,2),dataSKpred(:,2)-dataSKpred(:,3),dataSKpred(:,2)-dataSKpred(:,4),'k')

xticks(DayPred)
xticklabels(DatePred)
xtickangle(90)
xlabel('Date')
ylabel('Cases')
%legend('95% Confidence (Confirmed prediction)','95% Confidence (Total, 5 d shift)','Cumulative confirmed','New confirmed',['Poly. Approximation (Confirmed) R^2=',num2str(R2)],['Exp. Approximation (Total), 5 d shift'],'Location','northwest')
legend('Cumulative confirmed','New confirmed',['Poly. Approximation (Confirmed) R^2=',num2str(R2Poly)],['Poly. Approximation (Total), 5 d shift'],'Predictions (1-day)','Location','northwest')
title(['SARS-CoV-2 Cases in Slovakia: Polynomial model, ',datestr(dt)])
axis([0,length(Day)+7,0,NdPredictedPoly(max(Day)+7)])
text(0.5,0.9,'covid19.gergelytakacs.com','FontSize',10,'rotation',90,'Color',[0.7 0.7 0.7])
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 20 10];

cd out
print(['skCOVID19_Poly_Cases_',datestr(dt)],'-dpng','-r0')
%print(['skCOVID19_Poly_Cases_',datestr(dt)],'-dpdf','-r0')

% axis([0,length(Day)+1,0,dataSKpred(end,4)])
% 
% print(['skCOVID19_Poly_Cases_Detail',datestr(dt)],'-dpng','-r0')
% %print(['skCOVID19_Poly_Cases_Detail',datestr(dt)],'-dpdf','-r0')
cd ..
