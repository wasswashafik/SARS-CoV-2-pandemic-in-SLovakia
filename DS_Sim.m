% DS
N=22768;
lambda=195/dayYear;
mu=258/dayYear;

run=150

%% Simulate for generic data

%SEIR_VD.Parameters(3).Value=lambda;
% 
% SIR_VD.Parameters(4).Maximum =inf;
% SIR_VD.Parameters(4).Minimum = -inf;
% 
%SEIR_VD.Parameters(4).Value=mu;

udata = iddata([],zeros(run,0),1);
opt = simOptions('InitialCondition',[N-1;1;0;0]);
simdatads=sim(SEIR_VD,udata,opt);
IsimDS=simdatads.OutputData(:,2);
RsimDS=simdatads.OutputData(:,3);
SsimDS=simdatads.OutputData(:,1);
t=1:1:run;
Itot=cumsum(IsimDS);

%% Plot




hold on

subplot(3,1,1)

plot(t,SsimDS)
title('Susceptible')
grid on
xlabel('Days')
ylabel('Cases')

subplot(3,1,2)
grid on
plot(t,IsimDS)
hold on
%plot(t,Itot)
title('Infected')
grid on
xlabel('Days')
ylabel('Cases')

subplot(3,1,3)
grid on
plot(t,RsimDS)
title('Removed')
grid on
xlabel('Days')
ylabel('Cases')
