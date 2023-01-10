t=0:0.1:300;
gamma=6.391;


gammaTau=10;
gammaT2=1/10*(1-exp(-t/(gammaTau)));


plot(t,gammaT2)

legend('Original','Replacement')