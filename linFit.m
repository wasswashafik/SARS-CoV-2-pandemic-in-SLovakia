function [fitresult, gof] = linFit(testDay, newTest)


[xData, yData] = prepareCurveData( testDay, newTest );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData(end-7+1:end), xData(end-7+1:end), ft );




