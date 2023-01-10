function [fitresult, gof] = expFit(Day, Nd)

[xData, yData] = prepareCurveData( Day, Nd );

% Set up fittype and options.
ft = fittype( 'a^x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.45763626173227;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );



