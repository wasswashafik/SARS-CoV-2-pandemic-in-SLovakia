function [fitresult, gof] = polyFit(Day, Nd)

[xData, yData] = prepareCurveData( Day, Nd );

% Set up fittype and options.
ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [1.16532858718078 1.7021157845075];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );



