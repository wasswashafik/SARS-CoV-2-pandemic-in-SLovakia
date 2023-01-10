function [fitresult, gof] = expFitStart(Day, Nd)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( Day, Nd );

% Set up fittype and options.
ft = fittype( '(a^x)*b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
%opts.StartPoint = [0.0274677556889124 0.304255006873599];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );



