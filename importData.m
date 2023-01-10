%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Day", "Date", "Confirmed", "Deaths","negTest","Recovered"];
opts.VariableTypes = ["double", "string", "double", "double","double","double"];
opts = setvaropts(opts, 2, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 2, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable("data-sk.csv", opts);

%% Convert to output type
Day = tbl.Day;
Date = tbl.Date;
Confirmed = tbl.Confirmed;
Deaths = tbl.Deaths;
negTest = tbl.negTest;
Recovered = tbl.Recovered;

%% Clear temporary variables
clear opts tbl