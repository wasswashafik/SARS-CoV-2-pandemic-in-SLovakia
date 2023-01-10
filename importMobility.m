%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Day", "Retail", "Grocery", "Parks", "Transit", "Work", "Residential"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
%opts.EmptyLineRule = "read";

% Import the data
datamobility = readtable("C:\Users\Gergely\Dropbox\SK-COVID-19\data-mobility.csv", opts);


%% Clear temporary variables
clear opts

Retail=datamobility.Retail/100;
Grocery=datamobility.Grocery/100;
Parks=datamobility.Parks/100;
Transit=datamobility.Transit/100;
Work=datamobility.Work/100;
Residential=datamobility.Residential/100;

Mobility=(Residential-Retail-Work-Transit-Grocery)./(max(Residential-Retail-Work-Transit-Grocery));
Mobility(Mobility<=0)=eps;
save Mobility Mobility