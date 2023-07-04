function ETS1 = importSheet2ET2021h(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  ETS1 = IMPORTFILE(FILE) reads data from the first worksheet in the
%  Microsoft Excel spreadsheet file named FILE.  Returns the data as a
%  table.
%
%  ETS1 = IMPORTFILE(FILE, SHEET) reads from the specified worksheet.
%
%  ETS1 = IMPORTFILE(FILE, SHEET, DATALINES) reads from the specified
%  worksheet for the specified row interval(s). Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  ETS1 = importfile("G:\Mi unidad\Proyecto Fracking\Microsismica\ET.xp-2021 Pilot_Survey_Run10_8.5in_3330.07m.xlsx", "Sheet2", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 17-May-2023 18:12:02

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 3);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = dataLines(1, :);

% Specify column names and types
opts.VariableNames = ["MD", "DEVI", "HAZI"];
opts.VariableTypes = ["double", "double", "double"];

% Import the data
ETS1 = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = dataLines(idx, :);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    ETS1 = [ETS1; tb]; %#ok<AGROW>
end

end