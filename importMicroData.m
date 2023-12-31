function data = importMicroData(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  HFINALEVENTSCAMPOINCHAUSPEARGENTINA2849MZDATUM = IMPORTFILE(FILE)
%  reads data from the first worksheet in the Microsoft Excel
%  spreadsheet file named FILE.  Returns the data as a table.
%
%  HFINALEVENTSCAMPOINCHAUSPEARGENTINA2849MZDATUM = IMPORTFILE(FILE,
%  SHEET) reads from the specified worksheet.
%
%  HFINALEVENTSCAMPOINCHAUSPEARGENTINA2849MZDATUM = IMPORTFILE(FILE,
%  SHEET, DATALINES) reads from the specified worksheet for the
%  specified row interval(s). Specify DATALINES as a positive scalar
%  integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  hfinaleventsCampoInchauspeArgentina2849mZdatum = importfile("G:\Mi unidad\Proyecto Fracking\Microsismica\2022h_final-events_Campo-Inchauspe_Argentina2_849m-Zdatum.xlsx", "2022_final_all-times", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 23-Feb-2023 11:30:30

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
opts = spreadsheetImportOptions("NumVariables", 22);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = dataLines(1, :);

% Specify column names and types
opts.VariableNames = ["X", "Y", "Z", "DATE", "TIME", "AMP", "SNR", "SEMB", "FOCMECH", "STRIKE", "DIP", "RAKE", "UX", "UY", "UZ", "STAGE", "WELL", "AMTI", "Mw", "m0", "PostFrac", "Fault"];
opts.VariableTypes = ["double", "double", "double", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
data = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = dataLines(idx, :);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    data = [data; tb]; %#ok<AGROW>
end

end