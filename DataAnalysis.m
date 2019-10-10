% DataAnalysis
%
% Analyze the data used for the estimation
%
% See also:
% SetDSGE, GenSymVars, DataAnalysis, PriorAnalysis, GenPost, MaxPost,
% MaxPostFcn, MakeTableMaxPost, MCMC, MCMCFcn, MCMCSearchScaleFactor,
% MakePlotsMCMCDraws, MCMCInference, MakeTableMCMCInference, 
% MakePlotsMCMCTrace, MakePlotsMCMCPriorPost, MCMCConv, MakeTableMCMCConv,
% MakeReportMCMCPlots, MakeReportMCMC, TimeIdxCreate
%
% ..............................................................................
%
% Created: March 18, 2008 by Vasco Curdia
% Updated: August 6, 2015 by Vasco Curdia
%
% Copyright 2008-2015 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% display
fprintf('Loading and analyzing data...\n')

%% Set Timer
tt.start('DataAnalysis')

%% Detect file type
isCSV = strcmp(FileName.Data(end-3:end),'.csv');

%% Load data
if isCSV
    Data = importdata(FileName.Data);
    Data.Var = Data.textdata(1,2:end);
    Data.TimeIdx = Data.textdata(2:end,1)';
    Data.TimeStart = Data.TimeIdx{1};
    Data.TimeEnd = Data.TimeIdx{end};
    Data.nVar = length(Data.Var);
    Data.T = length(Data.TimeIdx);
    Data = rmfield(Data,'textdata');
    [tfDates,idxDates] = ismember(TimeIdx,Data.TimeIdx);
    if ~all(tfDates)
        error('DateLabels require time periods outside data file!')
    end
    [tfVar,idxVar] = ismember(ObsVar,Data.Var);
    if ~all(tfVar)
        error('Observable names do not match csv headers!')
    end
    Data = Data.data(idxDates,idxVar);
else
    % if CSV not specified in file name, assume that it is mat file
    Data = load(FileName.Data);
    if ~exist('DataVarName','var')
        Data.VarNames = fieldnames(Data);
        if length(Data.VarNames)==1
            DataVarName = Data.VarNames{1};
        else
            error('Too many variables in mat file. Need to specify which one to use.')
        end
    end
    Data = Data.(DataVarName);
end

%% analyze data
[nRow,nCol] = size(Data);
if nCol~=nObsVar
    if T~=nCol
        error('Data size does not match observables')
    end
    T = nCol; 
    Data = Data';
else
    T = nRow;
end
clear nRow nCol

%% Elapsed time
tt.stop('DataAnalysis')

%% -----------------------------------------------------------------------------
