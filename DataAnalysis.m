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
% Updated: August 31, 2009 by Vasco Curdia
%
% Copyright 2008-2011 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% display
fprintf('Loading and analyzing data...\n')

%% Set Timer
TimeElapsed.DataAnalysis = toc();

%% load the data
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
TimeElapsed.DataAnalysis = toc-TimeElapsed.DataAnalysis;

%% -----------------------------------------------------------------------------
