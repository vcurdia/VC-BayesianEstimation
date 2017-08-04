% MCMCAnalysis
%
% Runs all the tests to analyze the MCMCDraws
%
% See also:
% SetDSGE, GenSymVars, DataAnalysis, PriorAnalysis, GenPost, MaxPost,
% MaxPostFcn, MakeTableMaxPost, MCMC, MCMCFcn, MCMCSearchScaleFactor,
% MakePlotsMCMCDraws, MCMCInference, MakeTableMCMCInference, 
% MakePlotsMCMCTrace, MakePlotsMCMCPriorPost, MCMCConv, MakeTableMCMCConv,
% MakeReportMCMCPlots, MakeReportMCMC, MCMCAnalysis
%
% .........................................................................
%
% Created: February 4, 2009 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: April 26, 2012 by Vasco Curdia
%          Allows for the exclusion of some tests.
% 
% Copyright 2009-2012 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Options
if ~exist('isMakePlotsStates','var'),isMakePlotsStates=1;end
if ~exist('isMakeIRF','var'),isMakeIRF=1;end
if ~exist('isMakeVD','var'),isMakeVD=1;end
  
%% Run tests
Tests2Run = {...
    'MakeMCMCDrawsRedux';
    'MCMCInference';
    'MakePlotsPriorPost';
    'MakePlotsMCMCDraws';
    'MakePlotsMCMCTrace';
    'MCMCConv';
    };
if isMakePlotsStates
  Tests2Run{end+1} = 'MakePlotsStates';
end
if isMakeIRF
  Tests2Run{end+1} = 'MakeIRF';
end
if isMakeVD
  Tests2Run{end+1} = 'MakeVD';
end
for jT2R=1:length(Tests2Run)
    eval(Tests2Run{jT2R})
    save tmpMCMCAnalysis
end
MakeReportMCMC
delete tmpMCMCAnalysis.mat

