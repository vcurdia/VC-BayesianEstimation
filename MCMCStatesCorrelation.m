% MCMCStatesCorrelation
%
% Computes the correlations across states as implied by the data.
%
% NOTE: this is different from the model implied correlations, unless we
%       have a very large dataset.
%
% Required settings:
%
%   nChains
%   Sets the number of MCMC chains to analyze.
%
%   FileName.MCMCDraws
%   Name of the file with the draws.
%
% Optional settings:
%
%   BurnIn
%   Percentage of draws to burn.
%   Default: 0.5
%
%   nDrawsStates
%   Number of draws to use.
%   Default: 1000
%
%   Prc2Compute
%   Percentiles to compute.
%   Default: [5,50,95]
%
%   Vars2Show
%   List of variables for which to show the correlations.
%
% See also:
% SetDSGE, GenSymVars, DataAnalysis, PriorAnalysis, GenPost, MaxPost,
% MaxPostFcn, MakeTableMaxPost, MCMC, MCMCFcn, MCMCSearchScaleFactor,
% MakePlotsMCMCDraws, MCMCInference, MakeTableMCMCInference, 
% MakePlotsMCMCTrace, MakePlotsMCMCPriorPost, MCMCConv, MakeTableMCMCConv,
% MakeReportMCMCPlots, MakeReportMCMC
%
% .........................................................................
%
% Created: February 4, 2009 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2009-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Check required settings
isError = 0;
if ~exist('nChains','var'),fprintf('Warning: nChains not specified.\n'),isError=1;end
if ~isfield(FileName,'MCMCDraws'),fprintf('Warning: FileName.MCMCDraws not specified.\n'),isError=1;end
if isError, disp('Warning: Errors found. Cannot Continue.'),return,end

%% Settings
if ~exist('BurnIn','var'),BurnIn=0.5; end
if ~exist('nDrawsStates','var'),nDrawsStates=1000; end
if ~exist('Prc2Compute','var'),Prc2Compute=[5,50,95]; end
if ~exist('Vars2Show','var'),Vars2Show=VarStates; end

%% ------------------------------------------------------------------------

%% Display
fprintf('\n***************************')
fprintf('\n* MCMC State Correlations *')
fprintf('\n***************************\n')

%% Set Timer
TimeStr = sprintf('MCMCStateCorrelations');
TimeElapsed.(TimeStr) = toc;

%% load the mcmc draws
fprintf('\nLoading data...\n')

for jChain=1:nChains
    load([FileName.MCMCDraws,int2str(jChain)],'xDraws','nDraws')
    nThinningStates = floor(nDraws*nChains/nDrawsStates*(1-BurnIn));
    xd(:,:,jChain) = xDraws(:,BurnIn*nDraws+1:nThinningStates:end);
end
clear xDraws
nDrawsUsed = size(xd,2)*nChains;
xd = reshape(xd,np,nDrawsUsed);
fprintf('Total number of draws per chain: %.0f\n', nDraws)
fprintf('Thinning interval: %.0f\n', nThinningStates)
fprintf('Burn in: %.0f%%\n', 100*BurnIn)
fprintf('Total number of draws used: %.0f\n', nDrawsUsed)

%% draw states
StateVarCorrd = zeros(nStateVar,nStateVar,nDrawsUsed);
for j=1:nDrawsUsed
    [postj,StateVarttj,SIGttj,G1j,G2j,Hj,ObsVarBarj]=feval(FileName.Post,xd(:,j));
    StateVarj = DrawStates(StateVarttj,SIGttj,G1j,G2j,T,nStateVar);
    StateVarVarj = cov(StateVarj');
    StateVarCorrd(:,:,j) = diag(diag(StateVarVarj).^(-1/2))*StateVarVarj*diag(diag(StateVarVarj).^(-1/2));
end

%% Get percentiles
namelength = [cellfun('length',Vars2Show)];
namelengthmax = max(namelength);
nVars2Show = length(Vars2Show);
for j=1:length(Prc2Compute)
    PrcLabel = sprintf('Prc%03.0f',10*Prc2Compute(j));
    StateVarPosteriorCorr.(PrcLabel) = prctile(StateVarCorrd,j,3);
    fprintf('\nCorrelation matrix (percentile %.1f):\n\n',Prc2Compute(j))
    fprintf([repmat(' ',1,namelengthmax),...
        repmat(['  %',int2str(max(namelengthmax,5)),'s'],1,nVars2Show),'\n'],Vars2Show{:})
    for jStateVar=1:nVars2Show
        fprintf(['%',int2str(namelengthmax),'s',...
            repmat(['  %',int2str(max(namelengthmax,5)),'.3f'],1,nVars2Show),'\n'],...
            Vars2Show{jStateVar},StateVarPosteriorCorr.(PrcLabel)(jStateVar,:))
    end
end

%% Clean up
clear xd StateVarCorrd nThinningStates

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------
