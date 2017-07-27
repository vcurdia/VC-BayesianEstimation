% MakePlotsStates
%
% Plots the evolution of the different state variables in the model.
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
%   FigShapeStates (cell array)
%   shape of each figure. 1st element is the number of plots per row and
%   the 2nd is the number of plots per column.
%   Default: {1,1}
%
%   Bands2Show
%   Percent interevals to be shown in the plots, centered around the
%   median.
%   Default: [50,60,70,80,90]
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
% Created: October 21, 2008 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: July 29, 2011 by Vasco Curdia
%          Change the way parallel works
% Updated: September 26, 2011 by Vasco Curdia
%          Named States, instead of MCMCStates
%          Creates plots in folder.
% Updated: April 29, 2012 by Vasco Curdia
%          Allow for choice of whether to draw states or not, useful in case we
%          use the median or mode -- in this case it does not make sense to draw
%          the states, as it would give a single draw, hence potentially
%          misleading.
% 
% Copyright 2008-2012 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Settings
if ~exist('FigShapeStates','var'),FigShapeStates={1,1};end
if ~exist('ShowFig','var'),ShowFig=0; end
useDateLabels = (exist('DateLabels','var') && isfield(DateLabels,'XTick')...
    && isfield(DateLabels,'XTickLabels'));
if ~exist('nPreSample','var'), nPreSample = 0; end
if ~exist('UseMode','var'),UseMode=0; end
if ~exist('UseMedian','var'),UseMedian=0; end
if UseMode && UseMedian, error('Cannot use both median and mode.'), end
if ~exist('isDrawStates','var'),isDrawStates = UseMode || UseMedian; end
if ~(UseMode || UseMedian)
    if ~exist('nDrawsStates','var'),nDrawsStates=1000; end
    if ~exist('BurnIn','var'),BurnIn = 0.25; end
    if ~exist('nThinning','var'),nThinning = 1; end
    if ~exist('Bands2Show','var'),Bands2Show=[50,60,70,80,90]; end
end
if ~exist('PlotDir','var') || ~isfield(PlotDir,'States')
  PlotDir.States = 'Plots/States/';
end
if~isdir(PlotDir.States),mkdir(PlotDir.States),end
FileName.PlotsStates = sprintf('%sPlotsStatesMCMCUpdate%.0f',FileName.Output,nUpdate);

%% Check required settings
if ~(UseMode || UseMedian)
    isError = 0;
    if ~exist('nChains','var'),fprintf('Warning: nChains not specified.\n'),isError=1;end
    if ~isfield(FileName,'MCMCDraws'),fprintf('Warning: FileName.MCMCDraws not specified.\n'),isError=1;end
    if isError, disp('Warning: Errors found. Cannot Continue.'),return,end
end

%% ------------------------------------------------------------------------

%% Display
fprintf('\n****************')
fprintf('\n* States Plots *')
fprintf('\n****************\n')

%% Set Timer
TimeStr = strrep(FileName.PlotsStates,FileName.Output,'');
TimeElapsed.(TimeStr) = toc;

%% load the mcmc draws
if ~(UseMode || UseMedian)
    fprintf('\nLoading data...\n')
    nThinningStates = max(floor(nDraws*nChains/nDrawsStates*(1-BurnIn)),1);
    for jChain=1:nChains
        load([FileName.MCMCDraws,int2str(jChain)],'xDraws','nDraws')
        xd(:,:,jChain) = xDraws(:,BurnIn*nDraws+1:nThinningStates:end);
    end
    clear xDraws
    nDrawsUsed = size(xd,2)*nChains;
    xd = reshape(xd,np,nDrawsUsed);
    fprintf('Total number of draws per chain: %.0f\n', nDraws)
    fprintf('Thinning interval: %.0f\n', nThinningStates)
    fprintf('Burn in: %.0f%%\n', 100*BurnIn)
    fprintf('Total number of draws used: %.0f\n', nDrawsUsed)
end

%% Using mode
if UseMode
    xd = [Params(:).postmode]';
end

%% Using median
if UseMedian
    xd = [Params(:).postp500]';
end

%% draw states
nDrawsUsed = size(xd,2);
StateVard = zeros(nStateVar,T,nDrawsUsed);
MatFileName = FileName.Mats;
sd = sum(100*clock);
parfor j=1:nDrawsUsed
  StateVard(:,:,j) = MakePlotsStatesFcn(xd(:,j),MatFileName,Data,nStateVar,...
    nObsVar,nShockVar,T,sd+j,isDrawStates);
end

%% prepare reference variables
nPlotsStates = FigShapeStates{1}*FigShapeStates{2};
nFigStates = ceil(nStateVar/nPlotsStates);

%% Draw plots
for jF=1:nFigStates
    if ShowFig
        figure
    else
        figure('Visible','off')
    end
    for jf=1:nPlotsStates 
        jStateVar = (jF-1)*nPlotsStates+jf;
        if jStateVar>nStateVar, break, end
        subplot(FigShapeStates{:},jf)
        if UseMode || UseMedian
            plot(StateVard(jStateVar,nPreSample+1:end))
        else
            vcPlotDistBands(squeeze(StateVard(jStateVar,nPreSample+1:end,:))','Bands2Show',Bands2Show);
        end
        title(StateVar{jStateVar})
        if useDateLabels
            set(gca,...
                'XTick',DateLabels.XTick,'XTickLabel',DateLabels.XTickLabels,...
                'FontSize',8)
        end
    end
    if FigShapeStates{1}==1 && FigShapeStates{2}==1
        FigTxt = StateVar{jStateVar};
    else
        FigTxt = sprintf('%.0f',jF);
    end
    vcPrintPDF(sprintf('%s%sFig%s',PlotDir.States,FileName.PlotsStates,FigTxt))
    clear FigTxt
end

%% Clean up
clear xd StateVard nBands Bands2Show nStatesShow nPlotsStates Band BandPath BandColor StateVarjj BandPath

%% close figures
if ~ShowFig, close all, end

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------
