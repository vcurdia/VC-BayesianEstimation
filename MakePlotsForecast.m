% MakePlotsForecast
%
% Plots the evolution of the different state variables in the model beyond
% the last period.
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
%   nForecast
%   Number of periods to forecast. Default: 4
%
%   tStart
%   Period in which to start the plot. Default: 1
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
% Created: April 26, 2009 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: July 29, 2011 by Vasco Curdia
%          Change the way parallel works
% Updated: September 26, 2011 by Vasco Curdia
%          Named Forecast, instead of MCMCForecast
%          Creates plots in folder.
% 
% Copyright 2009-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Check required settings
isError = 0;
if ~exist('nChains','var'),fprintf('Warning: nChains not specified.\n'),isError=1;end
if ~isfield(FileName,'MCMCDraws'),fprintf('Warning: FileName.MCMCDraws not specified.\n'),isError=1;end
if isError, disp('Warning: Errors found. Cannot Continue.'),return,end

%% Settings
if ~exist('nForecast','var'),nForecast=4; end
if ~exist('tStart','var'),tStart=1; end
if ~exist('BurnIn','var'),BurnIn=0.5; end
if ~exist('nDrawsStates','var'),nDrawsStates=1000; end
if ~exist('ShowFig','var'),ShowFig=0; end
if ~exist('Bands2Show','var'),Bands2Show=[50,60,70,80,90]; end
if ~exist('useForecastLabels','var')
    useForecastLabels = (exist('DateLabels','var') && isfield(DateLabels,'XTick')...
        && isfield(DateLabels,'XTickLabels'));
end
if useForecastLabels && ~exist('ForecastLabels','var'),ForecastLabels=DateLabels;end
isNewData = exist('DataNew','var');
if ~exist('PlotDir','var') || ~isfield(PlotDir,'Forecast')
  PlotDir.Forecast = 'Plots/Forecast/';
end
if ~isdir(PlotDir.Forecast),mkdir(PlotDir.Forecast),end
FileName.PlotsForecast = sprintf('%sPlotsForecastMCMCUpdate%.0f',FileName.Output,nUpdate);


%% ------------------------------------------------------------------------

%% Display
fprintf('\n******************')
fprintf('\n* Forecast Plots *')
fprintf('\n******************\n')

%% Set Timer
TimeStr = strrep(FileName.PlotsForecast,FileName.Output,'');
TimeElapsed.(TimeStr) = toc;

%% load the mcmc draws
fprintf('\nLoading data...\n')
xd = zeros(np,nDrawsStates/nChains,nChains);
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
StateVard = zeros(nStateVar,T+nForecast,nDrawsStates);
ObsVard = zeros(nObsVar,T+nForecast,nDrawsStates);
FileNamePost = FileName.Post;
parfor j=1:nDrawsUsed
    [ObsVard(:,:,j),StateVard(:,:,j)] = MakePlotsForecastFcn(...
        FileNamePost,xd(:,j),nStateVar,T,nForecast,Mats.KF,nObsVar,Params);
end

%% Draw plots
for j=1:2
    for jF=1:((j==1)*nObsVar+(j==2)*nStateVar)
        if ShowFig
            figure
        else
            figure('Visible','off')
        end
        if j==1
            var2show = ObsVar{jF};
            vcPlotDistBands(squeeze(ObsVard(jF,nPreSample+tStart:end,:))','Bands2Show',Bands2Show);
        else
            var2show = StateVar{jF};
            vcPlotDistBands(squeeze(StateVard(jF,nPreSample+tStart:end,:))','Bands2Show',Bands2Show);
        end
        if j==1 && isNewData
            hold on
            plot(DataNew(nPreSample+1:end,jF),'-r','LineWidth',2)
            hold off
        end
        title(var2show)
        if useForecastLabels
            set(gca,...
                'XTick',ForecastLabels.XTick,'XTickLabel',ForecastLabels.XTickLabels,'FontSize',8)
        end
        vcPrintPDF(sprintf('%s%sFig_%s',PlotDir.Forecast,FileName.PlotsForecast,var2show))
    end
end

%% Clean up
clear xd StateVard ObsVard nBands Bands2Show nPlotsStates Band BandPath BandColor StateVarjj BandPath var2show

%% close figures
if ~ShowFig, close all, end

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------
