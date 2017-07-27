% MakePlotsMCMCTrace
%
% Plots draws and histograms to help checking for convergence
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
%   nThinning
%   Thinning desired. Plots will use every nThinning-th draw.
%   Default: 1
%
%   nTraceStep
%   Number of draws after which the trace is updated. Means and stdev will
%   be updated every nTraceStep-th draw.
%   Default: nDraws/100
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
% Created: April 6, 2008 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: September 26, 2011 by Vasco Curdia
%          Creates plots in folder.
% 
% Copyright 2008-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Check required settings
isError = 0;
if ~exist('nChains','var'),fprintf('Warning: nChains not specified.\n'),isError=1;end
if ~isfield(FileName,'MCMCDraws'),fprintf('Warning: FileName.MCMCDraws not specified.\n'),isError=1;end
if isError, disp('Warning: Errors found. Cannot Continue.'),return,end

%% Settings
if ~exist('BurnIn','var'),BurnIn=0.5; end
if ~exist('nThinning','var'),nThinning=1; end
if ~exist('nTraceStep','var'),nTraceStep=floor(max(nDraws/nThinning/100,1)); end
if ~exist('ShowFig','var'),ShowFig=0; end
if ~exist('PlotDir','var') || ~isfield(PlotDir,'MCMCTrace')
  PlotDir.MCMCTrace = 'Plots/MCMCTrace/';
end
if~isdir(PlotDir.MCMCTrace),mkdir(PlotDir.MCMCTrace),end
FileName.PlotsMCMCTrace = sprintf('%sPlotsMCMCTraceUpdate%.0f',FileName.Output,nUpdate);

%% ------------------------------------------------------------------------

%% Display
fprintf('\n********************')
fprintf('\n* MCMC Trace Plots *')
fprintf('\n********************\n')

%% Set Timer
TimeStr = strrep(FileName.PlotsMCMCTrace,FileName.Output,'');
TimeElapsed.(TimeStr) = toc;

%% load the mcmc draws
fprintf('\nLoading data...\n')
for jChain=1:nChains
    load([FileName.MCMCDraws,int2str(jChain)],'xDraws','postDraws','nDraws')
    xd(:,:,jChain) = xDraws(:,1:nThinning:end);
    postd(:,:,jChain) = postDraws(:,1:nThinning:end);
end
clear xDraws postDraws
nDrawsUsed = (1-BurnIn)*nDraws/nThinning;
xd = xd(:,BurnIn*nDraws/nThinning+1:end,:);
fprintf('Total number of draws per chain: %.0f\n',nDraws)
fprintf('Thinning interval: %.0f\n',nThinning)
fprintf('Burn in: %.0f%%\n\n',100*BurnIn)
fprintf('Trace step: %.0f\n',nTraceStep)

%% create trace plots
for jp=1:np
    Sample = squeeze(xd(jp,:,:));
    SampleID = [nTraceStep:nTraceStep:nDrawsUsed]';
    nSample = length(SampleID);
    clear RollingMean RollingSD
    for js=1:nSample;
        RollingMean(js,:) = mean(Sample(1:SampleID(js),:),1);
        RollingSD(js,:) = std(Sample(1:SampleID(js),:),0,1);
    end
    Mean = mean(Sample(:));
    SD = std(Sample(:));
    MeanBounds = [min(min(RollingMean(:)),Mean-2*SD) max(max(RollingMean(:)),Mean+2*SD)];
    MeanBounds = MeanBounds + (-1).^(1:-1:0)*0.01*(MeanBounds(2)-MeanBounds(1));
    SDBounds = [0 1.01*max(max(RollingSD(:)),2*SD)];
    if ShowFig
        figure
    else
        figure('Visible','off')
    end
    for jj=1:nChains
        subplot(nChains,2,(jj-1)*2+1)
        plot(SampleID,Mean*ones(size(SampleID)),'-','Color',[0,.5,0])
        hold on
        plot(SampleID,(Mean-SD)*ones(size(SampleID)),':r',...
             SampleID,(Mean+SD)*ones(size(SampleID)),':r')
        plot(SampleID,RollingMean(:,jj),'-b','LineWidth',2)
    title(['Parameter ',Params(jp).prettyname,' in each chain'])
        if jj==1, title(['Rolling mean of ',Params(jp).prettyname]), end
        ylim(MeanBounds)
        xlim(SampleID([1,end]))
        set(gca,'FontSize',8)
        subplot(nChains,2,jj*2)
        plot(SampleID,SD*ones(size(SampleID)),'-','Color',[0,.5,0])
        hold on
        plot(SampleID,RollingSD(:,jj),'-b','LineWidth',2)
        if jj==1, title(['Rolling std dev of ',Params(jp).prettyname]), end
        ylim(SDBounds)
        xlim(SampleID([1,end]))
        set(gca,'FontSize',8)
    end
    vcPrintPDF(sprintf('%s%s_%s',PlotDir.MCMCTrace,FileName.PlotsMCMCTrace,Params(jp).name))
end
clear xd Mean SD MeanBounds SDBounds Sample SampleID nSample RollingMean postd

%% close figures
if ~ShowFig, close all, end

%% Elapsed time
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------
