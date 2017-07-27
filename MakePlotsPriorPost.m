% MakePlotsPriorPost
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
%   nBin
%   number of bins for the histogram
%   Default: 50
%
%   FigShape (cell array)
%   shape of each figure. 1st element is the number of plots per row and
%   the 2nd is the number of plots per column.
%   Default: {3,3}
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
if ~exist('nBin','var'),nBin=50;end
if ~exist('FigShape','var'),FigShape={3,3};end
if ~exist('ShowFig','var'),ShowFig=0; end
if ~exist('PlotDir','var') || ~isfield(PlotDir,'PriorPost')
  PlotDir.PriorPost = 'Plots/PriorPost/';
end
if~isdir(PlotDir.PriorPost),mkdir(PlotDir.PriorPost),end
FileName.PlotsPriorPost = sprintf('%sPlotsPriorPostUpdate%.0f',FileName.Output,nUpdate);

%% ------------------------------------------------------------------------

%% Display
fprintf('\n*************************')
fprintf('\n* Prior-Posterior Plots *')
fprintf('\n*************************\n')

%% Set Timer
TimeStr = strrep(FileName.PlotsPriorPost,FileName.Output,'');
TimeElapsed.(TimeStr) = toc;

%% load the mcmc draws
fprintf('\nLoading data...\n')
for jChain=1:nChains
    load([FileName.MCMCDraws,int2str(jChain)],'xDraws','nDraws')
    xd(:,:,jChain) = xDraws(:,1:nThinning:end);
end
clear xDraws postDraws
nDrawsUsed = (1-BurnIn)*nDraws/nThinning*nChains;
xd = reshape(xd(:,BurnIn*nDraws/nThinning+1:end,:),np,nDrawsUsed);
fprintf('Total number of draws per chain: %.0f\n', nDraws)
fprintf('Thinning interval: %.0f\n', nThinning)
fprintf('Burn in: %.0f%%\n', 100*BurnIn)
fprintf('Total number of draws used: %.0f\n', nDrawsUsed)

%% Draw plots
nPlots = FigShape{1}*FigShape{2};
nFig = ceil(np/nPlots);
for jF=1:nFig
    if ShowFig
        figure
    else
        figure('Visible','off')
    end
    for jf=1:nPlots 
        jp = (jF-1)*nPlots+jf;
        if jp>np, break, end
        subplot(FigShape{:},jf)
        % posterior
        xPost = xd(jp,:);
        [xFreq,xOut] = hist(xPost,nBin);
        xStep = xOut(2)-xOut(1);
        xOutMin = min(xOut);
        xOutMax = max(xOut);
        % prior
        xCrit = [Params(jp).priorp025,Params(jp).priorp975];
        xPlot = [sort(xOutMin-xStep:-xStep:xCrit(1)) xOut xOutMax+xStep:xStep:xCrit(2)];
        nPlot = zeros(size(xPlot));
        nIdx = ismember(xPlot,xOut);
        nPlot(nIdx) = xFreq;
        for jx=1:length(xPlot)
            xPriorPdf(jx) = eval(sprintf(Params(jp).priorpdfcmd,xPlot(jx)));
        end
%         nPlot = nPlot*(max(xPriorPdf)-min(xPriorPdf))/(max(nPlot)-min(nPlot));
        % normalize hist to have unit area
        nPlot = nPlot/sum(nPlot*xStep);
        % plot
        bar(xPlot,nPlot)
        hold on
        plot(xPlot,xPriorPdf,'r','LineWidth',2)
        hold off
        title(Params(jp).prettyname)
        xBounds = xPlot([1,end]);
        xBounds = xBounds+(-1).^(1:-1:0)*0.01*(xBounds(2)-xBounds(1));
        xlim(xBounds)
        xTickLabels = xBounds(1):(xBounds(2)-xBounds(1))/8:xBounds(2);
        yBounds = max(max(nPlot),max(xPriorPdf));
        yBounds = [0,yBounds+0.01*yBounds];
        ylim(yBounds)
        set(gca,'YTick',[],'XTick',xTickLabels([2,5,8]),'FontSize',8)
        clear xPost xFreq xOut xStep xOutMin xOutMax xPlot nIdx xPriorPdf xBounds xTickLabels yBounds
    end
    vcPrintPDF(sprintf('%s%sFig%.0f',PlotDir.PriorPost,FileName.PlotsPriorPost,jF))
end

%% Clean up
clear xd nPlots nFig

%% close figures
if ~ShowFig, close all, end

%% Elapsed time
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------
