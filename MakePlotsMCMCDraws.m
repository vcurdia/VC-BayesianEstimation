% MakePlotsMCMCDraws
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
% See also:
% SetDSGE, GenSymVars, DataAnalysis, PriorAnalysis, GenPost, MaxPost,
% MaxPostFcn, MakeTableMaxPost, MCMC, MCMCFcn, MCMCSearchScaleFactor,
% MakePlotsMCMCDraws, MCMCInference, MakeTableMCMCInference, 
% MakePlotsMCMCTrace, MakePlotsMCMCPriorPost, MCMCConv, MakeTableMCMCConv,
% MakeReportMCMCPlots, MakeReportMCMC
%
% .........................................................................
%
% Created: April 1, 2008 by Vasco Curdia
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
if ~exist('ShowFig','var'),ShowFig=0; end
if ~exist('PlotDir','var') || ~isfield(PlotDir,'MCMCDraws')
  PlotDir.MCMCDraws = 'Plots/MCMCDraws/';
end
if~isdir(PlotDir.MCMCDraws),mkdir(PlotDir.MCMCDraws),end
FileName.PlotsMCMCDraws = sprintf('%sPlotsMCMCDrawsUpdate%.0f',FileName.Output,nUpdate);

%% ------------------------------------------------------------------------

%% Display
fprintf('\n********************')
fprintf('\n* MCMC Draws Plots *')
fprintf('\n********************\n')

%% Set Timer
TimeStr = strrep(FileName.PlotsMCMCDraws,FileName.Output,'');
TimeElapsed.(TimeStr) = toc;

%% load the mcmc draws
fprintf('\nLoading data...\n')
for jChain=1:nChains
    load([FileName.MCMCDraws,int2str(jChain)],'xDraws','postDraws','nDraws')
    xd(:,:,jChain) = xDraws(:,1:nThinning:end);
    postd(:,:,jChain) = postDraws(:,1:nThinning:end);
end
clear xDraws postDraws
nDrawsUsed = nDraws/nThinning;
fprintf('Total number of draws per chain: %.0f\n',nDraws)
fprintf('Thinning interval: %.0f\n',nThinning)
fprintf('Burn in: %.0f%%\n\n',100*BurnIn)

%% plot the posterior values
p_max = max(postd(:));
p_min = min(postd(:));
p_spread = p_max-p_min;
p_max = p_max+.01*p_spread;
p_min = p_min-.01*p_spread;
if ShowFig
    figure
else
    figure('Visible','off')
end
for jChain=1:nChains
    subplot(nChains,1,jChain)
    plot(postd(:,:,jChain))
    ylim([p_min p_max])
    set(gca,'FontSize',8)
end
subplot(nChains,1,1)
title('Posterior density (log) value in each chain')
vcPrintPDF(sprintf('%s%s_post',PlotDir.MCMCDraws,FileName.PlotsMCMCDraws))
clear postd

%% Plot scalar estimates
for jp=1:np
    p_max = max(max(squeeze(xd(jp,:,:))));
    p_min = min(min(squeeze(xd(jp,:,:))));
    p_spread = p_max-p_min;
    p_max = p_max+.01*p_spread;
    p_min = p_min-.01*p_spread;
    if ShowFig
        figure
    else
        figure('Visible','off')
    end
    for jChain=1:nChains
        subplot(nChains,2,(jChain-1)*2+1)
        plot(xd(jp,:,jChain))
        ylim([p_min p_max])
        set(gca,'XTick',BurnIn*nDrawsUsed,'XGrid','on','FontSize',8)
        subplot(nChains,2,jChain*2)
        hist(xd(jp,BurnIn*nDrawsUsed+1:end,jChain),100)
        xlim([p_min p_max])
        set(gca,'FontSize',8)
    end
    subplot(nChains,2,1)
    title(['Parameter ',Params(jp).prettyname,' in each chain'])
    subplot(nChains,2,2)
    title(['Hist excluding initial ',int2str(100*BurnIn),'% of obs'])
    vcPrintPDF(sprintf('%s%s_%s',PlotDir.MCMCDraws,FileName.PlotsMCMCDraws,Params(jp).name));
end
clear p_min p_max p_spread xd nDrawsUsed

%% close figures
if ~ShowFig, close all, end

%% Elapsed time
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------
