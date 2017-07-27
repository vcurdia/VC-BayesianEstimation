% MakePlotsHistDecomp
%
% Plots the historical decomposition of state and observation variables.
%
% NOTE: assumes that innovations show only once in the original system.
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
%   nDrawsHistDecomp
%   Number of draws to use.
%   Default: 1000
%
%   FigShapeStates (cell array)
%   shape of each figure. 1st element is the number of plots per row and
%   the 2nd is the number of plots per column.
%   Default: {1,1}
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
% Created: April 29, 2012 by Vasco Curdia
% Updated: May 10, 2012 by Vasco Curdia
%          Allow the use of MCMCDrawsRedux
% 
% Copyright 2012 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Settings
if ~exist('Vars2Show','var'),Vars2Show = ObsVar;end
if ~exist('ShockGroups','var')
  for jS=1:nShockVar
    ShockGroups{jS} = {ShockVar{jS}};
  end
end
if ~exist('ShockGroupsLabels','var')
  for jS=1:length(ShockGroups)
    ShockGroupsLabels{jS} = [ShockGroups{jS}{:}];
  end
end
if ~exist('ShockGroupsPrettyLabels','var')
  for jS=1:length(ShockGroups)
    ShockGroupsPrettyLabels{jS} = [ShockGroups{jS}{:}];
  end
end
if ~exist('ShockGroupsColors','var')
  ShockGroupsColors = vcColorScheme('nColors',length(ShockGroupsLabels));
end
if ~exist('FigShapeHistDecomp','var'),FigShapeHistDecomp={1,1};end
if ~exist('ShowFig','var'),ShowFig=0; end
if ~exist('ShowDates','var'),ShowDates = DateLabels;end
% useDateLabels = (exist('DateLabels','var') && isfield(DateLabels,'XTick')...
%     && isfield(DateLabels,'XTickLabels'));
if ~exist('nPreSample','var'), nPreSample = 0; end

if ~exist('UseMode','var'),UseMode=0; end
if ~exist('UseMedian','var'),UseMedian=0; end
if UseMode && UseMedian, error('Cannot use both median and mode.'), end
if ~exist('isDrawStates','var'),isDrawStates = UseMode || UseMedian; end
if ~(UseMode || UseMedian)
    if ~exist('nDrawsHistDecomp','var'),nDrawsHistDecomp=1000; end
    if ~exist('BurnIn','var'),BurnIn = 0.25; end
    if ~exist('nThinning','var'),nThinning = 1; end
end
if ~exist('UseMCMCDrawsRedux','var'),UseMCMCDrawsRedux=0; end

if ~exist('PlotDir','var') || ~isfield(PlotDir,'HistDecomp')
  PlotDir.HistDecomp = 'Plots/HistDecomp/';
end
if~isdir(PlotDir.HistDecomp),mkdir(PlotDir.HistDecomp),end
if ~exist('nUpdate','var')
  FileName.PlotsHistDecomp = sprintf('%sPlotsHistDecompMCMC',FileName.Output);
else
  FileName.PlotsHistDecomp = sprintf('%sPlotsHistDecompMCMCUpdate%.0f',FileName.Output,nUpdate);
end
if ~exist('LegPos','var'),LegPos = 'EO';end
if ~exist('LegOrientation','var'),LegOrientation = 'vertical';end

%% Check required settings
if ~(UseMode || UseMedian || UseMCMCDrawsRedux)
    isError = 0;
    if ~exist('nChains','var'),fprintf('Warning: nChains not specified.\n'),isError=1;end
    if ~isfield(FileName,'MCMCDraws'),fprintf('Warning: FileName.MCMCDraws not specified.\n'),isError=1;end
    if isError, disp('Warning: Errors found. Cannot Continue.'),return,end
end

%% ------------------------------------------------------------------------

%% Display
fprintf('\n***********************************')
fprintf('\n* Historical Decpomposition Plots *')
fprintf('\n***********************************\n')

%% Set Timer
TimeStr = strrep(FileName.PlotsHistDecomp,FileName.Output,'');
TimeElapsed.(TimeStr) = toc;

%% load draws
if UseMode
  fprintf('Using mode.\n')
  xd = [Params(:).postmode]';
elseif UseMedian
  fprintf('Using median.\n')
  xd = [Params(:).postp500]';
else
  fprintf('\nLoading data...\n')
  if UseMCMCDrawsRedux
    load(FileName.MCMCDrawsRedux)
    nThinningStates = nThinningRedux;
    clear postd
  else
    nThinningStates = max(floor(nDraws*nChains/nDrawsHistDecomp*(1-BurnIn)),1);
    for jChain=1:nChains
      load([FileName.MCMCDraws,int2str(jChain)],'xDraws','postDraws','nDraws')
      xd(:,:,jChain) = xDraws(:,1:nThinningStates:end);
    end
    clear xDraws postDraws
    nDrawsUsed = (1-BurnIn)*nDraws/nThinning*nChains;
    xd = reshape(xd(:,BurnIn*nDraws/nThinning+1:end,:),np,nDrawsUsed);
  end
  fprintf('Total number of draws per chain: %.0f\n', nDraws)
  fprintf('Thinning interval: %.0f\n', nThinningStates)
  fprintf('Burn in: %.0f%%\n', 100*BurnIn)
  fprintf('Total number of draws used: %.0f\n', nDrawsUsed)
end

%% Prepare indices
clear idxHD
idxHD.G2 = any(SymMats.StateEq.Gamma2~=0,2);
idxHD.G2 = any(SymMats.StateEq.Gamma1(idxHD.G2,:)~=0,1);
nShockGroups = length(ShockGroups);
idxHD.Shocks = zeros(nShockGroups,nShockVar);
DecompList = [ShockVar;{'Residual'}];
for jG=1:nShockGroups
  idxHD.Shocks(jG,:) = find(ismember(DecompList,ShockGroups{jG}));
end
nVars2Show = length(Vars2Show);
idxHD.StateVar = zeros(nVars2Show,1);
idxHD.ObsVar = zeros(nVars2Show,1);
for jV=1:nVars2Show
  [tf,idx] = ismember(Vars2Show{jV},StateVar);
  idxHD.StateVar(jV,:) = idx;
  [tf,idx] = ismember(Vars2Show{jV},ObsVar);
  idxHD.ObsVar(jV,:) = idx;
end


%% draw states
nDrawsUsed = size(xd,2);
HistDecompd = zeros(nVars2Show,T,nShockVar+1,nDrawsUsed);
FileNamePost = FileName.Post;
parfor j=1:nDrawsUsed
  HistDecompd(:,:,:,j) = MakePlotsHistDecompFcn(FileNamePost,xd(:,j),...
    nStateVar,nShockVar,nObsVar,nVars2Show,idxHD,T,j,isDrawStates);
end

%% Draw plots
nPlotsHistDecomp = FigShapeHistDecomp{1}*FigShapeHistDecomp{2};
nFigHistDecomp = ceil(nVars2Show/nPlotsHistDecomp);
TSample = T-nPreSample;
tid = find(ismember(TimeIdx,TimeIdxCreate(ShowDates.Start,ShowDates.End)));
TShow = length(tid);
ShowDates.XTick = find(ismember(TimeIdx,ShowDates.XTickLabels));
for jF=1:nFigHistDecomp
    if ShowFig
        figure
    else
        figure('Visible','off')
    end
    for jf=1:nPlotsHistDecomp 
      jVar = (jF-1)*nPlotsHistDecomp+jf;
      if jVar>nVars2Show, break, end
      subplot(FigShapeHistDecomp{:},jf)
      HistDecompj = squeeze(HistDecompd(jVar,nPreSample+tid,:,:));
      PlotData = zeros(TShow,nShockGroups,size(HistDecompj,3));
      for jG=1:nShockGroups
        PlotData(:,jG,:) = squeeze(sum(HistDecompj(:,idxHD.Shocks(jG,:),:),2));
      end
      PlotData(:,end+1,:) = squeeze(sum(HistDecompj,2));
      if ~(UseMode || UseMedian)
        PlotData = median(PlotData,3);
      end
%         PlotData = PlotData;
      bar(tid,PlotData(:,1:end-1),'BarWidth',1.5)
      colormap(ShockGroupsColors)
      hold on
      plot(tid,PlotData(:,end),'-k','LineWidth',2)
      hold off
      title(Vars2Show{jVar})
%         xlim(tid([1,end])+[-1,+1])
      set(gca,...
          'XTick',ShowDates.XTick,'XTickLabel',ShowDates.XTickLabels,...
          'FontSize',8)
%       legend(ShockGroupsLabels,'Location','SO','Orientation','horizontal')
      hleg = legend(ShockGroupsPrettyLabels,'Location',LegPos,'Orientation',LegOrientation);
      if strcmp(LegPos,'SO')
        legPos = get(hleg,'Position');
        legPos(2) = 0;
        set(hleg,'Position',legPos)
      end
    end
    if FigShapeHistDecomp{1}==1 && FigShapeHistDecomp{2}==1
        FigTxt = ['_',Vars2Show{jVar},'_',ShockGroupsLabels{:}];
    else
        FigTxt = sprintf('_%.0f',jF);
    end
    vcPrintPDF(sprintf('%s%sFig%s',PlotDir.HistDecomp,FileName.PlotsHistDecomp,FigTxt))
    clear FigTxt
end

%% Clean up
clear xd StateVard nBands Bands2Show nStatesShow nPlotsStates Band BandPath 
clear StateVarjj BandPath BandColor
clear idxHD ShockGroupsColors ShockGroupsLabels ShockGroupsPrettyLabels 
clear DecompList

%% close figures
if ~ShowFig, close all, end

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------
