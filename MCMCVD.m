% MCMCVD
%
% Generates Variance decomposition tables.
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
%   VDHorizons
%   Horizons at which to compute the VD.
%   Default: [1:32,T,inf]]
%
%   BurnIn
%   Percentage of draws to burn.
%   Default: 0.25
%
%   nThinning
%   Thinning desired. Plots will use every nThinning-th draw.
%   Default: 1
%
%   UseMode
%   computes VD at mode
%   Default: 0
%
%   UseMedian
%   computes VD at median
%   Default: 1
%
% See also:
% SetDSGE, GenSymVars, DataAnalysis, PriorAnalysis, GenPost, MaxPost, 
% MakeTableMaxPost, MCMC, MCMCSearchScaleFactor, MakePlotsMCMCConv, 
% MCMCInference, MakeTableMCMCInference, MakePlotsMCMCTrace, 
% MakePlotsMCMCPriorPost, MCMCConv, MakeTableMCMCConv, MCMCVD, MCMCVDFcn
%
% .........................................................................
%
% Created: April 15, 2010 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: July 29, 2011 by Vasco Curdia
%          Change the way parallel works
% Updated: May 10, 2012 by Vasco Curdia
%          - Allow the use of MCMCDrawsRedux
%          - option to keep result matrices
%          - allow for variable selection
%          - show figure with VD
% 
% Copyright 2010-2012 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Settings
if ~exist('Vars2Show','var'),Vars2Show = ObsVar;end
if ~exist('Vars2ShowPretty','var'),Vars2ShowPretty = Vars2Show;end

if ~exist('VDHorizons','var')
  VDHorizons = [1:32,T,inf];
  xTickIdx = [1,4:4:32,34];
  xTickLabels = {1,4:4:32,'   inf'};
end
if ~exist('VDPrctile','var'),VDPrctile = 50; end
if ~exist('ShowFig','var'),ShowFig=0; end
if ~exist('FigShapeVD','var'),FigShapeVD={1,1};end
if ~exist('ShockColors','var')
  ShockColors = vcColorScheme('nColors',nShockVar);
end

if ~exist('UseMode','var'),UseMode=0; end
if ~exist('UseMedian','var'),UseMedian=1; end
if UseMode && UseMedian, error('Cannot use both median and mode.'), end
if ~(UseMode || UseMedian)
    if ~exist('BurnIn','var'),BurnIn = 0.25; end
    if ~exist('nThinning','var'),nThinning = 1; end
    if ~exist('Bands2Show','var'),Bands2Show=[50,60,70,80,90]; end
end
if ~exist('UseMCMCDrawsRedux','var'),UseMCMCDrawsRedux=0; end
if ~exist('KeepResultMats','var'),KeepResultMats = 0; end
if ~exist('isSilent','var'),isSilent = 0; end

if ~exist('PlotDir','var') || ~isfield(PlotDir,'VD')
  PlotDir.VD = 'Plots/VD/';
end
if~isdir(PlotDir.VD),mkdir(PlotDir.VD),end
if ~exist('nUpdate','var')
  FileName.PlotsVD = sprintf('%sPlotsVDMCMC',FileName.Output);
else
  FileName.PlotsVD = sprintf('%sPlotsVDMCMCUpdate%.0f',FileName.Output,nUpdate);
end
if ~exist('ShowPlotTitle','var'),ShowPlotTitle=1;end
if ~exist('LegTxt','var'),LegTxt = ShockVar;end
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
fprintf('\n*******************************')
fprintf('\n* MCMC Variance Decomposition *')
fprintf('\n*******************************\n')

%% Set Timer
TimeStr = 'VD';
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
    nThinning = nThinningRedux;
    clear postd
  else
    for jChain=1:nChains
      load([FileName.MCMCDraws,int2str(jChain)],'xDraws','postDraws','nDraws')
      xd(:,:,jChain) = xDraws(:,1:nThinning:end);
    end
    clear xDraws postDraws
    nDrawsUsed = (1-BurnIn)*nDraws/nThinning*nChains;
    xd = reshape(xd(:,BurnIn*nDraws/nThinning+1:end,:),np,nDrawsUsed);
  end
  fprintf('Total number of draws per chain: %.0f\n', nDraws)
  fprintf('Thinning interval: %.0f\n', nThinning)
  fprintf('Burn in: %.0f%%\n', 100*BurnIn)
  fprintf('Total number of draws used: %.0f\n', nDrawsUsed)
end

%% Generate VD
if any(ismember(Vars2Show,ObsVar))
  isObsMats=2;
else
  isObsMats=0;
end
nVars2Show = length(Vars2Show);
idxVars = zeros(nVars2Show,1);
for ii = 1:length(Vars2Show)
  idxVars(ii) = find(ismember([StateVar;ObsVar],Vars2Show(ii)));
end
isInfHorizon = ismember(inf,VDHorizons);
VDHorizons = sort(VDHorizons);
nHorizons = length(VDHorizons);
MaxHorizon = VDHorizons(end-isInfHorizon);
idxMat = eye(nShockVar);
nDrawsUsed = size(xd,2);
VD = zeros(nVars2Show,nShockVar,nHorizons,nDrawsUsed);
FileNameMats = FileName.Mats;
parfor jd=1:nDrawsUsed
  VD(:,:,:,jd) = MCMCVDFcn(xd(:,jd),nShockVar,idxVars,FileNameMats,isObsMats,...
    nHorizons,idxMat,MaxHorizon,VDHorizons,isInfHorizon,isSilent);
end

%% Prepare variables
nVars2Show = length(Vars2Show);
AllVars = [StateVar;ObsVar];

%% Create table
fprintf('\nVariance decomposition:')
fprintf('\n=======================\n')
if UseMode
    fprintf('(at the mode)\n')
elseif UseMedian
    fprintf('(at the median)\n')
else
    fprintf('(Percentile: %.1f)\n',VDPrctile)
end
znamelength = [cellfun('length',Vars2Show)];
znamelengthmax = max(znamelength);
snamelength = [cellfun('length',ShockVar)];
snamelengthmax = max(max(snamelength),5);
for jh=xTickIdx
    if UseMode || UseMedian
        VDj = VD(:,:,jh);
    else
        VDj = prctile(VD(:,:,jh,:),VDPrctile,4);
    end
    fprintf('\nHorizon: %.0f\n',VDHorizons(jh))
    fprintf(['%-',int2str(znamelengthmax),'s',...
        repmat(['   %',int2str(snamelengthmax),'s'],1,nShockVar),'\n'],'',ShockVar{:})
    for jz=1:nVars2Show
      fprintf(['%-',int2str(znamelengthmax),'s',...
          repmat(['   %',int2str(snamelengthmax),'.3f'],1,nShockVar),'\n'],Vars2Show{jz},VDj(jz,:))
    end
end

%% Create plots
nPlotsVD = FigShapeVD{1}*FigShapeVD{2};
nFigVD = ceil(nVars2Show/nPlotsVD);
tid = 1:nHorizons;
if ~exist('xTickIdx','var'),xTickIdx=tid;end
if ~exist('xTickLabels','var'),xTickLabels=VDHorizons;end
for jF=1:nFigVD
  if ShowFig
    figure
  else
    figure('Visible','off')
  end
  for jf=1:nPlotsVD
    jVar = (jF-1)*nPlotsVD+jf;
    if jVar>nVars2Show, break, end
    hsubp(jf) = subplot(FigShapeVD{:},jf);
    if UseMode || UseMedian
        PlotData = squeeze(VD(jVar,:,:));
    else
        PlotData = squeeze(prctile(VD(jVar,:,:,:),VDPrctile,4));
    end
    bar(tid,PlotData','stacked','BarWidth',1,'EdgeColor','none')
    axis tight
    colormap(ShockColors)
    if ShowPlotTitle
      title(Vars2ShowPretty{jVar})
    end
    ylim([0,1])
    set(gca,'XTick',xTickIdx,'XTickLabel',xTickLabels,'FontSize',8)
%     legend(ShockVar,'Location','SO','Orientation','horizontal')
    if all([FigShapeVD{:}]==1)
      hleg = legend(LegTxt,'Location',LegPos,'Orientation',LegOrientation);
      if strcmp(LegPos,'SO')
        legPos = get(hleg,'Position');
        legPos(2) = 0;
        set(hleg,'Position',legPos)
      end
    elseif jf==nPlotsVD || jVar==nVars2Show
      hleg = legend(LegTxt,'Location',LegPos,'Orientation','horizontal');
      legPos = get(hleg,'Position');
      xL = get(hsubp((FigShapeVD{1}-1)*FigShapeVD{2}+1),'Position');
      xR = get(hsubp((FigShapeVD{1}-1)*FigShapeVD{2}),'Position');
      legPos(1) = xL(1)+(xR(1)-xL(1))/2+(xL(3)-legPos(3))/2;
      legPos(2) = 0;
      set(hleg,'Position',legPos)
    end
  end
  if FigShapeVD{1}==1 && FigShapeVD{2}==1
    FigTxt = ['_',Vars2Show{jVar}];
  else
    FigTxt = sprintf('_%.0f',jF);
  end
  vcPrintPDF(sprintf('%s%sFig%s',PlotDir.VD,FileName.PlotsVD,FigTxt))
  clear FigTxt
end


%% Clean up
clear xd postd Vs
if ~KeepResultMats
  clear VD
end
% for jp=1:np, eval(sprintf('syms %s',Params(jp).name)), end

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------

