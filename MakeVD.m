% MakeVD
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
% Updated: June 9, 2014 by Vasco Curdia
% 
% Copyright 2010-2014 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Settings

if ~exist('UseDist','var'), UseDist = 'postdraws'; end
if ~exist('nDrawsUsed','var'), nDrawsUsed = 1000; end
if ~exist('nThinning','var'), nThinning = 1; end
if ~exist('BurnIn','var'), BurnIn = 0.25; end
if ~exist('Bands2Show','var'),Bands2Show=[50,60,70,80,90]; end

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
%   ShockColors = vcColorScheme('nColors',nShockVar);
    ShockColors = vcColorScheme('nColors',nShockVar,...
                                'UseColors','Excel2010',...
                                'LightFactors',[0,0.4]);
end

if ~exist('KeepResultMats','var'),KeepResultMats = 0; end
if ~exist('isSilent','var'),isSilent = 0; end

if ~exist('PlotDir','var') || ~isfield(PlotDir,'VD')
  PlotDir.VD = 'Plots/VD/';
end
if~isdir(PlotDir.VD),mkdir(PlotDir.VD),end
FileName.PlotsVD = sprintf('%sVD',FileName.Output);

if ~exist('ShowPlotTitle','var'),ShowPlotTitle=1;end
if ~exist('LegTxt','var'),LegTxt = ShockVar;end
if ~exist('LegPos','var'),LegPos = 'EO';end
if ~exist('LegOrientation','var'),LegOrientation = 'vertical';end

%% Auxilliary Variables
nPanels = length(Panels);
if nPanels>1
    TempVars = union(Vars2Show{:});
else
    TempVars = Vars2Show{1};
end
if any(ismember(TempVars,ObsVar))
    isObsMats=2;
else
    isObsMats=0;
end


%% ------------------------------------------------------------------------

%% Display
fprintf('\n*******************************')
fprintf('\n* MCMC Variance Decomposition *')
fprintf('\n*******************************\n')

%% Set Timer
TimeStr = 'VD';
TimeElapsed.(TimeStr) = toc;

%% Prepare Draws
if strcmp(UseDist,'priordraws')
    if ~isfield(FileName,'GenPriorDraw')
        MakeGenPriorDraw
    end
    xd = feval(FileName.GenPriorDraw,nDrawsUsed);
elseif strcmp(UseDist,'postdraws')
    if ~isfield(FileName,'MCMCDrawsRedux')
        MakeMCMCDrawsRedux
    end
    load(FileName.MCMCDrawsRedux,'xd')
    xd = xd(:,1:nDrawsUsed);
else
    if ~isfield(Params,UseDist)
        fprintf(2,'Did not recognize distribution to use. Cannot proceed.\n');
        return
    end
    xd = [Params(:).(UseDist)]';
end
nDrawsUsed = size(xd,2);

%% Create Indices
nTempVars = length(TempVars);
idxVars = zeros(nTempVars,1);
AllVars = {StateVar{:},ObsVar{:}};
for ii = 1:nTempVars
    idxVars(ii) = find(ismember(AllVars,TempVars(ii)));
end

%% Generate VD
isInfHorizon = ismember(inf,VDHorizons);
VDHorizons = sort(VDHorizons);
nHorizons = length(VDHorizons);
MaxHorizon = VDHorizons(end-isInfHorizon);
idxMat = eye(nShockVar);
VD = zeros(nTempVars,nShockVar,nHorizons,nDrawsUsed);
FileNameMats = FileName.Mats;
parfor jd=1:nDrawsUsed
  VD(:,:,:,jd) = MakeVDFcn(xd(:,jd),nShockVar,idxVars,FileNameMats,isObsMats,...
    nHorizons,idxMat,MaxHorizon,VDHorizons,isInfHorizon,isSilent);
end

%% Reassign VD to Vars2Show
ShowVD = cell(1,nPanels);
for ii = 1:nPanels
    [tf,idx] = ismember(Vars2Show{ii},TempVars);
    ShowVD{ii} = VD(idx,:,:,:);
end

%% Create tables
fprintf('\nVariance decomposition:')
fprintf('\n=======================\n')
fprintf('Using %s',UseDist)
if nDrawsUsed>1
    fprintf(' (Percentile: %.1f)',VDPrctile)
end
fprintf('\n')
for jPanel=1:nPanels
    fprintf('\nPanel: %s\n',Panels{jPanel})
    Vars2Showj = Vars2Show{jPanel};
    nVars2Showj = length(Vars2Showj);
    ShowVDj = ShowVD{jPanel};
    znamelength = [cellfun('length',Vars2Showj)];
    znamelengthmax = max(znamelength);
    snamelength = [cellfun('length',ShockVar)];
    snamelengthmax = max(max(snamelength),5);
    for jh=xTickIdx
        if nDrawsUsed==1
            VDj = ShowVDj(:,:,jh);
        else
            VDj = prctile(ShowVDj(:,:,jh,:),VDPrctile,4);
        end
        fprintf('\nHorizon: %.0f\n',VDHorizons(jh))
        fprintf(['%-',int2str(znamelengthmax),'s',...
                 repmat(['   %',int2str(snamelengthmax),'s'],1,nShockVar),...
                 '\n'],'',ShockVar{:})
        for jz=1:nVars2Showj
            fprintf(['%-',int2str(znamelengthmax),'s',...
                     repmat(['   %',int2str(snamelengthmax),'.3f'],1,...
                            nShockVar),'\n'],Vars2Showj{jz},VDj(jz,:))
        end
    end
end

%% plot prep
if exist('Shape2Plot','var')
    if length(Shape2Plot)~=length(Vars2Show)
        error(...
            'Shape2Plot dimensions must be specified for each Vars2Show cell.');
    end
elseif ~exist('Shape2Plot','var')
    for ii = 1:nPanels
        Shape2Plot{ii} = [1,1];
        jDim = 1;
        while prod(Shape2Plot{ii})<length(Vars2Show{ii})
            Shape2Plot{ii}(jDim) = Shape2Plot{ii}(jDim)+1;
            jDim = ~(jDim-1)+1;
        end
    end
end

%% Create plots
tid = 1:nHorizons;
if ~exist('xTickIdx','var'),xTickIdx=tid;end
if ~exist('xTickLabels','var'),xTickLabels=VDHorizons;end
for ii=1:nPanels
    if ShowFig
        figure
    else
        figure('Visible','off')
    end
    ShowVDj = ShowVD{ii};
    for jj=1:length(Vars2Show{ii})
        hsubp(jj) = subplot(Shape2Plot{ii}(1),Shape2Plot{ii}(2),jj);
        if nDrawsUsed==1
            PlotData = squeeze(ShowVDj(jj,:,:));
        else
            PlotData = squeeze(prctile(ShowVDj(jj,:,:,:),VDPrctile,4));
        end
        bar(tid,PlotData','stacked','BarWidth',1,'EdgeColor','none')
        axis tight
        colormap(ShockColors)
        if ShowPlotTitle
            title(Vars2ShowPretty{ii}{jj})
        end
        ylim([0,1])
        set(gca,'XTick',xTickIdx,'XTickLabel',xTickLabels,'FontSize',8)
    end
    if all(Shape2Plot{ii}==1)
        hleg = legend(LegTxt,'Location',LegPos,'Orientation',LegOrientation);
        if strcmp(LegPos,'SO')
            legPos = get(hleg,'Position');
            legPos(2) = 0;
            set(hleg,'Position',legPos)
        end
    else
        hleg = legend(LegTxt,'Location',LegPos,'Orientation','horizontal');
        legPos = get(hleg,'Position');
        xL = get(hsubp((Shape2Plot{ii}(1)-1)*Shape2Plot{ii}(2)+1),'Position');
        xR = get(hsubp((Shape2Plot{ii}(1)-1)*Shape2Plot{ii}(2)),'Position');
        legPos(1) = xL(1)+(xR(1)-xL(1))/2+(xL(3)-legPos(3))/2;
        legPos(2) = 0;
        set(hleg,'Position',legPos)
    end
    vcPrintPDF(sprintf('%s%s_%s',PlotDir.VD,FileName.PlotsVD,Panels{ii}))
end


%% Clean up
clear xd postd Vs ShowVDj
if ~KeepResultMats
  clear VD ShowVD
end
% for jp=1:np, eval(sprintf('syms %s',Params(jp).name)), end

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------

