%% MakeIRF
% 
% Master file for making Impulse Response Functions.
%
% Settings:
% 
%   UseDist
%   String variable indicating whether to use draw(s) from the prior or
%   posterior distribution. 
%   Options: 'priordraws', 'postdraws', or any of the fields in Params
%   Default: 'postdraws'
% 
%   nDrawsUsed
%   Number of draws to plot for priordraws or postdraws.
%   Default: 1000
% 
%   Shocks2Show
%   Cell array containing shocks to plot.
%   Choose shock names from 'ShockVar'.
%   Default: ShockVar
%
%   Vars2Show
%   Cell array containing cells of states and/or observables to plot.
%   Choose variables from 'StateVar' and 'ObsVar'.
%   Default: ObsVar
%
%   Shape2Plot
%   Cell array of 2 x 1 matrices indicating figure panel dimensions
%   for each of the Vars2Plot cells.
%   Default: 3 x 3 panel of figures.
% 
%   nSteps
%   Number of steps used in IRFs.
%   Default: 25
%
%   TickStepIRF
%   X-axis tick step for figures.
%   Default: 4
%
%   nThinning
%   Thinning desired. Plots will use every nThinning-th draw.
%   Default: 1
%
%   BurnIn
%   Percentage of draws to burn.
%   Default: 0.25
%
%   ShowFig
%   If set to 1 then the figures are shown in the display other wise
%   figures will be printed to eps files only.
%   Default: 0
%
%   yPerSlack
%   Percentage slack on y-axis tightening.
%   Default: 0.05
%
%   yMaxSlack
%   Absolute slack about y-axis.
%   Default: []
% 
% 
% See also:
% MakeIRFFcn,MakeMats
%
% .............................................................................
% 
% Created: June 28, 2010 by Vasco Curdia and Ging Cee Ng
% Updated: June 9, 2014 by Vasco Curdia
% 
% Copyright 2008-2014 by Vasco Curdia

%-----------------------------------------------------------------------------%

%% Settings

if ~exist('UseDist','var'), UseDist = 'postdraws'; end
if ~exist('nDrawsUsed','var'), nDrawsUsed = 1000; end
if ~exist('nThinning','var'), nThinning = 1; end
if ~exist('BurnIn','var'), BurnIn = 0.25; end
if ~exist('Bands2Show','var'),Bands2Show=[50,60,70,80,90]; end

if ~exist('Vars2Show','var'), 
    Vars2Show = {ObsVar}; 
    Panels = {'ObsVar'};
end
if ~exist('Panels','var')
    for j=1:length(Vars2Show)
        Panels{j} = int2str(j);
    end
end
if ~exist('Vars2ShowPretty','var'), Vars2ShowPretty = Vars2Show; end
if ~exist('Shocks2Show','var'), Shocks2Show = ShockVar; end
if ~exist('nSteps','var'), nSteps = 25; end
if ~exist('ShowFig','var'), ShowFig = 0; end
if ~exist('TickStepIRF','var'), TickStepIRF = 4; end
if ~exist('yPerSlack','var'), yPerSlack=0.05; end
if ~exist('yMaxSlack','var'), yMaxSlack=[]; end
if ~exist('yMinScale','var'), yMinScale=0; end

if ~exist('KeepData','var'), KeepData=0; end

if ~exist('PlotDir','var') || ~isfield(PlotDir,'IRF')
    PlotDir.IRF = 'Plots/IRF/';
end
if ~isdir(PlotDir.IRF), mkdir(PlotDir.IRF), end
if ~isfield(FileName,'PlotsIRF'), 
    FileName.PlotsIRF = sprintf('%sIRF',FileName.Output); 
end

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
if ~exist('Vars2ShowScale','var')
  for j=1:nPanels
    Vars2ShowScale{j} = ones(1,length(Vars2Show{j}));
  end
end

%% ------------------------------------------------------------------------

%% Display
fprintf('\n*******')
fprintf('\n* IRF *')
fprintf('\n*******\n')

%% Set Timer
TimeStr = 'MakeIRF';
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
nShocks2Show = length(Shocks2Show);
idxShocks = zeros(nShocks2Show,1);
for ii = 1:nShocks2Show
    idxShocks(ii) = find(ismember(ShockVar,Shocks2Show(ii)));
end    
nTempVars = length(TempVars);
idxVars = zeros(nTempVars,1);
AllVars = {StateVar{:},ObsVar{:}};
for ii = 1:nTempVars
    idxVars(ii) = find(ismember(AllVars,TempVars(ii)));
end

%% Generate IRF
FileNameMats = FileName.Mats;
IRF = zeros(nTempVars,nSteps,nShocks2Show,nDrawsUsed);
fprintf('Generating IRFs...\n');
parfor j=1:nDrawsUsed
  IRF(:,:,:,j) = MakeIRFFcn(xd(:,j),nSteps,idxShocks,idxVars,FileNameMats,...
    isObsMats,nObsVar);
end

%% Reassign IRF to Vars2Show
PlotIRF = cell(1,nPanels);
for ii = 1:nPanels
    [tf,idx] = ismember(Vars2Show{ii},TempVars);
    PlotIRF{ii} = IRF(idx,:,:,:);
end

%% IRF plot prep
XTicks = 0:TickStepIRF:(nSteps-1);
for ii = 1:length(XTicks)
    XTickLabels{ii} = sprintf('%i',XTicks(ii));
end
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

%% Plot IRFs
fprintf('Plotting IRFs...\n');
tid = 0:(nSteps-1);
for js=1:length(Shocks2Show)
    for ii = 1:nPanels
        if ShowFig 
            figure('Name',['Responses to innovation in ',Shocks2Show{js}])
        else
            figure('Visible','off')
        end
        for jj=1:length(Vars2Show{ii})
            subplot(Shape2Plot{ii}(1),Shape2Plot{ii}(2),jj);
            if nDrawsUsed==1
                vcPlot(tid,Vars2ShowScale{ii}(jj)*PlotIRF{ii}(jj,:,js))
            else
                vcPlotDistBands(tid,squeeze(...
                    Vars2ShowScale{ii}(jj)*PlotIRF{ii}(jj,:,js,:))',...
                                'Bands2Show',Bands2Show);
            end
            hold on
            plot(tid,zeros(1,nSteps),':k')
            hold off
            title(Vars2ShowPretty{ii}{jj})
            set(gca,'Xtick',XTicks,'XTickLabel',XTickLabels);
            axis tight
            yl = get(gca,'YLim');
            ySlack = max([yPerSlack*(yl(2)-yl(1)),yMaxSlack]);
            set(gca,'YLim',[min(yl(1)-ySlack,-yMinScale),...
                            max(yMinScale,yl(2)+ySlack)])
        end
        vcPrintPDF(sprintf('%s%s_%s_%s',...
                           PlotDir.IRF,FileName.PlotsIRF,...
                           Panels{ii},Shocks2Show{js}))
    end
end

%% Clean up
if ~KeepData
  clear xd postd irf IRF TempVars XTicks XTickLabels 
  clear PlotIRF yl ySlack
end

%% close figures
if ~ShowFig, close all, end

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------
