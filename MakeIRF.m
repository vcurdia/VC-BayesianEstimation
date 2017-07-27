%% MakeIRF
% 
% Master file for making Impulse Response Functions.
%
% Settings:
% 
%   UseDist
%   String variable indicating whether to use draw(s) from the prior or
%   posterior distribution. 
%   Options: 'Prior' or 'Post'
%   Default: 'Post'
% 
%   UseStat
%   String variable indicating which draw(s) to use for generating the
%   IRFs. 
%   Options: 'Median','Mode','Mean', or 'FanChart' (uses nDrawsIRF)
%   Default: 'Median'
%
%   nDrawsInf
%   Number of draws to plot for FanChart setting.
%   Default: 1000
% 
%   Shocks2Plot
%   Cell array containing shocks to plot.
%   Choose shock names from 'ShockVar'.
%   Default: ShockVar
%
%   Vars2Plot
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
% .........................................................................
% 
% Date Created: June 28, 2010 by Vasco Curdia and Ging Cee Ng
% Updated: July 26, 2011 by Vasco Curdia
% Updated: July 28, 2011 by Ging Cee Ng (added dotted zero line)
% Updated: July 29, 2011 by Vasco Curdia
%          Change the way parallel works
% 
% Copyright 2008-2011 by Vasco Curdia and Ging Cee Ng

%% ------------------------------------------------------------------------

%% Fill in missing options
if ~isfield(FileName,'PlotsIRF'), FileName.PlotsIRF = sprintf('%sIRF',FileName.Output); end
if ~exist('UseDist','var'), UseDist = 'Post'; end
if ~exist('UseStat','var'), UseStat = 'Median'; end
if ~exist('nDrawsInf','var'), nDrawsInf = 1000; end
if ~exist('nSteps','var'), nSteps = 25; end
if ~exist('nThinning','var'), nThinning = 1; end
if ~exist('BurnIn','var'), BurnIn = 0.25; end
if ~exist('Vars2Plot','var'), Vars2Plot = {ObsVar}; end
if ~exist('Vars2PlotPretty','var'), Vars2PlotPretty = Vars2Plot; end
if ~exist('Shocks2Plot','var'), Shocks2Plot = ShockVar; end
if ~exist('ShowFig','var'), ShowFig = 0; end
if ~exist('TickStepIRF','var'), TickStepIRF = 4; end
if ~exist('Bands2Show','var'),Bands2Show=[50,60,70,80,90]; end
if ~exist('yPerSlack','var'), yPerSlack=0.05; end
if ~exist('yMaxSlack','var'), yMaxSlack=[]; end
if ~exist('yMinScale','var'), yMinScale=0; end
if ~exist('KeepData','var'), KeepData=0; end

%% Auxilliary Variables
nP = length(Vars2Plot);
if nP>1
    TempVars = union(Vars2Plot{:});
else
    TempVars = Vars2Plot{1};
end
if any(ismember(TempVars,ObsVar))
    isObsMats=2;
else
    isObsMats=0;
end
if ~exist('Vars2PlotScale','var')
  for j=1:nP
    Vars2PlotScale{j} = ones(1,length(Vars2Plot{j}));
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

%% Identify Parameter Vector
if strcmp('Prior',UseDist), distParams = 'prior';
elseif strcmp('Post',UseDist), distParams = 'post';
else error('Must use either Prior or Post.')
end 

if strcmp('Median',UseStat), distStat = 'p500';
elseif strcmp('Mode',UseStat), distStat = 'mode';
elseif strcmp('Mean',UseStat), distStat = 'mean';
elseif strcmp('FanChart',UseStat), distStat = 'Fan';
else error('Must use either Median, Mode, Mean, or FanChart.');
end

PlotParam = strcat(sprintf('%s%s',distParams,distStat));
clear distParams distStat

%% Construct Parameter Draws Matrix
if strcmp('postFan',PlotParam)
    fprintf('\nFanChart mode: Loading MCMC Draws...\n')
    for jChain=1:nChains
        load([FileName.MCMCDraws,int2str(jChain)],'xDraws','postDraws','nDraws')
        xd(:,:,jChain) = xDraws(:,1:nThinning:end);
    end
    clear xDraws postDraws
    nDrawsUsed = (1-BurnIn)*nDraws/nThinning*nChains;
    xd = reshape(xd(:,BurnIn*nDraws/nThinning+1:end,:),np,nDrawsUsed);
    fprintf('Total number of draws per chain: %.0f\n', nDraws)
    fprintf('Thinning interval: %.0f\n', nThinning)
    fprintf('Burn in: %.0f%%\n', 100*BurnIn)
    fprintf('Total number of draws used: %.0f\n', nDrawsUsed)
elseif strcmp('priorFan',PlotParam)
    rand('state',sum(clock)*1000);
    randn('state',sum(clock)*1000);
    xd = zeros(np,nDrawsInf);
    for ii = 1:nDrawsInf
        for pp = 1:length(Params)
            eval(sprintf('xd(pp,ii)=%s;',Params(pp).priorrndcmd));
        end
    end
else
    xd = zeros(length(Params),1);
    for ii = 1:length(Params)
        eval(sprintf('xd(%i) = Params(%i).%s;',ii,ii,PlotParam));
    end
end    

%% Create Indices
nShocks2Plot = length(Shocks2Plot);
idxShocks = zeros(nShocks2Plot,1);
for ii = 1:length(Shocks2Plot)
    idxShocks(ii) = find(ismember(ShockVar,Shocks2Plot(ii)));
end    
nTempVars = length(TempVars);
idxVars = zeros(nTempVars,1);
for ii = 1:length(TempVars)
    idxVars(ii) = find(ismember([StateVar;ObsVar],TempVars(ii)));
end

%% Generate IRF
nDrawsInf = size(xd,2);
FileNameMats = FileName.Mats;
IRF = zeros(nTempVars,nSteps,nShocks2Plot,nDrawsInf);
fprintf('Generating IRFs...\n');
parfor j=1:nDrawsInf
  IRF(:,:,:,j) = MakeIRFFcn(xd(:,j),nSteps,idxShocks,idxVars,FileNameMats,...
    isObsMats,nObsVar);
end

%% Reassign IRF to Vars2Plot
PlotIRF = cell(1,nP);
if strcmp('FanChart',UseStat)
    for ii = 1:nP
        PlotIRF{ii} = IRF(ismember(Vars2Plot{ii},TempVars),:,:,:);
    end
else
    for ii = 1:nP
        PlotIRF{ii} = IRF(ismember(Vars2Plot{ii},TempVars),:,:);
    end
end    

%% IRF plot prep
XTicks = 0:TickStepIRF:(nSteps-1);
for ii = 1:length(XTicks)
    XTickLabels{ii} = sprintf('%i',XTicks(ii));
end
Fdims = {[1,1];[2,1];[3,1];[2,2];[3,2];[3,2];[4,2];[4,2];[3,3]};
if exist('Shape2Plot','var')
    if length(Shape2Plot)~=length(Vars2Plot)
        error('Shape2Plot dimensions must be specified for each Vars2Plot cell.');
    end
elseif ~exist('Shape2Plot','var')
    for ii = 1:nP
        if length(Vars2Plot{ii}) >=9
            Shape2Plot{ii} = [3,3];
        else
            Shape2Plot{ii} = Fdims{length(Vars2Plot{ii})};
        end
    end
end

%% Plot IRFs
fprintf('Plotting IRFs...\n');
tid = 0:(nSteps-1);
for js=1:length(Shocks2Plot)
    cc = 1;
    for ii = 1:nP
        nF = Shape2Plot{ii}(1)*Shape2Plot{ii}(2);
        idF = 1:nF:length(Vars2Plot{ii});
        for jF = 1:length(idF)
           if ShowFig 
               figure('Name',['Responses to innovation in ',Shocks2Plot{js}])
           else
               figure('Visible','off')
           end
           for jj=1:nF
               if idF(jF)+jj-1 <=length(Vars2Plot{ii})
                   subplot(Shape2Plot{ii}(1),Shape2Plot{ii}(2),jj);
                   Vj = idF(jF)+jj-1;
                   if ~strcmp('FanChart',UseStat)
                       vcPlot(tid,Vars2PlotScale{ii}(Vj)*PlotIRF{ii}(Vj,:,js))
                   else
                       vcPlotDistBands(squeeze(Vars2PlotScale{ii}(Vj)*PlotIRF{ii}(Vj,:,js,:))','Bands2Show',Bands2Show);
                   end
                   hold on
                   plot(tid,zeros(1,nSteps),':k')
                   hold off
                   title(Vars2PlotPretty{ii}{idF(jF)+jj-1})
                   set(gca,'Xtick',XTicks,'XTickLabel',XTickLabels);
                   axis tight
                   yl = get(gca,'YLim');
                   ySlack = max([yPerSlack*(yl(2)-yl(1)),yMaxSlack]);
                   set(gca,'YLim',[min(yl(1)-ySlack,-yMinScale),max(yMinScale,yl(2)+ySlack)])
               end
           end
%            print('-depsc',sprintf('%s_%s_Figure%i.eps',FileName.PlotsIRF,Shocks2Plot{js},cc))
           vcPrintPDF(sprintf('%s_%s_Figure%i',FileName.PlotsIRF,Shocks2Plot{js},cc))
           cc = cc+1;
        end
    end
end

%% Clean up
if ~KeepData
  clear xd postd irf IRF TempVars Fdims nP nF idF cc XTicks XTickLabels PlotIRF yl ySlack
end

%% close figures
if ~ShowFig, close all, end

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------
