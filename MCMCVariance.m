% MCMCVariance
%
% Generates Variance tables.
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
%   VHorizons
%   Horizons at which to compute the V.
%   Default: [1,4,inf]
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
%   computes V at mode
%   Default: 0
%
%   UseMedian
%   computes V at median
%   Default: 0
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
% Updated: May 10, 2012 by Vasco Curdia
%          - Allow the use of MCMCDrawsRedux
%          - option to keep result matrices
% 
% Copyright 2010-2012 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Settings
if ~exist('VHorizons','var'),VHorizons = [1,4,inf]; end
if ~exist('VPrctile','var'),VPrctile = 50; end
if ~exist('UseMode','var'),UseMode=0; end
if ~exist('UseMedian','var'),UseMedian=0; end
if UseMode && UseMedian, error('Cannot use both median and mode.'), end
if ~(UseMode || UseMedian)
    if ~exist('BurnIn','var'),BurnIn = 0.25; end
    if ~exist('nThinning','var'),nThinning = 1; end
    if ~exist('Bands2Show','var'),Bands2Show=[50,60,70,80,90]; end
end
if ~exist('UseMCMCDrawsRedux','var'),UseMCMCDrawsRedux=0; end
if ~exist('KeepResultMats','var'),KeepResultMats = 0; end
if ~exist('isSilent','var'),isSilent = 0; end


%% Check required settings
if ~(UseMode || UseMedian || UseMCMCDrawsRedux)
    isError = 0;
    if ~exist('nChains','var'),fprintf('Warning: nChains not specified.\n'),isError=1;end
    if ~isfield(FileName,'MCMCDraws'),fprintf('Warning: FileName.MCMCDraws not specified.\n'),isError=1;end
    if isError, disp('Warning: Errors found. Cannot Continue.'),return,end
end

%% settings for parallel computing
if ~exist('UseParallel','var'),UseParallel = 0; end
if ~exist('nMaxWorkers','var'), nMaxWorkers = 4; end

%% ------------------------------------------------------------------------

%% Display
fprintf('\n*****************')
fprintf('\n* MCMC Variance *')
fprintf('\n*****************\n')

%% Set Timer
TimeStr = 'V';
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
      postd(:,:,jChain) = postDraws(:,1:nThinning:end);
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

%% Generate V
isInfHorizon = ismember(inf,VHorizons);
VHorizons = sort(VHorizons);
nHorizons = length(VHorizons);
MaxHorizon = VHorizons(end-isInfHorizon);
idxMat = eye(nShockVar);
nDrawsUsed = size(xd,2);
V = zeros(nStateVar,nHorizons,nDrawsUsed);
parfor jd=1:nDrawsUsed
    V(:,:,jd) = MCMCVarianceFcn(xd(:,jd),Params,SymMats,nStateVar,nShockVar,np,1,...
        nHorizons,idxMat,MaxHorizon,VHorizons,isInfHorizon,isSilent);
end
if ~(UseMode || UseMedian)
    V = prctile(V,VPrctile,3);
end

%% Create table
fprintf('\nVariance:')
fprintf('\n=========\n')
if UseMode
    fprintf('(at the mode)\n')
elseif UseMedian
    fprintf('(at the median)\n')
else
    fprintf('(Percentile: %.1f)\n',VPrctile)
end
znamelength = [cellfun('length',StateVar)];
znamelengthmax = max(znamelength);
fprintf(['%-',int2str(znamelengthmax),'s',repmat('   %7.0f',1,nHorizons),'\n'],'Horizon',VHorizons(:))
for jz=1:nStateVar
    fprintf(['%-',int2str(znamelengthmax),'s',repmat('   %7.3f',1,nHorizons),'\n'],StateVar{jz},100^2*V(jz,:))
end

%% Clean up
clear xd postd Vs
% for jp=1:np, eval(sprintf('syms %s',Params(jp).name)), end
if ~KeepResultMats
  clear V
end


%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------

