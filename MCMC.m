% MCMC
%
% Runs several chains of MCMC
%
% Required settings:
%
% Optional settings:
%
%   UseParallel
%   whether to use parallel or not. Default: false
% 
%   nChains
%   Sets the number of MCMC chains.
%   Default: 4
%
%   nDraws
%   number of MCMC draws per chain
%
% See also:
% SetDSGE, GenSymVars, DataAnalysis, PriorAnalysis, GenPost, MaxPost, 
% MakeTableMaxPost, MCMC, MCMCSearchScaleFactor, MakePlotsMCMCConv, 
% MCMCInference, MakeTableMCMCInference, MakePlotsMCMCTrace, 
% MakePlotsMCMCPriorPost, MCMCConv, MakeTableMCMCConv, MCMCVD
%
% .........................................................................
%
% Created: March 27, 2008 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: July 29, 2011 by Vasco Curdia
%          Allow for verbose options
% 
% Copyright 2008-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Settings
if ~exist('nChains','var'), nChains = 4; end
if ~exist('nDraws','var'), nDraws = 50000; end
if ~exist('nUpdate','var'), nUpdate = 0; end
if ~exist('KeepLogsMCMC','var'), KeepLogsMCMC = 1; end
if ~exist('MCMCOptions','var'),MCMCOptions = struct;end
if ~exist('UseOneInitDrawAtMode','var'), UseOneInitDrawAtMode=0; end
FileName.MCMCDraws = sprintf('%sMCMCDrawsUpdate%.0fChain',FileName.Output,nUpdate); 
MCMCLogFileName = sprintf('%sMCMCDrawsUpdate%.0fChain',FileName.Output,nUpdate);

%% ------------------------------------------------------------------------

%% run chains

%% Display
fprintf('\n*****************************')
fprintf('\n* Generating MCMC sample(s) *')
fprintf('\n*****************************\n\n')

%% Set Timer
TimeStr = sprintf('MCMCUpdate%.0f',nUpdate);
TimeElapsed.(TimeStr) = toc;

%% Prepare options
JobOptions = cell(nChains,1);
for jChain=1:nChains
  MCMCOptionsj = MCMCOptions;
  MCMCOptionsj.fn = sprintf('%s%.0f.log',MCMCLogFileName,jChain);
  if UseOneInitDrawAtMode && jChain==1
    MCMCOptionsj.x0 = [Params(:).postmode]';
  end
  JobOptions{jChain} = {FileName.Post,sprintf('%s%.0f',FileName.MCMCDraws,jChain),...
    Post,nDraws,jChain,MCMCOptionsj};
end
clear MCMCOptionsj

%% Run MCMCFcn
parfor jChain=1:nChains
  MCMCFcn(JobOptions{jChain}{:});
end
    
%% show rejection rates
for jChain=1:nChains
    eval(sprintf('load %s%.0f nDraws nRejections ScaleJumpFactor',FileName.MCMCDraws,jChain))
    fprintf('\nChain %02.0f: ScaleJumpFactor = %4.2f, rejection rate = %5.1f%%',...
        jChain,ScaleJumpFactor,nRejections/nDraws*100)
end
disp(' ')

%% ------------------------------------------------------------------------

%% Clean up
if ~KeepLogsMCMC
    for jChain=1:nChains
        delete(sprintf('%s%.0f.log',MCMCLogFileName,jChain));
    end
end
clear options

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------
