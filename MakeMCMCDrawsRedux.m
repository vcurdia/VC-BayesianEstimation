% MakeMCMCDrawsRedux
%
% -------------------------------------------------------------------------
%
% Created: May 1, 2012 by Vasco Curdia
%
% Copyright 2012 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Check required settings
isError = 0;
if ~exist('nChains','var'),fprintf('Warning: nChains not specified.\n'),isError=1;end
if ~isfield(FileName,'MCMCDraws'),fprintf('Warning: FileName.MCMCDraws not specified.\n'),isError=1;end
if isError, disp('Warning: Errors found. Cannot Continue.'),return,end

%% Settings
if ~exist('BurnIn','var'),BurnIn=0.5; end
if ~exist('nDrawsRedux','var'),nDrawsRedux=1000;end

%% Display
fprintf('\nGenerating MCMC Draws Redux...\n')

%% Set Timer
TimeStr = sprintf('MCMCMakeDrawsReduxUpdate%.0f',nUpdate);
TimeElapsed.(TimeStr) = toc;

%% load the mcmc draws
nThinningRedux = max(floor(nDraws*nChains*(1-BurnIn)/nDrawsRedux),1);
for jChain=1:nChains
    load([FileName.MCMCDraws,int2str(jChain)],'xDraws','postDraws')
    xd(:,:,jChain) = xDraws(:,(BurnIn*nDraws+1):nThinningRedux:end);
    postd(:,:,jChain) = postDraws(:,(BurnIn*nDraws+1):nThinningRedux:end);
end
clear xDraws postDraws
nDrawsUsed = size(xd,2)*nChains;
xd = reshape(xd,np,nDrawsUsed);
postd = reshape(postd,nDrawsUsed,1);
fprintf('Total number of draws per chain: %.0f\n', nDraws)
fprintf('Thinning interval: %.0f\n', nThinning)
fprintf('Burn in: %.0f%%\n', 100*BurnIn)
fprintf('Total number of draws used: %.0f\n', nDrawsUsed)

%% Save MCMC draws Redux
FileName.MCMCDrawsRedux = strrep(FileName.MCMCDraws,'Chain','Redux');
save(FileName.MCMCDrawsRedux,'xd','postd','nDraws','nThinningRedux','nChains','BurnIn','nDrawsUsed')
fprintf('Saved MCMC draws redux to: %s.mat\n',FileName.MCMCDrawsRedux)

%% Clean up
clear xd postd

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% -----------------------------------------------------------------------------

