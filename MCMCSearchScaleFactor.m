% MCMCSearchScaleFactor
%
% Searches for the good scale factor for the MCMC jump distribution
%
% Required settings:
%
% Optional settings:
%
%   UseParallel
%   whether to use parallel or not. Default: false
% 
%   nChains
%   Sets the number of MCMC chains to use in each search block.
%   Default: 4
%
%   nDrawsSearch
%   number of MCMC draws per chain
%
%   ScaleJumpFactor
%   Initial factor. 
%   Default: 2.4
%
%   nIRS
%   size of importance resample. Default: 2000
%
%   dscale
%   double array with alternative scale changes to be tried in succession 
%   until proper rejection rate is confirmed.
%   Default: [0.2,0.05,0.01]
%
%   nConfirm
%   number of times the rejection rate in range needs to be confirmed.
%   Default: 1
%
%   RejectionRateMax
%   Maximum rejection rate
%   Default: 0.80
%
%   RejectionRateMin
%   Minimum rejection rate
%   Default: 0.70
%
%   MaxReverseDirection
%   Maximum number of times that the direction of search can be reverted
%   before attempting to change scale or simply terminate search.
%   Default: 2
%
% See also:
% SetDSGE, GenSymVars, DataAnalysis, PriorAnalysis, GenPost, MaxPost,
% MaxPostFcn, MakeTableMaxPost, MCMC, MCMCFcn, MCMCSearchScaleFactor,
% MakePlotsMCMCDraws, MCMCInference, MakeTableMCMCInference, 
% MakePlotsMCMCTrace, MakePlotsMCMCPriorPost, MCMCConv, MakeTableMCMCConv,
% MakeReportMCMCPlots, MakeReportMCMC
%
% ..............................................................................
%
% Created: March 31, 2008 by Vasco Curdia
% Updated: August 25, 2010 by Vasco Curdia
% Updated: February 14, 2011 by Vasco Curdia
%          Added minimum level to the search scale factor, to prevent it from
%          going negative.
% Updated: July 26, 2011 by Vasco Curdia
% Updated: July 29, 2011 by Vasco Curdia
%          Allow for verbose options
% 
% Copyright 2008-2011 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Settings
if ~exist('nChains','var'), nChains = 4; end
if ~exist('nDrawsSearch','var'), nDrawsSearch = 2000; end
if ~exist('nUpdate','var'), nUpdate = 0; end
if ~exist('KeepLogsMCMCSearchScale','var'), KeepLogsMCMCSearchScale = 0; end
if ~exist('MCMCOptions','var'),MCMCOptions = struct;end
FileName.MCMCDrawsSSF = sprintf('%sMCMCDrawsSSFUpdate%.0fChain',FileName.Output,nUpdate); 
MCMCSSFLogFileName = sprintf('%sMCMCDrawsSSFUpdate%.0fChain',FileName.Output,nUpdate);

%% settings regarding the scale search
if ~isfield(MCMCOptions,'ScaleJumpFactor'), MCMCOptions.ScaleJumpFactor = 2.4; end
if ~exist('dscale','var'), dscale = [0.2,0.05,0.01]; end
if ~exist('nConfirm','var'), nConfirm = 1; end
if ~exist('RejectionRateMax','var'), RejectionRateMax = 0.80; end
if ~exist('RejectionRateMin','var'), RejectionRateMin = 0.70; end
if ~exist('MaxReverseDirection','var'), MaxReverseDirection = 2; end
if ~exist('MinSearchScale','var'), MinSearchScale = 0.1; end

%% ------------------------------------------------------------------------

%% Display
fprintf('\n**************************************************')
fprintf('\n* Search Scale Factor for MCMC jump distribution *')
fprintf('\n**************************************************\n\n')

%% Set Timer
TimeStr = sprintf('MCMCSearchScaleFactorUpdate%.0f',nUpdate);
TimeElapsed.(TimeStr) = toc;

%% Search for scale factor
jConfirm = 0;
jscale = 1;
nRevertDirection = 0;
currChange = 0;
nBlocks = 0;
FileNamePost = FileName.Post;
FileNameMCMCDrawsSSF = FileName.MCMCDrawsSSF;
while jConfirm<=nConfirm 
    nBlocks = nBlocks+1;
%     BlockFileName = sprintf('%sb%02.0fc',FileName.MCMCDrawsSSF,nBlocks);
%     % update ScaleJumpFactor
%     MCMCOptionsj = MCMCOptions;
%     MCMCOptions{idxOp} = ScaleJumpFactor;
    % simulate MCMC
    parfor jChain=1:nChains
      MCMCOptionsj = MCMCOptions;
      MCMCOptionsj.fn = sprintf('%s%0.0f.log',MCMCSSFLogFileName,jChain);
      MCMCFcn(FileNamePost,sprintf('%s%.0f',FileNameMCMCDrawsSSF,jChain),Post,...
        nDrawsSearch,jChain,MCMCOptionsj);
    end
    % display rejection rates
    fprintf('\nLast chain results:\n')
    nAbove = 0;
    nBelow = 0;
    RejectionRates(nBlocks,:) = zeros(1,nChains);
    for jChain=1:nChains
        eval(sprintf('load %s%.0f nRejections',FileName.MCMCDrawsSSF,jChain))
        eval(sprintf('delete %s%.0f.mat',FileName.MCMCDrawsSSF,jChain))
        RejectionRates(nBlocks,jChain) = nRejections/nDrawsSearch;
        ScaleJumpFactors(nBlocks) = MCMCOptions.ScaleJumpFactor;
        fprintf('Block %03.0f, Chain %02.0f: ScaleJumpFactor = %4.2f, rejection rate = %5.1f%%\n',...
            nBlocks,jChain,MCMCOptions.ScaleJumpFactor,RejectionRates(nBlocks,jChain)*100)
        nAbove = nAbove + (RejectionRates(nBlocks,jChain)>RejectionRateMax);
        nBelow = nBelow + (RejectionRates(nBlocks,jChain)<RejectionRateMin);
    end
    fprintf('\nNumber of chains with rejection rate too high: %.0f',nAbove)
    fprintf('\nNumber of chains with rejection rate too low: %.0f\n\n',nBelow)
    if nAbove==nBelow
        jConfirm = jConfirm+1;
        if jConfirm>nConfirm
            fprintf('Results confirmed!\n')
        else
            fprintf('Rejection rate in range. Confirming results...\n')
        end
    else
        lastChange = currChange;
        currChange = (nAbove<nBelow)-(nAbove>nBelow);
        fprintf('Direction: %.0f\n',currChange)
        if lastChange+currChange==0
            nRevertDirection = nRevertDirection+1;
            fprintf('Direction reverted %.0f times\n',nRevertDirection)
        end
        if nRevertDirection>MaxReverseDirection
            fprintf('Direction reverted too many times. Changing increments...\n')
            jscale = jscale+1;
            nRevertDirection=0;
        end
        if jscale>length(dscale)
            jConfirm = jConfirm+1;
            if jConfirm>nConfirm
                fprintf('Results confirmed!\n')
            else
                fprintf('Increments changed too many times. Confirming results...\n')
            end
        else
            ScaleJumpFactorOld = MCMCOptions.ScaleJumpFactor;
            MCMCOptions.ScaleJumpFactor = MCMCOptions.ScaleJumpFactor +currChange*dscale(jscale);
            if MCMCOptions.ScaleJumpFactor<MinSearchScale
                MCMCOptions.ScaleJumpFactor = MinSearchScale;
                fprintf('Minimum scale breached. Setting it to minimum level: %.2f\n',...
                  MCMCOptions.ScaleJumpFactor)
            else
                fprintf('Scale changed to %.2f\n',MCMCOptions.ScaleJumpFactor)
            end
            jConfirm = (ScaleJumpFactorOld==MCMCOptions.ScaleJumpFactor);
        end            
    end
end

%% Show rejection rates:
fprintf('\nResults for ScaleJumpSearch:\n\n')
for j=1:nBlocks
    for jj=1:nChains
        fprintf('Block %03.0f, Chain %02.0f: ScaleJumpFactor = %4.2f, rejection rate = %4.1f%%\n',...
            j,jj,ScaleJumpFactors(j),RejectionRates(j,jj)*100)
    end
    fprintf('\n')
end

%% show reason to stop
fprintf('\nReason to stop: ')
if jscale>length(dscale)
    fprintf('Increments changed too many times.\n\n')
else
    fprintf('Results confirmed.\n\n')
end

%% Clean up
if ~KeepLogsMCMCSearchScale
    for jChain=1:nChains
        delete(sprintf('%s%0.0f.log',MCMCSSFLogFileName,jChain));
    end
end
clear dscale jscale jBlock jChain idxOp

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------
