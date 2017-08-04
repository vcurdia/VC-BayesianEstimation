% MCMCInference
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
%   Prc2Compute
%   Percentiles to be computed for the posterior
%   Default: [1,2.5,5,10,15,20,25,30,40,50,60,70,75,80,85,90,95,97.5]
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
if ~exist('Prc2Compute','var'),Prc2Compute=[1,2.5,5,10,15,20,25,30,40,50,60,70,75,80,85,90,95,97.5,99]; end
FileName.TableMCMCInference = sprintf('%sTableMCMCInferenceUpdate%.0f',FileName.Output,nUpdate);

%% ------------------------------------------------------------------------

%% Display
fprintf('\n******************')
fprintf('\n* MCMC Inference *')
fprintf('\n******************\n')

%% Set Timer
TimeStr = sprintf('MCMCInferenceUpdate%.0f',nUpdate);
TimeElapsed.(TimeStr) = toc;

%% load the mcmc draws
fprintf('\nLoading data...\n')
for jChain=1:nChains
    load([FileName.MCMCDraws,int2str(jChain)],'xDraws','postDraws','nDraws')
    xd(:,:,jChain) = xDraws(:,1:nThinning:end);
    postd(:,:,jChain) = postDraws(:,1:nThinning:end);
end
clear xDraws postDraws
nDrawsUsed = (1-BurnIn)*nDraws/nThinning*nChains;
xd = reshape(xd(:,BurnIn*nDraws/nThinning+1:end,:),np,nDrawsUsed);
postd = reshape(postd(:,BurnIn*nDraws/nThinning+1:end,:),nDrawsUsed,1);
fprintf('Total number of draws per chain: %.0f\n', nDraws)
fprintf('Thinning interval: %.0f\n', nThinning)
fprintf('Burn in: %.0f%%\n', 100*BurnIn)
fprintf('Total number of draws used: %.0f\n', nDrawsUsed)

%% Calculate moments
Post.Mean = mean(xd,2);
Post.Var = xd-repmat(Post.Mean,1,nDrawsUsed);
Post.Var = Post.Var*Post.Var'/nDrawsUsed;
for jPrc=1:length(Prc2Compute)
    Post.(sprintf('p%03.0f',10*Prc2Compute(jPrc))) = prctile(xd',Prc2Compute(jPrc))';
end
CorrMat = zeros(np);
for jr=1:np
    for jc=1:np
        CorrMat(jr,jc) = Post.Var(jr,jc)/sqrt(Post.Var(jr,jr)*Post.Var(jc,jc));
    end
end
Post.Corr = CorrMat; clear CorrMat

%% check if we find new mode in mcmc
[NewModePost,idxMax] = max(postd);
fprintf('\nChecking if mode in the MCMC draws...\n')
fprintf('Previous mode posterior density: %.6f\n',Post.ModePost)
fprintf('Highest posterior density in MCMC draws: %.6f\n',NewModePost)
if NewModePost>Post.ModePost
    fprintf('Found MCMC draw with higher posterior density.\n')
    fprintf('Posterior mode updated!\n')
    fprintf(['%',int2str(namelengthmax),'s %-7s %-7s\n'],'','Old','New')
    namelength = [cellfun('length',{Params(:).name})];
    namelengthmax = max(namelength);
    for jp=1:np
        fprintf(['%',int2str(namelengthmax),'s %7.4f %7.4f\n'],...
            Params(jp).name,Post.Mode(jp),xd(jp,idxMax))
    end
    Post.ModePost = NewModePost;
    Post.Mode = xd(:,idxMax);
else
    fprintf('Did not find MCMC draw with higher posterior density.\n')
    fprintf('Previous posterior mode kept!\n')
end
clear NewModePost

%% update Params
for jp=1:np
    Params(jp).postmode = Post.Mode(jp);
    Params(jp).postmean = Post.Mean(jp);
    Params(jp).postse = Post.Var(jp,jp)^(1/2);
    for jPrc=1:length(Prc2Compute)
        Prcj = sprintf('p%03.0f',10*Prc2Compute(jPrc));
        Params(jp).(['post',Prcj]) = Post.(Prcj)(jp);
    end
end

%% Marginal likelihood
% fprintf('\nComputing marginal likelihood\n')
tau = 0.1:0.1:0.9;
ntau = length(tau);
% compute parameter moments in robust fashion
% InvVar = inv(Post.Var);
xdd = xd-repmat(Post.Mean,1,nDrawsUsed);
xdvar = xdd*xdd'/nDrawsUsed;
[xdvaru,xdvars,xdvarv] = svd(xdvar);
xdvarmd = min(size(xdvars));
bigev = find(diag(xdvars(1:xdvarmd,1:xdvarmd))>1e-6);
xdvardim = length(bigev);
% [rank(xdvar),np,xdvardim]
xdvarlndet = 0;
for j=1:np
  if j>xdvardim
    xdvars(j,j) = 0;
  else
    xdvarlndet = xdvarlndet+log(xdvars(j,j));
    xdvars(j,j) = 1/xdvars(j,j);
  end
end
InvVar = xdvaru*xdvars*xdvaru';
postMax = Post.ModePost;
postMean = mean(postd);
% [postMax,postMean]
% Constant terms
% lfConst = -log(tau)-np/2*log(2*pi)-1/2*log(det(Post.Var));
lfConst = -log(tau)-np/2*log(2*pi)-1/2*xdvarlndet;
chi2Crit = chi2inv(tau,np);
% Calculate the ratio of f(x)/post(x)
pw = zeros(ntau,nDrawsUsed);
nB = 10;
for jB=1:nB
%     fprintf('Computing %.0f of %.0f...\n',jB,nB)
    for jd=1:nDrawsUsed/nB
        idx = nDrawsUsed/nB*(jB-1)+jd;
        pratio = xdd(:,idx)'*InvVar*xdd(:,idx);
        pw(:,idx) = (pratio<=chi2Crit).*exp(lfConst-1/2*pratio-postd(idx)+postMean);
    end
end
% harmonic mean
% LogMgLikelihood = postMean-log(mean(pw,2));
LogMgLikelihood = zeros(ntau,1);
for jtau=1:ntau
  pwj = pw(jtau,:);
  pwj = pwj(~isinf(pwj));
  pwj = pwj(~isnan(pwj));
  LogMgLikelihood(jtau) = postMean-log(mean(pwj,2));
end
Post.LogMgLikelihood = LogMgLikelihood(tau==0.5);
% Post.LogMgLikelihood = mean(LogMgLikelihood);
for jt = 1:ntau
    Post.LogMgLikelihoodtau(jt).tau = tau(jt);
    Post.LogMgLikelihoodtau(jt).value = LogMgLikelihood(jt);
end

%% display results on screen
fprintf('\nResults from MCMC inference:')
fprintf('\n============================\n')
namelength = [cellfun('length',{Params(:).name})];
namelengthmax = max(namelength);
DispList = {'','','name';
            'Prior','dist','priordist';
            '','   5%','priorp050';
            '',' median','priorp500';
            '','   95%','priorp950';
            'Posterior',' mode ','postmode';
            '','  mean','postmean';
            '','  se','postse';
            '','   5%','postp050';
            '',' median','postp500';
            '','   95%','postp950';
            }';
nc = size(DispList,2);
for jr=1:2
    str2show = sprintf(['%-',int2str(namelengthmax),'s'],DispList{jr,1});
    str2show = sprintf('%s  %-4s',str2show,DispList{jr,2});
    for jc=3:nc
        str2show = sprintf('%s  %-7s',str2show,DispList{jr,jc});
    end
    disp(str2show)
end
for j=1:np
    str2show = sprintf(['%',int2str(namelengthmax),'s'],Params(j).(DispList{3,1}));
    str2show = sprintf('%s  %4s',str2show,Params(j).(DispList{3,2}));
    for jc=3:nc
        str2show = sprintf('%s  %7.4f',str2show,Params(j).(DispList{3,jc}));
    end
    disp(str2show)
end
disp(' ')

% show mg Likelihood
fprintf('\nMarginal likelihood:')
fprintf('\n====================\n')
for jt = 1:ntau
    fprintf('tau = %3.1f, log-marginal likelihood = %9.4f\n',tau(jt),LogMgLikelihood(jt))
end
disp(' ')

% show correlation matrix
fprintf('\nCorrelation matrix:')
fprintf('\n===================\n')
namelength = [cellfun('length',{Params(:).name})];
namelengthmax = max(namelength);
str2show = sprintf(['%',int2str(namelengthmax),'s'],'');
for jc=1:np
    str2show = sprintf(['%s  %',int2str(namelengthmax),'s'],str2show,Params(jc).name);
end
disp(str2show)
for jr=1:np
    str2show = sprintf(['%-',int2str(namelengthmax),'s'],Params(jr).name);
    for jc=1:np
        str2show = sprintf(['%s  %',int2str(namelengthmax),'.4f'],str2show,Post.Corr(jr,jc));
    end
    disp(str2show)
end
disp(' ')

%% Clean up
clear xd postd xdd postMax pw pwj

%% create tex table
% MakeTableMCMCInference(FileName.TableMCMCInference,Params,Post,...
%     nChains,nDraws,BurnIn,nThinning,nDrawsUsed)

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------

