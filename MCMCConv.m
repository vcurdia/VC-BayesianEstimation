% MCMCConv
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
%   nSPMmeans
%   Number of means for the SPM test.
%   Default: 4
%
%   DrawsFraction
%   Fraction of the draws used to compute the number of autocovariance
%   matrices for the spectral density needed.
%   Default: 0.04
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
% Created: April 6, 2008 by Vasco Curdia
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
FileName.TableMCMCConv = sprintf('%sTableMCMCConvergenceUpdate%.0f',FileName.Output,nUpdate);
if ~exist('nSPMmeans','var'),nSPMmeans = 4;end
if ~exist('DrawsFraction','var'),DrawsFraction = 0.04; end

%% ------------------------------------------------------------------------

%% Display
fprintf('\n********************')
fprintf('\n* MCMC Convergence *')
fprintf('\n********************\n')

%% Set Timer
TimeStr = sprintf('MCMCConvUpdate%.0f',nUpdate);
TimeElapsed.(TimeStr) = toc;

%% load the mcmc draws
fprintf('\nLoading data...\n')
for jChain=1:nChains
    load([FileName.MCMCDraws,int2str(jChain)],'xDraws','nDraws')
    xd(:,:,jChain) = xDraws(:,1:nThinning:end);
end
clear xDraws
nDrawsUsed = (1-BurnIn)*nDraws/nThinning;
xd = xd(:,BurnIn*nDraws/nThinning+1:end,:);
fprintf('Total number of draws per chain: %.0f\n',nDraws)
fprintf('Thinning interval: %.0f\n',nThinning)
fprintf('Burn in: %.0f%%\n\n',100*BurnIn)
namelength = [cellfun('length',{Params(:).name})];
namelengthmax = max(namelength);

%% check convergence with R statistic from Gelman et al
Mj = mean(xd(:,:,:),2);
M = mean(Mj,3);
s2j = sum((xd(:,:,:)-repmat(Mj,[1,nDrawsUsed,1])).^2,2)/(nDrawsUsed-1);
W = mean(s2j,3);
B = nDrawsUsed/(nChains-1)*sum((Mj-repmat(M,[1,1,nChains])).^2,3);
varplus = (nDrawsUsed-1)/nDrawsUsed*W+1/nDrawsUsed*B;
fprintf('\nR statistic for convergence analysis:')
fprintf('\n=====================================\n\n')
Post.ConvR = squeeze(sqrt(varplus./W));
for jp=1:np
    fprintf(['%',int2str(namelengthmax),'s %7.4f\n'],Params(jp).name,Post.ConvR(jp))
end
fprintf('\n')

%% Effective number of independent draws a la Gelman et al.
fprintf('\nEffective number of independent draws a la Gelman et al.:')
fprintf('\n=========================================================\n\n')
Post.Convmn_eff = squeeze(nChains*nDrawsUsed*varplus./B);
for jp=1:np
    fprintf(['%',int2str(namelengthmax),'s %10.0f\n'],Params(jp).name,Post.Convmn_eff(jp))
end
fprintf('\n')
clear omegadc Mj M s2j W B varplus

%% Geweke's separated partial means test
npm = floor(nDrawsUsed/2/nSPMmeans);
Lp = round(DrawsFraction*npm);
for jChain=1:nChains
    for jp=1:nSPMmeans
        xj = xd(:,(2*jp-1)*npm+1:2*jp*npm,jChain);
        mean_jp(:,jp) = mean(xj,2);
        xj = xj-repmat(mean_jp(:,jp),1,npm);
        s0 = mean(xj.^2,2);
        for jL=1:Lp-1
            cL = sum(xj(:,1+jL:npm).*xj(:,1:npm-jL),2)/npm;
            s0 = s0 + 2*(Lp-jL)/Lp*cL;
        end
        S0(:,jp) = s0; 
    end
    for j=1:np
        hp = mean_jp(j,2:end)'-mean_jp(j,1:end-1)';
        Vp = diag(S0(j,2:end))+diag(S0(j,1:end-1));
        Vp = Vp - [zeros(nSPMmeans-1,1) [diag(S0(j,2:end-1));zeros(1,nSPMmeans-2)]];
        Vp = Vp - [zeros(1,nSPMmeans-1);[diag(S0(j,2:end-1)) zeros(nSPMmeans-2,1)]];
        Vp = Vp/npm;
        Post.ConvSPM(j,jChain) = hp'*inv(Vp)*hp;
    end
end
fprintf('\nSPM test results for each chain:')
fprintf('\n================================')
fprintf('\ncritical value (95%%) for chi-square(%.0f) is %f\n\n',nSPMmeans-1,chi2inv(0.95,nSPMmeans-1))
for jp=1:np
    fprintf(['%',int2str(namelengthmax),'s %10.4f %10.4f %10.4f %10.4f\n'],...
        Params(jp).name,Post.ConvSPM(jp,:))
end
fprintf('\n')
clear npm Lp xj mean_jp s0 cL S0 hp Vp

%% Calculate the within chain number of effective sample size
L = round(DrawsFraction*nDrawsUsed);
for jChain=1:nChains
    xj = xd(:,:,jChain);
    xj = xj-repmat(mean(xj,2),1,nDrawsUsed);
    c0 = mean(xj.^2,2);
    S0 = c0;
    for jL=1:L-1
        cL = sum(xj(:,1+jL:nDrawsUsed).*xj(:,1:nDrawsUsed-jL),2)/nDrawsUsed;
        S0 = S0 + 2*(L-jL)/L*cL;
    end
    Post.Convn_eff(:,jChain) = nDrawsUsed*c0./S0;
end
fprintf('\nn_eff for each chain:')
fprintf('\n=====================\n\n')
for jp=1:np
    fprintf(['%',int2str(namelengthmax),'s %10.0f %10.0f %10.0f %10.0f\n'],...
        Params(jp).name,Post.Convn_eff(jp,:))
end
fprintf('\n')
clear L xj c0 S0

%% Create tables
% MakeTableMCMCConv(FileName.TableMCMCConv,Params,Post,nChains,nDraws,BurnIn,nThinning,nSPMmeans)

%% Clean up
clear xd

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\n%s %s\n\n',TimeStr,vctoc([],TimeElapsed.(TimeStr)))

%% ------------------------------------------------------------------------
