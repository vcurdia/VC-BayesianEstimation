function V=MCMCVarianceObsFcn(xd,Params,Mats,nObsVar,nStateVar,nShockVar,np,nDraws,...
    nHorizons,idxMat,MaxHorizon,VHorizons,isInfHorizon,isSilent)

% MCMCVarianceFcn
%
% Used by MCMCVariance
%
% See also:
% SetDSGE, GenSymVars, DataAnalysis, PriorAnalysis, GenPost, MaxPost, 
% MakeTableMaxPost, MCMC, MCMCSearchScaleFactor, MakePlotsMCMCConv, 
% MCMCInference, MakeTableMCMCInference, MakePlotsMCMCTrace, 
% MakePlotsMCMCPriorPost, MCMCConv, MakeTableMCMCConv, MCMCVD
% .........................................................................
%
% Created: April 15, 2010 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2010-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

V = zeros(nObsVar,nHorizons,nDraws);
for jd=1:nDraws
    for jp=1:np, eval(sprintf('%s = xd(jp,jd);',Params(jp).name)), end
%     Omegaj = diag(xd(end-nShockVar+1:end,jd));
    GammaBarj = eval(Mats.StateEq.GammaBar);
    Gamma0j = eval(Mats.StateEq.Gamma0);
    Gamma1j = eval(Mats.StateEq.Gamma1);
    Gamma2j = eval(Mats.StateEq.Gamma2);
    Gamma4j = eval(Mats.StateEq.Gamma4);
    Gamma3j = eye(nStateVar);
%     Gamma2j = Gamma2j*chol(Omegaj)';
    cv = (all(Gamma0j(1:nStateVar,:)==0,2)~=0);
    Gamma0j(cv,:) = -Gamma1j(cv,:);
    Gamma1j(cv,:) = Gamma4j(cv,:);
    Gamma3j(:,cv) = [];
    if ~all(all(Gamma4j(~cv,:)==0,2)),error('Incorrect system reduction'),end
    [G1j,GBarj,G2j,fmat,fwt,ywt,gev,eu]=gensys(Gamma0j,Gamma1j,GammaBarj,Gamma2j,Gamma3j);
    Hj = eval(Mats.ObsEq.H);
    for js=1:nShockVar
        Vs = G2j*idxMat(:,js)*idxMat(js,:)*G2j';
        Vo = Hj*Vs*Hj';
        for jh=2:MaxHorizon
            Vs(:,:,jh) = G1j*Vs(:,:,jh-1)*G1j'+Vs(:,:,1);
            Vo(:,:,jh) = Hj*Vs(:,:,jh)*Hj';
        end
        Vo = Vo(:,:,VHorizons(1:end-isInfHorizon));
        if isInfHorizon
            Vo(:,:,nHorizons) = Hj*real(lyapcsdSilent(G1j,Vs(:,:,1),isSilent))*Hj';
        end
        for jh=1:nHorizons
            V(:,jh,jd) = V(:,jh,jd)+diag(Vo(:,:,jh));
        end
    end
end