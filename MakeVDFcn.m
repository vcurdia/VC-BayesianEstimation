function VD=MakeVDFcn(xd,nShockVar,idxVars,FileName,isObsMats,...
  nHorizons,idxMat,MaxHorizon,VDHorizons,isInfHorizon,isSilent)

% MakeVDFcn
%
% Used by MakeVD
%
% See also:
% SetDSGE, GenSymVars, DataAnalysis, PriorAnalysis, GenPost, MaxPost, 
% MakeTableMaxPost, MCMC, MCMCSearchScaleFactor, MakePlotsMCMCConv, 
% MCMCInference, MakeTableMCMCInference, MakePlotsMCMCTrace, 
% MakePlotsMCMCPriorPost, MCMCConv, MakeTableMCMCConv, MCMCVD
% .........................................................................
%
% Created: April 15, 2010 by Vasco Curdia
% Updated: June 9, 2014 by Vasco Curdia
% 
% Copyright 2010-2014 by Vasco Curdia

%% ------------------------------------------------------------------------

nVars = length(idxVars);
nDraws = size(xd,2);
VD = zeros(length(idxVars),nShockVar,nHorizons,nDraws);

for jd=1:nDraws
  %% Obtain Matrices
  eval(sprintf('Mats = %s(xd(:,%i),isObsMats,1);',FileName,jd));
  for js=1:nShockVar
    Vs = Mats.REE.G2*idxMat(:,js)*idxMat(js,:)*Mats.REE.G2';
    for jh=2:MaxHorizon
      Vs(:,:,jh) = Mats.REE.G1*Vs(:,:,jh-1)*Mats.REE.G1'+Vs(:,:,1);
    end
    Vs = Vs(:,:,VDHorizons(1:end-isInfHorizon));
    if isInfHorizon
      Vs(:,:,nHorizons) = real(lyapcsdSilent(Mats.REE.G1,Vs(:,:,1),isSilent));
    end
    if isObsMats>0
      for jh=1:nHorizons
        Vo(:,:,jh) = Mats.ObsEq.H*Vs(:,:,jh)*Mats.ObsEq.H';
      end
    end
    for jh=1:nHorizons
      vd = diag(Vs(:,:,jh));
      if isObsMats>0
        vd = [vd;diag(Vo(:,:,jh))];
      end
      VD(:,js,jh,jd) = vd(idxVars);
    end
  end
  VD(:,:,:,jd) = abs(VD(:,:,:,jd)./...
                     repmat(sum(VD(:,:,:,jd),2),[1,nShockVar,1,1]));
end
