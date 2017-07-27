function HistDecomp=MakePlotsHistDecompFcn(PostFileName,xd,nStateVar,...
  nShockVar,nObsVar,nVars2Show,idxHD,T,jd,isDrawStates)

% MakePlotsHistDecompFcn
%
% Used by MakePlotsHistDecomp
%
% .........................................................................
%
% Created: April 29, 2012 by Vasco Curdia
% 
% Copyright 2012 by Vasco Curdia

%% ------------------------------------------------------------------------

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)*jd));

nDraws = size(xd,2);
HistDecomp = zeros(nVars2Show,T,nShockVar+1,nDraws);
isObs = any(idxHD.ObsVar(:));
for j=1:nDraws
    [postj,StateVarttj,Pttj,G1j,G2j]=feval(PostFileName,xd(:,j));
    StateVarj = DrawStates(StateVarttj,Pttj,G1j,G2j,T,nStateVar,isDrawStates);
    Shockj = [zeros(nShockVar,1),...
      G2j(idxHD.G2,:)\(StateVarj(idxHD.G2,2:end)-G1j(idxHD.G2,:)*StateVarj(:,1:end-1))];
    StateVarDecompj=zeros(nStateVar,nShockVar+1,T);
    StateVarDecompj(:,end,1) = StateVarj(:,1);
    for t=2:T
      StateVarDecompj(:,:,t) = G1j*StateVarDecompj(:,:,t-1)+...
        [G2j*diag(Shockj(:,t)),zeros(nStateVar,1)];
    end
    if isObs
      Hj = RPRatioMats(xd(:,j),1);
      Hj = Hj.ObsEq.H;
      ObsVarDecompj = zeros(nObsVar,nShockVar+1,T);
      for t=1:T
        ObsVarDecompj(:,:,t) = Hj*StateVarDecompj(:,:,t);
      end
      ObsVarDecompj = permute(ObsVarDecompj,[1,3,2]);
    end
    StateVarDecompj = permute(StateVarDecompj,[1,3,2]);
    for jV=1:nVars2Show
      if idxHD.StateVar(jV)
        HistDecomp(jV,:,:,j) = StateVarDecompj(idxHD.StateVar(jV),:,:);
      else
        HistDecomp(jV,:,:,j) = ObsVarDecompj(idxHD.ObsVar(jV),:,:);
      end
    end
end

