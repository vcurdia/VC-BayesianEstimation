function StateVard=MakePlotsStatesFcn(xd,MatFileName,Data,nStateVar,...
  nObsVar,nShockVar,T,sd,isDrawStates)

% MakePlotsStatesFcn
%
% Used by MakePlotsMCMCStates
%
% .........................................................................
%
% Created: February 13, 2009 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: April 29, 2012 by Vasco Curdia
%          Allow for choice of whether to draw states or not, useful in case we
%          use the median or mode -- in this case it does not make sense to draw
%          the states, as it would give a single draw, hence potentially
%          misleading.
% Updated: January 24, 2013 by Vasco Curdia
% 
% Copyright 2009-2013 by Vasco Curdia

%% ------------------------------------------------------------------------
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)*jd));
rng(sd)

nDraws = size(xd,2);
StateVard = zeros(nStateVar,T,nDraws);
for j=1:nDraws
%   [postj,StateVarttj,Pttj,Mats]=feval(PostFileName,xd(:,j));
%   StateVard(:,:,j) = DrawStates(StateVarttj,Pttj,Mats.REE.G1,Mats.REE.G2,...
%     T,nStateVar,isDrawStates);
  Mats = feval(MatFileName,xd(:,j),2,1);
  StateVard(:,:,j) = DrawStatesDK(Data'-repmat(Mats.KF.ObsVarBar,1,T),...
    Mats.KF.s00,Mats.KF.sig00,Mats.REE.G1,Mats.REE.G2,Mats.ObsEq.H,nShockVar,...
    nStateVar,nObsVar,T,isDrawStates);
end
