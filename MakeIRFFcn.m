function IRF=MakeIRFFcn(xd,nSteps,idxShocks,idxVars,FileName,isObsMats,nObsVar)
%% MakeIRFFcn
%
% Function file to make IRFs. Used by MakeIRF.
%
% Inputs:
%   xd
%   Matrix of parameters, where each row corresponds to a parameter and
%   each column corresponds to a vector of parameters.
%
%   nSteps
%   Number of IRF steps.
%
%   idxShocks
%   Index vector of shock positioning.
%
%   idxVars
%   Index vector of positioning for States/Observables to plot.
%
%   FileName
%   String corresponding to name of Mats function (usually FileName.Mats)
%
%   isObsMats
%   Indicates whether to output Observable matrices and/or Gensys matrices
%   in the [Specification]Mats() function.
%   0: No observable matrices.
%   1: Function outputs observable matrices ONLY--no gensys matrices.
%   2: Function outputs observable matrices and gensys matrices, according
%   to 'isEvalGensys' options.
%
%   nobs
%   Number of observables.
%
% Outputs:
%   IRF
%   4-D matrix of irfs, where the first three dimensions are the standard
%   irf (nvar x nSteps x nshocks) and the fourth dimension represents each
%   column(parameter draw) in xd. When xd is a vector, IRF will be 3-D.
%
% See also:
% SetDSGE, MakeIRF, MakeMats
%
% .........................................................................
%
% Created: June 28, 2010 by Vasco Curdia and Ging Cee Ng
% Updated: July 26, 2011 by Vasco Curdia
%
% Copyright 2010-2011 by Vasco Curdia and Ging Cee Ng

%% ------------------------------------------------------------------------

nShocks = length(idxShocks);
IRF = zeros(length(idxVars),nSteps,nShocks,size(xd,2));

for jd = 1:size(xd,2)
  %% Obtain Matrices
  eval(sprintf('Mats = %s(xd(:,%i),isObsMats,1,1,0);',FileName,jd));
  G2 = Mats.REE.G2(:,idxShocks);
  
  %% States IRF
  irf = G2;
  for ii = 2:nSteps
    irf(:,:,ii) = Mats.REE.G1*irf(:,:,ii-1);
  end
  
  %% Observable IRF
  if isObsMats>0
    YY=zeros(nObsVar,nShocks,nSteps);
    for ii=1:nSteps
      YY(:,:,ii)=Mats.ObsEq.H*irf(:,:,ii);
    end
    irf = [irf;YY];
  end
  
  %% Format
  IRF(:,:,:,jd) = permute(irf(idxVars,:,:),[1,3,2]);
end
