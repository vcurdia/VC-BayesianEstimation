function varargout=DrawStates(StateVartt,SIGtt,G1,G2,T,nStateVar,isDrawStates)

% CRDrawStates
%
% Draws states using the Carter-Kohn procedure.
%
% Usage:
%   StateVard = DrawStates(StateVartt,SIGtt,G1,G2,T,isDrawStates)
%   [StateVartTn,SIGtTn] = DrawStates(...)
%   [StateVard,StateVartTn,SIGtTn] = DrawStates(...)
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
% Created: October 21, 2008 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: April 29, 2012 by Vasco Curdia
%          Allow for choice of whether to draw states or not, useful in case we
%          use the median or mode -- in this case it does not make sense to draw
%          the states, as it would give a single draw, hence potentially
%          misleading.
% 
% Copyright 2008-2012 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Run smoother and draw state
StateVartTn = StateVartt';
SIGtTn = SIGtt;
Om = G2*G2';
StateVard = StateVartTn;
if isDrawStates
  StateVard(T,:)=DrawState(StateVartTn(T,:),SIGtTn(:,:,T),nStateVar);
end
for t=T-1:-1:1
    [StateVartTnt,SIGtTnt]=vcDistSmooth(StateVartTn(t,:),SIGtTn(:,:,t),StateVard(t+1,:),G1,Om);
    if isDrawStates
      StateVard(t,:) = DrawState(StateVartTnt,SIGtTnt,nStateVar);
    else
      StateVard(t,:) = StateVartTnt;
    end
    StateVartTn(t,:) = StateVartTnt;
    SIGtTn(:,:,t) = SIGtTnt;
end

%% prepare output
if nargout==1
    varargout = {StateVard'};
elseif nargout==2
    varargout = {StateVartTn',SIGtTn};
else
    varargout = {StateVard',StateVartTn',SIGtTn};
end

%% ------------------------------------------------------------------------

function StateVard=DrawState(StateVartTn,SIGtTn,nStateVar)

% NOTE: the trick below is intended at rounding errors that could make the
%       matrix look like it is not semi-definite positive.

% Old method
% % NumPrecision = 1e-16;
% [uu,dd,vv]=svd(SIGtTn);
% % [uu,dd]=eig(SIGtTn);
% % dd=round(dd/NumPrecision)*NumPrecision;
% StateVard=StateVartTn+mvnrnd(zeros(nStateVar,1),dd)*uu';

% ptT = vcCovDefPos(ptT);
% StateVard=mvnrnd(StateVartT,ptT);

[u,d,v]=svd(SIGtTn);
StateVard=StateVartTn+mvnrnd(zeros(nStateVar,1),d.^(1/2))*u';

%% ------------------------------------------------------------------------
