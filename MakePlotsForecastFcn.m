 function [ObsVard,StateVard]=MakePlotsForecastFcn(PostFileName,xd,nStateVar,T,nForecast,MatsKF,nObsVar,Params)

% MakePlotsForecastFcn
%
% Used by MakePlotsForecast
%
% .........................................................................
%
% Created: April 26, 2009 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2009-2011 by Vasco Curdia

%% ------------------------------------------------------------------------
nDraws = size(xd,2);
StateVard = zeros(nStateVar,T+nForecast,nDraws);
ObsVard = zeros(nObsVar,T+nForecast,nDraws);
for j=1:nDraws
    [postj,StateVarttj,Pttj,G1j,G2j]=feval(PostFileName,xd(:,j));
    StateVard(:,1:T,j) = DrawStates(StateVarttj,Pttj,G1j,G2j,T,nStateVar);
    for t=1:nForecast
        StateVard(:,T+t,j) = G1j*StateVard(:,T+t-1,j);
    end
    Hd = double(subs(MatsKF.H,{Params(:).name},{xd(:,j)},0));
    HBard = double(subs(MatsKF.HBar,{Params(:).name},xd(:,j),0));
    ObsVard(:,:,j) = repmat(HBard,1,T+nForecast)+Hd*StateVard(:,:,j);
end
