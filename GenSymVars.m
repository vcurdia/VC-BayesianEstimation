% GenSymVars
%
% Creates symbolic variables and vectors for Parameters, ObsVar, StateVar and
% ShockVar(including the leads and lags)
%
% See also:
% SetDSGE, GenSymVars, DataAnalysis, PriorAnalysis, GenPost, MaxPost, 
% MakeTableMaxPost, MCMC, MCMCSearchScaleFactor, MakePlotsMCMCConv, 
% MCMCInference, MakeTableMCMCInference, MakePlotsMCMCTrace, 
% MakePlotsMCMCPriorPost, MCMCConv, MakeTableMCMCConv, MCMCVD
%
% .........................................................................
%
% Created: March 17, 2008 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: May 23, 2012 by Vasco Curdia
%          - Now it generates ParLists
% 
% Copyright 2008-2012 by Vasco Curdia

%% ------------------------------------------------------------------------

%% display
fprintf('Generating required symbolic variables and structures...\n')

%% Set Timer
TimeElapsed.GenSymVars = toc();

%% Parameters
if size(Params,2)==4, Params(:,5) = Params(:,1); end
Params = cell2struct(Params,{'name','priordist','priormean','priorse','prettyname'},2);
np = length(Params);
for j=1:np, eval(['syms ',Params(j).name]), end
ParList = {Params(:).name};

%% Observation variables
nObsVar = length(ObsVar);
ObsVar_t = sym(zeros(1,nObsVar));
for j=1:nObsVar
    eval(sprintf('syms %s_t',ObsVar{j}))
    ObsVar_t(j) = eval([ObsVar{j},'_t']);
end

%% State space variables
nStateVar = length(StateVar);
StateVar_t = sym(zeros(1,nStateVar)); StateVar_tF = sym(zeros(1,nStateVar)); 
StateVar_tL = sym(zeros(1,nStateVar));
for j=1:nStateVar
    eval(sprintf('syms %1$s_t %1$s_tF %1$s_tL',StateVar{j}))
    StateVar_t(j) = eval([StateVar{j},'_t']);
    StateVar_tF(j) = eval([StateVar{j},'_tF']);
    StateVar_tL(j) = eval([StateVar{j},'_tL']);
end

%% Shocks
nShockVar = length(ShockVar);
ShockVar_t = sym(zeros(1,nShockVar));
for j=1:nShockVar
    eval(['syms ',ShockVar{j},'_t'])
    ShockVar_t(j) = eval([ShockVar{j},'_t']);
end

%% constant variable
syms one

%% Time Elapsed
TimeElapsed.GenSymVars = toc-TimeElapsed.GenSymVars;

