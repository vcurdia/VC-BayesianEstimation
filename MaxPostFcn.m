function Out=MaxPostFcn(PostFileName,Params,jMax,MinParams,fn,PostArgIn)

% MaxPostFcn
%
% Maximizes the posterior by minimizing the negative of it.
%
% Inputs: [REQUIRED]
%   
%   PostFileName
%   Name of posterior evaluation function
%
%   Params
%   Params structure
%
% See also:
% SetDSGE, GenSymVars, DataAnalysis, PriorAnalysis, GenPost, MaxPost,
% MaxPostFcn, MakeTableMaxPost, MCMC, MCMCFcn, MCMCSearchScaleFactor,
% MakePlotsMCMCDraws, MCMCInference, MakeTableMCMCInference, 
% MakePlotsMCMCTrace, MakePlotsMCMCPriorPost, MCMCConv, MakeTableMCMCConv,
% MakeReportMCMCPlots, MakeReportMCMC
%
% ..............................................................................
%
% Created: March 25, 2008 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: July 29, 2011 by Vasco Curdia
%          Changed codes to use new verbose options, and switch parallel part to
%          use matlabpool.
% Updated: October 7, 2011 by Vasco Curdia
%          save temporary results per maximization, so that it's possible to
%          recover maximization results in case of crash.
% Updated: January 2, 2014 by Vasco Curdia
%  - added varargin to pass optional arguments to vcgensys
% 
% Copyright 2008-2014 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Check output destination
if ~exist('fn','var') || isempty(fn) || ~ischar(fn)
  fid=1;
else
  fid=fopen(fn,'wt');
end
if ~exist('PostArgIn','var') || isempty(PostArgIn)
  PostArgIn={};
end
PostArgIn = [{fid},PostArgIn];
MinParams.varargin = PostArgIn;

%% Set name of matfile to save temporary results
MinParams.MatFn = sprintf('MaxPost_%03.0f',jMax);

%% reset random number generator
% s = clock;s = (s(6)-round(s(6)))*jMax*100000;rand('state',s);
% s = clock;s = (s(6)-round(s(6)))*jMax*100000;randn('state',s);
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)*jMax));

%% default options
nMaxDraws = 1000;
UsePriorDist = 0;
x0Mean = [Params(:).priormean]';
x0SE = [Params(:).priorse];
x0df = 4;

%% Update options
if isfield(MinParams,'Guess')
    GuessFields = fieldnames(MinParams.Guess);
    for jf=1:length(GuessFields)
        eval([GuessFields{jf},'=MinParams.Guess.',GuessFields{jf},';'])
    end
    MinParams = rmfield(MinParams,'Guess');
end

%% Guess values
np = length(Params);
if ~exist('x0','var')
    for jg=1:nMaxDraws
        if UsePriorDist
            fprintf(fid,'Guess using prior distribution.\n');
            for j=1:np
                x0(j,1) = eval(Params(j).priorrndcmd);
            end
        else
            fprintf(fid,'Guess using t-student.\n');
            for jp=1:np
                for jpg=1:nMaxDraws*10
                    x0j = x0Mean(jp)+x0SE(jp)*trnd(x0df);
                    if strcmp(Params(jp).priordist,'N')
                        break
                    elseif strcmp(Params(jp).priordist,'B') && x0j>0 && x0j<1
                        break
                    elseif strcmp(Params(jp).priordist,'G') && x0j>0
                        break
                    elseif strcmp(Params(jp).priordist,'IG1') && x0j>0
                        break
                    elseif strcmp(Params(jp).priordist,'IG2') && x0j>0
                        break
                    end
                end
                x0(jp,1) = x0j;
            end
        end
        f0 = feval(PostFileName,x0,PostArgIn{:});
        if f0>1e50
            fprintf(fid+(fid==1),'WARNING: Bad initial guess. Drawing a new one.\n');
        else
            break
        end
    end
end


%% initial inverse hessian
if ~isfield(MinParams,'H0')
    H0 = [Params(:).priorse].^2;
    idx = (H0==inf); H0(idx) = [Params(idx).priormean];
    MinParams.H0 = diag(H0);
end

%% maximize the posterior
Out = vcRobustMin(PostFileName,x0,MinParams,fid);

%% close output file
if fid~=1,fclose(fid);end

%% -----------------------------------------------------------------------------

