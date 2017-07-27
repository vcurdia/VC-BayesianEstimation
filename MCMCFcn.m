function MCMCFcn(PostName,SaveName,Post,nDraws,jChain,options)

% MaxPostFcn
%
% Maximizes the posterior by minimizing the negative of it.
%
% Usage:
%
% Out=MCMCFcn(PostName,SaveName,Post,nDraws,jChain,fn)
% Out=MCMCFcn(...,<OptionName>,<OptionValue>)
%
% Inputs: [REQUIRED]
%   
%   PostName
%   Name of posterior evaluation function
%
%   SaveName
%   Name of the mat file where to save the results.
%
%   Post
%   Structure with the current information about the posterior. Needs
%   to have fields:
%     Post.Mode
%     Post.ModePost
%     Post.Var
%   
%   nDraws
%   number of draws
%
% Optional settings should be set with a sequence of argument pairs, where the
% first is the option name (string), and the second, the option value.
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
% Created: March 27, 2008 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: July 29, 2011 by Vasco Curdia
%          Allow for verbose options
% Updated: August 3, 2011 by Vasco Curdia
%          Change the way options are used and updated
% 
% Copyright 2008-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

%% set default options
ScaleJumpFactor = 2.4;
nIRS = 1000;
nBlocks = 10;
ExplosionScale = [.4,.5,.6,.7];
ExplosionProb = [.3,.6,.8,1];
nChangednIRS = 2;
InitDrawTStudent = 1;
InitDrawItMax = 1000;
InitDrawTol = 20;
InitDrawVarFactor = 1;
InitDrawDF = 6;
SaveList = {'xDraws','postDraws','nRejections','nDraws','ScaleJumpFactor','ScaleJump','VarJump'};

%% Check options
% if ~isempty(varargin)
%   nOptions = length(varargin);
%   if mod(nOptions,2), error('Incorrect number of optional arguments.'), end
%   for jO=1:nOptions/2
%     eval(sprintf('%s = varargin{%.0f};',varargin{(jO-1)*2+1},jO*2))
%   end
% end
if exist('options','var')
  OpList = fieldnames(options);
  for j=1:length(OpList)
    eval(sprintf('%1$s=options.%1$s;',OpList{j}))
  end
end

%% Check output destination
if ~exist('fn','var') || isempty(fn) || ~ischar(fn)
  fid=1;
else
  fid=fopen(fn,'wt');
end

%% Check posterior arguments...
if ~exist('PostArgIn','var') || isempty(PostArgIn)
  PostArgIn={};
end
PostArgIn = [{fid},PostArgIn];

%% reset random number generator
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)*jChain));

%% preliminary calculations
np = length(Post.Mode);
if ~exist('isAugment','var'), isAugment = 0; end

%% ------------------------------------------------------------------------

%% generate initial draw
if ~isAugment && ~exist('x0','var')
    fprintf(fid,'Generating initial draw...\n');
    if InitDrawTStudent
        InitDrawVarChol = chol(InitDrawVarFactor^2*Post.Var)';
        InitDrawCorr = eye(np);
        for jd=1:InitDrawItMax
            x0 = Post.Mode+InitDrawVarChol*mvtrnd(InitDrawCorr,InitDrawDF)';
            post0=-feval(PostName,x0,PostArgIn{:});
            if post0>Post.ModePost-InitDrawTol
                break
            end
        end
        if post0<Post.ModePost-InitDrawTol
            fprintf(fid+(fid==1),'Warning: Initial draw has low posterior density.\n');
        end
    else
        nExplosion = length(ExplosionScale);
        jChangednIRS = 0;tic
        while jChangednIRS<=nChangednIRS
            ExplosionDraws = unifrnd(0,1,1,nIRS);
            x0c = zeros(np,nIRS);
            x0c1 = zeros(np,nIRS);
            post0c = zeros(1,nIRS);
            post0c1 = zeros(1,nIRS);
        %     VarCholT = chol(0.5^2*Post.Var)';
            for jd=1:nIRS
                for je=1:nExplosion
                    if ExplosionDraws(jd)<=ExplosionProb(je)
                        ExpV = Post.Var*ExplosionScale(je)^2;
                        break
                    end
                end
                x0c(:,jd) = mvnrnd(Post.Mode,ExpV)';
                post0c(jd) = -feval(PostName,x0c(:,jd),PostArgIn{:});
        %         x0c1(:,jd) = Post.Mode+VarCholT*mvtrnd(eye(np),x0DF)';
        %         post0c1(jd) = -feval(PostName,x0c1(:,jd));
            end
        %     post0cRaw = post0c;
            post0c = exp(post0c - repmat(Post.ModePost,1,nIRS));
            post0c = post0c./sum(post0c);
            [cpost0c,idx0]=sort(post0c);
            cpost0c = cumsum(cpost0c);
            pick0 = unifrnd(0,1);
            for j=1:nIRS
                if cpost0c(j)>=pick0
                    x0 = x0c(:,idx0(j));
                    % check whether the draw was ok...
                    WeightPick = post0c(idx0(j));
                    WeightMax = max(post0c);
                    fprintf(fid,'weight of the picked draw: %.6f\n',WeightPick);
                    fprintf(fid,'maximum weight for the draw: %.6f\n',WeightMax);
                    break
                end
            end
            if exist('x0','var')
                break
            else
                if jChangednIRS==nChangednIRS
                    x0 = Post.Mode;
                    WeightPick = NaN;
                    WeightMax = NaN;
                    fprintf(fid+(fid==1),...
                      'Warning: Did not find any suitable candidate after increasing nIRS. Using mode.\n');
                else
                    nIRS = nIRS*10;
                    jChangednIRS = jChangednIRS+1;
                    fprintf(fid+(fid==1),...
                      'Warning: Did not find any suitable candidate. Trying with increased nIRS.\n');
                end
            end
        end
    end
end

%% Prepare variables
if isAugment
    fprintf(fid,'Loading existing chain...\n');
    nDrawsNew = nDraws;
    load(SaveName)
    nDrawsOld = nDraws;
    nDraws = nDrawsNew;
    nDrawsToGen = nDraws-nDrawsOld;
    njj = ceil(nDrawsToGen/nBlocks);
    x0 = xDraws(:,end);
else
    nDrawsToGen = nDraws;
    njj = ceil(nDraws/nBlocks);
    ScaleJump = ScaleJumpFactor/np^(.5);
    VarJump = Post.Var*ScaleJump^2;
    xDraws=[];
    postDraws=[];
    nRejections = 0;
end

%% MCMC
% nj = nDraws/nBlocks;
% ScaleJump = ScaleJumpFactor/np^(.5);
% VarJump = Post.Var*ScaleJump^2;
% xDraws=[];
% postDraws=[];
% nRejections = 0;
post0=-feval(PostName,x0,PostArgIn{:});
% for jj=1:nBlocks
%     fprintf('Generating set %2.0f out of %.0f...\n',jj,nBlocks)
%     xjj=[x0,zeros(np,nj)];
%     postjj=zeros(1,nj);
%     postjj(1)=post0;
%     clear x0;
for jj=1:nBlocks
    fprintf(fid,'Generating set %2.0f out of %.0f...\n',jj,nBlocks);
    nj = min(njj,nDrawsToGen-(jj-1)*njj);
    xjj = zeros(np,nj);
    postjj = zeros(1,nj);
    for j=1:nj
        xc = mvnrnd(x0,VarJump,1)';
        postc=-feval(PostName,xc,PostArgIn{:});
        if unifrnd(0,1)<exp(postc-post0)
            x0 = xc;
            post0 = postc;
        else
            nRejections = nRejections+1;
        end
        xjj(:,j) = x0;
        postjj(j) = post0;
    end
    xDraws=[xDraws,xjj];
    postDraws=[postDraws,postjj];
    % save output
    save(SaveName,SaveList{:})
end

%% show number of rejections
% disp([int2str(nRejections),' rejections out of ',int2str(nDraws),' draws (',num2str(,'%).'])
fprintf(fid,'%.0f rejections out of %.0f draws (%.2f%%).\n',nRejections,nDraws,nRejections/nDraws*100);

%% ------------------------------------------------------------------------

%% save output
save(SaveName,SaveList{:})

%% close printed output file
if fid~=1,fclose(fid);end

%% ------------------------------------------------------------------------
