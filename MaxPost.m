% MaxPost
%
% Runs several maximizations of the posterior
%
% Optional settings:
%
%   nMax
%   Sets the number of maximizations, starting at different random points.
%   Default: 1
%
%   ShowRobustness
%   If set to 1 it still shows the robustness summary in the main output. If set 
%   to zero, none of these are shown/kept. Default: 1
%
%   KeepLogsMaxPost
%   If set to 1, individual logs are kept, otherwise they are deleted after
%   conclusion of the minimization. Default: 1
%
%   KeepMatsMaxPost
%   If set to 1, individual mat files with individual minimization outputs are 
%   kept, otherwise they are deleted after conclusion of the minimization.
%   Usually once it concludes all minimization outputs are stored in a big
%   mat file in structure called MaxPostOut, so in general no need to keep 
%   these, as it will be duplicated. Use this only if MaxPostOut doesn't seem
%   right. Default: 0
%
%   DrawAll
%   If set to one, aside from user defined x0 it will simply draw remainder
%   nMax-nx0. If set to zero (default) then it adds one column to x0 (if it
%   exists already or creates x0 if it doesn't exist).
%
%   PrcGuessUsePriorDist
%   Percentage of guess draws using the prior distribution. This applies to
%   nMax-nx0. Default: 0.5 (50%).
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
%          - Changed codes to use new verbose options, and switch parallel
%            part to use matlabpool.
% Updated: October 7, 2011 by Vasco Curdia
%          - save temporary results per maximization, so that it's possible
%            to recover maximization results in case of crash.
% Updated: April 16, 2012 by Vasco Curdia
%          - Use matlabpool if available. There seems to be something wrong
%            with vcSubmitQueuedJobs.
% Updated: May 23, 2012 by Vasco Curdia
%          - Allow for page breaks and lines in between parameter blocks
%            when generating the maximization table.
%
% Copyright 2008-2012 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Set options
if ~exist('nMax','var'),nMax = 1; end
if ~isfield(FileName,'TableMaxPost')
  FileName.TableMaxPost = [FileName.Output,'TableMaxPost'];
end
if ~exist('ShowRobustness','var'),ShowRobustness = 1;end
if ~exist('KeepLogsMaxPost','var'),KeepLogsMaxPost = 1;end
if ~exist('KeepMatMaxPost','var'),KeepMatMaxPost = 0;end
if ~exist('DrawAll','var'),DrawAll=0;end
if ~exist('PrcGuessUsePriorDist','var'),PrcGuessUsePriorDist=0.5;end
if ~exist('PostArgIn','var'),PostArgIn={};end
if ~exist('UseParallel','var'),UseParallel=0;end

%% settings for parallel computing
if UseParallel
    if ~exist('nMaxWorkers','var'), nMaxWorkers = 4; end
    if ~exist('nMaxTasks','var'), nMaxTasks = 1; end
end

%% check other optional inputs for MaxPostFcn
if ~exist('MinParams','var') || ~isfield(MinParams,'verbose')
  MinParams.verbose = 1;
end
if isfield(MinParams,'Guess') && isfield(MinParams.Guess,'x0')
  nx0 = size(MinParams.Guess.x0,2);
else
  nx0 = 0;
end
if ~DrawAll
  if nx0==0
    MinParams.Guess.x0 = [Params(:).priormean]';
  else
    MinParams.Guess.x0 = [MinParams.Guess.x0,[Params(:).priormean]'];
  end
  nx0 = nx0+1;
end
nMax = max(nMax,nx0);
JobOptions = cell(nMax,1);
for jm=1:nMax
  if exist('MinParams','var'), MinParamsj = MinParams; end
  if jm<=nx0
    MinParamsj.Guess.x0 = MinParamsj.Guess.x0(:,jm);
  else
    MinParamsj.Guess = rmfield(MinParamsj.Guess,'x0');
  end
  if jm<=nx0+floor((nMax-nx0)*PrcGuessUsePriorDist);
    MinParamsj.Guess.UsePriorDist = 1;
  end
  JobOptions{jm} = {FileName.Post,Params,jm,MinParamsj,...
    sprintf('MaxPost_%03.0f.log',jm),PostArgIn};
end

%% -----------------------------------------------------------------------------

%% run maximizations

%% Display
fprintf('\n****************************')
fprintf('\n* Maximizing the posterior *')
fprintf('\n****************************\n')

%% Set Timer
TimeStr = 'MaxPost';
TimeElapsed.(TimeStr) = toc;

%% Run minimizations
MaxPostOut = cell(1,nMax);
parfor jm=1:nMax
  MaxPostOut{jm} = MaxPostFcn(JobOptions{jm}{:});
end

% if ~UseParallel
%   for jm=1:nMax
%     MaxPostOut{jm} = MaxPostFcn(JobOptions{jm}{:});
%   end
% else
%   MaxPostOut = vcSubmitQueuedJobs(nMax,@MaxPostFcn,1,JobOptions,...
%     ListPathDependencies,nMaxWorkers,'nMaxTasks',nMaxTasks,...
%     'SaveTmpFileName',[FileName.Output,'_MaxPostTmp']);
% end
save([FileName.Output,'MaxPost'])
MaxPostOut = [MaxPostOut{:}];
if ~KeepMatMaxPost
  for jm=1:nMax
    delete(sprintf('MaxPost_%03.0f.mat',jm))
  end
end
if ~KeepLogsMaxPost
  for jm=1:nMax
    delete(JobOptions{jm}{end})
  end
end
clear JobOptions
nMax = length(MaxPostOut);

%% -----------------------------------------------------------------------------

%% If desired show history evolution of robustness
if ShowRobustness
  for jm=1:nMax
    Outj = MaxPostOut(jm).MinOutput;
    fprintf('\nRobustness analysis for minimization %3.0f',jm)
    fprintf('\n----------------------------------------\n')
    if isempty(Outj)
      fprintf('%s\n%s\n\n',MaxPostOut(jm).rcMsg,MaxPostOut(jm).RrcMsg)
      continue
    end
    fprintf('Iteration %3.0f: function value: %15.8f',0,Outj(1).fh)
    %         fprintf('         %15s','')
    fprintf(' change: %15.8f',Outj(1).fh-MaxPostOut(jm).f0)
    fprintf(' stopped at iteration %.0f\n',Outj(1).itct)
    itbest.fh = Outj(1).fh;
    itbest.idx = 1;
    for jr=2:length(Outj)
      fprintf('Iteration %3.0f: function value: %15.8f',jr-1,Outj(jr).fh)
      itchg = Outj(jr).fh-itbest.fh;
      fprintf(' change: %15.8f',itchg)
      fprintf(' stopped at iteration %.0f\n',Outj(jr).itct)
      if itchg<0
        itbest.fh = Outj(jr).fh;
        itbest.idx = jr;
      end
    end
    fprintf('Best iteration: %.0f\n',itbest.idx-1)
    fprintf('Best iteration function value: %.8f\n',itbest.fh)
    fprintf('Robustness change over initial min: %.8f\n',itbest.fh-Outj(1).fh)
    fprintf('Robustness message: %s\n',MaxPostOut(jm).RrcMsg)
    fprintf('Best iteration message: %s\n\n',Outj(itbest.idx).rcMsg)
  end
end
clear Outj jm itbest itchg jr

%% -----------------------------------------------------------------------------

%% Show Starting Values
fprintf('\nGuess Values used (same order as in MaxPostOut):')
fprintf('\n================================================\n\n')
namelength = [cellfun('length',{Params(:).name})];
% show headers
str2show = '    ';
for jp=1:np
  str2show = sprintf(['%s %',int2str(max(namelength(jp),9)),'s'],...
    str2show,Params(jp).name);
end
disp(str2show)
for jm=1:nMax
  str2show = sprintf('%4.0f',jm);
  for jp=1:np
    str2show = sprintf(['%s %',int2str(max(namelength(jp),9)),'.4f'],...
      str2show,MaxPostOut(jm).x0(jp));
  end
  disp(str2show)
end
disp(' ')

%% show results for each maximization
fprintf('\nIndividual maximization results (ordered from best to worst):')
fprintf('\n=============================================================\n\n')
namelength = [cellfun('length',{Params(:).name})];
% show headers
str2show = '     log-density';
for jp=1:np
  str2show = sprintf(['%s %',int2str(max(namelength(jp),9)),'s'],...
    str2show,Params(jp).name);
end
disp(str2show)
% order maximizations
[SortPost,idxSortPost] = sort([MaxPostOut(:).f]);
% show results
for jm=1:nMax
  jShow = idxSortPost(jm);
  str2show = sprintf('%4.0f %11.4f',jShow,-MaxPostOut(jShow).f);
  for jp=1:np
    str2show = sprintf(['%s %',int2str(max(namelength(jp),9)),'.4f'],...
      str2show,MaxPostOut(jShow).x(jp));
  end
  disp(str2show)
end
disp(' ')

%% extract the best one
[Post.ModePost, idxMax] = min([MaxPostOut(:).f]);

%% Update Post
Post.ModePost = -MaxPostOut(idxMax).f;
Post.Mode = MaxPostOut(idxMax).x;
Post.Var = MaxPostOut(idxMax).H;

%% update Params
for j=1:np
  Params(j).postmode = Post.Mode(j);
  Params(j).postse = Post.Var(j,j)^(1/2);
end

%% display results on screen
fprintf('\nResults from maximization of posterior:')
fprintf('\n=======================================\n')
namelength = [cellfun('length',{Params(:).name})];
namelengthmax = max(namelength);
DispList = {'','','name';
  'Prior','dist','priordist';
  '','  mode','priormode';
  '','  mean','priormean';
  '','   se','priorse';
  '','   5%','priorp050';
  '',' median','priorp500';
  '','   95%','priorp950';
  'Posterior',' mode ','postmode';
  '','  se','postse';
  }';
nc = size(DispList,2);
for jr=1:2
  str2show = sprintf(['%-',int2str(namelengthmax),'s'],DispList{jr,1});
  str2show = sprintf('%s  %-4s',str2show,DispList{jr,2});
  for jc=3:nc
    str2show = sprintf('%s  %-7s',str2show,DispList{jr,jc});
  end
  disp(str2show)
end
for j=1:np
  str2show = sprintf(['%',int2str(namelengthmax),'s'],Params(j).(DispList{3,1}));
  str2show = sprintf('%s  %4s',str2show,Params(j).(DispList{3,2}));
  for jc=3:nc
    str2show = sprintf('%s  %7.4f',str2show,Params(j).(DispList{3,jc}));
  end
  disp(str2show)
end
fprintf('\nposterior log density at mode: %.6f\n\n',Post.ModePost)

%% save minimization output
% no need to be saved with everything else
save([FileName.Output,'MaxPostOut'],'MaxPostOut','idxMax')
clear MaxPostOut idxMax

%% clean up
clear str2show options

%% create tex table
% MakeTableMaxPost(FileName.TableMaxPost,Params,Post,o)
MakeTableMaxPost

%% Show time taken
TimeElapsed.(TimeStr) = toc-TimeElapsed.(TimeStr);
fprintf('\nMaxPost %s\n\n',vctoc([],TimeElapsed.(TimeStr)))

%% -----------------------------------------------------------------------------

