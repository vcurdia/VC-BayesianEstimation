% GenPost
%
% This script generates the posterior function for the model
%
% Optional settings:
%
%   FileName.Post
%   Name of the posterior evaluation file.
%   Default: [FileName.Output,'Post']
%
%   KFinit.StateVar
%   Set the initial levels for Kalman Filter (in deviations from steady
%   state). Default: set them to zero
%
%   KFinit.sig
%   Set the initial covariance matrix for Kalman Filter.
%   Default: identity matrix
%
%   BigNumber
%   This is the value of the log-posterior (in absolute value) if either
%   the priors have value zero or if REE solution is not normal.
%   Default: 'inf'
%
%   nPreSample
%   Number of periods used to initiate the KF before start counting for the
%   likelihood. Default: 0
%
% See also:
% SetDSGE, GenSymVars, DataAnalysis, PriorAnalysis, GenPost, MaxPost,
% MaxPostFcn, MakeTableMaxPost, MCMC, MCMCFcn, MCMCSearchScaleFactor,
% MakePlotsMCMCDraws, MCMCInference, MakeTableMCMCInference,
% MakePlotsMCMCTrace, MakePlotsMCMCPriorPost, MCMCConv, MakeTableMCMCConv,
% MakeReportMCMCPlots, MakeReportMCMC, MakeMats
%
% ..............................................................................
%
% Created: March 18, 2008 by Vasco Curdia
% Updated: February 3, 2015 by Vasco Curdia
%
% Copyright 2008-2015 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% display
fprintf('Generating posterior function...\n')

%% Set Timer
TimeElapsed.GenPost = toc();

%% specify some options, if not yet specified
if ~isfield(FileName,'Post'),FileName.Post = [FileName.Output,'Post'];end
if ~exist('BigNumber','var'),BigNumber = 'inf'; end
if ~exist('nPreSample','var'),nPreSample = 0; end

%% ------------------------------------------------------------------------

%% Initiate file
fid=fopen([FileName.Post,'.m'],'wt');
fprintf(fid,'function [post,varargout]=%s(x,fid,verbose,div)\n\n',FileName.Post);
fprintf(fid,'%% Usage:\n');
fprintf(fid,'%%   post=%s(x)\n',FileName.Post);
fprintf(fid,'%%   post=%s(x,fid)\n',FileName.Post);
fprintf(fid,'%%   post=%s(x,fid,verbose)\n',FileName.Post);
fprintf(fid,'%%   post=%s(x,fid,verbose,div)\n',FileName.Post);
fprintf(fid,'%%   [post,StateVartt,SIGtt,G1,G2,H,ObsVarBar]=%s(x)\n',FileName.Post);
fprintf(fid,'%%\n%% See also:\n');
fprintf(fid,'%% MakeMats\n%%\n');
TodayDate = clock;
fprintf(fid,'%% Created: %.0f/%.0f/%.0f \n\n',TodayDate(1:3));

fprintf(fid,'%%%% Initialize  some variables\n');
fprintf(fid,'StateVartt=[];SIGtt=[];G1=[];G2=[];H=[];ObsVarBar=[];\n\n');
fprintf(fid,'if nargout>1,varargout = {StateVartt,SIGtt,G1,G2,H,ObsVarBar};end\n\n');
fprintf(fid,'OpList = {};\n');
fprintf(fid,'if nargin>=2, OpList{1} = fid; end\n');
fprintf(fid,'if nargin>=3, OpList{2} = verbose; end\n');
fprintf(fid,'if nargin==4, OpList{3} = div; end\n');

%% map parameters from x to names
fprintf(fid,'%% Map parameters\n\n');
for j=1:np
  fprintf(fid,'%s = x(%.0f);\n',Params(j).name,j);
end

%% construct prior part of the posterior
fprintf(fid,'\n%% Evaluate priors\n\n');
fprintf(fid,'post = 0;\n');
for j=1:np
  fprintf(fid,'post = post + %s;\n',Params(j).priorlpdfcmd);
end
%fprintf(fid,'post = post + %.16f;\n',Prior.LogTruncationCorrection);
fprintf(fid,'if post==-inf, post = %s; return, end;\n\n',BigNumber);

%% List the data
fprintf(fid,'\n%% Data\n\n');
fprintf(fid,'Data = [...\n');
for jeq=1:T
  fprintf(fid,'  ');
  for jc=1:nObsVar
    fprintf(fid,' %.16f',Data(jeq,jc));
    if jc==nObsVar
      fprintf(fid,';\n');
    else
      fprintf(fid,',');
    end
  end
end
fprintf(fid,'  ];\n\n');

%% Run MakeMats
if ~isfield(FileName,'Mats')
  MakeMats
end
fprintf(fid,'%%%% Get Mats\n');
fprintf(fid,'Mats=%s(x,2,1,OpList{:});\n',FileName.Mats);
% fprintf(fid,'if ~all(Mats.REE.eu==1),post = %s;return,end\n\n',BigNumber);
fprintf(fid,'if ~all(Mats.REE.eu==1)||Mats.KF.sig00rc~=0,post = %s;return,end\n\n',BigNumber);

%% ------------------------------------------------------------------------

%% Run kalman filter
% we have the system as:
%   ObsVar_t = HBar + H*StateVar_t
%   StateVar_t = GBar + G1*StateVar_tL + G2*ShockVar_t
% and we want
%   ObsVar_t = H*StateVar_t
%   StateVar_t = G1*StateVar_tL + G2*ShockVar_t
fprintf(fid,'%%%% Kalman Filter\n\n');

%% Run Kalman Filter
fprintf(fid,'stt = Mats.KF.s00;\n');
fprintf(fid,'sigtt = Mats.KF.sig00;\n');
fprintf(fid,'StateVartt = zeros(%.0f,%.0f);\n',nStateVar,T);
fprintf(fid,'SIGtt = zeros(%.0f,%.0f,%.0f);\n',nStateVar,nStateVar,T);
fprintf(fid,'for t=1:%.0f;\n',T);
fprintf(fid,'  idxNoNaN = ~isnan(Data(t,:));\n');
fprintf(fid,'  [stt,sigtt,lh,ObsVarhat]=kf(Data(t,idxNoNaN)''-Mats.KF.ObsVarBar(idxNoNaN),...\n');
fprintf(fid,'    Mats.ObsEq.H(idxNoNaN,:),stt,sigtt,Mats.REE.G1,Mats.REE.G2);\n');
fprintf(fid,'  if t>%.0f\n',nPreSample);
fprintf(fid,'    post=post+lh*[1;1];\n');
fprintf(fid,'  end\n');
fprintf(fid,'  StateVartt(:,t) = stt;\n');
fprintf(fid,'  SIGtt(:,:,t) = sigtt;\n');
fprintf(fid,'end\n\n');

%% ------------------------------------------------------------------------

%% add normalization
fprintf(fid,'%% Add normalization\n');
if nPreSample
  realdata=sum(~isnan(Data(:)))-sum(sum(~isnan(Data(1:nPreSample,:))));
else
  realdata=sum(~isnan(Data(:)));
end
fprintf(fid,'post = -( post - %.0f/2*log(2*pi) );\n\n',realdata);

%% provide additional optional output
fprintf(fid,'%%%% Provide additional optional output\n');
fprintf(fid,'if nargout>1\n');
fprintf(fid,'  varargout = {StateVartt,SIGtt,Mats};\n');
fprintf(fid,'end\n\n');

%% ------------------------------------------------------------------------

%% close file
fclose(fid);

%% Test Posterior
fprintf('Testing posterior function...\n');
post = -feval(FileName.Post,[Params(:).priormean]');
fprintf('The log-posterior value using the prior mean is %0.4f.\n\n',post);

%% Clean up
clear H MatNames nCols post varargout realdata

%% Elapsed time
TimeElapsed.GenPost = toc-TimeElapsed.GenPost;

%% ------------------------------------------------------------------------
