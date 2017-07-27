%% MakeMats
%
% Generate matrices for gensys.m
%
% gensys requires the following form:
%   Gamma0*z_t = GammaBar + Gamma1*z_t_1 + Gamma2*eps_t + Gamma3*eta_t
% and yields system in form:
%   z_t = GBar + G1*z_t_1 + G2*eps_t
%
% The observation equations are given by:
%   ObsVar_t = HBar + H*StateVar_t
%
% See also:
% GenPost,gensys,vcgensys
%
% ..............................................................................
%
% Created: June 30, 2010 by Vasco Curdia and Ging Cee Ng
% Updated: July 26, 2011 by Vasco Curdia
% Updated: July 28, 2011 by Vasco Curdia
%   - Allow for additional options to be fed into gensys
% Updated: October 17, 2013 by Vasco Curdia
%   - Generates KF matrices.
% Updated: January 2, 2014 by Vasco Curdia
%  - added varargin to pass optional arguments to vcgensys
%
% Copyright 2010-2014 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Settings
if ~exist('GensysAuthor','var'),GensysAuthor='CS';end

%% Rename preexisting Mats structure
if exist('Mats','var') && ~exist('MatsSym','var')
  MatsSym = Mats;
  clear Mats;
end

%% display
fprintf('Making Mats for Gensys...\n')
FileName.Mats = sprintf('%sMats',FileName.Output);

%% Set Timer
TimeElapsed.MakeMats = toc();

%% -----------------------------------------------------------------------------

%% Initiate file
fidMats=fopen([FileName.Output,'Mats.m'],'wt');
fprintf(fidMats,...
  'function Mats=%s(x,isObsMats,isEvalGensys,fid,verbose,varargin)\n\n',...
  FileName.Mats);
fprintf(fidMats,'%% Options:\n');
fprintf(fidMats,'%%   1. [Mats]=%s(x)\n',FileName.Mats);
fprintf(fidMats,'%%   2. [Mats]=%s(x,isObsMats)\n',FileName.Mats);
fprintf(fidMats,'%%   3. [Mats]=%s(x,isObsMats,isEvalGensys)\n',FileName.Mats);
fprintf(fidMats,'%%   4. [Mats]=%s(x,isObsMats,isEvalGensys,fid)\n',...
  FileName.Mats);
fprintf(fidMats,'%%   5. [Mats]=%s(...,varargin)\n',FileName.Mats);
fprintf(fidMats,['%%   Default: Option 1 where Mats is a structure with ',...
  'fields ''ObsEq'',''StateEq'', and ''REE''\n']);
fprintf(fidMats,'%%\n');
fprintf(fidMats,'%% Inputs:\n');
fprintf(fidMats,'%%   x\n');
fprintf(fidMats,'%%   Vector of parameters\n%%\n');
fprintf(fidMats,'%%   isObsMats (Optional)\n');
fprintf(fidMats,['%%   Indicates whether to output Observable matrices ',...
  'and/or Gensys matrices.\n']);
fprintf(fidMats,'%%   0: No observable matrices.\n');
fprintf(fidMats,'%%   1: Function outputs observable matrices ONLY--no gensys matrices.\n');
fprintf(fidMats,'%%   2: Function outputs observable matrices and gensys matrices, according to ''isEvalGensys'' options.\n');
fprintf(fidMats,'%%   Default: 2\n%%\n');
fprintf(fidMats,'%%   isEvalGensys (Optional)\n');
fprintf(fidMats,'%%   Flags whether to evalute gensys and determines what function will output from gensys.\n');
fprintf(fidMats,'%%   0: Do not evaluate gensys. Function outputs gensys input matrices only.\n');
fprintf(fidMats,'%%   1: Evalute gensys. Function outputs gensys output matrices only.\n');
fprintf(fidMats,'%%   2: Evalute gensys. Function outputs both gensys input and output matrices.\n');
fprintf(fidMats,'%%   Default: 1\n%%\n');
fprintf(fidMats,'%%   varargin (optional)\n');
fprintf(fidMats,'%%   Allows for arguments for verbose and flie id to pass to gensys if needed.\n');
fprintf(fidMats,'%% Outputs:\n');
fprintf(fidMats,'%%   Mats is a structure with fields dependent on the values of the input flags:\n');
fprintf(fidMats,'%%      ObsEq.HBar,ObsEq.H\n');
fprintf(fidMats,'%%          Transition matrices from observable equation\n%%\n');
fprintf(fidMats,'%%      REE.GBar,REE.G1,REE.G2,REE.eu\n');
fprintf(fidMats,'%%          Gensys output matrices, including the error code vector eu.\n%%\n');
fprintf(fidMats,'%%      StateEq.Gamma0,StateEq.Gamma1,StateEq.GammaBar,StateEq.Gamma2,StateEq.Gamma3\n');
fprintf(fidMats,'%%          Gensys inputs\n%%\n');
fprintf(fidMats,'%% Gensys requires the following form:\n');
fprintf(fidMats,'%%   Gamma0*StateVar_t = GammaBar + Gamma1*StateVar_t_1 + Gamma2*ShockVar_t + Gamma3*eta_t\n');
fprintf(fidMats,'%% and yields system in form:\n');
fprintf(fidMats,'%%   StateVar_t = GBar + G1*StateVar_t_1 + G2*ShockVar_t\n');
fprintf(fidMats,'%%\n');
fprintf(fidMats,'%% The observation equations are given by:\n');
fprintf(fidMats,'%%   ObsVar_t = HBar + H*StateVar_t\n');
TodayDate = clock;
fprintf(fidMats,'%%\n%% Created: %.0f/%.0f/%.0f with MakeMats.m \n\n',TodayDate(1:3));

%% Set default input values
fprintf(fidMats,'%%%% Set default input values\n');
fprintf(fidMats,'if nargin<2, isObsMats=2; end\n');
fprintf(fidMats,'if nargin<3, isEvalGensys=1; end\n');
fprintf(fidMats,'if nargin<4, fid = 1; end\n\n');
fprintf(fidMats,'if nargin<5, verbose = 0; end\n\n');

%% map parameters from x to names
fprintf(fidMats,'%%%% Map parameters\n');
for j=1:np
  fprintf(fidMats,'%s = x(%.0f);\n',Params(j).name,j);
end

%% generate matrices for gensys
% gensys requires the following form:
%   Gamma0*StateVar_t = GammaBar + Gamma1*StateVar_tL + Gamma2*ShockVar_t 
%                       + Gamma3*eta_t
% and yields system in form:
%   StateVar_t = GBar + G1*StateVar_t_1 + G2*ShockVar_t
fprintf(fidMats,'\n%%%% Create Gensys Inputs\n');

%% raw matrices
SymMats.StateEq.GammaBar = jacobian(StateEq,one);
SymMats.StateEq.Gamma0 = -jacobian(StateEq,StateVar_tF);
SymMats.StateEq.Gamma1 = jacobian(StateEq,StateVar_t);
SymMats.StateEq.Gamma4 = jacobian(StateEq,StateVar_tL);
SymMats.StateEq.Gamma2 = jacobian(StateEq,ShockVar_t);
MatNames = fieldnames(SymMats.StateEq);
nCols = [1,nStateVar,nStateVar,nStateVar,nShockVar];
for jM=1:length(MatNames)
  fprintf(fidMats,'%s = [...\n',MatNames{jM});
  for jeq=1:nStateVar
    fprintf(fidMats,'   ');
    for jc=1:nCols(jM)
      fprintf(fidMats,' %s',char(eval(sprintf('SymMats.StateEq.%s(jeq,jc)',...
        MatNames{jM}))));
      if jc==nCols(jM)
        fprintf(fidMats,';\n');
      else
        fprintf(fidMats,',');
      end
    end
  end
  fprintf(fidMats,'    ];\n\n');
end
fprintf(fidMats,'Gamma3 = eye(%.0f);\n\n',nStateVar);

%% Check matrices and reassign them
fprintf(fidMats,'cv = (all(Gamma0(1:%.0f,:)==0,2)~=0);\n',nStateVar);
fprintf(fidMats,'Gamma0(cv,:) = -Gamma1(cv,:);\n');
fprintf(fidMats,'Gamma1(cv,:) = Gamma4(cv,:);\n');
fprintf(fidMats,'Gamma3(:,cv) = [];\n');
fprintf(fidMats,'if ~all(all(Gamma4(~cv,:)==0,2))\n');
fprintf(fidMats,'  error(''Incorrect system reduction'')\n');
fprintf(fidMats,'end\n\n');

%% ------------------------------------------------------------------------

%% Generate raw matrices for obs eq
fprintf(fidMats,'%%%% Generate matrices for observation equations\n');
fprintf(fidMats,'if isObsMats>0\n');
H0 = -jacobian(ObsEq,ObsVar_t);
SymMats.ObsEq.HBar = jacobian(ObsEq,one);
SymMats.ObsEq.H = jacobian(ObsEq,StateVar_t);
SymMats.ObsEq.HBar = H0\SymMats.ObsEq.HBar;
SymMats.ObsEq.H = H0\SymMats.ObsEq.H;
MatNames = fieldnames(SymMats.ObsEq);
nCols = [1,nStateVar];
for jM=1:length(MatNames)
  fprintf(fidMats,'   %s = [...\n',MatNames{jM});
  for jeq=1:nObsVar
    fprintf(fidMats,'   ');
    for jc=1:nCols(jM)
      fprintf(fidMats,' %s',char(eval(sprintf('SymMats.ObsEq.%s(jeq,jc)',...
        MatNames{jM}))));
      if jc==nCols(jM)
        fprintf(fidMats,';\n');
      else
        fprintf(fidMats,',');
      end
    end
  end
  fprintf(fidMats,'    ];\n');
end
fprintf(fidMats,'end\n');

%% ------------------------------------------------------------------------

%% run gensys
fprintf(fidMats,'\n%%%% Run gensys\n');
fprintf(fidMats,'if isEvalGensys>0\n');
fprintf(fidMats,...
  '  [G1,GBar,G2,fmat,fwt,ywt,gev,eu]=vcgensys(Gamma0,Gamma1,GammaBar,...\n');
fprintf(fidMats,...
  '  Gamma2,Gamma3,''%s'',fid,verbose,varargin{:});\n',GensysAuthor);
fprintf(fidMats,'end\n');

%% ------------------------------------------------------------------------

%% Kalman Filter
fprintf(fidMats,'if isObsMats>0 && isEvalGensys>0\n');

%% get constants
% fprintf(fid,'StateVarBar = inv(eye(%.0f)-G1)*GBar;\n',nStateVar);
% fprintf(fid,'ObsVarBar = HBar + H*zBar;\n\n');
fprintf(fidMats,'if all(GBar(:)==0)\n');
fprintf(fidMats,'  StateVarBar = zeros(%.0f,1);\n',nStateVar);
fprintf(fidMats,'else\n');
fprintf(fidMats,'  StateVarBar = (eye(%.0f)-G1)\\GBar;\n',nStateVar);
fprintf(fidMats,'end\n');
fprintf(fidMats,'ObsVarBar = HBar + H*StateVarBar;\n\n');

%% Initial level for Kalman Filter
if exist('KFinit','var') && isfield(KFinit,'s')
  fprintf(fidMats,'s00 = [...\n');
  for jeq=1:nStateVar
    fprintf(fidMats,'    %.16f;\n',KFinit.s(jeq));
  end
  fprintf(fidMats,'    ];\n\n');
else
  fprintf(fidMats,'s00 = zeros(%.0f,1);\n\n',nStateVar);
end

%% Initial variance for Kalman Filter
if exist('KFinit','var') && isfield(KFinit,'sig')
  fprintf(fidMats,'sig00 = [...\n');
  for jeq=1:nStateVar
    fprintf(fidMats,'   ');
    for jc=1:nStateVar
      fprintf(fidMats,' %0.16f',KFinit.sig(jeq,jc));
      if jc==nStateVar
        fprintf(fidMats,';\n');
      else
        fprintf(fidMats,',');
      end
    end
  end
  fprintf(fidMats,'    ];\n\n');
  fprintf(fidMats,'sig00rc = 0;\n');
else
  fprintf(fidMats,'[sig00,sig00rc]=lyapcsd(G1,G2*G2'');\n');
  fprintf(fidMats,'sig00 = real(sig00);sig00 = (sig00+sig00'')/2;\n');
  fprintf(fidMats,'if sig00rc~=0\n');
  fprintf(fidMats,'  if verbose\n');
  fprintf(fidMats,['    fprintf(fid,''Warning: Could not find ',...
    'unconditional variance!\\n'');\n']);
  fprintf(fidMats,'  end\n\n');
  fprintf(fidMats,'end\n\n');
end
fprintf(fidMats,'end\n');

%% ------------------------------------------------------------------------

%% Assign Output Vars
fprintf(fidMats,'\n%%%% Assign Output\n');
fprintf(fidMats,'if isObsMats>0\n');
fprintf(fidMats,'  Mats.ObsEq.HBar = HBar;\n');
fprintf(fidMats,'  Mats.ObsEq.H = H;\n');
fprintf(fidMats,'end\n');
fprintf(fidMats,'if isObsMats>0 && isEvalGensys>0\n');
fprintf(fidMats,'  Mats.KF.StateVarBar = StateVarBar;\n');
fprintf(fidMats,'  Mats.KF.ObsVarBar = ObsVarBar;\n');
fprintf(fidMats,'  Mats.KF.s00 = s00;\n');
fprintf(fidMats,'  Mats.KF.sig00 = sig00;\n');
fprintf(fidMats,'  Mats.KF.sig00rc = sig00rc;\n');
fprintf(fidMats,'end\n');
fprintf(fidMats,'if isObsMats~=1 && isEvalGensys>0\n');
fprintf(fidMats,'  Mats.REE.GBar=GBar;\n');
fprintf(fidMats,'  Mats.REE.G1=G1;\n');
fprintf(fidMats,'  Mats.REE.G2=G2;\n');
fprintf(fidMats,'  Mats.REE.eu=eu;\n');
fprintf(fidMats,'end\n');
fprintf(fidMats,'if isObsMats~=1 && isEvalGensys~=1\n');
fprintf(fidMats,'   Mats.StateEq.Gamma0=Gamma0;\n');
fprintf(fidMats,'   Mats.StateEq.Gamma1=Gamma1;\n');
fprintf(fidMats,'   Mats.StateEq.GammaBar=GammaBar;\n');
fprintf(fidMats,'   Mats.StateEq.Gamma2=Gamma2;\n');
fprintf(fidMats,'   Mats.StateEq.Gamma3=Gamma3;\n');
fprintf(fidMats,'end\n');

%% ------------------------------------------------------------------------

%% close file
fclose(fidMats);

%% Clean up
clear MatNames nCols H0

%% Elapsed time
TimeElapsed.MakeMats = toc-TimeElapsed.MakeMats;

%% ------------------------------------------------------------------------
