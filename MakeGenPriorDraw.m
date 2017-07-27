%% MakeGenPriorDraw
%
% Generates a function that will generate draws from the prior.
% 
% See also:
% SetDSGE, PriorAnalysis
%
% ..............................................................................
% 
% Created: July 19, 2011 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2011 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Settings

%% display
fprintf('Generating function to draw from prior...\n')
if ~isfield(FileName,'GenPriorDraw')
  FileName.GenPriorDraw = sprintf('%sGenPriorDraw',FileName.Output);
end

%% Set Timer
TimeElapsed.MakeGenPriorDraw = toc();

%% -----------------------------------------------------------------------------

%% Initiate file
fid = fopen([FileName.Output,'GenPriorDraw.m'],'wt');
fprintf(fid,'function x=%s(nDraws)\n\n',FileName.GenPriorDraw);
TodayDate = clock;
fprintf(fid,'%%%% %s\n',FileName.GenPriorDraw);
fprintf(fid,'%%\n%% Draws from the prior.\n');
fprintf(fid,'%%\n%% Created: %.0f/%.0f/%.0f with MakeMats.m \n\n',TodayDate(1:3));

%% map parameters from x to names
fprintf(fid,'if ~exist(''nDraws'',''var''), nDraws = 1; end\n');
fprintf(fid,'x = zeros(%.0f,nDraws);\n',np);
fprintf(fid,'for jd=1:nDraws\n');
for j=1:np
    fprintf(fid,'  x(%.0f,jd) = %s;\n',j,Params(j).priorrndcmd);
end
fprintf(fid,'end\n');

%% close file
fclose(fid);

%% Clean up
clear fid

%% Elapsed time
TimeElapsed.MakeGenPriorDraw = toc-TimeElapsed.MakeGenPriorDraw;

%% ------------------------------------------------------------------------
