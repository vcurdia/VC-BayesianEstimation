% function MakeTableMaxPost(TableFilename,Params,Post,o)
% 
% MakeTableMaxPost
%
% Creates a table in TeX that can be loaded to SWP or compiled with LaTeX.
%
% Optional settings:
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
% Created: March 26, 2008 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: August 22, 2011 by Vasco Curdia
%          - Added PDFLaTeX copilation and auxiliary file names deletion.
% Updated: May 23, 2012 by Vasco Curdia
%          - shifted from function to script
%          - Allow for page breaks and lines in between parameter blocks
%            when generating the maximization table.
% 
% Copyright 2008-2012 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Check options
if ~exist('TableBreaks','var'), TableBreaks = 35:35:np; end
if ~exist('TableLines','var'), TableLines = []; end
TableBreaks(TableBreaks>np)=[];
if ~ismember(np,TableBreaks),TableBreaks(end+1)=np;end
nBreaks = length(TableBreaks);

%% Initiate file
fid=fopen([FileName.TableMaxPost,'.tex'],'wt');
fprintf(fid,'\n\\documentclass{article}\n');
fprintf(fid,'\\setlength{\\oddsidemargin}{0.0in}\n');
fprintf(fid,'\\setlength{\\evensidemargin}{0.0in}\n');
fprintf(fid,'\\setlength{\\topmargin}{-0.5cm}\n');
fprintf(fid,'\\setlength{\\textheight}{9in}\n');
fprintf(fid,'\\setlength{\\textwidth}{6.5in}\n');
fprintf(fid,'\\usepackage{fancyhdr}\n');
fprintf(fid,'\\renewcommand{\\headrulewidth}{0pt}\n');
fprintf(fid,'\\pagestyle{fancy}\n');
fprintf(fid,'\\fancyhf{}\n');
fprintf(fid,'\\fancyhead[C]{\\textsc{%s}}\n',FileName.TableMaxPost);
fprintf(fid,'\\begin{document}\n');

%% value of the posterior density
% fprintf(fid,'\\begin{center}\n');
fprintf(fid,'\\begin{eqnarray*} \n');
fprintf(fid,'\\begin{tabular}{rl} \n');
fprintf(fid,'posterior log density at mode: & %.4f\n',Post.ModePost);
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{eqnarray*}\n');
% fprintf(fid,'\\end{center}\n');

%% check if need page breaks
idxPar = 0;
for jBreak=1:nBreaks
  idxPar = (idxPar(end)+1):TableBreaks(jBreak);
  fprintf(fid,'\\begin{eqnarray*} \n');
  fprintf(fid,'\\begin{tabular}{cccccccccccc} \n');
  fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
  fprintf(fid,'& & \\multicolumn{7}{c}{Prior} & & \\multicolumn{2}{c}{Posterior} \\\\[0.5ex]\n');
  fprintf(fid,'& & Dist & Mode & Mean & SE & 5\\%% & Median & 95\\%% & & Mode & SE \n');
  fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
  % fprintf(fid,'& & & & & & & & & & & \\\\\n');
  DispList = {'priormode','priormean','priorse','priorp050','priorp500',...
      'priorp950','postmode','postse'};
  nc = length(DispList);
  idxPost = find(ismember(DispList,'postmode'));
  for jr=idxPar
      str2show = ['$',strrep(Params(jr).prettyname,'\','\\'),'$'];
      str2show = [str2show,' & & ',Params(jr).priordist];
      for jc=1:nc
          if jc==idxPost,str2show = [str2show,' &'];end
          str2show = sprintf('%s & %.4f',str2show,Params(jr).(DispList{jc}));
      end
      str2show=[str2show,' \\\\\n'];
      fprintf(fid,str2show);
      if ismember(jr,TableLines) && jr~=idxPar(end)
        fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
      end        
  end
  % fprintf(fid,'\\\\ \\multicolumn{12}{l}{} \\\\ \n');
  fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
  % fprintf(fid,['\\multicolumn{12}{l}{* For the igam ',...
  %     'distribution, the shape parameter is set to 2 when SE is $\\infty$.} \n']);
  fprintf(fid,'\\end{tabular} \n');
  fprintf(fid,'\\end{eqnarray*} \n');
  fprintf(fid,'\\newpage\n');
end

%% finish document and close file
fprintf(fid,'\\end{document}\n');
fclose(fid);

%% Compile and Cleanup
pdflatex(FileName.TableMaxPost)

%% ------------------------------------------------------------------------
