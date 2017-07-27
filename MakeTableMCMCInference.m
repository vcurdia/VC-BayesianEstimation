function MakeTableMCMCInference(TableFilename,Params,Post,nChains,nDraws,BurnIn,nThinning,nDrawsUsed)

% MakeTableMCMCInference
%
% Creates a table in TeX that can be loaded to SWP or compiled with LaTeX.
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
% Created: April 1, 2008 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2008-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

%% preamble
np = length(Params);

%% Initiate file
fid=fopen([TableFilename,'.tex'],'wt');
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
fprintf(fid,'\\fancyhead[C]{\\textsc{%s}}\n',TableFilename);
fprintf(fid,'\\begin{document}\n');

%% heading
% fprintf(fid,'\\begin{center}\n');
fprintf(fid,'\\begin{eqnarray*} \n');
fprintf(fid,'\\begin{tabular}{rl} \n');
fprintf(fid,'number of chains: & %.0f\\\\\n',nChains);
fprintf(fid,'size of each chain: & %.0f\\\\\n',nDraws);
fprintf(fid,'burn in used: & %.0f (%.0f\\%%)\\\\\n',BurnIn*nDraws,BurnIn*100);
fprintf(fid,'thinning used: & %.0f\\\\\n',nThinning);
fprintf(fid,'number of draws used: & %.0f\\\\\\\\\n',nDrawsUsed);
fprintf(fid,'log-marginal likelihood: & %.4f\n',Post.LogMgLikelihood);
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{eqnarray*}\n');
% fprintf(fid,'\\end{center}\n');

%% Table with marginal density
% ntau = length(Post.LogMgLikelihoodtau);
% fprintf(fid,'\\begin{eqnarray*} \n');
% fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('c',1,ntau));
% fprintf(fid,'\\hline\\hline\n');
% fprintf(fid,'$\\tau$');
% for jt=1:ntau,fprintf(fid,' & %.1f',Post.LogMgLikelihoodtau(jt).tau);end
% fprintf(fid,'\\\\ \n$\\ln p(Y|\\mathcal{M})$');
% for jt=1:ntau,fprintf(fid,' & %.2f',Post.LogMgLikelihoodtau(jt).value);end
% fprintf(fid,'\\\\ \n\\hline\\hline\n');
% fprintf(fid,'\\end{tabular}\n');
% fprintf(fid,'\\end{eqnarray*}\n');

%% Begin table
fprintf(fid,'\\begin{eqnarray*} \n');
fprintf(fid,'\\begin{tabular}{ccccccccccccc} \n');
fprintf(fid,'\\hline\\hline\n');
fprintf(fid,'& & \\multicolumn{4}{c}{Prior} & & \\multicolumn{6}{c}{Posterior} \\\\\n');
fprintf(fid,'& & Dist & 5\\%% & Median & 95\\%% & & Mode & Mean & SE & 5\\%% & Median & 95\\%% \\\\\n');
fprintf(fid,'\\hline\n');
% fprintf(fid,'& & & & & & & & & & & \\\\\n');

%% Body of table
DispList = {'priorp050','priorp500','priorp950','postmode','postmean','postse',...
    'postp050','postp500','postp950'};
nc = length(DispList);
idxPost = find(ismember(DispList,'postmode'));
for jr=1:np
    str2show = ['$',strrep(Params(jr).prettyname,'\','\\'),'$'];
    str2show = [str2show,' & & ',Params(jr).priordist];
    for jc=1:nc
        if jc==idxPost,str2show = [str2show,' &'];end
        str2show = sprintf('%s & %.4f',str2show,Params(jr).(DispList{jc}));
    end
    str2show=[str2show,' \\\\\n'];
    fprintf(fid,str2show);
end

%% finish table
% fprintf(fid,'\\\\ \\multicolumn{12}{l}{} \\\\ \n');
fprintf(fid,'\\hline\\hline\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{eqnarray*}\n');

%% finish document and close file
fprintf(fid,'\\end{document}\n');
fclose(fid);

%% ------------------------------------------------------------------------
