function MakeTableMCMCConv(TableFilename,Params,Post,nChains,nDraws,BurnIn,nThinning,nSPMmeans)

% MakeTableMCMCConv
%
% Creates two tables in TeX that can be loaded to SWP or compiled with 
% LaTeX.
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
% Created: April 6, 2008 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2008-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

%% preamble
np = length(Params);

%% ------------------------------------------------------------------------

%% Table 1

%% Initiate file
fid=fopen([TableFilename,'T1.tex'],'wt');
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
fprintf(fid,'\\fancyhead[C]{\\textsc{%sT1}}\n',TableFilename);
fprintf(fid,'\\begin{document}\n');

%% heading
% fprintf(fid,'\\begin{center}\n');
fprintf(fid,'\\begin{eqnarray*} \n');
fprintf(fid,'\\begin{tabular}{rl} \n');
fprintf(fid,'number of chains: & %.0f\\\\\n',nChains);
fprintf(fid,'size of each chain: & %.0f\\\\\n',nDraws);
fprintf(fid,'burn in used: & %.0f (%.0f\\%%)\\\\\n',BurnIn*nDraws,BurnIn*100);
fprintf(fid,'thinning used: & %.0f\\\\\n',nThinning);
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{eqnarray*}\n');
% fprintf(fid,'\\end{center}\n');

%% Begin table
fprintf(fid,'\\begin{eqnarray*} \n');
fprintf(fid,['\\begin{tabular}{crr',repmat('r',1,nChains),'} \n']);
fprintf(fid,'\\hline\\hline\n');
nc = 2+nChains;
str2show = '& \\multicolumn{1}{c}{$\\hat{R}$} & \\multicolumn{1}{c}{$mn_{eff}$} ';
for j=1:nChains
    str2show = [str2show,'& \\multicolumn{1}{c}{$n_{eff}(',int2str(j),')$} '];
end
str2show = [str2show,'\\\\ \n'];
fprintf(fid,str2show);
fprintf(fid,'\\hline\n');

%% Body of table
DispList = {'Post.ConvR(jr)','Post.Convmn_eff(jr)'};
for jChain=1:nChains, DispList{end+1} = sprintf('Post.Convn_eff(jr,%.0f)',jChain); end
for jr=1:np
    str2show = ['$',strrep(Params(jr).prettyname,'\','\\'),'$'];
    str2show = sprintf('%s & %.4f',str2show,eval(DispList{1}));
    for jc=2:nc
        str2show = sprintf('%s & %.0f',str2show,eval(DispList{jc}));
    end
    str2show=[str2show,' \\\\\n'];
    fprintf(fid,str2show);
end

%% finish table
fprintf(fid,'\\hline\\hline\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{eqnarray*}\n');

%% finish document and close file
fprintf(fid,'\\end{document}\n');
fclose(fid);

%% ------------------------------------------------------------------------

%% Table 2

%% Initiate file
fid=fopen([TableFilename,'T2.tex'],'wt');
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
fprintf(fid,'\\fancyhead[C]{\\textsc{%sT2}}\n',TableFilename);
fprintf(fid,'\\begin{document}\n');

%% heading
fprintf(fid,'\\begin{center}\n');
fprintf(fid,'number of chains: $%.0f$\\\\\n',nChains);
fprintf(fid,'size of each chain: $%.0f$\\\\\n',nDraws);
fprintf(fid,'burn in used: $%.0f$ ($%.0f\\%%$)\\\\\n',BurnIn*nDraws,BurnIn*100);
fprintf(fid,'thinning used: $%.0f$\\\\\n',nThinning);
fprintf(fid,'\\end{center}\n');

%% Begin table
fprintf(fid,'\\begin{eqnarray*} \n');
fprintf(fid,['\\begin{tabular}{c',repmat('rl',1,nChains),'} \n']);
fprintf(fid,'\\hline\\hline\n');
str2show = '';
for j=1:nChains
    str2show = [str2show,'& \\multicolumn{1}{c}{$SPM_{',int2str(nSPMmeans),'}(',int2str(j),')$} & '];
end
str2show = [str2show,'\\\\ \n'];
fprintf(fid,str2show);
fprintf(fid,'\\hline\n');

%% Body of table
SPMc95 = chi2inv(0.95,nSPMmeans-1);
SPMc99 = chi2inv(0.99,nSPMmeans-1);
for jr=1:np
    str2show = ['$',strrep(Params(jr).prettyname,'\','\\'),'$'];
    for jc=1:nChains
        SPMj = Post.ConvSPM(jr,jc);
        str2show = sprintf('%s & %.2f & ',str2show,SPMj);
        if SPMj>=SPMc99
            str2show = [str2show,'**'];
        elseif SPMj>=SPMc95
            str2show = [str2show,'*'];
        end
    end
    str2show=[str2show,' \\\\\n'];
    fprintf(fid,str2show);
end

%% finish table
fprintf(fid,'\\hline\\hline\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{eqnarray*}\n');

%% Legend
fprintf(fid,'The null hypothesis of the SPM test is that the mean in two\n');
fprintf(fid,'separate subsamples is the same. * indicates pvalue less than 5\\%%.\n');
fprintf(fid,'** indicates pvalue less than 1\\%%.\n');

%% finish document and close file
fprintf(fid,'\\end{document}\n');
fclose(fid);

%% ------------------------------------------------------------------------
