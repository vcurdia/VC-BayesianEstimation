% MakeReportMCMCPlots
%
% Creates a TeX document compiling plots about the MCMC draws
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
% Created: July 25, 2008 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2008-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

%% preamble
np = length(Params);
bookmarkslevel = 1;

%% Report Name
ReportFilename = strrep(FileName.PlotsMCMCConv,'Conv','');
ReportTitle = ['MCMC Plots:\\',FileName.Output];
[tf,idx] = ismember('Update',FileName.PlotsMCMCConv);
if all(tf)
    ReportTitle = [ReportTitle,', Update ',FileName.PlotsMCMCConv(idx(end)+1:end)];
end

%% Prepare sections
Section = {...
    'Draws','Draws','Conv';...
    'Trace Plots','Trace','Trace';...
    'Prior vs Posterior','PriorPost','PriorPost';...
    };
Section = cell2struct(Section,{'Title','Label','PlotID'},2);
nS = length(Section);

%% ------------------------------------------------------------------------

%% Begin tex file
fid=fopen(sprintf('%s.tex',ReportFilename),'wt');
fprintf(fid,'\n\\documentclass[12pt]{article}\n');
fprintf(fid,'\\usepackage{amsmath}\n');
% fprintf(fid,'\\usepackage{indentfirst}\n');
fprintf(fid,'\\usepackage{graphicx}\n');
fprintf(fid,'\\usepackage[dvipsnames,usenames]{color}\n');
fprintf(fid,'\\usepackage[pdftex,pdfstartview=Fit,pdfpagelayout=SinglePage,bookmarksopen,');
fprintf(fid,'bookmarksopenlevel=%.0f,colorlinks,linkcolor=MyDarkBlue,',bookmarkslevel);
fprintf(fid,'citecolor=MyDarkBlue,urlcolor=MyDarkBlue,filecolor=MyDarkBlue,');
fprintf(fid,'naturalnames]{hyperref}\n');
    
fprintf(fid,'\\usepackage[round,longnamesfirst]{natbib}\n');
fprintf(fid,'\\usepackage{fancyhdr}\n');

fprintf(fid,'%%TCIDATA{OutputFilter=LATEX.DLL}\n');
fprintf(fid,'%%TCIDATA{Version=5.50.0.2953}\n');
fprintf(fid,'%%TCIDATA{Codepage=65001}\n');
fprintf(fid,'%%TCIDATA{<META NAME="SaveForMode" CONTENT="2">}\n');
fprintf(fid,'%%TCIDATA{BibliographyScheme=BibTeX}\n');
fprintf(fid,'%%TCIDATA{Created=Tuesday, June 08, 2004 14:52:05}\n');
fprintf(fid,'%%TCIDATA{LastRevised=Wednesday, January 30, 2008 21:00:28}\n');
fprintf(fid,'%%TCIDATA{<META NAME="GraphicsSave" CONTENT="32">}\n');
fprintf(fid,'%%TCIDATA{<META NAME="DocumentShell" CONTENT="Standard LaTeX');
fprintf(fid,'\\Blank - Standard LaTeX Article">}\n');
fprintf(fid,'%%TCIDATA{Language=American English}\n');
fprintf(fid,'%%TCIDATA{CSTFile=40 LaTeX article.cst}\n');

fprintf(fid,'\\setlength{\\oddsidemargin}{0.0in}\n');
fprintf(fid,'\\setlength{\\evensidemargin}{0.0in}\n');
fprintf(fid,'\\setlength{\\topmargin}{0cm}\n');
fprintf(fid,'\\setlength{\\textheight}{8.5in}\n');
fprintf(fid,'\\setlength{\\textwidth}{6.5in}\n');
fprintf(fid,'\\setlength{\\hoffset}{0.0in}\n');
fprintf(fid,'\\setlength{\\headheight}{14.5pt}\n');
fprintf(fid,'\\pagestyle{fancy}\n');
fprintf(fid,'\\fancyhf{}\n');
% fprintf(fid,'\\fancyhead[L]{\\leftmark}\n');
fprintf(fid,'\\fancyhead[C]{\\textsc{%s}}\n',ReportFilename);
fprintf(fid,'\\fancyfoot[C]{\\thepage}\n');
fprintf(fid,'\\renewcommand{\\headrulewidth}{0pt}\n');
fprintf(fid,'\\definecolor{MyDarkBlue}{rgb}{0,0.08,0.45}\n');
fprintf(fid,'\\graphicspath{{%s/}}\n',FileName.Path);
fprintf(fid,'\\input{tcilatex}\n');

fprintf(fid,'\\begin{document}\n');

%% title page
fprintf(fid,'\\title{%s}\n',ReportTitle);
fprintf(fid,'\\maketitle\n');

fprintf(fid,'\\thispagestyle{empty}\n');
fprintf(fid,'\\newpage\n');

%% table of contents
fprintf(fid,'%%TCIMACRO{\\TeXButton{Table of Contents}{\\tableofcontents}}%%\n');
fprintf(fid,'%%BeginExpansion\n');
fprintf(fid,'\\tableofcontents%%\n');
fprintf(fid,'%%EndExpansion\n');
fprintf(fid,'\\newpage\n');

%% Insert plots in each section
for jS=1:nS
    fprintf(fid,'\\section{%s}\n',Section(jS).Title);
    if jS==1
        fprintf(fid,'\\subsection{Posterior density}\n');
        fprintf(fid,'%%TCIMACRO{\\TeXButton{B}{\\begin{figure}[htbp] \\centering}}%%\n');
        fprintf(fid,'%%BeginExpansion \n');
        fprintf(fid,'\\begin{figure}[htbp] \\centering%%\n');
        fprintf(fid,'%%EndExpansion \n');
        fprintf(fid,'\\label{%s_Post}%%\n',Section(jS).Label);
        fprintf(fid,'%%TCIMACRO{%%\n');
        fprintf(fid,'%%\\TeXButton{includegraphics}{\\includegraphics[scale=1]{%s_%s.pdf}}}%%\n',...
            FileName.(['PlotsMCMC',Section(jS).PlotID]),'post');
        fprintf(fid,'%%BeginExpansion \n');
        fprintf(fid,'\\includegraphics[scale=1]{%s_%s.pdf}%%\n',...
            FileName.(['PlotsMCMC',Section(jS).PlotID]),'post');
        fprintf(fid,'%%EndExpansion \n');
        fprintf(fid,'%%TCIMACRO{\\TeXButton{E}{\\end{figure}}}%%\n');
        fprintf(fid,'%%BeginExpansion \n');
        fprintf(fid,'\\end{figure}%%\n');
        fprintf(fid,'%%EndExpansion \n');
        fprintf(fid,'\\newpage \n');
    end
    if strcmp(Section(jS).PlotID,'PriorPost')
        for jF=1:ceil(np/(FigShape{1}*FigShape{2}))
            fprintf(fid,'\\subsection{Fig %.0f}\n',jF);
            fprintf(fid,'%%TCIMACRO{\\TeXButton{B}{\\begin{figure}[htbp] \\centering}}%%\n');
            fprintf(fid,'%%BeginExpansion \n');
            fprintf(fid,'\\begin{figure}[htbp] \\centering%%\n');
            fprintf(fid,'%%EndExpansion \n');
            fprintf(fid,'\\label{%s_Fig%.0f}%%\n',Section(jS).Label,jF);
            fprintf(fid,'%%TCIMACRO{%%\n');
            fprintf(fid,'%%\\TeXButton{includegraphics}{\\includegraphics[scale=1]{%sFig%.0f.pdf}}}%%\n',...
                FileName.(['PlotsMCMC',Section(jS).PlotID]),jF);
            fprintf(fid,'%%BeginExpansion \n');
            fprintf(fid,'\\includegraphics[scale=1]{%sFig%.0f.pdf}%%\n',...
                FileName.(['PlotsMCMC',Section(jS).PlotID]),jF);
            fprintf(fid,'%%EndExpansion \n');
            fprintf(fid,'%%TCIMACRO{\\TeXButton{E}{\\end{figure}}}%%\n');
            fprintf(fid,'%%BeginExpansion \n');
            fprintf(fid,'\\end{figure}%%\n');
            fprintf(fid,'%%EndExpansion \n');
            fprintf(fid,'\\newpage \n');
        end
    else
        for jP=1:np
            fprintf(fid,'\\subsection{%s}\n',Params(jP).name);
            fprintf(fid,'%%TCIMACRO{\\TeXButton{B}{\\begin{figure}[htbp] \\centering}}%%\n');
            fprintf(fid,'%%BeginExpansion \n');
            fprintf(fid,'\\begin{figure}[htbp] \\centering%%\n');
            fprintf(fid,'%%EndExpansion \n');
            fprintf(fid,'\\label{%s_%s}%%\n',Section(jS).Label,Params(jP).name);
            fprintf(fid,'%%TCIMACRO{%%\n');
            fprintf(fid,'%%\\TeXButton{includegraphics}{\\includegraphics[scale=1]{%s_%s.pdf}}}%%\n',...
                FileName.(['PlotsMCMC',Section(jS).PlotID]),Params(jP).name);
            fprintf(fid,'%%BeginExpansion \n');
            fprintf(fid,'\\includegraphics[scale=1]{%s_%s.pdf}%%\n',...
                FileName.(['PlotsMCMC',Section(jS).PlotID]),Params(jP).name);
            fprintf(fid,'%%EndExpansion \n');
            fprintf(fid,'%%TCIMACRO{\\TeXButton{E}{\\end{figure}}}%%\n');
            fprintf(fid,'%%BeginExpansion \n');
            fprintf(fid,'\\end{figure}%%\n');
            fprintf(fid,'%%EndExpansion \n');
            fprintf(fid,'\\newpage \n');
        end
    end
end

%% finish tex file
fprintf(fid,'\\end{document}\n');
fclose(fid);
