% MakeReportMCMC
%
% Creates a TeX document compiling plots and tables about the MCMC draws,
% including convvergence and inference analysis
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
% Updated: June 9, 2014 by Vasco Curdia
% 
% Copyright 2008-2014 by Vasco Curdia

%% ------------------------------------------------------------------------

%% preamble
np = length(Params);
if ~exist('bookmarkslevel','var'), bookmarkslevel = 1; end
if ~exist('isMoveLeft','var'), isMoveLeft = 1; end
if ~exist('ReportSections','var'), ReportSections = {'All'}; end
if ~exist('isMakePlotsStates','var'),isMakePlotsStates=1;end
if ~exist('isMakeIRF','var'),isMakeIRF = 1;end
if ~exist('isMakeVD','var'),isMakeVD = 1;end

%% Check options regarding tables
if ~exist('TableBreaks','var'), TableBreaks = 35:35:np; end
if ~exist('TableLines','var'), TableLines = []; end
TableBreaks(TableBreaks>np)=[];
if ~ismember(np,TableBreaks),TableBreaks(end+1)=np;end
nBreaks = length(TableBreaks);

%% Report Name
ReportFilename = sprintf('%sReportMCMCUpdate%.0f',FileName.Output,nUpdate);
ReportTitle = sprintf('MCMC Report:\\\\%s, Update %.0f',FileName.Output,nUpdate);

%% ------------------------------------------------------------------------

%% Begin tex file
fid=fopen(sprintf('%s.tex',ReportFilename),'wt');
fprintf(fid,'\n\\documentclass[12pt]{article}\n');
fprintf(fid,'\\usepackage{amsmath}\n');
fprintf(fid,'\\usepackage{graphicx}\n');
fprintf(fid,'\\usepackage[dvipsnames,usenames]{color}\n');
fprintf(fid,'\\usepackage[pdftex,pdfstartview=Fit,pdfpagelayout=SinglePage,bookmarksopen,');
fprintf(fid,'bookmarksopenlevel=%.0f,colorlinks,linkcolor=MyDarkBlue,',bookmarkslevel);
fprintf(fid,'citecolor=MyDarkBlue,urlcolor=MyDarkBlue,filecolor=MyDarkBlue,');
fprintf(fid,'naturalnames]{hyperref}\n');
    
fprintf(fid,'\\usepackage{fancyhdr}\n');


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

fprintf(fid,'\\begin{document}\n');

%% title page
fprintf(fid,'\\title{%s}\n',ReportTitle);
fprintf(fid,'\\maketitle\n');
fprintf(fid,'\\thispagestyle{empty}\n');

%% Description of MCMC sample
% fprintf(fid,'\\section{MCMC sample properties}\n');
fprintf(fid,'\\begin{eqnarray*} \n');
fprintf(fid,'\\begin{tabular}{rl} \n');
fprintf(fid,'number of chains: & %.0f\\\\\n',nChains);
fprintf(fid,'size of each chain: & %.0f\\\\\n',nDraws);
fprintf(fid,'burn in used: & %.0f (%.0f\\%%)\\\\\n',BurnIn*nDraws,BurnIn*100);
fprintf(fid,'thinning used: & %.0f\\\\\n',nThinning);
fprintf(fid,'\\\\log-marginal likelihood: & %.4f\n',Post.LogMgLikelihood);
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{eqnarray*}\n');
fprintf(fid,'\\newpage\n');

%% table of contents
fprintf(fid,'\\tableofcontents\n');
fprintf(fid,'\\newpage\n');

%% Parameter Moments Table
if any(ismember({'All','ParameterMoments'},ReportSections))
  fprintf(fid,'\\section{Parameter Moments}\n');
  idxPar = 0;
  for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):TableBreaks(jBreak);
    fprintf(fid,'\\begin{eqnarray*} \n');
    if isMoveLeft
      fprintf(fid,'\\hspace{-0.5in}\n');
    end
    fprintf(fid,'\\begin{tabular}{ccccccccccccc} \n');
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& & \\multicolumn{4}{c}{Prior} & & \\multicolumn{6}{c}{Posterior} \\\\[0.5ex]\n');
    fprintf(fid,'& & Dist & 5\\%% & Median & 95\\%% & & Mode & Mean & SE & 5\\%% & Median & 95\\%% \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    DispList = {'priorp050','priorp500','priorp950','postmode','postmean','postse',...
      'postp050','postp500','postp950'};
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
    fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\end{eqnarray*}\n');
    fprintf(fid,'\\newpage\n');
  end
end

%% Prior Post plots
if any(ismember({'All','PriorPost'},ReportSections))
  fprintf(fid,'\\section{Prior vs Posterior}\n');
  for jF=1:ceil(np/(FigShape{1}*FigShape{2}))
    fprintf(fid,'\\subsection{Fig %.0f}\n',jF);
    fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
    fprintf(fid,'\\label{PriorPost_Fig%.0f}\n',jF);
    fprintf(fid,'\\includegraphics[scale=1]{%s%sFig%.0f.pdf}\n',...
      PlotDir.PriorPost,FileName.PlotsPriorPost,jF);
    fprintf(fid,'\\end{figure}\n');
    fprintf(fid,'\\newpage \n');
  end
end

%% State plots
if any(ismember({'All','States'},ReportSections))&&isMakePlotsStates
  fprintf(fid,'\\section{States}\n');
  for jF=1:ceil(nStateVar/(FigShapeStates{1}*FigShapeStates{2}))
    if FigShapeStates{1}==1 && FigShapeStates{2}==1
      FigTxt = StateVar{jF};
    else
      FigTxt = sprintf('%.0f',jF);
    end
    fprintf(fid,'\\subsection{Fig %s}\n',strrep(FigTxt,'_',''));
    fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
    fprintf(fid,'\\label{States_Fig%s}\n',FigTxt);
    fprintf(fid,'\\includegraphics[scale=1]{%s%sFig%s.pdf}\n',...
      PlotDir.States,FileName.PlotsStates,FigTxt);
    fprintf(fid,'\\end{figure}\n');
    fprintf(fid,'\\newpage \n');
  end
end

%% IRF plots
if any(ismember({'All','IRF'},ReportSections))&&isMakeIRF
    fprintf(fid,'\\section{IRF}\n');
    for jPanel=1:nPanels
        for jShock=1:nShocks2Show
            fprintf(fid,'\\subsection{Response of  %s to %s}\n',...
                    Panels{jPanel},Shocks2Show{jShock});
            fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
            fprintf(fid,'\\label{IRF_%s_%s}\n',...
                    Panels{jPanel},Shocks2Show{jShock});
            fprintf(fid,'\\includegraphics[scale=1]{%s%s_%s_%s.pdf}\n',...
                    PlotDir.IRF,FileName.PlotsIRF,...
                    Panels{jPanel},Shocks2Show{jShock});
            fprintf(fid,'\\end{figure}\n');
            fprintf(fid,'\\newpage \n');
        end
    end
end

%% VD plots
if any(ismember({'All','VD'},ReportSections))&&isMakeVD
    fprintf(fid,'\\section{Variance Decomposition}\n');
    for jPanel=1:nPanels
        fprintf(fid,'\\subsection{%s}\n',Panels{jPanel});
        fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
        fprintf(fid,'\\label{VD_%s}\n',Panels{jPanel});
        fprintf(fid,'\\includegraphics[scale=1]{%s%s_%s.pdf}\n',...
                PlotDir.VD,FileName.PlotsVD,Panels{jPanel});
        fprintf(fid,'\\end{figure}\n');
        fprintf(fid,'\\newpage \n');
    end
end

%% Convergence table 1
if any(ismember({'All','ConvergenceT1'},ReportSections))
  fprintf(fid,'\\section{Convergence}\n');
  fprintf(fid,'\\subsection{Table 1}\n');
  idxPar = 0;
  for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):TableBreaks(jBreak);
    fprintf(fid,'\\begin{eqnarray*} \n');
    fprintf(fid,['\\begin{tabular}{crr',repmat('r',1,nChains),'} \n']);
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    nc = 2+nChains;
    str2show = '& \\multicolumn{1}{c}{$\\hat{R}$} & \\multicolumn{1}{c}{$mn_{eff}$} ';
    for j=1:nChains
      str2show = [str2show,'& \\multicolumn{1}{c}{$n_{eff}(',int2str(j),')$} '];
    end
    str2show = [str2show,'\\\\ \n'];
    fprintf(fid,str2show);
    fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
    % Body of table
    DispList = {'Post.ConvR(jr)','Post.Convmn_eff(jr)'};
    for jChain=1:nChains, DispList{end+1} = sprintf('Post.Convn_eff(jr,%.0f)',jChain); end
    for jr=idxPar
      str2show = ['$',strrep(Params(jr).prettyname,'\','\\'),'$'];
      str2show = sprintf('%s & %.4f',str2show,eval(DispList{1}));
      for jc=2:nc
        str2show = sprintf('%s & %.0f',str2show,eval(DispList{jc}));
      end
      str2show=[str2show,' \\\\\n'];
      fprintf(fid,str2show);
      if ismember(jr,TableLines) && jr~=idxPar(end)
        fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
      end        
    end
    % finish table
    fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\end{eqnarray*}\n');
    fprintf(fid,'\\newpage\n');
  end
end

%% Convergence table 2
if any(ismember({'All','ConvergenceT2'},ReportSections))
  fprintf(fid,'\\subsection{Table 2}\n');
  idxPar = 0;
  for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):TableBreaks(jBreak);
    fprintf(fid,'\\begin{eqnarray*} \n');
    fprintf(fid,['\\begin{tabular}{c',repmat('rl',1,nChains),'} \n']);
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    str2show = '';
    for j=1:nChains
      str2show = [str2show,'& \\multicolumn{1}{c}{$SPM_{',int2str(nSPMmeans),'}(',int2str(j),')$} & '];
    end
    str2show = [str2show,'\\\\ \n'];
    fprintf(fid,str2show);
    fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
    SPMc95 = chi2inv(0.95,nSPMmeans-1);
    SPMc99 = chi2inv(0.99,nSPMmeans-1);
    for jr=idxPar
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
      if ismember(jr,TableLines) && jr~=idxPar(end)
        fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
      end        
    end
    fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\end{eqnarray*}\n');
    fprintf(fid,'The null hypothesis of the SPM test is that the mean in two\n');
    fprintf(fid,'separate subsamples is the same. * indicates pvalue less than 5\\%%.\n');
    fprintf(fid,'** indicates pvalue less than 1\\%%.\n');
    fprintf(fid,'\\newpage\n');
  end
end

%% Convergence plots
if any(ismember({'All','ConvergencePlots'},ReportSections))
  Section = {...
    'Draws','Draws';...
    'Trace Plots','Trace';...
    };
  Section = cell2struct(Section,{'Title','Label'},2);
  nS = length(Section);
  for jS=1:nS
    fprintf(fid,'\\section{%s}\n',Section(jS).Title);
    if jS==1
      fprintf(fid,'\\subsection{Posterior density}\n');
      fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
      fprintf(fid,'\\label{%s_Post}\n',Section(jS).Label);
      fprintf(fid,'\\includegraphics[scale=1]{%s%s_%s.pdf}\n',...
        PlotDir.(['MCMC',Section(jS).Label]),...
        FileName.(['PlotsMCMC',Section(jS).Label]),'post');
      fprintf(fid,'\\end{figure}\n');
      fprintf(fid,'\\newpage \n');
    end
    for jP=1:np
      fprintf(fid,'\\subsection{%s}\n',Params(jP).name);
      fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
      fprintf(fid,'\\label{%s_%s}\n',Section(jS).Label,Params(jP).name);
      fprintf(fid,'\\includegraphics[scale=1]{%s%s_%s.pdf}\n',...
        PlotDir.(['MCMC',Section(jS).Label]),...
        FileName.(['PlotsMCMC',Section(jS).Label]),Params(jP).name);
      fprintf(fid,'\\end{figure}\n');
      fprintf(fid,'\\newpage \n');
    end
  end
end

%% finish tex file
fprintf(fid,'\\end{document}\n');
fclose(fid);

%% Compile and Cleanup
pdflatex(ReportFilename)

