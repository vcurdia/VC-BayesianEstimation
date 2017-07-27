% PriorAnalysis
%
% Analyzes the priors
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
% Created: March 18, 2008 by Vasco Curdia
% Updated: November 12, 2008 by Vasco Curdia
% Updated: July 19, 2011 by Vasco Curdia
%          Allows for calibrated parameters.
% Updated: July 26, 2011 by Vasco Curdia
% Updated: April 19, 2012 by Vasco Curdia
%          Added additional distribution: Truncated Normal, assuming
%          truncation at zero, so that x>0.
% Updated: April 26, 2012 by Vasco Curdia
%          Fixed TN pdf, cdf, inverse-cdf, and random number generator.
% 
% Copyright 2008-2012 by Vasco Curdia

%% ------------------------------------------------------------------------

%% display
fprintf('Analyzing priors...\n')

%% Set Timer
TimeElapsed.PriorAnalysis = toc;

%% check options/defaults
if ~exist('ShowPriorTable','var'), ShowPriorTable = 1; end

%% list percentiles to generate
Percentiles = [0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99];
nprc = length(Percentiles);

%% Analyze priors
np = length(Params);
for j=1:np
  if strcmp(Params(j).priordist,'N')
    pmean = Params(j).priormean;
    pse = Params(j).priorse;
    Params(j).priormode = pmean;
    for jprc=1:nprc
      eval(sprintf('Params(%.0f).priorp%03.0f = norminv(%f,%.16f,%.16f);',...
        j,[1000,1]*Percentiles(jprc),pmean,pse))
    end
    Params(j).priorparams = {pmean,pse};
    Params(j).priorlpdfcmd = sprintf('log(normpdf(%s,%.16f,%.16f))',Params(j).name,pmean,pse);
    Params(j).priorpdfcmd = sprintf('normpdf(%%.16f,%.16f,%.16f)',pmean,pse);
    Params(j).priorrndcmd = sprintf('normrnd(%.16f,%.16f)',pmean,pse);
  elseif strcmp(Params(j).priordist,'TN')
    % Assume x>=0
    pmean = Params(j).priormean;
    pse = Params(j).priorse;
    a = -pmean/pse;
    acdf = normcdf(a,0,1);
    aZ = 1-normcdf(a,0,1);
    alambda = normpdf(a,0,1)/aZ;
    adelta = alambda*(alambda-a);
    Params(j).priormean = pmean + pse*alambda;
    Params(j).priorse = pse*(1-adelta)^(1/2);
    Params(j).priormode = max(0,pmean);
    for jprc=1:nprc
      eval(sprintf('Params(%.0f).priorp%03.0f = norminv(%f,%.16f,%.16f);',...
        j,1000*Percentiles(jprc),Percentiles(jprc)*aZ+acdf,pmean,pse))
    end
    Params(j).priorparams = {pmean,pse};
    Params(j).priorlpdfcmd = sprintf(...
      'log((%1$s>=0)*normpdf((%1$s-%2$.16f)/%3$.16f,0,1)/%3$.16f/%4$.16f)',...
      Params(j).name,pmean,pse,aZ);
    Params(j).priorpdfcmd = sprintf(...
      '(%%1$.16f>=0)*normpdf((%%1$.16f-%1$.16f)/%2$.16f,0,1)/%2$.16f/%3$.16f',...
      pmean,pse,aZ);
    Params(j).priorrndcmd = sprintf('norminv(rand*%.16f+%.16f,%.16f,%.16f)',...
      aZ,acdf,pmean,pse);
    %         Params(j).priorrndcmd = sprintf('normrnd(%.16f,%.16f)',pmean,pse);
    clear alambda adelta aZ acdf
  elseif strcmp(Params(j).priordist,'B')
    pmean = Params(j).priormean;
    pse = Params(j).priorse;
    a = pmean*(pmean-pmean^2-pse^2)/pse^2;
    b = a*(1/pmean-1);
    Params(j).priormode = min(max(0,(a-1)/(a+b-2)),1);
    for jprc=1:nprc
      eval(sprintf('Params(%.0f).priorp%03.0f = betainv(%f,%.16f,%.16f);',...
        j,[1000,1]*Percentiles(jprc),a,b))
    end
    Params(j).priorparams = {a,b};
    Params(j).priorlpdfcmd = sprintf('log(betapdf(%s,%.16f,%.16f))',Params(j).name,a,b);
    Params(j).priorpdfcmd = sprintf('betapdf(%%.16f,%.16f,%.16f)',a,b);
    Params(j).priorrndcmd = sprintf('betarnd(%.16f,%.16f)',a,b);
  elseif strcmp(Params(j).priordist,'G')
    pmean = Params(j).priormean;
    pse = Params(j).priorse;
    a = (pmean/pse)^2;
    b = pmean/a;
    if a>=1
      Params(j).priormode = (a-1)*b;
    else
      Params(j).priormode = NaN;
    end
    for jprc=1:nprc
      eval(sprintf('Params(%.0f).priorp%03.0f = gaminv(%f,%.16f,%.16f);',...
        j,[1000,1]*Percentiles(jprc),a,b))
    end
    Params(j).priorparams = {a,b};
    Params(j).priorlpdfcmd = sprintf('log(gampdf(%s,%.16f,%.16f))',Params(j).name,a,b);
    Params(j).priorpdfcmd = sprintf('gampdf(%%.16f,%.16f,%.16f)',a,b);
    Params(j).priorrndcmd = sprintf('gamrnd(%.16f,%.16f)',a,b);
  elseif strcmp(Params(j).priordist,'IG1')
    pmean = Params(j).priormean;
    pse = Params(j).priorse;
    if pse==inf
      a = 1;
    else
      fname = sprintf('igamsolve%.0f',cputime*1e10);
      fid=fopen([fname,'.m'],'wt');
      fprintf(fid,'function f=%s(x)\n',fname);
      fprintf(fid,'pmean = %.16f;\n',pmean);
      fprintf(fid,'pvar = %.16f;\n',pse^2);
      fprintf(fid,'for j=1:length(x)\n');
      fprintf(fid,'    a = x(j);\n');
      fprintf(fid,'    f(j) = 1/(a-1)*(pmean*gamma(a)/gamma(a-1/2))^2-pmean^2-pvar;\n');
      fprintf(fid,'end\n');
      fclose(fid);
      [a,rc] = csolve(fname,5,[],1e-10,1000);
      if rc~=0, error('Search for iGam parameters failed, rc=%.0f',rc), end
      delete([fname,'.m'])
    end
    fname = sprintf('gamfcn%.0f',cputime*1e10);
    fid=fopen([fname,'.m'],'wt');
    fprintf(fid,'function f=%s(x)\n',fname);
    fprintf(fid,'f = gamma(x);\n');
    fclose(fid);
    b = (feval(fname,a-1/2)/pmean/feval(fname,a))^2;
    delete([fname,'.m'])
    Params(j).priormode = (1/b/(a+1/2))^(1/2);
    for jprc=1:nprc
      eval(sprintf('Params(%.0f).priorp%03.0f = gaminv(1-%f,%.16f,%.16f)^(-1/2);',...
        j,[1000,1]*Percentiles(jprc),a,b))
    end
    Params(j).priorparams = {a,b};
    Params(j).priorlpdfcmd = sprintf(...
      'log((%1$s>0)*(gampdf(%1$s^(-2),%2$.16f,%3$.16f)*2/%1$s^3))',Params(j).name,a,b);
    Params(j).priorpdfcmd = sprintf(...
      '(%%1$.16f>0)*(gampdf(%%1$.16f^(-2),%1$.16f,%2$.16f)*2/%%1$.16f^3)',a,b);
    Params(j).priorrndcmd = sprintf('gamrnd(%.16f,%.16f)^(-1/2)',a,b);
  elseif strcmp(Params(j).priordist,'IG2')
    pmean = Params(j).priormean;
    pse = Params(j).priorse;
    if pse==inf
      a = 2;
    else
      a = 2+pmean^2/pse^2;
    end
    b = 1/pmean/(a-1);
    Params(j).priormode = 1/b/(a+1);
    for jprc=1:nprc
      eval(sprintf('Params(%.0f).priorp%03.0f = gaminv(1-%f,%.16f,%.16f)^(-1);',...
        j,[1000,1]*Percentiles(jprc),a,b))
    end
    Params(j).priorparams = {a,b};
    Params(j).priorlpdfcmd = sprintf(...
      'log((%1$s>0)*(gampdf(%1$s^(-1),%2$.16f,%3$.16f)/%1$s^2))',Params(j).name,a,b);
    Params(j).priorpdfcmd = sprintf(...
      '(%%1$.16f>0)*(gampdf(%%1$.16f^(-1),%1$.16f,%2$.16f)/%%1$.16f^2)',a,b);
    Params(j).priorrndcmd = sprintf('gamrnd(%.16f,%.16f)^(-1)',a,b);
  elseif strcmp(Params(j).priordist,'C')
    pmean = Params(j).priormean;
    pse = 0;
    Params(j).priorse = 0;
    Params(j).priormode = pmean;
    for jprc=1:nprc
      eval(sprintf('Params(%.0f).priorp%03.0f = %.16f;',j,1000*Percentiles(jprc),pmean))
    end
    Params(j).priorparams = {pmean,pse};
    Params(j).priorlpdfcmd = sprintf('log(1*(%s==%.16f))',Params(j).name,pmean);
    Params(j).priorpdfcmd = sprintf('1*(%%.16f==%.16f)',pmean);
    Params(j).priorrndcmd = sprintf('%.16f',pmean);
  end
end

%% ------------------------------------------------------------------------

%% display results on screen
if ShowPriorTable
  fprintf('\nResults from Prior analysis:')
  fprintf('\n============================\n')
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
  disp(' ')
end

%% ------------------------------------------------------------------------

%% clean up
clear a b pmean pse jprc

%% Elapsed time
TimeElapsed.PriorAnalysis = toc-TimeElapsed.PriorAnalysis;

%% ------------------------------------------------------------------------

