function PostPlotCond(Par2Show,Params,FileName,Options)

% PostPlotCond
%
% Generates a surface plot for the value of the negative of the
% log-posterior value given a benchmark vector of parameters, conditional
% on a grid of values for the two parameters listed in Par2Show.
%
% ..............................................................................
%
% Created: September 12, 2011 by Vasco Curdia
% Copyright (2011) by Vasco Curdia

%% -----------------------------------------------------------------------------

%% defaults
nGrid = 20;
BaseMoment = 'postmode';
xBoundsMoment = 'postse';
NewFig = 1;
BoundsFactor = 1;
UseContour = 0;

%% Load new options
if exist('Options','var')
  OptionList = fieldnames(Options);
  for j=1:length(OptionList)
    eval(sprintf('%1$s = Options.%1$s;',OptionList{j}))
  end
end

%% Preliminary
if ~exist('xBase','var'),xBase = [Params(:).(BaseMoment)]';end
ParamNames = {Params(:).name};
np = length(Params);
idxPar2Show = find(ismember(ParamNames,Par2Show));
if ~exist('xBounds','var')
  xBounds = [Params(idxPar2Show).(xBoundsMoment)]'*BoundsFactor;
  xBounds = [xBase(idxPar2Show)-xBounds,xBase(idxPar2Show)+xBounds];
end
FileNamePost = FileName.Post;

%% Prepare Grid
SinglePar = (length(idxPar2Show)==1);
if SinglePar
  idxPar2Show = idxPar2Show(1);
  xGrid = zeros(nGrid,1);
  xRange = xBounds(1,2)-xBounds(1,1);
  xGrid(:,1) = xBounds(1,1):(xRange/(nGrid-1)):xBounds(1,2);
  XGrid = xGrid';
  nSim = nGrid;
  Xd = repmat(xBase,1,nGrid);
  for jC=1:nGrid
    Xd(idxPar2Show,jC) = XGrid(:,jC);
  end
  xd = Xd;
  postd = zeros(nGrid,1);
else
  xGrid = zeros(nGrid,2);
  for j=1:2
    xRange = xBounds(j,2)-xBounds(j,1);
    xGrid(:,j) = xBounds(j,1):(xRange/(nGrid-1)):xBounds(j,2);
  end
  XGrid = zeros(2,nGrid,nGrid);
  [XGrid(1,:,:),XGrid(2,:,:)] = meshgrid(xGrid(:,1),xGrid(:,2));
  nSim = nGrid*nGrid;
  Xd = repmat(xBase,[1,nGrid,nGrid]);
  for jR=1:nGrid
    for jC=1:nGrid
      Xd(idxPar2Show,jR,jC) = XGrid(:,jR,jC);
    end
  end
  xd = reshape(Xd,np,nSim);
  postd = zeros(nGrid,nGrid);
end

%% Evaluate Post
parfor j=1:nSim
  postd(j) = feval(FileNamePost,xd(:,j));
end

%% Show plot
if NewFig
  figure('Name',['logPost(',Par2Show{1},',',Par2Show{2},')'])
end
if SinglePar
  plot(xGrid(:,1),postd)
  hold on
  plot(xBase(idxPar2Show)*[1,1],[min(postd),max(postd)],':k')
%   plot(Params(idxPar2Show).priormean*[1,1],[min(postd),max(postd)],':r')
%   plot(Params(idxPar2Show).priormode*[1,1],[min(postd),max(postd)],':g')
  hold off
  axis tight
  xlabel(Params(idxPar2Show(1)).prettyname)
else
  if UseContour
    contour(xGrid(:,1),xGrid(:,2),postd)
  else
    surf(xGrid(:,1),xGrid(:,2),postd)
  end
  axis tight
  xlabel(Params(idxPar2Show(1)).prettyname)
  ylabel(Params(idxPar2Show(2)).prettyname)
end

%% -----------------------------------------------------------------------------


