function loadingsandleverages(data,f,varargin)
%
% <strong>Syntax</strong>
%   <strong>loadingsandleverages</strong>(data,f,it,fignames)
%
% <a href="matlab: doc loadingsandleverages">help for loadingsandleverages</a> <- click on the link

% Visualise the PARAFAC model loadings and leverages. Leverages show the
% impact of each sample and wavelength (excitation, emission) on the model.
%
% USEAGE:   
%       loadingsandleverages(data,f,it,fignames)
%
% INPUTS: 
%      data: dataset structure containing PARAFAC model results
%         f: Number of components in the model to be plotted, 
%            e.g. 6 to plot the 6-component model in data.Model6.
%       run: (optional) 
%            Run number of model to be plotted, when data contains
%            multiple runs of the same model, as from the output
%            of randinitanal.
%             []: the main model will be plotted (default)
%              n: plot run number x (the model plotted will be
%                data.Modeln_runx, e.g. data.Model6_run2)
%  fignames: (optional) 
%             []: by default the figure is called 
%               "Model f - scores, leverages and loading"
%       'mytext': prefix the figure name by the word(s) in mytext
%
% EXAMPLES:
%    loadingsandleverages(Test1,6)
%    loadingsandleverages(Test1,6,5)
%    loadingsandleverages(S.Split(1),6,5,'Split 1, Run 5')
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% loadingsandleverages: Copyright (C) 2019 Kathleen R. Murphy & Urban J. Wuensch
% Chalmers University of Technology
% Sven Hultins Gata 6
% 41296 Gothenburg
% Sweden
% murphyk@chalmers.se & wuensch@chalmers.se
%
% $ Version 0.2.0 $ March 2019     $ plotting reimplemented with cursor info
% $ Version 0.1.0 $ September 2013 $ First Release
%%%%%%

%Initialise
narginchk(2,4)
fignames=[];
R=[];

data.iseq=(1:data.nSample)';
data.iEm=(1:data.nEm)';
data.iEx=(1:data.nEx)';

if length(f)>1
    error('Specify one value of ''f'' at a time');
else
    modelf=['Model' num2str(f)];
end
if nargin>2
    R=varargin{1};
    if nargin>3
        fignames=varargin{2};
    end
end

if ~isempty(R)
    if ~isnumeric(R)
        error('Input for variable ''run'' is not understood.')
    else
        modelf=[modelf '_run' int2str(R)];
    end
end

if ~isfield(data,modelf)
   disp(modelf)
   disp(data)
   error('loadingsandleverages:fields',...
       'The dataset does not contain a model with the specified number of factors') 
end

h=dreemfig;
set(h,'units','normalized','pos',[0.2 0.25 0.6 0.4])
if isempty(fignames)
    set(h,'Name',[modelf ' - scores,leverage and loadings']);
else
    set(h,'Name',[fignames ' ' modelf ' - scores,leverage and loadings']);
end
%%
% Data extraction
factors=getfield(data,{1,1},modelf);
lev=cell(1,numel(factors));
for n=1:numel(factors)
    lev{n}=diag(factors{n}*(factors{n}'*factors{n})^-1*factors{n}');
end
xname={'iseq','Em','Ex'};      % x-axis scale

xaxl={'Sample # (data.i)','Em. (nm)','Ex. (nm)'};
yaxl={'Scores','Loadings','Loadings','Leverage'};
taxl={'Sample','Emission','Excitation'};
ax = gobjects(1,3);
hax= cell(3,1);
for n=1:3
    ax(n)=subplot(2,3,n);
    x=data.(xname{n});
    [x,idx]=sort(x);
    hax{n,1}=line(x,factors{n}(idx,:),'Marker','.','MarkerSize',7,'LineWidth',0.5,'DisplayName',['lod',num2str(n)]);
    xlabel(xaxl{n})
    ylabel(yaxl{n})
    axis tight
    box on
    title(taxl{n})
    
    ax(n+3)=subplot(2,3,n+3);
    hax{n,2}=line(x,lev{n}(idx,:),'Marker','+','Color','k',...
        'LineStyle','none','DisplayName',['lev',num2str(n)]);
    axis tight
    xlabel(xaxl{n})
    ylabel(yaxl{4})

    box on
end
dreemfig(h);
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@dispA,data});
datacursormode on

axes(ax(3))
h=leg(f,lines(f));
legend(h,num2str((1:f)'),'Location','Best')
disp(' ')
disp('Cursor active. Click on any leverage to obtain details (sample index wavelength).')
disp(' ')


end


function txt = dispA(~,event_obj,data)

pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
try 
    flist=data.filelist{I};
catch
    flist='no filename';
end
switch event_obj.Target.DisplayName
    case 'lod1'
        txt = {['i: ',num2str(data.i(I))],...
               ['current index: ',num2str(I)],...
               ['filename: ',flist],...
               ['Score: ',num2str(pos(2))]};
    case 'lod2'
        txt = {['current index: ',num2str(I)],...
               ['Em: ',num2str(data.Em(I)),'nm'],...
               ['loading: ',num2str(pos(2))]};

    case 'lod3'
        txt = {['current index: ',num2str(I)],...
               ['Ex: ',num2str(data.Ex(I)),'nm'],...
               ['loading: ',num2str(pos(2))]};        
    case 'lev1'
        txt = {['i: ',num2str(data.i(I))],...
               ['current index: ',num2str(I)],...
               ['filename: ',flist],...
               ['Leverage: ',num2str(pos(2))]};

    case 'lev2'
        txt = {['current index: ',num2str(I)],...
               ['Em: ',num2str(data.Em(I)),'nm'],...
               ['Leverage: ',num2str(pos(2))]};
    case 'lev3'
        txt = {['current index: ',num2str(I)],...
               ['Ex: ',num2str(data.Ex(I)),'nm'],...
               ['Leverage: ',num2str(pos(2))]};
    otherwise
        warning('''DisplayName'' not identified')
end


end

function [ h ] = leg( numPlots,col )
h = gobjects(numPlots, 1);
for n=1:numPlots
    hold on;
    h(n) = line(NaN,NaN,'Marker','o','MarkerSize',5,'MarkerFaceColor',col(n,:),'MarkerEdgeColor','k','LineStyle','none');
end
end