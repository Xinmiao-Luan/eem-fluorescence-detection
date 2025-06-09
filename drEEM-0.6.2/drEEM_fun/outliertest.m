function varargout=outliertest(data,varargin)
%
% <strong>Syntax</strong>
%   dataout=<strong>outliertest</strong>(data,wavestep,myfactors,constraints,CC)
%
% <a href="matlab: doc outliertest">help for outliertest</a> <- click on the link

% Fit fast PARAFAC models to identify outlier variables and / or samples
% Models are fit with a convergence criterion of 1e-2 for speed. The best
% solution out of two initializtions is picked for robuestness.
%
% If models already exist, the function will ask whether model diagnostics
% should simply be plotted instead.
%
% USEAGE:
%    Test=outliertest(data,wavestep,myfactors,constraints,CC,Stepwise)
%
% INPUT:
%       data: data structure containing EEMs to be modelled in data.X
%
%   wavestep: (optional wavelength selection - [ExInt,EmInt])
%         wavestep=[] or [1,1]: use all the data (default)
%         wavestep=[1,2]: use every other Em wavelength 
%         wavestep=[2,2]: every other Ex and Em wavelength
%
%  myfactors: (optional, numbers of factors in models)
%                     []: construct models with 2-6 components (default)
%               e.g. 2:5: construct models with 2,3,4 and 5 factors
%
%constraints: (optional, default = 'unconstrained')
%        'nonnegativity': Apply non-negativity constraints on all modes
%        'unconstrained': No constraints
%        'unimodnonneg' - positive solutions, unimodal emission
%
%         CC: (optional) convergence criterion for PARAFAC model
%                     []: default value used is 1e-2
%
%  Stepwise : (optional, default = 'stepwise')
%             In stepwise mode,  models are run one at a time 
%             and diagnostic plots  created immediately after each  model. 
%             Press enter to continue to the next PARAFAC model. 
%             'stepwise': run model in stepwise mode
%              'at once': run models at once
%
% Examples
%   Test=outliertest(Xs) %Defaults: 2-6 factors, unconstrained models,
%   Test=outliertest(Xs,[2,2],2:7)
%   Test=outliertest(Xs,[1,2],5,'nonnegativity')
%   Test=outliertest(Xs,[],4:5,'unconstrained',1e-6,'at once')
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
%Copyright (C) 2020- Urban Wuensch
%Copyright (C) 2013- Kate Murphy krm@unsw.edu.au
%Copyright (C) 2008- Colin Stedmon
%
% $ Version 0.1.0 $ September 2013 $ First Release
%% Check inputs
try
    narginchk(1,6)
catch ME
    help outliertest
    pause(1)
    rethrow(ME)
end



EmInt=1;ExInt=1;
myfactors=2:6;
convgcrit=1e-2;
constraints='nonnegativity';
if nargin>1
    wavestep=varargin{1};
    if isempty(wavestep)
        wavestep=[1,1];
    end
    ExInt=wavestep(1);
    EmInt=wavestep(2);
    if nargin>2
        myfactors=varargin{2};
        if nargin>3
            NN=varargin{3};
            switch NN
                case{'nonnegativity','Nonnegativity'}
                    constraints='nonnegativity';
                case{'unconstrained','Unconstrained'}
                    constraints='unconstrained';
                case{'unimodnonneg','Unimodnonneg'}                   
                    constraints='unimodnonneg';
                otherwise
                    error('Allowed constraints are nonnegativity, unconstrained, and unimodnonneg');
            end
            if nargin>4
                convgcrit=varargin{4};
                if isempty(convgcrit)
                    convgcrit=1e-2;
                end
                if ~isnumeric(convgcrit)
                    error('Unrecognised input for convergence criterion')
                elseif convgcrit>0.01
                    error(['Unexpected convergence criterion: ' num2str(convgcrit)]);
                end
                if nargin>5
                    disp('''Stepwise'' is an obsolete option since drEEM v0.6.0')
                end
            end
        end
    end
end

Test=data;
Test.X=data.X(:,1:EmInt:end,1:ExInt:end);
if isfield(Test,'Xf')
    Test.Xf=data.Xf(:,1:EmInt:end,1:ExInt:end);
end
Test.Em=data.Em(1:EmInt:end);
Test.Ex=data.Ex(1:ExInt:end);
Test.nEm=length(Test.Em);
Test.nEx=length(Test.Ex);
Test.nSample=size(Test.X,1); %Test.nSample=length(Test.X); %BUG FIXED 2013/12/20

%% Welcome message
disp('  ')
disp('-----')
disp(' outliertest - Fit fast PARAFAC models to identify outlier variables and / or samples')
disp('-----')

%% Too much information?
if mean(diff(Test.Em))<=3
    disp(' ')
    disp('Emission resolution was automatically reduced')
    Test=subdataset(Test,[],1:2:Test.nEm,[]);
end
disp(' ')
disp(' ')
disp(' ')

fthere=zeros(1,100);
for n=1:100
    if isfield(Test,['Model',num2str(n)])
        fthere(n)=n;
    else
        fthere(n)=0;
    end
end
fthere(fthere==0)=[];
%% Modeling
fitmodels=true;
if ~isempty(fthere)
    disp('Found existing models in ''data'' input. Do you want to plot these instead of calculating new models?')
    decision=input('y/n: ','s');
    if strcmp(decision,'y')
        fitmodels=false;
        dataout=data;
        f=fthere;
    elseif strcmp(decision,'n')
        disp('New outliertest-models will be fit to the EEMs')
    end
end

if fitmodels
    dataout=randinitanal(Test,myfactors,'starts',2,'convgcrit',convgcrit,'maxit',1000,'init','svd','mode','outliertest','constraints',constraints);
    close all
end

%% Show results

f_there=zeros(1,100);
for n=1:100
    if isfield(dataout,['Model',num2str(n)])&&size(dataout.(['Model',num2str(n)]),2)==3
        if ~isempty(dataout.(['Model',num2str(n)]){2})
            f_there(n)=n;
        else
            f_there(n)=0;
        end
    end
end
f_there(f_there==0)=[];



dims=[2 3;1 3;1 2];
lev=cell(numel(f_there),3);
ssq=cell(numel(f_there),3);
modall=cell(numel(f_there),1);
core=nan(numel(f_there),1);
pere=nan(numel(f_there),1);

for ii=1:numel(f_there)
    mod=dataout.(['Model',num2str(f_there(ii))]);
    if isempty(mod)
        modall{ii}={nan(data.nSample,f_there(ii)) nan(data.nEm,f_there(ii)) nan(data.nEx,f_there(ii))};
        continue
    end
    res=dataout.X-nmodel(mod);
    for n=1:3
        lev{ii,n}=diag(mod{n}*(mod{n}'*mod{n})^-1*mod{n}');
        ssq{ii,n}=rcvec(squeeze(nansum(nansum(res.^2,dims(n,1)),dims(n,2))),'column');
    end
    modall{ii}=mod;
    core(ii)=dataout.(['Model',num2str(f_there(ii)),'core']);
    pere(ii)=dataout.(['Model',num2str(f_there(ii)),'percentexpl']);
end

%%
% figHandles = get(0,'Children');
% if any(contains({figHandles.Name},'Outliertest:'))
%     close(figHandles(contains({figHandles.Name},'Outliertest:')))
% end
fhandle=dreemfig;
w=.4;h=.55;l=(1-w)/2;b=(1-h)/2;
set(fhandle,...
    'units','normalized','outerposition',[l b w h],...
    'Name','Outliertest: Model diagnostics');
ax = gobjects(4*3-1,1);
for n=1:4*3-1
    ax(n)=subplot(3,4,n);
end

uicontrol(fhandle,'Style','popupmenu','String',...
    cellstr(strcat(num2str(f_there'),repmat('-comp. model',numel(f_there),1))),...
    'Units','normalized',...
    'position',[0.7470,0.2477,0.1871,0.0475],...
    'Callback',{@plotdata,dataout,ax,lev,ssq,modall},...
    'Value',1);

plot(ax(4),f_there,pere,'o--','Color',[0.3467    0.5360    0.6907],...
    'MarkerSize',5,'MarkerFaceColor',[0.3467    0.5360    0.6907],'MarkerEdgeColor','k',...
    'LineWidth',1.4,'DisplayName','perexpl');

plot(ax(8),f_there,core,'o--','Color',[0.9153    0.2816    0.2878],...
    'MarkerSize',5,'MarkerFaceColor',[0.9153    0.2816    0.2878],'MarkerEdgeColor','k',...
    'LineWidth',1.4,'DisplayName','corecon');

title(ax(4),{'Other metrics',' ','% expl.'})
ylim(ax(4),[min(pere) 100])
xlabel(ax(4),'No. components')
ylabel(ax(4),'% expl.')

title(ax(8),{'Core consistency'})
ylim(ax(8),[-10 100])
xlabel(ax(8),'No. components')
ylabel(ax(8),'% Core consistency')

plotdata([],[],dataout,ax,lev,ssq,modall)

dcm_obj = datacursormode(fhandle);
set(dcm_obj,'UpdateFcn',{@dispA,dataout});
datacursormode on
dreemfig(fhandle);
if nargout==1
    varargout{1}=dataout;
end

end


%%
function [vout] = rcvec(v,rc)
% Make row or column vector
% v: vector
% rc: either 'row' ([1:5])or 'column' ([1:5]')
sz=size(v);
if ~any(sz==1)
    error('Input is not a vector')
end

switch rc
    case 'row'
        if ~(sz(1)<sz(2))
            vout=v';
        else
            vout=v;
        end
    case 'column'
        if ~(sz(1)>sz(2))
            vout=v';
        else
            vout=v;
        end
    otherwise
            error('Input ''rc'' not recognized. Options are: ''row'' and ''column''.')
end


end


%%
function plotdata(src,event,dataout,ax,lev,ssq,mod)

try
    f= src.Value;
catch
    f=1;
end
% model
try
    plot(ax(1),dataout.i,[mod{f}{1}],'DisplayName','load1');ax(1).Box='off';
    plot(ax(5),dataout.Em,[mod{f}{2}],'DisplayName','load2');ax(5).Box='off';
    plot(ax(9),dataout.Ex,[mod{f}{3}],'DisplayName','load3');ax(9).Box='off';
catch
    disp('No valid model calculated.')
    cla(ax(1),'reset')
    cla(ax(5),'reset')
    cla(ax(9),'reset')
end
for n=[1 5 9],axis(ax(n),'tight'),end

% Leverage
try
    h(1)=scatter(ax(2),dataout.i,[lev{f,1}],10,'filled','DisplayName','lev1');
    h(2)=scatter(ax(6),dataout.Em,[lev{f,2}],10,'filled','DisplayName','lev2');
    h(3)=scatter(ax(10),dataout.Ex,[lev{f,3}],10,'filled','DisplayName','lev3');
catch
    disp('No valid model calculated.')
    cla(ax(2),'reset')
    cla(ax(6),'reset')
    cla(ax(10),'reset')
end

% Residual
try
    h(4)=scatter(ax(3),dataout.i,[ssq{f,1}],10,'filled','DisplayName','res1');
    h(5)=plot(ax(7),dataout.Em,[ssq{f,2}],'k','DisplayName','res2');ax(5).Box='off';
    h(6)=plot(ax(11),dataout.Ex,[ssq{f,3}],'k','DisplayName','res3');ax(8).Box='off';
catch
    disp('No valid model calculated.')
    cla(ax(3),'reset')
    cla(ax(7),'reset')
    cla(ax(11),'reset')
end
for n=1:4
    alpha(ax(n),.5)
    h(n).MarkerEdgeColor='k';
    h(n).MarkerEdgeAlpha=.5;
end


title(ax(1),{'Loadings',' ','Sample'})
title(ax(5),{'Em'})
title(ax(9),{'Ex'})

title(ax(2),{'Leverage',' ','Sample'})
title(ax(6),{'Em'})
title(ax(10),{'Ex'})

title(ax(3),{'Residuals',' ','Sample'})
title(ax(7),{'Em'})
title(ax(11),{'Ex'})

for n=[5:7 9:11 ],xlabel(ax(n),'Wavelength (nm)'),end
for n=1:3,xlabel(ax(n),'Sample no.'),end
for n=[1 5 9],ylabel(ax(n),'Loadings'),end
for n=[2 6 10],ylabel(ax(n),'Lev.'),end
for n=[3 7 11],ylabel(ax(n),'Sum Sq. Res.'),end

legentries=strcat(repmat('C ',size(mod{f}{1},2),1),num2str((1:size(mod{f}{1},2))'));
legend1=legend(ax(1),legentries);
pos=legend1.Position;
legend1.Position=[0.8025 0.0523 pos(3) pos(4)];
try legend1.ItemTokenSize=[30/2 18/2];end

dreemfig(gcf);

end


function txt = dispA(~,event_obj,data)

pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
try 
    flist=data.filelist{I};
catch
    flist='no filename';
end
seccharTarget=str2double(event_obj.Target.DisplayName(2));
if ~isnan(seccharTarget)
    % New R2020b behaviour fix. Legend has changed DisplayName.
    % This affects only load1. Look for number in the second character.
    txt = {['data.i: ',num2str(data.i(I))],...
               ['i: ',num2str(I)],...
               ['fname: ',flist],...
               ['Score: ',num2str(pos(2))]};
else
    switch event_obj.Target.DisplayName
        case 'load1'
            txt = {['data.i: ',num2str(data.i(I))],...
                   ['i: ',num2str(I)],...
                   ['fname: ',flist],...
                   ['Score: ',num2str(pos(2))]};
        case 'load2'
            txt = {['i: ',num2str(I)],...
                   ['Em: ',num2str(data.Em(I)),'nm'],...
                   ['load: ',num2str(pos(2))]};

        case 'load3'
            txt = {['i: ',num2str(I)],...
                   ['Ex: ',num2str(data.Ex(I)),'nm'],...
                   ['load: ',num2str(pos(2))]};        
        case 'lev1'
            txt = {['data.i: ',num2str(data.i(I))],...
                   ['i: ',num2str(I)],...
                   ['fname: ',flist],...
                   ['Lev: ',num2str(pos(2))]};
        case 'lev2'
            txt = {['current index: ',num2str(I)],...
                   ['Em: ',num2str(data.Em(I)),'nm'],...
                   ['Lev: ',num2str(pos(2))]};
        case 'lev3'
            txt = {['i: ',num2str(I)],...
                   ['Ex: ',num2str(data.Ex(I)),'nm'],...
                   ['Lev: ',num2str(pos(2))]};
        case 'res1'
            txt = {['data.i: ',num2str(data.i(I))],...
                   ['i: ',num2str(I)],...
                   ['fname: ',flist],...
                   ['Res: ',num2str(pos(2))]};
        case 'res2'
            txt = {['current index: ',num2str(I)],...
                   ['Em: ',num2str(data.Em(I)),'nm'],...
                   ['Res: ',num2str(pos(2))]};
        case 'res3'
            txt = {['i: ',num2str(I)],...
                   ['Ex: ',num2str(data.Ex(I)),'nm'],...
                   ['Res: ',num2str(pos(2))]};
        case 'perexpl'
            txt = {['Components: ',num2str(event_obj.Target.XData(I))],...
                   ['%expl.: ',num2str(pos(2))]};
        case 'xcorecon'
            txt = {['Components: ',num2str(event_obj.Target.XData(I))],...
                   ['%core con.: ',num2str(pos(2))]};

        otherwise
            warning('''DisplayName'' not identified')
    end
end


end