function [best,all,details] = randinitanal(data,f,varargin)
%
% <strong>Syntax</strong>
%   [best,all,details] = <strong>randinitanal</strong>(data,f,Name,Value)
%   [best,all,details] = <strong>randinitanal</strong>(data,f,it,constraints,convgcrit)
%
% <a href="matlab: doc randinitanal">help for randinitanal</a> <- click on the link

% Run PARAFAC multiple times to identify the least squares solution
%
%USEAGE: 
%       [best,all,details] = randinitanal(data,f,varargin) OR
%       [best,all,details] = randinitanal(data,f,it,constraints,convgcrit)
%
%INPUTS:     data: data structure containing EEMs to be modelled in data.X.
%               f: Number of components in the model to be fitted, e.g. 2:5
%        varargin: Optional input arguments (see below)
%
%
%OPTIONAL INPUTS: (option name followed by value, e.g. "'starts',5")
%        constraints: Constrain the model loadings.
%                     'unconstrained', 'nonnegativity', or 'unimodnonneg' (positive solutions, unimodal emission)
%                     Default: 'nonnegativity'.
%             starts: Number of random starts per model. Default: 40
%          convgcrit: Convergence criterion. Default: 1e-6
%              maxit: Maximum number of iterations. Default: 2500
%               tbox: PARAFAC engine. 'nway' or 'pls'. Default: 'nway'
%   NOTE: Optional inputs can also be supplied in the sequence ->
%            it,constraints,convgcrit
%
%
%  OUTPUTS:
%               best: Least squares model chosen from all model runs.
%                all: Data structure containing all iterated models.
%            details: Table, summary of convergence results for each run.
%
% The function automatically tests if it can benefit from multiple CPUs
%
%  NOTES:
%      If multiple CPUs are used, DO NOT use Ctrl+C to cancel this function
%      Use the cancel button on the waitbar instead.
%
%
%Examples
%   model=randinitanal(data,5:8) 
%   [bestmodel,allmodel,moddetail]=randinitanal(data,5,'starts',10,'constraints','nonnegativity','convgcrit',1e-10,'maxit',5000)
%   [bestmodel,allmodel,moddetail]=randinitanal(data,2,10,'nonnegativity',1e-2)
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013,
%     DOI:10.1039/c3ay41160e.
%
% randinitanal: Copyright (C) 2019 Urban J. Wuensch
% Chalmers University of Technology
% Sven Hultins Gata 6
% 41296 Gothenburg
% Sweden
% wuensch@chalmers.se
% $ Version 0.1.0 $ August 2019 $ First Release

%% Input argument handling
if nargin==0
    help randinitanal
    return
end
% Plot or calculate?

if numel(varargin)==0&&~exist('f','var')
    plotmodels=true;
else
    plotmodels=false;
end
% Decide whether or not to use input parser
parseargs={'constraints','starts','convgcrit','maxit','tbox','init','mode'};
if numel(varargin)>0
    if any(contains(cellstr(string(varargin{1})),parseargs)) % Name,Value-scheme is used
        parser=true;
    else
        parser=false;
    end
else
    parser=true;
end
if parser
    try
        params = inputParser;
        params.addParameter('constraints', 'nonnegativity', @checkconstraints);
        params.addParameter('starts', 40, @isnumeric);
        params.addParameter('convgcrit', 1e-6, @isnumeric);
        params.addParameter('maxit', 2500, @isnumeric);
        params.addParameter('tbox', 'nway', @checktoolbox);
        params.addParameter('init', 'random', @checkinit);
        params.addParameter('mode', 'default', @ischar);
        params.addParameter('consoleoutput', 'all', @checkconsoleoutput);
        params.addParameter('initvalues', {}, @iscell);
        params.parse(varargin{:});
    catch ME
        help randinitanal
        error('Function input not understood. Please see text above for help on input specification.')
    end
    constraints = params.Results.constraints;
    convgcrit  = params.Results.convgcrit;
    maxit  = params.Results.maxit;
    tbox = params.Results.tbox;
    starts = params.Results.starts;
    initstyle = params.Results.init;
    initvalues = params.Results.initvalues;
    fmode = params.Results.mode;
    consoleoutput = params.Results.consoleoutput;
else %old-style inputs are supplied
    % Try to parse these inputs manually
    numel(varargin)
    if numel(varargin)>=1
        starts = varargin{1};
        if isempty(starts);starts=40;end
    else
        starts=40;
    end
    if numel(varargin)>=2
        constraints = varargin{2};
        if isempty(constraints);constraints='nonnegativity';end
    else
        constraints = 'nonnegativity';
    end
    if numel(varargin)>=3
        convgcrit = varargin{3};
        if isempty(convgcrit);convgcrit=1e-6;end
    else
        convgcrit = 1e-6;
    end
    maxit  = 2500;
    tbox = 'nway';
    initstyle = 'random';
    fmode = 'default';
    consoleoutput='all';
end
% Check that the inputs (if supplied) were correct
if ~isnumeric(starts)
    error('Input to ''it'' must be numeric')
end
if ~isnumeric(convgcrit)||convgcrit>1e-2||convgcrit<1e-15
    error('Input to ''convgcrit'' is either not numeric, too small, or too big.')
end
if ~any([contains(constraints,'nonnegativity') contains(constraints,'unconstrained') contains(constraints,'unimodnonneg')])
    error('Input to ''constraints'' must be ''unconstrained'', ''nonnegativity'', or ''unimodnonneg''')
end

if strcmp(initstyle,'given')
    if numel(initvalues)~=numel(f)
        error('''initvalues'' must have as many elements as input to ''f''')
    end
    for ii=1:numel(f)
        if ~sum(size(initvalues{ii}{1})==[data.nSample,f(ii)])==2
            error('''initvalues'' must have dimensions compatible with data.X')
        end
        if ~sum(size(initvalues{ii}{2})==[data.nEm,f(ii)])==2
            error('''initvalues'' must have dimensions compatible with data.X')
        end
        if ~sum(size(initvalues{ii}{3})==[data.nEx,f(ii)])==2
            error('''initvalues'' must have dimensions compatible with data.X')
        end
    end
end

%% Plot if f is not provided
if plotmodels
    ftoplot=zeros(1,100);
    cnt=1;
    for n=1:100
        if isfield(data,['Model',num2str(n)])&&size(data.(['Model',num2str(n)]),2)==3
            ftoplot(n)=n;
            mod{cnt}=data.(['Model',num2str(n)]);
            sz(cnt,1)=size(mod{cnt}{2},1);
            sz(cnt,2)=size(mod{cnt}{3},1);
            cnt=cnt+1;
            
        else
            ftoplot(n)=0;
        end

    end
    ftoplot(ftoplot==0)=[];
    if exist('mod','var')
        c1=sum( sz(:,1)==data.nEm )==numel(mod);
        c2=sum( sz(:,2)==data.nEx )==numel(mod);
        if c1&&c2
            plotfacs(mod,ftoplot,[],data.Em,data.Ex)
        else
            warning('Cannot plot the existing models. Em and / or Ex does not match the model')
        end
    else
        warning('Did not find any models to plot')
    end
    return
end
%% Setup
funmode=parallelcomp(consoleoutput);
[oldp,newp]=tboxinit(tbox);
opt = setoptions(tbox,constraints,convgcrit,maxit,initstyle);

fac=f;
numstarts=starts*numel(fac);
facCalls=reshape(repmat(fac,starts,1),numstarts,1);
if strcmp(initstyle,'given')
    ivalsCalls=reshape(repmat(initvalues,starts,1),numstarts,1);
else
    ivalsCalls=cell(1,numstarts);
end
%% Welcome messages
if strcmp(fmode,'default')&&~strcmp(consoleoutput,'none')
    disp('  ')
    disp('-----')
    disp(' randinitanal.m - Random initializations of PARAFAC models')
    disp('-----')
    switch tbox
        case 'pls'
            disp(['PARAFAC engine:                  ','PLS_toolbox'])
        case 'nway'
            disp(['PARAFAC engine:                  ','N-way toolbox'])
    end
    switch funmode
        case 'sequential'
            disp(['Parallelization:                 ','no'])
        case 'parallel'
            disp(['Parallelization:                 ','yes'])
    end
    disp(['Models with components:          [',num2str(f),']'])
    disp(['Number of random starts:         ',num2str(starts)])
    disp(['Convergence criterion:           ',num2str(convgcrit)])
    disp(['Maximum number of iterations:    ',num2str(maxit)])
    disp(['Constraint:                      ',constraints])
    disp(['% missing numbers:               ',num2str(round( (sum(isnan(data.X(:)))./numel(data.X))*100 ,3))])
    disp('-----')
    disp('  ')
end
%% Submit the models
ttot=tic;
switch funmode
    case 'sequential'
         modout = repelem(struct('model',{},'err',{},'iterations',{}),1,numstarts); % Preallocate output
         h=waitbar(0,'PARAFAC models are being calculated');
        for i=1:numstarts
            modout(i)=dreemparafac(data.X,facCalls(i),opt,tbox,ivalsCalls{i}); 
            try
            h=waitbar(i/numstarts,h,['PARAFAC models are being calculated [',num2str(i),'/',num2str(numstarts),']']);
            catch
                disp('Waitbar was closed, but function will keep running.')
            end
        end
        close(h)
        Model={modout.model}';
        Iter=arrayfun(@(x) x.iterations,modout,'UniformOutput',false)';
        Err=arrayfun(@(x) x.err,modout,'UniformOutput',false)';
        for ii=1:numel(Model);corecon{ii,1}=corcond(data.X,Model{ii},[],0);end

    case 'parallel'
        modout(1:numstarts) = parallel.FevalFuture; % Preallocate output
        
        for i=1:numstarts
            modout(i)=parfeval(@dreemparafac,1,data.X,facCalls(i),opt,tbox,ivalsCalls{i});
        end
        [Model,Iter,Err,corecon,ttc]=trackprogress(modout,numstarts,facCalls,tbox,data.X,data.Em,data.Ex,consoleoutput);
end

restoreoldpath(oldp,newp); % In case matlabpath was changed, restore it.
%% Retreive results

[best,all,details] = analyzeresults(data,fac,starts,maxit,facCalls,convgcrit,constraints,Model,Iter,Err,corecon,initstyle);

switch funmode
    case 'parallel'
    ttot=toc(ttot);
    ttc=sum(seconds(ttc));
    if ~strcmp(consoleoutput,'none')
        disp(' ')
        disp('Finished.')
        disp(['Time elaped: ',num2str(round(ttot./60,2)),'min. Parallel computing speedup: x',num2str(round(ttc/ttot,2))])
        disp(' ');
    end
 case 'sequential'
    ttot=toc(ttot);
    if ~strcmp(consoleoutput,'none')
        disp(' ')
        disp('Finished.')
        disp(['Time elaped: ',num2str(round(ttot./60,2)),'min.'])
    end
end


end

%% Internal functions used above (shared with splitanalysis)
%%
function checkconstraints(finput)
if ~contains(finput,{'unconstrained', 'nonnegativity', 'unimodnonneg'})
    error('Input to ''constraints'' must be ''unconstrained'', ''nonnegativity'', or ''unimodnonneg''')
end
end

%%
function checkinit(finput)
if ~contains(finput,{'svd', 'random','multiplesmall','given'})
    error('Input to ''init'' must be ''svd'', ''given'' or ''random''.')
end
end
%%
function checkconsoleoutput(finput)
if ~contains(finput,{'all', 'minimum','none'})
    error('Input to ''consoleoutput'' must be ''all'', ''minimum'' or ''none''.')
end
end

%%
function funmode=parallelcomp(consoleoutput)
test=ver;
funmode='sequential';
if any(contains({test.Name},'Parallel'))
    funmode='parallel';
    try
        initppool(consoleoutput)
    catch
        funmode='sequential';
    end
end
end

%%
function checktoolbox(finput)
if ~contains(finput,{'nway', 'pls'})
    error('Input to ''tbox'' must be ''nway'', or ''pls''. Other software products are currently not supported.')
end
end

%%
function initppool(consoleoutput)
warning off
try
    poolsize=feature('NumCores');
    p = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(p)
        parpool('local',poolsize);
    elseif p.NumWorkers~=poolsize
        if ~strcmp(consoleoutput,'none')
            disp('Found existing parpool with wrong number of workers.')
            disp(['Will now create pool with ',num2str(poolsize),' Workers.'])
        end
        delete(p);
        parpool('local',poolsize);
    else
        if ~strcmp(consoleoutput,'none')
            disp(['Existing parallel pool of ',num2str(p.NumWorkers),' workers found and used...'])
        end
    end
catch ME
    rethrow(ME)
end
warning on
end

%%
function [oldp,newp]=tboxinit(tbox)
warning off
oldp = path;
mlpath = path;mlpath=textscan(mlpath,'%s','delimiter',pathsep);mlpath=mlpath{:};
nway=contains(mlpath,'drEEM');
plsp=contains(mlpath,'PLS');

switch tbox
    case 'nway'
        if find(nway,1,'first') > find(plsp,1,'first')
            rmpath(mlpath{nway|plsp});
            addpath(mlpath{plsp});
            addpath(mlpath{nway});
        end
    case 'pls'
        if find(nway,1,'first') < find(plsp,1,'first')
            rmpath(mlpath{nway|plsp});
            addpath(mlpath{nway});
            addpath(mlpath{plsp});
        end
        clearvars pls
        try
            pls test
            disp('Testing PLS_toolbox'),evridebug
        catch
            error('PLS_toolbox either not installed or not functional.')
        end
    otherwise
        error('That''s a toolbox I don''t know about.')
end
newp = path;
warning on
end

%%
function opt = setoptions(tbox,constraints,convgcrit,maxIt,initstyle)
if contains(tbox,'pls')
    opt=parafac('options');
    opt.plots='none';
    opt.waitbar='off';
    opt.coreconsist='off';
    opt.display='on';

    if contains(constraints,'nonnegativity')
        % Set default: All dims nonnegative
        cdim=[1:3];
        % Check for custom input (deletes default)
        if ~isempty(erase(constraints,'nonnegativity'))
            t=(erase(constraints,'nonnegativity'));
            cdim=[];
            for n=1:numel(t)
                cdim=[cdim str2double(t(n))];
            end
            if numel(cdim)>3||any(cdim>3)
                error('numeric input to nonnegativity not understood')
            end
        end
        for i=cdim;opt.constraints{i}.type='nonnegativity';end
    elseif contains(constraints,'unimodality')
        for i=[1 3];opt.constraints{i}.type='nonnegativity';end
        opt.constraints{2}.type='unimodnonneg';
    end
    if contains(initstyle,'svd')
        opt.init=2;
    elseif contains(initstyle,'random')
        opt.init=3;
    elseif contains(initstyle,'given')
        opt.init=0;
    end
    opt.stopcriteria.iterations=maxIt;
    opt.stopcriteria.relativechange=convgcrit;
elseif contains(tbox,'nway')
    if contains(constraints,'nonnegativity')
        % Set default: All dims nonnegative
        cdim=[1:3];
        % Check for custom input (deletes default)
        if ~isempty(erase(constraints,'nonnegativity'))
            t=(erase(constraints,'nonnegativity'));
            cdim=[];
            for n=1:numel(t)
                cdim=[cdim str2double(t(n))];
            end
            if numel(cdim)>3||any(cdim>3)
                error('numeric input to nonnegativity not understood')
            end
        end
        for i=cdim;opt.constraints{i}.type=2;end
        for i=setdiff(1:3,cdim);opt.constraints{i}.type=0;end
        
    elseif contains(constraints,'unimodnonneg')
        for i=[1 3];opt.constraints{i}.type=2;end
        opt.constraints{2}.type=3;
    elseif contains(constraints,'unconstrained')
        for i=1:3;opt.constraints{i}.type=0;end
    end
    if contains(initstyle,'svd')
        opt.init=1;
    elseif contains(initstyle,'random')
        opt.init=2;
    elseif contains(initstyle,'multiplesmall')
        opt.init=10;
    elseif contains(initstyle,'given')
        opt.init=0;
    end
    opt.stopcriteria.iterations=maxIt;
    opt.stopcriteria.relativechange=convgcrit;
end

end

%%
function [out] = dreemparafac(tensor,f,opt,tbox,initvalues)
% Start the PARAFAC model. PLS_toolbox reads a pref-file, which is tricky if
% that's done in parallel. This here is to catch errors related to that.
if opt.init==0
    for n=1:3
        randvals=rand(size(initvalues{n},1),size(initvalues{n},2));
        randvals=randvals./vecnorm(randvals,2);
        idx=isnan(initvalues{n});
        initvalues{n}(idx)=randvals(idx);
    end
end
% plotfac(initvalues)
switch tbox
    case 'pls'
        tries=0;
        while tries<100
            try
                disp(datetime)
                if opt.init~=0
                    outlocal=parafac(tensor,f,opt);
                elseif opt.init==0
                    outlocal=parafac(tensor,initvalues,opt);
                    disp('Initialization with given values')
                end
                tries=102;
            catch ME
                tries = tries+1;
                pause(rand(1)/16);
                if tries>100
                    disp('PLS PARAFAC could not be initiated due to error')
                    rethrow(ME)
                end
            end
        end
        out.model = outlocal.loads;
        out.iterations = outlocal.detail.critfinal(3);
        out.err = outlocal.detail.ssq.residual;
    case 'nway'
        if opt.init~=0
            [Factors,it,err,~] = nwayparafac(tensor,f,...
                [opt.stopcriteria.relativechange opt.init 0 0 50 opt.stopcriteria.iterations],...
                [opt.constraints{1}.type opt.constraints{2}.type opt.constraints{3}.type]);
        elseif opt.init==0
            [Factors,it,err,~] = nwayparafac(tensor,f,...
                [opt.stopcriteria.relativechange opt.init 0 0 50 opt.stopcriteria.iterations],...
                [opt.constraints{1}.type opt.constraints{2}.type opt.constraints{3}.type],...
                initvalues,[0 0 0]);
        end
        out.model = Factors;
        out.iterations = it;
        out.err = err;
        clearvars Factors it err
        
end
end

%% 
function [hittype,printinfo] = scandiary(diary,tbox)            %hittype
switch tbox
    case 'nway'
        targets={'PRELIMINARY',2,[];...                                 %1
            'Sum-of-Squares   Iterations  Explained',2,1:5;...          %2
            'The algorithm converged',0,[];...                          %3
            'Warning: Matrix is close to singular ',0,[];...            %4
            'The misfit is approaching the machine uncertainty',0,[];...%5
            'WARNING, SOME FACTORS ARE HIGHLY CORRELATED.',0,[]};       %6
    case 'pls'
        targets={'Fitting new PARAFAC ...',2,[];...                     %1
            'plsupdate',0,[];};                                         %2
end

    
hittype=0;
printinfo='';
for n=1:size(targets,1)
    if ~strcmp(targets{n},'plsupdate')
        res=strfind(diary,targets{n});
        res=find(not(cellfun('isempty',res)));
        if ~isempty(res)
            try
                printinfo=diary{res+targets{n,2}};
                hittype=n;
            catch
                printinfo='';
                hittype=0;
            end
            return
        end
    elseif strcmp(targets{n},'plsupdate')
        idx=not(cellfun('isempty',cellfun(@str2num,diary,'UniformOutput',false)));
        if any(idx)
            printinfo=diary{idx};
            hittype=2;
        end
    end
end
end

%%
function [diariesOld] = printupdates(diaries,diariesOld,modelID,tbox)
% Function for reading futures Diary of parfeval futures
if numel(diariesOld)~=numel(diaries)
    sizeOld=numel(diariesOld); 
    sizeNew=numel(diaries); 
    dispIdx=(sizeOld+1):(sizeNew-1);
    [hittype,printinfo] = scandiary(diaries(dispIdx),tbox);
    switch hittype
        case 0 % Nothing to disply
        case 1 % New model
            disp(['mod-id ',sprintf('%03d',modelID),':  Model initialized. ',printinfo])
        case 2 % Display prelim. fit numbers
            switch tbox
                case 'nway'
                    fitinfo=str2num(printinfo);
                    if ~isempty(fitinfo)
                        disp(['mod-id: ',     sprintf('%03d',modelID),...
                            ':    it. ',    sprintf('%-*s',4,num2str(fitinfo(2))),...
                            '     err.: ', sprintf('%.3f',(round(fitinfo(1),3))),...
                            '     %exp.: ',  sprintf('%.2f',(round(fitinfo(3),2)))]);
                    end
                case 'pls'
                    fitinfo=str2num(printinfo);
                    if ~isempty(fitinfo)
                        disp(['mod-id: ',sprintf('%03d',modelID),...
                            ':    it.: ',sprintf('%-*s',4,num2str(fitinfo(1))),...
                            '     err.: ',sprintf('%.3f',(round(fitinfo(4),3))),...
                            '     chgfit.: ',sprintf('%.2u',fitinfo(3))]);
                    end
                    
            end
        case 4 % warning issued
        msgbox({'Warning during PARAFAC fitting: ',...
            'Warning: Matrix is close to singular or badly scaled.',...
            'Results may be inaccurate',' ',...
            'Acknowledge this message by clicking ''ok''';},...
            'WARNING','warn','replace');
        case 5 % warning issued
        msgbox({'Warning during PARAFAC fitting: ',...
            'The misfit is approaching the machine uncertainty.',...
            'If pure synthetic data is used this is OK. Otherwise',...
            'multiply the EEMs by a large number',' ',...
            'Acknowledge this message by clicking ''ok''';},...
            'WARNING: Machine Precision','warn','replace');
        case 6 % warning issued
        msgbox({'Warning during PARAFAC fitting: ',...
            'SOME FACTORS ARE HIGHLY CORRELATED.',...
            'Consider modeling fewer components.',...
            'Acknowledge this message by clicking ''ok''';},...
            'WARNING: Correlated factors','warn','replace');
    end
    diariesOld=diaries;
end
end

%%
function restoreoldpath(oldp,newp)
    if ~isequal(oldp,newp)
        path(oldp)
    end
end

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
        if ~[sz(1)<sz(2)]
            vout=v';
        else
            vout=v;
        end
    case 'column'
        if ~[sz(1)>sz(2)]
            vout=v';
        else
            vout=v;
        end
    otherwise
            error('Input ''rc'' not recognized. Options are: ''row'' and ''column''.')
end


end

%% Internal functions specifically for randinitanal
%%
function [cncl] = blockbar(names,state)
% Blockbar is similar to a multiwaitbar, but visualizes the state of
% individual jobs rather than a % completed bar.
% (c) Urban Wünsch, 2019
%% Default options
incompcol = [0.8         0.8         0.8];
compcol   = [0.204       0.569           1];

multi=numel(names);
init=true; % init state is true by default. When old figure is found, init = false.
cncl=false(numel(names),1);
%% Find existing blockbar and axes within it (not the buttons)
figHandles =  get(groot, 'Children');
for n=1:numel(figHandles)
    if strcmp(figHandles(n).Name,'Progress...')
        init=false; % Success, do not make another figure, use the old
        old=figHandles(n);
        break
    end
end
% Close blockbar if 'close' is requested by user
if strcmp(class(names),'char')&&sum(strcmp(names,'close'))==1
    try
       close(old)
       return
    catch
        warning('Coudn''t close the blockbar. Probably, there was none.')
        return
    end
end
% Check input
if numel(names)~=size(state,1)
    warning('Size of ''names'' and ''state'' are incompatible. Can'' draw the blockbar...')
    return
end
%% Draw the blockbar
stp=[0:1/size(state,2):1 1];
if init % This generates the figure and bars
    set(0, 'Units', 'pixel');
    screenSize = get(0,'ScreenSize');
    pointsPerPixel = 72/get(0,'ScreenPixelsPerInch');
    width = 360 * pointsPerPixel;
    height = multi * 75 * pointsPerPixel;
    fpos = [screenSize(3)/2-width/2 screenSize(4)/2-height/2 width height];
    figureh=figure('Position',fpos,...
        'MenuBar','none',...
        'Numbertitle','off',...
        'Name','Progress...',...
        'Resize','off');
    axeswidth = 172;
    axesheight = 172/14;
    axesbottom = fliplr(linspace(15,fpos(4)-axesheight*3,multi));
    addprop(figureh,'state');
    addprop(figureh,'sname');
    figureh.state=state;
    figureh.sname=names;
    for i=1:numel(names)
        axPos = [axesheight*2.5 axesbottom(i) axeswidth axesheight];
        bPos = [axPos(1)+axPos(3)+10 axPos(2) 50 axPos(4)*1.5];
        bPos(2) = bPos(2)-0.5*(bPos(4)-axPos(4));
        ax(i)=axes('units','pixel','pos',axPos);
        addprop(ax(i),'sname');
        set(ax(i),'units','pixel','Tag',names{i});
        uic=uicontrol( 'Style', 'togglebutton', ...
            'String', '', ...
            'FontWeight', 'Bold', ...
            'units','pixel',...
            'Position',bPos,...
            'Tag',num2str(i),...
            'String','Cancel','FontWeight','normal');
        addprop(uic,'snames');
        uic.snames=i;
        
        title(ax(i),names{i});
        set(ax(i),'YTickLabel','');
        set(ax(i),'XTickLabel','');
        for n=1:size(state,2)
            if state(i,n)
                col=compcol;
            else
                col=incompcol;
            end
            p=patch(ax(i),[stp(n) stp(n+1) stp(n+1) stp(n)],[0 0 1 1],col,...
                'EdgeColor','none');hold on
            addprop(p,'id');
            p.id=n;
        end
        set(ax(i),'YColor',[0 0 0 0.5],'XColor',[0 0 0 0.5],'Box','on','YTick','','XTick','')
    end
else % This updates  the figure and bars
    oldstate=old.state;
    if ~isequal(oldstate,state)
        ax=get(old,'Children');
        ax=ax(strcmp(get(ax,'Type'),'axes'));
        ni=multi:-1:1;
        for n=1:numel(names)
            newstate=state(n,:);
            if ~isequal(oldstate,newstate)
                h=get(ax(ni(n)),'Children');
                ii=size(state,2):-1:1;
                for i=1:size(state,2)
                    if state(n,i)&&~oldstate(n,i)
                        h(ii(i)).FaceColor=compcol;
                    end
                end
            end
        end
        old.state=state;
    end
    childs=get(old,'Children');
    tbutt=childs(strcmp(get(childs,'Type'),'uicontrol'));
    ni=multi:-1:1;    
    for n=1:numel(names)
        if logical(tbutt(ni(n)).Value)&any(~state(n,:))
            cncl(n)=true;
        else
            cncl(n)=false;
        end
    end
    
end
end

%%
function [best,all,details] = analyzeresults(data,factors,it,MaxIt,facCalls,CC,constraints,Model,Iter,Err,corecon,initstyle)

take=~cellfun(@isempty,Model);
notake=find(cellfun(@isempty,Model));
numcomplete=sum(take);
if numel(facCalls)~=numcomplete
    warning on
    warning('User interuption before all iterations were completed.')
    for ii=1:numel(notake)
        Iter{notake(ii)}=nan;
        Err{notake(ii)}=nan;
        corecon{notake(ii)}=nan;
    end
end
Model   =          reshape(Model,  it,numel(factors));
Iter    = cell2mat(reshape(Iter,   it,numel(factors)));
Err     = cell2mat(reshape(Err,    it,numel(factors)));
corecon = cell2mat(reshape(corecon,it,numel(factors)));

all=data;
best=data;
F=fieldnames(best);
best=rmfield(best,F(~cellfun(@isempty,strfind(F,'run'))));clearvars F
convergence=cell(1,numel(factors));
nonconvmessage='';
for i=1:numel(factors)
    try
        f=factors(i);
        % #1: All results
        for j=1:size(Model(:,i),1)
            all.(['Model',num2str(f),'_run',    int2str(j)])=Model{j,i};
            all.(['Model',num2str(f),'err_run', int2str(j)])=Err(j,i);
            all.(['Model',num2str(f),'it_run',  int2str(j)])=Iter(j,i);
            all.(['Model',num2str(f),'core_run',int2str(j)])=corecon(j,i);
        end
        
        % #2: Warn when critwarn % of models did not converge (currently 30%)
        critwarn=0.3;
        pnc=sum(Iter(:,i)>=MaxIt)./it;
        if pnc>=critwarn
            nonconvmessage=[nonconvmessage newline num2str(factors(i)),'-component models: ',num2str(pnc*100),' % of models did not converge. ']; %#ok<AGROW>
        end
        
        % #3: LSQ result
        niterr=[Err(:,i) Iter(:,i)]; % Matrix of model results
        niterr(niterr(:,2)>=MaxIt,:)=nan;    % 'NaN' all unconverged models
        
        % New: if unconstrained mode, check for 2FDG:
        tfd=false(it,1);
        % This is a feature that needs more work. Therefore, deactivated for now.
%         if ~any(strcmp(constraints,{'nonnegativity','unimodnonneg'}))
%             for iii=1:it
%                 x=Model{iii,i}{1};
%                 combs=nchoosek(1:i,2);
%                 for n=1:size(combs,1)
%                     r(n)=corr(x(:,combs(n,1)),x(:,combs(n,2)));
%                 end
%                 if any(abs(r)>0.99)
%                     tfd(iii)=true;
%                 else
%                     tfd(iii)=false;
%                 end
%             end
%         end
        niterr(tfd,:)=nan;    % 'NaN' degenerate models
        if sum(isnan(niterr(:,1)))==it
            error('randinitanal: No convergence.')
        end
        [~,min_i]=min(niterr(:,1),[],'omitnan');  % Find LSS
        %calculation for component size
        sizeF=nan(1,f);
        for ii=1:f
            modelledF=nmodel([{Model{min_i,i}{1}(:,ii)} {Model{min_i,i}{2}(:,ii)} {Model{min_i,i}{3}(:,ii)}]);
            sizeF(ii)=100 * (1 - (sum((data.X(:)-modelledF(:)).^2,'omitnan') / sum(data.X(:).^2,'omitnan')));

        end
        best.(['Model' int2str(f)])               = Model{min_i,i};
        best.(['Model' int2str(f) 'err' ])        = Err(min_i,i);
        best.(['Model' int2str(f) 'it' ])         = Iter(min_i,i);
        best.(['Model' int2str(f) 'core' ])       = corecon(min_i,i);
        best.(['Model' int2str(f) 'source' ])     = ['Model' int2str(f) 'it_' int2str(min_i)];
        best.(['Model' int2str(f) 'convgcrit'])   = CC;
        best.(['Model' int2str(f) 'constraints']) = constraints;
        best.(['Model' int2str(f) 'initialise'])  = initstyle;
        best.(['Model' int2str(f) 'percentexpl']) = 100 * (1 - Err(min_i,i) / sum(data.X(:).^2,'omitnan'));
        best.(['Model' int2str(f) 'compsize'])= sizeF;      
    catch
        warning(['No valid model found, no Model',num2str(f),' created.'])
        best.(['Model' int2str(f)])= {'No valid model found'};
        best.(['Model' int2str(f) 'err' ])= {};
        best.(['Model' int2str(f) 'it' ])=  {};
        best.(['Model' int2str(f) 'core' ])= {};
        best.(['Model' int2str(f) 'source' ])=  {};
        best.(['Model' int2str(f) 'convgcrit'])= CC;
        best.(['Model' int2str(f) 'constraints'])= constraints;
        best.(['Model' int2str(f) 'initialise'])= initstyle;
        best.(['Model' int2str(f) 'percentexpl'])=NaN;
        best.(['Model' int2str(f) 'compsize'])=NaN;
    end
    convergence{i}=[(1:it)' Err(:,i) Iter(:,i)];
end

if ~isempty(nonconvmessage)
    nonconvmessage=[nonconvmessage newline newline 'Increase of number of random starts recommended!' newline 'Please investigate if the issue persists.'];
    warning(nonconvmessage)
    msgbox(nonconvmessage,...
        'WARNING: Unconverged models','warn','replace');
end
details=array2table([facCalls cell2mat(convergence')],...
    'VariableNames',{'components','init','error','iterations'});
end

%%
function [Model,Iter,Err,corecon,ttc] = trackprogress(futures,numtry,facCalls,toolbox,data,Em,Ex,consoleoutput)
% Monitor parfeval progress supervision (c) Urban Wünsch, 2016-2019

% Allocation of diary variables (aka console output in parfeval)
diariesOld = cell(numel(futures)); % Old parfeval diaries
diaries    = cell(numel(futures)); % parfeval formatted diary

% Allocation of PARAFAC outputs
Model = cell(numtry,1);            % Allocate Model cell
Iter = cell(numtry,1);             % Allocate cell for #of interations
Err = cell(numtry,1);              % Allocate cell for SSE
corecon = cell(numtry,1);          % Allocate cell for core consistency

% Allocation of monitoring variables
fetchthese=true(1,numel(facCalls));% fetchNext: only these (not cancelled)
completed    = false(numtry,1);    % Complete / incomplete?
numCompleted = 0;                  % Total number of completed models
ttc=repmat(duration,numtry,1);     % time-to-convergence for each run

% Allocation of variables for prelim. plotfacs
[facs,iF]=unique(facCalls);
idx=nan(numel(facs),1);
idxOld=nan(numel(facs),1);

% 1st Blockbar
state=false(numel(facs),numtry/numel(facs));
mname=cellstr(strcat(repmat('Model',numel(facs),1),num2str((facs))))';
% Close any old blockbar
figHandles = get(0,'Children');
for n=1:numel(figHandles)
    if strcmp(figHandles(n).Name,'Progress...')
        close(figHandles(n));
    end
end
clearvars n figHandles
blockbar(mname,state);


while numCompleted < numtry
    ffetch=futures(fetchthese&~completed');
    [completedIdx,Modelres] = fetchNext(ffetch,0.1);
    % If fetchNext returned an output, let's extract it.
    if ~isempty(completedIdx)
        % Get actual index of future if ffetch-futures were subset
        [~,completedIdx]=intersect([futures.ID],ffetch(completedIdx).ID);
        numCompleted = numCompleted + 1;
        % Update list of completed futures.
        completed(completedIdx) = true;
        % Update state to account for completed futures
        state=reshape(completed,numtry/numel(facs),numel(facs))';
        Model{completedIdx,1} = rcvec(Modelres.model,'row');
        Iter{completedIdx,1} = Modelres.iterations;
        Err{completedIdx,1} = Modelres.err;
        try
            corecon{completedIdx,1} = corcond(data,Model{completedIdx,1},[],0);
        catch
            warning('Could not calculate core [Unknown error in N-way toolbox.]')
            corecon{completedIdx,1} = NaN;
        end
        
        ttc(completedIdx)=datetime(futures(completedIdx).FinishDateTime,'TimeZone','Europe/London')-...
            datetime(futures(completedIdx).StartDateTime,'TimeZone','Europe/London');
        if ~strcmp(consoleoutput,'none')||strcmp(consoleoutput,'minimum')
            disp(['mod-id ',sprintf('%03d',completedIdx),' done. | #it.: ',...
                sprintf('%-*s',4,num2str(Iter{completedIdx,1})),' | err.: ',num2str(Err{completedIdx,1}),...
                ' | core%: ',num2str(round(corecon{completedIdx,1})),...
                ' | time: ',char(ttc(completedIdx)),...
                ' | sec./it.: ',num2str(round(seconds(ttc(completedIdx))/Iter{completedIdx,1},2))]);
        end
    else % Analyze parfeval diary if future is still running.
        for i=1:numel(futures)
            if strcmp(futures(i).State,'running')
                if ~strcmp(consoleoutput,'none')||strcmp(consoleoutput,'minimum')
                    diaries{i}=strsplit(futures(i).Diary,'\n');
                    [diariesOld{i}]=printupdates(diaries{i},diariesOld{i},i,toolbox);
                end
            end
        end
    end
    % Check status of blockbar
    c=blockbar(mname,state);
    
    % Check if any models were cancelled and cancel the unfinished jobs among these models
    fetchthese=~repelem(c,numtry/numel(unique(facCalls)),1)';
    if any(~fetchthese)
        if ~isempty(futures(~fetchthese&~completed'))
            cancel(futures(~fetchthese&~completed'));
            Err(~fetchthese&~completed',1)=num2cell(nan(sum(~fetchthese&~completed'),1));
            numCompleted=numCompleted+numel(futures(~fetchthese&~completed'));
            completed(~fetchthese) = true;
            warning('User canceled some models prematurely. No output fetched for those.')
        end
    end
    
    
    for n=1:numel(facs)
        if ~all(isnan(cell2mat(Err(facCalls==facs(n)))))&&all(completed(facCalls==facs(n)))
            [~,idx(n)]=min(cell2mat(Err(facCalls==facs(n))));
            idx(n)=idx(n)+iF(n)-1;
        else
            idx(n)=nan;
        end
    end
    if ~isequal(isnan(idxOld),isnan(idx))&&~all(isnan(idx))
        try
            plotfacs(Model(idx(~isnan(idx))),facCalls(idx(~isnan(idx)))',[],Em,Ex);
        catch
            disp('Couln''t plot that')
        end
        idxOld=idx;
    end
    
end
% If complete, cancel the futures and delete the waitbar.
cancel(futures);
blockbar('close');
end

%%
function plotfacs(ModelFin,factors,ford,Em,Ex)
figHandles = get(0,'Children');
if isempty(figHandles)
    hf=dreemfig;
    w=.8;h=.4;l=(1-w)/2;b=(1-h)/2;
    set(hf, 'units','normalized','outerposition',[l b w h],...
        'Name','Scores / Spectral loadings plot');
end
figHandles = get(0,'Children');

if any(contains({figHandles.Name},'Spectral'))
    if sum(contains({figHandles.Name},'Spectral'))>1
        close(figHandles(contains({figHandles.Name},'Spectral')))
        hf=dreemfig;
        w=.8;h=.4;l=(1-w)/2;b=(1-h)/2;
        set(hf, 'units','normalized','outerposition',[l b w h],...
            'Name','Scores / Spectral loadings plot');
    else
        ax = (findobj(figHandles(contains({figHandles.Name},'Spectral')), 'type', 'axes'));
        for n=1:numel(ax)
            delete(ax(n))
        end
        hf=figure(figHandles(contains({figHandles.Name},'Spectral')));
    end
else
    hf=dreemfig;
    w=.8;h=.4;l=(1-w)/2;b=(1-h)/2;
    set(hf, 'units','normalized','outerposition',[l b w h],...
        'Name','Scores / Spectral loadings plot');
end
mfac=max(factors);
nfac=length(factors);
plotcount=1;
for j=1:nfac
    [A,B,C]=fac2let(ModelFin{j});
    
    if ~isempty(ford)
        if iscell(ford)
            fordbyn=ford{j};
        else
            fordbyn=ford;
        end
        
        if ~isempty(fordbyn)
            A=A(:,fordbyn);
            B=B(:,fordbyn);
            C=C(:,fordbyn);
        end
    end
    subplot(nfac,mfac+1,plotcount);
    plot(A,'LineWidth',1.2,'Marker','.','MarkerSize',7)
    if all(A>=0)
        ylim([0 max(A(:))])
    end
    title('Scores')
    xlabel('Sample')
    ylabel(['Scores: ' num2str(factors(j))])
    plotcount=plotcount+1;
    col=lines(mfac);
    for i=1:mfac
        try
            B(:,i);
            subplot(nfac,mfac+1,plotcount);
            plot(Em,B(:,i),'-','color',col(i,:),'LineWidth',1.2);
            if j==nfac
                xlabel('Wave. (nm)');
            end
            if i==1
                ylabel(['Loads: ' num2str(factors(j))]);
            end
            hold on
            plot(Ex,C(:,i),'-.','color',col(i,:),'LineWidth',1.2);
            axis tight
            grid on
            plotcount=plotcount+1;
        catch
            plotcount=plotcount+1;
        end
    end
end
dreemfig(hf);

if any(contains({figHandles.Name},'Progress...'))
    figure(figHandles(contains({figHandles.Name},'Progress...')))
end

end