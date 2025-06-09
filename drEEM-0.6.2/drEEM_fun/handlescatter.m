function [dataout] = handlescatter(data,varargin)
%
% <strong>Syntax</strong>
%
%   dataout = <strong>handlescatter</strong>(data,options)
% <a href="matlab:opt=handlescatter('options')">opt=handlescatter('options')</a>
%
%
% <a href="matlab: doc handlescatter">help for handlescatter</a> <- click on the link

% Excise EEM scatter and (optionally) interpolate between missing values.
% Primary and secondary Rayleigh and Raman are removed and interpolated if 
% requested, or left as NaNs. Zeros may be placed at a specified
% distance below the line Em=Ex. Optional plots can be shown that compare 
% the EEMs before and after excising/smoothing.
%
% USEAGE:
% [dataout] = handlescatter(data,options) OR
% [dataout] = handlescatter(Xin,Ray1,Ram1,Ray2,Ram2,NaNfilter,d2zero,freq,plotview)
%
%
% OPTIONS:
% Options can be conveniantly supplied as a single strucure. Default values are
% obtained by calling 'opt = handlescatter('options')'. The values in the
% fields of opt can then be modified as desired.
%
% IMPORTANT: handlescatter accepts inputs just like they used to be provided with 
% smootheem. If you are used to providing inputs the "old" way, handlescatter
% can also digest these inputs. Some newer options will however be left at
% their default values!
%
% Fields:
% cutout=[1 1 1 1];         % Cut scatter?  [Ray1 Ram1 Ray2 Ram2]
% interpolate=[0 1 1 1];    % Interpolate scatter? [Ray1 Ram1 Ray2 Ram2]
% ray1=[10 10];             % Rayleigh 1st: [below above]
% ram1=[15 15];             % Raman 1st:    [below above]
% ray2=[5  5];              % Rayleigh 2nd: [below above]
% ram2=[5  5];              % Raman 2nd:    [below above]
% d2zero=60;                % Distance in nm below Ray1. When Em<Ray1-d2zero, values are forced zero
% iopt='normal';            % Chose if overlapping NaN's should be interpolated {'normal'|'conservative'}
% imethod='inpaint';        % Interpolation method {'inpaint'|'fillmissing'}
% negativetozero='on';      % If 'on', all negative values will be set to zero
% iopt='normal';            % Should overlapping NaN's areas  be left uninterpolated (conserved)? {'normal'|'conservative'}
% plot='on';                % Plot raw, cut, and final EEMs, sample-by-sample
%                             NOTE: Plots will always be shown for all samples,
%                             but simply closing the window will terminate
%                             plotting and return the smoothed data.
% samples='all'             % 'all': all samples will be cut as specified.
%                              If numeric vector is supplied, only part of the dataset will be treated.
% plottype='mesh';          % If plot is 'on', which type of plot should be shown {'mesh'|'surface'|'contourf'}
%
%
% OUTPUT:
%
%   dataout: A data structure with the smoothed EEM in dataout.X
%
%
% handlescatter: Copyright (C) 2020 Urban J. Wuensch
% Chalmers University of Technology
% Sven Hultins Gata 6
% 41296 Gothenburg
% Sweden
% wuensch@chalmers.se
% $ Version 0.1.0 $ September 2019 $ First Release
% 
% inpaintn: Copyright (c) 2017, Damien Garcia
% All rights reserved.
%
%% Check input arguments and react to different scenarios

% Scenario: Return options
if strcmp(data,'options')
    options=defaultoptions;
    if strcmp(data,'options')
        dataout=options;
        return
    end
end

% Scenario: No options provided
if isstruct(data)&&nargin==1
    options=defaultoptions;
    disp(' ')
    warning(sprintf(['No options were provided, the defaults were assumed.\n'...
        '     Please inspect the result and see if adjustments are necessary.\n'...
        '     Options can be obtained by calling ''handlescatter(''options'')'''])) %#ok<SPWRN>
    disp(' ')
end

% Scenario: data & options provided
if isstruct(data)&&nargin==2
    options=varargin{1};
end

% Scenario: Old input provided.
if nargin>3
    disp('     Old smootheem input detected. Attempting to parse these inputs...')
    disp('     The remaining options were left at default values:')
    disp('     If you wish to modify these, call opt=handlescatter(''options'')')
    disp('     and modify the values within the structure ''opt''')
    disp('     Type ''help handlescatter'' for more information.')
    disp('     The options currently are: ')
    disp(' ')
    pause(0.5)
    options=defaultoptions;
    optin=varargin;
    if numel(optin)~=8
        error(sprintf(['When classic smootheem input is provided, all options MUST be explicitly set (provided).' ...
            '\n\n smootheem(Xin,Ray1,Ram1,Ray2,Ram2,NaNfilter,d2zero,freq,plotview)' ...
            '\n\n see <a href="matlab: doc smootheem">documentation for smootheem</a>']))
    end
    fieldnames={'ray1' 'ram1' 'ray2' 'ram2' 'interpolate' 'd2zero' 'freq' 'plot'};
    for n=1:numel(fieldnames)
        if ~isempty(optin{n})
            if numel(optin{n})==2
                options.(fieldnames{n})=optin{n}([2 1]);
            else
                options.(fieldnames{n})=optin{n};
            end
        end
    end
    if strcmp(options.plot,'pause')
        options.plot='on';
    else
        options.plot='off';
    end
    disp(options)
    disp(' ')
end

%% Check user input for settings that will cause issues

if options.ray1(1)>options.d2zero&&options.interpolate(1)&&options.cutout(1)
    disp(' ')
    disp(' ')
    disp(sprintf(['    You are planning to cut and interpolate Rayleigh 1st order scatter\n'...
        '    and then to force part of the interpolation to zero.'])); %#ok<DSPS>
    disp(sprintf(['    It is advisable to leave some room between interpolation and forced zeros.\n'...
        '    options.d2zero was reset to options.ray1(1)+5 '])); %#ok<DSPS>
    disp(' ')
    pause(1)
    options.d2zero=options.ray1(1)+5;
end

if ~(ischar(options.samples)||isnumeric(options.samples))
    error('options.samples must either be character or numeric.')
end
if ischar(options.samples)
    if ~strcmp(options.samples,'all')
        error('options.samples must  be ''all'' if character is provided.')
    end
end

if ischar(options.samples)
    if strcmp(options.samples,'all')
        allsamples=true;
    end
end

if isnumeric(options.samples)
    allsamples=false;
end
%% Preparation of matricies that will determine treatment
% Missing-number matrix (will be used to NaN the EEMs)
nanmat=false(data.nEm,data.nEx);
% Interpolation matrix (will be used to identify interpolation scenarios)
imat=ones(data.nEm,data.nEx); % imat: 1(don't handle) 2(interpolate) 3(do no interpolate) (2&3 are set by labeler)
types={'ray1';'ram1';'ray2';'ram2'};
for n=1:numel(types)
    [nanmat,imat]=labeler(nanmat,imat,data.Ex,data.Em,types{n},options.(types{n})(1),options.(types{n})(2),options.interpolate(n),options.cutout(n));
end
clearvars types

%% Data treatment
if allsamples
    X=data.X;  %X  -> raw data
    
else
    X=data.X(options.samples,:,:);  %X  -> raw data
end
Xc=X;      %Xc -> X (cut)
% (1) NaN the scatter diagonals
Xc(:,nanmat)=nan;



% (2) Zero negatives (they may otherwise impact interpolation)
if strcmp(options.negativetozero,'on')
    Xc(Xc<0)=0;
end
Xi=Xc;     %Xi -> X (interpolated)


% (3) Interpolate ALL scatter types, some may be NaN'ed later (depends on settings)
if any(options.interpolate)
    Xi=interpolatetheeem(Xi,options.imethod);
end
Xi(:,imat==3&nanmat)=nan; % NaN the types that should not have been interpolated.


% (4) Zero negatives again to account for bad interpolation
if strcmp(options.negativetozero,'on')
    Xi(Xi<0)=0;
end

% (5) Zero data below 1st order Rayleigh
if ~isempty(options.d2zero)
    d2zeromat=false(data.nEm,data.nEx);
    [d2zeromat,~]=labeler(d2zeromat,[],data.Ex,data.Em,'d2zero',options.d2zero);
    Xi(:,d2zeromat)=0;
end

% (6) NaN areas of intersecting scatter (if desired)
if strcmp(options.iopt,'conservative')
    ex=200:800;
    em=200:1000;
    imattemplate=ones(numel(em),numel(ex)); % imat: 1(don't handle) 2(interpolate) 3(do no interpolate) (2&3 are set by labeler)
    types={'ray1';'ram1';'ray2';'ram2'};
    for n=1:numel(types)
        if options.cutout(n)
            [~,imat2{n}]=labeler([],imattemplate,ex,em,types{n},options.(types{n})(1),options.(types{n})(2),1);
            [compmat{n},]=labeler(false(data.nEm,data.nEx),[],data.Ex,data.Em,types{n},options.(types{n})(1),options.(types{n})(2),1);
        end
    end

    firstorderexcl=imat2{1}+imat2{2};
    secondorderexcl=imat2{3}+imat2{4};

    clearvars imattemplate imat2 types n

    exlim=ex(find(any(firstorderexcl==4,1),1,'last'));
    emlim=em(find(any(firstorderexcl==4,2),1,'last'));

    nanmat2=false(data.nEm,data.nEx);
    if exlim>=min(data.Ex)&&emlim>=min(data.Em)
        nanmat2(:,1:knnsearch(data.Ex,exlim))=true;
        nanmat2(1:knnsearch(data.Ex,exlim),:)=true;
    end
    nanmat3=nanmat2&compmat{1}|nanmat2&compmat{2};
    exlim=ex(find(any(secondorderexcl==4,1),1,'last'));
    emlim=em(find(any(secondorderexcl==4,2),1,'last'));

    nanmat2=false(data.nEm,data.nEx);
    if exlim>=min(data.Ex)&&emlim>=min(data.Em)
        nanmat2(:,1:knnsearch(data.Ex,exlim))=true;
        nanmat2(1:knnsearch(data.Ex,exlim),:)=true;
    end
    nanmat4=nanmat2==true&compmat{1}==true|nanmat2==true&compmat{2}==true;
    
    nanmat5=nanmat3|nanmat4;
    clearvars nanmat4 nanmat3 nanmat2 exlim emlim firstorderexcl secondorderexcl
    Xi(:,nanmat5)=NaN;
end

%% Output variable definition
dataout=data;
if allsamples
    dataout.X=Xi;
else
    dataout.X(options.samples,:,:)=Xi;
end
if ~isfield(dataout,'Smooth')
    dataout.Smooth=options;
else
    try
        dataout.Smooth=[dataout.Smooth options];
    catch
        error('Could not add settngs information to the dataset. Delete field "data.Smooth". If you''re tyring to smooth one dataset with multiple settings, only use ''handlescatter''.')
    end
end
%% Plotting

if strcmp(options.plot,'on')
    fh=dreemfig;
    set(fh, 'units','normalized','pos',[0.05    0.1611    0.9    0.3700]);
    ax=gobjects(1,3);
    for n=1:3
        ax(n)=subplot(1,3,n);
    end
    huic = uicontrol('Style', 'pushbutton','String','Next',...
    'Units','normalized','Position', [0.9323 0.0240 0.0604 0.0500],...
    'Callback',{@pltnext});
    huic2 = uicontrol('Style', 'pushbutton','String','Close figure',...
    'Units','normalized','Position', [0.9323 0.0816 0.0604 0.0500],...
    'Callback',{@endfunc,fh});
    az=[83.2,83.2,83.2];
    el=[55.26,55.26,55.26];
    for n=1:dataout.nSample
        switch options.plottype
            case 'mesh'
                mesh(ax(1),dataout.Ex,dataout.Em,squeeze(X(n,:,:)))
                mesh(ax(2),dataout.Ex,dataout.Em,squeeze(Xc(n,:,:)))
                mesh(ax(3),dataout.Ex,dataout.Em,squeeze(Xi(n,:,:)))
                for k=1:3
                    view(ax(k),az(k),el(k))
                    colormap(ax(k),cmap)
                    ylabel(ax(k),'Emission (nm)')
                    xlabel(ax(k),'Excitation (nm)')
                    zlabel(ax(k),'Intensity')
                end
            case 'surface'
                for k=1:3,cla(ax(k),'reset');end
                surface(ax(1),dataout.Ex,dataout.Em,squeeze(X(n,:,:)),'EdgeAlpha',0.5)
                surface(ax(2),dataout.Ex,dataout.Em,squeeze(Xc(n,:,:)),'EdgeAlpha',0.5)
                surface(ax(3),dataout.Ex,dataout.Em,squeeze(Xi(n,:,:)),'EdgeAlpha',0.5)
                for k=1:3
                    view(ax(k),az(k),el(k))
                    colormap(ax(k),cmap)
                    ylabel(ax(k),'Emission (nm)')
                    xlabel(ax(k),'Excitation (nm)')
                    zlabel(ax(k),'Intensity')
                    hold(ax(k),'off')
                end
            case 'contourf'
                contourf(ax(1),dataout.Ex,dataout.Em,squeeze(X(n,:,:)),50,'LineStyle','none')
                contourf(ax(2),dataout.Ex,dataout.Em,squeeze(Xc(n,:,:)),50,'LineStyle','none')
                contourf(ax(3),dataout.Ex,dataout.Em,squeeze(Xi(n,:,:)),50,'LineStyle','none')
                for k=1:3
                    colormap(ax(k),cmap)
                    ylabel(ax(k),'Emission (nm)')
                    xlabel(ax(k),'Excitation (nm)')
                end
        end
        title(ax(1),'Raw data'),title(ax(2),'Cut'),title(ax(3),'Final')
        uicontrol(huic)
        uiwait(fh)
        if ~ishandle(fh); return; end % Ends function when plot is closed by user
        
        for k=1:3,[az(k),el(k)]=view(ax(k));end
    end
end
end
%% Internal functions used above
%%
function options=defaultoptions
options.cutout=[1 1 1 1];                   % Cut scatter?  [Ray1 Ram1 Ray2 Ram2]
options.interpolate=[0 1 1 1];              % Interpolate scatter? [Ray1 Ram1 Ray2 Ram2]
options.negativetozero='on';                % If 'on', all negative values will be set to zero
options.ray1=[10 10];                       % Rayleigh 1st: [below above]
options.ram1=[15 15];                       % Raman 1st:    [below above]
options.ray2=[5  5];                        % Rayleigh 2nd: [below above]
options.ram2=[5  5];                        % Raman 2nd:    [below above]
options.d2zero=60;                          % Distance in nm below Ray1. When Em<Ray1-d2zero, values are forced zero
options.imethod='inpaint';                  % Interpolation method {'inpaint'|'fillmissing'}
options.iopt='normal';                      % Chose if overlapping NaN's areas should be left uninterpolated (conserved) {'normal'|'conservative'}
options.plot='on';                          % Plot raw, cut, and final EEMs, sample-by-sample
options.plottype='mesh';                    % If plot is 'on', which type of plot should be shown {'mesh'|'surface'|'contourf'}
options.samples='all';                      % If 'all', all EEMs will be treated, if numeric vector, only the specified samples will be treated.
options.description='Options for handlescatter.m';
end

%%
function [nanout,iout]=labeler(nanin,iin,x,y,type,below,above,iswitch,cutswitch)
nanout=nanin;
iout=iin;
syncx=nan(numel(x),numel(y));
% switch-case for different types (d2zero immediately returns)
switch type
    case 'ray1'
        for n=1:numel(x)
            syncx(n,:) = y-x(n);
        end
    case 'ray2'
        for n=1:numel(x)
            syncx(n,:) = y-2*x(n);
        end
    case 'ram1'
        for n=1:numel(x)
            syncx(n,:)= y -(1e7*((1e7)/(x(n))-3382)^-1);
        end
    case 'ram2'
        for n=1:numel(x)
            syncx(n,:) = y - ((1e7*((1e7)/(x(n))-3382)^-1)*2);
        end
    case 'd2zero'
        for n=1:numel(x)
            syncx=y-x(n);
            nanout(syncx<=-below,n)=true;
        end
        return
    otherwise
        error('Input to ''type'' must be one of these: ray1 ray2 ram1 ram2 d2zero')
end

for n=1:numel(x)
    if cutswitch
        nanout(syncx(n,:)>=-below&syncx(n,:)<=above,n)=true;
    else
        nanout(syncx(n,:)>=-below&syncx(n,:)<=above,n)=false;
    end
    if iswitch
         iout(syncx(n,:)>=-below&syncx(n,:)<=above,n)=2;
    else
        iout(syncx(n,:)>=-below&syncx(n,:)<=above,n)=3;
    end
end

end

%%
function Xout=interpolatetheeem(Xin,method)
Xout=zeros(size(Xin,1),size(Xin,2),size(Xin,3));


switch method
    case 'inpaint'
        disp(['    inpaint-interpolation takes time. Please wait... (approx. ',num2str(round(1.6/40000*numel(Xin)*2/60,1)),'min)'])        
        parfor n=1:size(Xin,1)
            x=squeeze(Xin(n,:,:));
            x=inpaintn(x,500);
            Xout(n,:,:)=x;
        end
    case 'fillmissing'
        for n=1:size(Xin,1)
            x=squeeze(Xin(n,:,:));
            x=fillmissing(x,'pchip',1,'EndValues','nearest');
            Xout(n,:,:)=x;
        end
end
        
end

%%
function pltnext(sosurce,event) %#ok<INUSD>
uiresume
end

%%
function endfunc(~,~,hfig) 
close(hfig)
end

%%
function cols=cmap
cols=parula;
end


%% INPAINTN functions
function y = inpaintn(x,n,y0,m)
% INPAINTN Inpaint over missing data in N-D array
%   Y = INPAINTN(X) replaces the missing data in X by extra/interpolating
%   the non-missing elements. The non finite values (NaN or Inf) in X are
%   considered as missing data. X can be any N-D array.
%
%   INPAINTN (no input/output argument) runs the following 3-D example.
%
%   Important note:
%   --------------
%   INPAINTN uses an iterative process baased on DCT and IDCT.
%   Y = INPAINTN(X,N) uses N iterations. By default, N = 100. If you
%   estimate that INPAINTN did not totally converge, increase N:
%   Y = INPAINTN(X,1000)
%
%   Y = INPAINTN(X,N,Y0) uses Y0 as initial guess. This could be useful if
%   you want to run the process a second time or if you have a GOOD guess
%   of the final result. By default, INPAINTN makes a nearest neighbor
%   interpolation (by using BWDIST) to obtain a rough guess.
%
%   References (please refer to the two following references)
%   ---------- 
%   1) Garcia D, Robust smoothing of gridded data in one and higher
%   dimensions with missing values. Computational Statistics & Data
%   Analysis, 2010;54:1167-1178. 
%   <a
%   href="matlab:web('http://www.biomecardio.com/publis/csda10.pdf')">download PDF</a>
%   2) Wang G, Garcia D et al. A three-dimensional gap filling method for
%   large geophysical datasets: Application to global satellite soil
%   moisture observations. Environ Modell Softw, 2012;30:139-142.
%   <a
%   href="matlab:web('http://www.biomecardio.com/publis/envirmodellsoftw12.pdf')">download PDF</a>
%
%   Examples
%   --------
%
%     %% ---- RGB image ---- %%
%     onion = imread('onion.png');
%     I = randperm(numel(onion));
%     onionNaN = double(onion); onionNaN(I(1:round(numel(I)*0.5))) = NaN;
%     subplot(211), imshow(uint8(onionNaN)), title('Corrupted image - 50%')
%     for k=1:3, onion(:,:,k) = inpaintn(onionNaN(:,:,k)); end
%     subplot(212), imshow(uint8(onion)), title('Inpainted image')
%
%     %% ---- 2-D data ---- %%
%     n = 256;
%     y0 = peaks(n);
%     y = y0;
%     I = randperm(n^2);
%     y(I(1:n^2*0.5)) = NaN; % lose 1/2 of data
%     y(40:90,140:190) = NaN; % create a hole
%     z = inpaintn(y,200); % inpaint data
%     subplot(2,2,1:2), imagesc(y), axis equal off
%     title('Corrupted data')
%     subplot(223), imagesc(z), axis equal off
%     title('Recovered data ...')
%     subplot(224), imagesc(y0), axis equal off
%     title('... compared with original data')
%
%     %% ---- 3-D data ---- %%
%     load wind
%     xmin = min(x(:)); xmax = max(x(:));
%     zmin = min(z(:)); ymax = max(y(:));
%     %-- wind velocity
%     vel0 = interp3(sqrt(u.^2+v.^2+w.^2),1,'cubic');
%     x = interp3(x,1); y = interp3(y,1); z = interp3(z,1);
%     %-- remove randomly 90% of the data
%     I = randperm(numel(vel0));
%     velNaN = vel0;
%     velNaN(I(1:round(numel(I)*.9))) = NaN;
%     %-- inpaint using INPAINTN
%     vel = inpaintn(velNaN);
%     %-- display the results
%     subplot(221), imagesc(velNaN(:,:,15)), axis equal off
%     title('Corrupted plane, z = 15')
%     subplot(222), imagesc(vel(:,:,15)), axis equal off
%     title('Reconstructed plane, z = 15')    
%     subplot(223)
%     hsurfaces = slice(x,y,z,vel0,[xmin,100,xmax],ymax,zmin);
%     set(hsurfaces,'FaceColor','interp','EdgeColor','none')
%     hcont = contourslice(x,y,z,vel0,[xmin,100,xmax],ymax,zmin);
%     set(hcont,'EdgeColor',[.7,.7,.7],'LineWidth',.5)
%     view(3), daspect([2,2,1]), axis tight
%     title('Original data compared with...')
%     subplot(224)
%     hsurfaces = slice(x,y,z,vel,[xmin,100,xmax],ymax,zmin);
%     set(hsurfaces,'FaceColor','interp','EdgeColor','none')
%     hcont = contourslice(x,y,z,vel,[xmin,100,xmax],ymax,zmin);
%     set(hcont,'EdgeColor',[.7,.7,.7],'LineWidth',.5)
%     view(3), daspect([2,2,1]), axis tight
%     title('... reconstructed data')
%
%     %% --- 4-D data --- %%
%     [x1,x2,x3,x4] = ndgrid(-2:0.2:2);
%     z0 = x2.*exp(-x1.^2-x2.^2-x3.^2-x4.^2);
%     I = randperm(numel(z0));
%     % remove 50% of the data
%     zNaN = z0; zNaN(I(1:round(numel(I)*.5))) = NaN;
%     % reconstruct the data using INPAINTN
%     z = inpaintn(zNaN);
%     % display the results (for x4 = 0)
%     subplot(211)
%     zNaN(isnan(zNaN)) = 0.5;
%     slice(x2(:,:,:,1),x1(:,:,:,1),x3(:,:,:,1),zNaN(:,:,:,11),...
%        [-1.2 0.8 2],2,[-2 0.2])
%     title('Corrupt data, x4 = 0')
%     subplot(212)
%     slice(x2(:,:,:,1),x1(:,:,:,1),x3(:,:,:,1),z(:,:,:,11),...
%        [-1.2 0.8 2],2,[-2 0.2])
%     title('Reconstructed data')
%
%   See also SMOOTHN, GRIDDATAN
%
%   -- Damien Garcia -- 2010/06, last update 2017/08
%   website: <a
%   href="matlab:web('http://www.biomecardio.com/en')">www.BiomeCardio.com</a>
if nargin==0&&nargout==0, RunTheExample, return, end
class0 = class(x);
x = double(x);
if nargin==1 || isempty(n), n = 100; end
sizx = size(x);
d = ndims(x);
Lambda = zeros(sizx);
for i = 1:d
    siz0 = ones(1,d);
    siz0(i) = sizx(i);
    Lambda = bsxfun(@plus,Lambda,...
        cos(pi*(reshape(1:sizx(i),siz0)-1)/sizx(i)));
end
Lambda = 2*(d-Lambda);
% Initial condition
W = isfinite(x);
if nargin==3 && ~isempty(y0)
    y = y0;
    s0 = 3; % note: s = 10^s0
else
    if any(~W(:))
        [y,s0] = InitialGuess(x,isfinite(x));
    else
        y = x;
        return
    end
end
x(~W) = 0;
if isempty(n) || n<=0, n = 100; end
% Smoothness parameters: from high to negligible values
s = logspace(s0,-6,n); 
RF = 2; % relaxation factor
if nargin<4 || isempty(m), m = 2; end
Lambda = Lambda.^m;
% h = waitbar(0,'Inpainting...');
for i = 1:n
        Gamma = 1./(1+s(i)*Lambda);
        y = RF*idctn(Gamma.*dctn(W.*(x-y)+y)) + (1-RF)*y;
%         waitbar(i/n,h)
end
% close(h)
y(W) = x(W);
y = cast(y,class0);
end
%% Initial Guess
function [z,s0] = InitialGuess(y,I)
if license('test','image_toolbox')
    %-- nearest neighbor interpolation
    [~,L] = bwdist(I);
    z = y;
    z(~I) = y(L(~I));
    s0 = 3; % note: s = 10^s0
else
    warning('MATLAB:inpaintn:InitialGuess',...
        ['BWDIST (Image Processing Toolbox) does not exist. ',...
        'The initial guess may not be optimal; additional',...
        ' iterations can thus be required to ensure complete',...
        ' convergence. Increase N value if necessary.'])
    z = y;
    z(~I) = mean(y(I));
    s0 = 6; % note: s = 10^s0
end
end
%% DCTN
function y = dctn(y)
%DCTN N-D discrete cosine transform.
%   Y = DCTN(X) returns the discrete cosine transform of X. The array Y is
%   the same size as X and contains the discrete cosine transform
%   coefficients. This transform can be inverted using IDCTN.
%
%   Reference
%   ---------
%   Narasimha M. et al, On the computation of the discrete cosine
%   transform, IEEE Trans Comm, 26, 6, 1978, pp 934-936.
%
%   Example
%   -------
%       RGB = imread('autumn.tif');
%       I = rgb2gray(RGB);
%       J = dctn(I);
%       imshow(log(abs(J)),[]), colormap(jet), colorbar
%
%   The commands below set values less than magnitude 10 in the DCT matrix
%   to zero, then reconstruct the image using the inverse DCT.
%
%       J(abs(J)<10) = 0;
%       K = idctn(J);
%       figure, imshow(I)
%       figure, imshow(K,[0 255])
%
%   -- Damien Garcia -- 2008/06, revised 2011/11
%   -- www.BiomeCardio.com --
y = double(y);
sizy = size(y);
y = squeeze(y);
dimy = ndims(y);
% Some modifications are required if Y is a vector
if isvector(y)
    dimy = 1;
    if size(y,1)==1, y = y.'; end
end
% Weighting vectors
w = cell(1,dimy);
for dim = 1:dimy
    n = (dimy==1)*numel(y) + (dimy>1)*sizy(dim);
    w{dim} = exp(1i*(0:n-1)'*pi/2/n);
end
% --- DCT algorithm ---
if ~isreal(y)
    y = complex(dctn(real(y)),dctn(imag(y)));
else
    for dim = 1:dimy
        siz = size(y);
        n = siz(1);
        y = y([1:2:n 2*floor(n/2):-2:2],:);
        y = reshape(y,n,[]);
        y = y*sqrt(2*n);
        y = ifft(y,[],1);
        y = bsxfun(@times,y,w{dim});
        y = real(y);
        y(1,:) = y(1,:)/sqrt(2);
        y = reshape(y,siz);
        y = shiftdim(y,1);
    end
end
        
y = reshape(y,sizy);
end
%% IDCTN
function y = idctn(y)
%IDCTN N-D inverse discrete cosine transform.
%   X = IDCTN(Y) inverts the N-D DCT transform, returning the original
%   array if Y was obtained using Y = DCTN(X).
%
%   Reference
%   ---------
%   Narasimha M. et al, On the computation of the discrete cosine
%   transform, IEEE Trans Comm, 26, 6, 1978, pp 934-936.
%
%   Example
%   -------
%       RGB = imread('autumn.tif');
%       I = rgb2gray(RGB);
%       J = dctn(I);
%       imshow(log(abs(J)),[]), colormap(jet), colorbar
%
%   The commands below set values less than magnitude 10 in the DCT matrix
%   to zero, then reconstruct the image using the inverse DCT.
%
%       J(abs(J)<10) = 0;
%       K = idctn(J);
%       figure, imshow(I)
%       figure, imshow(K,[0 255])
%
%   See also DCTN, IDSTN, IDCT, IDCT2, IDCT3.
%
%   -- Damien Garcia -- 2009/04, revised 2011/11
%   -- www.BiomeCardio.com --
y = double(y);
sizy = size(y);
y = squeeze(y);
dimy = ndims(y);
% Some modifications are required if Y is a vector
if isvector(y)
    dimy = 1;
    if size(y,1)==1
        y = y.';
    end
end
% Weighing vectors
w = cell(1,dimy);
for dim = 1:dimy
    n = (dimy==1)*numel(y) + (dimy>1)*sizy(dim);
    w{dim} = exp(1i*(0:n-1)'*pi/2/n);
end
% --- IDCT algorithm ---
if ~isreal(y)
    y = complex(idctn(real(y)),idctn(imag(y)));
else
    for dim = 1:dimy
        siz = size(y);
        n = siz(1);
        y = reshape(y,n,[]);
        y = y.*w{dim};
        y(1,:) = y(1,:)/sqrt(2);
        y = ifft(y,[],1);
        y = real(y*sqrt(2*n));
        I = (1:n)*0.5+0.5;
        I(2:2:end) = n-I(1:2:end-1)+1;
        y = y(I,:);
        y = reshape(y,siz);
        y = shiftdim(y,1);            
    end
end
        
y = reshape(y,sizy);
end
