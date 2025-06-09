function [] = eemreview(DS,varargin)
%
% <strong>Syntax</strong>
%
%   <strong>eemreview</strong>(data,Name,Value)
%
% <strong>Name,Value arguments</strong>
% 
%    eem:       character
%    title:     character
%    samples:   numeric
%    nContours: numeric
%    hold:      true|false
%    LineStyle: see linespec
%    peaks:     true|false
%
% <a href="matlab: doc eemreview">help for eemreview</a> <- click on the link


% View EEMs for selection of samples with the option to plot Em and Ex scans for selected samples
%
% USEAGE:
%           [] = eemreview(DS,Name,value)
%
% INPUTS
%            DS:                drEEM-format dataset
%            (optional):        Name,value: Parameter name followed by option, e.g. 'PeakSwitch','on'
%            eem:               (char): Name of variable in DS that should be plotted.
%                                      Default: 'X'
%                                      Special: ModelX     -> Modeled data (X-component model)
%                                               ModelXresi -> Residuals of X-component model
%            title:             (cell) character cells providing title above EEM.
%                                      Default: DS.filelist
%            samples:           (numeric) Selection of samples for plotting
%                                      Default: all samples in the dataset
%            nContours:         (numeric) Number of contours
%                                      Default: 15. If > 50, LineStyle will be 'none' automatically
%            hold:              (logical)  if true, 2D spectra will be plotted superimposed.
%                                      Default: 'off'
%            LineStyle:         (char) LineStyle of contourplot. 'none' turns lines off.
%                                      Default: '-'. If ncontours >50 this option will be forced to 'none'.
%            peaks:             (logical) if true, displays predefined peak positions
%                                      Default: false.
%
%
% Examples:
%       	1. eemreview(EEMcor)
%           2. eemreview(EEMcor,'peaks',true)
%           3. eemreview(EEMcor,'nContours',25,'LineStyle','none')
%           4. eemreview(EEMcor,'nContours',15,'LineStyle','none','samples',[1 10 50])
%           5. eemreview(EEMcor,'nContours',15,'LineStyle','none','samples',[1 10 50],'title','treatment')
%           6. eemreview(EEMcor,'eem','Xnotscaled')
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013,
%     DOI:10.1039/c3ay41160e.
%
% eemreview: Copyright (C) 2019 Urban J. Wuensch
% Chalmers University of Technology
% Sven Hultins Gata 6
% 41296 Gothenburg
% Sweden
% wuensch@chalmers.se
% $ Version 0.1.0 $ March 2019 $ First Release

if nargin==0
    help eemreview
    return
end

%% Input Parser
params = inputParser;
params.addParameter('eem', 'X', @ischar);
params.addParameter('title', 'filelist', @ischar);
params.addParameter('samples', [], @isnumeric);
params.addParameter('nContours', 15, @isnumeric);
params.addParameter('LineStyle', '-', @ischar);
params.addParameter('peaks', false, @islogical);
params.addParameter('hold', false, @islogical);
params.parse(varargin{:});

Xname = params.Results.eem;
ContLine  = params.Results.LineStyle;
titlefield  = params.Results.title;
ncontours = params.Results.nContours;
samples   = params.Results.samples;
peakSwitch   = params.Results.peaks;
holdSwitch   = params.Results.hold;
if isempty(samples)
    samples=1:DS.nSample;
end

if peakSwitch
    ContLine='none';
end

if ncontours>=50
    ContLine='none';
end

%% Welcome messages
disp(' ')
disp(' ')
disp('eemreview.m')
disp('----------')
disp('Plot selected EEMs and Em and Ex scans')
disp('Close the figure to exit this function.')
disp('----------')
disp('For each sample, you have the following options:')
disp('1. Press ''Spectra'' to view both Ex and Em spectrum at the selected wavelengths.')
disp('2. Use the slider to adjust colorscale for the contour plot.')
disp('3. Missed something? Press ''Previous sample'' to go back one sample.')
disp('4. Press ''next'' to continue with the next sample. By default, hitting the ''spacebar'' will plot the next sample')
disp('----------')



% EEM extraction
if ~any(strfind(Xname,'Model'))&&~any(strfind(Xname,'resi'))
    X=DS.(Xname);
    %X=getfield(DS,Xname); % Field extraction
    X=X(samples,:,:); %Sample selection if no special cases are selected
elseif any(strfind(Xname,'Model'))&&~any(strfind(Xname,'resi')) % In case X is model, modeled data will be caluclated
    try
        X=DS.(Xname(1:6));
    catch
        error(['No field with the name ' Xname(1:6) ' found.'])
    end
    X=nmodel(X);
    X=X(samples,:,:);
elseif any(strfind(Xname,'Model'))&&any(strfind(Xname,'resi'))  % In case X is model and residuals are requested
    try
        X=DS.(Xname(1:6));
    catch
        error(['No field with the name ' Xname(1:6) ' found.'])
    end
    X=nmodel(X);
    X=X(samples,:,:);
    X=DS.X(samples,:,:)-X;
end


try
    filelist=DS.(titlefield);
    filelist=filelist(samples,1);
catch
%     warning('Could not field ''filelist'' or user-specified field given as option ''title''.')
    filelist=cellstr(num2str((1:DS.nSample)'));
end
% Figure business (close old eemreview figures, they might cause trouble)
figHandles = get(0,'Children');
for n=1:numel(figHandles);if strfind(figHandles(n).Name,'eemreview');close(figHandles(n));end;end

%% Define figure and axes
fig1=dreemfig;
set(fig1, 'units','normalized','pos',[0.2453    0.1611    0.5448    0.5000]);
if ~isfield(DS,'name')
    set(fig1,'Name','eemreview')
else
    try set(fig1,'Name',['eemreview - dataset:',DS.name{1}])
    catch
        set(fig1,'Name',['eemreview - dataset:',DS.name])
    end
end

pan(1) = uipanel(fig1,'Position',[0.001 0.16    0.60    0.8 ],...
    'BackgroundColor','w','BorderType','none'); % EEM axis
pan(2) = uipanel(fig1,'Position',[0.61  0.02    0.39    0.95],...
    'BackgroundColor','w','BorderType','none'); % 2D figures
pan(3) = uipanel(fig1,'Position',[0.001 0.0200  0.3     0.15],...
    'BackgroundColor','w','BorderType','none'); % Sample navigation
pan(4) = uipanel(fig1,'Position',[0.45   0.0200  0.15     0.15],...
    'BackgroundColor','w','BorderType','none'); % Spectra button
pan(5) = uipanel(fig1,'Position',[0.32   0.0200  0.15     0.15],...
    'BackgroundColor','w','BorderType','none'); % Control toggles

ax(1) = axes(pan(1),'pos',[0.15    0.15    0.6    0.75]);
ax(2) = axes(pan(2),'pos',[0.15    0.2    0.8    0.3]);
ax(3) = axes(pan(2),'pos',[0.15    0.65    0.8    0.3]);
ax(2).Title.String='Ex.spectra';
ax(3).Title.String='Em.spectra';
for n=1:3;set(ax(n),'Box','on');end

% Define uicontrol elements
positions=[0.2,0.68,0.75,0.3;... % popupmenu
    0.02,0.01,0.9,0.05;...       % RayRam text field
    0.2,0.2,0.3,0.4;...          % Previous
    0.1,0.4,0.8,0.4;...          % Spectra
    0.55,0.2,0.4,0.4;...         % Next
    0.9157 0.15 0.0661 0.75;...  % CLim slider
    0.05 0.6 0.9 0.2;...  % Checkbox 1
    0.05 0.2 0.9 0.2;...  % Checkbox 2
    ];

huic(1)=uicontrol(pan(3),'Style','popupmenu',...
    'Position',[100 20 100 20],...
    'String',filelist,...
    'Units','normalized',...
    'Position', positions(1,:));
huic(3) = uicontrol(pan(3),...
    'Style', 'pushbutton',...
    'String','Prev.',...
    'Units','normalized',...
    'Position', positions(3,:));
huic(5) = uicontrol(pan(3),...
    'Style', 'pushbutton',...
    'String','Next',...
    'FontWeight','bold',...
    'Units','normalized',...
    'Position', positions(5,:));
huic(4) = uicontrol(pan(4),...
    'Style', 'pushbutton',...
    'String','Spectra',...
    'FontWeight','bold',...
    'Units','normalized',...
    'Position', positions(4,:));

huic(2) = uicontrol(pan(2),...
    'Style','text',...
    'Position',[100 20 100 20],...
    'String','Black: Raman, red: Rayleigh scatter ',...
    'Units','normalized',...
    'Position', positions(2,:),...
    'BackgroundColor',[1 1 1],...
    'HorizontalAlignment','center');
huic(6) = uicontrol(pan(1),...
    'Style', 'slider',...
    'String','CLim',...
    'Min',0,'Max',1,'Value',0.5,...
    'Units','normalized',...
    'Position', positions(6,:));

huic(7) = uicontrol(pan(5),...
    'Style', 'checkbox',...
    'String','hold spectra',...
    'Min',false,'Max',true,'Value',holdSwitch,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position', positions(7,:));
huic(8) = uicontrol(pan(5),...
    'Style', 'checkbox',...
    'String','show peaks',...
    'Min',false,'Max',true,'Value',peakSwitch,...
    'Units','normalized',...
    'BackgroundColor',[1 1 1],...
    'Position', positions(8,:),...
     'Callback', {@overlaypeaks,DS,ax(1)});

for n=1:numel(huic)
    huic(n).FontSize=9;
    huic(n).FontName='Arial';
end

%% Plotting
sample_i=1;
while sample_i<=numel(samples)
    try
        if ~ishandle(fig1); return; end % Ends function when plot is closed by user
        ploteem([],[],DS,X,sample_i,filelist,samples,ncontours,ContLine,ax,huic)
        % uicontrol objects
        set(huic(1),'Callback',{@pulldownselect,DS,sample_i,ax,huic(1),sample_i},...
            'Value',sample_i);
        set(huic(3),'Callback',{@prevsample,DS,sample_i,ax});
        set(huic(5),'Callback',{@nextsample,DS,sample_i,ax});
        set(huic(4),'Callback',{@init_exemmview,DS,X,sample_i,ax,huic});
        
        ConLim=get(ax(1),'CLim');
        set(huic(6),'Min',ConLim(1)./2,'Max',ConLim(2)*2,'Value',ConLim(2),...
            'Callback', {@ContourCLim,ConLim,ax,huic});            
        
        
        %wait for uiinput
        datacursormode on
        uiwait(fig1)
        if ~ishandle(fig1); return; end % Ends function when plot is closed by user
        counter=str2double(ax(1).Children(end).DisplayName); % Callback functions cannot return values. User decision stored in DisplayName of EEM.
        if counter==-1&&sample_i==1
            %warning('You''ve reached the first sample. Try ''Next sample''.')
        else
            sample_i=sample_i+counter;
        end
    catch ME
        disp(' ');rethrow(ME); % Repeats last error message
    end
end
disp(' ')
disp('Done.')
end

function pulldownselect(~,~,~,~,ax,pdhdl,sample_i)
tit={'','Excitation','Emission'};
for n=2:3
    cla(ax(n))
    legend(ax(n),'off')
    title(ax(n),tit{n})
end
ax(1).Children(end).DisplayName=num2str(pdhdl.Value-sample_i);
uiresume
end


function nextsample(~,~,~,~,ax)
tit={'','Excitation','Emission'};
for n=2:3
    cla(ax(n))
    legend(ax(n),'off')
    title(ax(n),tit{n})
end
ax(1).Children(end).DisplayName='1';
uiresume
end

function prevsample(~,~,~,~,ax)
tit={'','Excitation','Emission'};
for n=2:3
    cla(ax(n))
    legend(ax(n),'off')
    title(ax(n),tit{n})
end
ax(1).Children(end).DisplayName='-1';
uiresume
end


function ContourCLim(source,~,ConLim,ax,huic)
val = source.Value;
try
    set(ax(1),'CLim',[ConLim(1) val])
catch
    set(ax(1),'CLim',[ConLim(1) ConLim(2)])
    warning('Could not set color limit to your preference.')
end
uicontrol(huic(5))
end

function ploteem(~,~,DS,X,i,filelist,samples,ncontours,ContLine,ax,huic)
% Data extraction
eem=(squeeze(X(i,:,:)));
% Plotting
warning('OFF','MATLAB:contourf:EmptyV6OutputArgument')
contourf(ax(1),DS.Ex,DS.Em,eem,ncontours,'Color',[.3 .3 .3],'LineStyle',ContLine,'DisplayName',num2str(i));
caxis([min(eem(~isnan(eem))) max(eem(~isnan(eem)))])
warning('ON','MATLAB:contourf:EmptyV6OutputArgument')
% Formatting
if isempty(findall(gcf,'type','ColorBar'))
    c=colorbar(ax(1));
    c.Position=[0.78 0.15 0.0470 0.75];
end
xlabel(ax(1),'Excitation (nm)');ylabel(ax(1),'Emission (nm)');
title(ax(1),[num2str(i),'/',num2str(numel(samples)),': ',char(filelist(i))], 'Interpreter', 'none');

overlaypeaks(huic(8),[],DS,ax(1))

uicontrol(huic(5))
end

function init_exemmview(~,~,DS,X,i,ax,huic)
[x,y]=ginput(1);

for n=2:3
    finddelete_scatterline(ax(n))
end

for n=2:3
    if huic(7).Value
        hold(ax(n),'on')
    else
        hold(ax(n),'off')
    end
end

exview(DS,X,i,y,ax(2))
emview(DS,X,i,x,ax(3))


axes(ax(1))
uicontrol(huic(5))
end


function [] = exview(DS,X,samSel,Em,ax)

Emsel = mindist(DS.Em,Em);
Em=DS.Em(Emsel,1);
h=plot(ax,DS.Ex,(squeeze(X(samSel,Emsel,:))),'Color','k','LineWidth',1.5,'DisplayName',['Em ',num2str(Em)]);
[numplot,~,exnam]=findexorem(ax);
h.Color=colmaker(numplot);
hold(ax,'on')
axis(ax,'tight')
physicatter_ex(Em,ax)
xlim(ax,[min(DS.Ex) max(DS.Ex)])

ylabel(ax,'Fl. intensity')
title(ax,['Ex. spectrum at em. ', num2str(DS.Em(Emsel)),' nm']);
if numplot>1
    hleg=leg(numplot,colmaker(1:numplot),ax);
    legend(hleg,fliplr(exnam))
end
end

function [] = emview(DS,X,samSel,Ex,ax)

Exsel=mindist(DS.Ex,Ex);
Ex=DS.Ex(Exsel,1);
h=plot(ax,DS.Em,(squeeze(X(samSel,:,Exsel))),'Color','k','LineWidth',1.5,'DisplayName',['Ex ',num2str(Ex)]);
[numplot,emnam,~]=findexorem(ax);
h.Color=colmaker(numplot);
hold(ax,'on')
axis(ax,'tight')
physicatter_em(Ex,ax)
xlabel(ax,'Emission (nm)')
xlim(ax,[min(DS.Em) max(DS.Em)])
ylabel(ax,'Fl. intensity')
title(ax,['Em. spectrum at ex. ', num2str(DS.Ex(Exsel)),' nm']);
if numplot>1
    hleg=leg(numplot,colmaker(1:numplot),ax);
    legend(hleg,fliplr(emnam))
end

end


function [idx,distance] = mindist( vec,value)
[distance,idx]=min(abs(vec-value));
end


function physicatter_ex(wl,ax)
%wl: Em
cylim=get(ax,'YLim');
cxlim=get(ax,'XLim');

y=[0 cylim(2)];
ray1st=[wl wl];
ray2nd=ray1st/2;

ram1=(1*10^7*((1*10^7)/(wl)+3382)^-1);
ram2=(1*10^7*((1*10^7)/(wl/2)+3382)^-1);

ram1st=[ram1 ram1];
ram2nd=[ram2 ram2];

line(ax,ram1st,y,...
    'DisplayName','Raman 1st','Color','k','LineWidth',0.5);
line(ax,ram2nd,y,...
    'DisplayName','Raman 2nd','Color','k','LineWidth',0.5);
line(ax,ray1st,y,...
    'DisplayName','Rayleigh 1st','Color','r','LineWidth',0.5);
line(ax,ray2nd,y,...
    'DisplayName','Rayleigh 2nd','Color','r','LineWidth',0.5);
set(ax,'XLim',cxlim);
set(ax,'YLim',cylim);

end

function physicatter_em(wl,ax)
cylim=get(ax,'YLim');
cxlim=get(ax,'XLim');

y=[0 cylim(2)];
ram1=1*10^7*((1*10^7)/(wl)-3382)^-1;
ram2=(1*10^7*((1*10^7)/(wl)-3382)^-1)*2;

ram1st=[ram1 ram1];
ram2nd=[ram2 ram2];
ray1st=[wl wl];
ray2nd=[wl*2 wl*2];

line(ax,ram1st,y,...
    'DisplayName','Raman 1st','Color','k','LineWidth',0.5);
line(ax,ram2nd,y,...
    'DisplayName','Raman 2nd','Color','k','LineWidth',0.5);
line(ax,ray1st,y,...
    'DisplayName','Rayleigh 1st','Color','r','LineWidth',0.5);
line(ax,ray2nd,y,...
    'DisplayName','Rayleigh 2nd','Color','r','LineWidth',0.5);

set(ax,'XLim',cxlim);
set(ax,'YLim',cylim);
end

function finddelete_scatterline(ax)
hc=get(ax,'Children');
for n=1:numel(hc)
    if ~isempty(strfind(hc(n).DisplayName,'Raman'))||~isempty(strfind(hc(n).DisplayName,'Rayleigh'))
        %contains(hc(n).DisplayName,{'Raman','Rayleigh'}) % If only everyone would use the latest software
        delete(hc(n))
    end
end
end


function overlaypeaks(source,~,data,ax)

% Formatting
lw=1.5;
OnOff=source.Value;
col=[0 0 0];
marker='+';
markersize=20;
% Definition of peaks
peaks(1).name=cellstr('B');
peaks(1).Em=num2cell(305);
peaks(1).Ex=num2cell(275);
peaks(2).name=cellstr('T');
peaks(2).Em=num2cell(340);
peaks(2).Ex=num2cell(275);
peaks(3).name=cellstr('A');
peaks(3).Em=num2cell(400:460);
peaks(3).Ex=num2cell(260);
peaks(4).name=cellstr('M');
peaks(4).Em=num2cell(370:410);
peaks(4).Ex=num2cell(290:310);
peaks(5).name=cellstr('C');
peaks(5).Em=num2cell(420:460);
peaks(5).Ex=num2cell(320:360);
peaks(6).name=cellstr('D');
peaks(6).Em=num2cell(509);
peaks(6).Ex=num2cell(390);
peaks(7).name=cellstr('N');
peaks(7).Em=num2cell(370);
peaks(7).Ex=num2cell(280);
if OnOff
    hold(ax,'on')
    for nPeak=1:size(peaks,2)
        if numel(peaks(nPeak).Em)==1&&numel(peaks(nPeak).Ex)==1
            plot(ax,data.Ex(mindist(data.Ex,peaks(nPeak).Ex{1}),1),data.Em(mindist(data.Em,peaks(nPeak).Em{1}),1),marker,'MarkerSize',markersize,'LineWidth',lw,'Color',col,'DisplayName','PeaksOverlayed'),hold on
            Labels = peaks(nPeak).name;    %'
            text(ax,data.Ex(mindist(data.Ex,peaks(nPeak).Ex{1}),1)+0.1,data.Em(mindist(data.Em,peaks(nPeak).Em{1}),1)+0.1, Labels, 'horizontal','left', 'vertical','bottom','Color',col,'DisplayName','PeaksOverlayed');
        elseif numel(peaks(nPeak).Ex)==1&&numel(peaks(nPeak).Em)~=1
            tempEm=data.Em(mindist(data.Em,peaks(nPeak).Em{1}):mindist(data.Em,peaks(nPeak).Em{end}),1);
            plot(ax,repmat(data.Ex(mindist(data.Ex,peaks(nPeak).Ex{1}),1),numel(tempEm)),tempEm,'LineWidth',lw,'Color',col,'DisplayName','PeaksOverlayed')
            Labels = peaks(nPeak).name;    %'
            text(ax,data.Ex(mindist(data.Ex,peaks(nPeak).Ex{1}),1)+0.1,data.Em(mindist(data.Em,peaks(nPeak).Em{1}),1)+0.1, Labels, 'horizontal','left', 'vertical','bottom','Color',col,'DisplayName','PeaksOverlayed');
        else
            x1=data.Ex(mindist(data.Ex,peaks(nPeak).Ex{1}),1);
            x2=data.Ex(mindist(data.Ex,peaks(nPeak).Ex{end}),1);
            y1=data.Em(mindist(data.Em,peaks(nPeak).Em{1}),1);
            y2=data.Em(mindist(data.Em,peaks(nPeak).Em{end}),1);
            x = [x1, x2, x2, x1, x1];
            y = [y1, y1, y2, y2, y1];
            plot(ax,x, y, 'b-', 'LineWidth', lw,'Color',col,'DisplayName','PeaksOverlayed');
            Labels = peaks(nPeak).name;
            text(ax,data.Ex(mindist(data.Ex,peaks(nPeak).Ex{1}),1)+0.1,data.Em(mindist(data.Em,peaks(nPeak).Em{1}),1)+0.1, Labels, 'horizontal','left', 'vertical','bottom','Color',col,'DisplayName','PeaksOverlayed');
        end
    end
    hold(ax,'off')
else
    hc=get(ax,'Children');
    for n=1:numel(hc)
        if ~isempty(strfind(hc(n).DisplayName,'PeaksOverlayed'))||~isempty(strfind(hc(n).DisplayName,'PeaksOverlayed'))
            %contains(hc(n).DisplayName,{'Raman','Rayleigh'}) % If only everyone would use the latest software
            delete(hc(n))
        end
    end
end
end


function [num,ex,em] = findexorem(ax)
hc=get(ax,'Children');
num=0;
ex=cell(1);
em=cell(1);
for n=1:numel(hc)
    if ~isempty(strfind(hc(n).DisplayName,'Ex'))||~isempty(strfind(hc(n).DisplayName,'Em'))
    %contains(hc(n).DisplayName,{'Ex','Em'}) % For 2017b and higher
        num=num+1;
        if ~isempty(strfind(hc(n).DisplayName,'Ex'))
            ex{num}=hc(n).DisplayName(4:end);
        elseif ~isempty(strfind(hc(n).DisplayName,'Em'))
            em{num}=hc(n).DisplayName(4:end);
        end
    end
end
end

function cols = colmaker(i)
collist= [0 0 0;...
    0.368627450980392,0.309803921568627,0.635294117647059;...
    0.197386603025884,0.512851710442123,0.740257214831838;...
    0.353920722469962,0.729472517932409,0.656177840632872;...
    0.597784406620190,0.840777034550949,0.644490274269565;...
    0.841701969675001,0.940885130292891,0.609611616709332;...
    0.931851609124017,0.957193786605551,0.661676168589149;...
    0.976267475434952,0.928364732010819,0.658313136694680;...
    0.995493878994122,0.822714483114813,0.482765825462557;...
    0.987654277338283,0.615414278604657,0.339114037801447;...
    0.943815034723900,0.390980258870390,0.266846392952373;...
    0.820562112287474,0.223926402795306,0.309367864350849;...
    0.619607843137255,0.00392156862745098,0.258823529411765];
if numel(i)==1
    try
        cols=collist(i,:);
    catch
        cols=collist(1,:);
        warning('No more unique Colors available. Plotting ''k''.')
    end
else
    cols=collist(i,:);
end
end

function [ h ] = leg( numPlots,col,ax)
h = gobjects(numPlots, 1);
for n=1:numPlots
    hold on;
    h(n) = line(ax,NaN,NaN,'Marker','o','MarkerSize',5,'MarkerFaceColor',col(n,:),'MarkerEdgeColor','k','LineStyle','none');
end
end