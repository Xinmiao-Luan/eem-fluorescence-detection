function out=diffeem(DS,varargin)
%
% <strong>Syntax</strong>
%   out=<strong>diffeem</strong>(DS,Name,Value)
%
% <a href="matlab: doc diffeem">help for diffeem</a> <- click on the link

% Calculate & visualize the difference between a reference EEM and other samples.
%
% USEAGE:
%          out=diffeem(DS,varargin)
%
% INPUTS
%            DS:                drEEM-format dataset
%            (optional):        Parameter name followed by option.
%            ref:               (numeric): Index of reference sample
%                                          Default: first sample.
%            samples:           (numeric): Index of samples to be subtracted from reference.
%                                          Default: all but the first sample.
%            plot:              (logical): If true, plots will be shown.
%                                          Default: true.
%            title:             (char):     Name of cell array containing sample names
%                                          Default: 'filelist'.
%            mode:              (char):    'absolute' or 'percentage'. Difference will be calculated 
%                                          as absolute of relative difference.
%                                          Default: 'absolute'.
%
% OUTPUTS
%            out:               drEEM-format dataset. Containing:
%                               out.Xsubtract: difference EEMs. The first consists of zeros, representing the 
%                                              reference subtraction. The remaining EEMs come in the sequence
%                                              supplied in the option "samples".
%                               out.Xsubtr_ref: Name of the reference, as supplied with the option "title".
%                               out.Xsubtr_titlelist: Names of the
%                               difference EEMs.
%
%           To visualize the function output, type the following
%               eemreview(out,'eem','Xsubtract','title','Xsubtr_titlelist')
%
% Examples:
%       	1. out=diffeem(DS)
%           2. out=diffeem(DS,'mode','percentage')
%           3. out=diffeem(DS,'ref',10,'samples',setdiff(1:DS.nSample,10))
%              (select 10th sample as reference and subtract from all other EEMs)
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
% $ Version 0.1.0 $ April 2019 $ First Release

%% Parse inputs
params = inputParser;
params.addParameter('ref', 1, @isnumeric);
params.addParameter('samples', 2:DS.nSample, @isnumeric);
params.addParameter('plot',true, @islogical);
params.addParameter('title', 'filelist', @ischar);
params.addParameter('mode', 'absolute', @diagmode);

params.parse(varargin{:});
init   = params.Results.ref;
subtr  = params.Results.samples;
doplot = params.Results.plot;
names  = params.Results.title;
fmode  = params.Results.mode;


%% Difference between EEMs in the original scale
adiff=(DS.X(subtr,:,:)-DS.X(init,:,:));           % absolute
pdiff=((DS.X(subtr,:,:)-DS.X(init,:,:))...
    ./DS.X(init,:,:))*100;                               % percentage

switch fmode
    case 'absolute'
        Y=adiff;
    case 'percentage'
        Y=pdiff;
    otherwise
        error('Input to ''mode'' not understood. Refer to help.')
end
%% Assign output
out=DS;
out.Xsubtract=[zeros(1,DS.nEm,DS.nEx);Y];

sname=DS.(names)(init);
if iscell(sname)
    sname=sname{1};
end
if ~ischar(sname)
    try
        sname=num2str(sname);
    catch
        sname='REF';
    end
end


out.Xsubtr_ref=['REF: ',sname];
out.Xsubtr_titlelist=[DS.(names)(init);DS.(names)(subtr)];

%% Plotting
if doplot
    dpeaks=pickdpeaks(out);
    dpeaks(1,:)={nan};
    fig1=dreemfig;
    set(fig1,'units','normalized',...
        'Name','diffeem: Difference between reference and samples',...
        'pos',[0.1589    0.1806    0.6026    0.3741])
    ax(1)=subplot(1,2,1);
    ax(2)=subplot(1,2,2);

    plot(ax(1),table2array(dpeaks),'LineWidth',0.5,'Marker','.','MarkerSize',15,'LineStyle','-.')
    axis(ax(1),'tight')
    if height(dpeaks)<50
        set(ax(1),'XTick',1:height(dpeaks))
        set(ax(1),'XTickLabel',out.Xsubtr_titlelist,'XTickLabelRotation',45)
        lab=get(ax(1),'XTickLabel');
        lab=cellstr(lab);
        lab{1}='REF';
        set(ax(1),'XTickLabel',lab)
    else
        set(ax(1),'XTick',floor(1:10:height(dpeaks)-1))
        set(ax(1),'XTickLabel',out.Xsubtr_titlelist(floor(1:10:height(dpeaks)-1)),'XTickLabelRotation',45)
    end
    
    switch fmode
        case 'absolute'
            title(ax(1),'Peaks: Absolute difference')
            ylabel(ax(1),'Absolute difference (original unit)')
        case 'percentage'
            title(ax(1),'Peaks: Relative difference')
            ylabel(ax(1),'% Difference from reference')
    end
    line(ax(1),get(ax(1),'XLim'),[0 0],'LineWidth',2,'LineStyle','--','Color','k')
    try
        legend(ax(1),dpeaks.Properties.VariableNames,'NumColumns',4,'location','best')
    catch
        legend(ax(1),dpeaks.Properties.VariableNames,'location','best')
    end
    
    disp('Press any key to continue or Ctrl + C to cancel.')
    for n=2:numel(subtr)+1
        switch fmode
            case 'percentage'
            llist=linspace(min(min(table2array(dpeaks(n,:))))-20, max(max(table2array(dpeaks(n,:))))+20,50);
            case 'absolute'
                mat=out.Xsubtract(n,:,:);mat=mat(:);
                llist=linspace(min(mat), max(mat),50);
        end
        [~,h]=contourf(ax(2),out.Ex,out.Em,squeeze(out.Xsubtract(n,:,:)),llist,'LineStyle','none');
        caxis([llist(1) llist(end)]);
        npos=numel(find(h.LevelList>0));
        nneg=numel(find(h.LevelList<0));
        cmap=mkclrs(npos,nneg);
        colormap(ax(2),cmap)
        c=colorbar(ax(2));
        
        
        sname=out.Xsubtr_titlelist(n);
        if iscell(sname)
            sname=sname{1};
        end
        if ischar(sname)
            sname=char(out.Xsubtr_titlelist(n));
        else
            try
                sname=num2str(out.Xsubtr_titlelist(n));
            catch
                sname=num2str(n);
            end
        end
        title({['EEM: ',num2str(n-1),'/',num2str(numel(subtr)),' :  '],[out.Xsubtr_ref,'-',sname]})
        switch fmode
            case 'absolute'
                ylabel(c,'Absolute difference (original unit)')
            case 'percentage'
                ylabel(c,'% Difference from reference')
        end
        hc=get(ax(1),'Children');
        for ii=1:numel(hc)
            if ~isempty(strfind(hc(ii).DisplayName,'progress'))
                delete(hc(ii));
            end
        end
        line(ax(1),[n n],get(ax(1),'YLim'),'DisplayName','progress',...
            'LineStyle','-','LineWidth',2,'Color','k')
        ylabel(ax(2),'Emission (nm)')
        xlabel(ax(2),'Excitation (nm)')
        dreemfig(fig1);
        pause
    end
    disp('Done.')
    close(fig1)
end
end


function dpeaks=pickdpeaks(out)
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


dpeaks=nan(size(out.Xsubtract,1),numel(peaks));
for i=1:size(out.Xsubtract,1)
    % Running through different scenarios
    % #1: Ex/Em pair
    for n=1:size(peaks,2)
        if numel(peaks(n).Em)==1&&numel(peaks(n).Ex)==1
            dpeaks(i,n)=out.('Xsubtract')(i,mindist(out.Em,peaks(n).Em{1}),mindist(out.Ex,peaks(n).Ex{1}));
            % #2: Specific Ex, but Em range
        elseif numel(peaks(n).Ex)==1&&numel(peaks(n).Em)~=1
            tempEm=out.('Xsubtract')(i,mindist(out.Em,peaks(n).Em{1}):mindist(out.Em,peaks(n).Em{end}),mindist(out.Ex,peaks(n).Ex{1}));
            dpeaks(i,n)=max(tempEm);
            % #3: Ex and Em range-peaks
        else
            x1=mindist(out.Ex,peaks(n).Ex{1});
            x2=mindist(out.Ex,peaks(n).Ex{end});
            y1=mindist(out.Em,peaks(n).Em{1});
            y2=mindist(out.Em,peaks(n).Em{end});
            dpeaks(i,n)=nanmax(nanmax(squeeze(out.('Xsubtract')(i,y1:y2,x1:x2))));
        end
    end
end
dpeaks=array2table(dpeaks,'VariableNames',[peaks.name]);
end

function [idx,distance] = mindist( vec,value)
[distance,idx]=min(abs(vec-value));
end


function cmap=mkclrs(npos,nneg)
nmap=[0 0.455 0.737;
0.01 0.46 0.74;
0.02 0.466 0.743;
0.03 0.471 0.745;
0.04 0.477 0.748;
0.051 0.482 0.751;
0.061 0.488 0.753;
0.071 0.493 0.756;
0.081 0.499 0.759;
0.091 0.504 0.761;
0.101 0.51 0.764;
0.111 0.515 0.766;
0.121 0.521 0.769;
0.131 0.526 0.772;
0.141 0.532 0.774;
0.152 0.537 0.777;
0.162 0.543 0.78;
0.172 0.549 0.782;
0.182 0.554 0.785;
0.192 0.56 0.788;
0.202 0.565 0.79;
0.212 0.571 0.793;
0.222 0.576 0.796;
0.232 0.582 0.798;
0.242 0.587 0.801;
0.253 0.593 0.804;
0.263 0.598 0.806;
0.273 0.604 0.809;
0.283 0.609 0.812;
0.293 0.615 0.814;
0.303 0.62 0.817;
0.313 0.626 0.82;
0.323 0.631 0.822;
0.333 0.637 0.825;
0.343 0.642 0.828;
0.354 0.648 0.83;
0.364 0.653 0.833;
0.374 0.659 0.835;
0.384 0.664 0.838;
0.394 0.67 0.841;
0.404 0.675 0.843;
0.414 0.681 0.846;
0.424 0.686 0.849;
0.434 0.692 0.851;
0.444 0.697 0.854;
0.455 0.703 0.857;
0.465 0.708 0.859;
0.475 0.714 0.862;
0.485 0.719 0.865;
0.495 0.725 0.867;
0.505 0.73 0.87;
0.515 0.736 0.873;
0.525 0.741 0.875;
0.535 0.747 0.878;
0.545 0.752 0.881;
0.556 0.758 0.883;
0.566 0.763 0.886;
0.576 0.769 0.889;
0.586 0.774 0.891;
0.596 0.78 0.894;
0.606 0.785 0.897;
0.616 0.791 0.899;
0.626 0.796 0.902;
0.636 0.802 0.904;
0.646 0.807 0.907;
0.657 0.813 0.91;
0.667 0.818 0.912;
0.677 0.824 0.915;
0.687 0.829 0.918;
0.697 0.835 0.92;
0.707 0.84 0.923;
0.717 0.846 0.926;
0.727 0.851 0.928;
0.737 0.857 0.931;
0.747 0.862 0.934;
0.758 0.868 0.936;
0.768 0.873 0.939;
0.778 0.879 0.942;
0.788 0.884 0.944;
0.798 0.89 0.947;
0.808 0.895 0.95;
0.818 0.901 0.952;
0.828 0.906 0.955;
0.838 0.912 0.958;
0.848 0.917 0.96;
0.859 0.923 0.963;
0.869 0.928 0.966;
0.879 0.934 0.968;
0.889 0.939 0.971;
0.899 0.945 0.973;
0.909 0.95 0.976;
0.919 0.956 0.979;
0.929 0.961 0.981;
0.939 0.967 0.984;
0.949 0.972 0.987;
0.96 0.978 0.989;
0.97 0.983 0.992;
0.98 0.989 0.995;
0.99 0.994 0.997;
1 1 1 ];
pmap=[1 1 1;
0.997 0.992 0.993;
0.994 0.984 0.985;
0.99 0.976 0.978;
0.987 0.968 0.971;
0.984 0.96 0.964;
0.981 0.952 0.956;
0.978 0.944 0.949;
0.975 0.936 0.942;
0.971 0.928 0.934;
0.968 0.92 0.927;
0.965 0.912 0.92;
0.962 0.904 0.913;
0.959 0.896 0.905;
0.956 0.888 0.898;
0.952 0.88 0.891;
0.949 0.872 0.883;
0.946 0.864 0.876;
0.943 0.856 0.869;
0.94 0.848 0.862;
0.937 0.84 0.854;
0.933 0.832 0.847;
0.93 0.824 0.84;
0.927 0.816 0.832;
0.924 0.808 0.825;
0.921 0.8 0.818;
0.918 0.792 0.81;
0.914 0.784 0.803;
0.911 0.776 0.796;
0.908 0.768 0.789;
0.905 0.76 0.781;
0.902 0.752 0.774;
0.899 0.744 0.767;
0.895 0.736 0.759;
0.892 0.728 0.752;
0.889 0.72 0.745;
0.886 0.712 0.738;
0.883 0.704 0.73;
0.88 0.696 0.723;
0.876 0.688 0.716;
0.873 0.68 0.708;
0.87 0.672 0.701;
0.867 0.664 0.694;
0.864 0.656 0.687;
0.861 0.648 0.679;
0.857 0.64 0.672;
0.854 0.632 0.665;
0.851 0.624 0.657;
0.848 0.616 0.65;
0.845 0.608 0.643;
0.842 0.6 0.636;
0.838 0.592 0.628;
0.835 0.584 0.621;
0.832 0.576 0.614;
0.829 0.568 0.606;
0.826 0.56 0.599;
0.823 0.552 0.592;
0.819 0.544 0.585;
0.816 0.536 0.577;
0.813 0.528 0.57;
0.81 0.52 0.563;
0.807 0.512 0.555;
0.804 0.504 0.548;
0.8 0.496 0.541;
0.797 0.488 0.534;
0.794 0.48 0.526;
0.791 0.472 0.519;
0.788 0.464 0.512;
0.785 0.456 0.504;
0.781 0.448 0.497;
0.778 0.44 0.49;
0.775 0.432 0.482;
0.772 0.424 0.475;
0.769 0.416 0.468;
0.766 0.408 0.461;
0.762 0.4 0.453;
0.759 0.392 0.446;
0.756 0.384 0.439;
0.753 0.376 0.431;
0.75 0.368 0.424;
0.747 0.36 0.417;
0.743 0.352 0.41;
0.74 0.344 0.402;
0.737 0.336 0.395;
0.734 0.328 0.388;
0.731 0.32 0.38;
0.727 0.312 0.373;
0.724 0.304 0.366;
0.721 0.296 0.359;
0.718 0.288 0.351;
0.715 0.28 0.344;
0.712 0.272 0.337;
0.708 0.264 0.329;
0.705 0.256 0.322;
0.702 0.248 0.315;
0.699 0.24 0.308;
0.696 0.232 0.3;
0.693 0.224 0.293;
0.689 0.216 0.286;
0.686 0.208 0.278];
cmap=[nmap(round(linspace(1,100,nneg)),:);pmap(round(linspace(1,100,npos)),:)];
end


function diagmode(inpt)
switch inpt
    case 'absolute'
        % ok
    case 'percentage'
        % ok
    otherwise
        error('Value for function option ''mode'' is invalid. Valid options: ''absolute'' or ''percentage''')
end
end