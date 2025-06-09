function fingerprint(data,f,varargin)
%
% <strong>Syntax</strong>
%   <strong>fingerprint</strong>(data,f,varargin)
%   <strong>fingerprint</strong>(data,f,subplots,viewscatter)
%
% <a href="matlab: doc fingerprint">help for fingerprint</a> <- click on the link

%Plot a PARAFAC model in fingerprint mode (contour plots).
%
%USEAGE: 
%       fingerprint(data,f,subplots)
%
%INPUTS: 
%        data: data structure containing a PARAFAC model in data. Modelf
%           f: Number of components in the model to be fitted.
%    subplots:(optional) plot arrangement e.g. [1,5] for 1 row x 5 col.
%    viewscatter: (optional). If true, scatter excision details will be
%    superimposed.
%
%OUTPUTS:
%    A figure with f plots will be produced, where each plot shows the 
%    spectra (in contours) for one of the components. The plots are shown
%    in the order that the components were resolved by PARAFAC
%     i.e. plot 1 -> component 1.

%Examples
%   fingerprint(LSmodel6,6)
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, (yr), 
%     DOI:10.1039/c3ay41160e. 
%
% fingerprint; Copyright (C) 2013 Kathleen R. Murphy
% The University of New South Wales
% Dept Civil and Environmental Engineering
% Water Research Center
% UNSW 2052
% Sydney
% krm@unsw.edu.au
%
% $ Version 0.2.0 $ March 2019 $ minor improvements.
% $ Version 0.1.0 $ September 2013 $ First Release

narginchk(2,4)
try
    M = data.(['Model' int2str(f)]);
catch
    error('A model with the stated number of components does not exist in this dataset.')
end
[~,B,C]=fac2let(M);
if ~isfield(data,'nEm'),data.nEm=size(B,1);end
if ~isfield(data,'nEx'),data.nEx=size(C,1);end

rc=      [1 1;1 2;1 3;2 2;  2 3;  2 3;  2 4;  2 4;  3 3;   2 5;     3 4;       3 4;     4 4;      4 4];
ylabpos={{1},{1},{1},{1,3},{1,4},{1,4},{1,5},{1,5},{1,4,7},{1,4,7},{1,4,7,11},{1,4,7,11},{1,5,9,12},{1,5,9,12}};
xlabpos={{1},{2},{2},{3:4},{3:6},{4:6},{5:8},{5:8},{7:9},  {7:10}, {9:12},    {9:12},     {12:16}, {12:16}};

if nargin==3
    rcin=varargin{1};
    rc(f,1)=rcin(1);
    rc(f,2)=rcin(2);
end
showscatter=false;
if nargin>3
    showscatter=varargin{2};
    if ~islogical(showscatter)
        error('Input to ''viewscatter'' must be logical.')
    end
end

if showscatter
    try
        data.Smooth.description;
        
        if numel(data.Smooth)>1
            warning('multiple smoothing settings used. Scatter for the first set of settings is shown.')
            data.Smooth=data.Smooth(1);
        end
        
    catch
        
        try
            sparam=data.Smooth;
            sparam=strsplit(sparam,',');
            for n=1:numel(sparam)
                sparam{n}=str2num(sparam{n});
            end
            opt=handlescatter('options');
            opt.ray1(1)=sparam{1}(2);opt.ray1(2)=sparam{1}(1);
            opt.ram1(1)=sparam{2}(2);opt.ram1(2)=sparam{2}(1);
            opt.ray2(1)=sparam{3}(2);opt.ray2(2)=sparam{3}(1);
            opt.ram2(1)=sparam{4}(2);opt.ram2(2)=sparam{4}(1);
            
            if ~all(opt.ray1==0)
                opt.cutout(1)=1;
            end
            if ~all(opt.ram1==0)
                opt.cutout(2)=1;
            end
            if ~all(opt.ray2==0)
                opt.cutout(3)=1;
            end
            if ~all(opt.ram2==0)
                opt.cutout(4)=1;
            end
            data.Smooth=opt;
        catch
            showscatter=false;
            warning('Something went wrong while parsing the output of smootheem. No scatter lines will be shown.')
        end
    end
end

hf=dreemfig;
set(hf,'units','normalized','pos',[0.2 0.25 0.5 0.4],'InvertHardcopy','off')
set(hf,'Name',['Contour plot for ' num2str(f) '-comp. model']);
for i=1:f
    subplot(rc(f,1),rc(f,2),i)
    Comp=reshape((krb(C(:,i),B(:,i))'),[1 data.nEm data.nEx]);
    contourf(data.Ex,data.Em,(squeeze(Comp(1,:,:))),'DisplayName','contourf');
    v=axis;
    handle=title(['Comp ' int2str(i)]);
    set(handle,'Position',[0.9*v(2) 1.05*v(3) 1],'FontWeight','bold','color',[1 1 1]);
    if ismember(i,cell2mat(ylabpos{f}))
        ylabel('Em. (nm)')
    end
    if ismember(i,cell2mat(xlabpos{f}))
        xlabel('Ex. (nm)')
    end
    if showscatter
        
        
        
        ram1=@(x) 1*10^7*((1*10^7)/(x)-3382)^-1;
        ram2=@(x) (1*10^7*((1*10^7)/(x)-3382)^-1)*2;
        
        ray1=@(x) x;
        ray2=@(x) x*2;
        
        
        if data.Smooth.cutout(2)
            line([min(data.Ex) max(data.Ex)],[ram1(min(data.Ex)) ram1(max(data.Ex))],'Color','r')
            line([min(data.Ex) max(data.Ex)],[ram1(min(data.Ex))+data.Smooth.ram1(2) ram1(max(data.Ex))+data.Smooth.ram1(2)],'Color','k')
            line([min(data.Ex) max(data.Ex)],[ram1(min(data.Ex))-data.Smooth.ram1(1) ram1(max(data.Ex))-data.Smooth.ram1(1)],'Color','k')
        end
        if data.Smooth.cutout(4)
            line([min(data.Ex) max(data.Ex)],[ram2(min(data.Ex)) ram2(max(data.Ex))],'Color','r')
            line([min(data.Ex) max(data.Ex)],[ram2(min(data.Ex))+data.Smooth.ram2(2) ram2(max(data.Ex))+data.Smooth.ram2(2)],'Color','k')
            line([min(data.Ex) max(data.Ex)],[ram2(min(data.Ex))-data.Smooth.ram2(1) ram2(max(data.Ex))-data.Smooth.ram2(1)],'Color','k')
        end
        
        
        if data.Smooth.cutout(1)
            line([min(data.Ex) max(data.Ex)],[ray1(min(data.Ex)) ray1(max(data.Ex))],'Color','r')
            line([min(data.Ex) max(data.Ex)],[ray1(min(data.Ex))+data.Smooth.ray1(2) ray1(max(data.Ex))+data.Smooth.ray1(2)],'Color','k')
            line([min(data.Ex) max(data.Ex)],[ray1(min(data.Ex))-data.Smooth.ray1(1) ray1(max(data.Ex))-data.Smooth.ray1(1)],'Color','k')
        end
        if data.Smooth.cutout(3)
            line([min(data.Ex) max(data.Ex)],[ray2(min(data.Ex)) ray2(max(data.Ex))],'Color','r')
            line([min(data.Ex) max(data.Ex)],[ray2(min(data.Ex))+data.Smooth.ray2(2) ray2(max(data.Ex))+data.Smooth.ray2(2)],'Color','k')
            line([min(data.Ex) max(data.Ex)],[ray2(min(data.Ex))-data.Smooth.ray2(1) ray2(max(data.Ex))-data.Smooth.ray2(1)],'Color','k')
        end
        
        
        ylim([min(data.Em) max(data.Em)])
        xlim([min(data.Ex) max(data.Ex)])
        

        
    end
end
dreemfig(hf);
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@parse4datacursor,data});
datacursormode on