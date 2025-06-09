function [ picklist ] = pickpeaks( DS,varargin)
%
% <strong>Syntax</strong>
%   picklist = <strong>pickpeaks</strong>(DS,varargin)
%   picklist = <strong>pickpeaks</strong>(DS,'Name','Value')
%
% <a href="matlab: doc pickpeaks">help for pickpeaks</a> <- click on the link

% Pick fluorescence intensities at predefined peaks and calculate
% fluorescence indicies. For references to the publications, see bottom of
% this help section
%
% [ picklist ] = pickpeaks( DS, Name, Value )
%
% INPUTS
%            DS:                drEEM-format dataset
%           (optional):         Name,value: Parameter name followed by option
%           eem:                (char): name of field containing EEMs. Default: 'X'
%           plot:               (logigal): option to trigger or surpress plotting of results
%                               default: true (plot results)
%           details:            (logigal): if true, plots for extracted scans are shown. 
%                                          This is to diagnose smoothing and general applicability.
%                                          default: false
%
% OUTPUTS
%           picklist:           table with the results. Check table header for description
%EXAMPLES
% 1.   [ picklist ] = pickpeaks( DS,'plot',false)
% 2.   [ picklist ] = pickpeaks( DS,'plot',false,'details',true)
%
% NOTICE:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013,
%     DOI:10.1039/c3ay41160e.
%
% REFERENCES:
%
% 1. Peaks A,B,C,T,N. Peak D is not named but specified in this publication:
%    Coble, P.G., 2007. Marine optical biogeochemistry: The chemistry of ocean color.
%    Chem. Rev. 107, 402–418. https://doi.org/10.1021/cr050350+
% 2. Fluorescence index:
%    ----Ratio of em wavelengths at 470 nm and 520 nm, obtained at ex370---
%    Maie, N., Parish, K.J., Watanabe, A., Knicker, H., Benner, R.,
%    Abe, T., Kaiser, K., Jaffé, R., 2006. Chemical characteristics of dissolved organic
%    nitrogen in an oligotrophic subtropical coastal ecosystem. Geochim. Cosmochim. Acta 70,
%    4491–4506. https://doi.org/10.1016/j.gca.2006.06.1554
% 3. Freshness index:
%    ----em 380 divided by the em maximum between 420 and 435, obtained at ex310----
%    Parlanti, E., Wörz, K., Geoffroy, L., Lamotte, M., 2000. Dissolved organic matter
%    fluorescence spectroscopy as a tool to estimate biological activity in a coastal
%    zone submitted to anthropogenic inputs. Org. Geochem. 31, 1765–1781.
%    https://doi.org/10.1016/S0146-6380(00)00124-8
% 4. Humification index:
%    ----area under the em spectra 435–480 divided by the peak area 300–345 + 435–480, at ex254---
%    Parlanti, E., Wörz, K., Geoffroy, L., Lamotte, M., 2000. Dissolved organic matter fluorescence
%    spectroscopy as a tool to estimate biological activity in a coastal zone submitted to
%    anthropogenic inputs. Org. Geochem. 31, 1765–1781. https://doi.org/10.1016/S0146-6380(00)00124-8
% 5. Biological index:
%    ---Ex 310, Em 380 divided by 430.
%    Huguet, A., Vacher, L., Relexans, S., Saubusse, S., Froidefond, J.M., Parlanti, E., 2009.
%    Properties of fluorescent dissolved organic matter in the Gironde Estuary. Org. Geochem. 40,
%    706–719. https://doi.org/10.1016/j.orggeochem.2009.03.002
%
% All extracted scans are smoothed with a Savitzky-Golay filter (2nd order,
% 21 step window length) to mitigate the influence of noise. If your not
% sure if that works for your data, make sure to set 'details' to true in
% the function input and inspect the output!
%
%
% pickpeaks: Copyright (C) 2019 Urban J. Wuensch
% Chalmers University of Technology
% Sven Hultins Gata 6
% 41296 Gothenburg
% Sweden
% $ Version 0.1.0 $ March 2019 $ First Release
if nargin==0
    help pickPeaks
    return
end

params = inputParser;
params.addParameter('eem', 'X', @ischar);
params.addParameter('plot', true, @islogical);
params.addParameter('details', false, @islogical);

params.parse(varargin{:});
Xname = params.Results.eem;
plt = params.Results.plot;
diagn = params.Results.details;
%% Check for products

mv=ver;
tbx=false;
for n=1:numel(mv)
    if strcmp(mv(n).Name,'MATLAB')
        mver=mv(n).Version;
    end
    if strfind(mv(n).Name,'Curve Fitting')
        tbx=true;
    end
end

if ~tbx
    error('pickpeaks requires the Curve Fitting Toolbox.')
end

mver=str2double(mver);
if mver<=8.4
    plt=false;
    diagn=false;
    warning('Matlab is too old for some features regarding plots of this function. All plots surpressed to avoid errors.')
end

%% Definition of peaks: B, T, A, M, C, D, E, N.
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
% peaks(7).name=cellstr('E');
% peaks(7).Em=num2cell(521);
% peaks(7).Ex=num2cell(455);
peaks(7).name=cellstr('N');
peaks(7).Em=num2cell(370);
peaks(7).Ex=num2cell(280);

%% Define anon. functs

vec=@(x) x(:);

%% Change resolution of EEMs to 1nm in order to properly calculate and extract peaks and indicies
Exint=1;
Emint=1;
Ex_i=round(DS.Ex(1)):Exint:round(DS.Ex(end));
Em_i=(round(DS.Em(1)):Emint:round(DS.Em(end)))';
nEx_i=size(Ex_i,2);
nEm_i=size(Em_i,1);
Xi=zeros(DS.nSample,size(Em_i,1),size(Ex_i,2));

if mean(diff(DS.Ex))>5
    warning('Excitation increment>5nm. The performed interpolation might introduce significant biases!')
    pause(2)
end
if mean(diff(DS.Em))>10
    warning('Emission increment>10nm. The performed interpolation might introduce significant biases!')
    pause(2)
end


for i=1:DS.nSample
    eem=squeeze(DS.(Xname)(i,:,:));
    Xi(i,:,:) = interp2(DS.Ex,DS.Em,eem,Ex_i,Em_i);
end
DS.(Xname)=Xi;
DS.nEm=nEm_i;
DS.nEx=nEx_i;
DS.Ex=Ex_i';
DS.Em=Em_i;


%% Check if distance between fluorescence peaks and dataset Ex and Em is too big

if ~any(DS.Ex==254)
    warning('Humification index cannot be calculated due to dataset limitations')
    pause(2)
    HIX_excl=true;
else
    HIX_excl=false;
end

distEx=nan(1,numel(peaks));
distEm=nan(1,numel(peaks));
for n=1:numel(peaks)
    [~,distEx(n)] = mindist( DS.Ex,peaks(n).Ex{1});
    [~,distEm(n)] = mindist( DS.Em,peaks(n).Em{1});
end

if any(distEx>=5)
    warning('Distance between some peak definitions and Dataset-wavelengths >=5nm: Excitation');
    pause(2)
end
if any(distEm>=5)
    warning('Distance between some peak definitions and Dataset-wavelengths >=10nm: Emission');
    pause(2)
end


%% Peak picking
Cpeak=nan(DS.nSample,numel(peaks));
for i=1:DS.nSample
    % Running through different scenarios
    % #1: Ex/Em pair
    for n=1:size(peaks,2)
        if numel(peaks(n).Em)==1&&numel(peaks(n).Ex)==1
            Cpeak(i,n)=DS.(Xname)(i,mindist(DS.Em,peaks(n).Em{1}),mindist(DS.Ex,peaks(n).Ex{1}));
            % #2: Specific Ex, but Em range
        elseif numel(peaks(n).Ex)==1&&numel(peaks(n).Em)~=1
            tempEm=DS.(Xname)(i,mindist(DS.Em,peaks(n).Em{1}):mindist(DS.Em,peaks(n).Em{end}),mindist(DS.Ex,peaks(n).Ex{1}));
            Cpeak(i,n)=max(tempEm,[],'omitnan');
            % #3: Ex and Em range-peaks
        else
            x1=mindist(DS.Ex,peaks(n).Ex{1});
            x2=mindist(DS.Ex,peaks(n).Ex{end});
            y1=mindist(DS.Em,peaks(n).Em{1});
            y2=mindist(DS.Em,peaks(n).Em{end});
            Cpeak(i,n)=max(vec(DS.(Xname)(i,y1:y2,x1:x2)),[],'omitnan');
        end
    end
end

%% Indices
FI=nan(1,DS.nSample);
FrI=nan(1,DS.nSample);
BIX=nan(1,DS.nSample);
HIX=nan(1,DS.nSample);
if diagn
    dfig=dreemfig;
    set(dfig,'units','normalized','pos',[0.1344    0.2537    0.7625    0.3194])
    set(dfig,'name','pickPeaks.m - raw vs. smoothed fluorescence for indicies')
end
for i=1:DS.nSample
    try
        % Fluorescence index
        EmScan370=DS.(Xname)(i,:,mindist(DS.Ex,370));
        EmScan370=naninterp(EmScan370,'PCHIP');
        if diagn
            subplot(1,3,1)
            hold on
            h1=plot(DS.Em,EmScan370,'LineStyle','none','Marker','+','Color','k');hold on;
            h2=plot(DS.Em,smooth(EmScan370,21,'sgolay',2),'LineWidth',1.7);
            axis tight
            cylim=get(gca,'YLim');
            ylim([0 cylim(2)])
            cylim=get(gca,'YLim');
            plot([470 470],cylim,'r','LineStyle','--')
            plot([520 520],cylim,'r','LineStyle','--')
            title('em at ex = 370 (Fl.index)')
            legend([h1 h2],{'raw','smoothed'},'location','best')
            box on
        end
        EmScan370 = smooth(EmScan370,21,'sgolay',2);
        
        Val1=EmScan370(mindist(DS.Em,470));
        Val2=EmScan370(mindist(DS.Em,520));
        FI(i)=Val1/Val2;
        
        % Freshness index
        EmScan310=DS.(Xname)(i,:,mindist(DS.Ex,310));
        EmScan310=naninterp(EmScan310,'PCHIP');
        if diagn
            subplot(1,3,2)
            hold on
            plot(DS.Em,EmScan310,'LineStyle','none','Marker','+','Color','k');
            hold on
            plot(DS.Em,smooth(EmScan310,21,'sgolay',2),'LineWidth',1.7)
            axis tight
            cylim=get(gca,'YLim');
            ylim([0 cylim(2)])
            cylim=get(gca,'YLim');
            plot([380 380],cylim,'r','LineStyle','--')
            y = [cylim(1) cylim(2) cylim(2) cylim(1)];
            x = [420 420 435 435];
            patch(x,y,'red','FaceAlpha',0.5)
            title('em at ex = 310 (Freshness index)')
            box on
        end
        EmScan310 = smooth(EmScan310,21,'sgolay',2);
        
        Val1=EmScan310(mindist(DS.Em,380));
        Val2=max(EmScan310(mindist(DS.Em,420):mindist(DS.Em,435)),[],'omitnan');
        FrI(i)=Val1/Val2;
        
        % Humification index
        if ~HIX_excl
            EmScan254=DS.(Xname)(i,:,mindist(DS.Ex,254));
            EmScan254=naninterp(EmScan254,'PCHIP');
            if diagn
                subplot(1,3,3)
                hold on
                plot(DS.Em,EmScan254,'LineStyle','none','Marker','+','Color','k');
                hold on
                plot(DS.Em,smooth(EmScan254,21,'sgolay',2),'LineWidth',1.7)
                axis tight
                cylim=get(gca,'YLim');
                ylim([0 cylim(2)])
                cylim=get(gca,'YLim');
                y = [cylim(1) cylim(2) cylim(2) cylim(1)];
                x = [435 435 480 480];
                patch(x,y,'red','FaceAlpha',0.5)
                
                y = [cylim(1) cylim(2) cylim(2) cylim(1)];
                x = [300 300 345 345];
                patch(x,y,'red','FaceAlpha',0.5)

                title('em at ex = 254 (HIX)')
                box on
            end
            EmScan254 = smooth(EmScan254,21,'sgolay',2);
            Val1=sum(EmScan254(mindist(DS.Em,435):mindist(DS.Em,480)),'omitnan');
            Val2=sum(EmScan254(mindist(DS.Em,300):mindist(DS.Em,345)),'omitnan')+...
                sum(EmScan254(mindist(DS.Em,435):mindist(DS.Em,480)),'omitnan');
            HIX(i)=Val1/Val2;
        end
        if diagn
            disp(['Spectrum ',num2str(i),' of ',num2str(DS.nSample),...
                '. Press any key to continue or Ctrl + C to cancel.'])
            pause
            subplot(1,3,1),cla,subplot(1,3,2),cla,subplot(1,3,3),cla,box on
        end
        % BIX
        Val1=EmScan310(mindist(DS.Em,380));
        Val2=EmScan310(mindist(DS.Em,430));
        BIX(i)=Val1/Val2;
    catch
        disp(['Issue with sample ',num2str(i),'. Returning NaN''s.'])
        FI(i)=nan;
        FrI(i)=nan;
        BIX(i)=nan;
    end
end


%% Results
VarName=[peaks.name]';
VarName{end+1}='FluI';VarName{end+1}='FrI'; VarName{end+1}='BIX'; VarName{end+1}='HIX';

if ~HIX_excl
    picklist=array2table([Cpeak FI' FrI' BIX' HIX'],...
        'VariableNames',VarName);
else
    picklist=array2table([Cpeak FI' FrI' BIX'],...
        'VariableNames',VarName(1:end-1));
end

switch plt
    case true
        fig1=dreemfig;
        set(fig1,'units','normalized','Name','pickPeaks: Extracted intensities of prefefined fluorescence peaks','pos',[0.2594    0.2296    0.4448    0.5130])
        subplot(2,1,1)
        for n=1:numel(peaks)
            plot(Cpeak(:,n),'LineWidth',1.5),hold on
        end
        legend([peaks.name]','location','best');
        title('Fluorescence peaks')
        xlabel('# of sample in dataset')
        axis tight
        
        subplot(2,1,2)
        plot(FI,'LineWidth',1.5),hold on
        plot(FrI,'LineWidth',1.5),hold on
        plot(BIX,'LineWidth',1.5),hold on
        if ~HIX_excl
            plot(HIX,'LineWidth',1.5),hold on
            legend('Fluorescence index','Freshness index' ,'Biological index','Humification index','location','best')
        else
            legend('Fluorescence index','Freshness index','Biological index','location','best')
        end
        title('Fluorescence indicies')
        xlabel('# of sample in dataset')
        hold off
        axis tight
        dreemfig(fig1);
    otherwise
        disp('No plots shown. Please check the function output to inspect the extracted values')
end
end


function [idx,distance] = mindist( vec,value)
[distance,idx]=min(abs(vec-value));
end

function X = naninterp(X,method)
% Interpolate over NaNs
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)),method);
end