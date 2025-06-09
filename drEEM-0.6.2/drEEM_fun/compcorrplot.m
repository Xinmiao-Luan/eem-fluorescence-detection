function compcorrplot(data,f,varargin)
%
% <strong>Syntax</strong>
%   <strong>compcorrplot</strong>(data,f,varargin)
%   <strong>compcorrplot</strong>(data,f,comparisons,codeby,loglog)
%
% <a href="matlab: doc compcorrplot">help for compcorrplot</a> <- click on the link

%Plot PARAFAC components against each other
%
%Useage:   compcorrplot(data,f,comparisons,codeby,loglog)
%Inputs:
%            data: data structure containing model results (e.g. data.Model4).
%              f:  Number of components in model
%    comparisons:  (optional)
%                  identify which components to plot. The default is to
%                  plot all possible combinations of components on one page.
%                  Separate the comparison pairs with a colon:
%                  e.g.
%                    [1 2; 3 4] plots only Comp 1 vs 2 and Comp 3 vs 4.
%        codeby:  (optional)
%                  color-code the data points
%                 'label' : if data.label is a field containing metadata,
%                           a separate plot symbol will be used for
%                           each value of label
%        loglog:  (optional)
%                  use log scale on both axes
%Examples
%   compcorrplot(data,6)
%   compcorrplot(data,6,[1 2; 1 5; 4 6])
%   compcorrplot(data,6,[],'site')
%   compcorrplot(data,6,[1 5],'site','loglog')
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013,
%     DOI:10.1039/c3ay41160e.
%
% compcorrplot: Copyright (C) 2013 Kathleen R. Murphy
% The University of New South Wales
% Dept Civil and Environmental Engineering
% Water Research Center
% UNSW 2052
% Sydney
% krm@unsw.edu.au
%
% $ Version 0.2.0 $ April 2019 $ Showing Rsq values for linear correlations $ UJW
% $ Version 0.1.0 $ September 2013 $ First Release

%Initialise and check inputs
narginchk(2,5)
comparisons=[];
sitedata=[];
LLplot=false;
model=data.(['Model' num2str(f)]);
scores=model{1};
if nargin>2
    comparisons=varargin{1};
    if nargin>3
        codeby=varargin{2};
        if ~isempty(codeby)
            sitedata=data.(codeby);
            sitelist=unique(sitedata);
            if length(sitelist)>50
                error(['Number of codes (' num2str(length(sitelist)) ') exceeds the maximum allowed (50)'])
            end
        end
        if nargin>4
            logplots=varargin{3};
            if strcmp(logplots,'loglog')
                LLplot=true;
            end
        end
    end
end
if isempty(comparisons)
    comparisons=nchoosek(1:f,2);
end
nc=size(comparisons,1);

%plotting
rc=      [1 1;1 2;1 3;2 2;2 3;2 3;2 4;2 4;3 3;2 5;3 4;3 4;4 4;4 4;4 4;4 4;4 5;4 5;4 5;4 5;4 6;4 6;4 6;4 6;5 5;5 6;5 6;5 6;5 6;5 6];
shapes= char('v','o','*','s','d','^','o','p','v','*','s','+','d','^','o','+','d','o','*','s','v','o','*','s','d','v','o','*','s','d','^','o','p','v','*','s','+','d','^','o','+','d','o','*','s','v','o','*','s','d');
colours=char('r','b','b','k','k','r','m','b','c','g','r','k','b','b','r','k','b','c','r','c','k','r','m','b','c','r','b','b','k','k','r','m','b','c','g','r','k','b','b','r','k','b','c','r','c','k','r','m','b','c');
fills=  char('r','g','b','c','y','k','r','g','b','g','y','k','r','g','c','k','w','k','r','b','y','k','r','g','m','r','g','b','c','y','k','r','g','b','g','y','k','r','g','c','k','w','k','r','b','y','k','r','g','m');

hf=dreemfig;
set(hf,'name',['Model ' int2str(f) ' correlations between components in ' num2str(nc) ' plots'])
set(hf,'units', 'normalized', 'pos', [0.1 0.2 0.8 0.4]);
for j=1:nc
    subplot(rc(size(comparisons,1),1),rc(size(comparisons,1),2),j)
    comp1=comparisons(j,1); comp2=comparisons(j,2);
    if isempty(sitedata)
        if LLplot
            h=loglog(scores(:,comp1),scores(:,comp2),...
                'ok','MarkerSize',6);
            addprop(h,'data');
            h.data=data;
            
        else
            h=plot(scores(:,comp1),scores(:,comp2),...
                'ok','MarkerSize',6,...
                'MarkerFaceColor',lines(1),...
                'MarkerEdgeColor','k');
            addprop(h,'data');
            h.data=data;
            
            lm=fitlm(scores(:,comp1),scores(:,comp2));
            pos=get(gca,'Pos');
            annotation(gcf,'textbox',...
                [pos(1) pos(2)-0.01 pos(3) pos(4)],...
                'FitBoxToText','on',...
                'String',{(['R^2 = ',num2str(round(lm.Rsquared.Adjusted,2))])},...
                'FontSize',8,'LineStyle','none');
            axis('tight')
        end
        xlabel(['Comp. ' num2str(comp1)],'fontsize',10);
        ylabel(['Comp. ' num2str(comp2)],'fontsize',10)
    else
        done=false;
        for i=1:length(sitelist)
            imat=strfind(sitedata,char(sitelist(i)));
            nonmatches = cellfun(@isempty,imat);
            matches=logical(ones(size(nonmatches))-nonmatches);
            if ~isempty(imat)
                if LLplot
                    h=loglog(scores(matches,comp1),scores(matches,comp2),...
                        'Marker',shapes(i),...
                        'Linestyle','none',...
                        'MarkerEdgeColor',colours(i),...
                        'MarkerFaceColor',fills(i),...
                        'MarkerSize',6);hold on
                    addprop(h,'data');
                    h.data=subdataset(data,nonmatches,[],[]);
                else
                    h=plot(scores(matches,comp1),scores(matches,comp2),...
                        'Marker',shapes(i),...
                        'Linestyle','none',...
                        'MarkerEdgeColor',colours(i),...
                        'MarkerFaceColor',fills(i),...
                        'MarkerSize',6);hold on
                    addprop(h,'data');
                    h.data=subdataset(data,nonmatches,[],[]);
                    if ~done
                        done=true;
                        lm=fitlm(scores(:,comp1),scores(:,comp2));
                        pos=get(gca,'Pos');
                        annotation(gcf,'textbox',...
                            [pos(1) pos(2)-0.01 pos(3) pos(4)],...
                            'FitBoxToText','on',...
                            'String',{(['R^2 = ',num2str(round(lm.Rsquared.Adjusted,2))])},...
                            'FontSize',8,'LineStyle','none');
                    end
                end
            end
        end
        hold off
        axis('tight')
        xlabel(['Comp. ' num2str(comp1)],'fontsize',10);
        ylabel(['Comp. ' num2str(comp2)],'fontsize',10)
        legend(sitelist)
        if j<nc
            legend off
        end
    end
end
dreemfig(hf);
dcm_obj = datacursormode(hf);
set(dcm_obj,'UpdateFcn',{@cursorinfo});
datacursormode on
end


function txt = cursorinfo(~,event_obj)

idx = get(event_obj, 'DataIndex');

try
    flist=event_obj.Target.data.filelist{idx};
catch
    flist='no filename';
end
try
    txt = {['data.i: ',num2str(event_obj.Target.data.i(idx))],...
        ['filename: ',flist]};
catch
end
end