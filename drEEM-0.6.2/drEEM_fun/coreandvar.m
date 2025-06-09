function [ ccve_tab ] = coreandvar(models)
%
% <strong>Syntax</strong>
%   <strong>coreandvar</strong>(models)
%
% <a href="matlab: doc coreandvar">help for coreandvar</a> <- click on the link

% Inspect core consistency and % explained variance of different models
%
%USEAGE: [ cc,ve ] = coreandvar( models)
%
%INPUT VARIABLES: 
% models: Data structures to be analyzed here.
%
%
%OUTPUT VARIABLES: 
% ccve_tab:  core consistencies and explained variances in a table

%EXAMPLES
% 1.   coreandvar(models)                (just plot cc and ev)
% 2.   ccve_tab = coreandvar( models)    (plots of cc and ev along with table)
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% coreandvar: Copyright (C) 2019 Urban J. Wuensch
% Chalmers University of Technology
% Sven Hultins Gata 6
% 41296 Gothenburg
% Sweden
% wuensch@chalmers.se
% $ Version 0.1.0 $ March 2019 $ First Release

%% Function init
if nargin==0
    help coreandvar
    return
end
narginchk(1,1)

f=zeros(1,100);
for n=1:100
    if isfield(models,['Model',num2str(n)])&&size(models.(['Model',num2str(n)]),2)==3
        f(n)=n;
    else
        f(n)=0;
    end

end
f(f==0)=[];
cc=nan(numel(f),numel(models));
ve=nan(numel(f),numel(models));
comp=nan(numel(f),numel(models));
vpc=nan(numel(f),f(end));
for n=f
    comp(n,1)=n;
    try
        cc(n,1) = models.(['Model',num2str(n),'core']);
    catch
        % When split models are analyzed, the CC was not calculated.
        cc(n,1) = corcond(models.X,models.(['Model',num2str(n)]));
    end
    try
        ve(n,1) = models.(['Model',num2str(n),'percentexpl']);
    catch
        % When split models are analyzed, the VE was not calculated.
        ve(n,1) = 100 * (1 - models.(['Model',num2str(n),'_err']) / sum(models.X(:).^2,'omitnan'));
    end
    try
        vpc(n,1:n)=round(models.(['Model',num2str(n),'compsize']),1);
    catch
        % When split models are analyzed, the compsize was not calculated.
        sizeF=nan(1,n);
        for ii=1:n
            modelledF=nmodel([{models.(['Model',num2str(n)]){1}(:,ii)} {models.(['Model',num2str(n)]){2}(:,ii)} {models.(['Model',num2str(n)]){3}(:,ii)}]);
            sizeF(ii)=100 * (1 - (sum((models.X(:)-modelledF(:)).^2,'omitnan') / sum(models.X(:).^2,'omitnan')));
            
        end
        vpc(n,1:n)=round(sizeF,1);
    end
end

sizename=erase(cellstr(strcat(repmat('Compsize_',f(end),1),num2str((1:f(end))')))',' ');
ccve_tab=array2table([f', cc(f) ve(f),vpc(f,:)],'VariableNames',[{'Factors','Core','ExplVar'},sizename]);
%%
figure1=dreemfig;
set(figure1,'units','normalized','pos',[0.2365    0.3898    0.5099    0.2000])
set(figure1,'Name',char(strcat('Core consistency and %expl. variance, ',{' '},'Model ',{' '},num2str(f))))
subplot(1,3,1)
plot(ccve_tab.Factors,ccve_tab.Core,'Color',lines(1),'LineWidth',1,'Marker','o','MarkerFaceColor',lines(1),'MarkerEdgeColor','k')
line(get(gca,'XLim'),[0 0],'LineStyle','-.','LineWidth',2,'Color','r')
try ylim([min(cc) 100])
catch; ylim([0 100])
end
title('Core Consistency')
ylabel('Core consistency (%)')
xlabel('# of factors')

subplot(1,3,2)
plot(ccve_tab.Factors,ccve_tab.ExplVar,'Color',lines(1),'LineWidth',1,'Marker','o','MarkerFaceColor',lines(1),'MarkerEdgeColor','k')
line(get(gca,'XLim'),[0 0],'LineStyle','-.','LineWidth',2,'Color','r')
try ylim([min(ve(ve~=0)) 100])
catch; ylim([0 100])
end
title('Variance explained')
ylabel('% variance explained')
xlabel('# of factors')


subplot(1,3,3)
col=lines(numel(f));
for n=1:height(ccve_tab)
    plot(1:f(end),table2array(ccve_tab(n,4:end)),'Color',col(n,:),'LineWidth',1,'Marker','o','MarkerFaceColor',col(n,:),'MarkerEdgeColor','k','MarkerSize',5)
    hold on
end
axis tight
try ylim([0 max(max(table2array(ccve_tab(:,4:end))))])
catch; ylim([0 100])
end
title('Component size')
ylabel({'Component size (%)','(does not add to 100%)'})
xlabel('Component no.')
legend(leg(height(ccve_tab),col),cellstr(num2str(ccve_tab.Factors)),'location','best')
dreemfig(figure1);

end

function [ h ] = leg( numPlots,col )
h = gobjects(numPlots, 1);
for n=1:numPlots
    hold on;
    h(n) = line(NaN,NaN,'Marker','o','MarkerSize',5,'MarkerFaceColor',col(n,:),'MarkerEdgeColor','k','LineStyle','none');
end
end