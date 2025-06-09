function [] = errorsandleverages(data,f)
%
% <strong>Syntax</strong>
%   <strong>errorsandleverages</strong>(data,f)
%
% <a href="matlab: doc errorsandleverages">help for errorsandleverages</a> <- click on the link

% Plot sum-of-squared-errors against leverages in all three dataset modes.
% A variable / sample with high leverage and low SSE has a high influence
% on the model. Low leverage / high residual variables or samples have
% little influence. This helps to identify problematic samples / variables
%
% errorandleverage(data,f)
%
%INPUT VARIABLES: 
% data:        One data structures to be plotted.
% f:           factor, number of components
%
%EXAMPLES
% 1.  errorsandleverages(LSmodel,6)
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% errorsandleverages: Copyright (C) 2019 Urban J. Wuensch
% Chalmers University of Technology
% Sven Hultins Gata 6
% 41296 Gothenburg
% Sweden
% $ Version 0.1.0 $ March 2019 $ First Release

%%
[A,B,C]=fac2let(data.(['Model',num2str(f)]));
levA=diag(A*(A'*A)^-1*A');
levB=diag(B*(B'*B)^-1*B');
levC=diag(C*(C'*C)^-1*C');

dat=data.X;
mod=nmodel(data.(['Model',num2str(f)]));
res=dat-mod;

resA=squeeze(nansum(nansum(res.^2,2),3));
resB=squeeze(nansum(nansum(res.^2,1),3))';
resC=squeeze(nansum(nansum(res.^2,1),2));

%%
hf=dreemfig;
set(hf,'Name',['Model',num2str(f),': Sum-of-squared-error (SSE) vs. leverage']);
set(hf,'units','normalized','pos',[0.2938    0.1120    0.4047    0.7009])
subplot(3,1,1)
scatter(levA,resA,30,data.i,'filled' ,'MarkerEdgeColor',[.5 .5 .5],'DisplayName','sample');
cmap=parula(data.nSample);
colormap(gca,cmap)
c=colorbar;
ylabel(c,'Sample index [field ''i'']')
hold on
ylabel('SSE'),xlabel('Leverage')
title('Sample mode')
ylim([0 inf])
xlim([0 inf])
box on

subplot(3,1,2)
scatter(levB,resB,30,data.Em,'filled' ,'MarkerEdgeColor',[.5 .5 .5],'DisplayName','emission');
cmap=parula(data.nEm);
colormap(gca,cmap)
c=colorbar;
ylabel(c,'Emission wavelength')
hold on
ylabel('SSE'),xlabel('Leverage')
title('Emission mode')
ylim([0 inf])
xlim([0 inf])
box on

subplot(3,1,3)
scatter(levC,resC,30,data.Ex,'filled' ,'MarkerEdgeColor',[.5 .5 .5],'DisplayName','excitation');
cmap=parula(data.nEx);
colormap(gca,cmap)
c=colorbar;
ylabel(c,'Excitation wavelength')
hold on
ylabel('SSE'),xlabel('Leverage')
title('Excitation mode')
ylim([0 inf])
xlim([0 inf])
box on

dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@parse4datacursor,data});
datacursormode on

dreemfig(hf);

end