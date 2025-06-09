function comparespectra(data,factors)
%
% <strong>Syntax</strong>
%   <strong>comparespectra</strong>(data,factors)
%
% <a href="matlab: doc comparespectra">help for comparespectra</a> <- click on the link

% Compare Ex and Em loadings for various PARAFAC models. A figure with
% n x 2 plots is produced, where n is the number of models. Ex and Em 
% spectra for each model are shown in overlaid plots. 
%
%USEAGE:   
%		 comparespectra(data,factors)
%INPUTS: 
%  data: dataset structure containing PARAFAC model results in Model3,
%        Model4, ... etc, where Modeln has n components
%     f: Number(s) of components in various models to be compared
%
%OUTPUTS: 
%    A figure is produced with separate plots of Ex and Em spectral 
%    loadings for each component in each model.
%
%EXAMPLE:
%  comparespectra(Test1,3:6)
%
%  See also SPECTRALLOADINGS
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% comparespectra: Copyright (C) 2013 Kathleen R. Murphy
% The University of New South Wales
% Dept Civil and Environmental Engineering
% Water Research Center
% UNSW 2052
% Sydney
% krm@unsw.edu.au
%
% $ Version 0.1.0 $ September 2013 $ First Release

narginchk(2,2)
plotlayout=[length(factors),2];

h=dreemfig;
set(h,'Name',['Compare Ex and Em loadings: ' num2str(factors) ' components']);
for j=1:plotlayout(1)
    modelf=['Model' num2str(factors(j))];
    if ~isfield(data,modelf)
        disp(data)
        disp(['Can not find ' modelf])
        error('CompareComponents:fields',...
            'The dataset does not contain a model with the specified number of factors')
    end
    M = getfield(data,{1,1},modelf);
    B=M{2};C=M{3};
    
    subplot(plotlayout(1),2,2*j-1);
    plot(data.Em,B);
    axis tight
    if j==plotlayout(1)
        xlabel('Em wavelength (nm)');
    end
    if ismember(j,1:length(factors))
        ylabel(['Loads: ' modelf]);
    end
    if j==plotlayout(1)
        legend(num2str((1:size(B,2))'),'location','northeast');
    end    
    subplot(plotlayout(1),2,2*j);
    plot(data.Ex,C);
    axis tight
    if j==plotlayout(1)
        xlabel('Ex wavelength (nm)');
    end
end
dreemfig(h);