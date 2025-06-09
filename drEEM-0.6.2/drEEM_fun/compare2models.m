function compare2models(data,fA,fB,varargin)
%
% <strong>Syntax</strong>
%   <strong>compare2models</strong>(data,fA,fB,varargin)
%   <strong>compare2models</strong>(data,fA,fB,erraxis)
%
% <a href="matlab: doc compare2models">help for compare2models</a> <- click on the link

% Plot measured, modeled and residual EEMs for two different models, one
% set above the other. One sample is shown per page. Press enter to
% progress to the next sample in the dataset. 
%
% USEAGE: compare2models(data,fA,fB,erraxis)
%
% INPUTS
%     data: Source of model e.g Test2
%       fA: Number of components in first model 
%       fB: Number of components in second model 
%  erraxis: (optional) set the relative scale of the error plots
%
% OUTPUTS
%    This function produces on-screen plots (2 rows x 3 columns) with
%    one sample per page. To stop before the last sample press
%    Ctrl-c. For greater control over plot appearance and other options,  
%    use eemview instead.
%
% EXAMPLES
%  compare2models(Test2,4,5) %compare 4 and 5 component models in Test2
%  compare2models(Test2,4,5,0.05) %set scale in residual plots
%                                  to ±5 percent of max. sample intensity
%
%  See also EEMVIEW
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% compare2models: Copyright (C) 2013 Kathleen R. Murphy
% The University of New South Wales
% Dept Civil and Environmental Engineering
% Water Research Center
% UNSW 2052
% Sydney
% krm@unsw.edu.au
%
% $ Version 0.1.0 $ September 2013 $ First Release
%%
narginchk(3,4)
erraxis=[];
if nargin>3
    erraxis=varargin{1};
end
datasource={data,data,data,data};
datafield={'X',['Model',int2str(fA)],'error_residuals','X',['Model',int2str(fB)],'error_residuals'};
istart=[];
plotview=[];
plotlayout=[2 3];
datalabel=[repmat(char('data ','model ','error '),[2 1]) [repmat(int2str(fA),[3,1]);repmat(int2str(fB),[3,1])]];
datalabel=cellstr(datalabel)';
axeslimits=[];
if ~isempty(erraxis)
    axeslimits=[1 1 erraxis 1 1 erraxis];
end
dorotate=[];
docolorbar=[];

eemview(datasource,datafield,plotlayout,istart,plotview,datalabel,axeslimits,dorotate,docolorbar)

