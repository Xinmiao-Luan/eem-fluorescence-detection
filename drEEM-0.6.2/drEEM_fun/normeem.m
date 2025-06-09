function Xnew=normeem(data,varargin)
%
% <strong>Syntax</strong>
%   Xnew=<strong>normeem</strong>(data,varargin)
%   Xnew=<strong>normeem</strong>(data,options,f,specs,excludelow)
%
% <a href="matlab: doc normeem">help for normeem</a> <- click on the link

% Normalise EEMs to reduce concentration effects, or reverse a normalisation 
% that was applied previously.
%
%USEAGE  Xnew=normeem(data,options,f,specs,excludelow)
%INPUTS
%        data:
%     options: (optional)
%              'reverse'- reverse normalisation if applied previously
%           f: (optional) number of components in model that is to be 
%              unscaled in reverse case.
%       specs: (optional) {convgcrit, constraint}
%              constraints and convergence criteria applied during
%              modelling. If not specified, this will be taken from data
%              in {data.Val_ConvgCrit, data.Val_Constraints}, or else
%              {data.Modelfconvgcrit, data.Modelfconstraints}, or else
%              {data.OutlierTest_convgcrit, data.OutlierTest_constraints}
%              Constraints recognised are:
%              'nonnegativity','unimodnonneg','unconstrained'
%       excludelow:(optional) New in v0.6.0. Automatically exclude or dont exclude low
%       intensity samples. Default: false, meaning the user has to chose.
%
%OUTPUTS
%       Xnew, which is the same as data except that 
%       if options =[], 
%             Xnew.X now contains the normalised dataset;
%             Xnew.Xnotscaled now contains the original dataset;
%       if options ='reverse', Xnew.Modelf contains the f-component model
%          with the true (unscaled) scores.
%             Xnew.X now contains the unscaled dataset;
%             Xnew.Xnorm now contains the normalised dataset;
%             Xnew.Xnotscaled has been removed;
%
%EXAMPLES
%
%       Xnew=normeem(data)   %normalise samples in data.X
%       Xnew=normeem(val6,'reverse',6) %reverse normalisation of 6 comp model
%       Xnew=normeem(data,'reverse',6,{1e-6,'nonnegativity'})
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% normeem: Copyright (C) 2013 Kathleen R. Murphy
% The University of New South Wales
% Dept Civil and Environmental Engineering
% Water Research Center
% UNSW 2052
% Sydney
% krm@unsw.edu.au
%
% $ Version 0.1.0 $ September 2013 $ First Release
% $ Version 0.3.0 $ September 2016 $ Third Release - bugfix constraints

narginchk(1,5)
Xnew=data;
cc=[];const=[];f=[];
excludechoice=false;
excludequery=true;
if nargin==1
    ncase='apply';
end
if nargin>1
    ncase=varargin{1};
    if isempty(ncase)
        ncase='apply';
    end
    if ~contains(ncase,{'reverse','apply'})
        error('When specifying input to ''options'', input has to be ''reverse'' (default is ''apply'').')
    end
    if nargin==2
        error('number of factors must be specified in the reverse case')
    end
    if nargin>2
        f=varargin{2};
        if isempty(f)&&strcmp(ncase,'reverse')
            error('Number of factors in model must be specified')
        end
        if nargin>3
            specs=varargin{3};
            if ~iscell(specs)&&strcmp(ncase,'reverse')
                error('Specify convergence criteria and constraints in curly brackets {}')
            end
            if iscell(specs)&&strcmp(ncase,'reverse')
                cc=specs{1};
                const=specs{2};
            end
        end
        if nargin>4
            if ~islogical(excludechoice)
                error('Input to ''excludelow'' must be either true or false')
            end
            excludechoice=varargin{3};
            excludequery=false;
        end
    end
end
if excludechoice
    excludechoice='y';
else
    excludechoice='n';
end

%%
switch ncase
    case 'apply'
        Xnew.Xnotscaled=data.X;
        [Xnew.X,~,scales]=nprocess(data.X,[0 0 0],[1 0 0]);
        Xnew.Preprocess='Normalised to unit variance in sample mode';
        delidx=find(scales{1}>5*median(scales{1}));
        if any(scales{1}>5*median(scales{1}))
            if excludequery
                disp('--- CAUTION ---')
                disp('Sample(s) with signals far below average signals identified')
                disp('Samples (data.i): ')
                disp(' ')
                try
                    for n=1:numel(delidx)
                        disp(strcat(num2str(data.i(delidx(n))),repmat(" : ",1,1),data.filelist(delidx(n))))
                    end
                catch
                    for n=1:numel(delidx)
                        disp(num2str(data.i(delidx(n))))
                    end
                end
                disp(' ')
                disp('These samples may negatively impact the PARAFAC analysis.')
                disp('Do you wish to subset the dataset to exclude them?')
                disp('  ')
                decision=input('--->  y/n: ','s');
            else
                disp('--- CAUTION ---')
                disp('Sample(s) with signals far below average signals identified')
                disp('Samples (data.i): ')
                disp(' ')
                try
                    for n=1:numel(delidx)
                        disp(strcat(num2str(data.i(delidx(n))),repmat(" : ",1,1),data.filelist(delidx(n))))
                    end
                catch
                    for n=1:numel(delidx)
                        disp(num2str(data.i(delidx(n))))
                    end
                end
                decision=excludechoice;
            end
           switch decision
               case 'y'
                   Xnew = subdataset(Xnew,delidx,[],[]);
                   disp('Samples excluded from dataset.')
               case 'n'
                   disp('Samples remained in dataset.')
               otherwise
                   error('keyboard input either ''y'' or ''n''')
           end
        end
    case 'reverse'
        if ~isfield(data,'Xnotscaled')
            error('data.Xnotscaled does not contain the unscaled X data')
        else
            if ~isequal(size(data.X),size(data.Xnotscaled))
                error('data.X and data.Xnotscaled must be of equal size')
            end
        end
        Db=size(data.Xnotscaled);
        for i=f
            if or(isempty(cc),isempty(const))
                if isfield(data,{'Val_ConvgCrit','Val_Constraints'})
                    cc=data.Val_ConvgCrit;
                    const=data.Val_Constraints;
                elseif isfield(data,{['Model' int2str(i) 'convgcrit'],['Model' int2str(i) 'constraints']})
                    cc=data.(['Model' int2str(i) 'convgcrit']);
                    const=data.(['Model' int2str(i) 'constraints']);
                elseif isfield(data,{'OutlierTest_convgcrit','OutlierTest_constraints'})
                    cc=data.OutlierTest_convgcrit;
                    const=data.OutlierTest_constraints;
                else
                    error('Specify convergence criteria and constraints in ''specs''')
                end
            end
            
            if strcmp(const,'unconstrained') %no constraints
                constr=[0 0 0];
            elseif strcmp(const,'nonnegativity') %all non-negative
                constr=[2 2 2];
            elseif strcmp(const,'unimodnonneg') %unimodal emission
                constr=[2 3 2];
            end

            modeln=data.(['Model' int2str(i)]);
            [~,B,C]=fac2let(modeln);
            forced=nwayparafac(data.Xnotscaled,i,cc,constr,{rand(Db(1),i);B;C},[0 1 1]);
            [Anew,B,C]=fac2let(forced);
            Xnew.(['Model' int2str(i)])={Anew;B;C};
            Xnew.(['Model' int2str(i) 'preprocess'])='Reversed normalisation to recover true scores';
        end
        Xnew.Xnorm=data.X;
        Xnew.X=data.Xnotscaled;
        Xnew=rmfield(Xnew,'Xnotscaled');
end

