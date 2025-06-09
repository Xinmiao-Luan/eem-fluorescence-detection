function DS = assembledataset(X, Ex, Em, unit, varargin)
%
% <strong>Syntax:</strong>
%
%   DS = <strong>assembledataset</strong>(X, Ex, Em, unit, varargin)
%
% <a href="matlab: doc assembledataset">help for assembledataset</a> <- click on the link

%Assemble EEMs and metadata into a single dataset structure for processing
%with the drEEM toolbox.
%
% USEAGE
% DS = assembledataset(X,Ex,Em, units,'data1name',data1,'data2name',data2,...., dataN, NUM)
%
%INPUTS
%  X = 3D dataset of EEMs
%  Ex = Ex wavelengths (a vector)
%  Em = Em wavelengths (a vector)
%  unit = Unit of measure for fluorescence intensity e.g. 'RU','QSE' or 'AU'
% (optional) sample metadata corresponding to the rows of X
%            Each field name enclosed in inverted commas is followed by the
%            name of the variable containing the field data.
%              e.g.'Site',site,'TempC',t,'pH',pHdata,...
%              produces data.Site containing the variable site
%                       data.TempC containing the variable t
%                       data.pH containing the variable pHdata
%
% (compulsary only if some optional metadata is formatted as number vectors) 
%            indicate which metadata are numeric rather than categorical
%             (e.g. pH,temp) 
%            NUM =[] => all  metadata are categorical
%            NUM =[2,3] indicates 2nd and 3rd metadata fields contain numeric data
%                 in the above example this corresponds to temperature and pH
%
%OUTPUTS
%  DS = a data structure containing the assembled datasets
%
%EXAMPLES
%  DS = assembledatset(X,Ex,Em,'QSE','site',sites,'date',dates,'ID',sampleID,[])
%  DS = assembledatset(X,Ex,Em,'RU','longID',filelist_eem,'pH',pH,[2])
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 
%     5, 6557-6566, 2013. DOI:10.1039/c3ay41160e. 
%
% assembledatset: Copyright (C) 2013 Kathleen R. Murphy
% Water Environment Technology
% Chalmers University of Technology
% Department of Civil and Environmental Engineering
% Water Environment Technology
% Sven Hultins gata 8
% 412 96 Gteborg
% murphyk@chalmers.se
%
% $ Version 0.1.0 $ September 2013 $ First Release


N_samples=size(X,1);
N_ex=max(size(Ex));
N_em=max(size(Em));
Em=CheckMatrixDim(Em,N_em,1,[],'Em');
Ex=CheckMatrixDim(Ex,N_ex,1,[],'Ex');
X=CheckMatrixDim(X,N_samples,N_em,N_ex,'X');
DS.Ex = Ex;    
DS.Em = Em;    
DS.X=X;
DS.IntensityUnit=unit;
DS.nEx=N_ex;     
DS.nEm=N_em;      
DS.nSample=N_samples;
DS.i=(1:N_samples)';
if ~isempty(varargin)
    metanames=varargin(1:2:end-1);
    metadata=varargin(2:2:end-1);
    metanum=varargin(end);
    metanum=metanum{:};
    if ~isequal(size(metadata,2),size(metanames,2))
        warning('Follow each metadata field name by the variable containing the metadata...')
        warning('...And follow that with a list of numeric fields, e.g.[1,4] or [] ')
        error('assembledataset:CheckMetaDat','Incomplete metadata');
    end
    for i=1:size(metanames,2)
        F=char(metanames(i));
        CheckMatrixDim(metadata{i},N_samples,1,[],[':''' F '''']);
        if ismember(i,metanum)
            if isnumeric(metadata{i})
                DS.(F)=metadata{i};
            else
                disp(F)
                warning(['Data contained in ''' F ''' are non-numeric'])
                error('assembledataset:nonnumeric',...
                    'Attempted to classify non-numeric metadata as numeric')
            end
        else
            if isnumeric(metadata{i})
                DS.(F)=cellstr(num2str(metadata{i}));
            elseif ischar(metadata{i})
                DS.(F)=cellstr(metadata{i});
            elseif iscell(metadata{i})
                DS.(F)=metadata{i};
            end
        end
    end
end
end

function M=CheckMatrixDim(M,Dim1,Dim2,Dim3,Mname)
% M=CheckMatrixDim(M,Dim1,Dim2,Dim3,Mname)
% Check the dimensions of matrix M (name = 'Mname') against expected dimensions Dim1, Dim2,and Dim3
% Transpose M if necessary to achieve correct dimensions or else generate error message.
% Copyright 2010 K.R. Murphy (FDOMcorr toolbox)

if isempty(Dim3); %2D data
    if size(M,1)~=Dim1;
        if and(size(M,2)==Dim1,size(M,1)==Dim2)
            M=M';
        elseif isempty(Dim2) %Number of columns is arbitrary
            if size(M,2)==Dim1
                M=M';
            else
                fprintf(['Check size of matrix ' Mname '.\n'])
                fprintf(['Expecting ' num2str(Dim1) ' rows.\n']),pause
                fprintf(['Current size is ' num2str(size(M,1)) ' rows and ' num2str(size(M,2)) ' columns.\n']),pause
                fprintf('Hit any key to continue.\n'),pause
                error('Unexpected Matrix Size')
            end
        else
            fprintf(['Check size of matrix ' Mname '.\n'])
            fprintf(['Expecting ' num2str(Dim1) ' rows and ' num2str(Dim2) ' columns.\n']),pause
            fprintf(['Current size is ' num2str(size(M,1)) ' rows and ' num2str(size(M,2)) ' columns.\n']),pause
            fprintf('Hit any key to continue.\n'),pause
            error('Unexpected Matrix Size')
        end
    end
else  %3D data
    if isequal(size(M),[Dim1 Dim2 Dim3]);
    else
        fprintf(['Check size of matrix ' Mname '.\n'])
        fprintf(['Expecting ' num2str(Dim1) ' x ' num2str(Dim2) ' x ' num2str(Dim3) '.\n']),pause
        fprintf(['Current size is ' num2str(size(M,1)) ' x ' num2str(size(M,2)) ' x ' num2str(size(M,3)) '.\n']),pause
        fprintf('Hit any key to continue.\n'),pause
        error('Unexpected Matrix Size')
    end
end
end
