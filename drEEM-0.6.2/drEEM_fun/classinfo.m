function classinfo(data)
%
% <strong>Syntax:</strong>
%
%   <strong>classinfo</strong>(data)
%
% <a href="matlab: doc classinfo">help for classinfo</a> <- click on the link

% Produce a summary of the metadata contained in an EEM dataset structure
%
% INPUT
%     data: A data structure containing Ex, Em, X and associated metadata
%           and/or models
% OUTPUT
%     A summary of metadata (information describing samples) is printed 
%     to screen, including the identity and number of unique values of each
%     metadata field.
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% classinfo: Copyright (C) 2013 Kathleen R. Murphy
% The University of New South Wales
% Dept Civil and Environmental Engineering
% Water Research Center
% UNSW 2052
% Sydney
% krm@unsw.edu.au
%
% $ Version 0.1.0 $ September 2013 $ First Release


%% Validate inputs
F=fieldnames(data);
%validation of EEM dataset structure
if sum(strcmp(F,'X'))==1,else error('Data structure is missing compulsary field: data.X'); end
if sum(strcmp(F,'Ex'))==1,else error('Data structure is missing compulsary field: data.Ex'); end
if sum(strcmp(F,'Em'))==1,else error('Data structure is missing compulsary field: data.Em'); end
if sum(strcmp(F,'nSample'))==1,else data.nSample=size(X,1); end
if sum(strcmp(F,'nEx'))==1,else data.nEx=length(data.Ex); end
if sum(strcmp(F,'nEm'))==1,else data.nEm=length(data.Em); end
    
%Class information
for i=1:size(F,1)
    fldnm=char(F(i));
    switch fldnm
        case {'X','Ex','Em','nSample','nEx','nEm','Xf','backupX','backupEx','backupEm','backupXf','i','Smooth','Zap','IntensityUnit'}
        otherwise
            disp(blanks(2)')
            f_i = data.(fldnm);
            fprintf(['Field: ' fldnm '\n'])
            try
                u_i=unique(f_i,'stable')';
                fprintf([fldnm ' # of rows = ' num2str(size(f_i,1)) '\n'])
                fprintf([fldnm ' # of unique values = ' num2str(length(u_i)) '\n'])
                fprintf('displaying first 20 values in order of appearance:\n')
                if length(u_i)>20
                    disp(u_i(:,1:20))
                else
                    disp(u_i)
                end
            catch ME
                if isstruct(f_i)
                    [l1,w1]=size(f_i);
                    FN2=fieldnames(data.(fldnm));
                    d2=size(FN2,1);
                    fprintf([fldnm ' is a ' int2str(l1) ' x ' int2str(w1) ...
                        ' data array with ' num2str(d2) ' fields as listed ' ...
                        'below\n'])
                    disp(FN2')
                elseif isequal(size(f_i),[1,3])
                    fprintf('This field may contain PARAFAC model loadings')
                else
                    disp(ME)
                end
            end
            fprintf('\n')
            disp(' %%%%%%%%%%% ')
            fprintf('\n')
    end
end
