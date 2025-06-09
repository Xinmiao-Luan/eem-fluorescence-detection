function C=lookforconflicts(f)
%
% <strong>Syntax</strong>
%   C=<strong>lookforconflicts</strong>(f)
%
% <a href="matlab: doc lookforconflicts">help for lookforconflicts</a> <- click on the link

% Check for potential conflicts between functions having the same names. 
% This function checks whether functions in the folder 
% defined in fn have the same names as functions in other folders 
% on the MATLAB path. If there are conflicting mfiles (i.e. functions that 
% have the same name and appear first on the MATLAB path, warnings and a
% final error are generated.
%
% Input: 
%    folder: name of a folder on the MATLAB path 
%            (enclose the name in inverted commas)
% Output: 
%         C: list of conflicting file names
%
% Examples: 
%         C=lookforconflicts('drEEM');
%         C=lookforconflicts('nway');
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% lookforconflicts: Copyright (C) 2013 Kathleen R. Murphy
% The University of New South Wales
% Dept Civil and Environmental Engineering
% Water Research Center
% UNSW 2052
% Sydney
% krm@unsw.edu.au
%
% $ Version 0.2.0 $ April 2019 $ Second Release
% $ Version 0.1.0 $ September 2013 $ First Release

% New for v0.5.0 (folder renamed)
if strcmp(f,'drEEM')
    f='drEEM_fun';
end

W = what(f);
if size(W,1)==0
    error([ f ' is not a recognised folder on the MATLAB path']);
else
    disp('----')
    disp(' ')
    disp(['checking the following files in ' f ':']);
    disp(W(1).m);
    disp('----')
    disp(' ')
end
s=filesep;
try
    folder=W.path;
catch ME
    disp(ME)
    disp(W)
end

[pathstr,name] = fileparts(folder);
conflicts='';
for i=1:size(W.m,1)
    mf=deblank(char(W.m(i,:)));
    if ~strcmp(mf,'contents.m')
        test=which(mf);
        tPath = fileparts(test);
        if ~strcmp(tPath,[pathstr s name])
            fprintf('\n')
            fprintf('\n')
            conflicts=char(conflicts,mf);
            warning(['An mfile named ' mf ' in ' tPath ' takes precedence over ' mf ' in ' name '. To ensure ' name ' functions work properly, place ' name ' ahead of conflicting folders on the MATLAB path list'])
            
        end
    end
end
if isempty(conflicts)
    C=[];
    fprintf('\n Result: No conflicting functions were found.\n');
else
    C=conflicts(2:end,:);
    error(['\n Result: ' num2str(size(C,1)) ' conflicts were found.\n']);
end

