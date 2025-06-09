function dreeminstall
%
% <strong>Syntax</strong>
%   <strong>dreeminstall</strong>
%
% <a href="matlab: doc dreeminstall">help for dreeminstall</a> <- click on the link

% Install the drEEM toolbox
% This function makes sure that the drEEM toolbox is added to the MATLAB
% path. It finds all current folders containing the string "drEEM", removes
% them, and adds the current version of the toolbox.
%
% The following checks are also performed:
%    1) Are old versions of drEEM in the MATLAB path?
%    2) Do conflicting functions take priority over drEEM-functions?
%    3) Is MATLAB newer than 9.0 and are necessary products installed?
%    4) Can parallel computing be used?
%
% IMPORTANT: Change the path to the path of the drEEM toolbox you want
% installed. Otherwise, this function will throw an error.
%
% EXAMPLE:
%         cd 'D:\Work\drEEM-v060';
%         dreeminstall
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013,
%     DOI:10.1039/c3ay41160e.
%
% dreeminstall: : Copyright (C) 2019 Urban J. Wuensch
% Chalmers University of Technology
% Sven Hultins Gata 6
% 41296 Gothenburg
% Sweden
% wuensch@chalmers.se
%
% $ Version 0.2.0 $ August 2020 $ Second Release for 0.6.0
% $ Version 0.1.0 $ April 2019 $ First Release
%% Initial check before anything else
% Check if the highest priority dreeminstall is in the current folder & if
% this specific function being executed is also the one in the current
% path.
if all([~strcmp(fileparts(which('dreeminstall.m')),pwd),~strcmp(fileparts(mfilename('fullpath')),pwd)])
    error('You must change the MATLAB directory to the one containing ''dreeminstall'' in order to correctly install drEEM')
end

% commithist=webread('https://gitlab.com/dreem/drEEM/commits/master?feed_token=_T6Ybe8-Fscb4Drzxqz-&format=atom');
% update=strfind(commithist,'<updated>');
% versdate=commithist(update(1)+9:update(1)+18);
% THIS NEXT LINE NEEDS CHANGING FOR EACH RELEASE!!!
vers='0.6.2';
versdate='October 2020';
clearvars commithist update

%% Welcome messages
disp(' ')
disp(' ')
disp('-----------------------------------------------------')
disp(['drEEM toolbox (',vers,', last update ',versdate,')'])
disp('-----------------------------------------------------')
disp(' ')
disp('Installation routine')
disp(['Installing drEEM stored in: ',pwd])
disp(' ')
pause(2)

%% STEP 1
disp('Step 1/6: Checking your MATLAB version and installed MATLAB products');pause(0.2)

mver=ver('MATLAB');
if str2double(mver.Version)<9.0
    warning(['Consider updating MATLAB in order for all functions to work as intended. drEEM v', vers ,' was only tested on R2019a. Current version: ',mver.Release(2:end-1)])
end
mallver=ver;
tbs={'Statistics and Machine Learning Toolbox' 'Parallel Computing Toolbox' 'Curve Fitting Toolbox'};
isthere=zeros(1,numel(tbs));
for n=1:numel(tbs)
    isthere(n)=any(~cellfun(@isempty,strfind({mallver.Name},tbs{n})));
end
if all(isthere)
    disp('          Success.')
else
    warning(['Missing toolbox(es):',tbs{~isthere},'. Some functions may not work as intended.'])
    disp('Press any key to acknowledge.')
    pause
end


%% STEP 2
disp('Step 2/6: Checking for old drEEM versions.');pause(0.2)
pbackup=path;


% Method:
% 1) Obtain matlabpath
% 2) Obtain list of functions contained in current drEEM toolbox
% 3) Cycle through all matlabpath entries, check for existing functions
%    that are part of the toolbox. If number of hits is > 50% of number of functions
%    flag that folder for removal.
% 4) Any flagged folder is removed.
CurPath = path;CurPath=textscan(CurPath,'%s','delimiter',pathsep);CurPath=CurPath{:};
queries={'drEEM_fun','nway','dreemroot','tutdemo'};
isold=false(numel(CurPath),1);
for n=1:numel(queries)
    funct=returnmfiles(queries{n});
    for i=1:numel(CurPath)
        list=dir(CurPath{i});
        list={list(:).name};
        sumhit=0;
        for k=1:numel(funct)
            if ~isempty(find(strcmp(list,funct{k})))
                sumhit=sumhit+1;
            end
        end
        if sumhit/numel(funct)>=.5
        isold(i)=true;
        end
    end
end

if any(isold)
    disp('          Other version of drEEM is currently installed.')
    disp('          These folders will now be removed from the MATLAB path.')
    disp('          They are NOT deleted, just not searchable in MATLAB anymore.')
    disp('          These folders are:')
    disp('  ')
    i=find(isold);
    for n=1:sum(isold)
        disp(['            ',CurPath{i(n)}])
    end
    disp('  ')
    pause(2)
    rmpath(CurPath{isold})  % This removes old drEEM folders
    disp('          Old drEEM folders removed from MATLAB path.')
else
    disp('          No previous version of drEEM found.')
end
%% STEP 3
disp('Step 3/6: Adding the new toolbox (make sure you changed the current directory)');pause(0.2)
fd=which('dreeminstall');
p=genpath(fd(1:end-14));
warning off
addpath(p)
warning on
disp('          Success.')
%% STEP 4

disp('Step 4/6: Checking for conflicting function names for drEEM...');pause(0.2)
try
    lookforconflictsinstall('drEEM');
    disp('          Success.')
catch ME
    matlabpath(pbackup)
    rethrow(ME)
end
%% STEP 5
disp('Step 5/6: Checking for conflicting function names for nway...');pause(0.2)
try
    lookforconflictsinstall('nway');
    disp('          Success.')
catch ME
    matlabpath(pbackup)
    rethrow(ME)
end



%% STEP 6
disp('Step 6/6: Initializing a few things (this could take a minute).');pause(0.2)
try
    gcp;
catch
end
% The first plots can be a bit slow, this is just calling and closing one.
fig1=figure;
plot(nan,nan)
close(fig1)
%% END
disp('-----------------------------------------------------')
disp(' ')
disp('Installation complete.')


end


function lookforconflictsinstall(f)
% (C) Urban Wuensch (2019)
% Quick rename of drEEM to drEEM_fun to reflect the toolbox reorganization
if strcmp(f,'drEEM')
    f='drEEM_fun';
end
%% Step 1: Check for folders.
% Check if folder is found, if multiple folders with the same name are
% present & if so: Is the installation prioritized?
W = what(f);
if size(W,1)==0
    error([ f ' is not a recognised folder on the MATLAB path']);
elseif size(W,1)>1
    if strcmp(W(1).path,[pwd,filesep,f])
        disp(['          ',f,': Multiple folders with the same name, but the drEEM toolbox is prioritized.'])
    else
        error([f,': Multiple folders with the same name & the installed one is not prioritized.'])
    end
end
%% Step 1: Check for functions.
% Check that each of the functions of drEEM are prioritized.
mfilesstock=returnmfiles(f);
for n=1:numel(mfilesstock)
    a=which(mfilesstock{n},'-all');
    if isempty(a)
        error(['Function ''',mfilesstock{n},''' should be there, but was not found.'])
    elseif iscell(a)&&numel(a)>1
        if strcmp(fileparts(a{1}),[pwd,filesep,f])
            %disp(['          ',mfilesstock{n},': Multiple functions with the same name, but the drEEM toolbox is prioritized.'])
            % Line is commented, because this message is kind of irrelevant
        else
            error([mfilesstock{n},': Multiple functions with the same name & the installed one is not prioritized.'])
        end
    end
end
end


function mlist=returnmfiles(tbx)

switch tbx
    case 'drEEM_fun'
        mlist={'assembledataset.m';'checkdataset.m';'classinfo.m';...
            'compare2models.m';'comparespectra.m';'compcorrplot.m';...
            'coreandvar.m';'describecomp.m';'diffeem.m';'dreemfig.m';...
            'eemreview.m';'eemview.m';'errorsandleverages.m';...
            'fdomcorrect.m';'fingerprint.m';'loadingsandleverages.m';...
            'lookforconflicts.m';'handlescatter.m';'matchsamples.m';'metadata.m';...
            'modelout.m';'modeloutmac.m';'normeem.m';'nwayparafac.m';...
            'openfluor.m';'openfluormatches.m';'outliertest.m';'pickpeaks.m';...
            'ramanintegrationrange.m';'randinitanal.m';'readineems.m';...
            'readinscans.m';'readlogfile.m';'relcomporder.m';...
            'scanview.m';'scores2fmax.m';'slopefit.m';'smootheem.m';...
            'specsse.m';'spectralloadings.m';'spectralvariance.m';...
            'splitanalysis.m';'splitds.m';'splitvalidation.m';...
            'subdataset.m';'undilute.m';'zap.m'};
    case 'nway'
        mlist={'calcore.m';'ckron.m';'cmatrep.m';'complpol.m';'contents.m';...
            'corcond.m';'coredian.m';'coreswdn.m';'corevarn.m';'demos.m';'derdia3.m';...
            'derswd3.m';'dervar3.m';'dtld.m';'eemtimize.m';'explcore.m';'fac2let.m';...
            'fastnnls.m';'fnipals.m';'fnnls.m';'getindxn.m';'gram.m';'gsm.m';'ini.m';...
            'inituck.m';'kr.m';'krb.m';'maxdia3.m';'maxswd3.m';'maxvar3.m';'missmean.m';...
            'missmult.m';'misssum.m';'monreg.m';'ncosine.m';'ncrossdecomp.m';'ncrossdecompn.m';...
            'ncrossreg.m';'neye.m';'nident.m';'nmodel.m';'nonneg.m';'normit.m';'npls.m';'npred.m';...
            'nprocess.m';'nsetdiff.m';'nshape.m';'ntimes.m';'parademo.m';'pfls.m';...
            'pfplot.m';'pftest.m';'plotfac.m';'ppp.m';'refold3.m';'setnans1.m';'setopts.m';...
            'stdnan.m';'t3core.m';'tuckdemo.m';'tucker.m';'tucker2.m';'tucktest.m';'two2n.m';...
            'ulsr.m';'unimodal.m';'unimodalcrossproducts.m'};
    case 'dreemroot'
        mlist={'dreeminstall.m';'getdreem.m';'Contents.m'};
    case 'tutdemo'
        mlist={'drEEM_dataImport.m';'drEEM_dataImport_AL.m';'drEEM_parafac_tutorial_portSurvey.m';'drEEM_parafac_tutorial_wuensch.m';'drEEM_scattertreatment_tutorial_wuensch.m';'readme.md'};
    otherwise
        error('Strange... This folder / toolbox is unknown.')
end
end