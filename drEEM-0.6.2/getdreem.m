function getdreem(developer) %#ok<*INUSD>
%
% <strong>Syntax</strong>
%   <strong>getdreem</strong>
%
% <a href="matlab: doc getdreem">help for getdreem</a> <- click on the link

% Download and install the drEEM toolbox. This function:
%  - fetches the version-list from 'https://gitlab.com/dreem/drEEM'
%  - downloads the latest stable release
%  - extracts & installs the toolbox with the routine provided 
%     in the respective version
%
% Notes: 
%   -  Developer option: If any input to ''developer'' is provided (any number or 
%      character), the user will be ased to specify the desired version.
%       - 'nightly beta' contains the latest (not yet released) features of drEEM.
%          This version can be considered generally stable, but the code may
%          change without anouncement. 
%          View 'CHANGELOG.md' for the list of ongoing changes.
%          It is not recommended for download.
%       - 'latest released' corresponds to the first entry. This is the
%          version intended for the normal user.
%
%
% Notice:
% This mfile is be part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013,
%     DOI:10.1039/c3ay41160e.
%
% getdreem: : Copyright (C) 2019 Urban J. Wuensch
% Chalmers University of Technology
% Sven Hultins Gata 6
% 41296 Gothenburg
% Sweden
% wuensch@chalmers.se
%
% $ Version 0.1.0 $ August 2019 $ First Release
%% Check input

% if ~ischar(reqver)
%     error('Input to ''reqver'' must be a character, e.g. ''0.5.0''')
% end
% 
olddir=pwd;

if ~checkonlinestatus
    error('Cannot reach ''google.com''. Your computer must have internet access for this function to work.')
end
disp(' ')
disp(' ')
disp(' ')
disp('----------------------------------------------------------------------------')
disp('   drEEM toolbox - decomposition routines for Emission-Excitation Matrices  ')
disp(' ')
disp('               Download and installation routine                            ')
disp(' ')
disp('----------------------------------------------------------------------------')
disp(' ')
disp('Fetching version directory ...')
try
    vhist=webread('https://gitlab.com/dreem/drEEM/raw/master/versions.txt');
catch
    try
        o = weboptions('CertificateFilename','');
        vhist=webread('https://gitlab.com/dreem/drEEM/raw/master/versions.txt',o);
    catch
        iaccess=checkonlinestatus;
        if ~iaccess
            error('Cannot reach ''google.com''. Your computer must have internet access for this function to work as intended.')
        elseif iaccess
            error('Cannot reach the version file, but internet connection is available. Please contact drEEM user support.')
        end
    end
end
disp('Success. ')
cellvhist=textscan(vhist,'%s%s','Delimiter',';');
tvhist=cell2table([cellvhist{:}],'VariableNames',{'Version','URL'});

if exist('developer','var')
    disp(' ')
    disp('--- Available versions -----------')
    message=strcat(num2str((1:height(tvhist))'),repmat(" : ",height(tvhist),1),tvhist.Version);
    for n=1:size(message,1)
        disp(message(n,:))
    end
    disp('----------------------------------')
    disp(' ')
    ui=input(['Enter number of the desired release (',num2str(1:2),' ... ):  '],'s');
    reqver=tvhist.Version{str2double(ui)}; 
else
    reqver=tvhist.Version{strcmp(tvhist.Version,'latest released (same as first entry)')};
end

urlselected=find(strcmp(tvhist.Version,reqver));


if isempty(urlselected)
    error('Something went wrong when trying to download drEEM \n%s',...
        '   Visit: https://gitlab.com/dreem/drEEM or: dreem.openfluor.org')
else
    temp=strsplit(tvhist.URL{urlselected},'/');temp=temp{end};temp=strsplit(temp,'.z');temp=temp{1};
    disp(['Downloading toolbox (',userpath,filesep,temp,filesep,')'])   
    clearvars temp
    try
        unzip(tvhist.URL{urlselected},userpath)
    catch
        try
            o = weboptions('CertificateFilename','');
            websave('dreemtemp',tvhist.URL{urlselected},o);
            unzip('dreemtemp.zip',userpath);
            delete('dreemtemp.zip')
        catch
            error('Could not download the requested version of drEEM. Contact support if the issue persists.')
        end
    end
    disp('Success. ')
    folder=strsplit(tvhist.URL{urlselected},{'/'});
    folder=folder{end};
    folder=strsplit(folder,'.zip');
    folder=folder{1};
    cd([userpath,filesep,folder])
    disp('Installing toolbox ... ')
    disp(' ')
    disp(' ')

    dreeminstall
    cd(olddir)
end

end

function access = checkonlinestatus
if ispc
    C = evalc('!ping -n 1 8.8.8.8');    
elseif isunix
    C = evalc('!ping -c 1 8.8.8.8');        
end

if ~contains(C,'unreachable')
    access=true;
else
    access=false;
end

end