function [X2m, X2m_list]=matchsamples(filelist1,filelist2,PairedList,X1,X2)
%
% <strong>Syntax:</strong>
%
%   [X2m, X2m_list]=<strong>matchsamples</strong>(filelist1,filelist2,PairedList,X1,X2)
%
% <a href="matlab: doc matchsamples">help for matchsamples</a> <- click on the link

% Create a 2 or 3 way dataset using data from X2 in a manner which matches
% with matrix X1 according to specified pairings between X1 and X2.
% Note that matches are case-insensitive (i.e. '1ABC.csv' = '1abc.csv')
%
% USEAGE:
% 	[X2m X2m_list]=matchsamples(filelist1,filelist2,PairedList,X1,X2);
%
% INPUTS:
%   filelist1:	list of filenames corresponding to X1
%   filelist2:	list of filenames corresponding to X2
%  PairedList:	2 column matrix with filenames from filelist1 in the first column and
%       		filenames from filelist2 in the second column, indicating pairings between
%       		the two datasets. Filenames in PairedList can appear in any order.
%  		   X1:	dataset corresponding to filelist 1 (e.g. Samples)
% 		   X2:	dataset corresponding to filelist 2 (e.g. Blanks, UV scans, or Absorbance scans)
%               for X1 and X2, samples appear in rows (1st dimension in 3-way case) 
%               and wavelengths are omitted
%
% OUTPUTS:
%  		 X2m:	dataset consisting of samples from X2 matched with filelist1
% 	X2m_list:	list of filenames corresponding to X2m
% 
% EXAMPLES:
%  Assume that
%   X1 and filelist1 are the sample dataset and filenames from ReadInEEMs.m
%   X2 and filelist2 are the blank dataset and filenames, OR
%   X2 and filelist2 are a matrix of Absorbance scans (wavelengths omitted) and filenames, 
%   X2 and filelist2 are a matrix of UV scans (wavelengths omitted) and filenames, 
%   PairedList has filelist1 in column 1, and matched blanks/scans from filelist2 in column 2.
%
% Then, to produce a matrix X2m to correspond with X1, use:
%   [X2m X2m_list]=matchsamples(filelist1,filelist2,PairedList,X1,X2);
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% readinscans: Copyright (C) 2013 Kathleen R. Murphy
% $ Version 0.1.0 $ September 2013 $ First Release
% $ Updated from FDOMcorr toolbox ver. 1.8
%
% Copyright (C) 2013 KR Murphy,
% Water Research Centre
% The University of New South Wales
% Department of Civil and Environmental Engineering
% Sydney 2052, Australia
% krm@unsw.edu.au
%%%%%%%%%%%%%%%%

%% Checks on the size and contents of file lists
%Check for blanks in first row of file lists, and remove if present
if strcmp(cellstr(filelist1(1,:)),{''})==true;
    filelist1=filelist1(2:end,:);
end
if strcmp(cellstr(filelist2(1,:)),{''})==true;
    filelist2=filelist2(2:end,:);
end

%Check that the size of X matches the size of Filenames
% if isequal(size(filelist2,1),size(X2,1))==false;
%     error('Number of samples in filelist 2 does not appear to match the size of X2');
% end
if isequal(size(filelist1,1),size(X1,1))==false;
    error('Number of samples in filelist 1 does not appear to match the size of X1');
end

%Check for duplicate samples in the sample log/paired list
%ContainsDuplicates=false;
nU=size(unique(PairedList(:,1)),1);
nP=size(PairedList(:,1),1);
if nU<nP
    warning('There are duplicate filenames in filelist1')
    %ContainsDuplicates=true;
end

%% Subset PairedList to include only samples in the Filelist1
%% Sort PairedList so that the samples are in the same order as filelist1
iii=NaN*ones(size(PairedList,1),1);
if ~isequal(filelist1,PairedList(:,1));
    if ~isequal(char(filelist1),char(PairedList(:,1)));
        % Subset PairedList
        for i=1:size(PairedList,1)
            CH=char(PairedList(i,1));
            %ii=strncmpi(filelist1,PairedList(i,1),length(CH))';
            %iii(i)=max(ii);
            ii=strncmpi(filelist1,PairedList(i,1),length(CH));
            iii(i)=max(ii);
        end
        PairedList=PairedList(iii==1,:);
        
        % Sort PairedList 
        iii=[];
        try
            for i=1:size(PairedList,1)
                ii=strcmpi(filelist1,PairedList(i,1));
                iii(i)=find(ii==1);
            end
        catch
            disp('filelist1')
            disp(filelist1)
            disp('PairedList')
            disp(PairedList(i,1))
            pause
            error('Sample Name not found. Could be a missing file extension e.g. .csv!')
        end
        [d_1, b]=sortrows(iii'); %#ok<ASGLU>
        PairedList=PairedList(b,:);
    end
end
%size(PairedList),pause

%% Check that all samples in filelist1 are in the paired list.
%Perform case insensitive match of filenames
chk1_case_sensitive=ismember(filelist1,PairedList(:,1));
chk1_case_ignore=ismember(upper(filelist1),upper(PairedList(:,1)));
chk1=max(chk1_case_sensitive,chk1_case_ignore);
%[chk1 chk1_case_sensitive chk1_case_ignore],pause
if min(chk1)==0;
    fprintf(' \n')
    fprintf(' Error! The following filenames appearing in Filelist1  \n')
    fprintf(' were not matched with any filenames in Paired List (1st column) \n')
    disp(filelist1(chk1==0,:))
    fprintf(' Press any key to continue... \n')
    in1=input(' To print the filenames in list 1 and 2 to Excel type "9", otherwise press enter: ');
    if in1==9;
        xlswrite('FileList_Check.xls',cellstr('list 1 from filelist'),'matchpairs1','A1')
        xlswrite('FileList_Check.xls',cellstr('list 2 from filelist'),'matchpairs1','B1')
        xlswrite('FileList_Check.xls',cellstr('matched=1'),'matchpairs1','B1')
        xlswrite('FileList_Check.xls',cellstr('list 1- matched'),'matchpairs1','E1')
        xlswrite('FileList_Check.xls',cellstr('list 2- matched'),'matchpairs1','F1')
        xlswrite('FileList_Check.xls',cellstr(filelist1),'matchpairs1','A2')
        xlswrite('FileList_Check.xls',cellstr(filelist2),'matchpairs1','B2')
        xlswrite('FileList_Check.xls',cellstr(num2str(chk1)),'matchpairs1','C2')
        xlswrite('FileList_Check.xls',cellstr(PairedList),'matchpairs1','E2')
        fprintf('Data have been saved to the Excel file named "FileList_Check.xls"...\n')
    end
    error('some samples could not be matched')
    % fprintf(' Press any key to continue ...\n')
    % pause
end

%% Check that samples in the new paired list are all in filelist2.
%Perform case insensitive match of filenames
chk2_case_sensitive=ismember(PairedList(:,2),filelist2);
chk2_case_ignore=ismember(upper(PairedList(:,2)),upper(filelist2));
chk2=max(chk2_case_sensitive,chk2_case_ignore);
if min(chk2)==0;
    fprintf(' Error!  The following filenames appearing in the Paired List (2nd column)  \n')
    fprintf(' were not matched with any filenames in Filelist2 \n')
    disp(unique(PairedList(chk2==0,2)))
    fprintf(' Press any key to continue... \n')
    in1=input(' To print the filenames in list 1 and 2 to Excel type "9", otherwise press enter: ');
    if in1==9;
        xlswrite('FileList_Check.xls',cellstr('list 1 from filelist'),'matchpairs2','A1')
        xlswrite('FileList_Check.xls',cellstr(filelist1),'matchpairs2','A2')
        xlswrite('FileList_Check.xls',cellstr('matched=1'),'matchpairs2','B1')
        xlswrite('FileList_Check.xls',cellstr(num2str(chk2)),'matchpairs2','B2')
        xlswrite('FileList_Check.xls',cellstr('list 2 from filelist'),'matchpairs2','D1')
        xlswrite('FileList_Check.xls',cellstr(filelist2),'matchpairs2','D2')
        xlswrite('FileList_Check.xls',cellstr('list 1- matched'),'matchpairs2','F1')
        xlswrite('FileList_Check.xls',cellstr('list 2- matched'),'matchpairs2','G1')
        xlswrite('FileList_Check.xls',cellstr(PairedList(chk2==1,:)),'matchpairs2','F2')
        fprintf('Data have been saved to the Excel file named "FileList_Check.xls"...\n')
    end
    error('some samples could not be matched')
end

%% Create a matrix the length of PairedList using matched samples from X2

NoSamples=size(PairedList,1);
SampleIndex=1:size(PairedList,1);

Xdim=max(size(size(X2)));
TFmat=[];
if Xdim==2; %Matrix, 2-way
        if ischar(X2)
            X2m=repmat(' ',[NoSamples 30]);
        elseif isnumeric(X2)
            X2m=NaN*ones(NoSamples,size(X2,2));
        end
    for i=1:size(filelist2,1)
        TF = strcmpi(deblank(filelist2(i,:)),deblank(PairedList(:,2)));
        TFmat=[TFmat TF]; %#ok<*AGROW>
        %sum(TF)
        %X2(i,:)
        %pause
        X2r=repmat(X2(i,:),[sum(TF) 1]);
        if ischar(X2r)
            X2m(SampleIndex(TF==1),1:length(X2r))=X2r;
        else
            X2m(SampleIndex(TF==1),:)=X2r;
        end
    end
elseif Xdim==3 %Cube, 3-way
    X2m=NaN*ones(NoSamples,size(X2,2),size(X2,3));
    for i=1:size(filelist2,1)
        %disp(i)
        TF = strcmpi(deblank(filelist2(i,:)),deblank(PairedList(:,2)));
        TFmat=[TFmat TF];
        X2r=repmat(X2(i,:,:),[sum(TF) 1 1]);
        X2m(SampleIndex(TF==1),:,:)=X2r;
    end
end
%FL1=filelist1;
%FL2=filelist2;
X2m_list=PairedList(:,2);
fprintf('A matched dataset has been created in your workspace. \n');
