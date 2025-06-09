function [A,c,v] = alignds(SampleLog,refinfo,newinfo,varargin)
%
% <strong>Syntax:</strong>
%
%   [A,c,v] = <strong>alignds</strong>(SampleLog,refinfo,newinfo,casesen)
%
% <a href="matlab: doc alignds">help for alignds</a> <- click on the link

% Align two datasets using a SampleLog
%
% [A,c,v] = alignds(SampleLog,refinfo,newinfo,casesen)
%
%INPUTS
%SampleLog : a data structure created with readlogfile.m 
%
%refinfo : a cell matrix, with {reflog,reffilelist}
%            reflog: the name of the field in the SampleLog that contains
%                    the file names for the reference dataset (this
%                    corresponds to a column header in the CSV log file)
%            reffilelist: the filenames corresponding to the reference
%                    dataset.
%
%newinfo : a cell matrix, with either of (a) or (b)
%           (a) {newlog,newfilelist,newdata}; or
%           (b) {newlog}
%           where:
%             newlog : the name of the field in the SampleLog that contains
%                      the file names for the dataset that will be aligned,
%                      which corresponds to a column header in the CSV log.
%             newfilelist: the list of files corresponding to newdata
%             newdata: 2 or 3-way dataset to be aligned e.g. EEMs or scans
%
%          choose (a) to match a loaded 2D or 3D dataset, 
%                e.g. a dataset of Abs scans or EEMs
%          choose (b) to extract numbers or text directly from the log, 
%                e.g. a list of SampleID or dilution factors.
%
%casesen (optional): default is a case-sensitive match between the log and the list
%         of filenames in matlab. To make it case-insenstive so that 
%         e.g. 'Samp1b' = 'samp1b' set casesen=0, otherwise casesen=1;
%
%OUTPUTS
%A : aligned data
%c : filelist for aligned data.
%v : indices of aligned data, so that c=newfilelist(v).
%
% examples case (a)
% Sabs= alignds(SampleLog,{'Log_EEMfile',filelist_eem},{'ABSfile',filelist_abs,S_abs});
% Sr= alignds(SampleLog,{'Log_EEMfile',filelist_eem},{'RamanFile',filelist_R,S_R});
% B= alignds(SampleLog,{'Log_EEMfile',filelist_eem},{'BlankFile',filelist_b,X_b});
% note: In the first of the above three examples, SampleLog.Log_EEMfile
%       contains the list of EEM files derived from the CSV Sample Log, and 
%       SampleLog.ABSfile contains the list of absorbance data files.
%
% examples case (b)
% SampleID= alignds(SampleLog,{'Log_EEMfile',filelist_eem},{'Log_ID'});
% date= alignds(SampleLog,{'Log_EEMfile',filelist_eem},{'Log_Date'});
% cruise= alignds(SampleLog,{'Log_EEMfile',filelist_eem},{'Log_Cruise'});
% note: In the first of the above three examples, SampleLog.Log_EEMfile
%       contains the list of EEM files derived from the CSV Sample Log, and 
%       SampleLog.Log_ID contains the list of ID numbers.      
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% alignds: Copyright (C) 2014 Kathleen R. Murphy
% $ Version 0.2.0 $ June 2014 $ Second Release
%
% Water Environment Technology
% Chalmers University of Technology
% Department of Civil and Environmental Engineering
% Water Environment Technology
% Sven Hultins gata 8
% 412 96 Göteborg
% murphyk@chalmers.se

narginchk(3,4)

casesen=true;
if nargin>3
    casesen=varargin{1};
end
if ~casesen
    disp('Performing case-insensitive search for matching file names')
elseif casesen
    disp('Performing case-sensitive search for matching file names')
end

if max(size(refinfo))==2
    refinlog=SampleLog.(refinfo{1});
    reffilelist=cellstr(refinfo{2});
else
    error('invalid input: refinfo')
end

if max(size(newinfo))==1
    newinlog=SampleLog.(newinfo{1});
    if ~isnumeric(newinlog)
        c=cell(size(reffilelist,1),1);
    else
        c=NaN*ones(size(newinlog));
    end
    newdata=newinlog;
    newfilelist=[];
elseif max(size(newinfo))==3
    c=cell(size(reffilelist,1),1);
    newinlog=SampleLog.(newinfo{1});
    newfilelist=cellstr(newinfo{2});
    newdata=newinfo{3};
else
    error('invalid input: newinfo')
end

%locate ref file in the log and get corresponding new file name/data
disp('reconciling reference data with log')
for i=1:length(reffilelist)
    fn=reffilelist(i);
    if ~casesen
        I = strcmpi(fn, refinlog);
    else
        I = strcmpi(fn, refinlog);
    end
    iin=find(I==1);
    if ~isempty(iin)
        try
        c(i)=newinlog(iin);
        catch
            c(i,:)=newinlog(iin);
        end
    else
        fprintf(['Could not locate sample at row ' int2str(i) ' of the sample log\n'])
        errmsg=[char(fn) ' is not found in the filelist.'];
        error(errmsg)
    end
end

%search for new file name in the new filelist
disp('reconciling new data with reference data')
if ~isempty(newfilelist)
    v=NaN*ones(size(newfilelist,1),1);
    for i=1:length(c)
        fn=c(i);
        if ~casesen
            I = strcmpi(fn, newfilelist);
        else
            I = strcmpi(fn, newfilelist);
        end
        iin=find(I==1);
        if ~isempty(iin)
            v(i)=iin;
        else
            fprintf(['Could not locate sample at row ' int2str(i) ' of the reference filelist\n'])
            errmsg=[char(fn) ' is not found in the new filelist.'];
            error(errmsg)
        end
    end
    
    %align
    if length(size(newdata))==2
        A=newdata(v,:);
    elseif length(size(newdata))==3
        A=newdata(v,:,:);
    end
else
    if isnumeric(newinlog)
        c=c(~isnan(c));
        A=c;
        c=cellstr(num2str(c));
    else
      A=c;
    end
end

%disp
disp('')
disp('Check that the lists of files match correctly!')
pause(1)
%size(reffilelist),size(c)
disp('    REFERENCE               ALIGNED')
disp([cellstr(reffilelist) c])

end
