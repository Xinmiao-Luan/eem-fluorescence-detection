function [X,Emmat,Exmat,filelist,outdata]=readineems(type,format,range,headers,display,outdat)
%
% <strong>Syntax</strong>
%   [X,Emmat,Exmat,filelist,outdata]=<strong>readineems</strong>(type,format,range,headers,display,outdat)
%
% <a href="matlab: doc readineems">help for readineems</a> <- click on the link

% Read in files of specified format (xls,xlsx,dat,txt or csv) in the current directory
% and convert to a 3D matrix of EEMs. Remove non-EEM files of specified format
% from the current directory before running.
% The files will be loaded in alphabetical order (see outputed file list).
% Input EEMs must be of the same size, type and format with same Ex and Em.
%
% USEAGE:
%    [X,Emmat,Exmat,filelist,outdata]=readineems(type,format,range,headers,display,outdat)
%
% INPUTS:
%	   type: 1-3, defined as below
%            1:  Ex in columns and Em in rows (usual output from Fluoromax, Hitachi etc.)
%            2:  Ex in columns and Em in rows, with Em headers displayed in every second column (Varian)
%            3:  Ex in columns and Em in rows, Ex in descending order (AquaLog)
%
%	 format: 'csv', 'xls', 'xlsx','dat', 'txt' reads in all files with the
%	         specified extension. To read in only the subset of these files
%	         conforming to a particular naming scheme, specify the pattern
%	         also. e.g. 'Sample*.csv' or 'QS*Waterfall Plot Sample.dat'
%             * by default, 'xls' and 'xlsx' data will be imported from the
%               first sheet in each excel workbook
%             * csv data exported from a Varian fluorometer contain 2 rows of text headers.
%               The Ex wavelengths will be automatically extracted from the first row
%               of the data file (do not include text rows in the range, see below).
%
% 	  range: Specify the range of cells in the raw EEMs that contain numbers (NOT TEXT).
%            For many EEM files, there are no text rows that need be excluded e.g. range ='A1..AA500'
%            For Varian files, exclude the first two rows plus any rows of text below the EEMs, e.g. range ='A3..CD113'
%            For AquaLog files, exclude the first row plus any rows of text below the EEMs, e.g. range ='A2..BV126'
%            The chosen range affects the designation of headers, see below.
%
%	headers: indicate presence/absence for Ex and Em, reflecting data excluded from the range
%           [1 1]:  present for Ex and Em  - this is usual if numerical Ex and Em headers are present and were not excluded from the range
%           [0 1]:  present for Em but not Ex   - this is usual for Varian and AquaLog files once text are excluded from range
%           [1 0]:  present for Ex but not Em
%           [0 0]:  absent for Ex and Em
%
%	display: length of time (in seconds) to show plots of scans on screen (e.g. 1 = display for 1 second)
%             if display=0, no plots will be created.
%
%	  outdat: Specify whether or not to save the data in the current directory
%               0: do not save or export data
%				1: save dataset and export outdata to MS Excel (*.xlsx, see below)
%               2: save dataset and export outdata to *.dat
%                  Saved data are:
%                 (a) outdata: intensities at common wavelength pairs (A,C,M,T,B)
%                 (b) X,Emmat,Exmat,filelist and outdata saved to a matlab file (RawData_FDOM.mat)
%
% OUTPUTS:
%   		X: 3D matrix of EEMs
%   	Emmat: matrix of Em corresponding with individual EEMs
%   	Exmat: matrix of Ex corresponding with individual EEMs (text data excluded and replaced with 1:N)
%    filelist: list of filenames
%     outdata: data that were saved to file, or [] otherwise
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 
%     5, 6557-6566, 2013. DOI:10.1039/c3ay41160e. 
%
% readineems: Copyright (C) 2013,2014 Kathleen R. Murphy
% $ Version 0.3.0 $ Apr 2019. Improved error handling (UJW)
% $ Version 0.2.0 $ June 2014 $ Second Release
% $ Version 0.1.1 $ May 2014 $ Second Release
% $ Updated from FDOMcorr toolbox ver. 1.8
%
% Copyright (C) 2014 KR Murphy,
% Water Environment Technology
% Chalmers University of Technology
% Department of Civil and Environmental Engineering
% Water Environment Technology
% Sven Hultins gata 8
% 412 96 Göteborg
% murphyk@chalmers.se

%% INITIALISE
filelist=[]; X=[];Emmat=[];Exmat=[]; x=[]; outdata=[];datelist=[]; 
EmDummy=0;ExDummy=0;
if display>0
    figure
end

%Obtain list of files conforming to specified format
period=strfind(format,'.');
if ~isempty(period)
    direc = dir(format);
    format=format(period+1:end);
else
    direc = dir(['*.' format]);
end


if size(direc,1)==0
    error(['No files with the specified pattern ''',format,''' were found in the directory'])
end

%Load EEMs
for i=1:size(direc,1)
    
    fprintf([direc(i).name '\n'])
    filelist=char(filelist,direc(i).name); %file list
    datelist(i)=direc(i).datenum; %date list
    %% LOAD FILES, excluding cells that contain text
    if strcmp(format,'csv')==1
        x = dlmread(direc(i).name,',',range);
    elseif strcmp(format,'xls')==1
        x = xlsread(direc(i).name,1,range);
    elseif strcmp(format,'xlsx')==1
        x = xlsread(direc(i).name,1,range);
    elseif strcmp(format,'dat')==1
        x = dlmread(direc(i).name,'\t',range);
    elseif strcmp(format,'txt')==1
        x = dlmread(direc(i).name,'\t',range);
    end
    
    if i==1
        xsize=size(x);
    else
        %check that the size of x is consistent with the first EEM
        if size(x)~=xsize
            fprintf('The size of the input EEMs has changed!\n')
            error('Check raw data file and resize as necessary.')
        end
    end
    
    %% REMOVE HEADERS
    %size(x),x(1:5,1:5)
    if headers(2)==1 %remove the Em headers
        if headers(1)==0 %no Ex headers
            Emmat(:,i)=x(:,1); %#ok<*AGROW>
            if ~isequal(type,2)
                x=x(:,2:end); % remove Em headers from first column
                ExDummy=1;Exmat(i,:)=1:size(x,2);
                if type==3 %AquaLog
                    exheaders=aqlhead(range);
                    try %try to automatically extract the headers from B1..Bn
                        Exmat(i,:)=fliplr(dlmread(direc(i).name,'\t',exheaders));
                        ExDummy=0;
                    catch
                        if i==1,fprintf('could not extract Ex headers automatically; add them manually'),end
                        Exmat(i,:)=1:(size(x,2));
                        %disp(Exmat(i,:))
                    end
                end
            elseif type==2 %varian
                x=x(:,2:2:end); % remove Em headers from every second column of varian file
                TxtLen=6; %Assume excitation wavelengths that are 6 characters long, e.g. 250.00
                try
                    TxtEx=ExtractVarianEx(direc(i).name,TxtLen);
                    %Exmat, pause
                    if isempty(Exmat)
                        Exmat(i,:)=TxtEx;
                    elseif isequal(length(TxtEx),size(Exmat,2))
                        Exmat(i,:)=TxtEx;
                    elseif ~isequal(length(TxtEx),size(Exmat,2))
                        if length(TxtEx)>size(Exmat,2)
                            Exmat(i,:)=TxtEx(1:size(Exmat,2))';
                            warning('These data have been truncated!'),
                            %pause(0.5)
                            %size(TxtEx),TxtEx',size(Exmat)
                            %disp(Exmat(i,:))
                        end
                    end
                catch 
                    ExDummy=1;
                    Exmat(i,:)=1:(size(x,2));
                    disp(Exmat(i,:))
                end
            end
        elseif headers(1)==1 %has Ex headers
            Emmat(:,i)=x(2:end,1);
            if type==1
                Exmat(i,:)=x(1,2:end);
                x=x(2:end,2:end); % remove Em headers from first column
            elseif type==2
                Exmat(i,:)=x(1,2:2:end);
                x=x(2:end,2:2:end); % remove Em headers from every second column of varian file
            end
        end
    elseif headers(2)==0 % no Em headers
        if type==2
            error('type=2 (alternate columns of Em) is not compatible with "no Em headers"');
        end
        EmDummy=1;
        if headers(1)==1 %remove the Ex headers
            Exmat(i,:)=x(1,:);
            x=x(2:end,:);
        elseif headers(1)==0 %
            ExDummy=1;
            Exmat(i,:)=1:(size(x,2));
        end
        Emmat(:,i)=(1:size(x,1))';
    end
    
    if type==3 %put Ex wavelengths in ascending order (AquaLog)
        x=fliplr(x);
    end
    
    X(i,:,:)=x;
    %size(Exmat),size(Emmat), size(x),x(1:5,1:5)
    
    if ~isequal(size(Emmat,1),size(x,1))
        fprintf('Em size is '), size(Emmat),
        fprintf('X size is '), size(x),
        error('Em and X matrix sizes not compatible'),
    end
    if ~isequal(size(Exmat,2),size(x,2))
        fprintf('Ex size is '), size(Exmat),
        fprintf('X size is '), size(x),
        error('Ex and X matrix sizes not compatible'),
    end
    
    %% DISPLAY a plot of each EEM as it is loaded
    if display>0
        %size(x),size(Emmat(:,i)), size(Exmat(i,:))
        contourf(Exmat(i,:),Emmat(:,i),x),
        title(direc(i).name,'interpreter','none')
        pause(display),
    end
end
filelist=filelist(2:end,:); %remove first row of blanks

%Create outdata matrix (assume Em and Ex as for the first file that was read in)
if EmDummy==0
    Em=round(Emmat(:,1)*2)/2; %round to nearest 0.5 nm
elseif EmDummy==1
    moveon=0;
    while moveon==0
        Em = input('Specify the emission wavelength range, e.g. 300:2:600:    ');
        if size(Em,2)==size(x,1)
            Emmat=repmat(Em',[1 size(X,1)]);
            moveon=1;
        else
            fprintf(['Size mismatch: expecting a header of size 1x' num2str(size(x,1)) '. Try again. \n']) ;
        end
    end
end
if ExDummy==0
    Ex=round(Exmat(1,:)*2)/2; %round to nearest 0.5 nm
elseif ExDummy==1
    moveon=0;
    while moveon==0
        Ex = input('Specify the excitation wavelength range, e.g. 250:5:450:    ');
        if size(Ex,2)==size(x,2)
            Exmat=repmat(Ex,[size(X,1) 1]);
            moveon=1;
        else
            fprintf(['Size mismatch: expecting a header of size 1x' num2str(size(x,2)) '. Try again. \n']) ;
        end
    end
end

%% EXPORT DATA
% Export EEM cube to MATLAB file "RawData_FDOM.mat"
% Export data for particular wavelength pairs to an excel spreadsheet "OutData_FDOM.xls"
% adjust or add to the following wavelength selection as required

if outdat>0
    inpairs=[350 450; ...    %C (humic-like)
        250 450; ...         %A (humic-like)
        290 350; ...         %T (tryptophan-like)
        270 304; ...         %B (tyrosine-like)
        320 412];            %M (marine/microbial-like)
    
    %modify the list above to correspond with Ex and Em in this dataset
    outwaves=nearestwave(inpairs,Ex,Em);
    
    outdata=NaN*ones(size(X,1)+2,size(outwaves,1)+1); %first two rows are headers
    outdata(1,2:end)=outwaves(:,1)';    %excitation headers
    outdata(2,2:end)=outwaves(:,2)';    %emission headers
    outdata(3:end,1)=(1:size(X,1))';  %number samples
    for i=1:size(outwaves,1)
        p=X(:,Em==outwaves(i,2),Ex==outwaves(i,1));
        if ~isempty(p)
            outdata(3:end,i+1)=nanmean(p,2);
        end
    end
    
    %% Export data to xlsx (not xls)
    FNout='OutData_FDOM.xlsx';
    if outdat==1
        
        try
            w = actxserver('Excel.Application');     % Fails if Excel not installed
        catch
            w = [];
        end
        %w=[] %test xlwrite
        
        if isempty(w) %Work around for Apple Macs and if MS Excel is not installed
            try
                xlworkaround(FNout,'TestXLwrite');
                %Write outdata matrix to excel file: OutData_FDOM.xls
                fprintf('writing OutData_FDOM.xlsx to your current directory...\n')
                %delete FNout; %Remove old versions of this spreadsheet.
                xlwrite(FNout,outdata,'Raw');
                xlwrite(FNout,cellstr('Ex wave'),'Raw','A1');
                xlwrite(FNout,cellstr('Em wave'),'Raw','A2');
                xlwrite(FNout,(1:size(X,1))','Filenames','A2');
                xlwrite(FNout,cellstr('List of Files'),'Filenames','B1');
                xlwrite(FNout,cellstr(filelist),'Filenames','B2');
            catch
                warning('Could not write to *.xls file. Use outdat=2 instead')
            end
        else
            %Write outdata matrix to excel file: OutData_FDOM.xls
            fprintf(['writing ' FNout ' to your current directory...\n'])
            %delete FNout; %Remove old versions of this spreadsheet.
            xlswrite(FNout,outdata,'Raw');
            xlswrite(FNout,cellstr('Ex wave'),'Raw','A1');
            xlswrite(FNout,cellstr('Em wave'),'Raw','A2');
            xlswrite(FNout,(1:size(X,1))','Filenames','A2');
            xlswrite(FNout,cellstr('List of Files'),'Filenames','B1');
            xlswrite(FNout,cellstr(filelist),'Filenames','B2');
        end
    elseif outdat==2 %write to .dat file
        fileID = fopen([FNout(1:end-5) '.dat'],'w');
        if fileID==-1
            error('readineems:fid',['Could not access  ' FNout '. It may be a permissions error. Close any files with the same name .'])
        end
        C = [cellstr(char('Ex wave','Em wave',filelist)) num2cell(outdata)];
        dline=repmat('%6.4f\t',[1,size(outdata,2)]);
        formatSpec = ['%s\t' dline  '\n'];
        [nrows,ncols] = size(C); %#ok<NASGU>
        for row = 1:nrows
            fprintf(fileID,formatSpec,C{row,:});
        end
        fclose(fileID);
        fprintf([FNout(1:end-5) '.dat has been written to your current directory...\n'])
    end
    
    fprintf('writing RawData_FDOM.mat to your current directory...\n')
    save RawData_FDOM.mat X Emmat Exmat filelist outdata datelist
end
end

%% Extract Ex from Varian text headers
function ExTXT=ExtractVarianEx(fn,TxtLen)
% function ExTXT=ExtractVarianEx(fn,TxtLen)
% Extract wavelengths embedded in text from the first row of a Varian data file
%
% INPUTS
%  fn       = filename
%  TxtLen   = number of characters in each Ex wavelength (typically six, for example, 220.00)
%
% OUTPUTS
%  ExTXT    = Excitation wavelengths extracted from the header
%
% Copyright (C) 2011 KR Murphy,
% Water Research Centre
% The University of New South Wales
% Department of Civil and Environmental Engineering
% Sydney 2052, Australia
% krm@unsw.edu.au
%%%%%%%%%%%%%%%%

fid=fopen(fn);          % Open the file
tline = fgetl(fid);     % Header line containing excitation wavelengths
fclose(fid);            % Close the file

CommaPos=strfind(tline,',');
C=[];
for i=1:TxtLen
    C=[(CommaPos(1:2:end)-i)' C]; %#ok<AGROW>
end
ExTXT=str2num(tline(C)); %#ok<ST2NM>
%%%%%%%%%%%%%%%%
end

function newpair=nearestwave(inpair,Ex,Em)
%find wavelength pairs in Ex and Em as similar
%to those listed in inpair as possible
newpair=inpair;
for i=1:size(inpair,1);
    [~, j1]=min(abs(Ex-inpair(i,1)));
    [~, j2]=min(abs(Em-inpair(i,2)));
    newpair(i,:)=[Ex(j1) Em(j2)];
end
end

function xlworkaround(FN,Sh)
%Work around for xlswrite, which does not work in the following situations:
% 1. The computer is a Mac
% 2. MS Excel is not installed
% 3. Communication between MATLAB and MS Excel fails for another reason.
%The work around is to use xlwrite instead of the built-in xlswrite.
% Download ''XLWRITE'' (File ID: #38591 by Alec de Zegher) from the MATHWORKS
% file exchange, at http://www.mathworks.com/matlabcentral/fileexchange/');
% Follow the instructions to install it, especially remembering to activate
% the POI files according to the help documentation.
warning('MATLAB can not access MS Excel. Checking for alternative route via xlwrite.');
S = which('xlwrite');
if isempty(S)
    warning('xlwrite not on the MATLAB path. Observe the following instructions.');
    disp('1. Download ''XLWRITE'' (File ID: #38591 by Alec de Zegher) ')
    disp('   from the MATHWORKS file exchange, at');
    disp('      http://www.mathworks.com/matlabcentral/fileexchange/');
    disp('2. Unzip the folder and place it somewhere general ')
    disp('      e.g C:/Program Files/MATLAB/R2012b/extra_toolboxes/xlwrite')
    disp('3. Add the xlwrite path with subfolders to your MATLAB path. Save.')
    disp('4. Next, make the XLWRITE folder your current directory, then follow  ');
    disp('   the XLWRITE help (notes) to install the POI.');
    disp('5. Finally, change the current directory back to the location')
    disp('   where you want this file to be written, for example')
    disp('   >> cd ''C:/Data/EEMs/Project1/ ''  ');
    disp('6. Now try executing the code again.');
    disp('');
    error('MATLAB:No_MSExcel','Export to MS Excel spreadsheets not currently supported.')
else
    disp('Located xlwrite. Performing a quick test to see if it works.')
    try
        disp('attempt 1')
        %xlwrite('xlwrite_test.xlsx',{'xlwrite: Success!'},Sh,'A1')
        xlwrite(FN,{'xlwrite: Success!'},Sh,'A1');
        disp('Passed.')
    catch
        warning('Test failed. Try adding POI library to the javapath...')
        currentdir=cd;
        trytoaddpath=true;
        try
            cd(S);
        catch ME
            disp(ME)
            trytoaddpath=false;
        end
        if trytoaddpath
            javaaddpath('poi_library/poi-3.8-20120326.jar');
            javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
            javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
            javaaddpath('poi_library/xmlbeans-2.3.0.jar');
            javaaddpath('poi_library/dom4j-1.6.1.jar');
            cd(currentdir)
            try
                disp('attempt 2')
                xlwrite(FN,{'xlwrite: Success!'},Sh,'A1');
                disp('Passed.')
            catch %#ok<*CTCH>
                error('MATLAB:xlwrite','Writing to an Excel Spreadsheet using xlwrite work-around failed for unknown reasons')
            end
        else
            error('Permisisons issue for xlswrite directory. Add the POI library manually')
        end
    end
end
end
function newrange=aqlhead(range)
%gets the range of AquaLog Ex headers, assuming they start in 'B1'
colno=regexp(range,'\d');
range(1:2)='B1';
newrange=[range(1:end+1-length(colno)) '1'];
end