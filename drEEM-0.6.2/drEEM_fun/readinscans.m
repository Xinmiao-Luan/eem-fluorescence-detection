function [S,W,wave,filelist]=readinscans(type,format,XLrange,display,outdat,pattern)

% Read files of specified format (csv, txt, dat, xls, xlsx) in the current directory
% and convert to a single 2D matrix of scans. Remove unwanted files of the same format
% from the current directory before running this program.
% The files will be loaded in alphabetical order (see output: filelist).
% All scans must be of the same size, type and format
%
% USEAGE:
%           [S,W,wave,filelist]=readinscans(type,format,XLrange,display,outdat,pattern)
% INPUTS:
%       type: Up to 10 characters describing the type of data scans, used to name the output file. For example: 
%            'Absorbance':  Absorbance scans
%            'UV350':  UV scans at 350 nm excitation
%
%     format: 'csv','xls','xlsx','txt','dat' or the special formats that are the same
%             as these but end in numbers e.g. 'dat_1_10', csv_2_5. See note 3 below.
%             Note 1. If importing from .txt or .dat files, first open them with Excel to determine
%                     what range of cells (XLrange) are appropriate.
%             Note 2. For xls or xls formats, data will be imported from the first worksheet of the excel file
%             Note 3. If the files have >2 columns and you wish to import
%                     only some colums, representing [wavelengths, data] then
%                     specify the indices of columns at the end of the
%                     format. Subtract columns that were excluded from the XLrange.
%                     e.g. 'dat_1_10': wavelengths from col. 1 and data from col. 10.
%                     e.g. 'csv_2_5': wavelengths from col. 2 and data from col. 5.
%
%    XLrange: Range of cells in the raw EEMs that contain numbers (see Note 1 above for 'txt' formats).
%             e.g. range ='A3..B500' or range ='A2..BB3' (or range = [] for .txt files only)
%
%    display: length of time to show plots of scans on screen (e.g. 1 = display for 1 second)
%             if display=0, no plots will be created.
%
%     outdat: 1 to write S, W, filelist and wave to an Excel file in the current directory; 
%                 Data are saved in "OutData_xx.xlsx": where xx is the contents of 'type' above.
%             2 to write S, W, filelist and wave to dat files
%               (this is faster than writing to excel)
%
%    pattern: a text string in the file name. E.g. 'QS*' means only filenames beginning with 'QS' will be read in.
%             If no pattern is specified the default will be to read in all files of the specified extension,
%             i.e. the same as reading in e.g. '*.csv' or '*.txt'.
%
% OUTPUTS:
%          S: A 2D matrix of scans
%          W: A 2D matrix of wavelengths corresponding with individual scans
%       wave: if all the rows of W are equal within 0.5 nm, wave is a single row of wavelengths rounded to 0.5nm
%   filelist: the list of filenames in the order they appear in S, W and wave.
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% readinscans: Copyright (C) 2013,2014 Kathleen R. Murphy
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


%% Initialize
narginchk(5,6)
nargoutchk(4,4)

%If files are in Excel format, delete old version of OutData.xls to avoid conflict.
if strcmp(format,'xls')==1;
    FN=['OutData_' type '.xls'];
    fprintf(['Any existing files of name ' FN ' in the current directory will be replaced. Press any key to continue...\n'])
    pause;
    eval(['delete ' FN])
end

%Check whether to reduce matrices to two columns
SubByCol=false;
if nargin>=5
    prescols = regexp(format,'_', 'once');
    if ~isempty(prescols)
        SubByCol=true;
        [format,col1,col2]=numext(format);
    end
    if nargin == 5
        direc = dir(['*.' format]);
    elseif nargin == 6
        direc = dir([pattern '.' format]);
    end
end

%initialise S and W from RangeIn
W=[];S=[];
filelist=[];frows=[];

%% Load Files
for i=1:size(direc,1)
    
    fprintf([direc(i).name '\n'])
    filelist=char(filelist,direc(i).name); %file list

    % LOAD FILES, excluding cells that contain text
    if strcmp(format,'csv')==1;
        x = dlmread(direc(i).name,',',XLrange);
    elseif strcmp(format,'txt')==1;
        x = dlmread(direc(i).name,'\t',XLrange);
    elseif strcmp(format,'xls')==1;
        x = xlsread(direc(i).name,1,XLrange);
    elseif strcmp(format,'xlsx')==1;
        x = xlsread(direc(i).name,1,XLrange);
    elseif strcmp(format,'dat')==1;
        x = dlmread(direc(i).name,'\t',XLrange);
    end
    
    if SubByCol
         x = x(:,[col1 col2]);
    end

    % convert to row vectors
    xsize=size(x);
    sizeerror=false;
    if min(xsize)==1 % single data pair
        if and(xsize(1)==1,xsize(2)==2); %x=[wave data]
            x=x';
        else
            sizeerror=true;
        end
    elseif min(xsize)==2 %data are in rows or in columns
        if xsize(1)>2; %data are in columns
            x=x';
        elseif xsize(1)==xsize(2); %data in rows or in columns
            if isempty(frows);
                fprintf('input data are 2 rows x 2 colums, more information needed... \n')
                while isempty(frows)
                    frows=input('type 1 if wavelengths are in rows, or 2 if wavelengths are in columns: ');
                    if ismember(frows,[1;2])==true
                    else frows=[];
                    end
                end
            end
            if frows==2;
              x=x';              
            end
        end
    elseif min(xsize)>2
        sizeerror=true;
    end
    if sizeerror==true;
        fprintf(['size of input data is ' num2str(xsize(1)) ' rows and ' num2str(xsize(2))  ' columns \n'])
        error(['Inappropriate data range. Scans must consist of a row or column ' ...
            'of wavelengths followed by 1 row or column of data']);
    end
    x=(sortrows(x',1))';    
    xsize=size(x);
    
    %check that x has no more than two rows (or columns)
    if i==1;
        xsize1=xsize;
    else
        %check that the size of x is consistent with the first scan
        if size(x)~=xsize1
            fprintf('The size of the input scans has changed!\n')
            fprintf(['size of first scan was ' num2str(xsize(1)) 'rows and ' num2str(xsize(2))  ' columns'])
            fprintf(['size of current scan is ' num2str(size(x,1)) 'rows and ' num2str(size(x,2))  ' columns'])
            error('Check raw data file and resize as necessary.')
        end
    end
    
    %% Assemble matrices
    %size(x),pause
    W(i,:)=x(1,:); %#ok<AGROW>
    x=x(2,:);
    S(i,:)=x; %#ok<AGROW>
    
    %% DISPLAY a plot of each scan as it is loaded
    if display>0;
        plot(W(i,:),S(i,:),'ko'),
        title(direc(i).name,'interpreter','none')
        if i>1
            hold on
            plot(W(1:i,:)',S(1:i,:)','c-'),
            legend('current scan','prior scans')
            hold off
        end
       pause(display),
       if i<size(direc,1)
       close
       end
    end
end
%remove leading row of blank filenames
filelist=filelist(2:end,:); 

%Check for consistency among wavelengths
wave=NaN*ones(1,size(S,2));
rndW=round(W*2)/2;
rndW1=repmat(rndW(1,:),[size(rndW,1) 1]);
Wdiff=rndW-rndW1;
if max(max(Wdiff))==0;
    wave=rndW(1,:);
else
    fprintf('\n')
    fprintf('Warning: not all input wavelengths were the same! \n')
    fprintf('Check the rows of output variable W for non-standard wavelengths.\n')
    fprintf('The output variable "wave" is a vector of NaNs.\n')
end

%% EXPORT DATA 
% Export scans to an excel spreadsheet
if outdat==1;
    xlver='.xlsx'; %if exporting to .xls note that no more than 256 columns can be exported; not Mac compatible. 
    FNout=['OutData_' type xlver];
    
    try
        w = actxserver('Excel.Application');     % Fails if Excel not installed
    catch
        w = [];
    end
    %w=[] %test xlwrite
    
    if isempty(w) %Work around for Apple Macs and if MS Excel is not installed
        warning('Workaround will be implemented for Excel/Mac incompatibility. To avoid this message use outdat=2!')
        pause(2);
        try
            xlworkaround(FNout,'TestXLwrite');
            xlwrite(FNout,wave,type,'A1');
            xlwrite(FNout,S,type,'A2');
            xlwrite(FNout,W,'Wavelengths');
            fprintf([FNout ' has been written to your current directory...\n'])
       catch
            warning('Could not write to *.xlsx file. Use outdat=2 instead')
        end
    else
        xlswrite(FNout,wave,type,'A1');
        xlswrite(FNout,S,type,'A2');
        xlswrite(FNout,W,'Wavelengths');
        xlswrite(FNout,cellstr(filelist),'Filenames','A2')
        fprintf([FNout ' has been written to your current directory...\n'])
    end

elseif outdat==2;
    fileID = fopen(['OutData_' type '_scan.dat'],'w');
    if fileID==-1;
        error('readinscans:fid',['Could not access  ' FNout '. It may be a permissions error. Close any files with the same name .'])
    end
            
    C1 = {'filenames',wave};
    firstline=repmat('%6.1f\t',[1,size(S,2)]);
    fprintf(fileID,['%s\t' firstline '\n'],C1{1,:});
    C = [cellstr(filelist) num2cell(S)];
    dline=repmat('%6.4f\t',[1,size(S,2)]);
    formatSpec = ['%s\t' dline  '\n'];
    [nrows,ncols] = size(C); %#ok<NASGU>
    for row = 1:nrows
        fprintf(fileID,formatSpec,C{row,:});
    end
    fclose(fileID);
    
    fileID = fopen(['OutData_' type '_inwave.dat'],'w');
    C1 = {'filenames','wavelengths read in'};
    fprintf(fileID,'%s\t %s\t \n',C1{1,:});    
    C = [cellstr(filelist) num2cell(W)];
    dline=repmat('%6.4f\t',[1,size(W,2)]);
    formatSpec = ['%s\t' dline  '\n'];
    [nrows,ncols] = size(C); %#ok<NASGU>
    for row = 1:nrows
        fprintf(fileID,formatSpec,C{row,:});
    end
    fclose(fileID);
    
    fprintf([['OutData_' type '_scan.dat'] ' has been written to your current directory...\n'])
    fprintf([['OutData_' type '_inwave.dat'] ' has been written to your current directory...\n'])
end
end

function varargout=numext(str)
%extract the column numbers from the text string
t=str;
i1 = regexp(str,'_');
i2 = regexp(str,'[0-9]');
if ~isempty(i1)
    t=t(1:i1-1);
    varargout{2}=str2double(str(i2(i2<i1(2))));
    varargout{3}=str2double(str(i2(i2>i1(2))));
elseif length(i1)>2
    error('NUMEXT:String','incompatible string format')
end
varargout{1}=t;
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