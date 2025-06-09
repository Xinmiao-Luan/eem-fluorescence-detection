function [SampleLog,H,C]=readlogfile(filename,strings)
%
% <strong>Syntax:</strong>
%
%   [SampleLog,H,C]=<strong>readlogfile</strong>(filename,strings)
%
% <a href="matlab: doc readlogfile">help for readlogfile</a> <- click on the link

%Reads in a log file (*.csv) consisting of columns of data with one row of headers.
%Outputs a structure called SampleLog, with each column of data imported 
%into a different field. The field names are taken from the column headers, 
%ignoring any blanks so that e.g. 'EEM file' becomes 'EEMfile' 
%Note that if columns contain numbers but are read in as strings, it can cause
%blank rows to appear in the SampleLog. Prevent this by
%converting numbers to string in the excel file before importing.
%
% Useage:  [SampleLog,H,C]=readlogfile(filename,strings)
% Inputs
%      filename: name of log file including the file extension
%      strings: specify which columns contain numbers versus text using
%               zero for numbers and 1 for text
% Outputs
%      SampleLog: Data structure containing the imported data from C and H
%              H: Row  1     of the log file (the column headers)
%              C: Rows 2:end of the log file (the data)
%
% Example
% readlogfile('SampleLog_PortSurveyDemo.csv',[0 1 1 1 0 0 1 1 1 1 1 1 0 1 0 1])
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% readlogfile: Copyright (C) 2014 Kathleen R. Murphy
% $ Version 0.2.0 $ May 2014 $ Second Release
%
% Copyright (C) 2014 KR Murphy,
% Water Environment Technology
% Chalmers University of Technology
% Department of Civil and Environmental Engineering
% Water Environment Technology
% Sven Hultins gata 8
% 412 96 Göteborg
% murphyk@chalmers.se


hformats=[repmat('%s',[1,size(strings,2)]) '\n'];
colformats=[];
for i=1:size(strings,2)
    icol=strings(i);
    switch icol
        case 0
            colformats=[colformats '%f ']; %#ok<*AGROW>
        case 1
            colformats=[colformats '%s '];
    end
end
           
fid = fopen(filename);
if fid==-1
    error(['Could not open file: ' filename '. Check the file name and path are correct, and close the file if it is open in another program.']) 
end
H = textscan(fid, hformats,1,'delimiter', ',');
C = textscan(fid, colformats,'delimiter', ',');
fclose(fid);

clear SampleLog
for i=1:size(C,2)
    try
        str=char(H{i});
        str(isspace(str))='';
        SampleLog.(str)=C{i};
    catch ME %#ok<NASGU>
        if isempty(char(H{i}))
            msg=['could not name column #' num2str(i) ' with header : (empty) '];
        else
        msg=['could not name column #' num2str(i) ' with header : ' char(H{i}) ' '];
        end
        warning(msg)
    end
end
%disp(SampleLog)
