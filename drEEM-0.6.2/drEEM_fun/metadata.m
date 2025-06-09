function [C,D]=metadata(DS,varargin)
%
% <strong>Syntax:</strong>
%
%   [C,D]= <strong>metadata</strong>(DS,varargin)
%
% <a href="matlab: doc metadata">help for metadata</a> <- click on the link

%Display metadata for a dataset or its splits. There is one row of metadata
%for each row of data.X (the EEMs). Assuming data.nSample = N, a field
%is assumed to contain metadata if it has N rows of text or numbers.
%Turn off fields you dont want displayed. 
%
%Inputs
%DS : Dataset containing metadata in fields DS.metadata1, Ds.metadata2,
%varargin (optional) 
%    text strings or numbers, where:
%     * Number indicates which splits to show, e.g. 1 or 1:4, 
%        instead of showing metadata for the main dataset
%     * Strings list fields to ignore, e.g. 'date','degC','salinity'
%
%Outputs are written to the screen and include
%     * the field names for the metadata 
%     * the metadata displayed as text
%
%Examples
% metadata(DS)
% metadata(DS,[1,3]) %metadata for splits 1 and 3
% metadata(DS,'ID','date') %don't display data in fields DS.ID, DS.date
% metadata(DS,1:4,'ID') %metadata for splits 1:4 (ignore DS.ID)
%
% metadata: Copyright (C) 2013,2014 Kathleen R. Murphy
% $ Version 0.3.0 $ August 2015 $ Third Release
%
% Water Environment Technology
% Chalmers University of Technology
% Department of Civil and Environmental Engineering
% Water Environment Technology
% Sven Hultins gata 8
% 412 96 Göteborg
% murphyk@chalmers.se
%
%%
splits=1;splitopt=false;
FN=fieldnames(DS);
if nargin>1
    for i=1:nargin-1
        if isnumeric(varargin{i})
            splitopt=true;
            splits=varargin{i};
        else
            FN=setxor(FN,varargin{i});
        end
    end
end

for s=1:length(splits)
    z=splits(s);
    if true(splitopt)
        data=DS.Split(z);
        FNs=fieldnames(data);
        FNi=intersect(FN,FNs);
    else
        data=DS;
        FNi=FN;
    end
    C=cellstr(num2str([1:data.nSample]'));D='#'; %#ok<NBRAK>
    for i=1:length(FNi)
        d=char(deblank(FNi(i,:)));
        c=data.(char(deblank(FNi(i,:))));
        if iscellstr(c)
            if size(c,1)==data.nSample
                C=[C c]; %#ok<*AGROW>
                D=char(D, d);
            end
        elseif isnumeric(c)
            if and(size(c,1)==data.nSample,size(c,2)==1)
                C=[C cellstr(char(num2str(c)))];
                D=char(D, d);
            end
        end
    end
    if true(splitopt)
        disp(['Split number ' num2str(z)])
    end
    D=cellstr(D);
    disp(cellstr(D))
    disp(C)
end
end