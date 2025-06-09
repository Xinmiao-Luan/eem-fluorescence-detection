function [FMax,B,C,varargout]=modeloutmac(data,f,excelname,varargin)
%
% <strong>Syntax:</strong>
%
%   [FMax,B,C,varargout]=<strong>modeloutmac</strong>(data,f,excelname,varargin)
%
% <a href="matlab: doc modeloutmac">help for modeloutmac</a> <- click on the link


%Export model loadings and validation results to an Excel file. Two to four
%worksheets will be created as follows:
%  "ModelfReport" contains the information to be reported about the model
%     (e.g. in publications),including preprocessing, split-validation 
%     design and results, and model loadings for the Ex and Em modes.
%  "ModelfLoadings" with  the fluorescence intensity of each
%     component in each sample ("FMax"), and the spectral loadings of each
%     component in the Ex and Em modes.
%  "ModelfMetadata"(optional) contains information about the samples.
%  "ModelfFmaxProjected" (optional) contains Fmax for the full dataset.
%
% USEAGE:
%           [FMax,B,C,FMaxFull,Proj]=modelout(data,f,excelname,fullds,metadata)
%
% INPUTS
%      data: data structure containing model for exporting. The model with
%            f components must be located in data.Modelf . Other compulsary
%            fields in data (data.field) are X,Ex,Em. Note that if the
%            model was generated externally of drEEM, then the contents of
%            fields that track dataset operations (e.g. sample removal,
%            splitting operations, validation methods) may be absent or may
%            contain inaccurate data.
%         f: number of components in model.
% excelname: name and location of the excel file to create. Note
%            that if a file of the same name with the above three sheets
%            already exists, data in these sheets will be overwritten.
%    fullds: (optional) the full dataset for projecting on the model.
%            if the full dataset is supplied , it will be projected
%            on the PARAFAC model. This will produce Fmax for
%            all samples including those that were excluded during model
%            building. Note it is critical that the two datasets
%            data.X and fullds.X differ only in numbers of samples, so
%            have identical Ex and Em wavelengths and were pretreated 
%            in the same way (including smoothing).
%            - if fullds has different Ex and Em dimensions,
%              an error will result.
%            - if data.Smooth is not the same as fullds.Smooth, a warning
%              will be given but the projection will be attempted. If this 
%              happens, check that Fmax is close to identical for samples
%              common to data.X and fullds.X. If the wrong dataset was given, 
%              if scatter was treated differently, or if the wrong model 
%              constraints are used, the model may fail to converge or else  
%              may produce wrong values for Fmax.
%              When projecting a dataset, model constraints and convergence
%              criteria will be extracted from data (if present), or 
%              else will can be input during execution of the function.
%  metadata: (optional)
%            list case-sensitive names of metadata fields for exporting
%            as {M1,M2,M3,...}
%
% OUTPUTS
%   outputs are matrices of size [f x nSample] containing for each sample
%   and each component:
%      Fmax: maximum fluorescence intensity of each PARAFAC component
%            calculated during model export.
%         B: emission spectra
%         C: excitation spectra
%  FmaxFull: (optional) Fmax for the full dataset (outliers reincluded)
%      Proj: (optional)  data structure with the projected scores
%            and loadings in Modelf_P
%
% Examples:
%  modelout(val5,5,'MyPARAFACresults.xls');
% [FMax,B,C]=modelout(val5,5,'5compmodel.xlsx',[],{'cruise','site','date'});
% [FMax,B,C,FMaxFull,Proj]=modelout(val6,6,'fullmodel.xlsx',Xs);
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% modeloutmac: Copyright (C) 2014 Kathleen R. Murphy
% $ Version 0.3.0 $ Sep 2015 $ Third Release
% $ Version 0.2.0 $ June 2014 $ Second Release
% $ Version 0.1.0 $ May 2013 $ First Release
%
% bug fix 2015/09- fixed formatting of exported data for excluded samples
%
% Water Environment Technology
% Chalmers University of Technology
% Department of Civil and Environmental Engineering
% Water Environment Technology
% Sven Hultins gata 8
% 412 96 Göteborg
% murphyk@chalmers.se
%
% $ Version 0.2.0 $ May 2014 $ Second Release


xlworkaround(excelname);            
sheetinlang='Sheet'; % EN: Sheet, DE: Tabelle, etc. (Language dependent)
try
    v=ver;
    toolbox=v(contains({v.Name},'drEEM Toolbox')).Version;
catch
    toolbox='drEEM 0.6.0 or higher';
end

narginchk(3,5);
projectmodel=false;
FMaxFull=[];metaflds=[];ndel=0;iremove=[];

%The model
[excelFilePath,fn,ext] = fileparts(excelname);
if isempty(excelFilePath)
    excelFilePath=pwd;
end
thefullfile=fullfile(excelFilePath, fn);
excelFileName=[thefullfile ext];
themodel=['Model' int2str(f)];
M=getfield(data,{1,1},themodel);
A=M{1};B=M{2};C=M{3};
nSample=size(A,1);
nEm=size(B,1);
nEx=size(C,1);
nComp=size(B,2);
BMax=max(B);
CMax=max(C);
report=[data.Ex C;data.Em B];

%Options
if nargin>3 
    errormsg1=false;
    fullds=varargin{1};
    if isstruct(fullds)
        Xb=fullds.X;
        Db=size(Xb);
        if ~isfield(fullds,{'X','Em','Ex'})
            error('The full dataset must at minimum contain the fields: X, Em and Ex')
        end
        if or(isfield(data,'Smooth'),isfield(fullds,'Smooth'))
            if ~and(isfield(data,'Smooth'),isfield(fullds,'Smooth'))
                disp('Incomplete records in data.Smooth and fullds.Smooth. Press any key to continue or CTRL+C to cancel.')
                warning('modelout:Smooth1','Projected PARAFAC scores may be inaccurate!');
                pause
            else
                if ~strcmp(fullds.Smooth,data.Smooth)
                    disp('Incompatible data about scatter removal in data.Smooth and fullds.Smooth. Press any key to continue or CTRL+C to cancel.')
                    warning('modelout:Smooth2','Projected PARAFAC scores may be inaccurate!');
                    pause
                end
            end
        end
        if ~and(isequal(Db(2),nEm),isequal(Db(3),nEx))
            warning('FullDS:X','The matrices in data.X and fullds.X are of incompatible sizes');
            errormsg1=true;
        end
        if ~isequal(fullds.Em,data.Em) %emission
            warning('FullDS:Em','Emission wavelengths for data.X and fullds.X are incompatible');
            errormsg1=true;
        end
        if ~isequal(fullds.Ex,data.Ex) %emission
            warning('FullDS:Ex','Excitation wavelengths for data.X and fullds.X are incompatible');
            errormsg1=true;
        end
        if errormsg1
            error('Can not recover scores for outlier samples due to incompatible wavelengths in the specified full dataset')
        else
            if isfield(data,'Val_Constraints')
                if strcmp(data.Val_Constraints,'nonnegativity')
                    const=[2 2 2];
                elseif strcmp(data.Val_Constraints,'unconstrained')
                    const=[0 0 0];
                end
            else
                const=[0 0 0];
                y=input('Type 1 if the PARAFAC model for exporting was derived using a non-negativity constraint');
                if y==1
                    const=[2 2 2];
                end
            end
            if isfield(data,'Val_ConvgCrit')
                cc=data.Val_ConvgCrit;
            else
                y=input('Input the convergence criterion used to develop the PARAFAC model, or press enter for the default value (1e-6) ');
                if ~isempty(y)
                    cc=y;
                else
                    cc=1e-6;
                end
            end
            projectmodel=true;
        end
    else
        if ~isempty(fullds)
            error('fullds must be a data structure containing X,Em,Ex')
        end
    end
    if nargin>4 %metadata
        metaflds=varargin{2};
        if isempty(metaflds)
            error('[] is not a valid input value for metafields')
        else
        compare=metaflds(isfield(data,metaflds)==0);
        end
        if  size(compare,2)>0
            for i=1:size(compare,2)
                warning(['Metadata field not found:   ' char(compare(i)) ])
            end
            error('One or more metadata field names not found. Note names are case-sensitive')
        else
            MD=cell(1,size(metaflds,2));
            for i=1:size(metaflds,2)
                mf=char(metaflds{i});
                MD{i}=data.(metaflds{i});
                if ~isequal(size(MD{i},1), nSample)
                    error([mf ' does not contain metadata (1 value per sample)']);
                end
            end
        end
    end
end

%Generate PARAFAC Loadings, Fmax
if sum(isfield(data,cellstr(char('backupX','i'))))==2
    nfull=size(data.backupX,1);
    ndel=nfull-nSample;
    iremove=setxor(1:nfull,data.i);
    if size(iremove,1)>size(iremove,2)
        iremove=iremove';
    end
    indices=data.i;
else
    nfull=nSample;
    indices=(1:nSample)';
end
FMax=NaN*ones(nSample,nComp);
for i=1:nSample
    FMax(i,:)=(A(i,:)).*(BMax.*CMax);
end
if projectmodel
    disp('Calculating projected Fmax values....')
    FMaxFull=rand(Db(1),nComp);
    forced=nwayparafac(Xb,nComp,cc,const,{FMaxFull;B;C},[0 1 1]);
    Afull=forced{1};
    for i=1:Db(1)
        FMaxFull(i,:)=(Afull(i,:)).*(BMax.*CMax);
    end
end

%Set up Excel Sheets
reportsheet=[themodel 'Report'];
loadsheet=[themodel 'Loading'];
metasheet=[themodel 'Metadata'];
optsheet=[themodel 'FmaxProjected'];

AtoZ=char(97:122);
removedefaultsheets=true;
sheets2go=cellstr(char([sheetinlang '1'],[sheetinlang '2'],[sheetinlang '3']));

%Warnings if overwriting existing data
if exist(excelFileName,'file')
    [~, desc, ~] = xlsfinfo(excelFileName);
    s2g=char(reportsheet,loadsheet,metasheet,optsheet,[sheetinlang '1'],[sheetinlang '2'],[sheetinlang '3']);
    sheets2go=(intersect(desc,cellstr(s2g)))';
    if ~isempty(sheets2go)
        warning('modelout:overwrite','This action may overwrite data in existing worksheet/s. Press any key to continue or CTRL+C to cancel.');
        pause;
        rmsheetwin(thefullfile,sheets2go,1)
        removedefaultsheets=false;
    end
end

valfields=cellstr(char('Split_Style','Split_NumBeforeCombine',...
    'Split_NumAfterCombine','Split_Combinations','Split_nSample','Split_AnalRuns',...
    'Split_PARAFAC_options','Split_PARAFAC_constraints','Split_PARAFAC_convgcrit',...
    'Split_PARAFAC_Initialise','Val_ModelName', 'Val_Source', 'Val_Err', ...
    'Val_It', 'Val_Result','Val_Splits', 'Val_Comparisons','Val_ConvgCrit',...
    'Val_Constraints','Val_Initialise','Val_Core',...
    'Val_PercentExpl','Val_CompSize','Val_Preprocess'));

exemhead=cellstr([repmat('Ex',[data.nEx,1]) ;repmat('Em',[data.nEm,1])]);
FmaxColHead=cellstr([repmat('Fmax',[nComp,1])  num2str((1:nComp)')])';
ExColHead=cellstr([repmat('Ex',[nComp,1])  num2str((1:nComp)')])';
EmColHead=cellstr([repmat('Em',[nComp,1])  num2str((1:nComp)')])';
CompColHead=cellstr([repmat('Comp',[nComp,1])  num2str((1:nComp)')])';

sepcol=2;
loadheads=cellstr(['i' FmaxColHead blanks(sepcol) 'Ex' ExColHead blanks(sepcol) 'Em' EmColHead]);

disp('Writing to Excel. This may take a few minutes......')

disp('Exporting Info....')
xlwrite(excelFileName,{'PARAFAC Model Report '},reportsheet,'A1');
xlwrite(excelFileName,{'Info'},reportsheet,'A3');
xlwrite(excelFileName,{'Toolbox'},reportsheet,'B4');
xlwrite(excelFileName,{toolbox},reportsheet,'C4');
xlwrite(excelFileName,{'Date'},reportsheet,'B5');
xlwrite(excelFileName,{datestr(now)},reportsheet,'C5');

disp('Exporting Dataset Description....')
PPHeaders=cellstr(char('nSample - full dataset','nSample - modeled dataset',...
    'No. excluded samples','Excluded samples -indices','Scatter Removal','Zapped (Samples,EmRange,ExRange)','Fluorescence unit','Scaling'));
PPData={nfull,nSample,ndel,iremove};
ppOptNames=cellstr(char('Smooth','Zap','IntensityUnit','Preprocess'));
xlwrite(excelFileName,{'Preprocessing'},reportsheet,'A7');
xlwrite(excelFileName,PPHeaders,reportsheet,'B8');
k=8;
for i=1:4
    if ~isempty(PPData{i})
        if ischar(PPData{i})
            xlwrite(excelFileName,cellstr(PPData{i}),reportsheet,['C' num2str(k)]);
        else
            xlwrite(excelFileName,PPData{i},reportsheet,['C' num2str(k)]);
        end
    end
    k=k+1;
end
for i=1:length(ppOptNames)
    ppvarin=ppOptNames{i};
    if isfield(data,ppvarin)
        if ischar(data.(ppvarin))
            xlwrite(excelFileName,cellstr(data.(ppvarin)),reportsheet,['C' num2str(k)]);
        else
            xlwrite(excelFileName,data.(ppvarin),reportsheet,['C' num2str(k)]);
        end
    end
    k=k+1;
end
k=k+1;

disp('Exporting Model Summary....')
MODHeaders=cellstr(char('No. PARAFAC components','No. Ex wavelengths','No. Em wavelengths'));
xlwrite(excelFileName,{'PARAFAC model'},reportsheet,['A' num2str(k)]);
xlwrite(excelFileName,MODHeaders,reportsheet,['B' num2str(k+1)]);
xlwrite(excelFileName,[nComp;nEx;nEm],reportsheet,['C' num2str(k+1)]);
k=k+3;
MOptNames=cellstr(char('OutlierTest_convgcrit','OutlierTest_constraints',...
    [themodel 'err'],[themodel 'it'],[themodel 'core'],[themodel 'source'],...
    [themodel 'convgcrit'],[themodel 'constraints'],[themodel 'initialise'],...
    [themodel 'percentexpl'],[themodel 'compsize'],[themodel 'preprocess']));
for i=1:length(MOptNames)
    Mvarin=MOptNames{i};
    if isfield(data,Mvarin)
        xlwrite(excelFileName,cellstr(Mvarin),reportsheet,['B' num2str(k)]);
        if ischar(data.(Mvarin))
            xlwrite(excelFileName,cellstr(data.(Mvarin)),reportsheet,['C' num2str(k)]);
        else
            xlwrite(excelFileName,data.(Mvarin),reportsheet,['C' num2str(k)]);
        end
        k=k+1;
    end
end
k=k+2;

disp('Exporting Validation Report....')
valheader=false;
kplus=0;
for i=1:size(valfields,1)
    outname=valfields{i};
    if isfield(data,outname)
        if ~valheader
            xlwrite(excelFileName,{'Validation'},reportsheet,['A' int2str(k)]);
            valheader=true;
        end
        k=k+kplus+1;
        kplus=0;
        outdat=data.(outname);
        if ischar(outdat)
            kplus=size(outdat,1)-1;
            outdat=cellstr(outdat);
        end
        namecell=['B' int2str(k)];
        xlwrite(excelFileName,cellstr(outname),reportsheet,namecell);
        datcell=['C' int2str(k)];
        xlwrite(excelFileName,outdat,reportsheet,datcell);
        if i==size(valfields,1)
            k=k+2;
        end
    else
        warning('modelout:ExportField',['No data located for: ' char(outname)])
    end
end

xlwrite(excelFileName,{'Spectra'},reportsheet,['A' int2str(k)]);
xlwrite(excelFileName,['mode' 'nm' CompColHead],reportsheet,['B' int2str(k+1)]);
xlwrite(excelFileName,exemhead,reportsheet,['B' int2str(k+2)]);
xlwrite(excelFileName,report,reportsheet,['C' int2str(k+2)]);

%Export Loadings
disp('Exporting Loadings....')
xlwrite(excelFileName,loadheads,loadsheet,'A1');
xlwrite(excelFileName,[indices FMax],loadsheet,[AtoZ(1) '2']);
xlwrite(excelFileName,[data.Ex C],loadsheet,[AtoZ(1+nComp+sepcol) '2']);
xlwrite(excelFileName,[data.Em B],loadsheet,[AtoZ(1+2*(nComp+sepcol)) '2']);

%Export Fmax for full dataset
if projectmodel
    disp('Exporting projected Fmax ....');
    xlwrite(excelFileName,{'Fmax for full dataset including samples excluded from model building'},optsheet,'A1');
    xlwrite(excelFileName,{'i'},optsheet,'A2');
    xlwrite(excelFileName,FmaxColHead,optsheet,'B2');
    xlwrite(excelFileName,[(1:Db(1))' FMaxFull],optsheet,'A3');
end
   
%Remove empty sheets   
if removedefaultsheets
    disp('Removing extra data sheets....')
    rmsheetwin(thefullfile,sheets2go,0)
end

%Export Metadata
if nargin>4
    disp('Exporting metadata....')
    for i=1:size(metaflds,2)
        xlwrite(excelFileName,metaflds,metasheet,'B1');
        xlwrite(excelFileName,MD{i},metasheet,[AtoZ(i+1) '2']);
    end
    xlwrite(excelFileName,(1:nSample)',metasheet,'A2');
end

if nargout>3
    varargout(1)={FMaxFull};
    if nargout>4
        proj=fullds;
        projname=[themodel '_P'];
        proj.(projname)=forced;
        varargout(2)={proj};
    end
end

disp('Finished!')
end

function rmsheetwin(thefullfile,sheetNames,doclear)
%Delete empty sheets or clear existing data (Windows only)

if ~or(ismac,isunix)
    % Open Excel file.
    objExcel = actxserver('Excel.Application');
    objExcel.Workbooks.Open(thefullfile);
    
    % Delete or clear sheets.
    if doclear
        cellfun(@(x) objExcel.ActiveWorkBook.Worksheets.Item(x).Cells.Clear, sheetNames);
    else
        try
            cellfun(@(x) objExcel.ActiveWorkbook.Worksheets.Item(x).Delete, sheetNames);
        catch ME
            disp(ME); % Do nothing.
        end
    end
    
    % Save, close and clean up.
    objExcel.ActiveWorkbook.Save;
    objExcel.ActiveWorkbook.Close;
    objExcel.Quit;
    objExcel.delete;
end
end

function xlworkaround(FN)
S = which('xlwrite');
if isempty(S)
    disp('       ');
    warning('xlwrite not on the MATLAB path. Observe the following instructions.');
    disp('       ');
    disp('1. Download ''XLWRITE'' (File ID: #38591 by Alec de Zegher) ')
    disp('   from the MATHWORKS file exchange, at');
    disp('      http://www.mathworks.com/matlabcentral/fileexchange/');
    disp('       ');
    disp('      *** do not confuse xlwrite with xlswrite (MATLAB built-in function)! ***');
    disp('       ');
    disp('2. Unzip the folder and place it somewhere general ')
    disp('      e.g C:/Program Files/MATLAB/R2012b/extra_toolboxes/xlwrite')
    disp('3. Add the xlwrite path with subfolders to your MATLAB path. Save.')
    disp('4. Next, make the XLWRITE folder your current directory, then follow  ');
    disp('   the XLWRITE help (notes) to install the POI.');
    disp('5. Finally, change the current directory back to the location')
    disp('   where you want this file to be written, for example')
    disp('   >> cd ''C:/Data/EEMs/Project1/ ''  ');
    disp('6. Now try executing the code again.');
    disp('   ** Note: The excel file must be closed or you will not be able to write to it.');
    disp('       ');
    error('MATLAB:No_MSExcel','Export to MS Excel spreadsheets not currently supported.')
else
    disp('Located xlwrite. Performing a quick test to see if it works.')
    try
        xlwrite(FN,{'xlwrite: Success!'},'Test_XLWRITE','A1');
        disp('Passed.')
    catch
        warning('Test failed, probably because the *.xlsx file is open or because xlwrite is not correctly installed')
        warning('Will try adding POI library to the javapath...')
        currentdir=cd;
        try
            cd(S);
            javaaddpath('poi_library/poi-3.8-20120326.jar');
            javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
            javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
            javaaddpath('poi_library/xmlbeans-2.3.0.jar');
            javaaddpath('poi_library/dom4j-1.6.1.jar');
            disp('POI library is installed')
            cd(currentdir)
        catch
            error('MATLAB:xlwrite','Could not access POI library directory and install files. Could not generate *.xlsx file')
        end
        try
            xlwrite(FN,{'xlwrite: Success!'},Sh,'A1');
            disp('Passed. Try executing the code again.')
        catch %#ok<*CTCH>
            error('MATLAB:xlwrite','Writing to an Excel Spreadsheet failed for undetermined reasons. See earlier warnings for clues.')
        end
    end
end
end

