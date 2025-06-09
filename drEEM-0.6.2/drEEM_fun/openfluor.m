function openfluor(data, f, filename)
%
% <strong>Syntax:</strong>
%
%   <strong>openfluor</strong>(data, f, filename)
%
% <a href="matlab: doc openfluor">help for openfluor</a> <- click on the link

%Write a PARAFAC model output to text file ready for importing in OpenFluor.
%
% USEAGE:
%     openfluor(data, f, filename)
%
% INPUTS:
%      data: data structure containing model for exporting. The model with
%            f components must be located in data.Modelf . Other compulsary
%            fields in data (data.field) are X,Ex,Em. 
%         f: number of components in model.
% filename:  name of the text file to create, 
%            e.g. 'Model6.txt'.
%
% OUTPUTS:
%         filename.txt will be created 
%         it will include place markers for metadata plus PARAFAC spectra
%         see www.openfluor.org
%         Ex and Em must be whole numbers; if not, the spectra are
%         interpolated.
%
% EXAMPLES:
%     openfluor(LSmodel5, 5, 'Model5.txt')
%     openfluor(LSmodel6, 6, 'C:/Data/MATLAB/PARAFAC/Model6.txt')
%
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% openfluor: Copyright (C) 2013 Kathleen R. Murphy
%
% $ Version 0.1.1 $ Jan 2014 
% $ Version 0.1.2 $ Oct 2014 - bug fix - rounding of Em spectrum

themodel=['Model' int2str(f)];
M=getfield(data,{1,1},themodel);
A=M{1};B=M{2};C=M{3};
nSample=size(A,1);
Ex=data.Ex;
Em=data.Em;
if ~isequal(round(Ex),Ex)
    disp('Ex must take whole number values')
    disp('Press any key to interpolate the spectra, or Cntrl-C to cancel')
    %pause(1)
    Exint=round((max(Ex)-min(Ex))/(length(Ex)-1));
    Exmin=round(min(Ex));
    Exmax=round(max(Ex));
    C=interp1(Ex,C,Exmin:Exint:Exmax,'spline','extrap');
    Ex=(Exmin:Exint:Exmax)';
    disp('Excitation wavelengths have been interpolated to whole numbers')
end
if ~isequal(round(Em),Em)
    disp('Em must take whole number values')
    disp('Press any key to interpolate the spectra, or Cntrl-C to cancel')
    %pause(1)
    Emint=round((max(Em)-min(Em))/(length(Em)-1));
    Emmin=round(min(Em));
    Emmax=round(max(Em));
    B=interp1(Em,B,Emmin:Emint:Emmax,'spline','extrap');
    Em=(Emmin:Emint:Emmax)';
    disp('Emission wavelengths have been interpolated to whole numbers')
end
report=[Ex C;Em B];
RowHead=char(repmat('Ex',[size(C,1),1]),repmat('Em',[size(B,1),1]));

metafields=cellstr(char(...
'name',...
'creator',...
'email',...
'doi',...
'reference',...
'unit',...
'toolbox',...
'date',...
'fluorometer',...
'nSample',...
'constraints',...
'validation',...
'methods',...
'preprocess',...
'sources',...
'ecozones',...
'description'));

metalocs=(char(...
'',...
'',...
'',...
'',...
'',...
'IntensityUnit',...
'',...
'date',...
'',...
'nSample',...
[themodel 'constraints'],...
'',...
'',...
'Preprocess',...
'',...
'',...
''));
d=datestr(now);
fid = fopen(filename, 'w');
fprintf(fid, '%s\t\n', '#');
fprintf(fid, '%s\t\n', '# Fluorescence Model');
fprintf(fid, '%s\t\n', '#');
for i=1:size(metafields,1)
    fprintf(fid, '%s\t', metafields{i,:});
    switch i
        case {6,11,14}
            try
                fprintf(fid, '%s\n', data.(deblank(metalocs(i,:))));
            catch %#ok<*CTCH>
                fprintf(fid, '%s\n', '?');
            end
        case 8
            fprintf(fid, '%s\n', d);
        case 10
            fprintf(fid, '%d\n', nSample);
        otherwise
            fprintf(fid, '%s\n', '?');
    end
end
fprintf(fid, '%s\t\n', '#');
fprintf(fid, '%s\t\n', '# Excitation/Emission (Ex, Em), wavelength (nm), component[n] (intensity)');
fprintf(fid, '%s\t\n', '#');

for i=1:size(report,1)
    fprintf(fid, '%s\t', RowHead(i,:));
    fprintf(fid, '%d\t', report(i,1));
    fprintf(fid, '%10.8f\t', report(i,2:end));
    if i<size(report,1);
        fprintf(fid, '\n');
    end
end
fclose(fid);
end