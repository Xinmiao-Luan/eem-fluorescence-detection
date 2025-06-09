function [IFC,K]=ABAife(Ex,Em,A,X)
%
% <strong>Syntax</strong>
%   [IFC,K]=<strong>ABAife</strong>(Ex,Em,A,X)
%
% <a href="matlab: doc ABAife">help for ABAife</a> <- click on the link

% USEAGE:
%  [IFC,K]=ABAife(Ex,Em,A,X)
%
% Absorbance (ABA) Method for correcting fluorescence inner filter effects.
% Assumes that EEMs and Absorbance were measured in a square cuvette with
% standard right-angle geometry.
% All samples must have the same wavelengths.
%
% Input Variables:
%   Ex      excitation wavelengths
%   Em      emission wavelengths
%   A       absorbance scans for n samples over the full range of Ex and Em, 
%               with wavelengths in 1st row, samples in rows 2:n+1.
% Optional
%   X       fluorescence EEM (n x Em rows x Ex columns) measured in a 1 cm cell
%
% Output Variables:
%   IFC     inner filter matrix (n x Em rows x Ex columns)
%   K       corrected EEM (n x Em rows x Ex columns), only if X was given
%
% Examples          
%   IFC=ABAife(Ex,Em,A)
%   [IFC,K]=ABAife(Ex,Em,A,X)
%
%Reference for ABA algorithm:
% Parker, C. A.; Barnes, W. J., Some experiments with spectrofluorimeters
% and filter fluorimeters. Analyst 1957, 82, 606-618.
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 
%     5, 6557-6566, 2013. DOI:10.1039/c3ay41160e. 

narginchk(3,4) %check input arguments
K=[];
N_samples=size(A,1)-1;
N_em=length(Em);
N_ex=length(Ex);
IFC=NaN*ones([N_samples,N_em,N_ex]);

for i=1:N_samples;
    Atot=zeros(N_em,N_ex);
    Ak=[A(1,:);A(i+1,:)];
    
    Aem=MatchWave(Ak,Em,0);
    Aex=MatchWave(Ak,Ex,0);
    
    for j=1:size(Aem,1)
        for k=1:size(Aex,1)
            Atot(j,k)=Aex(k,2)+Aem(j,2);
        end
    end
    ifc=10.^(1/2*Atot);
    IFC(i,:,:)=ifc;
end 
if nargin==4
K=IFC.*X;
end  
end

function Y=MatchWave(V1,V2,doplot)
%Y=MatchWave(V1,V2,doplot)
%Automatically match the size of two vectors, interpolating if necessary.
%Vector V1 is resized and interpolated to have the same wavelengths as V2;
%For example, V1 is a correction file (0.5nm intervals) and V2 is an Emission scan in 2 nm intervals.
%Errors are produced if 
%(1) there are wavelengths in V2 that are outside the upper or lower limit of V1.
%(2) Be aware of rounding errors. For example, 200.063 is NOT equivalent to 200 nm.
% Errors can often be resolved by restricting the wavelength range of V2
% Copyright K.R. Murphy
% July 2010

if size(V1,2)>2; V1=V1';end  %V1 should have multiple rows
if size(V2,2)>2; V2=V2';end  %V2 should have multiple rows

%Restrict the wavelength range of V1 so it is the same as V2
t=V1(find(V1(:,1)<=V2(1,1),1,'last'):find(V1(:,1)>=V2(end,1),1,'first'),:); 
if isempty(t)
    fprintf('\n')
    fprintf(['Error - Check that your VECTOR 1 (e.g. correction file) '...
        'encompasses \n the full range of wavelengths in VECTOR 2 (e.g. emission scan)\n'])
    fprintf('\n Hit any key to continue...')
    pause
    fprintf(['\n VECTOR 1 range is ' num2str(V1(1,1)) ' to ' num2str(V1(end,1)) ' in increments of ' num2str(V1(2,1)-V1(1,1)) '.']),pause
    fprintf(['\n VECTOR 2 range is ' num2str(V2(1,1)) ' to ' num2str(V2(end,1)) ' in increments of ' num2str(V2(2,1)-V2(1,1)) '.']),pause
    fprintf('\n\n This error can usually be resolved by restricting the wavelength range of VECTOR 2,');
    fprintf('\n also, be aware of rounding errors (e.g. 200.063 is NOT equivalent to 200 nm)\n');
    error('fdomcorrect:wavelength','Error - abandoning calculations.');
else
    Y=[V2(:,1) interp1(t(:,1),t(:,2),V2(:,1))]; %V1 corresponding with EEM wavelengths
end

if doplot==true;    
figure, 
plot(V1(:,1),V1(:,2)), hold on, plot(Y(:,1),Y(:,2),'ro')
legend('VECTOR 1', 'VECTOR 2')
end
end