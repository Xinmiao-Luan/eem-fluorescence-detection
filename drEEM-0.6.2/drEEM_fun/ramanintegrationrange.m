function [varargout] = ramanintegrationrange(RAMmat,flist,landa,varargin)
%
% <strong>Syntax</strong>
%
%   [varargout] = <strong>ramanintegrationrange</strong>(RAMmat,flist,landa,varargin)
%   [IR,IRmed,IRdiff]=<strong>ramanintegrationrange</strong>(RAMmat,flist,landa,halfwidth,sequence,tolerance,plotview)
%
% <a href="matlab: doc ramanintegrationrange">help for ramanintegrationrange</a> <- click on the link

% function [IR IRmed IRdiff] = ramanintegrationrange(RAMmat,flist,landa,halfwidth,sequence,tolerance,plotview);
% Determine the integration range [min max] for the water Raman peak at 
% specified excitation wavelengths for a set of water Raman scans
%
%INPUTS
% 		RAMmat:	matrix of Raman scans, with wavelength in 1st row and
%        	 	scans for n samples in n subsequent rows.
% 		flist:  list of filenames (text) or a list of sample numbers, e.g. 1:n
% 		landa:  excitation wavelength; either a single value for all scans,
%         		or a list of 1:n wavelengths matching each scan.
% 	halfwidth:	(optional input) maximum width (cm-1) of the scatter peak. If not
%           	specified, the default value of 1800 cm-1 will be used.
% 	 sequence: 	number of consecutive points identifying a sequence (default=8)
% 	tolerance: 	minimum threshold identifying a non-zero slope. Default = 1% of maximum slope 
% 	 plotview: 	length of time to show plots of scans on screen (e.g. 1 = display for 1 second)
%             	if display=0, plots will be created, displayed 0.5s, but not saved. 
%             	if display=[], the default is used, where the default is
%               to wait for user input before showing the next plot, then save plots to file
%
%OUTPUTS
% 		   IR:	integration range [min max] for the Raman peak at specified landa 
% 	    IRmed:	if a single landa was specified, IRmed is the median integration range across the dataset
% 	   IRdiff:	if a single landa was specified, IRdiff is the deviation of each sample
%             	from the median integration range.
%PLOTS
% 	Plots are automatically appended to a postscript file 
% 	and saved in the current directory to the file "RamanIntegrationRange.ps" 
% 	View saved plots using a postscript reader, e.g. Ghostview (freeware), Adobe Acrobat Pro, etc
%
%EXAMPLES
% [IR IRmed IRdiff] = ramanintegrationrange(RAMmat,1:10,350); %landa=350, default settings
% [IR] = ramanintegrationrange(RAMmat,1:3,[300 320 350]); %landa varies, default settings
% [IR] = ramanintegrationrange(RAMmat,['a','b'],[300 320],[],6); %landa varies, sequence =6
% [IR IRmed IRdiff] = ramanintegrationrange(RAMmat,['a','b','c'],300,[],[],0.03); %landa=300, tolerance=3%
% [IR IRmed IRdiff] = ramanintegrationrange(RAMmat,['a','b','c'],300,[],[],[],1); %display plots for 1 sec
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% ramanintegrationrange: Copyright (C) 2013 Kathleen R. Murphy
% $ Version 0.1.0 $ September 2013 $ First Release
% $ Updated from FDOMcorr toolbox ver. 1.8
%
% Copyright (C) 2013 KR Murphy,
% Water Environment Technology
% Chalmers University of Technology
% Department of Civil and Environmental Engineering
% Water Environment Technology
% Sven Hultins gata 8
% 412 96 Göteborg
% murphyk@chalmers.se
%%%%%%%%%%%%%%%%

%Check and assign input variables
%Assign default values to missing variables
narginchk(3,7)
halfwidth=1800;
sequence=8;
tolerance=0.01;
plotview=0.5;
if nargin>3;
    halfwidth=varargin{1};
    if isempty(halfwidth)
        halfwidth=1800;
    end
    if nargin>4;
        sequence=varargin{2};
        if isempty(sequence)
            sequence=8;
        end
        if nargin>5;
            tolerance=varargin{3};
            if isempty(tolerance)
                tolerance=0.01;
            end
        end
        if nargin>6;
            plotview=varargin{4};
        end
    end
end

if ispos(gradient(RAMmat(1,:)))==0
    error(['check that wavelengths are in the first row of RAMmat, and '...
    'that scans for different samples are in subsequent rows!'])
end

nargoutchk(1,3)
if nargout==3
    if max(size(landa))>1;
        error('specify one output argument only when Ex varies between scans')
    end
end


%Generate plots of the results and save them to the current directory
%Plots are appended to the existing file named RamanIntegrationRange.ps
%[NoRows d_1]=size(RAMmat);
NoRows=size(RAMmat,1);
NoSamples=NoRows-1;
IR=zeros(NoSamples,2);
delete RamanIntegrationRange.ps
figure
for i=2:NoRows
    if max(size(landa))==1;
        landa_ex=landa;
    else
        landa_ex=landa(i-1);
    end
    R=[RAMmat(1,:);RAMmat(i,:)];
    fn=flist(i-1,:);
    if isnumeric(fn)
        fn=num2str(fn);
    end
    disp(fn)
    IR(i-1,:)=RamanIR(R,fn,landa_ex,halfwidth,sequence,tolerance);
    orient landscape,
    if plotview==0
        pause(0.5)
        %close
    else
        print -dpsc -append RamanIntegrationRange  %save (append) figures to file
        if isempty(plotview)
            pause
            fprintf('Plot written to file. Press any key to continue \n\n\n')
            %close
        else
            pause(plotview)
            fprintf('Plot written to file.  \n\n\n')
            %close
        end
    end
end
if ~plotview==0
    wh=what;
    fprintf('plots saved to the current directory:  \n');
    wh.path
    pause(2)
end
IRmed=round(median(IR,1)*2)/2;
IRdif=IR-ones(NoSamples,1)*IRmed;
%fprintf('results=[(1:NoSamples) IR_sample IR_median difference]')
%[(1:NoSamples)' IR ones(NoSamples,1)*IRmed IR-ones(NoSamples,1)*IRmed]

%ASSIGN OUTPUT VARIABLES
varargout{1}=IR;
if max(size(landa))==1;
    varargout{2}=IRmed;
    varargout{3}=IRdif;
end

function IR=RamanIR(R,fn,landa_ex,halfwidth,sequence,tolerance)
% function IR=RamanIR(R,fn,landa_ex,halfwidth,sequence,tolerance)
% Calculates the integration range IR for a the raman peak in a scan
% based upon the gradient of the points in the scan. IR is bounded by the
% halfwidth. The peak is assumed to begin at the first defined-length
% sequence of positive slopes above a threshold value (defined by tolerance)
% and end at the last sequence of negative slopes above threshold.
%
% INPUTS
% R     =       Raman scan, first row wavelength, second row fluorescence
% fn    =       filename or sample number
% landa_ex =    excitation wavelength for the Raman scan
% halfwidth =   maximum spread of Raman scatter from the peak (cm-1)
% sequence =    number of consecutive points identifying a sequence 
% tolerance =   minimum threshold identifying a non-zero slope
%
% OUTPUTS
% IR =      integration range [min max] for the Raman peak at specified landa_ex 
% Copyright (C) 2010 KR Murphy,
% The University of New South Wales
% Department of Civil and Environmental Engineering
% Sydney 2052, Australia
% krm@unsw.edu.au

if size(R,1)<size(R,2)
    R=R';
end

[landa_em,Emmin,Emmax]=RamanPeakPos(landa_ex,halfwidth);
landa_em=round(landa_em*10)/10;

%Determine gradients
x=R(:,1);
y=R(:,2);
inc=round((max(x)-min(x))/(length(x)-1)*2)/2;
[~,fy] = gradient([x y],inc);

g=fy(:,2);
sg=smooth(g); %5 point moving average
%[x y sg];

%Identify sequences of positive and negative gradients
p=zeros(length(sg),1);n=p;
sq=sequence-1;
for i=1:length(sg)-sq
    sgsub_pos=sg(i:i+sq);
    p(i)=ispos(sgsub_pos);
end
for i=1+sq:length(sg)
    sgsub_neg=sg(i-sq:i);
    n(i)=isneg(sgsub_neg);
end
sgp=sg.*p;
sgn=sg.*n;
%[x sg n  p sgn sgp], pause
sgmax=max(sg(x>Emmin)); %avoid rayleigh scatter in determining the maximum gradient
xi=x(sgp-tolerance*sgmax>0);sgpi=sgp(sgp-tolerance*sgmax>0);
xf=x(sgn+tolerance*sgmax<0);sgnf=sgn(sgn+tolerance*sgmax<0);
gradpos=[xi(xi>Emmin) sgpi(xi>Emmin)]; %positive sequences within integration range
gradneg=[xf(xf<Emmax) sgnf(xf<Emmax)]; %negative sequences within integration range
if isempty(gradpos)
    %Emmin,Emmax,[x sg sgp],[xi sgpi],
    %figure,plot(x,sg),hold on, plot(x,sgp,'ro')
    fprintf('\n no positive sequences of desired length within the integration range')
    fprintf('\n check that excitation wavelengths are correctly specified')
    error('to trouble shoot this error, try reducing sequence length and setting tolerance = 0')
end
if isempty(gradneg)
    %Emmin,Emmax,[x sg sgn],[xi sgni],
    fprintf('\n')
    fprintf('\n no negative sequences of desired length within the integration range')
    fprintf('\n check that excitation wavelengths are correctly specified')
    error('to trouble shoot this error, try reducing sequence length and setting tolerance = 0')
end
fprintf('\n')
IR=[floor(gradpos(1,1)*2)/2 ceil(gradneg(end,1)*2)/2]; %integration range

%Figure with three subplots summarizing the results
%figure
subplot(1,3,1)
plot(x,y,'bo'),
hold on, 
plot(x,g,'cx')
plot(x,sg,'k-')
axis('tight');
legend('Raman','gradient','smooth gradient')
ylabel('Signal height (arbitrary units)')
hold off

subplot(1,3,2)
plot(x,sg,'k-');
hold on
plot(xi,sgpi,'go');
plot(xf,sgnf,'co');
plot(gradpos(1,1),gradpos(1,2),'r*');
plot(gradneg(end,1),gradneg(end,2),'r*');
axis('tight');
ylim([-.25*sgmax .25*sgmax])
legend('smooth gradient','gradient +','gradient -','integration boundary')
plot([min(x);max(x)],[0;0],'r--');
plot([Emmin;Emmin],[min(sg)/4;max(sg)/4],'r--');
plot([Emmax;Emmax],[min(sg)/4;max(sg)/4],'r--');
title(fn)
xlabel('Emission wavelength (nm)')
hold off

subplot(1,3,3)
xsub=x(intersect(find(x>IR(1)),find(x<IR(2))));
ysub=y(intersect(find(x>IR(1)),find(x<IR(2))));
plot(x,y,'bo');
hold on,
bar(xsub,ysub);
axis('tight');
ylim([0 1.02*max(ysub)])
text(landa_em+10,2.95/3*max(ysub),['Ex = ' num2str(landa_ex) 'nm'])
text(landa_em+10,2.85/3*max(ysub),['Em_p_e_a_k = ' num2str(landa_em) 'nm'])
text(landa_em+10,2.75/3*max(ysub),['IR_m_i_n = ' num2str(IR(1)) 'nm'])
text(landa_em+10,2.65/3*max(ysub),['IR_m_a_x = ' num2str(IR(2)) 'nm'])
text(landa_em+10,2.55/3*max(ysub),['peak width = ' num2str(IR(2)-IR(1)) 'nm'])
text(landa_em+10,2.45/3*max(ysub),['tol = ' num2str(tolerance)])
text(landa_em+10,2.35/3*max(ysub),['sequence = ' num2str(sequence)])
hold off

function [landa_em,Emmin,Emmax]=RamanPeakPos(landa_ex,halfwidth)
% function [landa_em,Emmin,Emmax]=RamanPeakPos(landa_ex,halfwidth)
% For a given excitation wavelength (landa_ex), calculate the water Raman
% peak position (landa_em) and an integration range (Emmin-Emmax) for 
% calculating Raman area.
% The water Raman peak appears at incident wavenumber minus 3380 cm-1
%
% INPUTS:
% landa_ex = excitation wavelength of the water Raman peak
% halfwidth=wave number distance (cm-1) from Raman peak maximum to its upper/lower extent
%
% OUTPUTS:
% landa_em = emission wavelength of the water Raman peak maximum
% Emmin = lower emission extent of the peak
% Emmax = upper emission extent of the peak
%
% Copyright (C) 2010 KR Murphy,
% The University of New South Wales
% Department of Civil and Environmental Engineering
% Sydney 2052, Australia
% krm@unsw.edu.au

landa_em=1/(1/landa_ex - 3380/10^7);         %nm  
Emmin=1/(1/landa_em + halfwidth*10^-7);      %nm   
Emmax=1/(1/landa_em - halfwidth*10^-7);      %nm   
Emmin=floor(Emmin);     %nm   
Emmax=ceil(Emmax);      %nm   

function TrueOrFalse=ispos(Xmatrix)
%  Returns true if Xmatrix contains only positive elements, 
%  and false otherwise. Ignore NaNs
if any(Xmatrix(isfinite(Xmatrix)) <= 0) 
	TrueOrFalse=0;
else
	TrueOrFalse=1;
end

function TrueOrFalse=isneg(Xmatrix)
%  Returns true if Xmatrix contains only negative elements
%  and false otherwise. Ignore NaNs
if any(Xmatrix(isfinite(Xmatrix)) >= 0) 
	TrueOrFalse=0;
else
	TrueOrFalse=1;
end

function z=smooth(x)
%5 point moving average
z=NaN*ones(size(x));
for i= 1:length(x),
    if i==1;
        z(i)=x(i);
    elseif i==2;
        z(i)=sum(x(i-1:i+1))/3;         
    elseif i>=3;
        if i<(length(x)-2)
            z(i)=sum(x(i-2:i+2))/5;
        elseif  i==(length(x)-1)
            z(i)=sum(x(i-1:i+1))/3;  
        elseif  i==length(x)
            z(i)=x(length(x));
        end
    end
end
