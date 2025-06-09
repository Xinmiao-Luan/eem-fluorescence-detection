function peaks=describecomp(data,f,varargin)
%
% <strong>Syntax</strong>
%   peaks=<strong>describecomp</strong>(data,f,varargin)
%
% <a href="matlab: doc describecomp">help for describecomp</a> <- click on the link

% Describe PARAFAC components in terms of local maxima. A list of peaks
% will be provided for each component's Ex and Em spectrum. Note that if 
% the spectrum does not include an inflection point (e.g. for a peak
% occurring below the range of the data) the peak will not appear in this list.
% Therefore it is always necessary to check the results of this function
% against the actual spectra.
% 
% Useage:   peaks=describecomp(data,f,run)
%
% Inputs: 
%      data: dataset structure containing PARAFAC model results
%         f: Number of components in the model to be plotted, 
%            e.g. 6 to plot the 6-component model in data.Model6.
%       run: (optional) 
%            Run number of model to be described, when data contains
%            multiple runs of the same model, as from the output
%            of randinitanal.
%             []: the main model will be described (default)
%              n: plot run number x (the model plotted will be
%                data.Modeln_runx, e.g. data.Model6_run2)
%
% Examples
%    describecomp(LSmodel6,6)
%    peaks=describecomp(DSit6,6,5)
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
%Copyright describecomp (C) 2013- KR Murphy
%Copyright findpeaksTH (C) 2009- Tom O'Haver 
%
% $ Version 0.1.0 $ September 2013 $ First Release

%%
%Initialise
narginchk(2,3)
R=[];
if length(f)>1
    error('Specify one value of ''f'' at a time');
else
    modelf=['Model' num2str(f)];
end
if nargin>2
    R=varargin{1};
end
if ~isempty(R)
    if ~isnumeric(R)
        error('Input for variable ''run'' is not understood.')
    else
        modelf=[modelf '_run' int2str(R)];
    end
end

if ~isfield(data,modelf)
   disp(modelf)
   disp(data)
   error('loadingsandleverages:fields',...
       'The dataset does not contain a model with the specified number of factors') 
end
M = getfield(data,{1,1},modelf);
Ex=data.Ex;
Em=data.Em;

v=ver;fp=0;
if any(strcmp('Signal Processing Toolbox', {v.Name}))
    fp=1;
else
    fprintf('Signal Processing Toolbox not located. Using open-source findpeaks from T.C.O''Haver\n')
end

peaks=cell(2,f);
for col=1:2
    switch col
        case 1
            w=Em;
        case 2
            w=Ex;
    end
    for i=1:f
        if fp==1
            [pks,locs]=findpeaks(M{col+1}(:,i));
            peaks(col,i)={[w(locs) (pks)]};
        else
            P=findpeaksTH(w,M{col+1}(:,i),0,-1,2,3);
            peaks(col,i)={P(:,2:3)};
        end
    end
end

disp(' ')
disp(' ')
disp('************************************')
disp(['Description of PARAFAC components - ' modelf])
for i=1:f
disp('************************************')
disp(['Component ' int2str(i)])
disp('Ex ')
t2=cell2mat(peaks(2,i));
disp(t2)
disp('Em ' )
t1=cell2mat(peaks(1,i));
disp(t1)
end

function P=findpeaksTH(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% function
% P=findpeaks(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype) 
% Function to locate the positive peaks in a noisy x-y time series
% data set.  Detects peaks by looking for downward zero-crossings in the
% first derivative that exceed SlopeThreshold. Returns list (P) containing
% peak number and position, height, width, and area of each peak, assuming
% a Gaussian peak shape. Arguments "slopeThreshold", "ampThreshold" and
% "smoothwidth" control peak sensitivity. Higher values will neglect
% smaller features. "Smoothwidth" is the width of the smooth applied before
% peak detection; larger values ignore narrow peaks. If smoothwidth=0, no
% smoothing is performed. "Peakgroup" is the number points around the top 
% part of the peak that are taken for measurement. If Peakgroup=0 the local
% maximum is takes as the peak height and position. The argument
% "smoothtype" determines the smooth algorithm:
%   If smoothtype=1, rectangular (sliding-average or boxcar) 
%   If smoothtype=2, triangular (2 passes of sliding-average)
%   If smoothtype=3, pseudo-Gaussian (3 passes of sliding-average)
%
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and 
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
%
% T. C. O'Haver, 1995.  Version 5.1, Last revised December, 2012
% Skip peaks if peak measurement results in NaN values
%
% Examples:
% findpeaks(0:.01:2,humps(0:.01:2),0,-1,5,5)
% x=[0:.01:50];findpeaks(x,cos(x),0,-1,5,5)
% x=[0:.01:5]';findpeaks(x,x.*sin(x.^2).^2,0,-1,5,5)
%
%
% Copyright (c) 2009, Tom O'Haver
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

if nargin~=7;smoothtype=1;end  % smoothtype=1 if not specified in argument
if smoothtype>3;smoothtype=3;end
if smoothtype<1;smoothtype=1;end 
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
if smoothwidth>1,
    d=fastsmooth(deriv(y),smoothwidth,smoothtype);
else
    d=y;
end
n=round(peakgroup/2+1);
P=[0 0 0 0 0];
vectorlength=length(y);
peak=1;
AmpTest=AmpThreshold;
for j=2*round(smoothwidth/2)-1:length(y)-smoothwidth,
    if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold*y(j), % if slope of derivative is larger than SlopeThreshold
            if y(j) > AmpTest,  % if height of peak is larger than AmpThreshold
                xx=zeros(size(peakgroup));yy=zeros(size(peakgroup));
                for k=1:peakgroup, % Create sub-group of points near peak
                    groupindex=j+k-n+2;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    xx(k)=x(groupindex);yy(k)=y(groupindex);
                end
                if peakgroup>3,
                    [coef,S,MU]=polyfit(xx,log(abs(yy)),2);  %#ok<ASGLU> % Fit parabola to log10 of sub-group with centering and scaling
                    c1=coef(3);c2=coef(2);c3=coef(1);
                    PeakX=-((MU(2).*c2/(2*c3))-MU(1));   % Compute peak position and height of fitted parabola
                    PeakY=exp(c1-c3*(c2/(2*c3))^2);
                    MeasuredWidth=norm(MU(2).*2.35482/(sqrt(2)*sqrt(-1*c3)));
                    % if the peak is too narrow for least-squares technique
                    % to work well, just use the max value of y in the
                    % sub-group of points near peak.
                else
                    PeakY=max(yy);
                    pindex=val2ind(yy,PeakY);
                    PeakX=xx(pindex(1));
                    MeasuredWidth=0;
                end
                % Construct matrix P. One row for each peak detected,
                % containing the peak number, peak position (x-value) and
                % peak height (y-value). If peak measurements fails and
                % results in NaN, skip this peak
                if isnan(PeakX) || isnan(PeakY) || PeakY<AmpThreshold,
                    % Skip any peak that gives a NaN or is less that the
                    % AmpThreshold
                else % Otherwiase count this as a valid peak
                    P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth  1.0646.*PeakY*MeasuredWidth];
                    peak=peak+1; % Move on to next peak
                end % if isnan(PeakX)...
            end % if y(j) > AmpTest...
        end  % if d(j)-d(j+1) > SlopeThreshold...
    end  % if sign(d(j)) > sign (d(j+1))...
end  % for j=....
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);

function d=deriv(a)
% First derivative of vector using 2-point central difference.
%  T. C. O'Haver, 1988.
n=length(a);
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
  d(j)=(a(j+1)-a(j-1)) ./ 2; %#ok<AGROW>
end

function SmoothY=fastsmooth(Y,w,type,ends)
% fastbsmooth(Y,w,type,ends) smooths vector Y with smooth 
%  of width w. Version 2.0, May 2008.
% The argument "type" determines the smooth type:
%   If type=1, rectangular (sliding-average or boxcar) 
%   If type=2, triangular (2 passes of sliding-average)
%   If type=3, pseudo-Gaussian (3 passes of sliding-average)
% The argument "ends" controls how the "ends" of the signal 
% (the first w/2 points and the last w/2 points) are handled.
%   If ends=0, the ends are zero.  (In this mode the elapsed 
%     time is independent of the smooth width). The fastest.
%   If ends=1, the ends are smoothed with progressively 
%     smaller smooths the closer to the end. (In this mode the  
%     elapsed time increases with increasing smooth widths).
% fastsmooth(Y,w,type) smooths with ends=0.
% fastsmooth(Y,w) smooths with type=1 and ends=0.
% Example:
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]
%  T. C. O'Haver, May, 2008.
if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
  switch type
    case 1
       SmoothY=sa(Y,w,ends);
    case 2   
       SmoothY=sa(sa(Y,w,ends),w,ends);
    case 3
       SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends);
  end

function SmoothY=sa(Y,smoothwidth,ends)
w=round(smoothwidth);
SumPoints=sum(Y(1:w));
s=zeros(size(Y));
halfw=round(w/2);
L=length(Y);
for k=1:L-w,
   s(k+halfw-1)=SumPoints;
   SumPoints=SumPoints-Y(k);
   SumPoints=SumPoints+Y(k+w);
end
s(k+halfw)=sum(Y(L-w+1:L));
SmoothY=s./w;
% Taper the ends of the signal if ends=1.
  if ends==1,
    startpoint=(smoothwidth + 1)/2;
    SmoothY(1)=(Y(1)+Y(2))./2;
    for k=2:startpoint,
       SmoothY(k)=mean(Y(1:(2*k-1)));
       SmoothY(L-k+1)=mean(Y(L-2*k+2:L));
    end
    SmoothY(L)=(Y(L)+Y(L-1))./2;
  end
% ----------------------------------------------------------------------


