function [L,H,S] = CDAife(Xu,Xd,p,e)
%
% <strong>Syntax</strong>
%   [L,H,S]=<strong>CDAife</strong>(Xu,Xd,p,e)
%
% <a href="matlab: doc CDAife">help for CDAife</a> <- click on the link

% Controlled Dilution (CDA) approach to IFE correction
%
% Input Variables:
%   Xu An undiluted fluorescence EEM (Em rows x Ex wavelengths)
%   Xd A diluted EEM corresponding to Xu (Em rows x Ex wavelengths)
%   p = dilution factor (undiluted/diluted, a scalar value with p>1)
% Optional input variables
%   e = assumed experimental error, use [] for the default e=0.05
%
% Output Variables:
%   L corrected sample EEM
%   H inner filter matrix
%   S = 1st order approximation of sensitivity
%
% Notes: Xu and Xd must be of identical size with identical excitation and emission wavelengths
%
% Reference for CDA algorithm: 
%    X. Luciani et al. / Chemometrics and Intelligent Laboratory Systems 96 (2009) 227–238
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 
%     5, 6557-6566, 2013. DOI:10.1039/c3ay41160e. 

narginchk(3,4) %check input arguments
if ~isscalar(p)
    error('drEEM:CDAife','The dilution factor must be a scalar value')
else
    if p<1
        warning('The dilution factor is less than one! Hit any key to continue...')
        pause
    end
end
if nargin<4;
    e=0.05;
end
if and(ismatrix(Xu),ismatrix(Xd))
    if isequal(size(Xu),size(Xd))
        L= (((p*Xd).^p)./Xu).^(1/(p-1));     %equ. 21
        H= (Xu./(p*Xd)).^(p/(p-1));          %equ. 22
        S=((p-1-log(p)-log(Xd)+log(Xu)))*e/(p-1)^2;
    else
        fprintf('\n Input matrices Xu and Xd are not the same size \n')
        dim1=size(Xu);dim2=size(Xd);
        fprintf(['\n Xu: ' num2str(dim1(1)) ' by ' num2str(dim1(2)) '\n'])
        fprintf(['\n Xd: ' num2str(dim2(1)) ' by ' num2str(dim2(2)) '\n'])
        error('Check input data')
    end
else
    error('Xu and Xd must be matrices of size m x n')
end
