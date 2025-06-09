function [r,c]=relcomporder(data,RefSplit)
%
% <strong>Syntax</strong>
%   [r,c]=<strong>relcomporder</strong>(data,RefSplit)
%
% <a href="matlab: doc relcomporder">help for relcomporder</a> <- click on the link

%determine the order of components for plotting overlaid spectra relative
%to a defined 'reference' split.
%
% USEAGE 
%        [r,c]=relcomporder(data,RefSplitNum)
%
% INPUTS
%        data: validation data structure obtained from SplitValidation.m
% RefSplitNum: Number of the reference split, to be plotted 'as it comes'
%              and which all other split models will be plotted in relation to.
%
% OUTPUTS
%     r:  list of matching components [nSplits x nComponents] relative to 
%         RefSplit. Row i shows the components in Split i that matched the 
%         reference split components 1:nComponents. NaNs indicate
%         components in the reference split that did not find a match.
%     c:  [nSplits x 2] matrix where the first column is the i'th split
%         and the second column is the reference split. 
%
% EXAMPLE
%         In the case that Val_Splits= {[1 3],[2 4],[1 4],[2 3]}
%         then:
%         RefSplit=1 produces the order for plotting against [1 3]
%         RefSplit=2 produces the order for plotting against [2 4]
%         RefSplit=3 produces the order for plotting against [1 4]
%         RefSplit=4 produces the order for plotting against [2 3]
%
%       Note:
%         since in this example Val_Splits does not include all possible
%         combinations (missing: [1 2],[3 4]) it is impossible to calculate
%         the order for plotting all components against each other, e.g.
%         split [1 2] vs other splits. To get orders for all possible
%         combinations, first run SplitValidation choosing the default
%         selection that compares all possible combinations of splits.
%         e.g.
%         compareall=SplitValidation(data,6); %to compare all the 6-comp models in data.Split 
%         [r,c]=relcomporder(compareall,4); %to match their components to data.Split(4).Model6
%
% Notice:
% This mfile is part of the drEEM toolbox. Please cite the toolbox
% as follows:
%
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 2013, 
%     DOI:10.1039/c3ay41160e. 
%
% relcomporder: Copyright (C) 2013 Kathleen R. Murphy
% The University of New South Wales
% Dept Civil and Environmental Engineering
% Water Research Center
% UNSW 2052
% Sydney
% krm@unsw.edu.au
%
% $ Version 0.1.0 $ September 2013 $ First Release

narginchk(2,2)
if ~isfield(data,'Val_Matches')
    error('Input data matrix does not contain split validation field: Val_Matches')
elseif ~isfield(data,'Val_Comparisons')
    error('Input data matrix does not contain split validation field: Val_Comparisons')
end

if isfield(data,'Val_Splits')
    m=cell2mat(data.Val_Matches(1:end-1));
else
    m=cell2mat(data.Val_Matches(1:end));
end
comparisons=data.Val_Comparisons_Num;
D=size(m);

insplits=unique(comparisons);
r=NaN*ones(length(insplits),D(2));
c=NaN*ones(length(insplits),2);
r(RefSplit,:)=1:D(2);
c(RefSplit,:)=[RefSplit RefSplit];

%forward matches
[i,~]=find(comparisons(:,2)==RefSplit);
rowsf=comparisons(i,:);
if ~isempty(rowsf)
    r(rowsf(:,1),:)=m(i,:);
    c(rowsf(:,1),:)=comparisons(i,:);
end

%backward matches
[j,~]=find(comparisons(:,1)==RefSplit);
rowsb=comparisons(j,:);
if ~isempty(rowsb)
    r(rowsb(:,2),:)=reverseorder(m(j,:));
    c(rowsb(:,2),:)=fliplr(comparisons(j,:));
end
end

function y=reverseorder(ord)
y=NaN*ones(size(ord));
D=size(y);
for j=1:D(1)
    for i=1:D(2)
        t=find(ord(j,:)==i);
        if ~isempty(t)
            y(j,i)=t;
        end
    end
end
end
