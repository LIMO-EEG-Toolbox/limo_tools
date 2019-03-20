function [cluster, num] = limo_findcluster(onoff, spatdimneighbstructmat, varargin)

% limo_findcluster returns all connected clusters in a 3 dimensional matrix
% with a connectivity of 6.
%
% FORMAT: [cluster, num] = limo_findcluster(onoff, spatdimneighbstructmat)
%         [cluster, num] = limo_findcluster(onoff, spatdimneighbstructmat, minnbchan)
%         [cluster, num] = limo_findcluster(onoff, spatdimneighbstructmat, spatdimneighbselmat, minnbchan)
%
% INPUT: onoff is a 3D boolean matrix with size N1xN2xN3 (typically T/F values > threshold)
%        spatdimneighbstructmat  defines the neighbouring channels/combinations, see below
%        spatdimneighbselmat is a special neighbourhood matrix that is used for selecting
%                            channels/combinations on the basis of the minnbchan criterium
%        minnbchan the minimum number of neighbouring channels/combinations
%
% The neighbourhood structure for the first dimension is specified using
% spatdimneighbstructmat, which is a 2D (N1xN1) matrix. Each row and each column corresponds
% to a channel (combination) along the first dimension and along that row/column, elements
% with "1" define the neighbouring channel(s) (combinations). The first dimension of
% onoff should correspond to the channel(s) (combinations). The lower triangle of 
% spatdimneighbstructmat, including the diagonal, is assumed to be zero.
%
% OUTPUT: cluster is the map of clusters with their respective label
%         num is the number of cluster
%
% This file is taken from FieldTrip, 
% see <http://www.ru.nl/neuroimaging/fieldtrip> for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%    % Copyright (C) 2004, Robert Oostenveld
%
% $Id: findcluster.m 952 2010-04-21 18:29:51Z roboos $
% reused as limo_ft_findluster for integration in LIMO toolbox: GAR, University of Glasgow, June 2010 
% Edit Marianne Latinus adding spm_bwlabel, June 2013
% limo_findcluster is a refactoring of limo_ft_findcluster with the same functionality
% and some optimisations for speed. spm_bwlabel is now used preferentially if availiable
% calls to ismember rearranged and replaced with ismembc
%
% Andrew Stewart, August 2014
% --------------------------------------
% Copyright (C) LIMO Team 2016


spatdimlength = size(onoff, 1);
nfreq = size(onoff, 2);
ntime = size(onoff, 3);

if length(size(spatdimneighbstructmat))~=2 || ~all(size(spatdimneighbstructmat)==spatdimlength)
    error('invalid dimension of spatdimneighbstructmat');
end

minnbchan=2;
if length(varargin)==1
    minnbchan=varargin{1};
end;
if length(varargin)==2
    spatdimneighbselmat=varargin{1};
    minnbchan=varargin{2};
end;

if minnbchan>0
    % For every (time,frequency)-element, it is calculated how many significant
    % neighbours this channel has. If a significant channel has less than minnbchan
    % significant neighbours, then this channel is removed from onoff.
    
    if length(varargin)==1
        selectmat = single(spatdimneighbstructmat | spatdimneighbstructmat');
    end;
    if length(varargin)==2
        selectmat = single(spatdimneighbselmat | spatdimneighbselmat');
    end;
    nremoved=1;
    while nremoved>0
        nsigneighb=reshape(selectmat*reshape(single(onoff),[spatdimlength (nfreq*ntime)]),[spatdimlength nfreq ntime]);
        remove=(onoff.*nsigneighb)<minnbchan;
        nremoved=length(find(remove.*onoff));
        onoff(remove)=0;
    end;
end;

% for each channel (combination), find the connected time-frequency clusters
labelmat = zeros(size(onoff));
total = 0;

% axs - Check for requisites once, rather than within loops
if exist('spm_bwlabel','file') == 3   % axs - preferentially use spm_bwlabel mex if it is in path
    
    for spatdimlev=1:spatdimlength
        [labelmat(spatdimlev, :, :), num] = spm_bwlabel(double(reshape(onoff(spatdimlev, :, :), nfreq, ntime)), 6);
        labelmat(spatdimlev, :, :) = labelmat(spatdimlev, :, :) + (labelmat(spatdimlev, :, :)~=0)*total;
        total = total + num;
    end
    
    
elseif exist('bwlabeln','file') == 2 && exist('spm_bwlabel','file') ~= 3
    
    for spatdimlev=1:spatdimlength
        
        [labelmat(spatdimlev, :, :), num] = bwlabeln(reshape(onoff(spatdimlev, :, :), nfreq, ntime), 4);
        
        labelmat(spatdimlev, :, :) = labelmat(spatdimlev, :, :) + (labelmat(spatdimlev, :, :)~=0)*total;
        total = total + num;
    end
    
    
else
    errordlg('You need either the Image Processing Toolbox or SPM in your path to do clustering');
end




% combine the time and frequency dimension for simplicity
labelmat = reshape(labelmat, spatdimlength, nfreq*ntime);

% combine clusters that are connected in neighbouring channel(s)
% (combinations).
replaceby=1:total;
for spatdimlev=1:spatdimlength
    neighbours=find(spatdimneighbstructmat(spatdimlev,:));
    for nbindx=neighbours
        indx = find((labelmat(spatdimlev,:)~=0) & (labelmat(nbindx,:)~=0));
        for i=1:length(indx)
            a = labelmat(spatdimlev, indx(i));
            b = labelmat(nbindx, indx(i));
            if replaceby(a)==replaceby(b)
                % do nothing
                continue;
            elseif replaceby(a)<replaceby(b)
                % replace all entries with content replaceby(b) by replaceby(a).
                replaceby(find(replaceby==replaceby(b))) = replaceby(a);
            elseif replaceby(b)<replaceby(a)
                % replace all entries with content replaceby(a) by replaceby(b).
                replaceby(find(replaceby==replaceby(a))) = replaceby(b);
            end
        end
    end
end



% renumber all the clusters
num = 0;
cluster = zeros(size(labelmat));

uniquelabel=unique(replaceby(:))';

for j = 1:length(uniquelabel)
    num = num+1;
    
    uniquelabel_here = find(replaceby==uniquelabel(j));
    %labelmat_is_unique_here = ismember(labelmat(:),uniquelabel_here);
    
    labelmat_is_unique_here = ismembc(labelmat(:),uniquelabel_here); % axs - Using the semi-documented ismembc is x3 faster than ismember
    % For ismembc to work, uniquelabel_here MUST be sorted
    % low-high. This should be the case.
    cluster(labelmat_is_unique_here) = num;
end


% reshape the output to the original format of the data
cluster = reshape(cluster, spatdimlength, nfreq, ntime);

