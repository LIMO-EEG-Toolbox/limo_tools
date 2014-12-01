function IDX = limo_clusterica(varargin)

% FORMAT IDX = eeg_clusterica(varargin)
%
% INPUT varargin is an array 'key', data, 'key;, data ...
%       the 'key' indicates which type of data are passed
%       these can be: 'erp' and data are of dim [ic time]
%                     'spect' and data are of dim [ic, freqs]
%                     'ersp' and data are of dim [ic, freqs, time]
%                     'itc' and data are of dim [ic, freqs, time]
%                     'smap' a scalp map of dim 64*64]
%                     'dipole' and data are of dim [x,y,z] for spatial lovcation
%                              or [x,y,z,xx,yy,zz] adding the orientation vector
%                              assuming a common origin (i.e. not [x y z])
%
% OUPUT IDX is the indices of the compoments after clustering
%
% The algorithm compute several similarity matrices based on (i) slightly
% lagged cross-correlations of erp, spec, ersp, itc (ii) correlation of scalp
% map (iii) euclidian distance of dipoles and (iv) angles of dipole orientation
% Once those measures are obtained, the clustering is performed using
% affinity propagation on the average of erp, ersp, spec, itc and the average
% of dipole location and angles - the intersectrion gives the final
% clusters. Note data are paased with all IC from all subjects - that is
% the clustering is blind of subject belonging and only relies on
% similarity measurements.
%
% example 
% IDX = limo_clusterica('erp',randn(75,200),'spec',randn(75,200), ...
%     'ersp',randn(75,125,125),'itc',randn(75,125,125))
%                 
%
% Cyril Pernet, The University of Edinburgh
% Arnaud Delorme, SCCN, INC, UCSD
% Ramon Martinez-Cancino,INC, SCCN
%% -----------------------------
% Copyright (C) LIMO Team 2014

% default lag
lag = 2;
index = 1;

% note one computes xcorr on transposed data, i.e. between ic
% -----------------------------------------------------------
for m=1:2:length(varargin)
    
    if strcmpi(varargin{m},'erp')
        % the absolute of xcorr to be insensitive to polarity
        disp('clustering components using ERPs')
        M{index} = local_xcorr(cell2mat(varargin(m+1))',lag,'abs');
        C{index} = apcluster(M{index},median(M{index}));
        varargin(m+1) = {[]}; % free memory as we go along
        index = index+1;
    end
    
    if strcmpi(varargin(m),'spec')
        % power is already positive, take max of xcorr
        disp('clustering components using Spectra')
        M{index} = local_xcorr(cell2mat(varargin(m+1)'),lag,'signed');
        C{index} = apcluster(M{index},median(M{index}));
        varargin(m+1) = {[]}; index = index+1;
    end
    
    if strcmpi(varargin(m),'ersp')
        % power is already positive, take max of xcorr
        % note the matrix M is the average of the best 5 frequency bands
        % no point looking at frequency bands with noise
        disp('clustering components using ERSP')
        data = cell2mat(varargin(m+1));
        for f = 1:size(data,2)
            S(:,:,f) = local_xcorr(squeeze(data(:,f,:))',lag,'signed');
            c = triu(squeeze(S(:,:,f)),1); score(f) = mean(c(:));            
        end
        [~,ranking]=sort(score);
        M{index} = mean(S(:,:,ranking(1:5)),3); clear S
        C{index} = apcluster(M{index},median(M{index}));
        varargin(m+1) = {[]}; index = index+1;
    end
    
    if strcmpi(varargin(m),'itc')
        % phase can be opposed for opposed polarities, so we take the abs
        % note the matrix M is the average of the best 5 frequency bands
        % no point looking at frequency bands with noise
        disp('clustering components using ITC')
        data = cell2mat(varargin(m+1));
        for f = 1:size(data,2)
            S(:,:,f) = local_xcorr(squeeze(data(:,f,:))',lag,'abs');
            c = triu(squeeze(S(:,:,f)),1); score(f) = mean(c(:));            
        end
        [~,ranking]=sort(score);
        M{index} = mean(S(:,:,ranking(1:5)),3); clear S
        C{index} = apcluster(M{index},median(M{index}));
        varargin(m+1) = {[]}; index = index+1; 
    end
    
    if strcmpi(varargin(m),'smap')
        % Scalp map [ic,64*64] and take abs(corr) for polarity inversion
        [n,d1,d2]=size(cell2mat(varargin(m+1)));
        if d1~=d2
            error('scalp maps must be square matrices');
        else
            data = cell2mat(varargin(m+1));
            disp('clustering components using scalp maps')
       end
        % for each ic compute the 2D correlation with other maps
        pairs = nchoosek([1:n],2);
        S = eye(n);
        for p=1:length(pairs)
            S(pairs(p,1),pairs(p,2)) = corr2(squeeze(data(pairs(p,1),:,:)),squeeze(data(pairs(p,2),:,:)));
            S(pairs(p,2),pairs(p,1)) = S(pairs(p,1),pairs(p,2));
        end
        M{index} = S; clear S data 
        C{index} = apcluster(M{index},median(M{index}));
        varargin(m+1) = {[]}; index = index+1; 
    end
    
    if strcmpi(varargin(m),'dip')
        % compute euclidian distance between diploles and normalize to 1
        data = cell2mat(varargin(m+1));
        D = squareform(pdist(data(:,[1 2 3],'euclidean')));
        M{index} = D ./ (max(D(:))); clear D
        C{index} = apcluster(M{index},median(M{index}));
        index = index+1;
        % compute the abs(angle) between orientations
        if size(data,2) == 6
            pairs = nchoosek([1:size(data,1)],2);
            S = eye(size(data,1));
            for p=1:length(pairs)
                v1 = data(pairs(p,1),[4 5 6]); v2 = data(pairs(p,2),[4 5 6]);
                S(pairs(p,1),pairs(p,2)) = abs(acosd(dot(v1,v2)/(norm(v1)*norm(v2))));
                S(pairs(p,2),pairs(p,1)) = S(pairs(p,1),pairs(p,2));
            end
            M{index} = S./max(S(:)); clear S
            C{index} = apcluster(M{index},median(M{index}));
            index = index+1;
        end
        clear data
        varargin(m+1) = {[]};
    end
end % closes the varargin loop

% we have a series of similarity matrices M and clustering matrices C
% the final clustering is the intersection of Cs

% similar components should have similar time courses, spectra, ersp, and
% itc so we can either take the intersection of these or take the mean
% similarity matrix and cluster this one

% similar components should also have similar origin, orientation and scalp
% topography so again we can take either the intersection of these or take
% the mean similarity matrix and cluster this one


end

function S = local_xcorr(data,lag,type)
% computes the cross correlation between columns with the specified
% lag and returns a correlation matrix S with the maximum value across
% all lags - S is of dim [size(data,2) size(data,2)]
% type is either 'abs' or 'signed' meaning one take either the max
% of absolute values or we take the max of signed values

% xcorr conputes the conjugate of matrices, which can be very large for
% matrices with many ic, leading to memory issue, instead we loop per lag
% creating lagged matrices and compute the correlation

% get max of abs cross corrrelation
if strcmp(type,'abs')
    out = max(abs(xcorr(data,lag,'coeff')));
elseif strcmp(type,'signed')
    out = max(xcorr(data,lag,'coeff'));
else
    error('error calling the subfunction local_xcorr, unidentified arg ''type''')
end

% reshape
S = reshape(out',size(data,2),size(data,2));
end
