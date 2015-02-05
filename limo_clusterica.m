function IDX = limo_clusterica(varargin)

% FORMAT IDX = eeg_clusterica(varargin)
%
% INPUT varargin is an array 'key', data, 'key;, data ...
%       the 'key' indicates which type of data are passed
%       these can be: 'erp' and data are of dim [ic time]
%                     'spect' and data are of dim [ic, freqs]
%                     'ersp' and data are of dim [ic, freqs, time]
%                     'itc' and data are of dim [ic, freqs, time]
%                     'smap' a scalp map of dim [64*64]
%                     'dipole' and data are of dim [x,y,z] for spatial lovcation
%                              or [x,y,z,xx,yy,zz] adding the orientation vector
%                              assuming a common origin (i.e. not [x y z])
%
% OUPUT IDX is the indices of the compoments after clustering
%       IDX{1} clustering based on mean corr of erp, ersp, itc ..
%       IDX{2} clustering based on dipole and scalp map
%       IDX{3} intersection of IX{1} and IDX{2} - the final result
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
% A = sin(randn(1500,200)); % 75 IC * 20 subject  200 time frames
% for erp = 1:1500
%     coef = fft(A(erp,:));
%     B(erp,:) = coef(1:100).*conj(coef(1:100))/ 100; % power spec
% end
% d = rand(75,3); D = d;
% for s=2:20
%     D = [D;d+ceil(randn(75,3))]; % create a bunch of dipoles
% end
% IDX = limo_clusterica('erp',A,'spec',B, 'dip',D);
%
% Cyril Pernet, The University of Edinburgh
% Arnaud Delorme, SCCN, INC, UCSD
% Ramon Martinez-Cancino,INC, SCCN
%% -----------------------------
% Copyright (C) LIMO Team 2015

% default lag - that is we take the max corr allowing a bit of lag in time
% frames or frequencies
lag = 2;

% note one computes xcorr on transposed data, i.e. between ic
% -----------------------------------------------------------
for m=1:2:length(varargin)
    
    if strcmpi(varargin{m},'erp')
        % the absolute of xcorr to be insensitive to polarity
        disp('clustering components using ERPs')
        M = local_xcorr(cell2mat(varargin(m+1))',lag,'abs');
        varargin{m+1} = M; clear M
    end
    
    if strcmpi(varargin(m),'spec')
        % power is already positive, take max of xcorr
        disp('clustering components using Spectra')
        M = local_xcorr(cell2mat(varargin(m+1))',lag,'signed');
        varargin{m+1} = M; clear M
    end
    
    if strcmpi(varargin(m),'ersp')
        % power is already positive, take max of xcorr
        % note the matrix M is the average of the best 5 frequency bands
        % no point looking at frequency bands with noise
        disp('clustering components using ERSP')
        data = cell2mat(varargin(m+1));
        S = NaN(size(data,1),size(data,1),size(data,2));
        for f = 1:size(data,2)
            S(:,:,f) = local_xcorr(squeeze(data(:,f,:))',lag,'signed');
            c = triu(squeeze(S(:,:,f)),1); score(f) = mean(c(:));
        end
        [~,ranking]=sort(score);
        M = mean(S(:,:,ranking(1:5)),3); clear S
        varargin{m+1} = M; clear M
    end
    
    if strcmpi(varargin(m),'itc')
        % phase can be opposed for opposed polarities, so we take the abs
        % note the matrix M is the average of the best 5 frequency bands
        % no point looking at frequency bands with noise
        disp('clustering components using ITC')
        data = cell2mat(varargin(m+1));
        S = NaN(size(data,1),size(data,1),size(data,2));
        for f = 1:size(data,2)
            S(:,:,f) = local_xcorr(squeeze(data(:,f,:))',lag,'abs');
            c = triu(squeeze(S(:,:,f)),1); score(f) = mean(c(:));
        end
        [~,ranking]=sort(score);
        M = mean(S(:,:,ranking(1:5)),3); clear S
        varargin{m+1} = M; clear M
    end
    
    if strcmpi(varargin(m),'smap')
        % Scalp map [ic,64*64] and take abs(corr) for polarity inversion
        [n,d1,d2]=size(cell2mat(varargin(m+1)));
        if d1~=d2
            error('Scalp Maps must be square matrices');
        else
            data = cell2mat(varargin(m+1));
            disp('clustering components using scalp maps')
        end
        % for each ic compute the 2D correlation with other maps
        pairs = nchoosek([1:n],2);
        M = eye(n);
        for p=1:length(pairs)
            M(pairs(p,1),pairs(p,2)) = corr2(squeeze(data(pairs(p,1),:,:)),squeeze(data(pairs(p,2),:,:)));
            M(pairs(p,2),pairs(p,1)) = M(pairs(p,1),pairs(p,2));
        end
        varargin{m+1} = M; clear M
    end
    
    if strcmpi(varargin(m),'dip')
        % compute euclidian distance between diploles and normalize to 1
        data = cell2mat(varargin(m+1));
        D = squareform(pdist(data(:,[1 2 3]),'euclidean'));
        % if only position, cluster now
        if size(data,2) == 3
            M = D ./ (max(D(:))); clear D
            % if orientation compute the angles
        elseif size(data,2) == 6
            pairs = nchoosek([1:size(data,1)],2);
            S = eye(size(data,1));
            for p=1:length(pairs)
                v1 = data(pairs(p,1),[4 5 6]); v2 = data(pairs(p,2),[4 5 6]);
                S(pairs(p,1),pairs(p,2)) = acosd(dot(v1,v2)/(norm(v1)*norm(v2)));
                S(pairs(p,2),pairs(p,1)) = S(pairs(p,1),pairs(p,2));
            end
            S = S./max(S); % normalize to 1
            D = (D+S)./2; % average distance and orientation
            M = D./max(D(:)); clear D S
        end
        varargin{m+1} = M; clear M data
    end
end % closes the varargin loop

disp('computing intersection of clusters')
% we have a series of similarity matrices M

% similar components should have similar time courses, spectra, ersp, and
% itc so we can take the mean similarity matrix and cluster this one
MM = [];
dim1 = size(cell2mat(varargin(2)),1);
for m=1:2:length(varargin)
    if any([strcmpi(varargin{m},'erp'), strcmpi(varargin{m},'spec'), strcmpi(varargin{m},'ersp'),strcmpi(varargin{m},'itc')])
        if isempty(MM)
            MM = varargin{m+1};
        else
            MM = (MM + varargin{m+1}) ./2;
        end
    end
end

if ~isempty(MM)
    IDX{1} = apcluster(MM,median(MM));
else
    IDX{1} = NaN(dim1,1);
end

% similar components should also have similar origin & orientation and scalp
% topography so again we can take the mean similarity matrix and cluster this one
MM = [];
for m=1:2:length(varargin)
    if strcmpi(varargin{m},'smap') || strcmpi(varargin{m},'dip')
        if isempty(MM)
            MM = varargin{m+1};
        else
            MM = (MM + varargin{m+1}) ./2;
        end
    end
end

if ~isempty(MM)
    IDX{2} = apcluster(MM,median(MM));
else
    IDX{2} = NaN(dim1,1);
end

% update output
% --------------
% IDX corresponds to the exemplar number, simply reindex starting at 1
% if the same IC is the exemplar for erp/spec/ersp/itc and for smap/dip
% then they have the same cluster number
for i=1:2
    out = unique(IDX{i});
    for v=1:length(out)
        IDX{i}(IDX{i}==(out(v))) = v;
    end
end

% the final clustering is thus
if ~isnan(IDX{1}(1)) && ~isnan(IDX{2}(1))
    common = IDX{1} == IDX{2};
    IDX{3} = NaN(dim1,1);
    IDX{3}(common) = IDX{1}(common);
else
    IDX{3} = NaN(dim1,1);
end
end

function S = local_xcorr(data,lag,type)
% computes the cross correlation between columns with the specified
% lag and returns a correlation matrix S with the maximum value across
% all lags - S is of dim [size(data,2) size(data,2)]
% type is either 'abs' or 'signed' meaning one take either the max
% of absolute values or we take the max of signed values (ie only positive
% correlations like for spec)

% xcorr conputes the conjugate of matrices, which can be very large for
% matrices with many ic, leading to memory issue, instead we loop per lag
% creating lagged matrices and compute the correlation using corrcoef_cell

% get max of abs cross corrrelation
index = 1;
C = NaN(size(data,2),size(data,2),length(-lag:lag));
for l=-lag:lag
    tmp = zeros(size(data,1)+abs(l),size(data,2));
    if l < 0
        for c=1:size(data,2)
            tmp(1:end+l,:) = data;
            tmp(:,c) = 0; tmp(1+abs(l):end,c) = data(:,c);
            all_but_c = setdiff([1:size(data,2)],c);
            C(c,all_but_c,index) = corr(tmp(:,c),tmp(:,all_but_c));
        end
    elseif l == 0
        C(:,:,index) = corr(data);
    elseif l > 0
        for c=1:size(data,2)
            tmp(1+l:end,:) = data;
            tmp(:,c) = 0; tmp(1:end-l,c) = data(:,c);
            all_but_c = setdiff([1:size(data,2)],c);
            C(c,all_but_c,index) = corr(tmp(:,c),tmp(:,all_but_c));
        end
    end
    index = index+1;
end

if strcmp(type,'abs')
    S = max(abs(C),[],3);
elseif strcmp(type,'signed')
    S = max(C,[],3);
else
    error('error calling the subfunction local_xcorr, unidentified arg ''type''')
end
end
