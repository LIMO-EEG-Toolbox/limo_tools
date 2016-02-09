function tfce_score = limo_tfce(varargin)

% implementation of the Threshold-free cluster enhancement method
% developped for fMRI by Smith & Nichols, NeuroImage 44(2009), 83-98
% tfce = sum(extent(h)^E*height^H*dh)
%
% INPUT tfce_score = limo_tfce(type,data,channeighbstructmat)
%       tfce_score = limo_tfce(type,data,channeighbstructmat,updatebar,E,H,dh)
%
%       type = 1 for 1D data (one channel ERP or Power),
%              2 for 2D data (ERP, Power, a single freq*time map),
%              3 for 3D data (ERSP)
%       data can be either a map of t/F values or a set of t/F maps computed under H0 (in last dim)
%       channeighbstructmat is the neighbourhood matrix for clustering - if empty for type 2, siwtch to bwlabel = freq*time map
%       updatebar is a flag (default = 1) to produce a waitbar
%       E, H and dh are the parameters of the tfce algorithm defaults are 0.5, 2, 0.1
%
%
% OUPUT tfce_score is a map of scores
%
% References
%
% Pernet, C., Latinus, M., Nichols, T.E., & Rousselet, G.A. (2015)
% Cluster-based computational methods for mass univariate analyses
% of event-related brain potentials/fields: a simulation study
% Journal Of Neuroscience Method 250, Pages 85–93
% <10.1016/j.jneumeth.2014.08.003>
%
% Pernet, Cyril; Rousselet, Guillaume (2014): Type 1 error rate using TFCE for ERP.
% figshare. http://dx.doi.org/10.6084/m9.figshare.1008325
%
% Cyril Pernet v4 28-07-2015
% fixed indices / got the loop faster /
% V5 20-08-2015
% use limo_findcluster which is faster (clustering speed x60)
% changed the integration from a loop to hist - thx to Bruno Giordano
% --------------------------------------
% Copyright (C) LIMO Team 2016

% precision max = 200; % define how many thresholds between min t/F map and
% max t/F map --> needed as sometime under H0 some values can be
% arbritrarily high due to low variance during resampling

%% check input

if nargin < 3
    error('not enough arguments')
elseif nargin == 3 || nargin == 4
    if nargin == 3
        updatebar = 1;
    else
        updatebar = varargin{4};
    end
    E = 0.5;
    H = 2;
    dh = 0.1;
elseif nargin == 7
    updatebar = varargin{4};
    E = varargin{5};
    H = varargin{6};
    dh = varargin{7};
elseif nargin > 8
    error('too many arguments')
end

type = varargin{1};
data = varargin{2};
channeighbstructmat = varargin{3};
clear varargin

%% start tcfe


switch type
    
    % ---------------------------------------------------------------------
    case{1}  % 1D data -- needs to be updated (doesn't work)
        % ---------------------------------------------------------------------
        
        if isvector(data)
            [~,x]=size(data);
            subtype = 1;
        else
            [x,b]=size(data);
            subtype = 2;
        end
        
        switch subtype
            
            case{1}                
                % ------- tfce real data -----------
                
                % define increment size forced by dh
                data_range = range(data(:));
                if data_range > 1
                    precision = round(data_range / dh);
                    if precision > 200 % arbitrary decision to limit precision to 200th of the data range - needed as sometime under H0 one value can be very wrong
                        increment = data_range / 200;
                    else
                        increment = data_range / precision;
                    end
                else
                    increment = data_range *dh;
                end
                
                % check negative values if so do negate and add scores
                if min(data(:)) >= 0
                    
                    % select a height, obtain cluster map, obtain extent map (=cluster
                    % map but with extent of cluster rather than number of the cluster)
                    % then tfce score for that height
                    index = 1;
                    tfce = NaN(1,x,length(min(data(:)):increment:max(data(:))));
                    if updatebar ==1;
                        f = waitbar(0,'Thresholding levels','name','TFCE');
                    end
                    nsteps = length(min(data(:)):increment:max(data(:)));
                    for h=min(data(:)):increment:max(data(:))
                        if updatebar ==1; waitbar(index/nsteps); end
                        [clustered_map, num] = bwlabel((data > h));
                        extent_map = zeros(1,x); % same as cluster map but contains extent value instead
                        extent_map = integrate(clustered_map,num,extent_map);
                        % for i=1:num
                        %    idx = clustered_map(:) == i;
                        %    extent_map(idx) = sum(idx);
                        % end
                        tfce(1,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    % compute final score
                    tfce_score = nansum(tfce,3);
                    try close(f); end
                    
                else
                    
                    pos_data = (data > 0).*data;
                    neg_data = abs((data < 0).*data);
                    
                    if updatebar ==1;
                        f = waitbar(0,'Thresholding levels','name','TFCE');
                    end
                    nsteps = length(min(data(:)):increment:max(data(:)));
                    clear data
                    
                    % select a height, obtain cluster map, obtain extent map
                    % then tfce score for that height
                    l = length(min(pos_data(:)):increment:max(pos_data(:)));
                    pos_increment = (max(pos_data(:)) - min(pos_data(:))) / l;
                    pos_tfce = NaN(1,x,l); index = 1;
                    for h=min(pos_data(:)):pos_increment:max(pos_data(:))
                        if updatebar ==1; waitbar(index/nsteps); end
                        [clustered_map, num] = bwlabel((pos_data > h));
                        extent_map = zeros(1,x);
                        extent_map = integrate(clustered_map,num,extent_map);
                        % for i=1:num
                        %     idx = clustered_map(:) == i;
                        %    extent_map(idx) = sum(idx);
                        % end
                        pos_tfce(1,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    hindex = index-1;
                    l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
                    neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
                    neg_tfce = NaN(1,x,l); index = 1;
                    for h=min(neg_data(:)):neg_increment:max(neg_data(:))
                        if updatebar ==1; waitbar((hindex+index)/nsteps); end
                        [clustered_map, num] = bwlabel((neg_data > h));
                        extent_map = zeros(1,x);
                        extent_map = integrate(clustered_map,num,extent_map);
                        % for i=1:num
                        %    idx = clustered_map(:) == i;
                        %    extent_map(idx) = sum(idx);
                        % end
                        neg_tfce(1,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    % compute final score
                    tfce_score = nansum(pos_tfce,3)+nansum(neg_tfce,3);
                    try close(f); end
                end
                
                
            case{2}
                % ------- tfce bootstrapped data under H0 --------------
                tfce_score = NaN(1,x,b);
                
                % check negative values if so do negate and add scores
                if min(data(:)) > 0
                    
                    % select a height, obtain cluster map, obtain extent map
                    % then tfce score for that height
                    if updatebar ==1;
                        f = waitbar(0,'percentage of bootstraps analyzed','name','TFCE');
                    end
                    
                    for boot=1:b
                        if updatebar ==1; waitbar(boot/b); end
                        tmp_data = squeeze(data(:,boot));
                        % define increment size forced by dh
                        data_range = range(tmp_data(:));
                        if data_range > 1
                            precision = round(data_range / dh);
                            if precision > 200
                                increment = data_range / 200;
                            else
                                increment = data_range / precision;
                            end
                        else
                            increment = data_range *dh;
                        end
                        
                        index = 1;
                        tfce = NaN(1,x,length(min(tmp_data(:)):increment:max(tmp_data(:))));
                        % fprintf('estimating tfce under H0 boot %g \n',boot)
                        
                        for h=min(tmp_data(:)):increment:max(tmp_data(:))
                            [clustered_map, num] = bwlabel((tmp_data > h));
                            extent_map = zeros(1,x);
                            extent_map = integrate(clustered_map,num,extent_map);
                            % for i=1:num
                            %     idx = clustered_map(:) == i;
                            %     extent_map(idx) = sum(idx);
                            % end
                            tfce(1,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        tfce_score(1,:,boot) = nansum(tfce,3);
                    end
                    try close(f); end
                    
                else
                    
                    if updatebar ==1;
                        f = waitbar(0,'percentage of bootstraps analyzed','name','TFCE');
                    end
                    
                    for boot=1:b
                        
                        if updatebar ==1; waitbar(boot/b); end
                        % fprintf('estimating tfce under H0 for boot %g \n',boot)
                        tmp_data = squeeze(data(:,boot));
                        
                        % define increment size
                        data_range = range(tmp_data(:));
                        if data_range > 1
                            precision = round(data_range / dh);
                            if precision > 200
                                increment = data_range / 200;
                            else
                                increment = data_range / precision;
                            end
                        else
                            increment = data_range *dh;
                        end
                        
                        pos_data = (tmp_data > 0).*tmp_data;
                        neg_data = abs((tmp_data < 0).*tmp_data);
                        clear tmp_data
                        
                        % select a height, obtain cluster map, obtain extent map
                        % then tfce score for that height
                        l = length(min(pos_data(:)):increment:max(pos_data(:)));
                        pos_increment = (max(pos_data(:)) - min(pos_data(:))) / l;
                        pos_tfce = NaN(1,x,l); index = 1;
                        for h=min(pos_data(:)):pos_increment:max(pos_data(:))
                            [clustered_map, num] = bwlabel((pos_data > h));
                            extent_map = zeros(1,x);
                            extent_map = integrate(clustered_map,num,extent_map);
                            % for i=1:num
                            %     idx = clustered_map(:) == i;
                            %     extent_map(idx) = sum(idx);
                            % end
                            pos_tfce(1,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        
                        l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
                        neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
                        neg_tfce = NaN(1,x,l); index = 1;
                        for h=min(neg_data(:)):neg_increment:max(neg_data(:))
                            [clustered_map, num] = bwlabel((neg_data > h));
                            extent_map = zeros(1,x);
                            extent_map = integrate(clustered_map,num,extent_map);
                            % for i=1:num
                            %     idx = clustered_map(:) == i;
                            %     extent_map(idx) = sum(idx);
                            % end
                            neg_tfce(1,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        
                        % compute final score
                        tfce_score(1,:,boot) = nansum(pos_tfce,3)+nansum(neg_tfce,3);
                    end
                    try close(f); end
                end
                
        end
        
        
        % ---------------------------------------------------------------------
    case{2} % 2D data
        % ---------------------------------------------------------------------
        
        
        [x,y,b]=size(data);
        if b == 1
            subtype = 1;
        else
            subtype = 2;
        end
        
        switch subtype
            
            case{1}
                
                % ------- tfce real data -----------
                
                % define increment size forced by dh
                data_range = range(data(:));
                if data_range > 1
                    precision = round(data_range / dh);
                    if precision > 200
                        increment = data_range / 200;
                    else
                        increment = data_range / precision;
                    end
                else
                    increment = data_range *dh;
                end
                
                % check negative values if so do negate and add scores
                if min(data(:)) > 0
                    
                    % select a height, obtain cluster map, obtain extent map (=cluster
                    % map but with extent of cluster rather than number of the cluster)
                    % then tfce score for that height
                    index = 1;
                    tfce = NaN(x,y,length(min(data(:)):increment:max(data(:))));
                    if updatebar ==1;
                        f = waitbar(0,'Thresholding levels','name','TFCE');
                    end
                    nsteps = length(min(data(:)):increment:max(data(:)));
                    for h=min(data(:)):increment:max(data(:))
                        if updatebar ==1; waitbar(index/nsteps); end
                        if isempty(channeighbstructmat)
                            [clustered_map, num] = bwlabel((data > h),4);
                        else
                            [clustered_map, num] = limo_findcluster((data > h), channeighbstructmat,2);
                        end
                        
                        extent_map = zeros(x,y); % same as cluster map but contains extent value instead
                        extent_map = integrate(clustered_map,num,extent_map);
                        % for i=1:num
                        %    idx = clustered_map(:) == i;
                        %     extent_map(idx) = sum(idx);
                        % end
                        tfce(:,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    % compute final score
                    tfce_score = nansum(tfce,3);
                    try close(f); end
                    
                else
                    
                    pos_data = (data > 0).*data;
                    neg_data = abs((data < 0).*data);
                    
                    if updatebar ==1;
                        f = waitbar(0,'Thresholding levels','name','TFCE');
                    end
                    nsteps = length(min(data(:)):increment:max(data(:)));
                    clear data
                    
                    % select a height, obtain cluster map, obtain extent map
                    % then tfce score for that height
                    l = length(min(pos_data(:)):increment:max(pos_data(:)));
                    pos_increment = (max(pos_data(:)) - min(pos_data(:))) / l;
                    pos_tfce = NaN(x,y,l); index = 1;
                    for h=min(pos_data(:)):pos_increment:max(pos_data(:))
                        if updatebar ==1; waitbar(index/nsteps); end
                        if isempty(channeighbstructmat)
                            [clustered_map, num] = bwlabel((pos_data > h),4);
                        else
                            [clustered_map, num] = limo_findcluster((pos_data > h), channeighbstructmat,2);
                        end
                        
                        extent_map = zeros(x,y);
                        extent_map = integrate(clustered_map,num,extent_map);
                        % for i=1:num
                        %     idx = clustered_map(:) == i;
                        %     extent_map(idx) = sum(idx);
                        % end
                        pos_tfce(:,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    hindex = index-1;
                    l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
                    neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
                    neg_tfce = NaN(x,y,l); index = 1;
                    for h=min(neg_data(:)):neg_increment:max(neg_data(:))
                        if updatebar ==1; waitbar((hindex+index)/nsteps); end
                        if isempty(channeighbstructmat)
                            [clustered_map, num] = bwlabel((neg_data > h),4);
                        else
                            [clustered_map, num] = limo_findcluster((neg_data > h), channeighbstructmat,2);
                        end
                        
                        extent_map = zeros(x,y);
                        extent_map = integrate(clustered_map,num,extent_map);
                        % for i=1:num
                        %     idx = clustered_map(:) == i;
                        %     extent_map(idx) = sum(idx);
                        % end
                        neg_tfce(:,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    % compute final score
                    tfce_score = nansum(pos_tfce,3)+nansum(neg_tfce,3);
                    try close(f); end
                end
                
                
            case{2}
                % ------- tfce bootstrapped data under H0 --------------
                tfce_score = NaN(x,y,b);
                
                % check negative values if so do negate and add scores
                if min(data(:)) > 0
                    
                    % select a height, obtain cluster map, obtain extent map
                    % then tfce score for that height
                    if updatebar ==1;
                        f = waitbar(0,'percentage of bootstraps analyzed','name','TFCE');
                    end
                    
                    for boot=1:b
                        if updatebar ==1; waitbar(boot/b); end
                        tmp_data = squeeze(data(:,:,boot));
                        % define increment size
                        data_range = range(tmp_data(:));
                        if data_range > 1
                            precision = round(data_range / dh);
                            if precision > 200
                                increment = data_range / 200;
                            else
                                increment = data_range / precision;
                            end
                        else
                            increment = data_range *dh;
                        end
                        
                        index = 1; tfce = NaN(x,y,length(min(tmp_data(:)):increment:max(tmp_data(:))));
                        % fprintf('estimating tfce under H0 boot %g \n',boot)
                        
                        for h=min(tmp_data(:)):increment:max(tmp_data(:))
                            if isempty(channeighbstructmat)
                                [clustered_map, num] = bwlabel((tmp_data > h),4);
                            else
                                [clustered_map, num] = limo_findcluster((tmp_data > h), channeighbstructmat,2);
                            end
                            
                            extent_map = zeros(x,y);
                            extent_map = integrate(clustered_map,num,extent_map);
                            % for i=1:num
                            %    idx = clustered_map(:) == i;
                            %    extent_map(idx) = sum(idx);
                            % end
                            tfce(:,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        tfce_score(:,:,boot) = nansum(tfce,3);
                    end
                    try close(f); end
                    
                else
                    
                    if updatebar ==1;
                        f = waitbar(0,'percentage of bootstraps analyzed','name','TFCE');
                    end
                    
                    for boot=1:b
                        
                        if updatebar ==1; waitbar(boot/b); end
                        % fprintf('estimating tfce under H0 for boot %g \n',boot)
                        tmp_data = squeeze(data(:,:,boot));
                        
                        % define increment size
                        data_range = range(tmp_data(:));
                        if data_range > 1
                            precision = round(data_range / dh);
                            if precision > 200
                                increment = data_range / 200;
                            else
                                increment = data_range / precision;
                            end
                        else
                            increment = data_range *dh;
                        end
                        
                        pos_data = (tmp_data > 0).*tmp_data;
                        neg_data = abs((tmp_data < 0).*tmp_data);
                        clear tmp_data
                        
                        % select a height, obtain cluster map, obtain extent map
                        % then tfce score for that height
                        l = length(min(pos_data(:)):increment:max(pos_data(:)));
                        pos_increment = (max(pos_data(:)) - min(pos_data(:))) / l;
                        pos_tfce = NaN(x,y,l); index = 1;
                        for h=min(pos_data(:)):pos_increment:max(pos_data(:))
                            if isempty(channeighbstructmat)
                                [clustered_map, num] = bwlabel((pos_data > h),4);
                            else
                                [clustered_map, num] = limo_findcluster((pos_data > h), channeighbstructmat,2);
                            end
                            
                            extent_map = zeros(x,y);
                            extent_map = integrate(clustered_map,num,extent_map);
                            % for i=1:num
                            %     idx = clustered_map(:) == i;
                            %     extent_map(idx) = sum(idx);
                            % end
                            pos_tfce(:,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        
                        l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
                        neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
                        neg_tfce = NaN(x,y,l); index = 1;
                        for h=min(neg_data(:)):neg_increment:max(neg_data(:))
                            if isempty(channeighbstructmat)
                                [clustered_map, num] = bwlabel((neg_data > h),4);
                            else
                                [clustered_map, num] = limo_findcluster((neg_data > h), channeighbstructmat,2);
                            end
                            
                            extent_map = zeros(x,y);
                            extent_map = integrate(clustered_map,num,extent_map);
                            % for i=1:num
                            %    idx = clustered_map(:) == i;
                            %    extent_map(idx) = sum(idx);
                            % end
                            neg_tfce(:,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        
                        % compute final score
                        tfce_score(:,:,boot) = nansum(pos_tfce,3)+nansum(neg_tfce,3);
                    end
                    try close(f); end
                end
                
        end
        
        
        % ---------------------------------------------------------------------
    case{3} % 3D data
        % ---------------------------------------------------------------------
        
        [x,y,z,b]=size(data);
        if b == 1
            subtype = 1;
        else
            subtype = 2;
        end
        
        switch subtype
            
            case{1}
                % ------- tfce real data -----------
                
                % define increment size forced by dh
                data_range = range(data(:));
                if data_range > 1
                    precision = round(data_range / dh);
                    if precision > 200
                        increment = data_range / 200;
                    else
                        increment = data_range / precision;
                    end
                else
                    increment = data_range *dh;
                end
                
                % check negative values if so do negate and add scores
                if min(data(:)) > 0
                    
                    % select a height, obtain cluster map, obtain extent map (=cluster
                    % map but with extent of cluster rather than number of the cluster)
                    % then tfce score for that height
                    index = 1;
                    tfce = NaN(x,y,z,length(min(data(:)):increment:max(data(:))));
                    if updatebar ==1;
                        f = waitbar(0,'Thresholding levels','name','TFCE');
                    end
                    
                    nsteps = length(min(data(:)):increment:max(data(:)));
                    for h=min(data(:)):increment:max(data(:))
                        if updatebar ==1; waitbar(index/nsteps); end
                        try
                            [clustered_map, num] = limo_findcluster((data > h), channeighbstructmat,2);
                        catch
                            [clustered_map,num] = bwlabel((data > h)); % this allow continuous mapping
                        end
                        
                        extent_map = zeros(x,y,z);
                        extent_map = integrate(clustered_map,num,extent_map);
                        % for i=1:num
                        %    idx = clustered_map(:) == i;
                        %    extent_map(idx) = sum(idx);
                        % end
                        tfce(:,:,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    % compute final score
                    tfce_score = nansum(tfce,4);
                    try close(f); end
                    
                else
                    
                    pos_data = (data > 0).*data;
                    neg_data = abs((data < 0).*data);
                    
                    if updatebar ==1;
                        f = waitbar(0,'Thresholding levels','name','TFCE');
                    end
                    
                    nsteps = length(min(data(:)):increment:max(data(:)));
                    clear data
                    
                    % select a height, obtain cluster map, obtain extent map
                    % then tfce score for that height
                    l = length(min(pos_data(:)):increment:max(pos_data(:)));
                    pos_increment = (max(pos_data(:)) - min(pos_data(:))) / l;
                    pos_tfce = NaN(x,y,z,l); index = 1;
                    for h=min(pos_data(:)):pos_increment:max(pos_data(:))
                        if updatebar ==1; waitbar(index/nsteps); end
                        try
                            [clustered_map, num] = limo_findcluster((pos_data > h), channeighbstructmat,2);
                        catch
                            [clustered_map,num] = bwlabel((pos_data > h));
                        end
                        
                        extent_map = zeros(x,y,z);
                        extent_map = integrate(clustered_map,num,extent_map);
                        % for i=1:num
                        %     idx = clustered_map(:) == i;
                        %     extent_map(idx) = sum(idx);
                        % end
                        pos_tfce(:,:,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    hindex = index-1;
                    l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
                    neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
                    neg_tfce = NaN(x,y,z,l); index = 1;
                    for h=min(neg_data(:)):neg_increment:max(neg_data(:))
                        if updatebar ==1; waitbar((hindex+index)/nsteps); end
                        try
                            [clustered_map, num] = limo_findcluster((neg_data > h), channeighbstructmat,2);
                        catch
                            [clustered_map,num] = bwlabel((neg_data > h));
                        end
                        
                        extent_map = zeros(x,y,z);
                        extent_map = integrate(clustered_map,num,extent_map);
                        % for i=1:num
                        %     idx = clustered_map(:) == i;
                        %     extent_map(idx) = sum(idx);
                        % end
                        neg_tfce(:,:,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    % compute final score
                    tfce_score = nansum(pos_tfce,4)+nansum(neg_tfce,4);
                    try close(f); end
                end
                
                
            case{2}
                % ------- tfce bootstrapped data under H0 --------------
                tfce_score = NaN(x,y,z,b);
                
                % check negative values if so do negate and add scores
                if min(data(:)) > 0
                    
                    % select a height, obtain cluster map, obtain extent map
                    % then tfce score for that height
                    if updatebar ==1;
                        f = waitbar(0,'percentage of bootstraps analyzed','name','TFCE');
                    end
                    
                    for boot=1:b
                        if updatebar ==1; waitbar(boot/b); end
                        tmp_data = squeeze(data(:,:,:,boot));
                        % define increment size
                        data_range = range(tmp_data(:));
                        if data_range > 1
                            precision = round(data_range / dh);
                            if precision > 200
                                increment = data_range / 200;
                            else
                                increment = data_range / precision;
                            end
                        else
                            increment = data_range *dh;
                        end
                        
                        index = 1; tfce = NaN(x,y,z,length(min(tmp_data(:)):increment:max(tmp_data(:))));
                        % fprintf('estimating tfce under H0 boot %g \n',boot)
                        
                        for h=min(tmp_data(:)):increment:max(tmp_data(:))
                            try
                                [clustered_map, num] = limo_findcluster((tmp_data > h), channeighbstructmat,2);
                            catch
                                [clustered_map,num] = bwlabel((tmp_data > h));
                            end
                            
                            extent_map = zeros(x,y,z);
                            extent_map = integrate(clustered_map,num,extent_map);
                            % for i=1:num
                            %     idx = clustered_map(:) == i;
                            %     extent_map(idx) = sum(idx);
                            % end
                            tfce(:,:,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        tfce_score(:,:,:,boot) = nansum(tfce,4);
                    end
                    try close(f); end
                    
                else
                    
                    if updatebar ==1;
                        f = waitbar(0,'percentage of bootstraps analyzed','name','TFCE');
                    end
                    
                    for boot=1:b
                        
                        if updatebar ==1; waitbar(boot/b); end
                        % fprintf('estimating tfce under H0 for boot %g \n',boot)
                        tmp_data = squeeze(data(:,:,:,boot));
                        
                        % define increment size
                        data_range = range(tmp_data(:));
                        if data_range > 1
                            precision = round(data_range / dh);
                            if precision > 200
                                increment = data_range / 200;
                            else
                                increment = data_range / precision;
                            end
                        else
                            increment = data_range *dh;
                        end
                        
                        pos_data = (tmp_data > 0).*tmp_data;
                        neg_data = abs((tmp_data < 0).*tmp_data);
                        clear tmp_data
                        
                        % select a height, obtain cluster map, obtain extent map
                        % then tfce score for that height
                        l = length(min(pos_data(:)):increment:max(pos_data(:)));
                        pos_increment = (max(pos_data(:)) - min(pos_data(:))) / l;
                        pos_tfce = NaN(x,y,z,l); index = 1;
                        for h=min(pos_data(:)):pos_increment:max(pos_data(:))
                            try
                                [clustered_map, num] = limo_findcluster((pos_data > h), channeighbstructmat,2);
                            catch
                                [clustered_map,num] = bwlabel((pos_data > h));
                            end
                            
                            extent_map = zeros(x,y,z);
                            extent_map = integrate(clustered_map,num,extent_map);
                            % for i=1:num
                            %     idx = clustered_map(:) == i;
                            %     extent_map(idx) = sum(idx);
                            % end
                            pos_tfce(:,:,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        
                        l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
                        neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
                        neg_tfce = NaN(x,y,z,l); index = 1;
                        for h=min(neg_data(:)):neg_increment:max(neg_data(:))
                            try
                                [clustered_map, num] = limo_findcluster((neg_data > h), channeighbstructmat,2);
                            catch
                                [clustered_map,num] = bwlabel((neg_data > h));
                            end
                            
                            extent_map = zeros(x,y,z);
                            extent_map = integrate(clustered_map,num,extent_map);
                            % for i=1:num
                            %     idx = clustered_map(:) == i;
                            %     extent_map(idx) = sum(idx);
                            % end
                            neg_tfce(:,:,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        
                        % compute final score
                        tfce_score(:,:,:,boot) = nansum(pos_tfce,4)+nansum(neg_tfce,4);
                    end
                    try close(f); end
                end
                
        end
end
end

%% faster integration a la Bruno Giordano
function extent_map = integrate(clustered_map,num,extent_map)

clustered_map=clustered_map(:);
nv=histc(clustered_map,0:num);
[~,idxall]=sort(clustered_map,'ascend');
idxall(1:nv(1))=[];
nv(1)=[];
ends=cumsum(nv);
inis=ends-nv+1;
for i=1:num
    idx=idxall(inis(i):ends(i));
    extent_map(idx)=nv(i);
end
end

