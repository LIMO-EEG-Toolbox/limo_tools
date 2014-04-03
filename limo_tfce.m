function tfce_score = limo_tfce(varargin)

% implementation of the Threshold-free cluster enhancement method
% developped for fMRI by Smith & Nichols, NeuroImage 44(2009), 83-98
%
% INPUT tfce_score = limo_tfce(type,data,channeighbstructmat)
%       tfce_score = limo_tfce(type,data,channeighbstructmat,E,H,dh)
%
%       type = 1 for 1D data (one channel), 2 for 2D data (ERP, Power), 3 for 3D data (ESP)
%       data can be either a map of t/F values or a set of t/F maps computed under H0 (in last dim)
%       E, H and dh are the parameters of the tfce algorithm defaults are 0.5, 2, 0.1
%       tfce = sum(extent(h)^E*height^H*dh)
%
% OUPUT tfce_score is a map of scores
%
% Ref
% Pernet, C., Latinus, M., Nichols, T.E., & Rousselet, G.A.
% Cluster-based computational methods for mass univariate analyses
% of event-related brain potentials/fields: a simulation study
% Journal Of Neuroscience Method - submitted
%
% Cyril Pernet v3 03-04-2014
% fixed indices / got the loop faster / 
% --------------------------------------
% Copyright (C) LIMO Team 2010


% precision max = 200; % define how many thresholds between min t/F map and
% max t/F map --> needed as sometime under H0 some values can be
% arbritrarily high due to low variance during resampling

%% check input

if nargin < 3
    error('not enough arguments')
elseif nargin == 3
    E = 0.5;
    H = 2;
    dh = 0.1;
elseif nargin == 6
    E = varargin{4};
    H = varargin{5};
    dh = varargin{6};
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
                    f = waitbar(0,'Thresholding levels','name','TFCE');
                    nsteps = length(min(data(:)):increment:max(data(:)));
                    for h=min(data(:)):increment:max(data(:))
                        waitbar(index/nsteps);
                        [clustered_map, num] = bwlabel((data > h));
                        extent_map = zeros(1,x); % same as cluster map but contains extent value instead
                        for i=1:num
                            idx = clustered_map(:) == i;
                            extent_map(idx) = sum(idx);
                        end
                        tfce(1,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    % compute final score
                    tfce_score = nansum(tfce,3);
                    close(f)
                    
                else
                    
                    pos_data = (data > 0).*data;
                    neg_data = abs((data < 0).*data);
                    
                    f = waitbar(0,'Thresholding levels','name','TFCE');
                    nsteps = length(min(data(:)):increment:max(data(:)));
                    clear data
                    
                    % select a height, obtain cluster map, obtain extent map
                    % then tfce score for that height
                    l = length(min(pos_data(:)):increment:max(pos_data(:)));
                    pos_increment = (max(pos_data(:)) - min(pos_data(:))) / l;
                    pos_tfce = NaN(1,x,l); index = 1;
                    for h=min(pos_data(:)):pos_increment:max(pos_data(:))
                        waitbar(index/nsteps);
                        [clustered_map, num] = bwlabel((pos_data > h));
                        extent_map = zeros(1,x); 
                        for i=1:num
                            idx = clustered_map(:) == i;
                            extent_map(idx) = sum(idx);
                        end
                        pos_tfce(1,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    hindex = index-1;
                    l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
                    neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
                    neg_tfce = NaN(1,x,l); index = 1;
                    for h=min(neg_data(:)):neg_increment:max(neg_data(:))
                        waitbar((hindex+index)/nsteps);
                        [clustered_map, num] = bwlabel((neg_data > h));
                        extent_map = zeros(1,x); 
                        for i=1:num
                            idx = clustered_map(:) == i;
                            extent_map(idx) = sum(idx);
                        end
                        neg_tfce(1,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    % compute final score
                    tfce_score = nansum(pos_tfce,3)+nansum(neg_tfce,3);
                    close(f)
                end
                
                
            case{2}
                % ------- tfce bootstrapped data under H0 --------------
                tfce_score = NaN(1,x,b);
                
                % check negative values if so do negate and add scores
                if min(data(:)) > 0
                    
                    % select a height, obtain cluster map, obtain extent map
                    % then tfce score for that height
                    f = waitbar(0,'percentage of bootstraps analyzed','name','TFCE');
                    
                    for boot=1:b
                        waitbar(boot/b);
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
                            for i=1:num
                                idx = clustered_map(:) == i;
                                extent_map(idx) = sum(idx);
                            end
                            tfce(1,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        tfce_score(1,:,boot) = nansum(tfce,3);
                    end
                    close(f)
                    
                else
                    
                    f = waitbar(0,'percentage of bootstraps analyzed','name','TFCE');
                    for boot=1:b
                        
                        waitbar(boot/b);
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
                            for i=1:num
                                idx = clustered_map(:) == i;
                                extent_map(idx) = sum(idx);
                            end
                            pos_tfce(1,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        
                        l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
                        neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
                        neg_tfce = NaN(1,x,l); index = 1;
                        for h=min(neg_data(:)):neg_increment:max(neg_data(:))
                            [clustered_map, num] = bwlabel((neg_data > h));
                            extent_map = zeros(1,x); 
                            for i=1:num
                                idx = clustered_map(:) == i;
                                extent_map(idx) = sum(idx);
                            end
                            neg_tfce(1,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        
                        % compute final score
                        tfce_score(1,:,boot) = nansum(pos_tfce,3)+nansum(neg_tfce,3);
                    end
                    close(f)
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
                    f = waitbar(0,'Thresholding levels','name','TFCE');
                    nsteps = length(min(data(:)):increment:max(data(:)));
                    for h=min(data(:)):increment:max(data(:))
                        waitbar(index/nsteps);
                        [clustered_map, num] = limo_ft_findcluster((data > h), channeighbstructmat,2);
                        extent_map = zeros(x,y); % same as cluster map but contains extent value instead
                        for i=1:num
                            idx = clustered_map(:) == i;
                            extent_map(idx) = sum(idx); 
                        end
                        tfce(:,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    % compute final score
                    tfce_score = nansum(tfce,3);
                    close(f)
                    
                else
                    
                    pos_data = (data > 0).*data;
                    neg_data = abs((data < 0).*data);
                    
                    f = waitbar(0,'Thresholding levels','name','TFCE');
                    nsteps = length(min(data(:)):increment:max(data(:)));
                    clear data
                    
                    % select a height, obtain cluster map, obtain extent map
                    % then tfce score for that height
                    l = length(min(pos_data(:)):increment:max(pos_data(:)));
                    pos_increment = (max(pos_data(:)) - min(pos_data(:))) / l;
                    pos_tfce = NaN(x,y,l); index = 1;
                    for h=min(pos_data(:)):pos_increment:max(pos_data(:))
                        waitbar(index/nsteps);
                        [clustered_map, num] = limo_ft_findcluster((pos_data > h), channeighbstructmat,2);
                        extent_map = zeros(x,y); 
                        for i=1:num
                            idx = clustered_map(:) == i;
                            extent_map(idx) = sum(idx);
                        end
                        pos_tfce(:,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    hindex = index-1;
                    l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
                    neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
                    neg_tfce = NaN(x,y,l); index = 1;
                    for h=min(neg_data(:)):neg_increment:max(neg_data(:))
                        waitbar((hindex+index)/nsteps);
                        [clustered_map, num] = limo_ft_findcluster((neg_data > h), channeighbstructmat,2);
                        extent_map = zeros(x,y); 
                        for i=1:num
                            idx = clustered_map(:) == i;
                            extent_map(idx) = sum(idx);
                        end
                        neg_tfce(:,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    % compute final score
                    tfce_score = nansum(pos_tfce,3)+nansum(neg_tfce,3);
                    close(f)
                end
                
                
            case{2}
                % ------- tfce bootstrapped data under H0 --------------
                tfce_score = NaN(x,y,b);
                
                % check negative values if so do negate and add scores
                if min(data(:)) > 0
                    
                    % select a height, obtain cluster map, obtain extent map
                    % then tfce score for that height
                    f = waitbar(0,'percentage of bootstraps analyzed','name','TFCE');
                    
                    for boot=1:b
                        waitbar(boot/b);
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
                            [clustered_map, num] = limo_ft_findcluster((tmp_data > h), channeighbstructmat,2);
                            extent_map = zeros(x,y); 
                            for i=1:num
                                idx = clustered_map(:) == i;
                                extent_map(idx) = sum(idx);
                            end
                            tfce(:,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        tfce_score(:,:,boot) = nansum(tfce,3);
                    end
                    close(f)
                    
                else
                    
                    f = waitbar(0,'percentage of bootstraps analyzed','name','TFCE');
                    for boot=1:b
                        
                        waitbar(boot/b);
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
                            [clustered_map, num] = limo_ft_findcluster((pos_data > h), channeighbstructmat,2);
                            extent_map = zeros(x,y); 
                            for i=1:num
                                idx = clustered_map(:) == i;
                                extent_map(idx) = sum(idx);
                            end
                            pos_tfce(:,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        
                        l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
                        neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
                        neg_tfce = NaN(x,y,l); index = 1;
                        for h=min(neg_data(:)):neg_increment:max(neg_data(:))
                            [clustered_map, num] = limo_ft_findcluster((neg_data > h), channeighbstructmat,2);
                            extent_map = zeros(x,y); 
                            for i=1:num
                                idx = clustered_map(:) == i;
                                extent_map(idx) = sum(idx);
                            end
                            neg_tfce(:,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        
                        % compute final score
                        tfce_score(:,:,boot) = nansum(pos_tfce,3)+nansum(neg_tfce,3);
                    end
                    close(f)
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
                    f = waitbar(0,'Thresholding levels','name','TFCE');
                    nsteps = length(min(data(:)):increment:max(data(:)));
                    for h=min(data(:)):increment:max(data(:))
                        waitbar(index/nsteps);
                        try
                            [clustered_map, num] = limo_ft_findcluster((data > h), channeighbstructmat,2);
                        catch
                            [clustered_map,num] = bwlabel((data > h)); % this allow continuous mapping
                        end
                        extent_map = zeros(x,y,z); 
                        for i=1:num
                            idx = clustered_map(:) == i;
                            extent_map(idx) = sum(idx);
                        end
                        tfce(:,:,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    % compute final score
                    tfce_score = nansum(tfce,4);
                    close(f)
                    
                else
                    
                    pos_data = (data > 0).*data;
                    neg_data = abs((data < 0).*data);
                    
                    f = waitbar(0,'Thresholding levels','name','TFCE');
                    nsteps = length(min(data(:)):increment:max(data(:)));
                    clear data
                    
                    % select a height, obtain cluster map, obtain extent map
                    % then tfce score for that height
                    l = length(min(pos_data(:)):increment:max(pos_data(:)));
                    pos_increment = (max(pos_data(:)) - min(pos_data(:))) / l;
                    pos_tfce = NaN(x,y,z,l); index = 1;
                    for h=min(pos_data(:)):pos_increment:max(pos_data(:))
                        waitbar(index/nsteps);
                        try
                            [clustered_map, num] = limo_ft_findcluster((pos_data > h), channeighbstructmat,2);
                        catch
                            [clustered_map,num] = bwlabel((pos_data > h));
                        end
                        extent_map = zeros(x,y,z); 
                        for i=1:num
                            idx = clustered_map(:) == i;
                            extent_map(idx) = sum(idx);
                        end
                        pos_tfce(:,:,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    hindex = index-1;
                    l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
                    neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
                    neg_tfce = NaN(x,y,z,l); index = 1;
                    for h=min(neg_data(:)):neg_increment:max(neg_data(:))
                        waitbar((hindex+index)/nsteps);
                        try
                            [clustered_map, num] = limo_ft_findcluster((neg_data > h), channeighbstructmat,2);
                        catch
                            [clustered_map,num] = bwlabel((neg_data > h));
                        end
                        extent_map = zeros(x,y,z); 
                        for i=1:num
                            idx = clustered_map(:) == i;
                            extent_map(idx) = sum(idx);
                        end
                        neg_tfce(:,:,:,index) = (extent_map.^E).*h^H.*increment;
                        index = index +1;
                    end
                    
                    % compute final score
                    tfce_score = nansum(pos_tfce,4)+nansum(neg_tfce,4);
                    close(f)
                end
                
                
            case{2}
                % ------- tfce bootstrapped data under H0 --------------
                tfce_score = NaN(x,y,z,b);
                
                % check negative values if so do negate and add scores
                if min(data(:)) > 0
                    
                    % select a height, obtain cluster map, obtain extent map
                    % then tfce score for that height
                    f = waitbar(0,'percentage of bootstraps analyzed','name','TFCE');
                    
                    for boot=1:b
                        waitbar(boot/b);
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
                                [clustered_map, num] = limo_ft_findcluster((tmp_data > h), channeighbstructmat,2);
                            catch
                                [clustered_map,num] = bwlabel((tmp_data > h));
                            end
                            extent_map = zeros(x,y,z); 
                            for i=1:num
                                idx = clustered_map(:) == i;
                                extent_map(idx) = sum(idx);
                            end
                            tfce(:,:,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        tfce_score(:,:,:,boot) = nansum(tfce,4);
                    end
                    close(f)
                    
                else
                    
                    f = waitbar(0,'percentage of bootstraps analyzed','name','TFCE');
                    for boot=1:b
                        
                        waitbar(boot/b);
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
                                [clustered_map, num] = limo_ft_findcluster((pos_data > h), channeighbstructmat,2);
                            catch
                                [clustered_map,num] = bwlabel((pos_data > h));
                            end
                            extent_map = zeros(x,y,z); 
                            for i=1:num
                                idx = clustered_map(:) == i;
                                extent_map(idx) = sum(idx);
                            end
                            pos_tfce(:,:,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        
                        l = length(min(neg_data(:)):increment:max(neg_data(:)))-1;
                        neg_increment = (max(neg_data(:)) - min(neg_data(:))) / l;
                        neg_tfce = NaN(x,y,z,l); index = 1;
                        for h=min(neg_data(:)):neg_increment:max(neg_data(:))
                            try
                                [clustered_map, num] = limo_ft_findcluster((neg_data > h), channeighbstructmat,2);
                            catch
                                [clustered_map,num] = bwlabel((neg_data > h));
                            end
                            extent_map = zeros(x,y,z); 
                            for i=1:num
                                idx = clustered_map(:) == i;
                                extent_map(idx) = sum(idx);
                            end
                            neg_tfce(:,:,:,index) = (extent_map.^E).*h^H.*increment;
                            index = index +1;
                        end
                        
                        % compute final score
                        tfce_score(:,:,:,boot) = nansum(pos_tfce,4)+nansum(neg_tfce,4);
                    end
                    close(f)
                end
                
        end
end
