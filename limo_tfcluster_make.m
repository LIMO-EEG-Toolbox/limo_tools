function [th,boot_values] = limo_tfcluster_make(bootf,bootp,alphav)

% limo_tfcluster_make creates a multivariate, cluster based, statistical threshold based
% on bootstrapped F tests. The function operates along the last dimension of the
% input data. For each bootstrap, independently at each electrode:
%       (1) significant F values are clustered in the time*frequency domain
%       (2) the sum of F values inside each cluster is computed
%       (3) the maximum sum across clusters is saved.
% The maximum sums of clusters of significant F values are then sorted to
% obtain a (1-alpha)% percentile threshold.
%
% FORMAT     th = limo_ecluster_make(bootf,bootp,alpha)
%
% INPUTS:
%           bootf = 4D matrix (elec * freq * time * bootstrap) or 3D (freq *
%                   * time * bootstrap) of observed F values
%           bootp = observed p values
%           alpha = type I error rate, default 0.05
%
% OUTPUT:
%           th is a structure with distributions of maxima
%           such distributions allows getting a threshold and p values
%           th.elec = cluster max at each electrode
%           th.max  = max cluster across electrodes
%                    this is a more conservative way to control for multiple
%                    comparisons than using a high dimension clustering
%                    when the full space is analyzed because thresholds are
%                    smaller
%          boot_values are all the cluster sums conmputed under H0
%
% Cyril Pernet - adaptation of limo_ecluster_make for time*freqency data - June 2014
% ----------------------------------------------------------------------------
%  Copyright (C) LIMO Team 2014
%
% See also LIMO_ECLUSTER_MAKE LIMO_ECLUSTER_TEST

%% inputs
if ndims(bootf) ~= ndims(bootp)
    error('imnput matrices are of different size')
end

if nargin < 3
    alphav = 0.05;
end

if ndims(bootf)==4 % electrode * freq * time *boot
    Ne = size(bootf,1);
    b = size(bootf,4);
    boot_values = zeros(Ne,b);
    
    fprintf('computing low dimensional thresholds \n')
    for E = 1:Ne % electrodes
        parfor kk=1:b % bootstrap samples
            % get cluster for each time frequency map
            if exist('spm_bwlabel','file') == 3
                [L,NUM] = spm_bwlabel(squeeze(bootp(E,:,:,kk))<=alphav,6);
            elseif exist('bwlabeln','file')==2 
                [L,NUM] = bwlabeln(squeeze(bootp(E,:,:,kk))<=alphav);
            else
                errordlg('You need either the Image Processing Toolbox or SPM in your path to execute this function');
            end
            
            if NUM~=0
                tmp=zeros(1,NUM);
                for C = 1:NUM % compute sum for each cluster
                    data = squeeze(bootf(E,:,:,kk));
                    tmp(C) = sum(data(L==C));
                end
                boot_values(E,kk) = max(tmp(:)); % save max across clusters
            else
                boot_values(E,kk) = 0;
            end
        end % bootstrap loop
    end % electrode loop
    
    th.elec = sort(boot_values,2); % thresholds at each electrode
    maxSC = max(boot_values,[],1); % max hresholds across electrodes
    th.max = sort(maxSC);
    
elseif ndims(bootf)==3 % time * freq * bootstrap
    
    b = size(bootf,3);
    boot_values = zeros(1,b);
    
    parfor kk=1:b % bootstrap samples
        if exist('spm_bwlabel','file') == 3
            [L,NUM] = spm_bwlabel(squeeze(bootp(:,:,kk))<=alphav,6);
        elseif exist('bwlabeln','file')== 2 
            [L,NUM] = bwlabeln(squeeze(bootp(:,:,kk))<=alphav);
        else
            errordlg('You need either the Image Processing Toolbox or SPM in your path to execute this function');
        end
        
        if NUM~=0
            data = squeeze(bootf(:,:,kk));
            tmp=zeros(1,NUM);
            for C = 1:NUM % compute sum for each cluster
                tmp(C) = sum(data(L==C));
            end
            boot_values(kk) = max(tmp(:)); % save max across clusters
        else
            boot_values(kk) = 0;
        end
    end % bootstrap loop
    th.max = sort(boot_values); % thresholds at the unique electrode
    
else
    error('tfcluster_make: bootf should be 2D or 3D')
end

