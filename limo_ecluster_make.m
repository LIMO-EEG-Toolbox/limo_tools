function [th,boot_values] = limo_ecluster_make(bootf,bootp,alphav)

% ECLUSTER_MAKE creates a multivariate, cluster based, statistical threshold based 
% on bootstrapped F tests. The function operates along the last dimension of the 
% input data. For each bootstrap, independently at each electrode:  
%       (1) significant F values are clustered 
%       (2) the sum of F values inside each cluster is computed
%       (3) the maximum sum across clusters is saved.
% The maximum sums of clusters of significant F values are then sorted to
% obtain a (1-alpha)% percentile threshold.
% 
% NOTE: for a two-tailed bootstrap t-test, enter a squared matrix of T values:
% th = ecluster_make(boott.^2,bootp,alpha)
%
% FORMAT     th = limo_ecluster_make(bootf,bootp,alpha)
% 
% INPUTS:
%           bootf = 2D (frames x bootstrap) or 3D matrix (eg. elec * time *
%                   bootstrap) of F values format electrode x
%           bootp = matrix of p values associated with bootf, with same dim
%           alpha = type I error rate, default 0.05
%
% OUTPUT:
%           th is a structure with percentile thresholds 
%           TH.ELEC = cluster threshold at each electrode
%           TH.MAX = max cluster threshold across the 1st dimension
%                    this is a more conservative way to control for multiple
%                    comparisons than using a spatial-temporal clustering 
%                    when the full space is analyzed.
%          boot_values are all the max cluster sums computed under H0
%
% v1 Guillaume Rousselet, University of Glasgow, August 2010
% Luisa Frei, 4 Nov 2011: fixed bug in all electrode ouput
% Marianne Latinus adding spm_bwlabel 2012
% Cyril Pernet - removed some useless computations to speed things up, June 2014
% GAR - commented out last line which crashed the function - September 2015
% CP added boot_values to then compute the p values
% ----------------------------------------------------------------------------
%  Copyright (C) LIMO Team 2016
%
% See also LIMO_TFCLUSTER_MAKE LIMO_ECLUSTER_TEST

cluster_val = [];
if nargin < 3
    alphav = 0.05;
end

if ndims(bootf)==3 % electrode*time/freq*boot 
    b = size(bootf,3);
    U = round((1-alphav)*b);
    Ne = size(bootf,1);
    boot_values = zeros(b,Ne);

    for E = 1:Ne % electrodes/freq
        for kk=1:b % bootstrap samples

            % get cluster along 1st dim
            if exist('spm_bwlabel','file')~=0
                [L,NUM] = spm_bwlabel(squeeze(bootp(E,:,kk))<=alphav,6);
            else
                if exist('bwlabeln','file')~=0
                    [L,NUM] = bwlabeln(squeeze(bootp(E,:,kk))<=alphav);
                else
                    errordlg('You need either the Image Processing Toolbox or SPM in your path to execute this function');
                end
            end
            
            if NUM~=0
                tmp=zeros(1,NUM);
                for C = 1:NUM % compute sum for each cluster
                    tmp(C) = sum(squeeze(bootf(E,L==C,kk)) );
                end
                boot_values(E,kk) = max(tmp(:)); % save max across clusters
            else
                boot_values(E,kk) = 0;
            end
        end % bootstrap loop
    end % electrode loop

    sortSC = sort(boot_values,2);
    th.elec = sortSC(:,U); % threshold at each electrode
    maxSC = max(boot_values,[],1); % max across electrodes
    sortmaxSC = sort(maxSC);
    th.max = sortmaxSC(U); % threshold of max across electrodes
    

elseif ndims(bootf)==2 % 1 electrode * time/freq

    b = size(bootf,2);
    U = round((1-alphav)*b);
    boot_values = zeros(b,1);

    for kk=1:b % bootstrap samples
        [L,NUM] = bwlabeln(squeeze(bootp(:,kk))<=alphav); % find clusters

        if NUM~=0
            tmp=zeros(1,NUM);
            for C = 1:NUM % compute sum for each cluster
                tmp(C) = sum( squeeze(bootf(L==C,kk)) );
            end
            boot_values(kk) = max(tmp(:)); % save max across clusters
        else
            boot_values(kk) = 0;
        end

    end % bootstrap loop

    sortSC = sort(boot_values);
    th.elec = sortSC(U); % threshold at the unique electrode

else
    error('ecluster_make: bootf should be 2D or 3D')
end

