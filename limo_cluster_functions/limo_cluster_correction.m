function [mask,cluster_p] = limo_cluster_correction(M,P,bootM,bootP,neighbouring_matrix,MCC,p)

% Correction for multiple testing using cluster mass. Since data are 
% correlated in space/time/frequency is makes sense to rely on clusters 
% to estimate if an an effect is probable under H0.
%
% INPUTS: M = matrix of observed F values (channel*[]*[])
%         P = matrix of observed p values
%         (note for a single electrode the format is 1*[]*[])
%         bootM = matrix of F values for data under H0 (channel*[]*[]*MC)
%         bootP = matrix of F values for data under H0
%         (note for a single electrode the format is 1*[]*[]*MC)
%         neighbouring_matrix is the bianry matrix that define how to cluster the 1st dim of M
%         MCC = 2 (spatial-temporal clustering) or 1 (temporal clustering)
%         p = threshold to apply (note this applied to create clusters and to  threshold the cluster map)
%
% OUTPUTS: maks = a binary matrix of significant/non-significant cells
%          p_vals = a matrix of cluster corrected p-values
%
% see also limo_findcluster, limo_cluster_test, limo_chancluster_test
%
% Cyril Pernet - outsourced from limo_stat_values
% ---------------------------------------------------------------------
% Copyright (C) LIMO Team 2016

if nargin == 5
    MCC = 2;
    p=0.05;
elseif nargin == 6
    p=0.05;
end

% quick dim check
if numel(size(M)) ~= numel(size(P)) ||  numel(size(bootM)) ~= numel(size(bootP))
    error('input matrices of stats and p values must be the same size')
end

if (numel(size(M))+1) ~= numel(size(bootM))
    error('dimension issue between observed and H0 input matrices')
end

if size(M,1) == 1
    MCC = 1;
end

% allocate memory
cluster_p = NaN(size(M));
mask = cluster_p;

%% compute cluster mass
ndim = numel(size(bootM));
if ndim ~=3 || ndim ~= 4
    error('H0 data must be 3D or 4D')
end

nboot = size(bootM,ndim);
U = round((1-p)*nboot); % bootstrap threshold

if MCC == 2
    minnbchan = 2;
    boot_maxclustersum=zeros(nboot,1); % compute bootstrap clusters
    disp('getting clusters under H0 boot ...');
    parfor boot = 1:nboot
        [clusterslabel,nclusters] = limo_findcluster((bootP(:,:,boot) <= p),neighbouring_matrix,minnbchan);
        bootM_b = bootM(:,:,boot);
        if nclusters~=0
            tmp=zeros(1,nclusters);
            for C = 1:nclusters % compute sum for each cluster
                tmp(C) = sum( bootM_b(clusterslabel==C) );
            end
            boot_maxclustersum(boot) = max(tmp(:)); % save max across clusters
        else
            boot_maxclustersum(boot) = 0;
        end
    end
    [mask, cluster_p] = limo_cluster_test(M,P,boot_maxclustersum,neighbouring_matrix,minnbchan,p);
    
elseif MCC == 1
    [mask, cluster_p] = limo_chancluster_test(squeeze(M),squeeze(P),squeeze(bootM),squeeze(bootP),p);
end





