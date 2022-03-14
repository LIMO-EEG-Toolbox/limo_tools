function [mask, pval, maxval, maxclustersum_th] = limo_cluster_test(ori_f,ori_p,boot_maxclustersum,channeighbstructmat,minnbchan,alphav)

% limo_cluster_test finds clusters of significant F values, computes the
% sum of F values inside each cluster, and compares that sum to a threshold
% sum of F values expected by chance (for t test - use t^2).
%
% FORMAT: [mask, pval, L, maxval, maxclustersum_th] = limo_cluster_test(ori_f,ori_p,boot_maxclustersum,channeighbstructmat,minnbchan,alphav)
%
% INPUTS: ori_f: 3D or 2D matrix of observed F values 
%         ori_p: 3D or 2D matrix of observed P values 
%         boot_maxclustersum = distribution of cluster maxima observed under H0 
%         channeighbstructmat = output of limo_ft_neighbourselection
%         minnbchan = minimum number of channels, default = 2
%         alpha level, default 0.05
%
% OUTPUTS: mask = 3D or 2D logical matrix of significant effects corrected
%                 for multiple comparisons by a cluster test
%          pval = corrected p values of significant clusters
%          maxval = maximum observed cluster mass
%          maxclustersum_th = max cluster sum (1-alpha) bootstrap threshold
%
% See also limo_clustering limo_getcluster_test limo_getclustersum
% 
% Guillaume Rousselet & Cyril Pernet
% ------------------------------
%  Copyright (C) LIMO Team 2019


if nargin<5; alphav=.05;  end
if nargin<4; minnbchan=2; end

% find clusters in the observed data
[posclusterslabelmat,nposclusters] = limo_findcluster(ori_p<=alphav,channeighbstructmat,minnbchan);

% get bootstrap parameters
nboot                              = length(boot_maxclustersum);
sort_clustermax                    = sort(boot_maxclustersum);
n                                  = sum(isnan(sort_clustermax)); 

% NaN present if there was no clusters under H0 - just warp around
% should not happen if using limo_getclustersum as it returns 0 in that case
if n == nboot
    error('no cluster max value found - check boot_maxclustersum parameter')
else
    sort_clustermax(isnan(sort_clustermax))=[];
    sort_clustermax = [NaN(n,1); sort_clustermax];
end
maxclustersum_th = sort_clustermax(round((1-alphav)*nboot));
fprintf('cluster mass threshold: %g\n',maxclustersum_th)

% compute the mask: for each cluster do the sum and set significant if > maxclustersum_th
mask = zeros(size(ori_f));
cluster_label = 1; 
if nposclusters~=0
    for C = nposclusters:-1:1 % compute cluster sums & compare to bootstrap threshold
        maxval(C) = sum(ori_f(posclusterslabelmat==C));
        if  maxval(C)>= maxclustersum_th
            mask(posclusterslabelmat==C)= cluster_label; % flag clusters above threshold
            cluster_label = cluster_label+1;
        end
    end
end

% compute corrected p-values: number of times observed mass > bootstrap
mask2 = logical(mask); % logical - faster for masking
pval  = ones(size(mask));
if any(mask2(:)) % sum(mask2(:))>0
    L       = posclusterslabelmat.*mask2; % remove non significant clusters
    CL_list = setdiff(unique(L),0);       % remove label 0
    for CL=1:length(CL_list)
        % sort_clustermax
        cluster_mass = sum(ori_f(L==CL_list(CL)));
        if any(cluster_mass == maxval) % double checking this is in the mask
            p = 1-sum(cluster_mass>=sort_clustermax)./nboot;
            if p ==0
                p = 1/nboot; % never 0
            end
            tmp = ones(size(mask));
            tmp(L==CL_list(CL)) = p; % set p-values for many cells
            pval = pval.*tmp; % tmp is never at the same location so we can just 'add' values
        else
            error('cannot find the cluster mass for p-value? while found for the mask?? something is seriously wrong')
        end
    end
end
pval(pval==1)=NaN;

% check again corrected p values are correct
if any(pval > alphav)
   error('some corrected p-values are above the set alpha value, which should not happen - this is a bug, please contact the LIMO team') 
end

% just for output
if exist('maxval','var')
     maxval = max(maxval);   % biggest cluster
else
     maxval = 0;
end

