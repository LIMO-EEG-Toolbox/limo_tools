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
if n~=0 
    sort_clustermax(isnan(sort_clustermax))=[];
    sort_clustermax = [NaN(n,1); sort_clustermax];
end
maxclustersum_th = sort_clustermax(round((1-alphav)*nboot));
fprintf('cluster mass threshold: %g\n',maxclustersum_th)

% compute the mask: for each cluster do the sum and set significant if > maxclustersum_th
mask = zeros(size(ori_f));
cluster_label = 1; 
if nposclusters~=0
    for C = 1:nposclusters % compute cluster sums & compare to bootstrap threshold
        maxval(C) = sum(ori_f(posclusterslabelmat==C));
        if  maxval(C)>= maxclustersum_th
            mask(posclusterslabelmat==C)= cluster_label; % flag clusters above threshold
            cluster_label = cluster_label+1;
        end
    end
end

if exist('maxval','var')
     maxval = max(maxval);   % biggest cluster
else
     maxval = 0;
end
mask2  = logical(mask); % logical - faster for masking

% compute corrected p-values: number of times observed mass > bootstrap
pval = NaN(size(mask));
if any(mask2(:))
    L       = posclusterslabelmat(mask2); % remove non significant clusters
    CL_list = setdiff(unique(L),0);
    for CL=1:length(CL_list)
        p = 1-(sum(sum(ori_f(L==CL_list(CL)))>=boot_maxclustersum)./nboot);
        if p ==0
            p = 1/nboot; % never 0
        end
        pval(L==CL_list(CL)) = p;
    end
end


