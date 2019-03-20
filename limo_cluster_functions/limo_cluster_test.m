function [mask, pval, maxval, maxclustersum_th] = limo_cluster_test(ori_f,ori_p,boot_maxclustersum,channeighbstructmat,minnbchan,alphav)

% limo_cluster_test finds clusters of significant F values, computes the
% sum of F values inside each cluster, and compares that sum to a threshold
% sum of F values expected by chance.
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
% See also limo~_getcluster_test limo_getclustersum
% 
% Guillaume Rousselet, University of Glasgow, June 2010
% optional L & NUM outputs: GAR, Feb 2012
% optional pval & maxclustersum_th outputs: GAR, Feb 2012
% Cyril Pernet changed pval to be a map with NaN or the cluster p value May 2013
% added a warping of NaN Mars 2014 
% -------------------------------------------------
%  Copyright (C) LIMO Team 2016


if nargin<5; alphav=.05;  end
if nargin<4; minnbchan=2; end

[posclusterslabelmat,nposclusters] = limo_findcluster(ori_p<=alphav,channeighbstructmat,minnbchan);
nboot           = length(boot_maxclustersum);
sort_clustermax = sort(boot_maxclustersum);
n               = sum(isnan(sort_clustermax)); 

% NaN present if there was no clusters under H0 - just warp around
% should not happen if using limo_getclustersum as it returns 0 in that case
if n~=0 
    sort_clustermax(isnan(sort_clustermax))=[];
    sort_clustermax = [NaN(n,1); sort_clustermax];
end
maxclustersum_th = sort_clustermax(round((1-alphav)*nboot));

mask = zeros(size(ori_f));
if nposclusters~=0
    for C = 1:nposclusters % compute cluster sums & compare to bootstrap threshold
        maxval(C) = sum(ori_f(posclusterslabelmat==C));
        if  maxval(C)>= maxclustersum_th;
            mask(posclusterslabelmat==C)=1; % flag clusters above threshold
        end
    end
end
maxval = max(maxval);
mask = logical(mask); 

if sum(mask(:))>0
    L=posclusterslabelmat.*mask;  % figure; imagesc(L)
    CL_list=setdiff(unique(L),0);
    pval = NaN(size(L));
    for CL=1:length(CL_list)
        L(L==CL_list(CL))=CL;
        p=1-(sum(sum(ori_f(L==CL))>=boot_maxclustersum)./nboot);
        if p ==0
            p = 1/nboot;
        end
        pval(L==CL_list(CL)) = p;
    end
else
    pval=NaN;
end



