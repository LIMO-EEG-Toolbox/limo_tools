function [mask L NUM] = limo_cluster_test(ori_f,ori_p,maxclustersum_th,channeighbstructmat,minnbchan,alpha)

%function [mask L NUM] = limo_cluster_test(ori_f,ori_p,maxclustersum_th,channeighbstructmat,minnbchan,alpha)
%
% LIMO_CLUSTER_TEST finds clusters of significant F values, computes the
% sum of F values inside each cluster, and compares that sum to a threshold
% sum of F values expected by chance.
%
% INPUTS:
%
% ori_f: 2D matrix of F values from an ANOVA or t-test - so if using T
%   values, enter T.^2
%
% ori_p: 2D matrix of P values from the same ANOVA or t-test
%
% maxclustersum_th = threshold sum of F values expected by chance
%
% channeighbstructmat = output of limo_ft_neighbourselection
%
% minnbchan = minimum number of channels, default = 2, see
%   LIMO_FT_FINDCLUSTER for details
%
% alpha level, default 0.05
%
% OUTPUTS:
%
% mask = 2D logical matrix showing electrodes and time points with
%   significant effects corrected for multiple comparisons by a cluster test
% Use:
% mask = limo_cluster_test(ori_f,ori_p,maxclustersum_th,channeighbstructmat,minnbchan,alpha)
%
% OPTIONAL OUTPUTS:
%
% L = 2D map of significant clusters identified by a number from 1 to NUM
%
% NUM = number of significant clusters
%
% See also LIMO_GETCLUSTERSUM

% -----------------------------
%  Copyright (C) LIMO Team 2010

% v1: GAR, University of Glasgow, June 2010
% optional L & NUM outputs: GAR, Feb 2012

if nargin<5;alpha=.05;end
if nargin<4;minnbchan=2;end

[posclusterslabelmat,nposclusters] = limo_ft_findcluster(ori_p<=alpha,channeighbstructmat,minnbchan);

mask = zeros(size(ori_f));

if nposclusters~=0

    for C = 1:nposclusters % compute cluster sums & compare to bootstrap threshold

        if sum(ori_f(posclusterslabelmat==C)) >= maxclustersum_th;
            mask(posclusterslabelmat==C)=1; % flag clusters above threshold
        end

    end

end

mask = logical(mask);

if nargout>1
    if sum(mask(:))>0
        L=posclusterslabelmat.*mask; % figure;imagesc(L)
        CL_list=setdiff(unique(L),0);
        NUM=length(CL_list);
        for CL=1:length(CL_list)
            L(L==CL_list(CL))=CL;
        end
    else
        NUM=0;
        L=zeros(size(mask));
    end
end

end


