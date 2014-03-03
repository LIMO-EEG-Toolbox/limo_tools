function maxclustersum = limo_getclustersum(f,p,channeighbstructmat,minnbchan,alpha)

%function maxclustersum = limo_getclustersum(f,p,channeighbstructmat,minnbchan,alpha)
%
% LIMO_GETCLUSTERSUM finds clusters of significant F values, computes the
% sum of each cluster, and return the maximum sum across all clusters.
%
% INPUTS:
%
% f = 2D/3D matrix of F values - so if using T values, enter T.^2
% p = 2D/3D matrix of P values
% channeighbstructmat = output of limo_ft_neighbourselection
% minnbchan = minimum number of channels, default = 2, see
%   LIMO_FT_FINDCLUSTER for details
% alpha level, default 0.05
%
% OUTPUTS:
%
% maxclustersum = max cluster sum across all clusters of significant F
%   values
%
% See also LIMO_CLUSTER_TEST LIMO_FT_FINDCLUSTER
% v1: GAR, University of Glasgow, June 2010
% -----------------------------
%  Copyright (C) LIMO Team 2010


if nargin<5;alpha=.05;end
if nargin<4;minnbchan=2;end

[posclusterslabelmat,nposclusters] = limo_ft_findcluster(p<=alpha,channeighbstructmat,minnbchan);
% figure;imagesc(posclusterslabelmat)

if nposclusters~=0
    
    tmp=zeros(1,nposclusters);
    
    for C = 1:nposclusters % compute sum for each cluster
        
        tmp(C) = sum( f(posclusterslabelmat==C) );
        
    end

    maxclustersum = max(tmp(:)); % save max across clusters

else

    maxclustersum = 0;

end

end




