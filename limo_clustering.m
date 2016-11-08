function [mask,cluster_p] = limo_clustering(M,P,bootM,bootP,LIMO,MCC,p,fig)

% FORMAT:  [mask,cluster_p] = limo_clustering(M,P,bootM,bootP,LIMO,MCC,p)
%
% INPUT
% M = 2D matrix of observed F values (for a single channel use 1 in the 1st dimension)
% P = 2D matrix of observed p values (for a single channel use 1 in the 1st dimension)
% bootM = 3D matrix of F values for data bootstrapped under H0
% bootP = 3D matrix of F values for data bootstrapped under H0
% LIMO = LIMO structure - the necessary information is
%                         LIMO.data.chanlocs: the structure describing channels
%                         LIMO.data.neighbouring_matrix: the binary matrix of neighbours
% MCC = 2 (spatial-temporal clustering) or 3 (temporal clustering)
% p   = threshold to apply (note this applied to create the excursion set,
%       and then use as the corrected p value)
% fig = 1/0 to plot the maximum stat un der H0
%
% OUTPUT
% mask is a binary matrix of the same size as M corresponding to a threshold
%      p corrected for multiple comparisons
% cluster_p are the p-values obtained via the matrix bootM (corrected)
% max_th is the cluster mass threshold controlling the type 1 FWER
%
% Cyril Pernet v1
% outsourced from limo_stat_values
% --------------------------------
% Copyright (C) LIMO Team 2016

if nargin == 7
    fig = 1;
end

if size(M,1) == 1
    MCC = 3;
end
cluster_p = [];
mask = [];

%% spatial-temporal clustering

if MCC == 2 && size(bootM,1)>1
    nboot = size(bootM,3);
    U = round((1-p)*nboot); % bootstrap threshold
    minnbchan = 2;
    expected_chanlocs = LIMO.data.chanlocs;
    channeighbstructmat = LIMO.data.neighbouring_matrix;
    boot_maxclustersum=zeros(nboot,1); % compute bootstrap clusters
    disp('getting clusters under H0 boot ...');
    parfor boot = 1:nboot
        % boot_maxclustersum(boot) = limo_getclustersum(bootM(:,:,boot),bootP(:,:,boot),channeighbstructmat,minnbchan,p);
        [posclusterslabelmat,nposclusters] = limo_findcluster((bootP(:,:,boot) <= p),channeighbstructmat,minnbchan);
        bootM_b = bootM(:,:,boot);
        if nposclusters~=0
            tmp=zeros(1,nposclusters);
            for C = 1:nposclusters % compute sum for each cluster
                tmp(C) = sum( bootM_b(posclusterslabelmat==C) );
            end
            boot_maxclustersum(boot) = max(tmp(:)); % save max across clusters
        else
            boot_maxclustersum(boot) = 0;
        end
    end
    [mask, pval, maxval, maxclustersum_th] = limo_cluster_test(M,P,boot_maxclustersum,channeighbstructmat,minnbchan,p);
    
end

%% temporal clustering

if MCC == 2 && size(bootM,1)==1 || MCC == 3
    nboot = size(bootM,3);
    U = round((1-p)*nboot); % bootstrap threshold
    [th,boot_maxclustersum] = limo_ecluster_make(squeeze(bootM),squeeze(bootP),p);   
    maxclustersum_th = th.elec;
    [sigcluster, maxval,pval] = limo_ecluster_test(squeeze(M),squeeze(P),th,p);
    mask = sigcluster.elec_mask; 
end

%% plot
if fig == 1
    mass= sort(boot_maxclustersum);
    figure('Name','Correction by clustering: results under H0')
    plot(mass,'LineWidth',3); grid on; hold on; 
    
    plot(min(find(mass==maxclustersum_th)),maxclustersum_th,'r*','LineWidth',5)
    txt = 'threashold \rightarrow';
    text(min(find(mass==maxclustersum_th)),maxclustersum_th,txt,'FontSize',12,'HorizontalAlignment','right');
    
    [val,loc]=min(abs(mass-maxval)); 
    plot(loc,maxval,'r*','LineWidth',5)
    txt = sprintf('maximum');
    txt = [txt '\rightarrow'];
    text(loc,maxval,txt,'FontSize',12,'HorizontalAlignment','right');
    
    title('Cluster-mass Maxima under H0','FontSize',12)
    xlabel('sorted bootstrap iterations','FontSize',12); 
    ylabel('Freq.','FontSize',12)
    box on; set(gca,'Layer','Top')
end

