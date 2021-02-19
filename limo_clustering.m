function [mask,cluster_pval,max_th] = limo_clustering(varargin)

% FORMAT:  [mask,cluster_p,max_th] = limo_clustering(M,P,bootM,bootP,LIMO,MCC,p,fig)
%
% INPUT
% M     = 2D matrix of observed F values (for a single channel use 1 in the 1st dimension)
% P     = 2D matrix of observed p values (for a single channel use 1 in the 1st dimension)
% bootM = 3D matrix of F values for data bootstrapped under H0
% bootP = 3D matrix of F values for data bootstrapped under H0
% LIMO  = LIMO structure - the necessary information is
%                         LIMO.data.chanlocs: the structure describing channels
%                         LIMO.data.neighbouring_matrix: the binary matrix of neighbours
% MCC   = 2 (spatial-temporal clustering) or 3 (temporal clustering)
% p     = threshold to apply (note this applied to create the excursion set,
%         and then use as the corrected p value)
% fig   = 1/0 to plot the maximum stat under H0 (if empty, turned on if no
%         significant values)
%
% OUTPUT
% mask is a labelled matrix of the same size as M corresponding to a threshold
%      p corrected for multiple comparisons
% cluster_p are the p-values obtained via the matrix bootM (corrected)
% max_th is the cluster mass threshold controlling the type 1 FWER
%
% Cyril Pernet 
% outsourced from limo_stat_values
% ------------------------------
%  Copyright (C) LIMO Team 2019

% check inputs 
M       = varargin{1};
P       = varargin{2};
bootM   = varargin{3};
bootP   = varargin{4};
LIMO    = varargin{5};
MCC     = varargin{6};
p       = varargin{7};
if nargin == 7
    fig = [];
else
    fig = varargin{8};
end
clear varargin

% switch behavioural to 1D clustering if one channel, no matter user choice
if size(M,1) == 1
    MCC = 3; 
end

% boostrap parameters
nboot = size(bootM,3);      % nb of boostrap performed
    
% set outputs empty as default
cluster_pval = [];
mask         = [];
max_th       = [];

%% spatial-temporal clustering

if MCC == 2 && size(bootM,1)>1
    minnbchan           = 2;
    if isfield(LIMO,'data')
        channeighbstructmat = LIMO.data.neighbouring_matrix;
    else
        channeighbstructmat = LIMO.channeighbstructmat;
    end
    boot_maxclustersum  = zeros(nboot,1);     % maximum cluster mass at each bootstrap
    
    disp('getting clusters under H0 boot ...');
    parfor boot = 1:nboot
        % 1st find the cluster, thresholding H0 pvalues <= threshold p
        [posclusterslabelmat,nposclusters] = limo_findcluster((bootP(:,:,boot) <= p),channeighbstructmat,minnbchan);
        
        % 2nd compute the mass for each cluster
        bootM_b = bootM(:,:,boot);
        if nposclusters~=0
            tmp = zeros(1,nposclusters);
            for C = 1:nposclusters 
                tmp(C) = sum(bootM_b(posclusterslabelmat==C)); % sum stat value in a cluster label
            end
            boot_maxclustersum(boot) = max(tmp(:)); % save max value only
        else
            boot_maxclustersum(boot) = 0;
        end
    end
    
    % 3rd threshold observed cluster mass by the distribution of cluster
    % max computed in step 2
    [mask, cluster_pval, maxval, max_th] = limo_cluster_test(M,P,...
        boot_maxclustersum,channeighbstructmat,minnbchan,p);
end

%% temporal clustering

if MCC == 2 && size(bootM,1)==1 || MCC == 3
    % 1st get the distribution of maxima under H0
    [th,boot_maxclustersum]           = limo_ecluster_make(squeeze(bootM),squeeze(bootP),p);   
    max_th                            = th.elec;
    % 2nd threshold observed data
    [sigcluster, cluster_pval,maxval] = limo_ecluster_test(squeeze(M),squeeze(P),th,p, boot_maxclustersum);
    mask                              = sigcluster.elec_mask; 
end

%% plot
% when nothing is significant, always show why
if sum(mask(:)) == 0 && isempty(fig)
    fig = 1 ;
end

if fig == 1 
    figure('Name','Correction by clustering: results under H0')
    mass = sort(boot_maxclustersum);
    plot(mass,'LineWidth',3); grid on; hold on; 
    
    plot(min(find(mass==max_th)),max_th,'r*','LineWidth',5)
    txt = ['bootstrap threashold ' num2str(max_th) '\rightarrow'];
    text(min(find(mass==max_th)),max_th,txt,'FontSize',10,'HorizontalAlignment','right');
    
    [val,loc]=min(abs(mass-maxval)); 
    plot(loc,maxval,'r*','LineWidth',5)
    txt = ['biggest observed cluster mass: ' num2str(maxval) '\rightarrow'];
    text(loc,double(maxval),txt,'FontSize',10,'HorizontalAlignment','right');
    
    title('Cluster-mass Maxima under H0','FontSize',12)
    xlabel('sorted bootstrap iterations','FontSize',12); 
    ylabel('Freq.','FontSize',12)
    box on; set(gca,'Layer','Top')
end

