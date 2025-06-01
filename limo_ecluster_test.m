function [sigcluster,pval,maxval] = limo_ecluster_test(orif,orip,th,alpha_value,boot_maxclustersum)
% function sigcluster = limo_ecluster_test(orif,orip,th,alpha_value)
%
% ECLUSTER_TEST computes sums of temporal clusters of significant F values and 
% compares them to a threshold obtained from ECLUSTER_MAKE.
% The function operates along the last dimension.
% For each electrode:
%       significant F values are clustered in time;
%       the sum of F values inside each cluster is computed;
%       this sum is compared to the threshold sum stored in TH.
% NOTE: for a two-tailed bootstrap t-test, enter a squared matrix of T values:
% th = ecluster_test(boott.^2,bootp,th,alpha_value)
%
% INPUTS:
%           ORIF = 2D or 1D matrix of F values with format electrode x
%               frames or 1 x frames;
%           ORIP = matrix of p values associated with orif, with same
%               format;
%           TH = output from ECLUSTER_MAKE
%           alpha_value = type I error rate, default 0.05
%           boot_maxclustersum = max cluster F value for each boot
%
% OUTPUTS:
%           SIGCLUSTER.ELEC = significant [0/1] data points at each
%               electrode, based on a cluster threshold defined
%               independently at each electrode;
%           SIGCLUSTER.MAX = significant [0/1] data points at each
%               electrode, based on a cluster threshold defined by 
%               taking the max cluster across electrodes; this is a
%               more conservative way to control for multiple comparisons 
%               than using a spatial-temporal clustering technique.
%           pval is the pvalue computed per channel cluster
%           maxval_channel is the maximum temporal cluster value
%
% Guillaume Rousselet, Marianne Latinus, Cyril Pernet 
% -------------------------------------------------------
%  Copyright (C) LIMO Team 2019
%
% See also limo_ecluster_make limo_tfcluster_make

if nargin < 4
    alpha_value = 0.05;
end

if nargin < 5
    boot_maxclustersum = [];
end

pval   = [];
maxval = [];
Ne     = size(orif,1); % electrodes or frequencies
Nf     = size(orif,2); % time frames

%% threshold the data based on the maximum cluster sum obtained over the 1st dimension
% even if temporal clustering is done initially, we still just use the max across the 1st dim
if isfield(th, 'max')

    sigcluster.max_mask = zeros(Ne,Nf);
    for E = 1:Ne % for each electrode or frequency
        if exist('spm_bwlabel','file') == 3
            [L,NUM] = spm_bwlabel(double(orip(E,:)<=alpha_value), 6); % find clusters
        elseif exist('bwlabel','file') == 2
            [L,NUM] = bwlabeln(orip(E,:)<=alpha_value);
        else
            errordlg('You need either the Image Processing Toolbox or SPM in your path'); return
        end

        maxval = zeros(1,NUM);
        cluster_label = 1;
        for C = 1:NUM % for each cluster compute cluster sums & compare to bootstrap threshold
            maxval(C) = sum(abs(orif(E,L==C)));
            if  maxval(C) >= th.max
                sigcluster.max_mask(E,L==C)= cluster_label; % flag clusters above threshold
                cluster_label =  cluster_label + 1;
            end
        end
        maxval= max(maxval);
    end
end

%% threshold the data based on the maximum of cluster sum obtained over each electrode
if isfield(th, 'elec') 

    sigcluster.elec_mask = zeros(Ne,Nf);
    maxval_channel       = zeros(Ne,1);
    pval                 = nan(Ne,Nf);
    cluster_label        = 1;
    
    % we still to control over the entire space, so the max across channels
    if size(boot_maxclustersum,1) > 1
        boot_maxclustersum = max(boot_maxclustersum,[],1);
    end
    
    % but cluster one channel at a time
    for channel = 1:Ne
        if exist('bwlabeln','file')~=0
            [L,NUM] = bwlabeln(orip(channel,:)<=alpha_value); % find clusters
        else
            errordlg('You need either the Image Processing Toolbox or SPM in your path'); return
        end

        maxval = zeros(1,NUM);
        for C = 1:NUM % compute cluster sums & compare to bootstrap threshold
            maxval(C) = sum(abs(orif(channel,L==C)));
            if maxval(C) >= th.elec(channel)
                sigcluster.elec_mask(channel,L==C)= cluster_label; % flag clusters above threshold
                cluster_label =  cluster_label+1;
                p = 1 - sum(maxval(C) >= boot_maxclustersum)./length(boot_maxclustersum);
                if p==0
                    p = 1/length(boot_maxclustersum);
                end
                pval(channel,L==C) = p;
            end
        end
        if ~isempty(maxval)
            maxval_channel(channel) = max(maxval);
        else
            maxval_channel(channel) = 0;
        end
    end
    maxval = maxval_channel;
end



