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
%
% v1 Guillaume Rousselet, University of Glasgow, August 2010
% edit Marianne Latinus adding spm_bwlabel
% Cyril Pernet Edited for return p values + compatible time frequency  
% ---------------------------------------------------------------------
%  Copyright (C) LIMO Team 2016
%
% See also limo_ecluster_make limo_tfcluster_make

if nargin < 4
    alpha_value = 0.05;
end

if nargin < 5
    boot_maxclustersum = [];
end

Ne = size(orif,1); % electrodes or frequencies
Nf = size(orif,2); % time frames

%% threshold the data base on the maximum cluster sum obtained over the 1st dimension
if isfield(th, 'max') 

    sigcluster.max_mask = zeros(Ne,Nf);
    ME = [];
    
    for E = 1:Ne % for each electrode or frequency
        if exist('spm_bwlabel','file') == 3
            [L,NUM] = spm_bwlabel(double(orip(E,:)<=alpha_value), 6); % find clusters
        elseif exist('bwlabel','file') == 2
            [L,NUM] = bwlabeln(orip(E,:)<=alpha_value);
        else
            errordlg('You need either the Image Processing Toolbox or SPM in your path'); return
        end
        
        maxval = zeros(1,NUM);
        for C = 1:NUM % for each cluster compute cluster sums & compare to bootstrap threshold
            maxval(C) = sum(abs(orif(E,L==C)));
            if  maxval(C) >= th.max;
                sigcluster.max_mask(E,L==C)=1; % flag clusters above threshold
            end
        end
        maxval = max(maxval);
    end 
end

%% threshold the data base on the maximum of cluster sum obtained over each electrode
if isfield(th, 'elec') 

    sigcluster.elec_mask = zeros(Ne,Nf);
    pval = nan(Ne,Nf);
    ME = [];
    try
        [L,NUM] = bwlabeln(orip<=alpha_value); % find clusters
    catch ME
        try
            [L,NUM] = spm_bwlabel(double(orip<=alpha_value), 6);
        catch ME
            errordlg('You need either the Image Processing Toolbox or SPM in your path'); return
        end
    end

    maxval = zeros(1,NUM);
    for C = 1:NUM % compute cluster sums & compare to bootstrap threshold
        maxval(C) = sum(abs(orif(L==C)));
        if maxval(C) >= th.elec
            sigcluster.elec_mask(L==C)=1; % flag clusters above threshold

            if ~isempty(boot_maxclustersum)
                p = 1 - sum(maxval(C) >= boot_maxclustersum)/length(boot_maxclustersum);            
                if p ==0
                    p = 1/length(boot_maxclustersum);
                end
                pval(L==C) = p;
            end
        end
    end
    maxval = max(maxval);
end




