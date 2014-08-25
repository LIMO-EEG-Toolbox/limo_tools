function sigcluster = limo_tfcluster_test(orif,orip,th,alphav)

% limo_tfcluster_test computes sums of time*frequency clusters of 
% significant F values and compares them to a threshold obtained 
% from limo_tfcluster_make. The function operates along the last dimension.
% For each electrode:
% (1) significant F values are clustered in time;
% (2) the sum of F values inside each cluster is computed;
% (3) this sum is compared to the threshold obtained from tf distribution(s)
%
% FORMAT: sigcluster = limo_tfcluster_test(orif,orip,th,alpha)
%
% INPUTS: orif = 3D (electrode * freq * time) matrix of observed F values 
%         orip = 3D matrix of observed p values 
%         th = output from limo_ecluster_make (sorted distribution of max under H0)
%         alpha = type I error rate, default 0.05
%
% OUTPUTS: sigcluster.max_mask significant [0/1] data points 
%                              based on a cluster threshold defined by 
%                              taking the max cluster for the entire space
%          sigcluster.max_pvalues corrected p_values
%          sigcluster.elec_mask significant [0/1] data points for each
%                                element of the 1st dimension (should be electrodes)
%                                - the correction is computed independently for each
%                                of these elements (ie each time*freq map is corrected
%                                independently) 
%         sigcluster.elec_mask corrected p_values
%
% See also limo_tfcluster_make 
% -------------------------------------------------------------------
%  Copyright (C) LIMO Team 2014

% Cyril Pernet v1 June 2014 (adapted from limo_ecluster_test)

%% Inputs
if ndims(orif) ~= ndims(orip)
    error('dimension issues between observed F and p values');
end

Ne = size(orif,1); % electrodes 
Nf = size(orif,2); % frequencies 
Nt = size(orif,3); % time

if nargin < 3
    alphav = 0.05;
end
    
%% threshold the data base on the maximum cluster sum obtained over the whole space
if isfield(th, 'max') 

    sigcluster.max_mask = zeros(Ne,Nf,Nt);
    sigcluster.max_pvalues = NaN(Ne,Nf,Nt);
    b = size(th.max,2);
    U = round((1-alphav)*b);
    
    for e = 1:Ne % for each electrode or frequency
        if exist('bwlabel','file') ~=0
            [L,NUM] = bwlabeln(squeeze(orip(e,:,:))<=alphav); % find clusters
        elseif exist('spm_bwlabel','file') ~=0
            [L,NUM] = spm_bwlabel(double(squeeze(orip(e,:,:))<=alphav), 6);
        else
            errordlg('You need either the Image Processing Toolbox or SPM in your path'); return
        end
        
        for C = 1:NUM % for each cluster compute cluster sums & compare to bootstrap threshold
            tmp = squeeze(orif(e,:,:));
            corrected_p = NaN(Nf,Nt);
            corrected_p((L==C)) = 1-sum(sum(tmp(L==C))>=th.max)./b;
            sigcluster.max_pvalues(e,:,:) = corrected_p;
            if sum(tmp(L==C)) >= th.max(U)
                mask = zeros(Nf,Nt);
                mask(L==C)=1; % flag clusters above threshold
                sigcluster.max_mask(e,:,:) = mask;
            end
        end
        
    end 
end

%% threshold the data base on the maximum of cluster sum obtained over each electrode
if isfield(th, 'elec') && ndims(orif) == 3

    sigcluster.elec_mask = zeros(Ne,Nf,Nt);
    sigcluster.elec_pvalues = NaN(Ne,Nf,Nt);
    b = size(th.elec,2);
    U = round((1-alphav)*b);
    
    for e=1:Ne
        if exist('bwlabel','file') ~=0
            [L,NUM] = bwlabeln(squeeze(orip(e,:,:))<=alphav); % find clusters
        elseif exist('spm_bwlabel','file') ~=0
            [L,NUM] = spm_bwlabel(double(squeeze(orip(e,:,:))<=alphav), 6);
        else
            errordlg('You need either the Image Processing Toolbox or SPM in your path'); return
        end
        
        for C = 1:NUM % compute cluster sums & compare to bootstrap threshold
            tmp = squeeze(orif(e,:,:));
            corrected_p = NaN(Nf,Nt);
            corrected_p((L==C)) = 1-sum(sum(tmp(L==C))>=th.elec(e,:))./b;
            sigcluster.elec_pvalues(e,:,:) = corrected_p;
            if sum(tmp(L==C)) >= th.elec(e,U)
                mask = zeros(Nf,Nt);
                mask(L==C)=1; % flag clusters above threshold
                sigcluster.elec_mask(e,:,:) = mask;
            end
        end
    end
end




