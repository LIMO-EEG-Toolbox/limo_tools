function [t,p,dfe]=limo_ttest_permute(Data,n_perm)

% pseudo one sample t-test using sign-test
%
% FORMAT: [t,p,dfe]=limo_ttest_permute(Data,n_perm)
%
% INPUTS:
% data = a matrix of data to be used in the one sample t-test
% n_perm  = number of permutations to do.
%
% OUTPUTS:
% [t,p,dfe] = t and p values under H0, dfe
%
% compute that the mean of the data is not different from 0 (H0)
% using permutation of the sign of the data - if not different from 0 then
% splitting the data as positive or negative should still give 0
%
% taken from mult_comp_perm_t1.m by David Groppe
%
% Cyril Pernet 10-01-2013
% GAR 18-02-13: embedded electrode loop inside permutation loop so that the
%               same permutation is applied to every electrode
% -----------------------------
% Copyright (C) LIMO Team 2015


%% check inputs
if numel(size(Data)) == 3
    [n_channels,n_var,n_obs]=size(Data);
else
    n_channels = 1;
    [n_var,n_obs]=size(Data);
    temp = Data; Data = NaN(1,n_var,n_obs);
    Data(1,:,:) = temp; clear temp
end

if n_obs<7,
    n_psbl_prms=2^n_obs;
    fprintf(' Due to the very limited number of observations, \n the total number of possible permutations is small. \n Thus only a limited number of p-values (at most %d) are possible \n and the test might be overly conservative.\n',n_psbl_prms);
end


%% Set up permutation test
if n_obs<=12,
    n_perm=2^n_obs; % total number of possible permutations
    exact=1;
    fprintf('Due to the limited number of observations, all possible permutations of the data will be computed instead of random permutations.\n');
else
    exact=0;
    if nargin == 1
        n_perm = 1000;
    end
end

fprintf('Executing permutation test with %d permutations...\n',n_perm);

%% Compute permutations

t = NaN(n_channels,n_var,n_perm); %p = t;
% Constant factor for computing t
sqrt_nXnM1=sqrt(n_obs*(n_obs-1));
dfe = n_obs-1;

if exact %Use all possible permutations
    
    for perm=1:n_perm
        % set sign of each trial / participant's data
        temp=perm-1;
        n_temp=length(temp);
        sn=-ones(n_obs,1);
        sn(1:n_temp,1)=2*temp-1; % -1/1
        sn_mtrx=repmat(sn,1,n_var);
        
        for c = 1:n_channels
            data = squeeze(Data(c,:,:));
            d_perm=data.*sn_mtrx';
            
            % computes t-score of permuted data across all channels and time points
            sm=sum(d_perm,2);
            mn=sm/n_obs;
            sm_sqrs=sum(d_perm.^2,2)-(sm.^2)/n_obs;
            stder=sqrt(sm_sqrs)/sqrt_nXnM1;
            t(c,:,perm)=mn./stder;
        end
    end
    
else %Use random permutations
    
    for perm=1:n_perm
        % randomly set sign of each trial / participant's data
        sn=(rand(n_obs,1)>.5)*2-1;
        sn_mtrx=repmat(sn,1,n_var);
        
        for c = 1:n_channels
            data = squeeze(Data(c,:,:));
            d_perm=data.*sn_mtrx';
            
            % computes t-score of permuted data across all channels and time points
            sm=sum(d_perm,2);
            mn=sm/n_obs;
            sm_sqrs=sum(d_perm.^2,2)-(sm.^2)/n_obs;
            stder=sqrt(sm_sqrs)/sqrt_nXnM1;
            t(c,:,perm)=mn./stder;
        end
    end
    
end

%     p(c,:,:)= 2*tcdf(-abs(t(c,:,:)), dfe);
p=2*tcdf(-abs(t), dfe);

%% thresholding
% for univariate stats t obs > 95 percentile of the t distribution computed
% here // in mult_comp_perm_t1.m, one takes the max of each permutation and
% threshold all the t obs based on the max distribution

