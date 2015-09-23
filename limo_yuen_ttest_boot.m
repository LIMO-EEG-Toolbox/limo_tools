function [Ty,diff,CI,p,tcrit,df,ori_mask,boot_mask,cluster_mask]=limo_yuen_ttest_boot(a,b,percent,alpha,nboot,minnbchan,channeighbstructmat)

% function [Ty,diff,CI,p,tcrit,df,ori_mask,boot_mask,cluster_mask]=limo_yuen_ttest_boot(a,b,percent,alpha,nboot,minnbchan,channeighbstructmat)
%
% Computes Ty (Yuen's T statistic) to compare the trimmed means of two
% INDEPENDENT groups (single-trials, or subjects with the same number of electrodes).
% Ty=(tma-tmb) / sqrt(da+db), where tma & tmb are trimmed means of a & b,
% and da & db are yuen's estimate of the standard errors of tma & tmb.
% Ty is distributed approximately as Student's t with estimated degrees of freedom, df.
% The p-value (p) is the probability of obtaining a t-value whose absolute value
% is greater than Ty if the null hypothesis (i.e., tma-tmb = mu) is true.
% In other words, p is the p-value for a two-tailed test of H0: tma-tmb=0;
%
% INPUTS:
%
% a & b are matrices EEGLAB format with trials as the third dimension; 
%       percent must be a number between 0 & 100.
% minnbchan = minimum number of channels, default = 2, see
%       LIMO_FT_FINDCLUSTER for details
% channeighbstructmat = output of low-level function
%       limo_ft_neighbourselection or higher-level function limo_get_channeighbstructmat
% 
%   Default values:
%   percent = 20;
%   alpha = 0.05;
%   nboot = 1000;
%   minnbchan = 2;
%
% OUTPUTS:
%
% Ty:           Yuen T statistics
% diff:         difference between trimmed means of a and b
% CI:           confidence interval around the difference
% p:            p value
% tcrit:        1-alpha/2 quantile of the Student's t distribution with adjusted
%               degrees of freedom
% df:           degrees of freedom
% ori_mask:     logical mask of significant data points based on the original Yuen t-test
% boot_mask:    logical mask of significant data points based on the bootstrapped T values
% cluster_mask: logical mask of significant data points based on the
%               clustered bootstrapped and original T values
%
% See Wilcox (2005), Introduction to Robust Estimation and Hypothesis
% Testing (2nd Edition), page 159-161 for a description of the Yuen
% procedure.
% You can also check David C. Howell, Statistical Methods for Psychology,
% sixth edition, p.43, 323, 362-363.
%
% Original Yuen t-test code by Prof. Patrick J. Bennett, McMaster University
% Various editing, GAR, University of Glasgow, Dec 2007
% Vector version for EEG: GAR, University of Glasgow, Nov 2009
% CI output: GAR, University of Glasgow, June 2010
% Bootstrap, cluster, LIMO version: GAR, University of Glasgow, July 2010
% Debugging of various output variables: Luisa Frei, University of Glasgow, October 2011
% Updated calls to limo_yuen_ttest: GAR, Sept 2015 
%
% -----------------------------
%  Copyright (C) LIMO Team 2010
%
% See also LIMO_YUEN_TTEST LIMO_GETCLUSTERSUM LIMO_CLUSTER_TEST

if isempty(a) || isempty(b) 
    error('limo_boot_yuen_ttest:InvalidInput', 'data vectors cannot have length=0');
end

if nargin<5;nboot=1000;end
if nargin<4;alpha=.05;end
if nargin<3;percent=20;end

if (percent >= 100) || (percent < 0)
    error('limo_boot_yuen_ttest:InvalidPercent', 'PERCENT must be between 0 and 100.');
end

if nargin==7 % check cluster function
    if nargin<6;minnbchan=2;end
    if size(channeighbstructmat,1) ~= size(a,1)
        error('limo_boot_yuen_ttest:InvalidInput', 'channeighbstructmat must have the same number of electrodes as your data\n');
    end
end



U = round((1-alpha)*nboot); % univariate bootstrap threshold
lo = round(nboot.*(alpha./2));
hi = nboot-lo;

%% t-tests on original data

fprintf('limo_boot_yuen_ttest: t-tests on original data...\n')
[ori_Ty,ori_diff,se,ori_CI,ori_p,ori_tcrit,ori_df]=limo_yuen_ttest(a,b,percent,alpha);
ori_mask = ori_p<=alpha;

%% centre data

ca=a-repmat(limo_trimmed_mean(a,percent),[1 1 size(a,3)]);
cb=b-repmat(limo_trimmed_mean(b,percent),[1 1 size(b,3)]);

%% bootstrap

boot_index_a =  ceil(rand(size(a,3),nboot)*size(a,3)); % sample 3rd dimension WITH replacement
boot_index_b =  ceil(rand(size(b,3),nboot)*size(b,3));

boot_t = zeros(size(a,1),size(a,2),nboot);
boot_p = boot_t;

for bb = 1:nboot % bootstrap loop
    
    fprintf('limo_boot_yuen_ttest: bootstrap %i / %i...\n',bb,nboot)
    [boot_t(:,:,bb),diff,se,CI,boot_p(:,:,bb),tcrit,df] = ...
        limo_yuen_ttest(ca(:,:,boot_index_a(:,bb)),cb(:,:,boot_index_b(:,bb)),percent,alpha);
    
end % bootstrap loop

%% bootstrap univariate test

sort_boot_t = sort(boot_t,3);
boot_mask = ori_Ty < sort_boot_t(:,:,lo+1) | ori_Ty > sort_boot_t(:,:,hi);

%% bootstrap cluster test

if nargin==7 % cluster F values
    
    boot_maxclustersum = zeros(nboot,1);
    for bb = 1:nboot % get cluster for each bootstrap
        boot_maxclustersum(bb) = limo_getclustersum(squeeze(boot_t(:,:,bb)).^2,squeeze(boot_p(:,:,bb)),channeighbstructmat,minnbchan,alpha);
    end

    sort_boot_maxclustersum = sort(boot_maxclustersum,1);
    boot_maxclustersum_th = sort_boot_maxclustersum(U);

    % cluster test on original data
    
    cluster_mask = logical(limo_cluster_test(ori_Ty.^2,ori_p,boot_maxclustersum_th,channeighbstructmat,minnbchan,alpha));
               
end

Ty=ori_Ty;
p=ori_p;
diff=ori_diff;
tcrit=ori_tcrit;
CI=ori_CI;
df=ori_df;


end




