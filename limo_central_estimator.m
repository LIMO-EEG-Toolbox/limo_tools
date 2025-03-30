function [est,HDI,bb] = limo_central_estimator(Y,estimator,prob_coverage,legacy_mode)

% Compute a data estimator and its highest density intervals (HDI) based 
% on bayesian bootstrap estimates.
%
% FORMAT [est,HDI,bb] = limo_central_estimator(Y,'estimator',prob_coverage,legacy_mode);
%
% INPUTS Y is a 2D matrix, e.g. frames x participants or trials
%        estimator is 'Mean', 'Trimmed mean', (default) 'HD' (Mid-decile Harell-Davis) or 'Median'
%        prob_coverage is the probability coverage- default 0.95%
%        legacy_mode is 'false' (default) to use gamma and quantiles
%               if 'true' it uses exponential smapling and HDI 
%
% OUTPUT est is the estimator
%        ci is the quantile interval or highest density interval
%        bb is a vector of (Bayes) bootstraped estimators 
%
% Bayesian bootstrap implementation based on orignal R code from Rasmus Baath:
% http://www.sumsar.net/blog/2015/07/easy-bayesian-bootstrap-in-r/
% HDI implementation based on original R code HDIofMCMC from John K. Kruschke:
% https://github.com/boboppie/kruschke-doing_bayesian_data_analysis/blob/master/1e/HDIofMCMC.R
% The default follows Etz A, Chávez de la Peña AF, Baroja L, Medriano K, Vandekerckhove J. 
% The HDI + ROPE decision rule is logically incoherent but we can fix it. 
% Psychol Methods. 2024 May 23. doi: 10.1037/met0000660. https://pubmed.ncbi.nlm.nih.gov/38780591/
%
% see also mean, median, limo_trimmed_mean, limo_harrell_davis
%
% Guillaume Rousselet, Cedric Cannard & Cyril Pernet February 2025
% ---------------------------------------------------------------
%  Copyright (C) LIMO Team 2025

Nb = 1000; % number of bootstrap samples 

if nargin == 3
    legacy_mode = false; 
elseif nargin == 2
    prob_coverage = 0.95;
    legacy_mode = false; 
elseif nargin == 1
    prob_coverage = 0.95;
    estimator = 'Trimmed mean';
    legacy_mode = false; 
end
alphav = (1-prob_coverage)/2;

% compute the estimator
if strcmpi(estimator,'Mean')
    est = mean(Y,2);
elseif strcmpi(estimator,'Trimmed mean')
    est = limo_trimmed_mean(Y,20); % default 20% trimmed mean
elseif strcmpi(estimator,'HD')
    est = limo_harrell_davis(Y,.5); % default to estimation of the 5th decile
elseif strcmpi(estimator,'Median')
    est = median(Y,2);
end

% sample with replacement from Dirichlet
% sampling = number of observations, e.g. participants
n = size(Y,2); 
bb = zeros(size(Y,1),Nb);
parfor boot=1:Nb % bootstrap loop
    if legacy_mode
        theta    = exprnd(1,[n,1]);
    else
        theta    = gamrnd(ones(n,1),1); 
    end
    weigths  = theta ./ repmat(sum(theta,1),n,1); % normalize to sum to 1 to align with Dirichlet distribution properties
    resample = (datasample(Y',n,'Replace',true,'Weights',weigths))';
    
    % compute the estimator
    if strcmpi(estimator,'Mean')
        bb(:,boot) = mean(resample,2);
    elseif strcmpi(estimator,'Trimmed Mean')
        bb(:,boot) = limo_trimmed_mean(resample,20); 
    elseif strcmpi(estimator,'HD')
        bb(:,boot) = limo_harrell_davis(resample,.5); 
    elseif strcmpi(estimator,'Median')
        bb(:,boot) = median(resample,2);
    end
end

sorted_data   = sort(bb,2); % sort bootstrap estimates
if legacy_mode
    upper_centile = floor(prob_coverage*size(sorted_data,2)); % upper bound
    nCIs          = size(sorted_data,2) - upper_centile;
    HDI           = zeros(2,size(Y,1));

    % for frame = 1:size(Y,1)
    %     tmp = sorted_data(frame,:);
    %     ci = 1:nCIs; ciWidth = tmp(ci+upper_centile) - tmp(ci); % all centile distances
    %     [~,index]=find(ciWidth == min(ciWidth)); % densest centile
    %     if length(index) > 1; index = index(1); end % many similar values
    %     HDI(1,frame) = tmp(index);
    %     HDI(2,frame) = tmp(index+upper_centile);
    % end

    % vectorized version of the loop above
    ci       = 1:nCIs;
    ciWidth  = sorted_data(:,ci+upper_centile) - sorted_data(:,ci); % all centile distances
    [~,J]    = min(ciWidth,[],2);
    r        = size(sorted_data,1);
    I        = (1:r)';
    index    = I+r.*(J-1); % linear index
    HDI(1,:) = sorted_data(index);
    index    = I+r.*(J+upper_centile-1); % linear index
    HDI(2,:) = sorted_data(index);
else
    % lim = quantile(oneframe, [alphav, 1 - alphav]);
    HDI(:, 1) = bb(:,round(alphav*Nb));
    HDI(:, 2) = bb(:,round((1-alphav)*Nb));
end

