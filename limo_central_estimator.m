function [est,HDI,bb] = limo_central_estimator(Y,estimator,prob_coverage)

% Compute highest density intervals (HDI) based on bayesian bootstrap
% estimates.
%
% FORMAT [est,HDI,bb] = limo_central_estimator(Y,estimator,prob_coverage);
%
% INPUTS:
%       Y is a 2D matrix, e.g. time frames x participants
%       estimator is 'mean', 'tm', 'hd' or 'median' - default if 20% trimmed mean (tm)
%       prob_coverage is the probability coverage - default 0.95
%
% OUTPUT est is the estimator
%        ci is the high density interval
%        bb is a vector of (Bayes) bootstraped estimators 
%
% Guillaume Rousselet & Cyril Pernet January 2016
%
% Bayesian bootstrap implementation based on orignal R code from Rasmus Baath:
% http://www.sumsar.net/blog/2015/07/easy-bayesian-bootstrap-in-r/
% HDI implementation based on original R code HDIofMCMC from John K. Kruschke:
% https://github.com/boboppie/kruschke-doing_bayesian_data_analysis/blob/master/1e/HDIofMCMC.R
% ------------------------------------------
% Copyright (C) LIMO Team 2016

Nb = 1000; % number of bootstrap samples 

if nargin == 2
    prob_coverage = 0.95;
elseif nargin == 1
    prob_coverage = 0.95;
    estimator = 'tm'; % trimmed mean
end

% compute the estimator
if strcmpi(estimator,'mean')
    est = mean(Y,2);
elseif strcmpi(estimator,'tm')
    est = limo_trimmed_mean(Y,20);
elseif strcmpi(estimator,'hd')
    est = limo_harrell_davis(Y,.5);
elseif strcmpi(estimator,'median')
    est = median(Y,2);
end

% sample with replacement from Dirichlet
n = size(Y,2); % Dirichlet sampling = number of observations, e.g. participants
% n = 1000; % in some situations it might be beneficial to sample more observations
bb = zeros(size(Y,1),Nb);
parfor boot = 1:Nb % bootstrap loop
    theta = exprnd(1,[n,1]);
    weigths = theta ./ repmat(sum(theta,1),n,1);
    resample = (datasample(Y',n,'Replace',true,'Weights',weigths))';
    
    % compute the estimator
    if strcmpi(estimator,'mean')
        bb(:,boot) = mean(resample,2);
    elseif strcmpi(estimator,'trimmed mean')
        bb(:,boot) = limo_trimmed_mean(resample,20); % default 20% trimmed mean
    elseif strcmpi(estimator,'HD')
        bb(:,boot) = limo_harrell_davis(resample,.5); % default to estimation of the 5th decile
    elseif strcmpi(estimator,'median')
        bb(:,boot) = median(resample,2);
    end
end

sorted_data = sort(bb,2); % sort bootstrap estimates
upper_centile = floor(prob_coverage*size(sorted_data,2)); % upper bound
nCIs = size(sorted_data,2) - upper_centile;
HDI = zeros(2,size(Y,1));
ci = 1:nCIs;
ciWidth = sorted_data(:,ci+upper_centile) - sorted_data(:,ci); % all centile distances
[~,J] = min(ciWidth,[],2);
r = size(sorted_data,1);
I = (1:r)';
index = I+r.*(J-1); % linear index
HDI(1,:) = sorted_data(index);
index = I+r.*(J+upper_centile-1); % linear index
HDI(2,:) = sorted_data(index);

