function [est,HDI,bb] = limo_central_estimator(Y,estimator,prob_coverage)

% Compute a data estimator and its highest density intervals (HDI) based 
% on bayesian bootstrap estimates.
%
% FORMAT [est,HDI,bb] = limo_central_estimator(Y,'estimator',prob_coverage);
%
% INPUTS Y is a 2D matrix, e.g. frames x participants or trials
%        estimator is 'Mean', 'Trimmed mean', (default) 'HD' (Mid-decile Harell-Davis) or 'Median'
%        prob_coverage is the probability coverage- default 0.95%
%
% OUTPUT est is the estimator
%        ci is the high density interval
%        bb is a vector of (Bayes) bootstraped estimators 
%
% Bayesian bootstrap implementation based on orignal R code from Rasmus Baath:
% http://www.sumsar.net/blog/2015/07/easy-bayesian-bootstrap-in-r/
% HDI implementation based on original R code HDIofMCMC from John K. Kruschke:
% https://github.com/boboppie/kruschke-doing_bayesian_data_analysis/blob/master/1e/HDIofMCMC.R
% 
% Modifications by Cedric Cannard (June 17, 2024):
% 1) Weights: align with Dirichlet distribution properties more closely, 
%   using Gamma distribution normalized to sum 1, instead of the exponential 
%   distribution. 
%
% 2) use quantile intervals instead HDIs, which can lead to different 
%   conclusions based on the parameterization of the model. Quantile 
%   intervals are derived directly from the cumulative distribution function
%   (from the sorted posterior bootstrap samples) and are based on 
%   probabilities rather than densities, ensuring coherence across different 
%   parameterizations. Implemeted following new guidelines by: 
%       Etz A, Chávez de la Peña AF, Baroja L, Medriano K, Vandekerckhove J. 
%       The HDI + ROPE decision rule is logically incoherent but we can fix it. 
%       Psychol Methods. 2024 May 23. doi: 10.1037/met0000660. https://pubmed.ncbi.nlm.nih.gov/38780591/
% 
% see also mean, median, limo_trimmed_mean, limo_harrell_davis
%
% Guillaume Rousselet & Cyril Pernet February 2016
% ------------------------------
%  Copyright (C) LIMO Team 2019

Nb = 1000; % number of bootstrap samples 

if nargin == 2
    prob_coverage = 0.95;
elseif nargin == 1
    prob_coverage = 0.95;
    estimator = 'Trimmed mean';
end

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
    % theta    = exprnd(1,[n,1]);   % generate weights using exponential distribution (related to Bayesian approach but not fully conform)
    % weights  = theta ./ repmat(sum(theta,1),n,1);
    theta    = gamrnd(ones(n,1),1); % generate weights using Gamma distribution
    weights  = theta / sum(theta);  % normalize to sum to 1 to align with Dirichlet distribution properties
    resample = (datasample(Y',n,'Replace',true,'Weights',weights))';
    
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

% % vectorized version of the loop above
% ci       = 1:nCIs;
% ciWidth  = sorted_data(:,ci+upper_centile) - sorted_data(:,ci); % all centile distances
% [~,J]    = min(ciWidth,[],2);
% r        = size(sorted_data,1);
% I        = (1:r)';
% index    = I+r.*(J-1); % linear index
% HDI(1,:) = sorted_data(index);
% index    = I+r.*(J+upper_centile-1); % linear index
% HDI(2,:) = sorted_data(index);

% Compute quantile interval following guidelines from Etz et al. (2024)
% instead of HDIs. See function description. 
HDI = zeros(2, size(Y, 1));
for i = 1:size(Y, 1)
    sampleVec = bb(i, :);
    alp = (1 - prob_cov) / 2;
    lim = quantile(sampleVec, [alp, 1 - alp]);
    HDI(:, i) = lim';
end
