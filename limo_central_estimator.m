function [est,HDI,bb] = limo_central_estimator(Y,estimator,prob_coverage)

% routine to compute high density intervals based on bayesian bootstrap
% estimates ie bayes bootstrap the data, compute the estimator, get the CI
%
% FORMAT [est,ci] = limo_central_estimator(Y,'estimator',prob);
%
% INPUT Y is a 2D matrix eg time frames * subjects
%       estimator is 'Mean', 'Trimmed mean', 'HD' or 'Median'
%       prob is the probability coverage eg 95%
%
% OUTPUT est is the estimator
%        ci is the high density interval
%        bb is a vector of (Bayes) bootstraped estimators 
%
% Guillaume ROusselet & Cyril Pernet January 2016
% ------------------------------------------
% Copyright (C) LIMO Team 2016

Nb = 1000; % number of bootstrap samples 

if nargin == 2
    prob_coverage = 0.95;
elseif nargin == 1
    prob_coverage = 0.95;
    estimator = 'Trimmed mean';
end

% compute the estimator
if strcmpi(estimator,'mean')
    est = mean(Y,2);
elseif strcmpi(estimator,'trimmed mean')
    est = limo_trimmed_mean(Y,20);
elseif strcmpi(estimator,'HD')
    est = limo_harrell_davis(Y,.5);
elseif strcmpi(estimator,'median')
    est = median(Y,2);
end

% sampling with replcaement from Dirichlet
n=size(Y,2); 
bb = zeros(size(Y,1),Nb);
parfor boot=1:Nb
    theta = exprnd(1,[n,1]);
    weigths = theta ./ repmat(sum(theta,1),n,1);
    resample= (datasample(Y',n,'Replace',true,'Weights',weigths))';
    
    % compute the estimator
    if strcmpi(estimator,'mean')
        bb(:,boot) = mean(resample,2);
    elseif strcmpi(estimator,'trimmed mean')
        bb(:,boot) = limo_trimmed_mean(resample,20);
    elseif strcmpi(estimator,'HD')
        bb(:,boot) = limo_harrell_davis(resample,.5);
    elseif strcmpi(estimator,'median')
        bb(:,boot) = median(resample,2);
    end
end

sorted_data = sort(bb,2); % sort bootstrap estimates
upper_centile = floor(prob_coverage*size(sorted_data,2)); % upper bound
nCIs = size(sorted_data,2) - upper_centile;
HDI = zeros(2,size(Y,1));
for frame = 1:size(Y,1)
    tmp = sorted_data(frame,:);
    ci = 1:nCIs; ciWidth = tmp(ci+upper_centile) - tmp(ci); % all centile distances
    [~,index]=find(ciWidth == min(ciWidth)); % densest centile
    if length(index) > 1; index = index(1); end % many similar values
    HDI(1,frame) = tmp(index);
    HDI(2,frame) = tmp(index+upper_centile);
end

