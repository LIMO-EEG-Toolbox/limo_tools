function [EST,CI,p] = limo_pbci(x,nboot,alpha,mu,est,q)

% function [EST,CI,p] = limo_pbci(x,nboot,alpha,mu,est,q)
%
% Computes a one-sample percentile bootstrap confidence interval
%
%   INPUTS:
%           x (EEG matrix electrodes x frames x trials/subjects)
%           nboot (resamples, default 1000)
%           alpha (default 5%)
%           est (estimator, e.g. 'mean' - default trimmed mean 'tm')
%           q (argument for estimator, default 20% for trimmed mean
%           mu (null hypothesis, necessary to get a p-value in
%           one-sample hypothesis testing procedure, default 0)
% 
%   OUTPUTS:
%           EST = estimator
%           CI  = confidence interval of the estimator
%           p   = associated p value
%
% See also LIMO_ROBUST_CI

% -----------------------------
%  Copyright (C) LIMO Team 2010

% first version: GAR - University of Glasgow - Dec 2007
% added EST output: GAR - Nov 2008
% EEG version: GAR - May 2010
% LIMO version: GAR, June 2010
% GAR, April 2013: replaced randsample with randi; updated help

if nargin<2 || isempty(nboot)
    nboot=1000;
end
if nargin<3 || isempty(alpha)
    alpha=.05;
end
if nargin<5 || isempty(est)
    est='tm';
end

if nargin<6 || isempty(q)
    switch lower(est)
        case {'tm'} % build in trimmed mean function
            q=20; % 20% trimming of both tails, 40% trimming in total
        case {'trimmean'} % Matlab trimmed mean function
            q=40; % Matlab function performs q/2 on both sides
    end
end

switch lower(est)
    case {'trimmean'} % trimmed mean - Matlab function
        eval(['EST=',est,'(x,q,3);'])
    case {'tm'}
        eval(['EST=',est,'(x,q);'])
    otherwise
        eval(['EST=',est,'(x,3);'])
end

n = size(x,3);
lo = round(nboot.*alpha./2);
hi = nboot - lo;
bootx = zeros(nboot,size(x,1),size(x,2));

for kk = 1:nboot % bootstrap with replacement loop
    switch lower(est)
        case {'trimmean'} % trimmed mean - Matlab function
            eval(['bootx(kk,:,:)=',est,'(x(:,:,randi(n,1,n)),q,3);'])
        case {'tm'} % trimmed mean - Matlab function
            eval(['bootx(kk,:,:)=',est,'(x(:,:,randi(n,1,n)),q);'])
        otherwise
            eval(['bootx(kk,:,:)=',est,'(x(:,:,randi(n,1,n)),3);'])
    end
end

diffsort = sort(bootx,1); % sort in ascending order
CI(:,:,1) = diffsort(lo+1,:,:);
CI(:,:,2) = diffsort(hi,:,:);

if nargout > 2
    if nargin < 4 || isempty(mu)
        mu=0;
    end
    p=sum(bootx>mu,1)./nboot;
    p = 2.*min(p,1-p);
end

end

%% compute trimmed mean using tm function

function tma = tm(a,percent)

na=size(a,3); % number of trials
ga=floor((percent/100)*na); % number of items to trim
asort=sort(a,3);
tma=mean(asort(:,:,(ga+1):(na-ga)),3); % trimmed means

end
