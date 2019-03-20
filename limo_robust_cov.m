function robust_cov = limo_robust_cov(data,percent)

% Computes the winsorized covariance matrix of Data
%
% FORMAT: limo_robust_cov(data,percent)
% INPUTS: data is a 2D matrix of trials/subjects * measures
%         percent is the amount of trimming (default 20%)
% OUTPUT: robust_cov is the covariance of winsorized data
%
% Cyril Pernet & Guillaume Rousselet April 2014
% -----------------------------------------------------------------
% Copyright (C) LIMO Team 2014

if numel(size(data)) ~=2
    error('data size must be 2D')
end

if nargin == 1
    percent = 20/100;
end

[n,measures]=size(data);

% winsorise data
g=floor(percent*n); % number of items to winsorize
[xsort,indices]=sort(data,1);
windata = data;
for j=1:measures
    windata(indices(1:g+1,j),j) = xsort(g+1,j);
    windata(indices(n-g:end,j),j) = xsort(n-g,j);
end

% get the covariance
h = n-2*floor(percent*n);
robust_cov = zeros(measures,measures);
robust_cov(1,1) = (n-1)*var(windata(:,1))/(h*(h-1));
for j = 2:measures
    jk = j-1;
    robust_cov(j,j) = (n-1)*var(windata(:,j))/(h*(h-1));
    for k = 1:jk
        wincov = cov(windata(:,j),windata(:,k));
        wincov = wincov(1,2);
        robust_cov(j,k) = (n-1)*wincov/(h*(h-1));
        robust_cov(k,j) = robust_cov(j,k);
    end
end

