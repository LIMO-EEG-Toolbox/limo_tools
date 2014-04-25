function robust_cov = limo_robust_cov(data,percent)

% computes the winsorized covariance matrix of Data
% INPUT data is a 2D matrix of trials/subjects * measures
% OUTPUT robust_cov is the covariance of winsorized data
%
% Cyril Pernet & Guillaume Rousselet April 2014
% -----------------------------------------------------------------
%  Copyright (C) LIMO Team 2014

if nargin == 1
    percent = 20/100;
end

[n,measures]=size(data);

% number of items to winsorize and trim
g=floor(percent*n);

% winsorise data
[xsort,indices]=sort(data,1);
windata = data;
for j=1:measures
    windata(indices(1:g+1,j),j) = xsort(g+1,j);
    windata(indices(n-g:end,j),j) = xsort(n-g,j);
end

robust_cov = cov(windata);


