function [b,W] = limo_WLS(X,Y)

% LIMO_WLS Limo Weighted Least Squares (WLS)
% WLS is used to find the maximum likelihood estimates of a generalized
% linear model, and in robust regression to find an M-estimator, as a way 
% of mitigating the influence of outliers in an otherwise 
% normally-distributed data set. 
%
% Weights are obtained using a Principal Components Projection
%
% FORMAT: [b w] = limo_WLS(X,Y)
%
% INPUTS:
%   X             = the design matrix 
%   Y             = 2D matrix of EEG data (dim trials x frames)
%
% OUTPUTS:
%   b             = betas (dim parameters * time frames)
%   w             = weights (dim trials)
%
% References:
%   P. Filzmoser, R. Maronna, M. Werner (2007). Outlier identification in 
%       high dimensions, Computational Statistics and Data Analysis. 
%    
% see also LIMO_PCOUT LIMO_IRLS
%
% Cyril Pernet v1 Octibre 2012
% -----------------------------
% Copyright (C) LIMO Team 2012

%% input check
if  nargin < 1      
    error(message('Too Few Inputs'));      
end 

[rows,cols] = size(X);
if (rows <= cols)
   error('WLS cannot be computed, there is not enough trials for this design');     
end

Y = Y(:,mad(Y,1) > 1e-6); 
if isempty(Y)
    error('WLS cannot be computed, for at least 1 condition, all trials have the same values')
end


%% get the weights from PC on adjusted residuals

% Hat matrix; Leverages for each observation
H = diag(X*pinv(X'*X)*X');

% Adjustment factor
adjfactor = 1 ./ sqrt(1-H);

% Get residuals from OLS
b = pinv(X)*Y;
res = Y - X*b;
resadj = res .* repmat(adjfactor, 1, size(Y,2));

% Robust Estimator
% 0.6745 is the 0.75 quantile of the standard normal distribution
% (makes the estimate unbiased)
re = median(abs(resadj)) ./ 0.6745;
re(find(re < 1e-5)) = 1e-5;
tune = 4.685; % tuning function for the bisquare see limo_IRLS
r = resadj ./ repmat(tune.*re, size(Y,1),1);

% projection and weights 
[W,out] = limo_pcout(r);
WY = Y .* repmat(W,1,size(Y,2));
WX = X .* repmat(W,1,size(X,2));
b = pinv(WX)*WY;

