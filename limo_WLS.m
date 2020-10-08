function [b,W,rf] = limo_WLS(X,Y)

% LIMO_WLS Limo Weighted Least Squares (WLS)
% WLS is used to find the maximum likelihood estimates of a generalized
% linear model, and in robust regression to find an M-estimator, as a way 
% of mitigating the influence of outliers in an otherwise 
% normally-distributed data set. 
%
% Weights are obtained using a Principal Components Projection
%
% FORMAT: [b,W,rf] = limo_WLS(X,Y)
%
% INPUTS:
%   X             = the design matrix 
%   Y             = 2D matrix of EEG data (dim trials x frames)
%
% OUTPUTS:
%   b             = betas (dim parameters * time frames) the parameters of the GLM
%   W             = weights (dim trials) to apply for each trial
%   rf            = the number of dimensions (compoments) removed to get W
%
% References:
%   P. Filzmoser, R. Maronna, M. Werner (2007). Outlier identification in 
%       high dimensions, Computational Statistics and Data Analysis. 
%    
% see also LIMO_PCOUT LIMO_IRLS
%
% Cyril Pernet v2 January 2014
% ------------------------------
%  Copyright (C) LIMO Team 2019

%% input check
disp('check x')
X
disp('check y')
Y
input('ok?')
if  nargin < 2      
    error(message('Too Few Inputs'));      
end 

[rows,cols] = size(X);
if (rows <= cols)
   error('WLS cannot be computed, there is not enough trials for this design');     
end

if isempty(Y(:,mad(Y,1) > 1e-6))
    error('WLS cannot be computed, for at least 1 condition, all trials have the same values')
end


%% get the weights from PC on adjusted residuals

% Hat matrix; Leverages for each observation
H = diag(X*pinv(X'*X)*X');

% Adjustment factor
adjfactor = 1 ./ sqrt(1-H);

% OLS solution
b = pinv(X)*Y;

% tuning function for the bisquare see also limo_IRLS
tune = 4.685; 

% Get residuals from previous fit
res                 = Y - X*b;
resadj              = res .* repmat(adjfactor, 1, size(Y,2));

% re - Robust Estimator
% 0.6745 is the 0.75- quantile of the standard normal distribution
% (makes the estimate unbiased)
re                  = median(abs(resadj)) ./ 0.6745;
re(find(re < 1e-5)) = 1e-5;
r                   = resadj ./ repmat(tune.*re, size(Y,1),1);

%% do the computation
[W,~,rf]            = limo_pcout(r);
WY                  = Y .* repmat(W,1,size(Y,2));
WX                  = X .* repmat(W,1,size(X,2));
b                   = pinv(WX)*WY;

