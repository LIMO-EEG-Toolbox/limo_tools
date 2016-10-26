function [b,W] = limo_WLS(X,Y,varargin)

% LIMO_WLS Limo Weighted Least Squares (WLS)
% WLS is used to find the maximum likelihood estimates of a generalized
% linear model, and in robust regression to find an M-estimator, as a way 
% of mitigating the influence of outliers in an otherwise 
% normally-distributed data set. 
%
% Weights are obtained using a Principal Components Projection
%
% FORMAT: [b w] = limo_WLS(X,Y,method)
%
% INPUTS:
%   X             = the design matrix 
%   Y             = 2D matrix of EEG data (dim trials x frames)
%   method        is either 'simple' (default) or 'iterative'
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
% Cyril Pernet v2 January 2014
% v3 July 2015 incliude an iterative framework (not validated)
% -----------------------------
% Copyright (C) LIMO Team 2015

%% input check
if  nargin < 2      
    error(message('Too Few Inputs'));      
end 

method = 'simple';
if nargin == 3;
    method = cell2mat(varargin{1});
    if sum([strcmpi(method,'simple') strcmpi(method,'iterative')]) == 0
        error('imput method must be either ''simple'' or ''iterative''')
    end
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
res = Y - X*b;
resadj = res .* repmat(adjfactor, 1, size(Y,2));

% re - Robust Estimator
% 0.6745 is the 0.75- quantile of the standard normal distribution
% (makes the estimate unbiased)
re = median(abs(resadj)) ./ 0.6745;
re(find(re < 1e-5)) = 1e-5;
r= resadj ./ repmat(tune.*re, size(Y,1),1);

%% do the computation

if strcmpi(method,'simple')
    
    [W,out] = limo_pcout(r);
    WY = Y .* repmat(W,1,size(Y,2));
    WX = X .* repmat(W,1,size(X,2));
    b = pinv(WX)*WY;
    
elseif strcmpi(method,'iterative')
    
    % iterate as to min res.
    % set a 100 iteration max
    numiter = 0; iterlim = 100;
    oldRes=1; newRes=10;
    
    while(sum(abs(oldRes-newRes)) > cols*(1E-4)) % on average it is small
        
        numiter = numiter+1;
        oldRes = newRes;
        
        if (numiter>iterlim)
            warning('limo_WLS could not converge');
            break;
        end
        
        % Get residuals from previous fit
        res = Y - X*b;
        resadj = res .* repmat(adjfactor, 1, size(Y,2));
        
        %re - Robust Estimator
        % 0.6745 is the 0.75- quantile of the standard normal distribution
        % (makes the estimate unbiased)
        re = median(abs(resadj)) ./ 0.6745;
        re(find(re < 1e-5)) = 1e-5;
        r= resadj ./ repmat(tune.*re, size(Y,1),1);
        
        % Compute new weights using Principal Component projection
        [W,out] = limo_pcout(r);
        WY = Y .* repmat(W,1,size(Y,2));
        WX = X .* repmat(W,1,size(X,2));
        b = pinv(WX)*WY;
        
        % newRes= sum(sum(res.^2));
        newRes= sum(res(:).^2);
    end
    
end
