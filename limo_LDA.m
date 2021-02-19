function [predicted, confmat, accuracy] = limo_LDA(varargin)
% Performs linear classification following Rencher 2002. 
%
% FORMAT [predicted, confusionmat, accuracy] = limo_LDA(trainingset,traininglabel, testset, method)
%
% INPUT
% Y       = trainingset [trials/subjects x electrodes]
% class   = traininglabel (1D)
% testset = [trials/subjects x electrodes]
% method  = regularized or pseudo inverse covariance
%
% OUTPUT
% predicted labels, confusion matrice, accuracy
%
% compares to https://www.mathworks.com/help/stats/discriminant-analysis.html
% see also https://en.wikipedia.org/wiki/Linear_discriminant_analysis
%
% Iege Bassez v1 May 2018
% C Pernet, some code optimization
% --------------------------------
% Copyright (C) LIMO Team 2018

%% Inputs
Y       = varargin{1}; 
class   = varargin{2}; 
testset = varargin{3};

if nargin < 4
    method = 'pseudo'; % default
else
    method = varargin{4}; 
end

meanFeaturesGroup =  grpstats(Y, class, {'mean'}); % get the means

%% Calculate covariance matrix
if strcmpi(method, 'pseudo')
    % transform class first to X matrix to calculate E same way as mglm
    X                 = dummyvar(class);
    X                 = [X, ones(size(X,1),1)];
    R                 = eye(size(Y,1)) - (X*pinv(X));                                     
    Spooled           = (Y'*R*Y)./(size(Y,1) - length(unique(class)));
    
    % get linear classification functions with PINV
    % calculate coefficients and intercept on trainingset
    coeff     = meanFeaturesGroup * pinv(Spooled);
    intercept = diag(meanFeaturesGroup * pinv(Spooled) * meanFeaturesGroup')/2;
    
elseif strcmpi(method, 'regularized')
    % get regularized covariance matrix
    [RegularizedCovariance, ~] = cov1para(Y);
    % get linear classification functions with Regularized Covariance
    % calculate coefficients and intercept on trainingset
    coeff     = meanFeaturesGroup * inv(RegularizedCovariance);
    intercept = diag(meanFeaturesGroup * inv(RegularizedCovariance) * meanFeaturesGroup')/2;
end

% calculate classification functions on testset
L = (coeff * testset') - repmat(intercept,1,size(testset,1)); % page 306
Logical = L == repmat(max(L),k,1);
[predicted, ~] = find(Logical == 1);

%% additional outputs
if nargout == 2
    confmat = confusionmat(class, predicted);
end

if nargout == 3
    accuracy = sum(diag(confmat/sum(sum(confmat))));
end
