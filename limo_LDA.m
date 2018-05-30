function [predicted] = limo_LDA(varargin)
% Performs linear classification following Rencher 2002. 
%
% FORMAT [accuracy, confusionmat] = limo_LDA(trainingset,traininglabel, testset, method)
%
% INPUT
% Y = trainingset [trials/subjects x electrodes]
% class = traininglabel (1D)
% testset [trials/subjects x electrodes]
% method = regularized or pseudo inverse 
%
% OUTPUT
% predicted labels 
%
% Iege Bassez v1 May 2018


Y = varargin{1}; 
class = varargin{2}; 
testset = varargin{3};

if nargin < 4
    method = 'pseudo'; % default
else
    method = varargin{4}; 
end
% Calculate covariance matrix
% transform class first to X matrix to calculate E same way as mglm
X     = dummyvar(class);  
X     = [X, ones(size(X,1),1)];
T     = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));  % SS Total
R     = eye(size(Y,1)) - (X*pinv(X));                                      % Projection on E
E     = (Y'*R*Y); % SS Error
n = size(Y,1);
k = length(unique(class)); 
ve = n - k;  
Spooled = E./ve;
meanFeaturesGroup =  grpstats(Y, class, {'mean'});

if strcmp(lower(method), 'pseudo')
    % get linear classification functions with PINV 
    % calculate coefficients and intercept on trainingset
    coeff = meanFeaturesGroup * pinv(Spooled);
    intercept = diag(meanFeaturesGroup * pinv(Spooled) * meanFeaturesGroup')/2;
    % calculate classification functions on testset
    L = (coeff * testset') - repmat(intercept,1,size(testset,1)); % page 306
    Logical = L == repmat(max(L),k,1);
    [predicted, ~] = find(Logical == 1);

elseif strcmp(lower(method), 'regularized')
    % get linear classification functions with Regularized Covariance
    % get regularized covariance matrix 
    [RegularizedCovariance, ~] = cov1para(Y);
    % calculate coefficients and intercept on trainingset
    coeff = meanFeaturesGroup * inv(RegularizedCovariance);
    intercept = diag(meanFeaturesGroup * inv(RegularizedCovariance) * meanFeaturesGroup')/2;
    % calculate classification functions on testset
    L = (coeff * testset') - repmat(intercept,1,size(testset,1)); % page 306
    Logical = L == repmat(max(L),k,1);
    [predicted, ~] = find(Logical == 1);
end
end % end function