function errors = limo_prederror(LIMO,k)

% routine to compute the k-fold prediction error of the model LIMO.design.X
% take N-k data, compute the predicted values for k and the error (Y-Yhat).^2
%
% FORMAT: errors = limo_prederror(LIMO,k)
%
% INPUT: LIMO the variable or variable name for the data to analyse
%        k the number of folds for the cross-validation (k=5 by default)
%
% OUTPUTS: errors is a matrix channels*frames*3
%          the last dimension contains 1. incompressible error    = var(residuals) training set / k
%                                      2. the prediction bias     = (mean(Y test-Yhat predicted).^2 - incompressible error) / k
%                                      3. the prediction variance = var(Yhat predicted) / k
%
% interpreation: underfitting means high bias and low variance
%                overfitting means low bias but high variance
%                sum of errors [1 2] in the last dimension is the mean prediction error
%
% see Error and Validation - Advanced Methods for Data Analysis (36-402/36-608)
%     http://www.stat.cmu.edu/~ryantibs/advmethods/notes/errval.pdf
%
% Cyril Pernet - October 2020
% ------------------------------
%  Copyright (C) LIMO Team 2020

% LIMO = 'F:\WakemanHenson_Faces\eeg\derivatives\sub-002\eeg\FaceRepAll_GLM_Channels_Time_OLS\LIMO.mat';

% check LIMO file
if ischar(LIMO)
    LIMO = load(LIMO);
    LIMO = LIMO.LIMO;
end

% load the data
cd(LIMO.dir)
Yr = load(fullfile(LIMO.dir,'Yr.mat'));
Yr = Yr.Yr;

% check fold - default is 5 folds
if nargin == 1
    N = size(Yr,ndims(Yr));
    k = 5;
end
Nfold = round(N/k);

% create place holder variables
if strcmpi(LIMO.Analysis,'Time-Frequency')
    errors = NaN(size(Yr,1),size(Yr,2),size(Yr,3),3,k);
else
    errors = NaN(size(Yr,1),size(Yr,2),3,k);
end

%% compute
array           = find(~isnan(Yr(:,1,1)));          % all good channels
foldindex       = 1:Nfold:N;                        % all non overlapping indices
foldindex(k+1)  = N;                                
shuffling_index = randperm(N);                      % since Yr is organized per condition
Data            = Yr(:,:,shuffling_index);          % mix it up
Design          = LIMO.design.X(shuffling_index,:); % apply to X

for fold = k:-1:1
    testindex  = foldindex(fold)+1:foldindex(fold+1);
    trainindex = setdiff(1:N,testindex);
    for c = 1:length(array)
        channel    = array(c);
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            for freq = size(Data,2):-1:1
                test       = squeeze(Data(channel,freq,:,testindex));
                train      = squeeze(Data(channel,freq,:,trainindex));
                [errors(channel,freq,:,1,fold),errors(channel,freq,:,2,fold), ...
                    errors(channel,freq,:,3,fold)] = geterrors(Design,LIMO.design.method,train,trainindex,test,testindex);
            end
        else
            test       = squeeze(Data(channel,:,testindex));
            train      = squeeze(Data(channel,:,trainindex));
                [errors(channel,:,1,fold),errors(channel,:,2,fold), ...
                    errors(channel,:,3,fold)] = geterrors(Design,LIMO.design.method,train,trainindex,test,testindex);
        end
    end
end

% average folds
if strcmpi(LIMO.Analysis,'Time-Frequency')
    errors           = mean(errors,5); 
    errors(:,:,:,2)  = errors(:,:,:,2)-errors(:,:,:,1); % mean predicted - incompressible = bias
else
    errors           = mean(errors,4); 
    errors(:,:,2)    = errors(:,:,2)-errors(:,:,1); 
end

% figure;
% subplot(1,3,1); imagesc(prediction_error(:,:,1)); title('irreducible error')
% subplot(1,3,2); imagesc(prediction_error(:,:,2)); title('Bias (prediction error)')
% subplot(1,3,3); imagesc(prediction_error(:,:,3)); title('Variance (prediction var)')
end

function [incompressible_error,prediction_error,prediciton_variance] = ...
    geterrors(Design,method,train,trainindex,test,testindex)
% train set
X = Design(trainindex,:);
if strcmpi(method,'WLS')
    [Betas,W,rf]     = limo_WLS(X,train');
    WX               = X.*repmat(W,1,size(X,2));
    HM               = WX*pinv(WX); % Hat matrix, projection onto X
    dfe              = trace((eye(size(HM))-HM)'*(eye(size(HM))-HM)) - (rf-1);
    Yhat             = WX*Betas; % modelled data
else
    Betas            = pinv(X)*train';
    Yhat             = X*Betas;
    dfe              = size(train,2)-rank(X);
end
residuals            = train'-Yhat;
incompressible_error = diag((residuals'*residuals) / dfe);
% test set
X                    = Design(testindex,:);
Yhat                 = X*Betas;
prediction_error     = mean((test'-Yhat).^2);
prediciton_variance  = var(Yhat);
end