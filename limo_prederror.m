function prediction_error = limo_prederror(LIMO,k)

% routine to compute the k-fold prediction error of the model LIMO.design.X
% take N-k data, compute the predicted values for k and the error (Y-Yhat).^2
%
% FORMAT: prediction_error = limo_prederror(LIMO,k)
%
% INPUT: LIMO the variable or variable name for the data to analyse
%        k the folding for the cross-validation
%
% OUTPUT: prediction_error is a matrix channels*frames*3
%         the last dimension contains 1. the error from training set = var(residuals)
%                                    2. the bias estimates = (Y-Yhat).^2
%                                    3. the variance estimates = var(Yhat)
%
% interpreation: underfitting means high bias and low variance
%                overfitting means low bias but high variance
%                sum(predicted_error,ndims(predicted_error)) is the full prediction error
% 
% see Error and Validation - Advanced Methods for Data Analysis (36-402/36-608)
%     http://www.stat.cmu.edu/~ryantibs/advmethods/notes/errval.pdf
%
% Cyril Pernet - October 2020
% ------------------------------
%  Copyright (C) LIMO Team 2020

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
if vargin == 1
    N = size(Yr,ndims(Yr));
    k = floor(N/5);
end
Nfold = round(N/k);

% create place holder variables
if strcmpi(LIMO.Analysis,'Time-Frequency')
    Yr = limo_tf_4d_reshape(Yr,LIMO.data.size3D);    
end
prediction_error = NaN(size(Yr,1),size(Yr,2),3);

% compute
for fold = Nfold:-1:1

end



