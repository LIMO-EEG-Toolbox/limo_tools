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


