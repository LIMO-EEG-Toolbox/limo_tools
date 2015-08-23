function limo_CheckWeight(varargin)

% general function designed to look at the weights computed for eacc trial
% at the each channel for each subject
%
% FORMAT limo_CheckWeight(list_of_LIMO.mat,options)
%
% INPUT list_of_LIMO.mat simple txt file listing were the LIMO.mat to llok
%                        at are located (or actual variable for that list)
%
%       options 'plot deciles' since weights are between 0 and 1,
%                it computes the average response for each decile
%               'test difference' compute an OLS between the good trials
%               weights = 1/0.9 and the bad ones weights 0/0.1
%               'check bias' check that the weights are distributed across
%               trials in a uniform manner, i.e. that not one conditions is
%               more affected than another which would bias the results, but
%               also indicate that something is going on ion the data
% 
% OUTPUT creates a folder called 'Weights_checking' with the different
%        results in it
%
% Cyril Pernet 21-08-2015
% -----------------------------
% Copyright (C) LIMO Team 2015


