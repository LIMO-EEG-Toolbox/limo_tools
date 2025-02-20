function [all_weights, channel_weights, all_errors, channel_errors,  ...
    recon_betas, recon_betas_all_weighted, recon_betas_channel_weighted] = ...
         limo_group_outliers(Beta_files, expected_chanlocs,framestart,frameend,adjacency_matrix)

% -------------------------------------------------------------------------
% LIMO_GROUP_OUTLIERS
%
%   This function:
%   1) Loads or matches channel information using 'expected_chanlocs'.
%   2) Calls a Python script (NiPyAEoutliers.py) via pyrunfile to get 
%      reconstructed beta values (learned_betas).
%   3) Compares these reconstructed betas against the 
%      original betas, computes error metrics, and normalizes them.
%   4) Derives weighting factors at the subject (global) and channel 
%      (local) levels, then applies these weights to the reconstructed data.
%   5) Identifies outliers based on a user-defined threshold.
%
%   INPUTS:
%       Beta_files        - A string or structure pointing to the 
%                          location(s) of .mat or .txt with Beta data.
%                          (Currently placeholder in this code.)
%       expected_chanlocs - A file name or struct with expected channel
%                          locations -- allows having same channels for all
%                          subhects
%       framestart        - 1st frame in time or freq, or [freq time]
%       frameend          - last frame in time or freq, or [freq time]
%       adjacency_matrix  - binary neighbouring matrix fo the graph
%
%   OUTPUTS:
%       all_weights                  : 1 x nSubj array containing the global
%                                      (subject-level) weighting factors.
%                                      exp(-zscore(all_errors))
%       channel_weights              : nChan x nSubj array containing the 
%                                      channel-level weighting factors.
%                                      exp(-zscore( channel_errors ))
%       all_errors                   : 1 x nSubj array with the global 
%                                      (subject-level) error metric.
%                                      Average of abs(Yhat - Yhat reconstructed)
%       channel_errors               : nChan x nSubj array with the channel-
%                                      level error metric.
%                                      Average of abs(Yhat - Yhat reconstructed)
%                                      over trials and time
%       recon_beta                   : The reconstructed betas,
%                                      size(betas) as matching inputs.
%
%
% Renxiang Qiu  & Cyril Pernet 
% see limo_random_select.m, limo_random_robust.m
% ------------------------------
%  Copyright (C) LIMO Team 2025

%% organize inputs
% to do check if framestart,frameend can be vectors from limo_randon_select

for subject = 1:max(size(Beta_files))
    tmp = load(Beta_files{subject}); S.Betas = tmp.Betas;
    tmp = load(fullfile(fileparts(Beta_files{subject}),'LIMO.mat')); S.LIMO = tmp.LIMO;
    % load LIMO.mat of the subject use LIMO.data.chanloc
    out = limo_match_elec(S.LIMO.data.chanlocs, ...
        expected_chanlocs,framestart,frameend,S.Betas);
    for nb = size(out,3):-1:1
        data(nb,:,:,subject) = squeeze(out(:,:,nb))';
    end
    clear tmp S out
end

%% 1) Autoencoder -- arguments in: matrix of beta values, neighbourgh_matrix
%             -- argument out: learned matrix of beta values
% needed beta_values shape: [nBeta, nTime, nChan, nSubject], neighbourgh_matrix shape:[nChan,nChan] 
warning('be patient - running the GAE ... ')
% pyenv('Version',"your python path")

learned_betas = pyrunfile("NiPyAEoutliers.py", "learned_betas", ...
    datain = data, binatry_matrix = adjacency_matrix);
learned_betas = struct(learned_betas);

%% 2) Extract or reorder 'learned_betas' from the Python struct

%    with shape [nBeta, nSubj, nChan, nTime].
if isfield(learned_betas, 'reconstructed')
    recon_all = learned_betas.reconstructed; % shape: (nBeta, nSubj, nChan, nTime)
else
    error('learned_betas struct has no field "reconstructed".');
end

% Reorder to match (nBeta, nChan, nTime, nSubj)
% --> we permute dimensions 2 <-> 4:
% original: (1=nBeta, 2=nSubj, 3=nChan, 4=nTime)
% desired:  (1=nBeta, 2=nChan, 3=nTime, 4=nSubj)
recon_betas = double(recon_all);
recon_betas = permute(recon_betas, [1 3 4 2]);  % now shape = (nBeta, nChan, nTime, nSubj)
[~, nChan, ~, nSubj] = size(recon_betas);
% save allinfo


%% 3) Initialize error containers
all_errors                   = zeros(1, nSubj);          % per-subject global error
channel_errors               = zeros(nChan, nSubj);      % per-channel, per-subject error
recon_betas_all_weighted     = zeros(size(recon_betas)); % Weighted versions:
recon_betas_channel_weighted = zeros(size(recon_betas));

%% 4) Loop over subjects to compute errors
for iSubj = 1:nSubj
    fprintf('computing errors and weights subject %d\n',iSubj)
    Yhat     = load(fullfile(fileparts(Beta_files{iSubj}),'Yhat.mat')); 
    Yhat     = Yhat.Yhat;
    tmp      = load(fullfile(fileparts(Beta_files{iSubj}),'LIMO.mat')); 
    X_matrix = tmp.LIMO.design.X;

    % -------------------------------------
    % a) Compute "sub_error" = (original - reconstructed) 
    %    for all betas simultaneously.
    %    shape: (nBeta, nChan, nTime)
    % -------------------------------------
    recon_subj = recon_betas(:,:,:,iSubj);
     for channel = size(recon_subj,2):-1:1
        Yhat_recon(channel,:,:) =  (X_matrix*squeeze(recon_subj(:,channel,:)))';
    end
    deltaYhat  = Yhat - Yhat_recon; 
    
    % absolute value
    deltaYhat = abs(deltaYhat);              % shape (betas, nChan, nTime)

    % b) "Compute mean errors"
    all_errors(iSubj) = mean(deltaYhat(:));  % global average
        
    % c) "sub_error_per_channel" = average over frames (time) and someDim
    tmp_mean                 = mean(deltaYhat, 3);   % average over trials
    channel_errors(:, iSubj) = mean(tmp_mean, 2);    % mean over time dimension => (nChan,1)

    clear Yhat X_matrix Yhat_recon % subject specific stuff
end

%% 5) Normalize errors globally
%    i.e., Z-score across subjects
all_errors_normalized = (all_errors - mean(all_errors)) ./ std(all_errors);
% For channel_errors, we have shape (nChan, nSubj). 
% We can normalize across *all values* (the simplest approach),
channel_mean = mean(channel_errors(:));
channel_std  = std(channel_errors(:));
channel_errors_normalized = (channel_errors - channel_mean) ./ channel_std;

%% 6) Compute weights -- when error goes up we want small weight so -error
%    all_weights => exp(-z) for each subject
all_weights = exp(-all_errors_normalized);            % size = (1, nSubj)
%    channel_weights => exp(-z) for each channel & subject
channel_weights = exp(-channel_errors_normalized);    % size = (nChan, nSubj)

% reformat for output
recon_betas = permute(recon_betas,[1 3 2 4]);

end
