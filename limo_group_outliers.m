function [all_weights, channel_weights, all_errors, channel_errors, outliers, ...
    recon_betas, recon_betas_all_weighted, recon_betas_channel_weighted] = ...
         limo_group_outliers(Beta_files, expected_chanlocs,framestart,frameend,adjacency_matrix,X_matrix)

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
%       channel_weights              : nChan x nSubj array containing the 
%                                      channel-level weighting factors.
%       all_errors                   : 1 x nSubj array with the global 
%                                      (subject-level) error metric.
%       channel_errors               : nChan x nSubj array with the channel-
%                                      level error metric.
%       outliers                     : Indices of subjects whose global
%                                      weight is below a specified threshold.
%       recon_beta                   : The unweighted reconstructed betas,
%                                      size (nBeta, nChan, nTime, nSubj).
%       recon_beta_all_weighted      : Subject-level (global) weighted
%                                      reconstruction, same size as recon_beta.
%       recon_beta_channel_weighted  : Channel-level weighted reconstruction,
%                                      same size as recon_beta.

%% organize inputs

for subject = 1:max(size(Beta_files))
    tmp = load(Beta_files{subject}); S.Betas = tmp.Betas;
    tmp = load(fullfile(fileparts(Beta_files{subject}),'LIMO.mat')); S.LIMO = tmp.LIMO;
    % load LIMO.mat of the subject use LIMO.data.chanloc
    out = limo_match_elec(S.LIMO.data.chanlocs, ...
        expected_chanlocs,framestart,frameend,S.Betas);
    for nb = 1:size(out,3) 
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
% save learned_betas

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
recon_all_mat = double(recon_all);
recon_all_mat = permute(recon_all_mat, [1 3 4 2]);  % now shape = (nBeta, nChan, nTime, nSubj)
[nBeta, nChan, nTime, nSubj] = size(recon_all_mat);


%% 3) Initialize error containers
all_errors = zeros(1, nSubj);            % per-subject global error
channel_errors = zeros(nChan, nSubj);    % per-channel, per-subject error



% We will also store the final reconstructions for reference
recon_betas  = recon_all_mat;   % same shape as betas
% Weighted versions:
recon_betas_all_weighted      = zeros(size(recon_betas));
recon_betas_channel_weighted  = zeros(size(recon_betas));

%% 4) Loop over subjects to compute errors
for iSubj = 1:nSubj
    % -------------------------------------
    % a) Compute "sub_error" = (original - reconstructed) 
    %    for all betas simultaneously.
    %    shape: (nBeta, nChan, nTime)
    % -------------------------------------
    betas = permute(data, [1 3 2 4]);
    orig_subj = betas(:,:,:,iSubj);    % (nBeta, nChan, nTime)
    recon_subj = recon_betas(:,:,:,iSubj); 
    sub_error = orig_subj - recon_subj; 
    
    
    % If X_matrix is size (nTrials, nBeta)
    X_transposed = X_matrix(1:nBeta, :)';  
    % X_transposed has shape (someDim, nBeta).
    % sub_error we want to multiply along the "beta" dimension. 
    
    
    % sub_error is [nBeta, nChan, nTime]
    % We want deltaYhat = [someDim, nChan, nTime]
    % i.e., for each (chan, time), we multiply X_transposed * sub_error(:,chan,time)
    [someDim, ~] = size(X_transposed);
    deltaYhat = zeros(someDim, nChan, nTime);
    
    for ch = 1:nChan
        for t = 1:nTime
            % vector of length nBeta
            err_vec = sub_error(:, ch, t); 
            % multiply: [someDim x nBeta] * [nBeta x 1] => [someDim x 1]
            deltaYhat(:, ch, t) = X_transposed * err_vec;
        end
    end
    
    % absolute value
    deltaYhat = abs(deltaYhat);  % shape (someDim, nChan, nTime)
    
    % b) "Compute mean errors"
    tmp_mean = mean(deltaYhat, 1);             % average over dimension=1 (someDim)
    tmp_mean = squeeze(tmp_mean);              % shape -> (nChan, nTime)
    sub_error_all = mean(tmp_mean(:));         % global average
    all_errors(iSubj) = sub_error_all;
    
    % c) "sub_error_per_channel" = average over frames (time) and someDim
    tmp_mean_over_time = mean(tmp_mean, 2);    % mean over time dimension => (nChan,1)
    sub_errors_per_channel = tmp_mean_over_time; 
    channel_errors(:, iSubj) = sub_errors_per_channel;
end

%% 5) Normalize errors globally
%    i.e., Z-score across subjects
all_errors_normalized = (all_errors - mean(all_errors)) ./ std(all_errors);
% For channel_errors, we have shape (nChan, nSubj). 
% We can normalize across *all values* (the simplest approach),
allvals = channel_errors(:);
channel_mean = mean(allvals);
channel_std  = std(allvals);
channel_errors_normalized = (channel_errors - channel_mean) ./ channel_std;

%% 6) Compute weights
%    all_weights => exp(-z) for each subject
all_weights = exp(-all_errors_normalized);            % size = (1, nSubj)
%    channel_weights => exp(-z) for each channel & subject
channel_weights = exp(-channel_errors_normalized);    % size = (nChan, nSubj)

%% 7) Second loop: Apply weights to reconstructions
%  recon_beta is already the unweighted reconstruction = recon_all
recon_betas = recon_all_mat;  % shape (nBeta, nChan, nTime, nSubj)

% Weighted arrays (same shape)
for iSubj = 1:nSubj
    w_subj = all_weights(iSubj);            % scalar
    w_chan = channel_weights(:, iSubj);     % (nChan x 1)
    
    for iBeta = 1:nBeta
        % Extract the reconstruction for subject iSubj, beta iBeta
        recon_iBeta = recon_betas(iBeta, :,:, iSubj);  % shape (1, nChan, nTime)
        
        % -- (a) Subject-level (global) weighting:
        recon_betas_all_weighted(iBeta, :,:, iSubj) = w_subj * recon_iBeta;
        
        % -- (b) Channel-wise weighting:
        % Multiply each channel by w_chan(ch).
        for ch = 1:nChan
            recon_betas_channel_weighted(iBeta, ch, :, iSubj) = w_chan(ch) * recon_iBeta(:,ch, :);
        end
    end
end

%% 8) Optionally define outliers
thresh = 2; 
outliers = find(all_errors_normalized > thresh);

end
