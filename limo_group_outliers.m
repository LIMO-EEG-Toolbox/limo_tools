function [all_weights, channel_weights, all_errors, channel_errors,  ...
    recon_betas, recon_betas_all_weighted, recon_betas_channel_weighted] = ...
         limo_group_outliers(Beta_files, expected_chanlocs,framestart,...
         frameend,adjacency_matrix,saveGAE)

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
%       expected_chanlocs - A file name or struct with expected channel
%                          locations -- allows having same channels for all
%                          subhects
%       framestart        - 1st frame in time or freq, or [freq time]
%       frameend          - last frame in time or freq, or [freq time]
%       adjacency_matrix  - binary neighbouring matrix fo the graph
%       saveGAE           - 'no' (default) or 'yes' or 'xai' -- 'yes' only 
%                           stores basic data, 'xai' performs explainable AI 
%                           operations, which is more time-consuming.
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
        expected_chanlocs,max(framestart),min(frameend),S.Betas);
    for nb = size(out,3):-1:1
        data(nb,:,:,subject) = squeeze(out(:,:,nb))';
    end
    clear tmp S out
end

%% 1) Autoencoder -- arguments in: matrix of beta values, neighbourgh_matrix
%             -- argument out: learned matrix of beta values
% needed beta_values shape: [nBeta, nTime, nChan, nSubject], neighbourgh_matrix shape:[nChan,nChan] 
% use pyenv('Version',"your python path") to change matlab calls
PYTHONENVIRONMENT = pyenv;
warning('MATLAB is now calling python %s\n',PYTHONENVIRONMENT.Library)
warning('be patient - running the GAE ... ')

% Processing saveGAE parameter
if isempty(saveGAE)
    saveGAE = 'no';
else
    saveGAE = lower(saveGAE);
    % Add python file path to python sys.path
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);
    if count(py.sys.path, thisDir) == 0
        insert(py.sys.path, int32(0), thisDir);
    end
end

learned_betas = pyrunfile("NiPyAEoutliers.py", "learned_betas", ...
    datain = data, binatry_matrix = adjacency_matrix, saveGAE = saveGAE);
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

%% 3) Initialize error containers
all_errors                   = zeros(1, nSubj);          % per-subject global error
channel_errors               = zeros(nChan, nSubj);      % per-channel, per-subject error
recon_betas_all_weighted     = zeros(size(recon_betas)); % Weighted versions:
recon_betas_channel_weighted = zeros(size(recon_betas));

%% 4) Loop over subjects to compute errors
for iSubj = 1:nSubj
    fprintf('computing errors and weights subject %d\n',iSubj)
    fpath    = fileparts(Beta_files{iSubj});
    name     = limo_get_subname(Beta_files{iSubj});
    Yhat     = load(fullfile(fpath,[name '_desc-Yhat.mat']));
    Yhat     = Yhat.Yhat;
    tmp      = load(fullfile(fpath,'LIMO.mat')); 
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

for iSubj = 1:nSubj
    fpath                   = fileparts(Beta_files{iSubj});
    tmp                     = load(fullfile(fpath,'LIMO.mat')); 
    LIMO                    = tmp.LIMO;    
    LIMO.weighting.global   = all_weights(iSubj);
    LIMO.weighting.channels = channel_weights(:,iSubj);
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO',"-v7.3");
end


% reformat for output
recon_betas = permute(recon_betas,[1 3 2 4]);

% save all_weights and channel_weights
if strcmpi(saveGAE,'yes') || strcmpi(saveGAE,'xai')
    newdir = 'Group_outlier_parametrization';
    if ~exist(newdir,'dir')
        mkdir(newdir);
    end
    cd(newdir)
    save('GAE-based_all_weights',"all_weights");
    save('GAE-based_channel_weights',"channel_weights");
    cd ..
end

% plot channel and region masking and channel ablation results on scalp
if strcmpi(saveGAE, 'xai')
   thisDir = pwd();
   cd('Group_outlier_parametrization'); cd('masking')
   imp_data = load("imp_chan_all.mat").imp_chan_all(1:end-1, :, :); % remove the last beta
   plot_importance_on_scalp(imp_data, expected_chanlocs, 'Channel Masking')
   imp_data = load("imp_region_all.mat").imp_region_all(1:end-1, :, :);
   plot_importance_on_scalp(imp_data, expected_chanlocs, 'Region Masking')

   cd ..; cd('ablation')
   imp_data = load("imp_chan_drop_all.mat").imp_chan_drop_all(1:end-1, :, :);
   plot_importance_on_scalp(imp_data, expected_chanlocs, 'Channel Ablation')

   cd(thisDir)
end   

function plot_importance_on_scalp(imp_data, expected_chanlocs, type)
    plotSaveFolder = [type ' Scalp Plots'];
    if ~exist(plotSaveFolder, 'dir')
        mkdir(plotSaveFolder);
    end

    [nBeta, ~, ~] = size(imp_data);
    chan_mean_norm = zeros(size(imp_data));
    eps_val = 1e-12;

    % Normalize per beta
    for iBeta = 1:nBeta
        X = squeeze(imp_data(iBeta,:,:));
        denom = max(X(:));
        denom = max(denom, eps_val);
        chan_mean_norm(iBeta,:,:) = X ./ denom;
    end

    % mean across subjects
    beta_chan_mean = squeeze(mean(chan_mean_norm, 2)); % [nBeta, nChan]
    % aggregated across betas
    final_chan_imp = mean(beta_chan_mean, 1)'; % [nChan,1]

    % helper function - plot one scalp
    function plot_one(val, title_str, filename)
        val_centered = val - mean(val);
        c = max(abs(val_centered));
        figure('Color','w','NumberTitle','off','Name','limo_best_electrodes.m');
        val_plot = val;
        cmax = max(val_plot);
        opt = {'electrodes','on', ...
           'maplimits',[0 cmax], ...
           'verbose','off', ...
           'colormap', limo_color_images(val_plot)};
        
        topoplot(val_centered, expected_chanlocs, opt{:});
        title(title_str);
        exportgraphics(gcf, fullfile(plotSaveFolder, filename));
        close(gcf);
        fprintf('Saved scalp plot: %s\n', fullfile(plotSaveFolder, filename));
    end

    % plot aggregated scalp
    filenameAgg  = 'Aggregated_Imp_scalp.png';
    titleAgg     = 'Aggregated Importance';
    plot_one(final_chan_imp, titleAgg, filenameAgg);

    % plot each beta scalp
    for iBeta = 1:nBeta
        vals_beta = beta_chan_mean(iBeta, :)';  % [nChan,1]
        filenameB = sprintf('Beta_%d_Imp_scalp.png', iBeta);
        titleB    = sprintf('Beta %d Importance',iBeta);
        plot_one(vals_beta, titleB, filenameB);
    end

end

end
