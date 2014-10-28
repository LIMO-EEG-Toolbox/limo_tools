function [limo_psd,limo_psd_freqlist] = limo_power_spec_from_erp(EEG,options)

% FORMAT [limo_psd limo_psd_freqlist] = limo_power_spec_from_erp(EEG,maxfreq,winsize,power_dataset_savename)
% 
% INPUT    EEG input from EEGLAB EEG.set
%          options is a structure with the following (optional) fields
%          option.maxfreq maximum of the frequency to be kept (default 50Hz)
%          option.winsize is length of the window in data points (default EEG.srate)
%          option.power_dataset_savename [EEG.setname '_power.set']
%          option.ica
%
% OUTPUTS  limo_psd the actual PSD (3D array of electrodes x psd freq bins x epochs)
%          limo_psd_freqlist list of frequencies (bin centered) computed 
% 
% Takes epoched EEG data and generates the spectral power (power spectrum density)
% at the specified frequencies for each channel or comp[onent and at each epoch.
% This allows this power data to be examined in LIMO EEG.
% 
% Expects EEG to be the EEGLAB structure, including sample rate, and EEG data to be epoched.
% 
% See also EEGLAB's spectopo and std_spec
%
% Andrew Stewart, v1 23 sep 2013
% Cyril Pernet, minor changes and edits Jan 2014
% Cyril Pernet Nov 20104 - changed input structure to option adding ica
%                and call spectopo using log + other minor changes
% -------------------------------------------------------------------
%  Copyright (C) LIMO Team 2014


%% Set defaults if not set

if EEG.nbchan == 0 || isstruct(EEG)~=1     % if the EEG structure is empty or absent, fill it    
    warndlg2(' EEG is empty. On the next window, please select the EEGLAB .set file for which you wish to generate PSD data');
    EEG=pop_loadset;
end

if nargin<4;power_dataset_savename = [EEG.setname,'_power.set'];end
if nargin<3;maxfreq=50;end
if nargin<2;winsize=EEG.srate;end

if maxfreq > (EEG.srate/2)
    error('maximum frequency cannot be above the Nyquist limit')
end
    
% Go through each trial, generating PSD for each trial
    for t=1:EEG.trials
        fprintf('computing single trial PSD for trial %g/%g\n',t,EEG.trials)
        %if data
         [powerdata_here,limo_psd_freqlist] = spectopo(EEG.data(:,:,t),0,EEG.srate,'chanloc',EEG.chanlocs,'winsize',winsize,'plot','off');
        %if ica
        
        % note - spectopo runs on all freqs up to sr/2, not just our range
        EEG.etc.psd_max_freq=find(limo_psd_freqlist<=maxfreq, 1, 'last' );
        
        if t==1   % On first loop, set size of psd data
            limo_psd = zeros(EEG.nbchan,EEG.etc.psd_max_freq,EEG.trials);
        end
        limo_psd(:,:,t) = powerdata_here(:,1:EEG.etc.psd_max_freq);  % Save only up desired freq
    end
    
    % if save
    % EEG.etc.limo_psd = limo_psd;
    % EEG.etc.limo_psd_freqlist = limo_psd_freqlist(1:EEG.etc.psd_max_freq);
    % EEG = pop_saveset(EEG, 'filename',power_dataset_savename);
end
    
    
    
    
    
        