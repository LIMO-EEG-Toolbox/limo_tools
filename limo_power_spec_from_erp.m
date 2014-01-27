function [limo_psd,limo_psd_freqlist] = limo_power_spec_from_erp(EEG,freq_range,winsize,power_dataset_savename)

% function [limo_psd limo_psd_freqlist] = limo_power_spec_from_erp(EEG,freq_range,winsize,power_dataset_savename)
% 
% Takes epoched EEG data and generates the spectral power (power spectrum density)
% at the specified frequencies for each electrode and at each epoch. This
% allows this power data to be examined in LIMO.
% 
% Expects EEG to be the EEGLAB structure, including sample rate, and EEG data to be epoched.
% 
% Default values:
% freq_range = [1:80];
% winsize = EEG.pnts; - length of the window (in data points)
% power_dataset_savename = [EEG.setname '_power.set']
% 
% OUTPUTS:
% 
% limo_spec_power is a 3D array of electrodes x psd freq bins x epochs
% psdfreqs is a list of the center of each of the frequency 
%
% See also EEGLAB's spectopo and std_spec

% -----------------------------
%  Copyright (C) LIMO Team 2013

% Andrew Stewart, University of Edinburgh 2013
% v1 23 sep 2013

% Set defaults if not set

if EEG.nbchan == 0 || isstruct(EEG)~=1     % if the EEG structure is empty or absent, fill it

warndlg(' EEG is empty. On the next window, please select the EEGLAB .set file for which you wish to generate PSD data');
pause(2);

EEG=pop_loadset

end

if nargin<4;power_dataset_savename = [EEG.setname,'_power.set'];end
if nargin<3;freq_range=[1:80];end
if nargin<2;winsize=EEG.pnts;end




% Go through each trial, generating PSD for each trial
    for t=1:EEG.trials
        
        [powerdata_here,limo_psd_freqlist] = spectopo(EEG.data(:,:,t),0,EEG.srate,'chanloc',EEG.chanlocs,'winsize',winsize,'plot','off');
        
        % note - spectopo runs on all freqs up to sr/2, not just our range, as this is easy to
        % calculate
        
        
        
        EEG.etc.psd_max_freq=find(limo_psd_freqlist<max(freq_range), 1, 'last' );
        
        if t==1   % On first loop, set size of psd data
            nfreqbins=numel(limo_psd_freqlist);        
            EEG.etc.limo_psd = zeros(EEG.nbchan,EEG.etc.psd_max_freq,EEG.trials);
            
            freq_step=limo_psd_freqlist(2);
        end
        
        
        EEG.etc.limo_psd(:,:,t) = powerdata_here(:,1:EEG.etc.psd_max_freq);  % Save only up desired freq
        
        out=['Currently on ',EEG.setname,'_and trial_',num2str(t),'_with freq_step_',num2str(freq_step),'_top frequency_',num2str(limo_psd_freqlist(EEG.etc.psd_max_freq))]
        
    end
    
    limo_psd = EEG.etc.limo_psd;
    EEG.etc.limo_psd_freqlist = limo_psd_freqlist(1:EEG.etc.psd_max_freq);
    
    EEG = pop_saveset(EEG, 'filename',power_dataset_savename);
    
end
    
    
    
    
    
        