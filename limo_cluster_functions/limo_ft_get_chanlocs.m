function [eeg_limo] = ft2limo_get_chanlocs(eeg_limo, defaults)

% This function is used to create the chanlocs structure required by LIMO
% from FieldTrip data.

default_label = defaults.template_elec.label;
default_chanpos = defaults.template_elec.chanpos;


if isfield(eeg_limo,'elec')
    eeg_limo.chanpos = eeg_limo.elec.chanpos;
    eeg_limo.label = eeg_limo.elec.label;
end
    
if ~isfield(eeg_limo,'chanpos')
    disp('ERROR in limo_batch_import_data: Channel positions not defined in EEG_DATA.chanpos')
    disp('Lets use the default electrode set-up');
    eeg_labels = ft_channelselection('eeg', default_label);
    eeg_limo.chanpos = zeros(length(eeg_labels),3);
    for i = 1:length(eeg_labels)
        lab = eeg_labels{i};
        eeg_limo.chanpos(i,:) = default_chanpos(strcmp(default_label,lab),:);
    end
end
if ~isfield(eeg_limo,'label')
    disp('ERROR in limo_batch_import_data: Channel labels not defined in EEG_DATA.label')
    disp('Lets use the default electrode set-up');
    if ~exist(eeg_labels,'var')
        eeg_limo.label = ft_channelselection('eeg', default_label);
    else
        eeg_limo.label = eeg_labels;
    end
end

eeg_limo.chanlocs = struct('labels', eeg_limo.label);
for i = 1:length(eeg_limo.chanpos)
    eeg_limo.chanlocs(i).X = eeg_limo.chanpos(i,1);
    eeg_limo.chanlocs(i).Y = eeg_limo.chanpos(i,2);
    eeg_limo.chanlocs(i).Z = eeg_limo.chanpos(i,3);
end
eeg_limo.chanlocs = convertlocs(eeg_limo.chanlocs,'cart2all');
    
end