function eeg_limo = limo_get_ft_chanlocs(eeg_limo, defaults)

% This function creates the chanlocs structure required by LIMO
% from FieldTrip data (channel or source).
%
% FORMAT: eeg_limo = limo_get_ft_chanlocs(eeg_limo, defaults)
%
% INPUTS: eeg_limo structure with at least the fields 'filepath' and 'filename'
%         defaults is a structure specifying all the parameters
%                to use in the GLM (ie set in LIMO.mat)
%
% OUTPUT: eeg_limo structure updated with the fields
%         - elec.chanpos and elec.label from fieldtrip
%         - chanlocs as per EEGLAB (for channels and sources)
%
% see also limo_batch, limo_batch_import_data
% -----------------------------------------
%  Copyright (C) LIMO Team 2021

israw      = ft_datatype(eeg_limo, 'raw');
istimelock = ft_datatype(eeg_limo, 'timelock');
isfreq     = ft_datatype(eeg_limo, 'freq');
issens     = ft_datatype(eeg_limo, 'sens');
issource   = ft_datatype(eeg_limo, 'source');
if issens
    eeg_limo.chanlocs = struct('labels', eeg_limo.label);
    for i = 1:length(eeg_limo.chanpos)
        eeg_limo.chanlocs(i).X = eeg_limo.chanpos(i,1);
        eeg_limo.chanlocs(i).Y = eeg_limo.chanpos(i,2);
        eeg_limo.chanlocs(i).Z = eeg_limo.chanpos(i,3);
    end
elseif israw||istimelock||isfreq||istimelock
    
    hasopto = isfield(eeg_limo, 'opto');
    hasgrad = isfield(eeg_limo, 'grad');
    haselec = isfield(eeg_limo, 'elec');
    chanlocs = struct([]);
    if hasopto
        tmp      = limo_get_ft_chanlocs(eeg_limo.opto, defaults);
        chanlocs = cat(2, chanlocs, tmp.chanlocs);
    end
    if hasgrad
        tmp      = limo_get_ft_chanlocs(eeg_limo.grad, defaults);
        chanlocs = cat(2, chanlocs, tmp.chanlocs);
    end
    if haselec
        tmp      = limo_get_ft_chanlocs(eeg_limo.elec, defaults);
        chanlocs = cat(2, chanlocs, tmp.chanlocs);
    end
    eeg_limo.chanlocs = chanlocs;
    
    if isempty(eeg_limo.chanlocs)  
        disp('Channel positions not found in the input data structure')
        disp('Lets use the default electrode set-up, assuming EEG');
        eeg_limo.chanlocs = limo_get_ft_chanlocs(defaults.template_elec, defaults);
    end
        
elseif issource
    
    if ~isfield(eeg_limo,'pos')
        error('ERROR: non existing "pos" field')
    elseif ~isfield(eeg_limo,'label')
        error('ERROR: non existing "label" field')
    end

    eeg_limo.chanlocs = struct('labels', eeg_limo.label); %FIXME this probably does not work
    for i = 1:length(eeg_limo.pos)
        eeg_limo.chanlocs(i).X = eeg_limo.pos(i,1);
        eeg_limo.chanlocs(i).Y = eeg_limo.pos(i,2);
        eeg_limo.chanlocs(i).Z = eeg_limo.pos(i,3);
    end
        
end
