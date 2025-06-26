function [neighbours,channeighbstructmat] = limo_get_channeighbstructmat(EEG,neighbourdist)

% This function creates a neighbourhood distance matrix used to control for 
% clustering/multiple comparisons correction, by associating neighbiurgh channels
% together. For more details on the comutation see help limo_ft_neighbourselection
%
% FORMAT [neighbours,channeighbstructmat] = limo_get_channeighbstructmat(EEG,neighbourdist)
%
% INPUTS EEG is the EEGLAB data structure
%            the function works on the current EEGLAB dataset, so the first thing 
%            you need to do is to load an EEGLAB dataset and use this as input. 
%        neighbourdist is the threshold distance that defines  neighbour electrodes. 
%            For instance, 0.37 is a good distance for Biosemi standard 128 electrodes 
%            configuration. You can check the accuracy of neighbourdist for your 
%            electrode montage using the output neighbours
%
% OUTPUTS neighbours structure that lists all the electrodes with their neighbours.
%         channeighbstructmat a matrix of electrode neighbourhood used in cluster analyses.
%
% Reference: Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of 
%            EEG- and MEG-data. J Neurosci Methods, 164(1), 177-190.
%            http://fieldtrip.fcdonders.nl/tutorial/cluster_permutation_timelock
%
% See also LIMO_EEGLAB2FIELDTRIP LIMO_FT_NEIGHBOURSELECTION
% LIMO_FT_PREPARE_LAYOUT LIMO_EXPECTED_CHANLOCS
%
% G.A. Rousselet & C.R. Pernet
% ------------------------------
%  Copyright (C) LIMO Team 2019

if nargin < 2
    error('missing inputs')
end

if isempty(EEG.chanlocs) || isfield(EEG,'chanlocs') == 0
    error('no channel locations found in the EEG/LIMO file');
end

disp('computing neighbourhood using fieldtrip tools .. ') 
tmpcfg               = limo_eeglab2fieldtrip(EEG, 'preprocessing', 'none'); % convert structure
lay                  = limo_ft_prepare_layout(tmpcfg, tmpcfg); % fieldtrip function for electrode layout
tmpcfg.layout        = lay;
tmpcfg.neighbourdist = neighbourdist;

[neighbours,channeighbstructmat] = limo_ft_neighbourselection(tmpcfg, []); % get the matrices

