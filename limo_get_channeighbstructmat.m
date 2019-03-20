function [neighbours,channeighbstructmat] = limo_get_channeighbstructmat(EEG,neighbourdist)
% function [neighbours,channeighbstructmat] = limo_get_channeighbstructmat(EEG,neighbourdist)
%
% This function creates a neighbourhood 
% distance matrix used to control for multiple comparisons.
%
% For explanations:
% >>help limo_ft_neighbourselection
%
% Maris, E., & Oostenveld, R. (2007). 
% Nonparametric statistical testing of EEG- and MEG-data. 
% J Neurosci Methods, 164(1), 177-190.
%
% http://fieldtrip.fcdonders.nl/tutorial/cluster_permutation_timelock
% 
% INPUTS:
% -------
%
% The function works on the current EEGLAB dataset, so the first thing you need to do 
%       is to load an EEGLAB dataset. 
% neighbourdist is the threshold distance that defines  neighbour
%   electrodes. For instance, 0.37 is a good distance for Biosemi standard
%   128 electrodes configuration. You can check the accuracy of
%   neighbourdist for your electrode montage using the output NEIGHBOURS
%   structure, which list all the electrodes with their neighbours.
%
% OUTPUTS:
% --------
% neighbours: see description in limo_ft_neighbourselection
%
% channeighbstructmat, a matrix of electrode neighbourhood
% used in cluster analyses. The matrix is calculated using the current EEGLAB dataset.
% -----------------------------
%  Copyright (C) LIMO Team 2014
%
% See also LIMO_EEGLAB2FIELDTRIP LIMO_FT_NEIGHBOURSELECTION
% LIMO_FT_PREPARE_LAYOUT LIMO_EXPECTED_CHANLOCS

disp('computing neighbourhood using fieldtrip tools .. ') 
if isempty(EEG.chanlocs) || isfield(EEG,'chanlocs') == 0
    error('no channel locations found in the EEG/LIMO file');
end
tmpcfg = limo_eeglab2fieldtrip(EEG, 'preprocessing', 'none');
lay = limo_ft_prepare_layout(tmpcfg, tmpcfg); % fieldtrip function 
tmpcfg.layout        = lay;
tmpcfg.neighbourdist = neighbourdist;
[neighbours,channeighbstructmat] = limo_ft_neighbourselection(tmpcfg, []); % fieldtrip function

