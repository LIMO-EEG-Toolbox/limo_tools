function channeighbstructmat = limo_neighbourdist(EEG,neighbourdist)

% This function creates a neighbourhood 
% distance matrix used to control for multiple comparisons.
% For explanations see limo_ft_neighbourselection
% http://fieldtrip.fcdonders.nl/tutorial/cluster_permutation_timelock
%
% FORMAT: channeighbstructmat = limo_neighbourdist(EEG,neighbourdist)
%
% INPUT EEG the eeglab structure
%       neighbourdist (optional) the distance between electrodes
%
% OUTOUT: channeighbstructmat a structure with the channel info and the
%         neigbouring matrix
%
% Guillaume Rousselet v1 11 July 2010
% adapted from code by Arnaud Delorme
% ------------------------------
%  Copyright (C) LIMO Team 2019

% reduce data array (thx Benedikt Ehinger)
EEG.data   = EEG.data(:,1,1);
EEG.trials = 1;
EEG.pnts   = 1;

% use fieldtrip for the layout and make the matrix
tmpcfg        = limo_eeglab2fieldtrip(EEG, 'preprocessing', 'none');
lay           = limo_ft_prepare_layout(tmpcfg, tmpcfg); % fieldtrip function
tmpcfg.layout = lay;
if nargin == 1
    neighbourdist = eval(cell2mat(inputdlg('enter neighbourhood distance','neighbourhood distance'))); % 0.37 for biosemi 128;
end
tmpcfg.neighbourdist             = neighbourdist;
[neighbours,channeighbstructmat] = limo_ft_neighbourselection(tmpcfg, []); % fieldtrip functio
