function channeighbstructmat = limo_neighbourdist(EEG)

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
% See also LIMO_FT_NEIGHBOURSELECTION
%
% Guillaume Rousselet v1 11 July 2010
% adapted from code by Arnaud Delorme
% -----------------------------
%  Copyright (C) LIMO Team 2010

tmpcfg = limo_eeglab2fieldtrip(EEG, 'preprocessing', 'none');
lay = limo_ft_prepare_layout(tmpcfg, tmpcfg); % fieldtrip function
tmpcfg.layout        = lay;
neighbourdist = eval(cell2mat(inputdlg('enter neighbourhood distance','neighbourhood distance'))); % 0.37 for biosemi 128;
tmpcfg.neighbourdist = neighbourdist;
[neighbours,channeighbstructmat] = limo_ft_neighbourselection(tmpcfg, []); % fieldtrip function
