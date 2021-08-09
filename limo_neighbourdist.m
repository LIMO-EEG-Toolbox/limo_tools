function channeighbstructmat = limo_neighbourdist(EEG, neighbourdist)

% This function takes as input an EEGLAB dataset (with channel locations)
% and creates a neighbourhood distance matrix used to control for 
% multiple comparisons. A second optional argument indicate the neighbourhood
% distance. If it is not provided, the function pops up a GUI to query for it.
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
% ------------------------------
%  Copyright (C) LIMO Team 2019

tmpcfg = limo_eeglab2fieldtrip(EEG, 'preprocessing', 'none');
lay = limo_ft_prepare_layout(tmpcfg, tmpcfg); % fieldtrip function
tmpcfg.layout        = lay;
if nargin < 2
  neighbourdist = eval(cell2mat(inputdlg('enter neighbourhood distance','neighbourhood distance'))); % 0.37 for biosemi 128;
end
tmpcfg.neighbourdist = neighbourdist;
[neighbours,channeighbstructmat] = limo_ft_neighbourselection(tmpcfg, []); % fieldtrip function
