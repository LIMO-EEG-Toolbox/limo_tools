function version = eegplugin_limo(fig,try_strings,catch_strings)

% simple plugin function for eeglab
% call the LInear MOdeling of EEG data toolbox
% Multiple choices are available depending on the type
% of analysis
% 
% Nicolas Chauveau - Cyril Pernet v1 13/03/2009
% Cyril Pernet v2 22-06-2013
% -----------------------------
%  Copyright (C) LIMO Team 2010

version ='LIMO EEG 2.0';
if nargin < 3
    error('eegplugin_limo requires 3 arguments');
end;


% add folder to path
% ------------------
    if ~exist('pop_loadeep')
        p = which('eegplugin_limo.m');
        p = p(1:findstr(p,'eegplugin_limo.m')-1);
        addpath( p );
    end;
    
% create menus
% ------------
menu=findobj(fig,'tag','tools');
submenu = uimenu(menu, 'Label', version, 'separator', 'on', 'userdata','startup:on,study:on');
uimenu(submenu,'Label', 'GUI','callback','limo_eeg', 'userdata','study:on');
uimenu(submenu,'Label', '1st level analysis','callback','limo_batch', 'userdata','study:off');
uimenu(submenu,'Label', '2nd level analysis','callback','limo_random_effect','userdata','study:on,epoch:off,continuous:off');
uimenu(submenu,'Label', 'LIMO EEG results','callback','limo_results','userdata','study:on');
uimenu(submenu,'Label', 'LIMO EEG tools','callback','limo_tools','userdata','study:on');  