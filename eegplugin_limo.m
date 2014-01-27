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

global EEG

version ='LInear MOdeling 1.5';
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
submenu = uimenu(menu, 'Label', version, 'separator', 'on');
uimenu(submenu,'Label', 'ERP analyzes','callback','limo_eeg');
uimenu(submenu,'Label', 'Time Frequency analyzes','callback','limo_eeg_tf');
uimenu(submenu,'Label', 'Component Analyzes','callback','limo_egg_cp');


    