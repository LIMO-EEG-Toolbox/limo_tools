function version = eegplugin_limo(fig,try_strings,catch_strings)

% simple plugin function for eeglab
% call the LInear MOdeling of EEG data toolbox
% Nicolas Chauveau - Cyril Pernet v1 13/03/2009
% -----------------------------
%  Copyright (C) LIMO Team 2010

global EEG

version ='LInear MOdeling 1.0';
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
    
% menu callback commands
%----------------------
FunctionImport=['limo_import(CURRENTSET)'];
FunctionModeling=0;
FunctionResults=0;

% create menus
% ------------
menu=findobj(fig,'tag','tools');
submenu = uimenu(menu, 'Label', version, 'separator', 'on');
uimenu(submenu,'Label', 'Import data','callback','limo_eeg(2)');
uimenu(submenu,'Label', 'Analyze','callback','limo_eeg(4)');
uimenu(submenu,'Label', 'Results','callback','limo_results');
uimenu(submenu,'Label', 'Gp effects','callback','limo_eeg(7)');
% uimenu(submenu,'Label', 'Model Selection','callback','limo_eeg(8)'); not implemented yet

