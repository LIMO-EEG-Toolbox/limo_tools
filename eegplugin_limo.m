function version = eegplugin_limo(fig,try_strings,catch_strings)

% simple plugin function for eeglab
% call the LInear MOdeling of EEG data toolbox
% Multiple choices are available depending on the type
% of analysis
% 
% Nicolas Chauveau - Cyril Pernet v1 13/03/2009
% Cyril Pernet v2 22-06-2013
% -----------------------------
%  Copyright (C) LIMO Team 2019

version ='LIMO EEG 2.0';
if nargin < 3
    error('eegplugin_limo requires 3 arguments');
end

% add folder to path
% ------------------
if ~exist('pop_loadeep')
    p = which('eegplugin_limo.m');
    p = p(1:findstr(p,'eegplugin_limo.m')-1);
    addpath( p );
    addpath( fullfile(p, 'external') );
end

% create menus
% ------------
menu=findobj(fig,'tag','tools');
submenu = uimenu(menu, 'Label', version, 'separator', 'on', 'userdata','startup:on,study:on');
uimenu(submenu,'Label', 'GUI','callback','limo_eeg', 'userdata','study:on');
uimenu(submenu,'Label', '1st level analysis','callback','limo_batch', 'userdata','study:off');
uimenu(submenu,'Label', '2nd level analysis','callback','limo_random_effect','userdata','study:on,epoch:off,continuous:off');
uimenu(submenu,'Label', 'LIMO EEG results','callback','limo_results','userdata','study:on');
uimenu(submenu,'Label', 'LIMO EEG tools','callback','limo_tools','userdata','study:on');  

% create STUDY menus for EEGLAB 15
% --------------------------------
[~,vnum] = eeg_getversion;
if isempty(vnum) || vnum >= 15
    menu=findobj(fig,'tag','study');
    onstudy      = 'startup:off;epoch:off;continuous:off;study:on';
    nocheck      = 'try,';
    e_catch      = 'catch, eeglab_error; LASTCOM= ''''; clear EEGTMP ALLEEGTMP STUDYTMP; end;';
    e_hist       = [e_catch 'EEG = eegh(LASTCOM, EEG);'];
    e_plot_study = [e_catch 'if ~isempty(LASTCOM), STUDY = STUDYTMP; STUDY = eegh(LASTCOM, STUDY); disp(''Done.''); end; clear ALLEEGTMP STUDYTMP; eeglab(''redraw'');']; % ALLEEG not modified
    
    cb_limorunchan = [ nocheck '[STUDYTMP LASTCOM]= pop_limo(STUDY, ALLEEG,''dat'');' e_plot_study];
    cb_limoruncomp = [ nocheck '[STUDYTMP LASTCOM]= pop_limo(STUDY, ALLEEG,''ica'');' e_plot_study];
%     cb_limoreschan = [ nocheck 'pop_limoresults(STUDY,''dat'');' e_hist];
%     cb_limorescomp = [ nocheck 'pop_limoresults(STUDY,''ica'');' e_hist];
    limo_eeg = uimenu(menu,'Label','LInear MOdeling of EEG Data (BETA)' , 'userdata', onstudy,'Tag','limoeeg','separator', 'on');
    uimenu( limo_eeg, 'Label', 'Estimate Model Parameters (channnel)'        , 'userdata', onstudy, 'CallBack', cb_limorunchan);
    uimenu( limo_eeg, 'Label', 'Estimate Model Parameters (components)'     , 'userdata', onstudy, 'CallBack', cb_limoruncomp);
    uimenu( limo_eeg, 'Label', '2nd level analysis','callback','limo_random_effect', 'userdata', onstudy,'separator', 'on');
    uimenu( limo_eeg, 'Label', 'LIMO EEG results (1st and 2nd level)','callback','limo_results', 'userdata', onstudy);
    
end