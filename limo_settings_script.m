% Script allowing to custimize LIMO tools behaviour
%
% Arnaud Delorme & Cyril Pernet

% new GUI available for plotting and some selection menus
limo_settings.newgui = true; 

% default place where data are written is into 'derivatives'
% set to empty to work in the current folder
limo_settings.workdir = 'derivatives'; 

% toggle the use of the PSOM library, use to batch over subjects as a
% pipeline, ie capture errors and restarts where it failed
limo_settings.psom = true; 

% overwrite using your own script
if exist('limo_settings_script_user')
    eval('limo_settings_script_user')
end

%% --------------------------------------------------
%                  do not edit
% ---------------------------------------------------
if isequal(limo_settings.workdir, 'derivatives')
    try
        STUDY=evalin('base','STUDY');
        limo_settings.workdir = fullfile(STUDY.filepath, 'derivatives');
    catch
        disp('Failed to get study');
        limo_settings.workdir = '';
        if ~exist('STUDY', 'var')
            STUDY = [];
        end
    end
else
    STUDY = [];
end
