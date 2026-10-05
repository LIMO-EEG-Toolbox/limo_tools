% Script allowing to custimize LIMO tools behaviour
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
        % A caller-supplied STUDY takes precedence over the GUI workspace.
        if ~exist('STUDY','var') || isempty(STUDY)
            settings_study = evalin('base','STUDY');
            if ~isstruct(settings_study) || isempty(settings_study)
                error('LIMO:MissingStudy', 'No valid STUDY available for the default output directory.');
            end
            STUDY = settings_study;
        end
        default_folder = [STUDY.filepath filesep 'derivatives'];
        if exist(fullfile(default_folder,['LIMO_' extractBefore(STUDY.filename,'.study')]),"dir")
            limo_settings.workdir = fullfile(default_folder,['LIMO_' extractBefore(STUDY.filename,'.study')]);
        else
            limo_settings.workdir = default_folder;
        end
    catch
        disp('Failed to find STUDY variable');
        limo_settings.workdir = '';
    end
end
clear settings_study
