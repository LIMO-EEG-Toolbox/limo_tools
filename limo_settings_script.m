limo_settings.psom = true; % toggle the use of the PSOM library
limo_settings.workdir = 'derivatives'; % empty is current folder
limo_settings.newgui = true; % empty is current folder

if exist('limo_settings_script_user')
    eval('limo_settings_script_user')
end

if isequal(limo_settings.workdir, 'derivatives')
    try
        STUDY=evalin('base','STUDY');
        limo_settings.workdir = fullfile(STUDY.filepath, 'derivatives');
        %limo_settings.workdir = STUDY.filepath;
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
