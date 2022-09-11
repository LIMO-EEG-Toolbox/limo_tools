limo_settings.psom = true; % toggle the use of the PSOM library
limo_settings.workdir = ''; % empty is current folder
limo_settings.newgui = false; % empty is current folder

if exist('limo_settings_script_user')
    eval('limo_settings_script_user')
end

STUDY = [];
if isequal(limo_settings.workdir, 'derivatives')
    try
        STUDY=evalin('base','STUDY');
        %limo_settings.workdir = fullfile(S.filepath, 'derivatives');
        limo_settings.workdir = STUDY.filepath;
    catch
    end
end
