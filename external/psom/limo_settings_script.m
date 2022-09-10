limo_settings.psom = false; % toggle the use of the PSOM library
limo_settings.workdir = ''; % empty is current folder

if exist('limo_settings_script_user')
    eval('limo_settings_script_user')
end

if isequal(limo_settings.workdir, 'derivatives')
    try
        S=evalin('base','STUDY');
        %limo_settings.workdir = fullfile(S.filepath, 'derivatives');
        limo_settings.workdir = S.filepath;
    catch
    end
end
