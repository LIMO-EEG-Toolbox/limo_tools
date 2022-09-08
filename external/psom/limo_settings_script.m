limo_settings.psom = true; % toggle the use of the PSOM library
limo_settings.workdir = ''; %'derivatives'; % empty is current folder

if isequal(limo_settings.workdir, 'derivatives')
    try
        S=evalin('base','STUDY');
        limo_settings.workdir = fullfile(S.filepath, 'derivatives');
    catch
    end
end
