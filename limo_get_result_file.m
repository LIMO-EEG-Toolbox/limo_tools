function [FileName,PathName,FilterIndex]= limo_get_result_file

% routine to get result file
%
% FORMAT [FileName,PathName,FilterIndex]= limo_get_result_files
%
% OUTPUT FileName, PathName, Full File names are returned as string

Names = {};
Paths = {};
Files = {};
txtFile = '';
FilterIndex = 0;
FileName = 0;
PathName = 0;

limo_settings_script;
if limo_settings.newgui
    cd(limo_settings.workdir);
end
dirContent1 = dir('AN(C)OVA*/*.mat');
dirContent2 = dir('one_sample_ttest*/*.mat');
dirContent3 = dir('one_sample_ttest*/*/*.mat');
dirContent4 = dir('paired_ttest*/*.mat');
dirContent5 = dir('two_samples_ttest*/*.mat');
dirContent6 = dir('regression*/*.mat');
dirContent  = [dirContent1;dirContent2;dirContent3;dirContent4;dirContent5;dirContent6];

% remove Yr and LIMO files
for iFile = length(dirContent):-1:1
    if contains(dirContent(iFile).name, 'LIMO.mat') || contains(dirContent(iFile).name, 'Yr.mat') || ...
            contains(dirContent(iFile).name, 'Yhat.mat') || contains(dirContent(iFile).name, 'R2.mat') || ...
            contains(dirContent(iFile).name, 'Y1r.mat') || contains(dirContent(iFile).name, 'Y2r.mat') || ...
            contains(dirContent(iFile).name, 'Betas.mat')
        dirContent(iFile) = [];
    else
        relPath = strrep(dirContent(iFile).folder, pwd, './');
        dirContent(iFile).fullname = fullfile(relPath, dirContent(iFile).name);
        if contains(dirContent(iFile).fullname, 'ess_')
            dirContent(iFile).fullname = [ dirContent(iFile).fullname ' (contrast)' ];
        end
    end
end

if isempty(dirContent)
    txt = [ 'No result files were found. Do you want' 10 ...
        'to browse folders for a .mat result file?' ];
    options = { 'Cancel', 'Browse'};
    res = limo_questdlg(txt, 'Result file', options{:}, options{end});

    if contains(res, options{1}) % cancel
        return
    else
        [FileName,PathName,FilterIndex]=uigetfile('*.*','Select Result to plot');
    end
else
    uiList = { {'style' 'text' 'string' 'Pick a result file' } ...
        { 'style' 'popupmenu' 'string' {dirContent.fullname}  'value', length(dirContent) }};
    res = inputgui('uilist', uiList, 'geometry', { [1] [1] }, 'cancel', 'Browse');
    if ~isempty(res)
        FileName = dirContent(res{1}).name;
        PathName = dirContent(res{1}).folder;
        FilterIndex = 1;
    else
        [FileName,PathName,FilterIndex]=uigetfile('*.*','Select Result to plot');
    end
end    
