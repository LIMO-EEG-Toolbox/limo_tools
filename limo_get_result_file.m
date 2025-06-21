function [FileName,PathName,FilterIndex]= limo_get_result_file

% routine to get result file(s) either from folder names or within a folder
% from file names
%
% FORMAT [FileName,PathName,FilterIndex]= limo_get_result_files
% OUTPUT FileName, PathName, Full File names are returned as string
%
% Arnaud Delorme & Cyril Pernet
% ------------------------------
%  Copyright (C) LIMO Team 2023

FilterIndex = 0;
FileName    = 0;
PathName    = 0;

limo_settings_script;
if limo_settings.newgui && ~isempty(limo_settings.workdir)
    if exist(limo_settings.workdir,'dir')
        cd(limo_settings.workdir);
    end
end

% level 1 or level 2
options         = { 'Level 1', 'Level 2'};
options_default = options{end};
if exist(fullfile(pwd,'LIMO.mat'),'file')
    LIMO = load(fullfile(pwd,'LIMO.mat'));
    if LIMO.LIMO.Level == 1
        options_default = options{1};
    end
end

res = limo_questdlg('Plot level 1 (subject) or level 2 (group) analysis file?', 'Result file', options{:}, options_default);
if isempty(res)
    return
else
    if contains(res, options{1}) % cancel
        level = 1;
        strGUI = 'Pick a result file for level 1 analysis (individual subject)';
    else
        level = 2;
        strGUI = 'pick a result file for level 2 analysis (group analysis)';
    end
end

if level == 1
    dirContent = dir('sub*/eeg/ses*/*/*.mat');
else
    % assuming default directory filenames
    % if scripted, also look for actual filenames, 
    % when the user if inside a directory
    dirContent1 = dir('One_Sample_Ttest*/*.mat');
    if isempty(dirContent1)
        dirContent1 = dir('One_Sample_Ttest*.mat');
    end
    if isempty(dirContent1)
        dirContent1 = dir('One_Sample_Ttest*/*/*.mat'); % since we can run many one samples at once these get nested
    end

    dirContent2 = dir('Paired_Samples_Ttest*/*.mat');
    if isempty(dirContent2)
        dirContent2 = dir('Paired_Samples_Ttest*.mat');
    end

    dirContent3 = dir('Two_Samples_Ttest*/*.mat');
    if isempty(dirContent3)
        dirContent3 = dir('Two_Samples_Ttest*.mat');
    end

    dirContent4 = dir('Regression*/*.mat');
    dirContent5 = [];
    if isempty(dirContent4)
        dirContent4 = dir('Regression*Covariate_effect_*.mat');
        dirContent5 = dir('Regression*con*.mat');
    end

    dirContent6 = dir('ANOVA*/*.mat');
    dirContent7 = [];
    if isempty(dirContent6)
        dirContent7 = dir('ANOVA*Condition_effect*.mat');
    end

    dirContent8 = dir('ANCOVA*/*.mat');
    dirContent9 = [];
    dirContent10 = [];
    if isempty(dirContent8)
        dirContent9 = dir('ANCOVA*Condition_effect*.mat');
        dirContent10 = dir('ANCOVA*Covariate_effect_*.mat');
    end

    dirContent11 = dir('Rep_Meas_ANOVA*/*.mat');
    dirContent12 = [];
    if isempty(dirContent11)
        dirContent12 = dir('Rep_Meas_ANOVA*Rep_ANOVA*.mat');
    end

    dirContent = [dirContent1;dirContent2;dirContent3;dirContent4;...
        dirContent5;dirContent6;dirContent7;dirContent8,...
        dirContent9,dirContent10,dirContent11,dirContent12];
end

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
    uiList = { {'style' 'text'       'string'  strGUI } ...
               { 'style' 'popupmenu' 'string' {dirContent.fullname}  'value', 1 } };
    res = inputgui('uilist', uiList, 'geometry', { [1] [1] }, 'cancel', 'Browse'); %#ok<NBRAK2>
    if ~isempty(res)
        FileName = dirContent(res{1}).name;
        PathName = dirContent(res{1}).folder;
        FilterIndex = 1;
    else
        [FileName,PathName,FilterIndex]=uigetfile('*.*','Select Result to plot');
    end
end    
