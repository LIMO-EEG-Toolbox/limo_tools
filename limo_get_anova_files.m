function [Names,Paths,Files,txtFile] = limo_get_anova_files(varargin)

% routine to get multifiles from different directories
%
% FORMAT [Names,Paths,Files] = limo_get_anova_files(gp,filter,title)
%
% OUTPUT Names , Paths, Full File names are returned as cells
%

Names = {};
Paths = {};
Files = {};
txtFile = '';

limo_settings_script;
if ~isempty(limo_settings.workdir)
    cd(limo_settings.workdir);
end
dirContent = dir('AN(C)OVA*');

if isempty(dirContent)
    txt = 'For 2nd-level contrast, you need to run an ANOVA first. What do you want to do instead?'
    options = { 'Cancel     ', 'Run 1st-level contrast', 'Pick ANOVA folder' };
    res = limo_questdlg(txt, '2nd-level contrast', options{:}, options{end});

    if contains(res, options{1}) % cancel
        return
    elseif contains(res, options{2}) % first level
        [Names,Paths,Files,txtFile] = limo_get_result_files(varargin);
        return;
    else
        path = uigetdir('*.*','Pick the ANOVA folder');
        if isequal(path, 0)
            return;
        end
    end
end
if length(dirContent) > 1
    uiList = { {'style' 'text' 'string' 'Pick an ANOVA folder' } ...
        { 'style' 'popupmenu' 'string' {dirContent.name}  'value', length(dirContent) }};
    res = inputgui('uilist', uiList, 'geometry', { [1] [1] }, 'cancel', 'Browse');
    if ~isempty(res)
        path = fullfile(dirContent(res{1}).folder, dirContent(res{1}).name);
    else
        path = uigetdir('*.*','Pick the ANOVA folder');
        if isequal(path, 0)
            return;
        end
    end
end    
res = dir(fullfile(path, 'LIMO.mat'));
if isempty(res)
    limo_errordlg('LIMO.mat file not found in folder')
else
    Names = { 'LIMO.mat' };
    Paths = { path };
end