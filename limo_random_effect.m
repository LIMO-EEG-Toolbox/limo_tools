function varargout = limo_random_effect(varargin)

% Result GUI for the LIMO_eeg toolbox
% Created using GUIDE 
% Cyril Pernet 25-08-2009 v1
% ------------------------------
%  Copyright (C) LIMO Team 2019

%% GUI stuffs
% -------------------------
% Begin initialization code
% -------------------------
warning off

limo_settings_script;
if limo_settings.newgui
    guiName = [mfilename '_new'];
else
    guiName = mfilename;
end

gui_Singleton = 1;
gui_State = struct('gui_Name',       guiName, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @limo_random_effect_OpeningFcn, ...
                   'gui_OutputFcn',  @limo_random_effect_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% -----------------------
% End initialization code
% -----------------------


% --------------------------------------------------
%   Executes just before the menu is made visible
% --------------------------------------------------
function limo_random_effect_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

% define handles used for the save callback
handles.b    = 0;
handles.tfce = 0;
handles.ica  = 0;
handles.dir  = [];
handles = get_chan_loc(handles);
guidata(hObject, handles);
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = limo_random_effect_OutputFcn(~, ~, ~) 
varargout{1} = 'LIMO random effect terminated';


%% Callbacks

%---------------------------
%   BASIC_STATS_PANEL
%---------------------------

% Robust estimates and CI
% ---------------------------------------------------------------
function Central_tendency_and_CI_Callback(hObject, ~, handles)

if handles.ica == 1
    disp('IC not supported yet')
elseif test_chan_loc(handles)
    if isempty(handles.dir)
        savedir = uigetdir('pwd','where to save data?');
        if savedir == 0
            return
        else
            cd(savedir)
        end
    end
    limo_central_tendency_and_ci(handles.chan_file);
    guidata(hObject, handles);
end

% --- Executes on button press in data_plot.
function data_plot_Callback(hObject, ~, handles)

if handles.ica == 1
    disp('IC not supported yet')
elseif test_chan_loc(handles)
    limo_add_plots;
    guidata(hObject, handles);
end


% --- Executes on button press in differences.
function differences_Callback(hObject, ~, handles)

if handles.ica == 1
    disp('IC not supported yet')
elseif test_chan_loc(handles)
    limo_plot_difference
    guidata(hObject, handles);
end

% Parameters_box_plot
% ---------------------------------------------------------------
function Parameters_box_plot_Callback(hObject, ~, handles)

if handles.ica == 1
    disp('IC not supported yet')
elseif test_chan_loc(handles)
    limo_plots(handles.chan_file)
    guidata(hObject, handles);
end

 

%-------------------------------
%    Parameters 
%------------------------------

% get the number of bootstraps
% ---------------------------------------------------------------
function bootstrap_Callback(hObject, ~, handles)

handles.b = str2double(get(hObject,'String'));
if isempty(handles.b)
    handles.b = 1000;
elseif handles.b == 0
    disp('bootstrap set to null')
    set(handles.TFCE,'Value',0)
else
    fprintf('bootstrap changed to %g \n',handles.b)
    if handles.b > 0 && handles.b < 1000
        limo_warndlg(['Our simulations suggest that you need at leat 1000 bootstraps, consider changing your value: current boot ' num2str(handles.b)],'Bootstrap issue');
    end
end
guidata(hObject, handles);

function bootstrap_CreateFcn(hObject, ~, ~)
 
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in TFCE.
function TFCE_Callback(hObject, ~, handles)
M = get(hObject,'Value');
if M == 1
    handles.tfce = 1;
    disp('tfce is on');
elseif M == 0
    handles.tfce = 0;
    disp('tfce is off');
end
guidata(hObject, handles);


% --- Executes on button press in IC_analysis.
function IC_analysis_Callback(hObject, ~, handles)
M = get(hObject,enamValue');
if M == 1
    handles.ica = 1;
    disp('Analysis of ICs is on');
elseif M == 0
    handles.ica = 0;
    disp('Analysis of ICs is off');
end
guidata(hObject, handles);


%-------------------------
%         TESTS_PANEL
%------------------------


% One_Sample_t_test
% ---------------------------------------------------------------
function One_Sample_t_test_Callback(~, ~, handles)

go = update_dir(handles,'one_sample_ttest');
if go == 0
    return
else
    if handles.ica == 1
        limo_random_select('one sample t-test',handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Components');
    elseif test_chan_loc(handles)
        limo_random_select('one sample t-test',handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Channels');
    end
end

% Two_Samples_t_test
% ---------------------------------------------------------------
function Two_Samples_t_test_Callback(~, ~, handles)

go = update_dir(handles,'two_samples_ttest');
if go == 0
    return
else
if handles.ica == 1
    limo_random_select('two-samples t-test',handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Components');
elseif test_chan_loc(handles)
    limo_random_select('two-samples t-test',handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Channels');
end
end


% Paired_t_test
% ---------------------------------------------------------------
% --- Executes on button press in Paired_t_test.
function Paired_t_test_Callback(~, ~, handles)

go = update_dir(handles,'paired_ttest');
if go == 0
    return
else
if handles.ica == 1
    limo_random_select('paired t-test',handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Components');
elseif test_chan_loc(handles)
    limo_random_select('paired t-test',handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Channels');
end
end

% Regression
% ---------------------------------------------------------------
% --- Executes on button press in Regression.
function Regression_Callback(~, ~, handles)

go = update_dir(handles,'regression');
if go == 0
    return
else
if handles.ica == 1
    limo_random_select('regression',handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Components');
elseif test_chan_loc(handles)
    limo_random_select('regression',handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Channels');
end
end

% ANOVA/ANCOVA
% ---------------------------------------------------------------

% --- Executes on button press in ANOVA.
function ANOVA_Callback(~, ~, handles)

go = update_dir(handles,'AN(C)OVA');
if go == 0
    return
else
    answer = limo_questdlg('Which of the following ANOVA models following do you want to apply to the data (bold value is the default)?', 'Model selection', ...
        '     N-Ways ANOVA     ','ANCOVA','Repeated Measures ANOVA','Repeated Measures ANOVA');
    if ~isempty(answer)
        if handles.ica == 1
            limo_random_select(answer,handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Components');
        elseif test_chan_loc(handles)
            limo_random_select(answer,handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Channels');
        end
    end
end

%------------------------
%         OTHERS
%------------------------

 % --- Executes on button press in CD.
function CD_Callback(hObject, ~, handles)

PathName=uigetdir(pwd,'select LIMO working directory');
if PathName ~= 0
    list = dir(PathName);
    for i=1:length(list)
        if strcmp(list(i).name,'LIMO.mat')
            stop = limo_questdlg('a LIMO file already exists, are you sure you want to move to this directory', ...
                'WARNING','Yes','No','No');
            if strcmp(stop,'Yes')
                cd(PathName);
                handles.dir = PathName;
            else
                return
            end
        else
            cd(PathName);
            handles.dir = PathName;
        end
    end
end
guidata(hObject, handles);


% --- Executes on button press in chan_cluster_neighbours.
function chan_cluster_neighbours_Callback(hObject, ~, handles)

[chan_file,chan_path,sts]=uigetfile('expected_chanlocs.mat','Select channel location file');
if sts == 1
    test = load([chan_path chan_file]);
    if isfield(test,'expected_chanlocs')
        test = test.expected_chanlocs;
    end
    
    if isstruct(test) && ~isempty(test(1).labels) && ~isempty(test(1).theta) && ~isempty(test(1).radius) ...
            && ~isempty(test(1).X) && ~isempty(test(1).Y) && ~isempty(test(1).Z) && ~isempty(test(1).sph_theta) ...
             && ~isempty(test(1).sph_phi) && ~isempty(test(1).sph_radius) 
            
        handles.chan_file = [chan_path chan_file];
        disp('channel location loaded');
        assignin('base', 'gp_level_chanlocs', [chan_path chan_file])
         
    else
        limo_warndlg('this file is not recognize as a channel location file or informations are missing','file error')
    end
end
guidata(hObject, handles);


% --- Executes on button press in Help.
% ---------------------------------------------------------------
function Help_Callback(~, ~, ~)

web(['file://' which('limo_random_effect.html')]);


% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Quit_Callback(hObject, ~, handles)

uiresume
guidata(hObject, handles);
delete(handles.figure1)
%limo_gui

% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Plotting_Callback(hObject, ~, handles)

uiresume
guidata(hObject, handles);
[FileName,PathName,FilterIndex]=limo_get_result_file;
if FilterIndex ~= 0
    tmp          = load([PathName filesep 'LIMO.mat']);
    limo_display_results(1,FileName,PathName,0.05,1,tmp.LIMO);
end

%delete(handles.figure1)
limo_results

% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Bootstrap_Checkbox_Callback(hObject, ~, handles)

if get(hObject, 'value')
    handles.b = str2double(get(findobj(hObject.Parent, 'tag', 'bootstrap'),'String'));
    set(findobj(hObject.Parent, 'tag', 'TFCE'),'enable','on');
else
    handles.b = 0;
    set(findobj(hObject.Parent, 'tag', 'TFCE'),'value',0);
    set(findobj(hObject.Parent, 'tag', 'TFCE'),'enable','off');
end
disp(handles.b);

% Bootstrap_Checkbox_Callback

% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Contrast_Callback(hObject, ~, handles)

if 0
    limo_batch('contrast only');
else
    % below the code allows to use bootstrap and TFCE
    [Names,Paths,allFiles,txtFile] = limo_get_anova_files;
  
    if isempty(Names)
        disp('No file selected, abording')
        return
    end
    disp('Looking up the corresponding LIMO file');
    LIMOfile = fullfile(Paths{1}, Names{1});
    handles = limo_contrast_manager(LIMOfile);
    if isempty(handles) || isempty(handles.C)
        disp('No contrast, abording')
        return
    end
    nBoot = str2double(get(findobj(hObject.Parent, 'tag', 'bootstrap'),'String'));
    TFCE  = get(findobj(hObject.Parent, 'tag', 'TFCE'),'value');
    
    if isempty(txtFile) % second level
        levels = 2;
        options = {'1st and 2nd level', '2nd level only'};
        res = limo_questdlg( [ 'This is a 2nd (group) level contrast. Do you want to calculate the same' 10 ...
                               'contrast at the 1st-level to be able to plot it for individual subjects? ' ], '1st-level contrast', options{:}, options{end});
        if isempty(res)
            return;
        elseif strcmpi(res, options{1})
            levels = [1 2];
        end
    else
        levels = 1;
    end
    if any(levels == 2)
        % 2nd level
        limo_contrast_execute(LIMOfile, handles);
        
        if any(levels == 1)
            % 1st level from second level
            LIMO = load('-mat', LIMOfile);
            betaFiles = LIMO.LIMO.data.data;
            anovaVars = LIMO.LIMO.design.labels;
            for gp = 1:length(betaFiles) % mutliple groups
                allLIMOfiles = cellfun(@(x)strrep(x, 'Betas.mat', 'LIMO.mat'), betaFiles{gp}, 'uniformoutput', false);

                % read first file to find indices of ANOVA betas in Subject
                % betas (could be the same)
                LIMOsubject = load('-mat', allLIMOfiles{1});
                subjectVars = LIMOsubject.LIMO.design.labels;
                indList = [];
                for iBeta = 1:length(anovaVars)
                    indTmp = strmatch(anovaVars(iBeta).description, { subjectVars.description }, 'exact');
                    if length(indTmp) ~= 1
                        error('Issue looking up variables for first level')
                    end
                    indList = [indList indTmp]; 
                end

                % run contrast
                contrast.LIMO_files = allLIMOfiles;
                contrast.mat = zeros(1, length(subjectVars))
                contrast.mat(indList) = handles.C;
                limo_settings_script; % for STUDY var
                limo_batch('contrast only', [], contrast, STUDY);
            end
        end
    elseif any(levels == 1)
        % first level only
        contrast.LIMO_files = strrep(txtFile, 'Beta_', 'LIMO_');
        contrast.mat = handles.C;
        limo_settings_script; % for STUDY var
        limo_batch('contrast only', [], contrast, STUDY);
    end
end

% ----------------------
% subfunction to find channel locations
% ----------------------
function handles = get_chan_loc(handles)

try 
    S=evalin('base','STUDY');
    handles.chan_file = S.design.limo.chanloc; clear S
catch
    try
        S=evalin('base','gp_level_chanlocs');
        handles.chan_file = S;
    catch
        go_to_working_dir;
        tmpChanFile = fullfile(pwd, 'limo_gp_level_chanlocs.mat');
        if exist(tmpChanFile, 'file')
            handles.chan_file = tmpChanFile;
        else
            handles.chan_file = [];
        end
    end
end
if ~isempty(handles.chan_file)
    fprintf('using study channel location file \n%s\n',handles.chan_file);
end

% ----------------------
% subfunction called before calling the others 
% to test chanlocs is loaded
% ----------------------
function go = test_chan_loc(handles)

if isempty(handles.chan_file)
    go = 0;
    limo_warndlg('chanloc not specified, please load one','missing file')
else
    go = 1;
end

% ----------------------
% create folder if necessary
% ----------------------
function go = update_dir(handles,test)
go = 0;
if isempty(handles.dir)
    go_to_working_dir; % no effect if limo_settings.workdir is empty
    if exist(fullfile(pwd, test),'dir')
        count = 2;
        while exist([ test num2str(count)],'dir')
            count = count + 1;
        end
        newtest = [ test num2str(count)];
        res = limo_questdlg(sprintf('The directory "%s" already exist, do you want to overwrite it or\ncreate a new one named "%s"?\nIf you overwrite, previous results will be lost.',test,newtest), ...
            'Directory containing results', 'Cancel', 'Create new', 'Overwrite', 'Overwrite');
        if strcmpi(res, 'cancel')
            go = 0;
            return;
        elseif strcmpi(res, 'Create new')
            test = newtest;
        end

        if exist(test,'dir')
            disp('Removing folder');
            rmdir(test, 's');
        end
    end

    %limo_warndlg(sprintf('Creating "%s" directory.\nRemember it so you can plot results.\nNow you will select the type of analysis and the single trial analysis result file.',test),'Directory containing results')
    mkdir(test); cd(test); handles.dir = pwd; go = 1;
else
    go = 1;
end

% ----------------------
% go to the right folder
% ----------------------
function go_to_working_dir()
if ~exist('limo_settings_script', 'file')
    pathTmp = fileparts(which('limo_eeg'));
    addpath(fullfile(pathTmp, 'external', 'psom'));
end
limo_settings_script;
if ~isempty(limo_settings.workdir)
    cd(limo_settings.workdir);
end


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11
