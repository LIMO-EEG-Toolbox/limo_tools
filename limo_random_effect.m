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
gui_Singleton = 1;
gui_State = struct('gui_Name',       'limo_random_effect', ...
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
set(hObject, 'DockControls', 'off', 'WindowStyle', 'normal');
handles.output = hObject;
guidata(hObject, handles);

% define handles used for the save callback
handles.b         = 0;
handles.tfce      = 0;
handles.type      = 'Channels';
handles.dir       = [];
handles.chan_file = [];
handles           = get_chan_loc(handles);
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

if strcmpi(handles.type,'Channels')
    go = test_chan_loc(handles);
    if go
        if isempty(handles.dir)
            go = update_dir(handles,'Summary_stat');
            if go == 0
                return
            end
        end
        uiresume
        guidata(hObject, handles);
        delete(handles.figure1)
        limo_central_tendency_and_ci(handles.chan_file);

    end
else
    disp('Currently only supporting channel anaysis type')
end


% Robust estimates and CI of differences
% ---------------------------------------------------------------
function differences_Callback(hObject, ~, handles)

if strcmpi(handles.type,'Channels')
    go = test_chan_loc(handles);
    if go
        if isempty(handles.dir)
            go = update_dir(handles,'Difference_stat');
            if go == 0
                return
            end
        end
        uiresume
        guidata(hObject, handles);
        delete(handles.figure1)
        limo_central_tendency_and_ci(handles.chan_file);
    end
    set(hObject, 'Visible', 'off'); 
    limo_plot_difference;
    set(hObject, 'Visible', 'on'); 
    guidata(hObject, handles);
else
    disp('Currently only supporting channel anaysis type')
end
guidata(hObject, handles);

%-------------------------------
%    Parameters 
%------------------------------

function Bootstrap_Checkbox_Callback(hObject, ~, handles)

if get(hObject, 'value') % tick box on
    handles.b = str2double(get(findobj(hObject.Parent, 'tag', 'bootstrap'),'String'));
    set(findobj(hObject.Parent, 'tag', 'TFCE'),'enable','on');
    set(findobj(hObject.Parent, 'tag', 'bootstrap'),'enable','on');
else
    handles.b = 0;
    set(findobj(hObject.Parent, 'tag', 'TFCE'),'value',0);
    set(findobj(hObject.Parent, 'tag', 'TFCE'),'enable','off');
    set(findobj(hObject.Parent, 'tag', 'bootstrap'),'enable','off');
end
guidata(hObject, handles);

% get the number of bootstraps
% ---------------------------------------------------------------
function handles = bootstrap_Callback(hObject, ~, handles)

handles.b = str2double(get(hObject,'String'));
if isempty(handles.b)
    handles.b = 0;
elseif handles.b == 0
    disp('bootstrap set to null')
    set(handles.TFCE,'Value',0)
else
    fprintf('bootstrap changed to %g \n',handles.b)
    if handles.b > 0 && handles.b < 1000
        if handles.b ~= 101
            limo_warndlg(['Our simulations suggest that you need at least 1000 bootstraps, your current value ' num2bouble(handles.b) ' may be overwritten'],'Bootstrap issue');
        else
            limo_warndlg('Using boostrap 101 - a trick to test bootstrap functionality, but remember it is not valid')
        end
    end
end
guidata(hObject, handles);

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

%-------------------------
%         TESTS_PANEL
%------------------------

% Analysis type
% ---------------------------------------------------------------
function Analysis_type_Callback(hObject, ~, handles)

M = get(hObject,'Value');
if M == 1
    handles.type = 'Channels';
elseif M == 2
    handles.type = 'Components';
elseif M == 3
    handles.type = 'Sources';
end
guidata(hObject, handles);

% One_Sample_t_test
% ---------------------------------------------------------------
function One_Sample_t_test_Callback(~, ~, handles)

limo_settings_script; % set STUDY and limo_settings
if  ~isempty(handles.dir)
    cd(handles.dir)
else
    if ~isempty(limo_settings.workdir)
        cd(limo_settings.workdir);
    end
end

go = update_dir(handles,'One_Sample_Ttest');
if go == 0
    return
else
    if strcmpi(handles.type,'Channels')
        go = test_chan_loc(handles);
        if go
            limo_random_select('one sample t-test',handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type',handles.type);
        end
    else
        limo_random_select('one sample t-test',[],'nboot',handles.b,'tfce',handles.tfce,'type',handles.type);
    end
end
if ~isempty(limo_settings.workdir)
    cd(limo_settings.workdir);
end

% Two_Samples_t_test
% ---------------------------------------------------------------
function Two_Samples_t_test_Callback(~, ~, handles)

limo_settings_script; % set STUDY and limo_settings
if  ~isempty(handles.dir)
    cd(handles.dir)
else
    if ~isempty(limo_settings.workdir)
        cd(limo_settings.workdir);
    end
end

go = update_dir(handles,'Two_Samples_Ttest');
if go == 0
    return
else
    if strcmpi(handles.type,'Channels')
        go = test_chan_loc(handles);
        if go
            limo_random_select('two-samples t-test',handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type',handles.type);
        end
    else
        limo_random_select('two-samples t-test',[],'nboot',handles.b,'tfce',handles.tfce,'type',handles.type);
    end
end


% Paired_t_test
% ---------------------------------------------------------------
% --- Executes on button press in Paired_t_test.
function Paired_t_test_Callback(~, ~, handles)

limo_settings_script; % set STUDY and limo_settings
if  ~isempty(handles.dir)
    cd(handles.dir)
else
    if ~isempty(limo_settings.workdir)
        cd(limo_settings.workdir);
    end
end

go = update_dir(handles,'Paired_Samples_Ttest');
if go == 0
    return
else
    if strcmpi(handles.type,'Channels')
        go = test_chan_loc(handles);
        if go
            limo_random_select('paired t-test',handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type',handles.type);
        end
    else
        limo_random_select('paired t-test',[],'nboot',handles.b,'tfce',handles.tfce,'type',handles.type);
    end
end

% Regression
% ---------------------------------------------------------------
% --- Executes on button press in Regression.
function Regression_Callback(~, ~, handles)

limo_settings_script; % set STUDY and limo_settings
if  ~isempty(handles.dir)
    cd(handles.dir)
else
    if ~isempty(limo_settings.workdir)
        cd(limo_settings.workdir);
    end
end

go = update_dir(handles,'Regression');
if go == 0
    return
else
    if strcmpi(handles.type,'Channels')
        go = test_chan_loc(handles);
        if go
            uiresume
            limo_random_select('regression',handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type',handles.type);
        end
    else
        uiresume
        limo_random_select('regression',[],'nboot',handles.b,'tfce',handles.tfce,'type',handles.type);
    end
end


% ANOVA/ANCOVA
% ---------------------------------------------------------------

% --- Executes on button press in ANOVA.
function ANOVA_Callback(~, ~, handles)

limo_settings_script; % set STUDY and limo_settings
if  ~isempty(handles.dir)
    cd(handles.dir)
else
    if ~isempty(limo_settings.workdir)
        cd(limo_settings.workdir);
    end
end

answer = limo_questdlg('Which of the following ANOVA models following do you want to apply to the data (bold value is the default)?', 'Model selection', ...
    '     N-Ways ANOVA     ','ANCOVA','Repeated Measures ANOVA','Repeated Measures ANOVA');

if strcmpi(answer,'Repeated Measures ANOVA')
    go = update_dir(handles,'Rep_Meas_ANOVA');
elseif strcmpi(answer,'     N-Ways ANOVA     ')
    go = update_dir(handles,'ANOVA');
elseif strcmpi(answer,'ANCOVA')
    go = update_dir(handles,'ANCOVA');
else
    return
end

if go == 0
    return
else
    if strcmpi(handles.type,'Channels')
        if strcmpi(handles.type,'Channels')
            go = test_chan_loc(handles);
            if go
                limo_random_select(answer,handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type',handles.type);
            end
        else
            limo_random_select(answer,[],'nboot',handles.b,'tfce',handles.tfce,'type',handles.type);
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
function handles = chan_cluster_neighbours_Callback(hObject, ~, handles)

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


% --- Executes on button press in close
% ---------------------------------------------------------------
function Quit_Callback(hObject, ~, handles)

uiresume
guidata(hObject, handles);
delete(handles.figure1)
limo_gui

% --- Executes on button press plot
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

% ----------------------
% subfunction to find channel locations
% ----------------------
function handles = get_chan_loc(handles)

limo_settings_script; % set STUDY and limo_settings

if strcmpi(handles.type,'channels')
    chan_file = fullfile(limo_settings.workdir, 'limo_gp_level_chanlocs.mat');
    if exist(chan_file, 'file')
        handles.chan_file = chan_file;
    end
    if ~isempty(STUDY) && isfield(STUDY, 'limo') && isfield(STUDY.limo, 'chanloc')
        if ~isempty(STUDY.limo.chanloc)
            handles.chan_file = STUDY.limo.chanloc;
            warning('using study channel location from STUDY');
        end
    end
end


% -----------------------------------------------------------------------
% subfunction called before calling the others to test chanlocs is loaded
% -----------------------------------------------------------------------
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

if isempty(handles.dir) && ...
        strcmpi(fullfile(fileparts(pwd), test),pwd) % already at the right place
    if exist(fullfile(pwd,'LIMO.mat'),'file')
        res = limo_questdlg(sprintf('The directory "%s" already contains stat files, do you want to overwrite them?',test), ...
            'Directory containing results', 'Cancel', 'Overwrite', 'Overwrite');
        if strcmpi(res, 'Cancel')
            go = 0;
            return;
        else
            handles.dir = pwd; go = 1; %#ok<*NASGU>
        end
    else
        handles.dir = pwd; go = 1;
    end
    return
end

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
    mkdir(test); cd(test); handles.dir = pwd; go = 1;

else
    if exist(fullfile(handles.dir,'LIMO.mat'),'file')
        res = limo_questdlg(sprintf('The directory "%s" already contains stat files, do you want to overwrite them?',handles.dir), ...
            'Directory containing results', 'Cancel', 'Overwrite', 'Overwrite');
        if strcmpi(res, 'cancel')
            go = 0;
            return;
        else
            go = 1;
        end
    else
        go = 1;
    end
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
if exist(limo_settings.workdir,'dir')
    cd(limo_settings.workdir);
end

