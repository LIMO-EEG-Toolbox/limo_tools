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
gui_State = struct('gui_Name',       mfilename, ...
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
handles.b    = 1000;
handles.tfce = 0;
handles.ica  = 0;
handles.dir  = [];
try S=evalin('base','STUDY');
    handles.chan_file = S.design.limo.chanloc; clear S
    fprintf('using study channel location file \n%s\n',handles.chan_file);
catch
    try
        S=evalin('base','gp_level_chanlocs');
        handles.chan_file = S;
    catch
        handles.chan_file = [];
    end
end
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
        warndlg(['Our simulations suggest that you need at leat 1000 bootstraps, consider changing your value: current boot ' num2str(handles.b)],'Bootstrap issue');
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
M = get(hObject,'Value');
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
    answer = questdlg('What ANOVA model do you want to run?', 'Model selection', 'Repeated Measures ANOVA', ...
        'N-Ways ANOVA','ANCOVA','Repeated Measures ANOVA');
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
            stop = questdlg('a LIMO file already exists, are you sure you want to move to this directory', ...
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
    test = test.(cell2mat(fieldnames(test)));
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
        warndlg('this file is not recognize as a channel location file or informations are missing','file error')
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

clc
uiresume
guidata(hObject, handles);
delete(handles.figure1)
limo_gui


% subfunction called before calling the others 
% to test chanlocs is loaded
function go = test_chan_loc(handles)

if isempty(handles.chan_file)
    go = 0;
    warndlg('chanloc not specified, please load one','missing file')
else
    go = 1;
end

function go = update_dir(handles,test)

go = 0;
if isempty(handles.dir)
    if ~exist(test,'dir')
        disp('directory not specified, creating one')
        mkdir(test); cd(test); handles.dir = pwd; go = 1;
    else
        warndlg2(sprintf('directory not specified, %s already exists \n please create and select a directory',test),'directory issue')
    end
else
    go = 1;
end



