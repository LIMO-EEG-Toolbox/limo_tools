function varargout = limo_random_effect(varargin)

% Result GUI for the LIMO_eeg toolbox
% Created using GUIDE 
% Cyril Pernet 25-08-2009 v1
% -----------------------------
%  Copyright (C) LIMO Team 2015

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
function limo_random_effect_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

% define handles used for the save callback
handles.b = 1000;
handles.tfce = 0;
handles.ica  = 0;
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
function varargout = limo_random_effect_OutputFcn(hObject, eventdata, handles) 
varargout{1} = 'LIMO random effect terminated';


%% Callbacks

%-------------------------
%         BASIC_STATS_PANEL
%------------------------

% Robust estimates and CI
% ---------------------------------------------------------------
function Central_tendency_and_CI_Callback(hObject, eventdata, handles)

if handles.ica == 1
    disp('IC not supported yet')
elseif test_chan_loc(handles)
    limo_central_tendency_and_ci(handles.chan_file);
    guidata(hObject, handles);
end


% --- Executes on button press in data_plot.
function data_plot_Callback(hObject, eventdata, handles)

if handles.ica == 1
    disp('IC not supported yet')
elseif test_chan_loc(handles)
    limo_add_plots;
    guidata(hObject, handles);
end


% --- Executes on button press in differences.
function differences_Callback(hObject, eventdata, handles)

if handles.ica == 1
    disp('IC not supported yet')
elseif test_chan_loc(handles)
    limo_plot_difference
    guidata(hObject, handles);
end

% Parameters_box_plot
% ---------------------------------------------------------------
function Parameters_box_plot_Callback(hObject, eventdata, handles)

if handles.ica == 1
    disp('IC not supported yet')
elseif test_chan_loc(handles)
    limo_plots(handles.chan_file)
    guidata(hObject, handles);
end

 

%-------------------------------
%         Parameters parameter
%------------------------------

% get the number of bootstraaps
% ---------------------------------------------------------------
function bootstrap_Callback(hObject, eventdata, handles)

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

function bootstrap_CreateFcn(hObject, eventdata, handles)
 
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in TFCE.
function TFCE_Callback(hObject, eventdata, handles)
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
function IC_analysis_Callback(hObject, eventdata, handles)
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
function One_Sample_t_test_Callback(hObject, eventdata, handles)

if handles.ica == 1
    limo_random_select(1,handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Components');
elseif test_chan_loc(handles)
    limo_random_select(1,handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Channels');
end


% Two_Samples_t_test
% ---------------------------------------------------------------
function Two_Samples_t_test_Callback(hObject, eventdata, handles)

if handles.ica == 1
    limo_random_select(2,handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Components');
elseif test_chan_loc(handles)
    limo_random_select(2,handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Channels');
end



% Paired_t_test
% ---------------------------------------------------------------
% --- Executes on button press in Paired_t_test.
function Paired_t_test_Callback(hObject, eventdata, handles)

if handles.ica == 1
    limo_random_select(3,handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Components');
elseif test_chan_loc(handles)
    limo_random_select(3,handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Channels');
end


% Regression
% ---------------------------------------------------------------
% --- Executes on button press in Regression.
function Regression_Callback(hObject, eventdata, handles)

if handles.ica == 1
    limo_random_select(4,handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Components');
elseif test_chan_loc(handles)
    limo_random_select(4,handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Channels');
end


% ANOVA/ANCOVA
% ---------------------------------------------------------------

% --- Executes on button press in ANOVA.
function ANOVA_Callback(hObject, eventdata, handles)

if handles.ica == 1
    limo_random_select(5,handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Components');
elseif test_chan_loc(handles)
    limo_random_select(5,handles.chan_file,'nboot',handles.b,'tfce',handles.tfce,'type','Channels');
end


%------------------------
%         OTHERS
%------------------------

 % --- Executes on button press in CD.
function CD_Callback(hObject, eventdata, handles)

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
function chan_cluster_neighbours_Callback(hObject, eventdata, handles)

[chan_file,chan_path,sts]=uigetfile('expected_chanlocs.mat','Select channel location file');
if sts == 1
    load ([chan_path chan_file])
    if exist('expected_chanlocs','var') == 1
        test = expected_chanlocs;
    else
        test = eval(chan_file(1:end-4));
    end
    
    if isstruct(test) && ~isempty(test(1).labels) && ~isempty(test(1).theta) && ~isempty(test(1).radius) ...
            && ~isempty(test(1).X) && ~isempty(test(1).Y) && ~isempty(test(1).Z) && ~isempty(test(1).sph_theta) ...
             && ~isempty(test(1).sph_phi) && ~isempty(test(1).sph_radius) && sum(channeighbstructmat(:)) ~= 0
             % && ~isempty(test(1).urchan) % urchan should not be needed
            
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
function Help_Callback(hObject, eventdata, handles)

origin = which('limo_eeg'); origin = origin(1:end-10); 
origin = sprintf('%shelp',origin); cd(origin)
web(['file://' which('limo_random_effect.html')]);


% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)

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



