function varargout = limo_tools(varargin)

% Tools GUI for the LIMO_eeg toolbox
% Created using GUIDE 
% Cyril Pernet 25-08-2009 v1
% -----------------------------
%  Copyright (C) LIMO Team 2010

%% GUI stuffs
% -------------------------
% Begin initialization code
% -------------------------
warning off

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @limo_tools_OpeningFcn, ...
                   'gui_OutputFcn',  @limo_tools_OutputFcn, ...
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
function limo_tools_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

% define handles used for the save callback
handles.b = 1000;
handles.chan_file = [];
guidata(hObject, handles);
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = limo_tools_OutputFcn(hObject, eventdata, handles) 
varargout{1} = 'LIMO random effect terminated';


%% Callbacks



% --- Executes on button press in Review_Design.
% ---------------------------------------------------------------
function Review_Design_Callback(hObject, eventdata, handles)
limo_review
guidata(hObject, handles);

% --- Executes on button press in File_creation.
% ---------------------------------------------------------------
function File_creation_Callback(hObject, eventdata, handles)
limo_create_files
guidata(hObject, handles);

% --- Executes on button press in path_updates.
function path_updates_Callback(hObject, eventdata, handles)
% ---------------------------------------------------------------
limo_path_update
guidata(hObject, handles);

% --- Executes on button press in Expected_chanlocs.
% ---------------------------------------------------------------
function Expected_chanlocs_Callback(hObject, eventdata, handles)
limo_expected_chanlocs
guidata(hObject, handles);

% --- Executes on button press in electrode_optimization.
% ---------------------------------------------------------------
function electrode_optimization_Callback(hObject, eventdata, handles)
limo_best_electrodes
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function batch_CreateFcn(hObject, eventdata, handles)
% --- Executes on button press in electrode_optimization.
% ---------------------------------------------------------------
function batch_Callback(hObject, eventdata, handles)
uiresume
delete(handles.figure1)
limo_batch;
%guidata(hObject, handles);

% --- Executes on button press in CD.
% ---------------------------------------------------------------
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

% --- Executes on button press in Help.
% ---------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)

origin = which('limo_eeg'); origin = origin(1:end-10); 
origin = sprintf('%shelp',origin); cd(origin)
web(['file://' which('limo_tools.html')]);

% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)

clc
uiresume
guidata(hObject, handles);
delete(handles.figure1)
limo_gui






