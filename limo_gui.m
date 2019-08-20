function varargout = limo_gui(varargin)

% GUI of LIMO_eeg toolbox
% created using GUIDE
% cyril pernet 12-03-2010 v1
% ------------------------------------------
% Copyright (C) LIMO Team 2014


%% GUI stuffs
% -------------------------
% Begin initialization code
% -------------------------
warning off

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @limo_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @limo_gui_OutputFcn, ...
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
function limo_gui_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.dir = pwd;
guidata(hObject, handles);
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = limo_gui_OutputFcn(hObject, eventdata, handles) 
varargout{1} = 'LIMO terminated';


%% Callbacks

% --- Executes on selection change in import_menu.
function import_menu_Callback(hObject, eventdata, handles)
handles.import = get(hObject,'Value');
uiresume
guidata(hObject, handles);
delete(handles.figure1)
limo_eeg(2,handles.import);


% --- Executes during object creation, after setting all properties.
function import_menu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in analyse.
function analyse_Callback(hObject, eventdata, handles)
[file,dir_path] = uigetfile('LIMO.mat','select a LIMO.mat file');
if file ==0
    return
    guidata(hObject, handles);
else
    load([dir_path file]);
    if strcmp(LIMO.design.status,'done');
        answer=questdlg('This design was already estimated, redo it?','Estimation','Yes','No','Yes');
    else
        answer = 'Yes';
    end
    
    if strcmp(answer ,'Yes')
        uiresume; delete(handles.figure1); 
        cd (dir_path); LIMO.design.status ='to do';
        save LIMO LIMO; clear LIMO; limo_eeg(4); limo_results;
    end
end

% --- Executes on button press in results.
function results_Callback(hObject, eventdata, handles)
uiresume
delete(handles.figure1)
limo_results;

% --- Executes on button press in contrasts.
function contrasts_Callback(hObject, eventdata, handles)
[file,dir_path] = uigetfile('*.mat','select a LIMO.mat file');
if file ==0
    return
    guidata(hObject, handles);
else
    cd(dir_path); handles.LIMO = load('LIMO.mat');
    uiresume
    delete(handles.figure1)
    limo_contrast_manager(handles.LIMO.LIMO);
end

% --- Executes on button press in Rdx.
function Rdx_Callback(hObject, eventdata, handles)
uiresume
delete(handles.figure1)
limo_random_effect;


% --- Executes on button press in limo_tools.
function limo_tools_Callback(hObject, eventdata, handles)
uiresume
delete(handles.figure1)
limo_tools;



% --- Executes on button press in Help.
% ---------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
origin = which('limo_eeg'); origin = origin(1:end-10); 
origin = sprintf('%shelp',origin); cd(origin)
web(['file://' which('limo_gui.html')]);
cd (handles.dir)


% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)
clc
uiresume
try
    matlabpool('close');
end
delete(handles.figure1)
