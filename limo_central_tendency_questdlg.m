function varargout = limo_central_tendency_questdlg(varargin)
% LIMO_CENTRAL_TENDENCY_QUESTDLG M-file for limo_central_tendency_questdlg.fig

% GUI of LIMO_eeg toolbox 
% Called by limo_central_tendency_and)ci.m
% created using GUIDE
% cyril pernet 12-03-2010 v1
% -----------------------------
%  Copyright (C) LIMO Team 2010

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @limo_central_tendency_questdlg_OpeningFcn, ...
                   'gui_OutputFcn',  @limo_central_tendency_questdlg_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before limo_central_tendency_questdlg is made visible.
function limo_central_tendency_questdlg_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for limo_central_tendency_questdlg
handles.output = hObject;
handles.Estimator1 = [];
handles.Estimator2 = [];

% Update handles structure
guidata(hObject, handles);
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = limo_central_tendency_questdlg_OutputFcn(hObject, eventdata, handles) 

if isempty(handles.Estimator1) && isempty(handles.Estimator2)
    varargout{1} = [] ;
    varargout{2} = [] ;
else
    varargout{1} = handles.Estimator1;
    varargout{2} = handles.Estimator2;
end
delete(handles.figure1)


% -------------------------------------------------------------------------
% Intra-subject panel
% -------------------------------------------------------------------------

% --- Executes on button press in mean_1.
function mean1_Callback(hObject, eventdata, handles)

status = get(hObject,'Value');
if status == 1
    handles.Estimator1 = 'Mean';
end
guidata(hObject, handles);


% --- Executes on button press in trimmed_mean1.
function trimmed_mean1_Callback(hObject, eventdata, handles)

status = get(hObject,'Value');
if status == 1
    handles.Estimator1 = 'Trimmed mean';
end
guidata(hObject, handles);


% --- Executes on button press in harrell_davis1.
function harrell_davis1_Callback(hObject, eventdata, handles)

status = get(hObject,'Value');
if status == 1
    handles.Estimator1 = 'HD';
end
guidata(hObject, handles);


% --- Executes on button press in median1.
function median1_Callback(hObject, eventdata, handles)

status = get(hObject,'Value');
if status == 1
    handles.Estimator1 = 'Median';
end
guidata(hObject, handles);


% -------------------------------------------------------------------------
% Inter-subject panel
% -------------------------------------------------------------------------

% --- Executes on button press in mean2.
function mean2_Callback(hObject, eventdata, handles)
status = get(hObject,'Value');
if status == 1
    handles.Estimator2 = 'Mean';
end
guidata(hObject, handles);

% --- Executes on button press in trimmed_mean2.
function trimmed_mean2_Callback(hObject, eventdata, handles)
status = get(hObject,'Value');
if status == 1
    handles.Estimator2 = 'Trimmed mean';
end
guidata(hObject, handles);


% --- Executes on button press in harrell_davis2.
function harrell_davis2_Callback(hObject, eventdata, handles)
status = get(hObject,'Value');
if status == 1
    handles.Estimator2 = 'HD';
end
guidata(hObject, handles);


% --- Executes on button press in median2.
function median2_Callback(hObject, eventdata, handles)
status = get(hObject,'Value');
if status == 1
    handles.Estimator2 = 'Median';
end
guidata(hObject, handles);



% -------------------------------------------------------------------------
% Other buttons
% -------------------------------------------------------------------------


% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)
uiresume
guidata(hObject, handles);
clear all


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
origin = which('limo_eeg'); origin = origin(1:end-10); 
origin = sprintf('%shelp',origin); cd(origin)
web(['file://' which('limo_central_tendency_questdlg.html')]);
cd (handles.dir)


% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)

if ~isempty(handles.Estimator1) && ~isempty(handles.Estimator2)
    varargout{1} = handles.Estimator1;
    varargout{2} = handles.Estimator2;
    uiresume
    guidata(hObject, handles);
end
    
 
