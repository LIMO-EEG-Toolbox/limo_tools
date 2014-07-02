function varargout = limo_itc_gui2(varargin)

% ITC INTERFACE 2 - Test selection
% created using GUIDE 
% Based on limo_batch_gui
%
% Andrew Stewart, May 2014
% -----------------------------
% Copyright (C) LIMO Team 2014


%% GUI stuff
% -------------------------
% Begin initialization code
% -------------------------
warning off

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @limo_itc_gui2_OpeningFcn, ...
                   'gui_OutputFcn',  @limo_itc_gui2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% -----------------------
% End initialization code
% -----------------------


% --------------------------------------------------
%   Executes just before the menu is made visible
% --------------------------------------------------
function limo_itc_gui2_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% define handles used for the save callback
handles.test_select = 0;

guidata(hObject, handles);
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = limo_itc_gui2_OutputFcn(hObject, eventdata, handles) 

guidata(hObject, handles);

varargout{1} = handles.test_select;

delete(handles.figure1)

%% Callbacks


%-------------------------
%         OTHERS
%------------------------

% --- Executes on button press in Done.
% ---------------------------------------------------------------
function Done_Callback(hObject, eventdata, handles)
  
defaults.test_select       = handles.test_select;

uiresume
guidata(hObject, handles);




% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)


clc
uiresume

guidata(hObject, handles);

limo_gui


% --- Executes on button press in One_Sample_t_test.
function One_Sample_t_test_Callback(hObject, eventdata, handles)
% hObject    handle to One_Sample_t_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.test_select = {1,'One_Sample_t_test'};

uiresume
guidata(hObject, handles);


% --- Executes on button press in Two_Samples_t_test.
function Two_Samples_t_test_Callback(hObject, eventdata, handles)
% hObject    handle to Two_Samples_t_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.test_select = {2,'Two_Samples_t_test'};

uiresume
guidata(hObject, handles);


% --- Executes on button press in Paired_t_test.
function Paired_t_test_Callback(hObject, eventdata, handles)
% hObject    handle to Paired_t_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.test_select = {3,'Paired_t_test'};

uiresume
guidata(hObject, handles);


% --- Executes on button press in Regression.
function Regression_Callback(hObject, eventdata, handles)
% hObject    handle to Regression (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.test_select = {4,'Regression'};

uiresume
guidata(hObject, handles);


% --- Executes on button press in ANOVA.
function ANOVA_Callback(hObject, eventdata, handles)
% hObject    handle to ANOVA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.test_select = {5,'ANOVA'};

uiresume
guidata(hObject, handles);
