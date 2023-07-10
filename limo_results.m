function varargout = limo_results(varargin)

% Result GUI for the LIMO_eeg toolbox
% Created using GUIDE
% Cyril Pernet 20-03-2009 v1
% ------------------------------
%  Copyright (C) LIMO Team 2019


%% GUI stuffs
% -------------------------
% Begin initialization code
% -------------------------
warning off
gui_Singleton = 1;
gui_State = struct('gui_Name','limo_results', ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @limo_results_OpeningFcn, ...
    'gui_OutputFcn',  @limo_results_OutputFcn, ...
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
function limo_results_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

% define handles used for the save callback
handles.p         = 0.05;
handles.MCC       = 1;
handles.dir       = pwd;
handles.bootstrap = 0;
handles.tfce      = 0;
handles.filter    = {'*.mat'};
guidata(hObject, handles);
%uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = limo_results_OutputFcn(hObject, eventdata, handles)
varargout{1} = 'LIMO result terminated';


%% Callbacks

%------------------------
%         VISUALIZE
%------------------------

% show results as image
% ---------------------------------------------------------------
function Image_results_Callback(hObject, eventdata, handles)

[FileName,PathName,FilterIndex]=limo_get_result_file;
if FilterIndex ~= 0
    cd(PathName);
    if contains(FileName,'tfce')
        if exist(fullfile(fileparts(PathName(1:end-1)),FileName(6:end)),'file')
            [PathName,FileName] = fileparts(fullfile(fileparts(PathName(1:end-1)),FileName(6:end)));
            FileName = [FileName '.mat'];
        else
           error('parent file %s not found',fullfile(fileparts(PathName(1:end-1)),FileName(6:end))) 
        end
    end
    tmp          = load([PathName filesep 'LIMO.mat']);
    handles.LIMO = tmp.LIMO; clear tmp
    limo_display_results(1,FileName,PathName,handles.p,handles.MCC,handles.LIMO);
end
guidata(hObject, handles);

% Topoplot
% ---------------------------------------------------------------
function Topoplot_Callback(hObject, eventdata, handles)

[FileName,PathName,FilterIndex]=limo_get_result_file;
if FilterIndex == 1
    cd(PathName);
    if contains(FileName,'tfce')
        if exist(fullfile(fileparts(PathName(1:end-1)),FileName(6:end)),'file')
            [PathName,FileName] = fileparts(fullfile(fileparts(PathName(1:end-1)),FileName(6:end)));
            FileName = [FileName '.mat'];
        else
           error('parent file %s not found',fullfile(fileparts(PathName(1:end-1)),FileName(6:end))) 
        end
    end
    tmp          = load([PathName filesep 'LIMO.mat']);
    handles.LIMO = tmp.LIMO; clear tmp
    limo_display_results(2,FileName,PathName,handles.p,handles.MCC,handles.LIMO);
end
guidata(hObject, handles);

% course plots
% ---------------------------------------------------------------
function ERP_Callback(hObject, eventdata, handles)

[FileName,PathName,FilterIndex]=uigetfile('*.mat','Select LIMO file or summary stat file');
if FilterIndex == 1
    if ~strcmpi(FileName,'LIMO.mat')
        if contains(FileName,'tfce')
            FileName = FileName(6:end);
            PathName = fileparts(fileparts(PathName)); % trailing \
        end
        cd(PathName);
    end
    tmp = load([PathName filesep 'LIMO.mat']);
    handles.LIMO = tmp.LIMO; clear tmp
    limo_display_results(3,FileName,PathName,handles.p,handles.MCC,handles.LIMO);
end
guidata(hObject, handles);

% Review Design
% --------------
function review_design_Callback(hObject, eventdata, handles)
limo_review
guidata(hObject, handles);

function review_design_CreateFcn(hObject, eventdata, handles)

% get the p value threshold for display
% ---------------------------------------------------------------
function p_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function p_value_Callback(hObject, eventdata, handles)

handles.p = str2double(get(hObject,'String'));
test = isempty(handles.p);
if test == 1
    handles.p = 0.05;
end
guidata(hObject, handles);

% get the multiple comparisons correction method
% ---------------------------------------------------------------
function MCC_Choice_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MCC_Choice_Callback(hObject, eventdata, handles)

handles.MCC = get(hObject,'Value');  % 1 = None, 2 = Cluster, 3 = TFCE, 4 = T max
test        = isempty(handles.MCC);
if test == 1
    handles.MCC = 1;
end
guidata(hObject, handles);

        
%-------------------------
%         NEW ANALYSES
%------------------------

% --- New contrast.
% ---------------------------------------------------------------
function New_Contrast_Callback(hObject, eventdata, handles)

clc; uiresume
guidata(hObject, handles);
delete(handles.figure1)
limo_contrast_manager

% Semi-Partial coef
% ---------------------------------------------------------------
function Partial_coef_Callback(hObject, eventdata, handles)
[FileName,PathName,FilterIndex]=uigetfile('LIMO.mat','Select a LIMO file');
if FileName ~=0
    cd(PathName); handles.LIMO = load('LIMO.mat');
    limo_semi_partial_coef(handles.LIMO.LIMO)
end
guidata(hObject, handles);

%------------------------
%         OTHERS
%------------------------

% --- Executes on button press in Help.
% ---------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
web(fullfile(fileparts(which('limo_eeg')),[filesep 'help' filesep 'limo_results.html']));
guidata(hObject, handles);

% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)
clc; uiresume
guidata(hObject, handles);
delete(handles.figure1)
limo_gui

