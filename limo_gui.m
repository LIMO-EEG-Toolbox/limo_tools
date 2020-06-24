function varargout = limo_gui(varargin)

% GUI of LIMO_eeg toolbox
% created using GUIDE
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

% --- Executes on button press in import
function import_Callback(hObject, eventdata, handles)
uiresume
delete(handles.figure1)
limo_batch('model specification');

% --- Executes on button press in analyse.
function analyse_Callback(hObject, eventdata, handles)
[file,dir_path] = uigetfile('LIMO.mat','select a LIMO.mat file');
[Names,Paths,Files] = limo_get_files([],{'*.mat;*.txt'},'Select a LIMO.mat file or a list of such files');
if isempty(Names)
    return
    guidata(hObject, handles);
else
    for f=1:size(Files,2)
        load(Files{f});
        if LIMO.Level ~=1
            fprintf('cannot reprocess 2nd level LIMO.mat\n %s\n',Paths{f})
        else
            cd(LIMO.dir)
            LIMO.design.status = 'to do';
            save('LIMO','LIMO');
            fprintf('reprocessing \n %s\n',Paths{f})
            limo_eeg(4);
            if isfield(LIMO,'contrast')
                limo_eeg(6);
            end
        end
        uiresume; delete(handles.figure1);
        cd (dir_path); limo_results;
    end
end

% --- Executes on button press in results.
function results_Callback(hObject, eventdata, handles)
uiresume
delete(handles.figure1)
limo_results;

% --- Executes on button press in contrasts.
function contrasts_Callback(hObject, eventdata, handles)

clc; uiresume
guidata(hObject, handles);
opt = questdlg('run constrasts for all subjects or open one subject?','option','all','one','all');
delete(handles.figure1)
if strcmp(opt,'one')
    limo_contrast_manager
else
    limo_batch('contrast only');
end


% [file,dir_path] = uigetfile('*.mat','select a LIMO.mat file');
% if file ==0
%     return
%     guidata(hObject, handles);
% else
%     cd(dir_path); handles.LIMO = load('LIMO.mat');
%     uiresume
%     delete(handles.figure1)
%     limo_contrast_manager(handles.LIMO.LIMO);
% end

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
    parpool('close');
end
delete(handles.figure1)
