function varargout = limo_create_single_trials_gui(varargin)

% GUI for limo_create_single_trials
% created using GUIDE
% cyril pernet 21-01-2015 v1
% ------------------------------------------
% Copyright (C) LIMO Team 2015


%% GUI stuffs
% -------------------------
% Begin initialization code
% -------------------------
warning off

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @limo_create_single_trials_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @limo_create_single_trials_gui_OutputFcn, ...
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
function limo_create_single_trials_gui_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output     = hObject;
handles.dir        = pwd;
handles.scalp_erp  = 'off';
handles.scalp_spec = 'off';
handles.scalp_ersp = 'off';
handles.scalp_itc  = 'off';
handles.ica_erp    = 'off';
handles.ica_spec   = 'off';
handles.ica_ersp   = 'off';
handles.ica_itc    = 'off';
handles.subjects   = [];
handles.MAT        = 0;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = limo_create_single_trials_gui_OutputFcn(hObject, eventdata, handles) 
varargout{1} = 'LIMO EEG create single trials terminated';


%% Callbacks

% --- Executes on button press in scalp_erp.
function scalp_erp_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
if h == 1; handles.scalp_erp = 'on'; else  handles.scalp_erp = 'off'; end
guidata(hObject, handles);

% --- Executes on button press in scalp_spec.
function scalp_spec_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
if h == 1; handles.scalp_spec = 'on'; else  handles.scalp_spec = 'off'; end
guidata(hObject, handles);

% --- Executes on button press in scalp_ersp.
function scalp_ersp_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
if h == 1; handles.scalp_ersp = 'on'; else  handles.scalp_ersp = 'off'; end
guidata(hObject, handles);

% --- Executes on button press in scalp_itc.
function scalp_itc_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
if h == 1; handles.scalp_itc = 'on'; else  handles.scalp_itc = 'off'; end
guidata(hObject, handles);

% --- Executes on button press in ica_erp.
function ica_erp_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
if h == 1; handles.ica_erp = 'on'; else  handles.scalp_erp = 'off'; end
guidata(hObject, handles);

% --- Executes on button press in ica_spec.
function ica_spec_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
if h == 1; handles.ica_spec = 'on'; else  handles.scalp_spec = 'off'; end
guidata(hObject, handles);

% --- Executes on button press in ica_ersp.
function ica_ersp_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
if h == 1; handles.ica_ersp = 'on'; else  handles.scalp_ersp = 'off'; end
guidata(hObject, handles);

% --- Executes on button press in ica_itc.
function ica_itc_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
if h == 1; handles.ica_itc = 'on'; else  handles.scalp_itc = 'off'; end
guidata(hObject, handles);

% -------------------------------------------------------------------------
% --- Executes on button press in MATsave.
function MATsave_Callback(hObject, eventdata, handles)
handles.MAT = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in load_data.
function load_data_Callback(hObject, eventdata, handles)
[~,~,handles.subjects]=limo_get_files([],{'*.txt;*.set;*.study'});
guidata(hObject, handles);

% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)

if isnumeric(handles.subjects) || isempty(handles.subjects)
    [~,~,handles.subjects]=limo_get_files([],{'*.txt;*.set;*.study'});
    guidata(hObject, handles);
end

if isnumeric(handles.scalp_erp);  handles.scalp_erp  = 'off'; end
if isnumeric(handles.scalp_spec); handles.scalp_spec = 'off'; end
if isnumeric(handles.scalp_ersp); handles.scalp_ersp = 'off'; end
if isnumeric(handles.scalp_itc);  handles.scalp_itc  = 'off'; end
if isnumeric(handles.ica_erp);    handles.ica_erp    = 'off'; end
if isnumeric(handles.ica_spec);   handles.ica_spec   = 'off'; end
if isnumeric(handles.ica_ersp);   handles.ica_ersp   = 'off'; end
if isnumeric(handles.ica_itc);    handles.ica_itc    = 'off'; end

if handles.MAT == 0
    format = 'cell';
else
    format = 'matrix';
end;

N = size(handles.subjects,2);
for s=1:N
    subject_set = handles.subjects{s};
    fprintf('processing subject %g/%g: %s\n',s,N,subject_set);
    
    if strcmp(handles.scalp_erp,'on') || strcmp(handles.scalp_spec,'on') || ...
            strcmp(handles.scalp_ersp,'on') || strcmp(handles.scalp_itc,'on')
        options = {'format',format, 'datatype', 'channels', ...
            'erp',handles.scalp_erp,'spec',handles.scalp_spec, ...
            'ersp',handles.scalp_ersp,'itc',handles.scalp_itc,...
            'rmicacomps','on', 'erpparams',[],'specparams',[],'erspparams',[], ...
            'interp','off', 'scalp','off','recompute','on','savetrials','on'};
        disp('creating single trials scalp data');
        limo_create_single_trials(subject_set,options{:})
    end
    
    if strcmp(handles.ica_erp,'on') || strcmp(handles.scalp_spec,'on') || ...
            strcmp(handles.scalp_ersp,'on') || strcmp(handles.scalp_itc,'on')
        options = {'format',format, 'datatype', 'ica', ...
            'erp',handles.ica_erp,'spec',handles.ica_spec, ...
            'ersp',handles.ica_ersp,'itc',handles.ica_itc,...
            'rmicacomps','on', 'erpparams',[],'specparams',[],'erspparams',[], ...
            'interp','off', 'scalp','off','recompute','on','savetrials','on'};
        disp('creating single trials ica data');
        limo_create_single_trials(subject_set,options{:})
    end
end
guidata(hObject, handles);
uiresume; delete(handles.figure1)

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
delete(handles.figure1)
