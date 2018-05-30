function varargout = limo_batch_gui(varargin)

% BATCH INTERFACE
% created using GUIDE 
% Based on limo_import_tf
% Cyril Pernet v1. May 2014
% -----------------------------
% Copyright (C) LIMO Team 2015


%% GUI stuffs
% -------------------------
% Begin initialization code
% -------------------------
warning off

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @limo_batch_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @limo_batch_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:4}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% -----------------------
% End initialization code
% -----------------------


% --------------------------------------------------
%   Executes just before the menu is made visible
% --------------------------------------------------
function limo_batch_gui_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% define handles used for the save callback
handles.FileName = [];
handles.CatName             = [];
handles.fullfactorial       = 0;
handles.ContName            = [];
handles.zscore              = 1;
handles.start               = [];
handles.end                 = [];
handles.lowf                = [];
handles.highf               = [];
handles.Analysis            = [];
handles.type                = 'Channels';
handles.type_of_analysis    = 'Mass-univariate';
handles.method              = 'WLS';
handles.bootstrap           = 0;
handles.tfce                = 0;
handles.quit                = 0;

guidata(hObject, handles);
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = limo_batch_gui_OutputFcn(hObject, eventdata, handles) 

if handles.quit == 1
    varargout{1} = [] ;
    varargout{2} = [] ;
    varargout{3} = [] ;
    varargout{4} = [] ;
else
    varargout{1} = handles.FileName;
    varargout{2} = handles.CatName;
    varargout{3} = handles.ContName;
    varargout{4} = handles.defaults ;
end
delete(handles.figure1); return


%% Callbacks

%-------------------------
%         IMPORT_DATA
%------------------------

% load a data set -- EEG 
% ---------------------------------------------------------------
function Import_data_Callback(hObject, eventdata, handles)
[FileName,PathName,FilterIndex]=uigetfile({'*.txt; *.mat; *.set'}, 'Pick sets or list', 'MultiSelect', 'on');

if FilterIndex ~= 0
    if iscell(FileName) % multiselect for .sets
        for f=1:size(FileName,2)
            if ~strcmp(FileName{f}(end-3:end),'.set')
                errordlg('only set files are allowed or use a single .mat/.txt list')
                return
            end
        end
        
    elseif ischar(FileName) && strcmp(FileName(end-3:end),'.set') % single subject
        tmp = FileName; clear FileName;
        FileName{1} = [PathName tmp]; clear tmp; % a single subject
        
    elseif strcmp(FileName(end-3:end),'.txt') ||  strcmp(FileName(end-3:end),'.mat') % use .mat or .txt
        
        if strcmp(FileName(end-3:end),'.txt')
            FileName = importdata(FileName);
        elseif strcmp(FileName(end-3:end),'.mat')
            FileName = load([PathName FileName]);
            FileName = getfield(FileName,cell2mat(fieldnames(FileName)));
        end
    
        disp('checking files .. ')
        for f=1:size(FileName,1)
            if ~exist(FileName{f},'file')
                errordlg(sprintf('%s \n file not found',FileName{f}));
                return
            else
                fprintf('%s found \n',FileName{f}); 
            end
        end
    end
    handles.FileName = FileName;
    
end
guidata(hObject, handles);


% --- Executes on selection change in analysis_type.
% ---------------------------------------------------------------
function analysis_type_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function analysis_type_Callback(hObject, eventdata, handles)

content = get(hObject,'Value');
if content == 1
    handles.Analysis = 'Time';
    set(handles.Starting_point,'Enable','on')
    set(handles.ending_point,'Enable','on')
    set(handles.low_freq,'Enable','off')
    set(handles.high_freq,'Enable','off')
elseif content == 2
    handles.Analysis = 'Frequency';
    set(handles.Starting_point,'Enable','off')
    set(handles.ending_point,'Enable','off')
    set(handles.low_freq,'Enable','on')
    set(handles.high_freq,'Enable','on')
elseif content == 3
    handles.Analysis = 'Time-Frequency';
    set(handles.Starting_point,'Enable','on')
    set(handles.ending_point,'Enable','on')
    set(handles.low_freq,'Enable','on')
    set(handles.high_freq,'Enable','on')
end
guidata(hObject, handles);


% get the starting time point of the analysis
% ---------------------------------------------------------------
function Starting_point_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Starting_point_Callback(hObject, eventdata, handles)

handles.start = str2double(get(hObject,'String'));
guidata(hObject, handles);


% get the ending time point of the analysis
% ---------------------------------------------------------------
function ending_point_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ending_point_Callback(hObject, eventdata, handles)

handles.end = str2double(get(hObject,'String'));
guidata(hObject, handles);


% get the starting frequency point of the analysis
% ---------------------------------------------------------------
function low_freq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function low_freq_Callback(hObject, eventdata, handles)

handles.lowf = str2double(get(hObject,'String'));
guidata(hObject, handles);


% get the ending frequency point of the analysis
% ---------------------------------------------------------------
function high_freq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function high_freq_Callback(hObject, eventdata, handles)

handles.highf = str2double(get(hObject,'String'));
guidata(hObject, handles);

%---------------------------
%      ANALYSIS
% --------------------------
% --- Executes on button press in scalp_data.
function scalp_data_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
if h == 0
    handles.type = '[]';
    set(handles.component_data,'Enable','on')
elseif h == 1
    handles.type = 'Channels';
    set(handles.component_data,'Enable','off')
end
guidata(hObject, handles);

% --- Executes on button press in component_data.
function component_data_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
if h == 0
    handles.type = '[]';
    set(handles.scalp_data,'Enable','on')
elseif h == 1
    handles.type = 'Components';
    set(handles.scalp_data,'Enable','off')
end
guidata(hObject, handles);


% type of analysis
% --- Executes on selection change in type_of_analysis.
% --- Executes during object creation, after setting all properties.
function type_of_analysis_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function type_of_analysis_Callback(hObject, eventdata, handles)

contents{1} = 'Mass-univariate';  
contents{2} = 'Multivariate'; 
handles.type_of_analysis = contents{get(hObject,'Value')};
if isempty(handles.type_of_analysis)
    handles.type_of_analysis = 'Mass-univariate';
end
fprintf('analysis selected %s \n',handles.type_of_analysis);
guidata(hObject, handles);


% method
function method_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function method_Callback(hObject, eventdata, handles)

contents{1} = 'WLS'; contents{2} = 'IRLS'; contents{3} = 'OLS';
handles.method = contents{get(hObject,'Value')};
if isempty(handles.method)
    handles.method = 'WLS';
end
fprintf('method selected %s \n',handles.method);
guidata(hObject, handles);


% bootstrap
% --- Executes on button press in boostrap_check_box.
function boostrap_check_box_Callback(hObject, eventdata, handles)
M = get(hObject,'Value');
if M == 1
    handles.bootstrap = 1;
    disp('bootstrap is on');
    set(handles.TFCE,'Enable','on')
elseif M == 0
    handles.bootstrap = 0;
    disp('boostrap is off');
    set(handles.TFCE,'Enable','off')
end
guidata(hObject, handles);


% TFCE
% --- Executes on button press in TFCE.
function TFCE_CreateFcn(hObject, eventdata, handles)
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



%-------------------------
%         SPECIFY
%------------------------

% --- Executes on button press in categorical_variable_input.
% ---------------------------------------------------------------
function categorical_variable_input_Callback(hObject, eventdata, handles)

[CatName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','LIMO categorical data','Multiselect','on');
if FilterIndex == 1
    if ischar(CatName) % NOT multiselect
        if strcmp(CatName(end-3:end),'.txt')
            if ~isnumeric(importdata([PathName CatName]))
                CatName = importdata([PathName CatName]);
            else
                tmp = CatName; clear CatName;
                CatName{1} = [PathName tmp]; clear tmp; % a single subject
            end
        elseif strcmp(CatName(end-3:end),'.mat')
            CatName = load([PathName CatName]);
            CatName = getfield(CatName,cell2mat(fieldnames(CatName)));
        end
        
        for f=1:size(CatName,1)
            if ~exist(CatName{f},'file')
                errordlg(sprintf('%s \n file not found',CatName{f}));
                return
            else
                fprintf('%s found \n',CatName{f});
            end
        end
    else
        errordlg('file selection not supported'); return
    end
    handles.CatName = CatName;
end
guidata(hObject, handles);


% --- Executes on button press in full_factorial.
% ---------------------------------------------------------------
function full_factorial_Callback(hObject, eventdata, handles)
M = get(hObject,'Value');
if M == 1
    handles.fullfactorial = 1;
    disp('full factorial on');
elseif M == 0
    handles.fullfactorial = 0;
    disp('full factorial off');
end
guidata(hObject, handles);


% --- Executes on button press in continuous_variable_input.
% ---------------------------------------------------------------
function continuous_variable_input_Callback(hObject, eventdata, handles)

[ContName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','LIMO continuous data','Multiselect','on');
if FilterIndex == 1
    if ischar(ContName) % NOT multiselect
        if strcmp(ContName(end-3:end),'.txt')
            if ~isnumeric(importdata([PathName ContName]))
                ContName = importdata([PathName ContName]);
            else
                tmp = ContName; clear ContName;
                ContName{1} = [PathName tmp]; clear tmp; % a single subject
            end
        elseif strcmp(ContName(end-3:end),'.mat')
            ContName = load([PathName ContName]);
            ContName = getfield(ContName,cell2mat(fieldnames(ContName)));
        end
        
        for f=1:size(ContName,1)
            if ~exist(ContName{f},'file')
                errordlg(sprintf('%s \n file not found',ContName{f}));
                return
            else
                fprintf('%s found \n',ContName{f});
            end
        end
    else
        errordlg('file selection not supported'); return
    end
    handles.ContName = ContName;
end
guidata(hObject, handles);


% --- Executes on button press in z_score.
% ---------------------------------------------------------------
function z_score_Callback(hObject, eventdata, handles)
M = get(hObject,'Value');
if M == 0
    handles.zscore = 1;
    disp('zscoring on');
elseif M == 1
    handles.zscore = 0;
    disp('zscoring off');
end
guidata(hObject, handles);


%-------------------------
%         OTHERS
%------------------------

% --- Executes on button press in Done.
% ---------------------------------------------------------------
function Done_Callback(hObject, eventdata, handles)
  
if isempty(handles.Analysis)
        errordlg('choose a type of analysis to perfom','error')
    return
end

defaults.analysis          = handles.Analysis;  
defaults.fullfactorial     = handles.fullfactorial;
defaults.zscore            = handles.zscore;
defaults.start             = handles.start;
defaults.end               = handles.end ;
defaults.lowf              = handles.lowf;
defaults.highf             = handles.highf;
defaults.method            = handles.method;
defaults.type_of_analysis  = handles.type_of_analysis;  
defaults.bootstrap         = handles.bootstrap;  
defaults.tfce              = handles.tfce;  
defaults.type              = handles.type; % model.defaults.type


% -----------------------------------------
% load the expected channel locations
% -----------------------------------------
if handles.bootstrap == 1 && ~strcmp(handles.type,'Components') 
    [chan_file,chan_path,whatsup]=uigetfile('expected_chanlocs.mat','Select channel location file');
    if whatsup == 1
        load (sprintf('%s%s',chan_path,chan_file))
        test = eval(chan_file(1:end-4));
        if isstruct(test) && ~isempty(test(1).labels) && ~isempty(test(1).theta) && ~isempty(test(1).radius) ...
                && ~isempty(test(1).X) && ~isempty(test(1).Y) && ~isempty(test(1).Z) && ~isempty(test(1).sph_theta) ...
                && ~isempty(test(1).sph_phi) && ~isempty(test(1).sph_radius) 
            defaults.chanloc = test; disp('channel location loaded');
        else
            warndlg('this file is not recognize as a channel location file or informations are missing','file error')
        end
    else
        disp('exiting batch mode'); limo_gui; return
    end
end
handles.defaults = defaults;

if isempty(handles.CatName) && isempty(handles.ContName);
    errordlg('no regressors were loaded','error')
    return
else
    uiresume
    guidata(hObject, handles);
end

% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)
persistent out

clc
uiresume
handles.quit = 1;
guidata(hObject, handles);
limo_gui
