function varargout = limo_itc_gui(varargin)

% BATCH ITC INTERFACE
% created using GUIDE 
% Based on limo_batch_gui
%
% Andrew Stewart, May 2014
% -----------------------------
% Copyright (C) LIMO Team 2014


%% GUI stuffs
% -------------------------
% Begin initialization code
% -------------------------
warning off

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @limo_itc_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @limo_itc_gui_OutputFcn, ...
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
function limo_itc_gui_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% define handles used for the save callback
handles.FileName = [];
handles.CatName             = [];
handles.fullfactorial       = 0;
handles.ContName            = [];
handles.zscore              = 1;
handles.Cat                 = [];
handles.Cont                = [];
handles.start               = [];
handles.end                 = [];
handles.lowf                = [];
handles.highf               = [];
handles.Analysis            = [];
handles.type_of_analysis    = 'Mass-univariate';
handles.method              = 'OLS';
handles.bootstrap           = 0;
handles.tfce                = 0;
handles.quit                = 0;

handles.Analysis = 'Time-Frequency';

guidata(hObject, handles);
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = limo_itc_gui_OutputFcn(hObject, eventdata, handles) 

varargout{1} = 'ITC out';

if handles.quit == 1
    varargout{1} = [] ;
    varargout{2} = [] ;
    varargout{3} = [] ;
    varargout{4} = [] ;
else
    varargout{1} = handles.FileName;
    varargout{2} = handles.Cat;
    varargout{3} = handles.Cont;
    varargout{4} = handles.defaults ;
end
delete(handles.figure1)


%% Callbacks

%-------------------------
%         IMPORT
%------------------------

% load a data set -- EEG 
% ---------------------------------------------------------------
function Import_data_set_Callback(hObject, eventdata, handles)

[FileName,PathName,FilterIndex]=uigetfile({'*.txt;*.mat','List of dataset files (*.txt or *.mat)'; ...
        '*set*', 'EEGLAB (EEG.set)'}, ...
        'Pick sets or list', 'MultiSelect', 'on');
if FilterIndex ~= 0
    
    if iscell(FileName) % multiselect for .sets
        for f=1:size(FileName,2)
            if ~strcmp(FileName{f}(end-3:end),'.set')
                errordlg('only multiple set files are allowed or a single .mat/.txt list')
                return
            end
        end
        
    else % use .mat or .txt
        
        if strcmp(FileName(end-3:end),'.txt')
            %FileName = importdata(FileName);
            FileName = importdata([PathName FileName]);  %axs
        elseif strcmp(FileName(end-3:end),'.mat')
            FileName = load([PathName FileName]);
            FileName = getfield(FileName,cell2mat(fieldnames(FileName)));
        end
    
        for f=1:size(FileName,1)
            if ~exist(FileName{f},'file')
                errordlg(sprintf('%s \n file not found',FileName{f}));
                return
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

% For ITC, set as 'Time-Frequency'

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






%-------------------------
%         OTHERS
%------------------------

% --- Executes on button press in Done.
% ---------------------------------------------------------------
function Done_Callback(hObject, eventdata, handles)
  
defaults.fullfactorial     = handles.fullfactorial;
defaults.zscore            = handles.zscore;
defaults.start             = handles.start;
defaults.end               = handles.end ;
defaults.lowf              = handles.lowf;
defaults.highf             = handles.highf;
defaults.method            = handles.method;
defaults.type_of_analysis  = handles.type_of_analysis;  
defaults.analysis          = 'ITC';                  % Set to ITC
defaults.bootstrap         = handles.bootstrap;  
defaults.tfce              = handles.tfce;

defaults.chanlocs          = handles.chanlocs;

handles.defaults = defaults;


    uiresume
    guidata(hObject, handles);




% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)


clc
uiresume
handles.quit = 1;
guidata(hObject, handles);
% delete(handles.figure1)
limo_gui


% --- Executes on button press in Select_chanlocs.
function Select_chanlocs_Callback(hObject, eventdata, handles)
% hObject    handle to Select_chanlocs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[FileName,PathName,FilterIndex]=uigetfile({'*.mat'}, ...
    'Select chanloc mat', 'MultiSelect', 'on');

FileName = importdata([PathName FileName]);  %axs

handles.chanlocs = FileName;

guidata(hObject, handles);


% --- Executes on selection change in popupmenu9.
function popupmenu9_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu9 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu9


% --- Executes during object creation, after setting all properties.
function popupmenu9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bootstrapchb.
function bootstrapchb_Callback(hObject, eventdata, handles)
% hObject    handle to bootstrapchb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bootstrapchb


% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox13


% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10


% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11


% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8


% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on bootstrapchb and none of its controls.
function bootstrapchb_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to bootstrapchb (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

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


% --- Executes on selection change in type_of_analysis.
function popupmenu11_Callback(hObject, eventdata, handles)
% hObject    handle to type_of_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns type_of_analysis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from type_of_analysis


% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to type_of_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in boostrap_check_box.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to boostrap_check_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of boostrap_check_box


% --- Executes on button press in TFCE.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to TFCE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TFCE


% --- Executes on selection change in method.
function popupmenu12_Callback(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from method


% --- Executes during object creation, after setting all properties.
function popupmenu12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function boostrap_check_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to boostrap_check_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in Help.
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dir = pwd;
origin = which('limo_eeg'); origin = origin(1:end-10); 
origin = sprintf('%shelp',origin); cd(origin)
web(['file://' which('limo_itc.html')]);
cd (handles.dir)


%-------------------------
%         SPECIFY
%------------------------

% --- Executes on button press in categorical_variable_input.
% ---------------------------------------------------------------
function categorical_variable_input_Callback(hObject, eventdata, handles)

[FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','LIMO categorical data');
if FilterIndex == 1 
    cd(PathName); 
    if strcmp(FileName(end-3:end),'.txt')
        handles.Cat = load(FileName);
    else
        load(FileName)
        handles.Cat = eval(FileName(1:end-4));
    end
    
    % if there is more than one factor, allow factorial design
    if size(handles.Cat,2) > 1
        set(handles.full_factorial,'Enable','on')
        handles.fullfactorial = 0;
    else
        handles.fullfactorial = 0;
    end
    disp('Categorical data loaded');
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

[FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','LIMO continuous data');
if FilterIndex == 1
    cd(PathName); 
    if strcmp(FileName(end-3:end),'.txt')
        handles.Cont = load(FileName);
    else
        load(FileName)
        handles.Cont = eval(FileName(1:end-4));
    end
    
    % if the regressors are not zscored, allow option to leave it as such 
    % test mean = 0 with a margin of 10^-5
    M = mean(mean(handles.Cont));
    centered = M>-0.00001 && M<0.00001;
    % if mean = 0 also test std = 1
    if centered == 1
        S = mean(std(handles.Cont));
        reducted = S>0.99999 && S<1.00001;
    else
        reducted = 0;        
    end
    
    if  centered~=1 && reducted~= 1
        set(handles.full_factorial,'Enable','on')
        handles.zscore = 1;
    else
        handles.zscore = 1;
    end
    disp('Continuous data loaded');
end
guidata(hObject, handles);


% --- Executes on button press in full_factorial.
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
