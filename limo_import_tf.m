function varargout = limo_import_tf(varargin)

% import function for getting timefreq data into limo
% created using GUIDE -- import the various
% information needed to process the data
% Based on Cyril's limo_import_t
% Andrew X Stewart, January 2014
% Fixed folder locations, time/freq components in .data
% updated data format Cyril v2. Feb 2014
% -----------------------------
%  Copyright (C) LIMO Team 2013


%% GUI stuffs
% -------------------------
% Begin initialization code
% -------------------------
warning off

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @limo_import_tf_OpeningFcn, ...
                   'gui_OutputFcn',  @limo_import_tf_OutputFcn, ...
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
function limo_import_tf_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% define handles used for the save callback
try
    clear LIMO
    LIMO     = [];
catch
    LIMO     = [];    
end

handles.data_dir            = [];
handles.data                = [];
handles.chanlocs            = [];
handles.type_of_analysis    = 'Mass-univariate';
handles.method              = 'WLS';
handles.rate                = [];
handles.trim1               = [];
handles.trim2               = [];
handles.trim_lowf           = [];
handles.trim_highf          = [];
handles.Cat                 = [];
handles.Cont                = [];
handles.bootstrap           = 0;
handles.start               = 0;
handles.end                 = 0;
handles.lowf                = 0;
handles.highf               = 0;
handles.tf_times            = 0;
handles.tf_freqs            = 0;
handles.dir                 = [];
handles.zscore              = 1;
handles.fullfactorial       = 0;
handles.dir                 = pwd;
handles.bootstrap           = 0;
handles.tfce                = 0;
handles.out                 = [];

guidata(hObject, handles);
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = limo_import_tf_OutputFcn(hObject, eventdata, handles) 
if isempty(handles.out)
    varargout{1} = 'LIMO import terminated';
else
    varargout{1} = handles.out;
end
delete(handles.figure1)

%% Callbacks

%-------------------------
%         IMPORT
%------------------------

% load a data set -- EEG 
% ---------------------------------------------------------------
function Import_data_set_Callback(hObject, eventdata, handles)
global EEGLIMO 

[FileName,PathName,FilterIndex]=uigetfile('*.set','EEGLAB EEG epoch data with TF data');
if FilterIndex ~= 0
    try
        disp('loading TF EEGLAB dataset. Please wait ...');
        EEGLIMO=pop_loadset([PathName FileName]);
        handles.data_dir = PathName;
        handles.data     = FileName;
        handles.chanlocs = EEGLIMO.chanlocs;
        handles.rate     = EEGLIMO.srate;
        
        if isfield(EEGLIMO.etc,'timeersp') == 1 && isfield(EEGLIMO.etc,'freqersp') == 1
            handles.start      = EEGLIMO.etc.timeersp(1);
            handles.end        = EEGLIMO.etc.timeersp(end);
            handles.lowf       = EEGLIMO.etc.freqersp(1);
            handles.highf      = EEGLIMO.etc.freqersp(end);
            handles.tf_times   = EEGLIMO.etc.timeersp;
            handles.tf_freqs   = EEGLIMO.etc.freqersp;
            cd(handles.dir)
            fprintf('Data set %s loaded \n',FileName);
        else
            errordlg('Can''t load the data. Ensure that time-frequency information is stored in EEG.etc - see help.'); return
        end
        
    catch
        errordlg('pop_loadset eeglab function error / not found','error'); return
    end
end
guidata(hObject, handles);


% get the starting time point of the analysis
% ---------------------------------------------------------------
function Starting_point_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Starting_point_Callback(hObject, eventdata, handles)
global EEGLIMO 

start = str2double(get(hObject,'String'));
if start < EEGLIMO.etc.timeersp(1)
    errordlg(['The earliest time possible is:',num2str(EEGLIMO.etc.timeersp(1)),'ms']);
else
    % Find a possible frequency bin close to the requested one
    [a1 ind] = min(abs(EEGLIMO.etc.timeersp-start));
    closest_start = EEGLIMO.etc.timeersp(ind);
    if start ~= closest_start
        warndlg2(['this will be adjusted to sampling rate, start at:',num2str(closest_start),'ms']);
    end
    handles.start   = closest_start;
    handles.trim1   = ind; % gives the 1st column to start the analysis
end
guidata(hObject, handles);


% get the ending time point of the analysis
% ---------------------------------------------------------------
function ending_point_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ending_point_Callback(hObject, eventdata, handles)
global EEGLIMO 

ending = str2double(get(hObject,'String'));
if ending > EEGLIMO.etc.timeersp(end)
    errordlg(['The latest time possible is:',num2str(EEGLIMO.etc.timeersp(end)),'ms']);
else
    % Find a possible frequency bin close to the requested one
    if isnan(ending)
        [a1 ind] = max(EEGLIMO.etc.timeersp);
    else
        [a1 ind] = min(abs(EEGLIMO.etc.timeersp-ending));
    end
    closest_ending = EEGLIMO.etc.timeersp(ind);
    if ending ~= closest_ending
        warndlg2(['this will be adjusted to sampling rate, end at:',num2str(closest_ending),'ms']);
    end
    handles.end   = closest_ending;
    handles.trim2 = ind; % gives the 1st column to start the analysis
end

guidata(hObject, handles);


% get the starting frequency point of the analysis
% ---------------------------------------------------------------
function low_freq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function low_freq_Callback(hObject, eventdata, handles)
global EEGLIMO 

lowf = str2double(get(hObject,'String'));
if lowf < EEGLIMO.etc.freqersp(1)
    errordlg(['The lowest frequency possible is:',num2str(EEGLIMO.etc.freqersp(1))]);
else
    % Find a possible frequency bin close to the requested one
    [a1 ind] = min(abs(EEGLIMO.etc.freqersp-lowf));
    closest_lowf = EEGLIMO.etc.freqersp(ind);
    if lowf ~= closest_lowf
        warndlg2(['this will be adjusted to the closest frequency bin:',num2str(closest_lowf),'Hz']);
    end
    handles.lowf    = closest_lowf;
    handles.trim_lowf    = ind; % gives the 1st column to start the analysis
end

guidata(hObject, handles);


% get the ending frequency point of the analysis
% ---------------------------------------------------------------
function high_freq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function high_freq_Callback(hObject, eventdata, handles)
global EEGLIMO 

highf = str2double(get(hObject,'String'));
if highf > EEGLIMO.etc.freqersp(end)
    errordlg(['The highest frequency possible is:',num2str(handles.tf_freqs(end))]);
else
    % Find a possible frequency bin close to the requested one
    if isnan(highf)
        [a1 ind] = max(EEGLIMO.etc.freqersp);
    else
        [a1 ind] = min(abs(EEGLIMO.etc.freqersp-highf));
    end
    closest_highf = EEGLIMO.etc.freqersp(ind);
    if highf ~= closest_highf
        warndlg2(['this will be adjusted to the closest frequency bin:',num2str(closest_highf),'Hz']);
    end
    handles.highf    = closest_highf;
    handles.trim_highf    = ind; % gives the 1st column to start the analysis
end

guidata(hObject, handles);

%---------------------------
%      ANALYSIS
% --------------------------

% type of analysis
% --- Executes on selection change in type_of_analysis.
function type_of_analysis_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4

contents{1} = 'Mass-univariate';  
contents{2} = 'Multivariate'; 
handles.type_of_analysis = contents{get(hObject,'Value')};
if isempty(handles.type_of_analysis)
    handles.type_of_analysis = 'Mass-univariate';
end
fprintf('analysis selected %s \n',handles.type_of_analysis);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function type_of_analysis_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% method
function method_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function method_Callback(hObject, eventdata, handles)

contents{1} = 'OLS'; contents{2} = 'WLS'; contents{3} = 'IRLS';
handles.method = contents{get(hObject,'Value')};
if isempty(handles.method)
    handles.method = 'OLS';
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

[FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','LIMO categorical data');
if FilterIndex == 1 
    cd(PathName); 
    if strcmp(FileName(end-3:end),'.txt')
        handles.Cat = load(FileName);
    else
        cat = load(FileName);
        handles.Cat = getfield(cat,cell2mat(fieldnames(cat)));
        clear cat
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
        cont = load(FileName);
        handles.Cont = getfield(cont,cell2mat(fieldnames(cont)));
        clear cont
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
        set(handles.z_score,'Enable','on')
        handles.zscore = 1;
    else
        handles.zscore = 1;
    end
    disp('Continuous data loaded');
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

% --- Executes on button press in Directory.
% ---------------------------------------------------------------
function Directory_Callback(hObject, eventdata, handles)

PathName=uigetdir(handles.dir,'select LIMO working directory');
if PathName ~= 0
    cd(PathName); 
    handles.dir = PathName;
end
guidata(hObject, handles);



% --- Executes on button press in Help.
% ---------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
global EEGLIMO LIMO 

origin = which('limo_eeg'); origin = origin(1:end-10); 
origin = sprintf('%shelp',origin); cd(origin)
web(['file://' which('limo_importtf.html')]);
cd (handles.dir)


 
% --- Executes on button press in Done.
% ---------------------------------------------------------------
function Done_Callback(hObject, eventdata, handles)
global EEGLIMO LIMO 
  
LIMO.data.data_dir            = handles.data_dir;
LIMO.data.data                = handles.data;
LIMO.data.chanlocs            = handles.chanlocs;
LIMO.data.sampling_rate       = handles.rate;
LIMO.data.Cat                 = handles.Cat;      
LIMO.data.Cont                = handles.Cont; 

LIMO.design.fullfactorial     = handles.fullfactorial;
LIMO.design.zscore            = handles.zscore;
LIMO.design.method            = 'WLS';
LIMO.design.type_of_analysis  = handles.type_of_analysis;  
LIMO.design.bootstrap         = handles.bootstrap;  
LIMO.design.tfce              = handles.tfce;  
LIMO.Level                    = 1;
LIMO.Analysis                 = 'Time-Frequency';

% set defaults - take from 1:numel if trim not set
if isempty(handles.trim_lowf)
    LIMO.data.trim_low_f = 1;
else
    LIMO.data.trim_low_f = handles.trim_lowf;
end

if isempty(handles.trim_highf)
    LIMO.data.trim_high_f = numel(handles.tf_freqs);
else
    LIMO.data.trim_high_f = handles.trim_highf;
end

LIMO.data.tf_freqs = handles.tf_freqs(LIMO.data.trim_low_f:LIMO.data.trim_high_f);
LIMO.data.lowf    = handles.lowf ;
LIMO.data.highf   = handles.highf ;

if isempty(handles.trim1)
    LIMO.data.trim1 = 1;
else
    LIMO.data.trim1 = handles.trim1;
end

if isempty(handles.trim2)
    LIMO.data.trim2 = numel(handles.tf_times);
else
    LIMO.data.trim2 = handles.trim2;
end

LIMO.data.tf_times = handles.tf_times(LIMO.data.trim1:LIMO.data.trim2);
LIMO.data.start    = handles.start;
LIMO.data.end      = handles.end ;

if isempty(handles.dir)
    LIMO.dir = handles.data_dir;
else
    LIMO.dir = handles.dir;
end

test = isempty(handles.Cat) + isempty(handles.Cont);
if test == 2
    errordlg('no regressors were loaded','error')
else
    cd (LIMO.dir);
    save LIMO LIMO
    uiresume
    guidata(hObject, handles);
end

% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)

clc
uiresume
handles.out = 'LIMO import aborded';
guidata(hObject, handles);
limo_gui
