function varargout = limo_tools(varargin)

% GUI for the LIMO EEG Tools
% Created using GUIDE 
% Cyril Pernet 25-08-2009 v1
% update 21-01-2015 v2
% -----------------------------
%  Copyright (C) LIMO Team 2015

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
varargout{1} = 'LIMO tools off';


%% Callbacks

% ---------------------------------------------------------------

% --- Executes on button press in create_single_trials.
function create_single_trials_Callback(hObject, eventdata, handles)
limo_create_single_trials_gui
guidata(hObject, handles);

% --- Executes on button press in Split_continuous.
function Split_continuous_Callback(hObject, eventdata, handles)
limo_split_continuous
guidata(hObject, handles);

% --- Executes on button press in File_creation.
function File_creation_Callback(hObject, eventdata, handles)
limo_create_files
guidata(hObject, handles);

% --- Executes on button press in path_updates.
function path_updates_Callback(hObject, eventdata, handles)
limo_path_update
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function batch_Callback(hObject, eventdata, handles)
limo_batch;
guidata(hObject, handles);


% ---------------------------------------------------------------

% --- Executes on button press in electrode_optimization.
function electrode_optimization_Callback(hObject, eventdata, handles)
limo_best_electrodes
guidata(hObject, handles);

% --- Executes on button press in Expected_chanlocs.
function Expected_chanlocs_Callback(hObject, eventdata, handles)
[expected_chanlocs, channeighbstructmat] = limo_expected_chanlocs;
figure;  topoplot(zeros(1,71), expected_chanlocs,'style','blank','electrodes','labelpoint','chaninfo',expected_chanlocs);
figure; imagesc(channeighbstructmat); 
label = {}; colormap(jet); 
for i=1:2:length(expected_chanlocs); label{i}= expected_chanlocs(i).labels; end
set(gca,'YTick',[1:2:length(expected_chanlocs)],'YTickLabel', label)
label = {};
for i=2:2:length(expected_chanlocs); label{i}= expected_chanlocs(i).labels; end
set(gca,'XTick',[2:2:length(expected_chanlocs)],'XTickLabel', label)
title('Connectivity matrix between channels','FontSize',14)
positive = 1; while positive == 1
    [y,x]=ginput(1);
    if x<0 || y<0
        positive = 0;
    else
        
        if channeighbstructmat(round(x),round(y)) == 0
            channeighbstructmat(round(x),round(y)) = 1;
            imagesc(channeighbstructmat)
        else
            channeighbstructmat(round(x),round(y)) = 0;
            imagesc(channeighbstructmat)
        end
        label = {}; colormap(hot);
        for i=1:2:length(expected_chanlocs); label{i}= expected_chanlocs(i).labels; end
        set(gca,'YTick',[1:2:length(expected_chanlocs)],'YTickLabel', label)
        label = {};
        for i=2:2:length(expected_chanlocs); label{i}= expected_chanlocs(i).labels; end
        set(gca,'XTick',[2:2:length(expected_chanlocs)],'XTickLabel', label)
        title('Connectivity matrix between channels','FontSize',14)
    end
end
D=uigetdir(pwd,'Save file in directory');
if D == 0
    disp('data not saved'); return
else
    cd(D); save expected_chanlocs expected_chanlocs channeighbstructmat % save all in one file
    fprintf('expected_chanlocs & channeighbstructmatfile saved\n');
end
guidata(hObject, handles);

% ---------------------------------------------------------------

% --- Executes on button press in Help.
function Help_Callback(hObject, eventdata, handles)

origin = which('limo_eeg'); origin = origin(1:end-10); 
origin = sprintf('%shelp',origin); cd(origin)
web(['file://' which('limo_tools.html')]);

% --- Executes on button press in Quit.
function Quit_Callback(hObject, eventdata, handles)

clc
uiresume
guidata(hObject, handles);
delete(handles.figure1)
limo_gui




