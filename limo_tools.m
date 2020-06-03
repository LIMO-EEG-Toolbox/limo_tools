function varargout = limo_tools(varargin)

% GUI for the LIMO EEG Tools
% Created using GUIDE
% Cyril Pernet
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


% ---------------------------------------------------------------
% --- Executes on button press in check_weights.
function check_weights_Callback(hObject, eventdata, handles)
limo_CheckWeight
guidata(hObject, handles);


% --- Executes on button press in electrode_optimization.
function electrode_optimization_Callback(hObject, eventdata, handles)
limo_best_electrodes
guidata(hObject, handles);

% --- Executes on button press in Expected_chanlocs.
function Expected_chanlocs_Callback(hObject, eventdata, handles)

choice =  questdlg2('Do you want to Create or Edit a groupo level file?', ...
    'Choice', 'Create', 'Edit', 'Edit');

if strcmp(choice,'Create')

    [expected_chanlocs, channeighbstructmat] = limo_expected_chanlocs;
    figure('Name','Channel locations')
    topoplot([], expected_chanlocs,'style','blank','electrodes','labelpoint','chaninfo',expected_chanlocs);
    figure('Name','Neighbouring matrix')
    imagesc(channeighbstructmat); colormap(gray);
    for i=length(expected_chanlocs):-1:1
        if isfield(expected_chanlocs,'urchan')
            label{i}= expected_chanlocs(i).urchan;
        else
            label{i}= expected_chanlocs(i).labels;
        end
    end
    set(gca,'YTick',1:3:length(expected_chanlocs),'YTickLabel', label(1:3:length(expected_chanlocs)))
    set(gca,'XTick',2:3:length(expected_chanlocs),'XTickLabel', label(2:3:length(expected_chanlocs)))
    axis([1 length(expected_chanlocs) 1 length(expected_chanlocs)]); axis square
    title(sprintf('Connectivity matrix between channels \n (click outside the matrix or right click when done)'),'FontSize',14)
    cmap = gray; cmap(1,:) = [0.25 0.25 0.25]; colormap(cmap)
    
elseif strcmp(choice,'Edit')
    
    [expected_chanlocs, channeighbstructmat] = limo_edit_expected_chanlocs;
end

D = uigetdir(pwd,'Save file in directory');
if D == 0
    disp('data not saved'); return
else
    if strcmp(choice,'Create')
        save([D filesep 'gp_level_expected_channel'],'expected_chanlocs','channeighbstructmat') % save all in one file
    elseif strcmp(choice,'Edit')
        save([D filesep 'edited_gp_level_expected_channel'],'expected_chanlocs','channeighbstructmat') 
    end
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






