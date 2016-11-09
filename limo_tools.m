function varargout = limo_tools(varargin)

% GUI for the LIMO EEG Tools
% Created using GUIDE
% Cyril Pernet
% -----------------------------
%  Copyright (C) LIMO Team 2016

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

choice =  questdlg2('Do you you to Create or Edit a groupo level file?', ...
    'Choice', ...
    'Create', 'Edit', 'Edit');
if strcmp(choice,'Create')
    [expected_chanlocs, channeighbstructmat] = limo_expected_chanlocs;
else
    [gp_level_file,filepath,sts]=uigetfile('*.mat','select gp level channel file');
    if sts ==0
        return
    else
        load([filepath filesep gp_level_file])
    end
end

% show channels
figure
topoplot(zeros(1,71), expected_chanlocs,'style','blank','electrodes','labelpoint','chaninfo',expected_chanlocs);

% show connectivity matrix
figure
imagesc(channeighbstructmat); colormap(gray);
for i=1:length(expected_chanlocs);
    try
        label{i}= expected_chanlocs(i).urchan;
    catch
        label{i}= expected_chanlocs(i).labels;
    end
end
set(gca,'YTick',[1:3:length(expected_chanlocs)],'YTickLabel', label(1:3:length(expected_chanlocs)))
set(gca,'XTick',[2:3:length(expected_chanlocs)],'XTickLabel', label(2:3:length(expected_chanlocs)))
axis([1 length(expected_chanlocs) 1 length(expected_chanlocs)]); axis square
title(sprintf('Connectivity matrix between channels \n'),'FontSize',14)

% interactive editing
positive = 1;
while positive == 1
    [y,x]=ginput(1);
    if x<0 || y<0
        positive = 0;
    else
        if channeighbstructmat(round(x),round(y)) == 0
            channeighbstructmat(round(x),round(y)) = 1;
            channeighbstructmat(round(y),round(x)) = 1;
            imagesc(channeighbstructmat); v = 'on';
        else
            channeighbstructmat(round(x),round(y)) = 0;
            channeighbstructmat(round(y),round(x)) = 0;
            imagesc(channeighbstructmat);  v = 'off';
        end
        colormap(gray);
        set(gca,'YTick',[1:3:length(expected_chanlocs)],'YTickLabel', label(1:3:length(expected_chanlocs)))
        set(gca,'XTick',[2:3:length(expected_chanlocs)],'XTickLabel', label(2:3:length(expected_chanlocs)))
        axis([1 length(expected_chanlocs) 1 length(expected_chanlocs)]); axis square
        title(sprintf('Connectivity matrix between channels \nconnection %g %g %s',round(x),round(y),v),'FontSize',14)
    end
end

% save
D=uigetdir(pwd,'Save file in directory');
if D == 0
    disp('data not saved'); return
else
    save([D filesep 'edited_' gp_level_file],'expected_chanlocs','channeighbstructmat') % save all in one file
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






