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
limo_settings_script;
if limo_settings.newgui
    guiName = [mfilename '_new'];
else
    guiName = mfilename;
end

gui_Singleton = 1;
gui_State = struct('gui_Name',guiName, ...
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

% set(handles.add_bootstrap,'Enable','off')
% set(handles.add_tfce,'Enable','off')

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
    handles      = check_boot_and_tfce(handles,fullfile(PathName,FileName));
    limo_display_results(1,FileName,PathName,handles.p,handles.MCC,handles.LIMO);
end

% reset selection - but annoying behaviour really
% uiresume
% guidata(hObject, handles);
% delete(handles.figure1)
% limo_results

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
    check_boot_and_tfce(handles,fullfile(PathName,FileName))
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
    check_boot_and_tfce(handles,fullfile(PathName,FileName));
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

% --- Executes during object creation, after setting all properties.
function add_bootstrap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to add_bootstrap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in add_bootstrap.
function add_bootstrap_Callback(hObject, eventdata, handles)
M = get(hObject,'Value');
if M == 1
    handles.bootstrap = 1000;
    disp('bootstrap is on');
elseif M == 0
    handles.bootstrap = 0;
    disp('boostrap is off');
end

guidata(hObject, handles);

% --- Executes on button press in add_tfce.
function add_tfce_Callback(hObject, eventdata, handles)
M = get(hObject,'Value');
if M == 1
    handles.tfce = 1;
    disp('tfce is on');
elseif M == 0
    handles.tfce = 0;
    disp('tfce is off');
end
guidata(hObject, handles);

% Model Selection (multivariate).
% ---------------------------------------------------------------
function Model_Selection_Callback(hObject, eventdata, handles)
[FileName,PathName,FilterIndex]=uigetfile('LIMO.mat','Select a LIMO file');
cd(PathName); handles.LIMO = load('LIMO.mat');
limo_model_selection(handles.LIMO.LIMO,1);
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

% -----------------------------------------------------------
% subroutine to check bootstrap and tfce on any button clicks
% ------------------------------------------------------------
function handles = check_boot_and_tfce(handles,currentfile)
if handles.bootstrap ~= 0 && handles.MCC ~= 1 || ...
        handles.tfce == 1 && handles.MCC == 3    
   
    handles.LIMO.design.bootstrap = handles.bootstrap;
    handles.LIMO.design.tfce      = handles.tfce;
    
    [filepath,filename,ext] = fileparts(currentfile);

    % deal with bootstrap
    if handles.bootstrap ~= 0
        if ~exist([filepath filesep 'H0' filesep 'H0_' filename ext],'file')
            if handles.LIMO.Level == 1
                if strncmp(filename,'con',3) || strncmp(filename,'ess',3)
                    if exist('warndlg2','file')
                        warndlg2(sprintf('This contrast cannot be bootstrapped now, \nbootstrap the model and recompute the contrast'))
                    else
                        warndlg(sprintf('This contrast cannot be bootstrapped now, \nbootstrap the model and recompute the contrast'))
                    end
                else
                    if strcmp(questdlg('Level 1: are you sure to compute all bootstraps for that subject?','bootstrap turned on','Yes','No','No'),'Yes')
                        LIMO                  = handles.LIMO;
                        LIMO.design.bootstrap = 800;
                        if handles.tfce == 1
                            LIMO.design.tfce  = 1;
                        end
                        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
                        limo_eeg(4);
                    end
                end
            else % handles.LIMO.Level == 2
                if contains(filename,'one_sample')
                        limo_random_robust(1,fullfile(handles.LIMO.dir,'Yr.mat'),...
                            str2num(filename(max(strfind(filename,'_'))+1:end)),handles.LIMO);
                elseif contains(filename,'two_samples')
                        limo_random_robust(2,fullfile(handles.LIMO.dir,'Y1r.mat'),...
                            fullfile(handles.LIMO.dir,'Y1r.mat'), str2num(filename(max(strfind(filename,'_'))+1:end)),handles.LIMO);
                elseif contains(filename,'paired_samples')
                        limo_random_robust(3,fullfile(handles.LIMO.dir,'Y1r.mat'),...
                            fullfile(handles.LIMO.dir,'Y1r.mat'), str2num(filename(max(strfind(filename,'_'))+1:end)),handles.LIMO);
                elseif contains(filename,'Covariate_effect') && contains(handles.LIMO.design.name,'Regression') 
                    LIMO = handles.LIMO; LIMO.design.bootstrap = 1000;
                    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO'); 
                    handles.LIMO = LIMO; limo_eeg(4,handles.LIMO.dir); clear LIMO
                elseif contains(filename,'ANOVA') && ~strncmpi(filename,'Rep_ANOVA',9)
                        limo_random_robust(5,fullfile(handles.LIMO.dir,'Yr.mat'),...
                            handles.LIMO.data.Cat,handles.LIMO.data.Cont,handles.LIMO,'go','yes');
                elseif contains(filename,'Rep_ANOVA')
                    if strncmp(filename,'con',3)
                        if exist([filepath filesep 'H0' filesep 'H0_' filesep 'H0_Betas.mat'],'file')
                            limo_contrast([filepath filesep 'Yr.mat'], ...
                                [filepath filesep 'H0' filesep 'H0_' filesep 'H0_Betas.mat'], handles.LIMO, 0,3);
                        else
                            errordlg2('there is no GLM bootstrap file for this contrast file')
                        end
                    elseif strncmp(filename,'ess',3)
                        if exist([filepath filesep 'H0' filesep 'H0_' filesep 'H0_Betas.mat'],'file')
                            limo_contrast([filepath filesep 'Yr.mat'], ...
                                [filepath filesep 'H0' filesep 'H0_' filesep 'H0_Betas.mat'], handles.LIMO, 1,3);
                        else
                            errordlg2('there is no bootstrap file for this contrast file')
                        end
                    else
                        disp('Bootstraping Repeated Measure ANOVA')
                        limo_random_robust(6,fullfile(filepath,'Yr.mat'),handles.LIMO.data.Cat, ...
                            handles.LIMO.design.repeated_measure, handles.LIMO, 'go','yes')
                    end
                end
            end
        end
    end
    
    % deal with tfce
    if handles.tfce == 1
        if ~exist([filepath filesep 'tfce' filesep 'tfce_' filename ext],'file')
            limo_tfce_handling(currentfile,'checkfile','yes')
        end
    end
end            

