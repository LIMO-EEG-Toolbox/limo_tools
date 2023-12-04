function varargout = limo_contrast_manager(varargin)

% GUI to create / review contrasts
% created using GUIDE
%
% Usage:
% 1- limo_display_results   (Call the GUI. The file LIMO.mat file must be selected afterwards)
% 2- limo_display_results(PathtoFile) Will call the main GUI and load the LIMO.mat
%    file. PathtoFile must be the full or the relative path to the LIMO.mat
%    file. File name must be provided in the path too.
%    i.e PathtoFile = '/Users/username/WORK/LIMO.mat'
%
% In display_matrix_CreateFcn --> load LIMO.mat and display X
% In New_Contrast_Callback --> test new Contrastmatrix
% In Done_Callback --> update LIMO.Contrastmatrix and run new Contrastmatrix
%
% see also limo_contrast_checking.m limo_contrast.m
%
% Nicolas Chauveau, Arnaud Delorme & Cyril Pernet 29-04-2009 v2
% updated for 1st level bootstrap + some fix 20-06-2013
% ------------------------------
%  Copyright (C) LIMO Team 2019

%% GUI stuffs
% -------------------------
% Begin initialization code
% -------------------------
warning on
global limofile
global result

limo_settings_script;
if limo_settings.newgui
    guiName = [mfilename '_new'];
else
    guiName = mfilename;
end

gui_Singleton = 1;
gui_State = struct('gui_Name', guiName, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @limo_contrast_manager_OpeningFcn, ...
    'gui_OutputFcn',  @limo_contrast_manager_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin == 0
    limofile = [];
    result = [];
end
if nargin && ischar(varargin{1})
    if exist(varargin{1},'file') == 2 && ~isempty(strfind(varargin{1},'LIMO.mat'))
        limofile = varargin{1};
        result = [];
    else 
        gui_State.gui_Callback = str2func(varargin{1}); % only 1 parameter
    end
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% -----------------------
% End initialization code
% -----------------------


%   Executes just before the menu is made visible
% --------------------------------------------------
function limo_contrast_manager_OpeningFcn(hObject, eventdata, handles, varargin)
global limofile

variables = '';
if ~isempty(limofile)
    limofile = load('-mat', limofile);
    try
        variables = { limofile.LIMO.design.labels.description };
    catch
        variables = '';
    end
    if length(variables) == 1
        uiresume
        guidata(hObject, handles);
        close(hObject);
        clearvars LIMO limofile
        limo_errordlg('LIMO design has one variable or less so a contrast cannot be defined')
    end
end

% define handles used for the save callback
handles.dir      = pwd;
handles.go       = 0;
handles.C        = [];
handles.F        = 0;
handles.X        = [];
handles.limofile = limofile;
handles.variables = variables;
handles.output    = hObject;
handles.Name      = '';
guidata(hObject,handles);
listfactors1 = findobj(hObject, 'tag', 'Factorlist1');
listfactors2 = findobj(hObject, 'tag', 'Factorlist2');
if ~isempty(listfactors1) && ~isempty(listfactors2)
    if ~isempty(variables)
        set(listfactors1,'string', variables, 'value', 1, 'max', 2);
        set(listfactors2,'string', variables, 'value', 2, 'max', 2);
        handles.C = zeros(1, length(variables));    
        handles.C(1:2) = [1 -1];
        contrast_CreateFcn(hObject, eventdata, handles)
        handles.go = 1;
    else
        set(listfactors1,'string', {'LIMO file contains no variable' '(might be an old file, please recompute)'}, 'value', [], 'max', 2, 'enable', 'off');
        set(listfactors2,'string', {'LIMO file contains no variable' '(might be an old file, please recompute)'}, 'value', [], 'max', 2, 'enable', 'off');
    end
end
guidata(hObject, handles);

set(hObject,'Tag','figure_limo_contrast_manager');
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
% ---------------------------------------------------------------
function varargout = limo_contrast_manager_OutputFcn(hObject, eventdata, handles)

waitfor(findobj(hObject, 'string', 'Done'), 'userdata', 'done');
varargout{1} = handles;

%% Callbacks
% --- Display the design matrix
% ---------------------------------------------------------------
function update_contrast_CreateFcn(hObject, eventdata, handles)

selection1 = get(findobj(hObject.Parent,'tag','Facorlist1'), 'value');
selection2 = get(findobj(hObject.Parent,'tag','Facorlist2'), 'value');

% --- Display the design matrix
% ---------------------------------------------------------------
function display_matrix_CreateFcn(hObject, eventdata, handles)
global LIMO 
global limofile

if ~isfield(handles,'limofile') || isempty(handles.limofile) && ~isempty(limofile) 
    handles.limofile = limofile;
end

if ~isfield(handles,'limofile') || isempty(handles.limofile)
    [FileName,PathName,FilterIndex]=uigetfile('LIMO.mat','Select a LIMO file');
    if FilterIndex == 0
        guidata(hObject, handles);
        limo_results
    end
else
    [PathName,FileName] = fileparts(handles.limofile);
    if isempty(PathName), PathName = pwd; end
    FilterIndex = 1;
end

if FilterIndex ==0
    clear variables;
    varargout{1} = 'contrast cancelled';
else
    LIMO = load(fullfile(PathName,'LIMO.mat')); LIMO = LIMO.LIMO;
    if LIMO.Level == 2 && ~isempty(strfind(LIMO.design.name,'t-test'))
        warndlg('no contrast can be performed for a t-test','modal');
        Xdisplay = [];
    elseif LIMO.Level == 2 && contains(LIMO.design.name,'ANOVA')  % always F
        if LIMO.design.nb_conditions == 1 || ~contains(LIMO.design.name,'repeated')
            Xdisplay = LIMO.design.X;
        else
            Xdisplay = LIMO.design.X(1:prod(LIMO.design.repeated_measure),1:prod(LIMO.design.repeated_measure));
        end
    else
        if ~isfield(LIMO.design,'nb_conditions') && ~isfield(LIMO.design,'nb_continuous')
            Xdisplay = [];
        elseif sum(LIMO.design.nb_conditions) ~=0 && LIMO.design.nb_continuous == 0
            Xdisplay = LIMO.design.X;
        elseif sum(LIMO.design.nb_conditions) == 0
            Xdisplay = [zscore(LIMO.design.X(:,1:end-1)) LIMO.design.X(:,end)];
        else
            N = sum(LIMO.design.nb_conditions)+sum(LIMO.design.nb_interactions)+1;
            REGdisplay = LIMO.design.X(:,N:end-1); % covariates
            REGdisplay = REGdisplay + max(abs(min(REGdisplay)));
            Xdisplay = LIMO.design.X;
            Xdisplay(:,N:end-1) = REGdisplay ./ max(max(REGdisplay));
        end
    end
    
    if isempty(Xdisplay)
        uiresume; guidata(hObject, handles);
    else
        allhandles = get(get(get(hObject,'Parent'),'Parent'),'Children');
        %allhandles = findobj(hObject.Parent.Parent,'tag', 'LIMOmatrix');
        axes(allhandles(end));
        imagesc(Xdisplay); colormap('gray');
        title('design matrix'); drawnow
        handles.output = hObject;
        guidata(hObject,handles)
    end
end

% --- Display Contrastmatrix matrix
% ---------------------------------------------------------------
function contrast_CreateFcn(hObject, eventdata, handles)

allhandles = findobj(hObject, 'Tag','Contrastmatrix');
if isempty(allhandles)
    allhandles = findobj(hObject.Parent, 'Tag','Contrastmatrix');
end
axes(allhandles);
try
    imagesc(handles.C);
catch
    imagesc([]);
end
set(gca, 'clim', [-1 1]);
handles.output = hObject;
str = num2str(handles.C);
str = strrep(str, '        ', ' ');
str = strrep(str, '   ', ' ');
str = strrep(str, ' ', ' ');
set(findobj(hObject.Parent, 'Tag','New_Contrast'), 'string', str);
guidata(hObject,handles)

function New_Contrast_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- GUI Name
% ---------------------------------------------------------------
function ContrastName_Callback(hObject, eventdata, handles)
handles.Name = get(hObject,'string');
guidata(hObject, handles);

% --- Evaluate Listbox 1
% ---------------------------------------------------------------
function Factorlist1_Callback(hObject, eventdata, handles)

selection1 = get(hObject,'value')
Factorlist2  = findobj(hObject.Parent,'tag','Factorlist2');
selection2 = get(Factorlist2, 'value');
selection2 = setdiff(selection2, selection1);
set(Factorlist2, 'value', selection2);

newConstrast = zeros(1,length(handles.variables));
if length(selection1) ~= length(selection2)
    newConstrast(selection1) = 1/length(selection1);
    newConstrast(selection2) = -1/length(selection2);
else
    newConstrast(selection1) = 1;
    newConstrast(selection2) = -1;
end
handles.C = newConstrast;%$ remove the T: or F: then eval string

cName = findobj(hObject.Parent,'tag','ContrastName');
if length(selection1) == 1 && length(selection2) == 1
    handles.Name = [ handles.variables{selection1} ' vs ' handles.variables{selection2} ];
    handles.Name = strrep(handles.Name, 'type - ', '');
    set(cName, 'string', handles.Name, 'userdata', 'auto');
else
    if strcmpi(char(get(cName, 'userdata')), 'auto')
        set(cName, 'string', '', 'userdata', '');
    end
end

guidata(hObject, handles);
contrast_CreateFcn(hObject, eventdata, handles)
guidata(hObject, handles);
    
% set(findobj(hObject.Parent,'tag','New_Contrast'), 'string', num2str([weights{:}]));
% update_contrast_CreateFcn(hObject, eventdata, handles)

% --- Evaluate Listbox 1
% ---------------------------------------------------------------
function Factorlist2_Callback(hObject, eventdata, handles)

selection2 = get(hObject,'value');

Factorlist1 = findobj(hObject.Parent,'tag','Factorlist1');
selection1 = get(Factorlist1,'value');
selection1 = setdiff(selection1, selection2);
set(Factorlist1, 'value', selection1);

newConstrast = zeros(1,length(handles.variables));
if length(selection1) ~= length(selection2)
    newConstrast(selection1) = 1/length(selection1);
    newConstrast(selection2) = -1/length(selection2);
else
    newConstrast(selection1) = 1;
    newConstrast(selection2) = -1;
end
handles.C = newConstrast;%$ remove the T: or F: then eval string

cName = findobj(hObject.Parent,'tag','ContrastName');
if length(selection1) == 1 && length(selection2) == 1
    handles.Name = [ handles.variables{selection1} ' vs ' handles.variables{selection2} ];
    handles.Name = strrep(handles.Name, 'type - ', '');
    set(cName, 'string', handles.Name, 'userdata', 'auto');
else
    if strcmpi(char(get(cName, 'userdata')), 'auto')
        set(cName, 'string', '', 'userdata', '');
    end
end

guidata(hObject, handles);
contrast_CreateFcn(hObject, eventdata, handles)
guidata(hObject, handles);

% --- Evaluate New contrastmatrix
% ---------------------------------------------------------------
function New_Contrast_Callback(hObject, eventdata, handles)
global LIMO 

handles.C = str2num(get(hObject,'String'));
set(findobj(hObject.Parent,'tag','Factorlist1'), 'value', []);
set(findobj(hObject.Parent,'tag','Factorlist2'), 'value', []);

if LIMO.Level == 2 && contains(LIMO.design.name,'Repeated') % always T2 in fact (contrast T but F stat)
    if contains(LIMO.design.name,'') && handles.F ~= 0
        uiwait(warndlg('Only T contrasts are used for this design','Hotelling T^2 test','modal'));
    end
    handles.C  = limo_contrast_checking(LIMO.dir, LIMO.design.X, handles.C);
    handles.go = limo_contrast_checking(handles.C,LIMO.design.X);
else % other things than repeated measure
    if handles.F == 0
        handles.C  = limo_contrast_checking(LIMO.dir, LIMO.design.X, handles.C);
        if size(handles.C,1) > 1
            handles.F = 1;
        end
        handles.go = limo_contrast_checking(handles.C,LIMO.design.X);
    else
        L = zeros(1,size(LIMO.design.X,2));
        L(1:length(handles.C)) = handles.C;
        C1 = zeros(size(LIMO.design.X,2));
        for i=1:size(C1,1)
            C1(i,i) = L(i);
        end
        handles.C  = C1;
        handles.C  = limo_contrast_checking(LIMO.dir, LIMO.design.X, handles.C);
        handles.go = limo_contrast_checking(handles.C,LIMO.design.X);
    end
end

guidata(hObject, handles);
if handles.go == 0
    contrast_CreateFcn(hObject, eventdata, handles)
    warndlg('invalid contrast')
else
    contrast_CreateFcn(hObject, eventdata, handles)
    disp('press Done to evaluate')
end
guidata(hObject, handles);

allhandles = findobj('Tag','Contrastmatrix');
axes(allhandles);
try
    imagesc(handles.C);
catch
    imagesc([]);
end
handles.output = hObject;
guidata(hObject,handles)

% --- Executes on button press in F_contrast.
% ---------------------------------------------------------------
function F_contrast_Callback(hObject, eventdata, handles)

do_F_contrast = get(hObject,'Value');
if do_F_contrast == 0
    handles.F = 0;
elseif do_F_contrast == 1
    handles.F = 1;
end

handles.output = hObject;
guidata(hObject,handles)

if ~isempty(handles.C)
    contrast_CreateFcn(hObject, eventdata, handles)
end

% --- Previous contrastmatrix
% ---------------------------------------------------------------
function Pop_up_previous_contrasts_CreateFcn(hObject, eventdata, handles)
global LIMO 

if isempty(LIMO)
    display_matrix_CreateFcn(hObject, eventdata, handles)
    handles.limofile = [pwd filesep 'LIMO.mat'];
end

if isfield(LIMO,'contrast')
    handles.C = LIMO.contrast;
    previous_con = length(LIMO.contrast);
else
    previous_con = 0;
end

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

if previous_con ~=0
    for i=1:previous_con
        if LIMO.contrast{i}.V == 'T'
            contrasts{i} = ['T: ' num2str(LIMO.contrast{i}.C)];
        else
            s = [];
            for j=1:size(LIMO.contrast{i}.C,1)-1
                s = [s num2str(LIMO.contrast{i}.C(j,:)) ';'];
            end
            s = [s num2str(LIMO.contrast{i}.C(end,:))];
            contrasts{i} = ['F: ' s];
        end
    end
    set(hObject,'String',contrasts);
    display_matrix_CreateFcn(hObject, eventdata, handles)
    if isfield(LIMO,'contrast')
        handles.C = LIMO.contrast{1}.C;
    end
    % contrast_CreateFcn(hObject, eventdata, handles)
else
    contrasts = {'none'};
    set(hObject,'String',contrasts);
    display_matrix_CreateFcn(hObject, eventdata, handles)
    % display_matrix_CreateFcn(hObject, eventdata, handles)
end

allhandles = findobj('Tag','Contrastmatrix');
axes(allhandles);
imagesc([]);

handles.output = hObject;
guidata(hObject,handles)

function Pop_up_previous_contrasts_Callback(hObject, eventdata, handles)

contents = get(hObject,'String');
if ~strcmp(contents(1),'none') % not 'none'
    index = get(hObject,'Value');
    handles.C = str2num(contents{index}(4:end));%$ remove the T: or F: then eval string
    contrast_CreateFcn(hObject, eventdata, handles)
    % if no associated bootstrap set handles.go to 1 allowing to evaluate
    if strcmpi(contents{index}(1),'F')
        if ~exist(fullfile(handles.dir,['H0' filesep 'H0_ess_' num2str(index) '.mat']),'file')
            handles.go = 1;
        end
    elseif strcmpi(contents{index}(1),'T')
        if ~exist(fullfile(handles.dir,['H0' filesep 'H0_con_' num2str(index) '.mat']),'file')
            handles.go = 1;
        end
    end
end
handles.output = hObject;
guidata(hObject,handles)

% --- Executes on button press in Done.
% ---------------------------------------------------------------
function Done_Callback(hObject, eventdata, handles)
global LIMO 
global result 
global limofile 

set(findobj(gcbf, 'string', 'Done'), 'userdata', 'done');
result = handles.C;

if ~isempty(handles.C) && isempty(limofile)
    if handles.go == 1
        limo_contrast_execute(LIMO, handles);
    else
        warndlg2('no new contrast to evaluate')
    end
    
    uiresume
    guidata(hObject, handles);
    close(get(hObject,'Parent'));
    if isempty(handles.limofile)
        limo_results;
    end
else
    if ~isempty(limofile)
        guidata(hObject, handles);
        close(get(hObject,'Parent'));
    end
end

% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)
uiresume
guidata(hObject, handles);
close(get(hObject,'Parent'));
if isempty(handles.limofile)
    limo_results;
end

clearvars LIMO limofile
return
