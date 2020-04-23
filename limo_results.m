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

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
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
handles.p = 0.05;
handles.MCC = 1;
handles.dir = pwd;
handles.bootstrap = 0;
handles.tfce = 0;
handles.filter = {'*.mat'};

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

% show results with handles.dirmagesc
% ---------------------------------------------------------------
function Image_results_Callback(hObject, eventdata, handles)


[FileName,PathName,FilterIndex]=uigetfile(handles.filter,'Select Univariate Results to display','*.mat');
if FilterIndex ~= 0
    cd(PathName); handles.LIMO = load('LIMO.mat');
    
    % check if bootstrap or tfce should be computed
    % ---------------------------------------------
    % 1st level
    if handles.LIMO.LIMO.Level == 1
        if handles.bootstrap == 1 && ~exist(sprintf('H0%sH0_%s', filesep, FileName), 'file') ...
                && strncmp(FileName,'con',3) == 0 && strncmp(FileName,'ess',3) ==0
            if strcmp(questdlg('Level 1: compute all bootstraps?','bootstrap turned on','Yes','No','No'),'Yes')
                LIMO = handles.LIMO.LIMO;
                LIMO.design.bootstrap = 1;
                if handles.tfce == 1
                    LIMO.design.tfce = 1;
                end
                save LIMO LIMO
                limo_eeg(4);
            end
        end
        
        if handles.tfce == 1 && ~exist(sprintf('TFCE%stfce_%s', filesep, FileName), 'file') ...
                && exist(sprintf('H0%sH0_%s', filesep, FileName), 'file') && strncmp(FileName,'con',3) == 0 ...
                && strncmp(FileName,'ess',3) ==0
            if strcmp(questdlg('Level 1: compute all tfce?','tfce turned on','Yes','No','No'),'Yes')
                LIMO = handles.LIMO.LIMO;
                LIMO.design.tfce = 1;
                save LIMO LIMO
                limo_eeg(4);
            end
        end
    end
    
    % contrasts stuff
    if handles.bootstrap == 1 && ~exist(sprintf('H0%sH0_%s', filesep, FileName), 'file')
        Yr       = load('Yr.mat'); Yr = Yr.Yr; 
        H0_Betas = load(['H0' filesep 'H0_Betas.mat']); H0_Betas = H0_Betas.H0_Betas;
        if strncmp(FileName,'con',3)
            limo_contrast(Yr, H0_Betas, handles.LIMO.LIMO, 0,3); clear Yr H0_Betas
        elseif strncmp(FileName,'ess',3)
            limo_contrast(Yr, H0_Betas, handles.LIMO.LIMO, 1,3); clear Yr H0_Betas
        end
    end
    
    if handles.tfce == 1 && ~exist(sprintf('TFCE%stfce_%s', filesep, FileName), 'file') ...
            && exist(sprintf('H0%sH0_%s', filesep, FileName), 'file')
        if strncmp(FileName,'con',3)
            load(FileName);
            if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency')
                x = 3;
            else
                x = size(con,1);
                if x~=1
                    x=2;
                end
            end
            tfce_score = limo_tfce(x,squeeze(con(:,:,2)),handles.LIMO.LIMO.data.neighbouring_matrix);
            cd TFCE;
            filename2 = sprintf('tfce_%s',FileName);
            save ([filename2], 'tfce_score');
            clear con tfce_score
            cd ..; cd H0;
            filename = sprintf('H0_%s',FileName);
            load(filename);
            tfce_H0_score = limo_tfce(x,squeeze(H0_ess(:,:,2,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
            filename2 = sprintf('tfce_%s',filename);
            save ([filename2], 'tfce_H0_score');
            clear H0_con tfce_score
        elseif strncmp(FileName,'ess',3)
            load(FileName);
            if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency')
                x = 3;
            else
                x = size(ess,1);
                if x~=1
                    x=2;
                end
            end
            tfce_score = limo_tfce(x,squeeze(ess(:,:,2)),handles.LIMO.LIMO.data.neighbouring_matrix);
            cd TFCE;
            filename2 = sprintf('tfce_%s',FileName);
            save ([filename2], 'tfce_score');
            clear ess tfce_score
            cd ..; cd H0;
            filename = sprintf('H0_%s',FileName); load(filename);
            tfce_H0_score = limo_tfce(x,squeeze(H0_ess(:,:,2,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
            filename2 = sprintf('tfce_%s',filename);
            save ([filename2], 'tfce_H0_score');
            clear H0_ess tfce_score
        end
    end
    
    % 2nd level
    % ------------
    nboot = 1000;
    if handles.LIMO.LIMO.Level == 2
        if handles.bootstrap == 1 && ~exist(sprintf('H0%sH0_%s', filesep, FileName), 'file')
            if strncmp(FileName,'one_sample',10)
                load Yr;
                limo_random_robust(1,Yr,eval(FileName(28:end-4)),nboot,handles.tfce); clear Yr;
            elseif strncmp(FileName,'two_samples',11)
                load Y1r; load Y2r;
                limo_random_robust(2,Y1r,Y2r,eval(FileName(29:end-4)),nboot,handles.tfce); clear Y1r Y2r;
            elseif strncmp(FileName,'paired_samples',14)
                load Y1r; load Y2r;
                limo_random_robust(3,Y1r,Y2r,eval(FileName(32:end-4)),nboot,handles.tfce); clear Y1r Y2r;
            elseif strncmp(FileName,'Rep_ANOVA',9)
                LIMO = handles.LIMO.LIMO;
                LIMO.design.bootstrap = 1000;
                LIMO.design.tfce = handles.tfce;
                save LIMO LIMO
                limo_random_robust(6,[],LIMO.data.Cat, LIMO.design.repeated_measure, ...
                    LIMO,LIMO.design.bootstrap, LIMO.design.tfce)
            else
                limo_eeg(4);
            end
            LIMO = handles.LIMO.LIMO; LIMO.design.bootstrap = 1; save LIMO LIMO
            handles.LIMO.LIMO.design.bootstrap = 1;
        end
        
        if handles.tfce == 1 && ~exist(sprintf('TFCE%stfce_%s', filesep, FileName), 'file') ...
                && exist(sprintf('H0%sH0_%s', filesep, FileName), 'file')
            mkdir tfce; load(FileName); load(sprintf('H0%sH0_%s', filesep, FileName));
            if strncmp(FileName,'one_sample',10)
                parameter = eval(FileName(28:end-4));
                tfce_name = sprintf('tfce_one_sample_ttest_parameter_%g',parameter);
                tfce_H0_name = sprintf('tfce_H0_one_sample_ttest_parameter_%g',parameter);
                if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency');
                    x = size(one_sample,1);
                    if x==1
                        x=2; LIMO.LIMO.data.neighbouring_matrix = [];
                    else
                        x=3;
                    end
                    tfce_one_sample = limo_tfce(x,squeeze(one_sample(:,:,:,4)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    save(['tfce', filesep, tfce_name], 'tfce_one_sample'); clear tfce_one_sample;
                    tfce_H0_one_sample = limo_tfce(x,squeeze(H0_one_sample(:,:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    save(['H0', filesep, tfce_H0_name],'tfce_H0_one_sample'); clear tfce_H0_one_sample;
                else
                    x = size(one_sample,1); if x~=1; x=2; end
                    tfce_one_sample = limo_tfce(x,squeeze(one_sample(:,:,4)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    save(['tfce', filesep, tfce_name], 'tfce_one_sample'); clear tfce_one_sample;
                    tfce_H0_one_sample = limo_tfce(x,squeeze(H0_one_sample(:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    save(['H0', filesep, tfce_H0_name],'tfce_H0_one_sample'); clear tfce_H0_one_sample;
                end
            elseif strncmp(FileName,'two_samples',11)
                parameter = eval(FileName(29:end-4));
                tfce_name = sprintf('tfce_two_samples_ttest_parameter_%g',parameter);
                tfce_H0_name = sprintf('tfce_H0_two_samples_ttest_parameter_%g',parameter);
                if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency');
                    x = size(one_sample,1);
                    if x==1
                        x=2; LIMO.LIMO.data.neighbouring_matrix = [];
                    else
                        x=3;
                    end
                    tfce_two_samples = limo_tfce(x,squeeze(two_samples(:,:,:,4)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    save(['tfce', filesep, tfce_name], 'tfce_two_samples'); clear tfce_two_samples;
                    tfce_H0_two_samples = limo_tfce(x,squeeze(H0_two_samples(:,:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    save(['H0', filesep, tfce_H0_name],'tfce_H0_two_samples'); clear tfce_H0_two_samples;
                else
                    x = size(two_samples,1); if x~=1; x=2; end
                    tfce_two_samples = limo_tfce(x,squeeze(two_samples(:,:,4)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    save(['tfce', filesep, tfce_name], 'tfce_two_samples'); clear tfce_two_samples;
                    tfce_H0_two_samples = limo_tfce(x,squeeze(H0_two_samples(:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    save(['H0', filesep, tfce_H0_name],'tfce_H0_two_samples'); clear tfce_H0_two_samples;
                end
            elseif strncmp(FileName,'paired_samples',14)
                parameter = eval(FileName(32:end-4));
                tfce_name = sprintf('tfce_paired_samples_ttest_parameter_%g',parameter);
                tfce_H0_name = sprintf('tfce_H0_paired_samples_ttest_parameter_%g',parameter);
                if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency')
                    x = size(one_sample,1);
                    if x==1
                        x=2; LIMO.LIMO.data.neighbouring_matrix = [];
                    else
                        x=3;
                    end
                    tfce_paired_samples = limo_tfce(x,squeeze(paired_samples(:,:,:,4)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    save(['tfce', filesep, tfce_name], 'tfce_paired_samples'); clear tfce_paired_samples;
                    tfce_H0_paired_samples = limo_tfce(x,squeeze(H0_paired_samples(:,:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    save(['H0', filesep, tfce_H0_name],'tfce_H0_paired_samples'); clear tfce_H0_paired_samples;
                else
                    x = size(paired_samples,1); if x~=1; x=2; end
                    tfce_paired_samples = limo_tfce(x,squeeze(paired_samples(:,:,4)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    save(['tfce', filesep, tfce_name], 'tfce_paired_samples'); clear tfce_paired_samples;
                    tfce_H0_paired_samples = limo_tfce(x,squeeze(H0_paired_samples(:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    save(['H0', filesep, tfce_H0_name],'tfce_H0_paired_samples'); clear tfce_H0_paired_samples;
                end
            elseif strncmp(FileName,'Covariate_effect',16)
                if size(Covariate_effect,1) == 1
                    tfce_score(1,:) = limo_tfce(1,squeeze(Covariate_effect(:,:,1)),handles.LIMO.LIMO.data.neighbouring_matrix);
                else
                    tfce_score = limo_tfce(2,squeeze(Covariate_effect(:,:,1)),handles.LIMO.LIMO.data.neighbouring_matrix);
                end
                tfce_name = sprintf('tfce_%s',FileName); save(['tfce', filesep, tfce_name],'tfce_score');
                clear Covariate_effect tfce_score;
                
                cd('H0'); fprintf('Creating H0 Covariate TFCE scores \n');
                name = sprintf('H0_Covariate_effect_%s.mat',FileName(18:end-4));
                H0_Covariate_effect = load(name); 
                H0_Covariate_effect = H0_Covariate_effect.(cell2mat(fieldnames(H0_Covariate_effect))); 
                PCT_test = ver('distcomp');
                if size(H0_Covariate_effect,1) == 1
                    if ~isempty(PCT_test)
                        tfce_H0_score = NaN(1,size(H0_Covariate_effect,2),handles.LIMO.LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Covariate_effect(:,:,:,1,b)),handles.LIMO.LIMO.data.neighbouring_matrix,0);
                        end
                    else
                        tfce_H0_score(1,:,:) = limo_tfce(1,squeeze(H0_Covariate_effect(:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    end
                else
                    if ~isempty(PCT_test)
                        tfce_H0_score = NaN(size(H0_Covariate_effect,1),size(H0_Covariate_effect,2),handles.LIMO.LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Covariate_effect(:,:,1,b)),handles.LIMO.LIMO.data.neighbouring_matrix,0);
                        end
                    else
                        tfce_H0_score = limo_tfce(2,squeeze(H0_Covariate_effect(:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                    end
                end
                full_name = sprintf('tfce_%s',name); save(full_name,'tfce_H0_score'); cd ..
                clear H0_Covariate_effect tfce_H0_score;
                
            elseif strncmp(FileName,'Rep_ANOVA',9)
                LIMO = handles.LIMO.LIMO;
                LIMO.design.tfce = handles.tfce;
                save LIMO LIMO
                limo_random_robust(6,[],LIMO.data.Cat, LIMO.design.repeated_measure, ...
                    LIMO,LIMO.design.bootstrap, LIMO.design.tfce)
            else
                limo_eeg(4);
            end
        end
        LIMO = handles.LIMO.LIMO; LIMO.design.tfce = 1; save LIMO LIMO
        handles.LIMO.LIMO.design.tfce = 1;
    end
    
    % do the figure
    % -------------
    limo_display_results(1,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO);
end
guidata(hObject, handles);

% Topoplot
% ---------------------------------------------------------------
function Topoplot_Callback(hObject, eventdata, handles)

[FileName,PathName,FilterIndex]=uigetfile(handles.filter,'Select Result to plot');
if FilterIndex == 1
    cd(PathName); handles.LIMO = load('LIMO.mat');
    limo_display_results(2,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO);
end
guidata(hObject, handles);

% course plots
% ---------------------------------------------------------------
function ERP_Callback(hObject, eventdata, handles)

supported_tests = {'one_sample','two samples','paired','Rep_ANOVA'};
[FileName,PathName,FilterIndex]=uigetfile('*.mat','Select LIMO file or summary stat file');
if FilterIndex == 1  
    if strcmp(FileName,'LIMO.mat') && exist(fullfile(PathName,FileName),'file')
        cd(PathName); handles.LIMO = load('LIMO.mat'); files = dir('*.mat'); FileName = [];
        for i=1:length(files)
            if contains(files(i).name,supported_tests)
                FileName = files(i).name;
            end
        end
        
        if isempty(FileName)
            errordlg2('Associated LIMO file not found/supported'); return
        else
            limo_display_results(3,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO);
        end
        
    elseif contains(FileName,supported_tests)
        if exist(fullfile(PathName,'LIMO.mat'),'file')
             handles.LIMO = load(fullfile(PathName,'LIMO.mat'));
             limo_display_results(3,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO);
        else
             errordlg2('No associated LIMO file was found'); return
        end
        
    else
        test = load(fullfile(PathName,FileName));
        if isfield(test,'data')
            try
                limo_add_plots(fullfile(PathName,FileName))
            catch no_plot
                error('the file fits with variable ''data'' yet there was an error while plotting it: \n%s',no_plot.message)
            end
        else
            errordlg2('this file doesn''t seem supported'); return
        end
    end
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
    handles.bootstrap = 1;
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

% [file,dir_path] = uigetfile('*.mat','select a LIMO.mat file');
% if file ==0
%     return
%     clc; uiresume
%     guidata(hObject, handles);
% else
%     cd(dir_path); handles.LIMO = load('LIMO.mat');
%     clc; uiresume
%     guidata(hObject, handles);
%     delete(handles.figure1)
%     limo_contrast_manager(handles.LIMO.LIMO);
% end

% Semi-Partial coef
% ---------------------------------------------------------------
function Partial_coef_Callback(hObject, eventdata, handles)
[FileName,PathName,FilterIndex]=uigetfile('LIMO.mat','Select a LIMO file');
if FileName ~=0
    cd(PathName); handles.LIMO = load('LIMO.mat');
    limo_semi_partial_coef(handles.LIMO.LIMO)
end
guidata(hObject, handles);

% Model Selection (multivariate).
% ---------------------------------------------------------------
function Model_Selection_Callback(hObject, eventdata, handles)
[FileName,PathName,FilterIndex]=uigetfile('LIMO.mat','Select a LIMO file');
cd(PathName); handles.LIMO = load('LIMO.mat');
limo_model_selection(handles.LIMO.LIMO,1);
guidata(hObject, handles);

%------------------------
%         OTHERS
%------------------------

% --- Executes on button press in Help.
% ---------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)

origin = which('limo_eeg'); origin = origin(1:end-10);
origin = sprintf('%shelp',origin); cd(origin)
web(['file://' which('limo_results.html')]);
cd (handles.dir)

% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)

clc; uiresume
guidata(hObject, handles);
delete(handles.figure1)
limo_gui
