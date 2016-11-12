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
% In New_Contrast_Callback --> test new contrast
% In Done_Callback --> update LIMO.contrast and run new contrast
%
% see also limo_contrast_checking.m limo_contrast.m
%
% Nicolas Chauveau & Cyril Pernet 29-04-2009 v2
% updated for 1st level bootstrap + some fix 20-06-2013
% ----------------------------------------------------
%  Copyright (C) LIMO Team 2010

%% GUI stuffs
% -------------------------
% Begin initialization code
% -------------------------
warning off
global limofile

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @limo_contrast_manager_OpeningFcn, ...
                   'gui_OutputFcn',  @limo_contrast_manager_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    if exist(varargin{1},'file') == 2 && ~isempty(strfind(varargin{1},'LIMO.mat'))
        limofile = varargin{1};
    else
        gui_State.gui_Callback = str2func(varargin{1});
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


% --------------------------------------------------
%   Executes just before the menu is made visible
% --------------------------------------------------
function limo_contrast_manager_OpeningFcn(hObject, eventdata, handles, varargin)
global handles

% define handles used for the save callback
handles.dir = pwd;
handles.go  = 0;
handles.C   = [];
handles.F   = 0;
handles.X   = [];
handles.output = hObject;
guidata(hObject,handles);
set(hObject,'Tag','figure_limo_contrast_manager');
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = limo_contrast_manager_OutputFcn(hObject, eventdata, handles) 
varargout{1} = 'contrast done';


%% Callbacks

% --- Display the design matrix
% ---------------------------------------------------------------
function display_matrix_CreateFcn(hObject, eventdata, handles)
global LIMO handles

if isempty(handles.limofile)
    [FileName,PathName,FilterIndex]=uigetfile('LIMO.mat','Select a LIMO file');
else
    [PathName,FileName] = fileparts(handles.limofile);
    if isempty(PathName), PathName = pwd; end
    FilterIndex = 1;
end
if FilterIndex ==0
    clear all; 
    varargout{1} = 'contrast cancelled';
else
    cd(PathName); load('LIMO.mat');
    if LIMO.Level == 2 && strcmp(LIMO.design.name,'one sample t-test all electrodes') || ...
            LIMO.Level == 2 && strcmp(LIMO.design.name,'one sample t-test one electrodes') || ...
            LIMO.Level == 2 && strcmp(LIMO.design.name,'paired t-test all electrodes') || ...
            LIMO.Level == 2 && strcmp(LIMO.design.name,'paired t-test one electrodes') || ...
            LIMO.Level == 2 && strcmp(LIMO.design.name,'two sample t-test all electrodes') || ...
            LIMO.Level == 2 && strcmp(LIMO.design.name,'two sample t-test one electrodes')
        warndlg('no contrast can be performed for a t-test'); guidata(hObject, handles);
    end
    
    if LIMO.Level == 2 && strcmp(LIMO.design.name(1:11),'Mixed ANOVA') || ...
            LIMO.Level == 2 && strcmp(LIMO.design.name(1:8),'Repeated') % always F
        if LIMO.design.nb_conditions == 1 % only one gp
            Xdisplay = LIMO.design.X;
        else
            Xdisplay = LIMO.design.X(1:prod(LIMO.design.repeated_measure),1:prod(LIMO.design.repeated_measure));
        end
    else
        if sum(LIMO.design.nb_conditions) ~=0 && LIMO.design.nb_continuous == 0
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
    
    imagesc(Xdisplay); colormap('gray');
    title('design matrix'); drawnow
    handles.output = hObject;
    guidata(hObject,handles)
end

% --- Display selected contrasts
% ---------------------------------------------------------------
function contrast_CreateFcn(hObject, eventdata, handles)
global handles

try
    imagesc(handles.C);
catch
    imagesc([]);
end
handles.output = hObject;
guidata(hObject,handles)


% --- Evaluate New Contrast
% ---------------------------------------------------------------
function New_Contrast_Callback(hObject, eventdata, handles)
global LIMO handles

handles.C = str2num(get(hObject,'String'));
  
if LIMO.Level == 2 && strncmp(LIMO.design.name,'Repeated',8) % always T2 in fact (contrast T but F stat)
    if handles.F ~= 0
        uiwait(warndlg('Only T contrasts are used for this design','Hotelling T^2 test','modal'));
    end
    % adjust X for a single group
    if LIMO.design.nb_conditions > 1
        index = find(LIMO.data.Cat==1);
        handles.C = limo_contrast_checking(LIMO.dir, LIMO.design.X(index,[1 2]), handles.C);
        handles.go = limo_contrast_checking(handles.C,LIMO.design.X(index,[1 2]));
    else
        handles.C = limo_contrast_checking(LIMO.dir, LIMO.design.X, handles.C);
        handles.go = limo_contrast_checking(handles.C,LIMO.design.X);
    end
else % other things than repeated measure
    if handles.F == 0
        handles.C = limo_contrast_checking(LIMO.dir, LIMO.design.X, handles.C);
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
        handles.C = C1; handles.C = limo_contrast_checking(LIMO.dir, LIMO.design.X, handles.C);
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
    
    
function New_Contrast_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in F_contrast.
% ---------------------------------------------------------------
function F_contrast_Callback(hObject, eventdata, handles)
global LIMO handles

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
 

% --- Previous Contrast
% ---------------------------------------------------------------
function Pop_up_previous_contrasts_CreateFcn(hObject, eventdata, handles)
global LIMO handles limofile

if exist('limofile','var') == 1 && ~isempty(limofile)
    handles.limofile = limofile;
else
    handles.limofile = [];
end
clear global limofile

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

previous_con = 0;
try
    previous_con = size(LIMO.contrast,2);
catch
    previous_con = 0;
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
    handles.C = LIMO.contrast{1}.C;
    % contrast_CreateFcn(hObject, eventdata, handles)
else
    contrasts = {'none'};
    set(hObject,'String',contrasts);
    display_matrix_CreateFcn(hObject, eventdata, handles)
    % display_matrix_CreateFcn(hObject, eventdata, handles)
end

handles.output = hObject;
guidata(hObject,handles)


  
function Pop_up_previous_contrasts_Callback(hObject, eventdata, handles)
global handles

contents = get(hObject,'String');
if ~strcmp(contents(1),'none') % not 'none'
    index = get(hObject,'Value');
    contents = contents{index}(4:end); %$ remove the T: or F: then eval string
    contrast=str2num(contents);
    handles.C=contrast;
    contrast_CreateFcn(hObject, eventdata, handles)
end

handles.output = hObject;
guidata(hObject,handles)



% --- Executes on button press in Done.
% ---------------------------------------------------------------
function Done_Callback(hObject, eventdata, handles)
global LIMO handles

if ~isempty(handles.C);
    if handles.go == 1
        
        if LIMO.design.bootstrap ==1
            choice = questdlg('an associated bootstrap should be present','bootstrap choice','compute bootstrap contrast','don''t compute any bootstraps','compute bootstrap contrast');
        else
            choice = 'don''t compute any bootstraps';
        end
        
        % ------------------------------------------------
        % 1st level contrat & 2st level ANOVA/ANCOVA/Regression
        % ------------------------------------------------
        if LIMO.Level == 1 || ...
                LIMO.Level == 2 && strncmp(LIMO.design.name,'regression',10) || ...
                LIMO.Level == 2 && strncmp(LIMO.design.name,'N-ways',6) || ...
                LIMO.Level == 2 && strncmp(LIMO.design.name,'ANCOVA',6)
            
            try
                previous_con = size(LIMO.contrast,2);
            catch ME
                previous_con = 0;
            end
            index = previous_con+1;
            
            % update LIMO.mat
            LIMO.contrast{index}.C = handles.C;
            if handles.F == 0
                LIMO.contrast{index}.V = 'T';
            else
                LIMO.contrast{index}.V = 'F';
            end
            save LIMO LIMO
            
            load Yr; load Betas;
            % -------------------------------------------------------
            if strcmp(LIMO.design.type_of_analysis,'Mass-univariate')
                % -------------------------------------------------------
                result = limo_contrast(Yr, Betas, LIMO, handles.F,1);
                
                if LIMO.design.bootstrap ~= 0
                    if strcmp(LIMO.Analysis ,'Time-Frequency')
                        disp('preparing data matrix');
                        clear Yr Betas % release some memory
                        cd H0; load H0_Betas.mat;
                        tmp = zeros(size(H0_Betas,1), size(H0_Betas,2)*size(H0_Betas,3), size(H0_Betas,4), size(H0_Betas,5));
                        for boot = 1:size(H0_Betas,5)
                           tmp(:,:,:,boot)= limo_tf_4d_reshape(squeeze(H0_Betas(:,:,:,:,boot)));
                        end
                        clear H0_Betas; cd ..; load Yr; cd(H0);
                        result = limo_contrast(limo_tf_4d_reshape(Yr), tmp, LIMO, handles.F,2);
                        clear tmp
                    else
                        clear Betas; cd H0; load H0_Betas
                        result = limo_contrast(Yr, H0_Betas, LIMO, handles.F,2);
                    end
                    clear Yr ; cd ..
                end
                
                if LIMO.design.tfce == 1
                    if handles.F == 0
                        filename = sprintf('con_%g.mat',index); load(filename);
                        if strcmp(LIMO.Analysis ,'Time-Frequency')
                            tfce_score = limo_tfce(3,squeeze(con(:,:,4)),LIMO.data.neighbouring_matrix);
                        else
                            tfce_score = limo_tfce(2,squeeze(con(:,:,4)),LIMO.data.neighbouring_matrix);
                        end
                        cd TFCE; filename2 = sprintf('tfce_%s',filename); save ([filename2], 'tfce_score'); clear con tfce_score
                        cd ..; cd H0; filename = sprintf('H0_%s',filename); load(filename);
                        if strcmp(LIMO.Analysis ,'Time-Frequency')
                            tfce_H0_score = limo_tfce(3,squeeze(H0_con(:,:,2,:)),LIMO.data.neighbouring_matrix);
                        else
                            tfce_H0_score = limo_tfce(2,squeeze(H0_con(:,:,2,:)),LIMO.data.neighbouring_matrix);
                        end
                        filename2 = sprintf('tfce_%s',filename); save ([filename2], 'tfce_H0_score'); clear H0_con tfce_score
                    else
                        filename = sprintf('ess_%g.mat',index); load(filename);
                        if strcmp(LIMO.Analysis ,'Time-Frequency')
                            tfce_score = limo_tfce(3,squeeze(ess(:,:,end-1)),LIMO.data.neighbouring_matrix);
                        else
                            tfce_score = limo_tfce(2,squeeze(ess(:,:,end-1)),LIMO.data.neighbouring_matrix);
                        end
                        cd TFCE; filename2 = sprintf('tfce_%s',filename); save ([filename2], 'tfce_score'); clear ess tfce_score
                        cd ..; cd H0; filename = sprintf('H0_%s',filename); load(filename);
                        if strcmp(LIMO.Analysis ,'Time-Frequency')
                            tfce_H0_score = limo_tfce(3,squeeze(H0_ess(:,:,end-1,:)),LIMO.data.neighbouring_matrix);
                        else
                            tfce_H0_score = limo_tfce(2,squeeze(H0_ess(:,:,end-1,:)),LIMO.data.neighbouring_matrix);
                        end
                        filename2 = sprintf('tfce_%s',filename); save ([filename2], 'tfce_H0_score'); clear H0_ess tfce_score
                    end
                end
                
            elseif strcmp(LIMO.design.type_of_analysis,'Multivariate')
                
                LIMO.contrast = handles.F;
                save LIMO LIMO
                result = limo_contrast(squeeze(Yr(:,time,:))', squeeze(Betas(:,time,:))', [], LIMO, handles.F,1);
                
            end
            
            clear Yr Betas
            disp('boostrapped contrasts done ...')
            
            
            % -------------------------------------------
            %          2nd level Repeated measure ANOVA
            % -------------------------------------------
        elseif LIMO.Level == 2 && strncmp(LIMO.design.name,'Repeated',8)
            
            try
                previous_con = size(LIMO.contrast,2);
            catch ME
                previous_con = 0;
            end
            index = previous_con+1;
            
            % update LIMO.mat
            LIMO.contrast{index}.C = handles.C;
            LIMO.contrast{index}.V = 'F';
            C = handles.C;
            
            % create ess files and call limo_rep_anova adding C
            load Yr; save LIMO LIMO
            result = limo_contrast(Yr,LIMO,3);
            
            if LIMO.design.bootstrap ~= 0
                cd H0; limo_contrast(Yr, LIMO, 4,); cd ..
            end
            
            if LIMO.design.tfce == 1
                filename = sprintf('ess_repeated_measure_%g.mat',index); load(filename);
                tfce_score = limo_tfce(squeeze(ess(:,:,1)),LIMO.data.neighbouring_matrix);
                cd TFCE; filename2 = sprintf('tfce_%s',filename); save ([filename2], 'tfce_score'); clear ess tfce_score
                cd ..; cd H0; filename = sprintf('H0_%s',filename); load(filename);
                tfce_H0_score = limo_tfce(squeeze(H0_ess(:,:,1,:)),LIMO.data.neighbouring_matrix);
                filename2 = sprintf('tfce_%s',filename); save ([filename2], 'tfce_H0_score'); clear H0_ess tfce_score
            end
            
            clear Yr LIMO
            disp('contrast evaluation done ...')
               
                
        end 
        
    else
        errordlg('no new contrast to evaluate')
    end
    
    uiresume
    guidata(hObject, handles);
    close(get(hObject,'Parent'));
    if isempty(handles.limofile)
        limo_results;
    end
end


% --- Executes on button press in Quit.
% ---------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)
clc; 
uiresume
guidata(hObject, handles);
close(get(hObject,'Parent')); 
if isempty(handles.limofile)
    limo_results;
end

clearvars LIMO handles limofile
return
