function LIMO_files = limo_batch(varargin)

% interactive function to run several 1st level analyses
% select directories and files - possibly enter contrasts of
% interests and let it run. The batch relies on PSOM (see Ref)
% see opt.mode for parallel computing on grid using qsub or msub
% <https://code.google.com/p/psom/wiki/ConfigurationPsom>
%
% FORMAT limo_batch
% limo_batch(option,model,contrast)
% limo_batch(option,model,contrast,eeglab_study)
%
% INPUT if empty uses GUI
%       option should be 'model specification' 'contrast only' or 'both'
%       model is a structure that specifiy information to build a model
%       model.set_files: a cell array of EEG.set (full path) for the different subjects
%       model.cat_files: a cell array of categorial variable or variable files
%       model.cont_files: a cell array of continuous variable or variable files
%       model.defaults: specifiy the parameters to use for each subject
%       model.defaults.type = 'Channels' or 'Components'
%       model.defaults.analysis 'Time' 'Frequency' or 'Time-Frequency'
%       model.defaults.method 'WLS' 'IRLS' 'OLS'
%       model.defaults.fullfactorial 0/1
%       model.defaults.zscore 0/1
%       model.defaults.start starting time in ms
%       model.defaults.end ending time in ms
%       model.defaults.lowf starting point in Hz
%       model.defaults.highf ending point in Hz
%       model.defaults.bootstrap 0/1
%       model.defaults.tfce 0/1
%       model.defaults.channloc common channel locations (necessary if bootstrap = 1)
%       contrast is a structure that specify which contrasts to run for which subject
%       contrast.LIMO_files: a list of LIMO.mat (full path) for the different subjects
%                            this is optional if option 'both' is selected
%       contrast.mat: a matrix of contrasts to run (assumes the same for all subjects)
%       eeglab_study is the STUDY structure allowing to create multiple design with consistant names etc ... 
%
% OUTPUT  LIMO a cell array of LIMO.mat (info about subjects' GLM)
%         also generate a directory per subject with GLM results in it
%
% see also limo_eeg limo_import_t limo_import_f limo_import_tf and psom in external folder
%
% Reference for pipeline engine
% Bellec P, Lavoie-Courchesne S, Dickinson P, Lerch JP, Zijdenbos AP and Evans AC (2012)
% The pipeline system for Octave and Matlab (PSOM): a lightweight scripting framework and
% execution engine for scientific workflows. Front. Neuroinform. 6:7.
% doi: 10.3389/fninf.2012.00007
%
% Cyril Pernet and Nicolas Chauveau 2012 wrote the version 1
% CP 24-06-2013 updated to be even more automatic + fix for new designs
% Cyril Pernet May 2014 - fully redesigned with a GUI and using psom
% Cyril Pernet and Ramon Martinez-Cancino, October 2014 updates for EEGLAB STUDY
% ----------------------------------------------------------------------
% Copyright (C) LIMO Team 2015

% programmer help
% ---------------
% we build a pipeline to import, buid the design and run the glm
% import - calls limo_batch_import_data
% design - calls limo_batch_design_matrix
% glm calls limo_eeg(4) or limo_eeg_tf(4)



opt.mode = 'session'; % run in the current session -- see psom for other options
opt.max_queued = Inf; % with a maximum of possible sessions
opt.time_between_checks = 2; % and 2sec between job submission
opt.flag_pause = false; % don't bother asking to start jobs
opt.flag_debug = true; % report a bit more of issues
psom_gb_vars


%% what to do

if nargin == 0
    option = questdlg('batch mode','option','model specification','contrast only','both','model specification');
    if isempty(option)
        return
    end
    
    % model
    if strcmp(option,'model specification') || strcmp(option,'both')
        [model.set_files,model.cat_files,model.cont_files,model.defaults]=limo_batch_gui;
        if isempty(model.set_files)
            return
        end
    end
    % contrast
    if strcmp(option,'both')
        [FileName,PathName,FilterIndex]=uigetfile({'*.mat','MAT-files (*.mat)'; ...
            '*.txt','Text (*.txt)'}, 'Pick a matrix of contrasts');
        if FilterIndex ~=0
            if strcmp(FileName(end-3:end),'.txt')
                batch_contrast.mat = importdata(FileName);
            elseif strcmp(FileName(end-3:end),'.mat')
                FileName = load([PathName FileName]);
                batch_contrast.mat = getfield(FileName,cell2mat(fieldnames(FileName)));
            end
        else
            disp('limo batch aborded')
        end
        
        % update paths
        for f=1:size(model.set_files,1)
            [root,~,~] = fileparts(model.set_files{f});
            folder = ['GLM_' model.defaults.analysis];
            batch_contrast.LIMO_files{f} = [root filesep folder filesep 'LIMO.mat'];
        end
    end
    
    if strcmp(option,'contrast only')
        
        % get paths
        [FileName,PathName,FilterIndex]=uigetfile({'*.txt','Text (*.txt)'; ...
            '*.mat','MAT-files (*.mat)'}, 'Pick a list of LIMO.mat files');
        if FilterIndex ~=0
            if strcmp(FileName(end-3:end),'.txt')
                batch_contrast.LIMO_files = importdata(fullfile(PathName, FileName));
            elseif strcmp(FileName(end-3:end),'.mat')
                FileName = load([PathName FileName]);
                batch_contrast.LIMO_files = getfield(FileName,cell2mat(fieldnames(FileName)));
            end
        else
            disp('limo batch aborded')
        end
        
        % get the constrasts
        [FileName,PathName,FilterIndex]=uigetfile({'*.mat','MAT-files (*.mat)'; ...
            '*.txt','Text (*.txt)'}, 'Pick a matrix of contrasts');
        if FilterIndex ~=0
            if strcmp(FileName(end-3:end),'.txt')
                batch_contrast.mat = importdata(FileName);
            elseif strcmp(FileName(end-3:end),'.mat')
                FileName = load([PathName FileName]);
                batch_contrast.mat = getfield(FileName,cell2mat(fieldnames(FileName)));
            end
        else
            disp('limo batch aborded')
        end
    end
else
    option = varargin{1};
    
    % model
    if strcmp(option,'model specification') || strcmp(option,'both')
        model = varargin{2};
    end
    
    % batch_contrast
    if strcmp(option,'contrast only') || strcmp(option,'both')
        batch_contrast = varargin{3};
        if ~isfield(contrast,'mat')
            errordlg('the field contrast.mat is missing'); return
        end
        
        if strcmp(option,'both') && ~isfield(contrast,'LIMO_files')
            for f=1:size(model.set_files,1)
                [root,~,~] = fileparts(model.set_files{f});
                folder = ['GLM_' model.defaults.analysis];
                batch_contrast.LIMO_files{f} = [root filesep folder filesep 'LIMO.mat'];
            end
            batch_contrast.LIMO_files = batch_contrast.LIMO_files';
        end
    end
end

if nargin == 4
    STUDY = varargin{4}; clear varargin{4};
    cd(STUDY.filepath); current =pwd;
    mkdir('limo_batch_report'); 
    mkdir(['LIMO_' STUDY.filename(1:end-6)]);
    study_root = [STUDY.filepath filesep ['LIMO_' STUDY.filename(1:end-6)]];
    LIMO_files.LIMO = study_root;
else
    current =pwd;
    mkdir('limo_batch_report')
end

%% -------------------------------------
%% build pipelines
%% -------------------------------------

if strcmp(option,'model specification') || strcmp(option,'both')
    % quick check
    if ~isempty(model.cat_files)
        if size(model.cat_files,1) ~= size(model.set_files,1)
            error('the number of set and cat files disagree')
        end
    end
    
    if ~isempty(model.cont_files)
        if size(model.cont_files,1) ~= size(model.set_files,1)
            error('the number of set and cat files disagree')
        end
    end
    
    % build the pipelines
    for subject = 1:size(model.set_files,1)
        
        % build LIMO.mat files from import
        command = 'limo_batch_import_data(files_in,opt.cat,opt.cont,opt.defaults)';
        pipeline(subject).import.command = command;
        pipeline(subject).import.files_in = model.set_files{subject};
        pipeline(subject).import.opt.defaults = model.defaults;

        if isfield(model.defaults,'type')
            pipeline(subject).import.opt.defaults.type = model.defaults.type;
        else
            pipeline(subject).import.opt.defaults.type = 'Channels';
        end
        
        if isfield(model.defaults,'method')
            pipeline(subject).import.opt.defaults.method = model.defaults.method;
        else
            pipeline(subject).import.opt.defaults.method = 'WLS';
        end
        
        if nargin == 4
            mkdir([study_root filesep cell2mat(STUDY.names(subject))]);
            root = [study_root filesep cell2mat(STUDY.names(subject))];
            glm_name = ['GLM' num2str(STUDY.design_index) model.defaults.method '_' model.defaults.analysis '_' model.defaults.type];
            batch_contrast.LIMO_files{subject} = [root filesep glm_name filesep 'LIMO.mat']; 
            pipeline(subject).import.opt.defaults.studyinfo = STUDY.design_info;
        else
            [root,~,~] = fileparts(model.set_files{subject});
            glm_name = ['GLM_' model.defaults.method '_' model.defaults.analysis '_' model.defaults.type];    
        end
        pipeline(subject).import.files_out = [root filesep glm_name filesep 'LIMO.mat'];
        
        if ~isempty(model.cat_files)
            pipeline(subject).import.opt.cat = model.cat_files{subject};
        else
            pipeline(subject).import.opt.cat = [];
        end
        if ~isempty(model.cont_files)
            pipeline(subject).import.opt.cont = model.cont_files{subject};
        else
            pipeline(subject).import.opt.cont = [];
        end
        pipeline(subject).import.opt.defaults.name = fileparts(pipeline(subject).import.files_out);
        LIMO_files.mat{subject}  = [root filesep glm_name filesep 'LIMO.mat'];
        LIMO_files.Beta{subject} = [root filesep glm_name filesep 'Betas.mat'];
       
        % make design and evaluate
        command = 'limo_batch_design_matrix(files_in)';
        pipeline(subject).design.command = command;
        pipeline(subject).design.files_in = pipeline(subject).import.files_out;
        pipeline(subject).design.files_out = [root filesep glm_name filesep 'Yr.mat'];
        
        % run GLM
        if strcmp(model.defaults.analysis,'Time') || strcmp(model.defaults.analysis,'Frequency');
            command = 'cd(fileparts(files_in)), limo_eeg(4)';
        else strcmp(model.defaults.analysis,'Time-Frequency');
            command = 'cd(fileparts(files_in)), limo_eeg_tf(4)';
        end
        pipeline(subject).glm.command = command;
        pipeline(subject).glm.files_in = pipeline(subject).import.files_out;
        pipeline(subject).glm.files_out = [root filesep glm_name filesep 'Betas.mat'];
    end
    
end

if strcmp(option,'contrast only') || strcmp(option,'both')
    
    for subject = 1:size(batch_contrast.LIMO_files,1)
        command = 'limo_batch_contrast(files_in,opt.C)';
        pipeline(subject).batch_contrast.command = command;
        pipeline(subject).batch_contrast.files_in = batch_contrast.LIMO_files{subject};
        if iscell(batch_contrast.mat)
            pipeline(subject).batch_contrast.opt.C = cell2mat(batch_contrast.mat);
        else
            pipeline(subject).batch_contrast.opt.C = batch_contrast.mat;
        end
        
        for c=1:size(batch_contrast.mat,1)
            name{c} = [PathName filesep 'con_' num2str(c) '.mat'];
        end
        pipeline(subject).batch_contrast.files_out = name{1};
        LIMO_files.con{subject} = name;
    end
end

%% -------------------------------------
%% run the analyses
%% -------------------------------------

% run pipelines and report
try
    N = size(model.set_files,1);
    LIMO_files.mat = LIMO_files.mat';
    LIMO_files.Beta = LIMO_files.Beta';
    remove_limo = zeros(1,N);
catch
    N = size(batch_contrast.LIMO_files,1);
end

if isfield(LIMO_files,'con')
    LIMO_files.con = LIMO_files.con';
    remove_con = zeros(1,N);
end

for subject = 1:N
    disp('--------------------------------')
    fprintf('processing subject %g/%g \n',subject,N)
    disp('--------------------------------')
    try
        opt.path_logs = [current filesep 'limo_batch_report' filesep 'subject' num2str(subject)];
        psom_run_pipeline(pipeline(subject),opt)
        report{subject} = ['subject ' num2str(subject) ' processed'];
    catch ME
        report{subject} = ['subject ' num2str(subject) ' failed'];
        if strcmp(option,'model specification') 
            remove_limo(subject) = 1;
        elseif strcmp(option,'both')
            remove_limo(subject) = 1;
            remove_con(subject) = 1;
        elseif strcmp(option,'contrast only')
            remove_con(subject) = 1;
        end
    end
end

% save as txt file the list of .set, Betas, LIMO and con
if exist('STUDY','var')
    cd(LIMO_files.LIMO)
    cell2csv('EEGLAB_set.txt',model.set_files)
else
    cd(current)
end

if strcmp(option,'model specification') || strcmp(option,'both')
    cell2csv('LIMO_files.txt', LIMO_files.mat(find(~remove_limo),:))
    cell2csv('Beta_files.txt', LIMO_files.Beta(find(~remove_limo),:))
end

if isfield(LIMO_files,'con')
    for c=1:size(batch_contrast.mat,1)
        index = 1;
        for subject = 1:N
            name{index} = [fileparts(pipeline(subject).glm.files_out) filesep 'con_' num2str(c) '.mat'];
            index = index + 1;
        end
        name = name';
        cell2csv('con_files.txt', name(find(~remove_con),:));
    end   
end

% save the report from psom
cd([current filesep 'limo_batch_report'])
cell2csv('batch_report.txt', report')

cd(current); 
failed = 0;
for subject=1:N; 
    if strfind(report{subject},'failed')
        failed = 1;
    end
end

if failed == 0
    disp('LIMO batch processing finished succesfully')
else
    disp('LIMO batch done, some errors where detected see report')
end


