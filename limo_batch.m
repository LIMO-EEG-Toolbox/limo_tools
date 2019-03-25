function [LIMO_files, procstatus] = limo_batch(varargin)

% interactive function to run several 1st level analyses
% select directories and files - possibly enter contrasts of
% interests and let it run. The batch relies on PSOM (see Ref)
% see opt.mode for parallel computing on grid using qsub or msub
% <https://github.com/PSOM>
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
%       model.defaults.type_of_analysis 'univariate' or 'multivariate'
%       model.defaults.fullfactorial 0/1
%       model.defaults.zscore 0/1
%       model.defaults.start starting time in ms
%       model.defaults.end ending time in ms
%       model.defaults.lowf starting point in Hz
%       model.defaults.highf ending point in Hz
%       model.defaults.bootstrap 0/1
%       model.defaults.tfce 0/1
%       model.defaults.neighbouring_matrix neighbouring matrix use for clustering (necessary if bootstrap = 1)
%
%       contrast is a structure that specify which contrasts to run for which subject
%       contrast.LIMO_files: a list of LIMO.mat (full path) for the different subjects
%                            this is optional if option 'both' is selected
%       contrast.mat: a matrix of contrasts to run (assumes the same for all subjects)
%       eeglab_study is the STUDY structure allowing to create multiple design with consistant names etc ... 
%
% OUTPUT  
% LIMO_files  - A cell array of LIMO.mat (info about subjects' GLM)
%               create a directory per subject with GLM results in it
%               create a log file directory with the pipleine and logs
% procstatus  - [1 x Number of subjects] binary vector. Status of the LIMO computations for each of the N subjects.
%               [0] Failed, [1] Processed.
%
% see also limo_eeg limo_import_t limo_import_f limo_import_tf 
% see also psom in external folder
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
% Copyright (C) LIMO Team 2018

% programmer help
% ---------------
% we build a pipeline to import, buid the design and run the glm
% import - calls limo_batch_import_data
% design - calls limo_batch_design_matrix
% glm calls limo_eeg(4) or limo_eeg_tf(4)

opt.mode                = 'session'; % run in the current session -- see psom for other options // in batch we use parfor
opt.max_queued          = Inf; % with a maximum of possible sessions
opt.time_between_checks = 2; % and 2sec between job submission
opt.flag_pause          = false; % don't bother asking to start jobs
opt.flag_debug          = true; % report a bit more of issues
psom_gb_vars

% Initializing Outputs
LIMO_files = [];
procstatus = [];

%% what to do

if nargin <= 1
    
    if nargin == 0
        option = questdlg('batch mode','option','model specification','contrast only','both','model specification');
        if isempty(option)
            return
        end
    else
        option = varargin{1};
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
            disp('limo batch aborded'); return
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
            disp('limo batch aborded'); return
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
            disp('limo batch aborded'); return
        end
    end
elseif nargin > 1
    option = varargin{1};
    
    % model
    if strcmp(option,'model specification') || strcmp(option,'both')
        model = varargin{2};
    end
    
    % batch_contrast
    if strcmp(option,'contrast only') || strcmp(option,'both')
        batch_contrast = varargin{3};
        if ~isfield(batch_contrast,'mat')
            errordlg('the field batch_contrast.mat is missing'); return
        end
    end
end

if nargin == 4
    STUDY = varargin{4}; 
    if isempty(STUDY.filepath)
        STUDY.filepath =pwd;
    end
    cd(STUDY.filepath); % go to study

    % make derivatives directory
    if isempty(findstr(STUDY.filepath,'derivatives'))
        if ~exist([STUDY.filepath filesep 'derivatives'],'dir')
            mkdir([STUDY.filepath filesep 'derivatives']);
        end
        cd('derivatives'); % if study is not in 'derivatives' go to it
    end
    current = pwd; 
    
    % make the batch report and study directories
    if exist('limo_batch_report','dir')               ~= 7, mkdir('limo_batch_report'); end
    if exist(['LIMO_' STUDY.filename(1:end-6)],'dir') ~= 7, mkdir(['LIMO_' STUDY.filename(1:end-6)]); end
    study_root = [current filesep ['LIMO_' STUDY.filename(1:end-6)]];
    LIMO_files.LIMO = study_root;
    
    % if clustering is used, check subjects' set ordering
    if isfield(model.defaults, 'icaclustering')
        unique_subjects  = STUDY.design(STUDY.currentdesign).cases.value'; % all designs have the same cases
        for s = 1:length(unique_subjects)
            order{s} = eval(unique_subjects{s}); % find(strcmp(unique_subjects{s},STUDY.subject));
        end
    end
    
else
    current = pwd;
    mkdir('limo_batch_report')
end

clear varargin

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
    for s = 1:size(model.set_files,1)
        if isfield(model.defaults, 'icaclustering')
            subject = order{s}; % allows picking up the order in set_files based in STUDY
        else
            subject = s;
        end
        
        % build LIMO.mat files from import
        command = 'limo_batch_import_data(files_in,opt.cat,opt.cont,opt.defaults)';
        pipeline(s).import.command = command;
        pipeline(s).import.files_in = model.set_files{s};
        pipeline(s).import.opt.defaults = model.defaults;

        if isfield(model.defaults,'type')
            pipeline(s).import.opt.defaults.type = model.defaults.type;
        else
            pipeline(s).import.opt.defaults.type = 'Channels';
        end
        
        if isfield(model.defaults,'method')
            pipeline(s).import.opt.defaults.method = model.defaults.method;
        else
            pipeline(s).import.opt.defaults.method = 'WLS';
        end
        
        if isfield(model.defaults,'type_of_analysis')
            pipeline(s).import.opt.defaults.type_of_analysis = model.defaults.type_of_analysis;
        else
            pipeline(s).import.opt.defaults.type_of_analysis = 'Mass-univariate';
        end
        
        if nargin == 4
            if ~isempty(findstr(model.set_files{s},'sub'))
                root = [study_root filesep 'sub-' num2str(subject)];
            else
                root = [study_root filesep num2str(subject)];
            end
            % root = [study_root filesep 'sub-' num2str(subject)];
            if exist(root,'dir') ~= 7; mkdir(root); end
            design_name = STUDY.design(STUDY.currentdesign).name; 
            design_name(isspace(design_name)) = [];
            if findstr(design_name,'STUDY.')
                design_name = design_name(7:end);
            end
            glm_name = [design_name '_GLM_' model.defaults.type '_' model.defaults.analysis '_' model.defaults.method];
            batch_contrast.LIMO_files{s} = [root filesep glm_name filesep 'LIMO.mat']; 
            % pipeline(subject).import.opt.defaults.studyinfo = STUDY.design_info;
        else
            [root,~,~] = fileparts(model.set_files{subject});
            glm_name = ['GLM_' model.defaults.method '_' model.defaults.analysis '_' model.defaults.type];    
        end
        pipeline(s).import.files_out = [root filesep glm_name filesep 'LIMO.mat'];
        
        if strcmp(option,'both') && ~isfield(batch_contrast,'LIMO_files')
            batch_contrast.LIMO_files{s} = [root filesep glm_name filesep 'LIMO.mat'];
            batch_contrast.LIMO_files = batch_contrast.LIMO_files';
        end

        if ~isempty(model.cat_files)
            pipeline(s).import.opt.cat = model.cat_files{s};
        else
            pipeline().import.opt.cat = [];
        end
        
        if ~isempty(model.cont_files)
            pipeline(s).import.opt.cont = model.cont_files{s};
        else
            pipeline(s).import.opt.cont = [];
        end
        
        pipeline(s).import.opt.defaults.name = fileparts(pipeline(s).import.files_out);
        LIMO_files.mat{s}  = [root filesep glm_name filesep 'LIMO.mat'];
        LIMO_files.Beta{s} = [root filesep glm_name filesep 'Betas.mat'];
       
        % make design and evaluate
        command = 'limo_batch_design_matrix(files_in)';
        pipeline(s).design.command = command;
        pipeline(s).design.files_in = pipeline(s).import.files_out;
        pipeline(s).design.files_out = [root filesep glm_name filesep 'Yr.mat'];
        
        % run GLM
        if strcmp(model.defaults.analysis,'Time') || strcmp(model.defaults.analysis,'Frequency')
            command = 'cd(fileparts(files_in)), limo_eeg(4)';
        elseif strcmp(model.defaults.analysis,'Time-Frequency')
            command = 'cd(fileparts(files_in)), limo_eeg_tf(4)';
        end
        pipeline(s).glm.command = command;
        pipeline(s).glm.files_in = pipeline(s).import.files_out;
        pipeline(s).glm.files_out = [root filesep glm_name filesep 'Betas.mat'];
    end
    
end

if strcmp(option,'contrast only') || strcmp(option,'both')
    for s = 1:length(batch_contrast.LIMO_files)
        command = 'limo_batch_contrast(files_in,opt.C)';
        pipeline(s).n_contrast.command = command;
        pipeline(s).n_contrast.files_in = batch_contrast.LIMO_files{s};
        if iscell(batch_contrast.mat)
            pipeline(s).n_contrast.opt.C = cell2mat(batch_contrast.mat);
        else
            pipeline(s).n_contrast.opt.C = batch_contrast.mat;
        end
        
        if strcmp(option,'both') % we can only be sure of the number if it's a new model
            for c=1:size(batch_contrast.mat,1)
                name{c} = [fileparts(batch_contrast.LIMO_files{s}) filesep 'con_' num2str(c) '.mat'];
            end
            pipeline(s).n_contrast.files_out = name; % name{1};
            LIMO_files.con{s} = name;
        end
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
procstatus = zeros(1,N);

if isfield(LIMO_files,'con')
    LIMO_files.con = LIMO_files.con';
    remove_con = zeros(1,N);
end


% ----------------------
%% Save pipeline
% useful to re-run, simply calling psom_run_pipeline
if ~exist('glm_name','var') && strcmp(option,'contrast only') 
    [~,glm_name]=fileparts(fileparts(pipeline(1).n_contrast.files_in));
end
save([current filesep 'limo_pipeline_' glm_name '.mat'],'pipeline')

% allocate names
for subject = 1:N
    limopt{subject}= opt;
    limopt{subject}.path_logs = [current filesep 'limo_batch_report' filesep glm_name filesep 'subject' num2str(order{s})];
end
    
parfor subject = 1:N
    disp('--------------------------------')
    fprintf('processing subject %g/%g \n',order{subject},N)
    disp('--------------------------------')
    try
        psom_run_pipeline(pipeline(subject),limopt{subject})
        % example of debugging 
        % ---------------------
        % psom reported with function failed, eg limo_batch_import
        % pipeline(subject).import tells you the command line to test
        % put the point brack where needed and call e.g.
        % limo_batch_import_data(pipeline(subject).import.files_in,pipeline(subject).import.opt.cat,pipeline(subject).import.opt.cont,pipeline(subject).import.opt.defaults)
        % limo_batch_design_matrix(pipeline(subject).design.files_in)
        report{subject} = ['subject ' order{subject} ' processed'];
        procstatus(subject) = 1;
    catch ME
        report{subject} = ['subject ' order{subject} ' failed'];
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

%% Save txt files
% save as txt file the list of .set, Betas, LIMO and con
% these lists can then be used in second level analyses

if exist('STUDY','var')
    cd(LIMO_files.LIMO)
    cell2csv(['EEGLAB_set_' glm_name '.txt'],model.set_files)
else
    cd(current)
end

if strcmp(option,'model specification') || strcmp(option,'both')
    cell2csv(['LIMO_files_' glm_name '.txt'], LIMO_files.mat(find(~remove_limo),:))
    cell2csv(['Beta_files_' glm_name '.txt'], LIMO_files.Beta(find(~remove_limo),:))
end

if strcmp(option,'contrast only') || strcmp(option,'both')
    for c=1:size(batch_contrast.mat,1)
        index = 1; clear name
        for subject = 1:N
            if strcmp(option,'contrast only')
                load([fileparts(pipeline(subject).n_contrast.files_in) filesep 'LIMO.mat']);
                for l=1:size(LIMO.contrast,2)
                    if isequal(LIMO.contrast{l}.C,batch_contrast.mat(c,:))
                        con_num = l; break
                    end
                end
                name{index} = [fileparts(pipeline(subject).n_contrast.files_in) filesep 'con_' num2str(con_num) '.mat'];
            else
                name{index} = [fileparts(pipeline(subject).glm.files_out) filesep 'con_' num2str(c) '.mat'];
                con_num = c;
            end
            index = index + 1;
        end
        name = name';
        
        if sum(remove_con) ~= 0
            cell2csv(['con' num2str(con_num) '_files_' glm_name '.txt'], name(find(~remove_con),:));
        else
            cell2csv(['con' num2str(con_num) '_files_'  glm_name '.txt'], name);
        end
    end
end
% save the report from psom
cd([current filesep 'limo_batch_report'])
cell2csv(['batch_report_' glm_name '.txt'], report')

cd(current); 
failed = 0;
for subject=1:N
    if strfind(report{subject},'failed')
        failed = 1;
    end
end

if failed == 0
    disp('LIMO batch processing finished succesfully')
else
    disp('LIMO batch done, some errors where detected see report')
end
disp('LIMO batch works thanks to PSOM by Bellec et al. (2012)')
disp('The Pipeline System for Octave and Matlab. Front. Neuroinform. 6:7')


