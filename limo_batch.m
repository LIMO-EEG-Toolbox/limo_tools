function [LIMO_files, procstatus] = limo_batch(varargin)

% interactive function to run several 1st level analyses
% select directories and files - possibly enter contrasts of
% interests and let it run. The batch relies on PSOM (see Ref)
% see opt.mode for parallel computing on grid using qsub or msub
% <https://github.com/PSOM>
%
% FORMAT limo_batch
% [LIMO_files, procstatus] = limo_batch(option,model,contrast)
% [LIMO_files, procstatus] = limo_batch(option,model,contrast,eeglab_study)
%
% INPUT if empty uses GUI
%       - option should be 'model specification', 'contrast only' or 'both'
%       - model is a structure that specifiy information to build a model
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
%       - contrast is a structure that specify which contrasts to run for which subject
%       contrast.LIMO_files: a list of LIMO.mat (full path) for the different subjects
%                            this is optional if option 'both' is selected
%       contrast.mat: a matrix of contrasts to run (assumes the same for all subjects, rows are contrasts,
%                     columns are variables in the GLM including the constant)
%       - eeglab_study is the STUDY structure allowing to create multiple design with consistant names etc ...
%
% OUTPUT
% LIMO_files  - A cell array of LIMO.mat (info about subjects' GLM)
%               create a directory per subject with GLM results in it
%               create a log file directory with the pipleine and logs
% procstatus  - [1 x Number of subjects] binary vector. Status of the LIMO computations for each of the N subjects.
%               [0] Failed, [1] Processed.
%
%
% Example: limo_batch('both',model,contrast,STUDY);
%                            model.defaults.datatype= 'erp'
%                            model.defaults.type= 'Channels'
%                            model.defaults.analysis= 'Time'
%                            model.defaults.start= -1000
%                            model.defaults.end= 1996
%                            model.defaults.lowf= []
%                            model.defaults.highf= []
%                            model.defaults.fullfactorial= 0
%                            model.defaults.zscore= 1
%                            model.defaults.bootstrap= 0
%                            model.defaults.tfce= 0
%                            model.defaults.method= 'WLS'
%                            model.defaults.Level= 1
%                            model.defaults.type_of_analysis= 'Mass-univariate'
%                            model.cat_files: {n×1 cell};
%                                  model.cat_files{n}' = [1 1 1 2 2 3 4 4 2 3 ....];
%                            model.cont_files: []
%                            model.set_files: {n×1 cell}
%                                  model.set_files{n} = 'D:\EEG\mysuperdataset\sub-001\sub-001_task_dostuff.set';
%                            contrast.mat = [1 0 -1 0 0 ; 0 1 0 -1 0];
%
%          limo_batch('contrast only',[],contrast,STUDY);
%                            contrast.LIMO_files = 'D:\EEG\mysuperdataset\derivatives\LIMOstuff\LIMO_files.txt';
%                            contrast.mat = [0 1 0 1 ; 1 0 1 0; 1 -1 1 -1];
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
%  Copyright (C) LIMO Team 2022

% programmer help
% ---------------
% we build a pipeline to import, buid the design and run the glm
% import - calls limo_batch_import_data
% design - calls limo_batch_design_matrix
% glm calls limo_eeg(4) or limo_eeg_tf(4)

opt.mode                = 'session'; % run in the current session -- see psom for other options // in batch we use parfor
opt.max_queued          = Inf; % with a maximum of possible sessions
opt.time_between_checks = 3; % and x sec between job submission
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
            procstatus = 'import aborded';
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
                % batch_contrast.mat = getfield(FileName,cell2mat(fieldnames(FileName)));
                batch_contrast.mat = FileName.(cell2mat(fieldnames(FileName)));
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
        limo_settings_script;
        FileName = '';
        if ~isempty(limo_settings.workdir)
            fileList = dir(fullfile(limo_settings.workdir, 'LIMO_*', 'LIMO_*.txt'));
            if ~isempty(fileList)
                for iFile = 1:length(fileList)
                    fileList(iFile).fullname = fullfile(fileList(iFile).folder,fileList(iFile).name);
                end
                uiList = { {'style' 'text' 'string' 'Pick a 1st level analysis file' } ...
                           { 'style' 'popupmenu' 'string' {fileList.name} } };
                res = inputgui('uilist', uiList, 'geometry', {1 1}, 'cancel', 'Browse');
                if ~isempty(res)
                    FileName = fileList(res{1}).name;
                    PathName = fileList(res{1}).folder;
                    FilterIndex = 1;
                end
            end
        end
        if isempty(FileName)
            [FileName,PathName,FilterIndex]=uigetfile({'*.txt','Text (*.txt)'; ...
                '*.mat','MAT-files (*.mat)'}, 'Pick a list of LIMO.mat files');
        end
        

        if FilterIndex ~=0
            if strcmp(FileName(end-3:end),'.txt')
                batch_contrast.LIMO_files = importdata(fullfile(PathName, FileName));
            elseif strcmp(FileName(end-3:end),'.mat')
                FileName = load([PathName FileName]);
                % batch_contrast.LIMO_files = getfield(FileName,cell2mat(fieldnames(FileName)));
                batch_contrast.LIMO_files = FileName.(cell2mat(fieldnames(FileName)));
            end
            LIMO_files.LIMO = PathName;
        else
            disp('limo batch aborded'); return
        end
        
        % get the constrasts
        limo_settings_script;
        if isempty(limo_settings.workdir)
            [FileName,PathName,FilterIndex]=uigetfile('*.*', 'Pick a matrix of contrasts');
            if FilterIndex ~=0
                if strcmp(FileName(end-3:end),'.txt')
                    batch_contrast.mat = importdata(FileName);
                elseif strcmp(FileName(end-3:end),'.mat')
                    FileName = load([PathName FileName]);
                    % batch_contrast.mat = getfield(FileName,cell2mat(fieldnames(FileName)));
                    batch_contrast.mat = FileName.(cell2mat(fieldnames(FileName)));
                end
            else
                disp('limo batch aborded'); return
            end
        else
            batch_contrast.mat = limo_contrast_manager(batch_contrast.LIMO_files{1});
            if isempty(batch_contrast.mat)
                disp('limo batch aborded'); return
            end
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
        
        if isfield(batch_contrast,'LIMO_files')
            if ~iscell(batch_contrast.LIMO_files)
                if strcmpi(batch_contrast.LIMO_files(end-3:end),'.txt')
                    batch_contrast.LIMO_files = importdata(batch_contrast.LIMO_files);
                else
                   error('contrast.LIMO_files must be either the path to a txt file or a cell array of file(s)') 
                end
            end
        end

        if ~isfield(batch_contrast,'mat')
            errordlg('the field batch_contrast.mat is missing'); return
        end
    end
end

% check EEGLAB STUDY
if nargin == 4
    STUDY = varargin{4}; clear varargin{4};
end

% not passed but in base workspace (case of batching contrast from GUI)
if ~exist('STUDY','var') && evalin('base', 'exist(''STUDY'',''var'')')
    STUDY = evalin('base', 'STUDY');
    if ~isstruct(STUDY); clear STUDY; end
end

if isempty(STUDY)
    clear STUDY
end

if exist('STUDY','var')
    if isempty(STUDY.filepath)
        STUDY.filepath =pwd;
    end
    cd(STUDY.filepath); % go to study
    current = pwd;

    [~,foldername] = fileparts(STUDY.filepath);
    if ~strcmp(foldername,'derivatives')
        % derivatives should have been created by std_limo if not already in the path
        if exist(['derivatives' filesep 'LIMO_' STUDY.filename(1:end-6)],'dir') ~= 7
            mkdir(['derivatives' filesep 'LIMO_' STUDY.filename(1:end-6)]);
        end
        if exist(['derivatives' filesep 'LIMO_' STUDY.filename(1:end-6) filesep 'limo_batch_report'],'dir') ~= 7
            mkdir(['derivatives' filesep 'LIMO_' STUDY.filename(1:end-6) filesep 'limo_batch_report']);
        end
        LIMO_files.LIMO = [current filesep ['derivatives' filesep 'LIMO_' STUDY.filename(1:end-6)]];
    else
        if exist(['LIMO_' STUDY.filename(1:end-6)],'dir') ~= 7
            mkdir(['LIMO_' STUDY.filename(1:end-6)]);
        end
        if exist(['LIMO_' STUDY.filename(1:end-6) filesep 'limo_batch_report'],'dir') ~= 7
            mkdir(['LIMO_' STUDY.filename(1:end-6) filesep 'limo_batch_report']);
        end
        LIMO_files.LIMO = [current filesep ['LIMO_' STUDY.filename(1:end-6)]];
    end
else % if not part of a EEGLAB STUDY - e.g. run locally or FieldTrip
    [~,foldername] = fileparts(pwd);
    if ~strcmp(foldername,'derivatives') % make a derivatives folder
        mkdir('derivatives'); cd derivatives
    end
    current = pwd;
    mkdir('limo_batch_report')
    if isempty(LIMO_files)
        LIMO_files.LIMO = current;
    end
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
        pipeline(subject).import.command = command; %#ok<*AGROW>
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
        
        if isfield(model.defaults,'type_of_analysis')
            pipeline(subject).import.opt.defaults.type_of_analysis = model.defaults.type_of_analysis;
        else
            pipeline(subject).import.opt.defaults.type_of_analysis = 'Mass-univariate';
        end
        
        if exist('STUDY','var')
            if ~contains(STUDY.datasetinfo(subject).filename,{'sub-'}) && ...
                    ~contains(STUDY.datasetinfo(subject).filename,{'_task-'}) % not bids
                root = [fileparts(LIMO_files.LIMO) filesep 'sub-' STUDY.datasetinfo(subject).subject];
            else
                subname = STUDY.datasetinfo(subject).subject;
                extra   = STUDY.datasetinfo(subject).filepath(strfind(STUDY.datasetinfo(subject).filepath,subname)+length(subname):end);
                root    = [fileparts(LIMO_files.LIMO) filesep subname extra]; % still in derivatives via LIMO_files.LIMO
            end
            
            % if session and data are not in a derivatives/sess, make subdir
            if ~isempty(STUDY.datasetinfo(subject).session)
                nsess = sum(strcmp(STUDY.datasetinfo(subject).subject,{STUDY.datasetinfo.subject}));
                if ~contains(root,'ses-') && nsess>=1
                    if ischar(STUDY.datasetinfo(subject).session)
                        reuse = dir(fullfile(root,['ses-*' STUDY.datasetinfo(subject).session]));
                        if ~isempty(reuse)
                            index = find(arrayfun(@(x) STUDY.datasetinfo(subject).session == eval(x.name(5:end)), reuse));
                            root = fullfile(reuse(index).folder,reuse(index).name);
                        else
                            root  = fullfile(root,['ses-' STUDY.datasetinfo(subject).session]);
                        end
                    else
                        reuse = dir(fullfile(root,['ses-*' num2str(STUDY.datasetinfo(subject).session)]));
                        if ~isempty(reuse)
                            index = find(arrayfun(@(x) STUDY.datasetinfo(subject).session == eval(x.name(5:end)), reuse));
                            root = fullfile(reuse(index).folder,reuse(index).name);
                        else
                            root = fullfile(root,['ses-' num2str(STUDY.datasetinfo(subject).session)]);
                        end
                    end
                end
            end
            
            % [root filesep eeg] - case of bids without ses-
            if exist(fullfile(root,'eeg'),'dir')
                root = fullfile(root,'eeg');
            end

            if exist(root,'dir') ~= 7
                mkdir(root);
            end
            design_name = STUDY.design(STUDY.currentdesign).name;
            design_name(isspace(design_name)) = [];
            if strfind(design_name,'STUDY.') %#ok<STRIFCND>
                design_name = design_name(7:end);
            end
            glm_name = [STUDY.filename(1:end-6) '_' design_name '_GLM_' model.defaults.type '_' model.defaults.analysis '_' model.defaults.method];
            batch_contrast.LIMO_files{subject} = [root filesep glm_name filesep 'LIMO.mat'];
            % pipeline(subject).import.opt.defaults.studyinfo = STUDY.design_info;
        else
            [root,~,~] = fileparts(model.set_files{subject});
            for l=min(length(LIMO_files.LIMO),length(root)):-1:1
                common(l) = root(l) == LIMO_files.LIMO(l);
            end
            root = fullfile(LIMO_files.LIMO,root(min(find(diff(common))):end)); 
            glm_name = ['GLM_' model.defaults.method '_' model.defaults.analysis '_' model.defaults.type];
        end
        pipeline(subject).import.files_out = [root filesep glm_name filesep 'LIMO.mat'];
        subname = limo_get_subname(pipeline(subject).import.files_in);
        if ~isempty(subname)
            subname = [subname '_desc-'];
        end

        if strcmp(option,'both') && ~isfield(batch_contrast,'LIMO_files')
            batch_contrast.LIMO_files{subject} = [root filesep glm_name filesep 'LIMO.mat'];
            batch_contrast.LIMO_files = batch_contrast.LIMO_files';
        end
        
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
        LIMO_files.Beta{subject} = [root filesep glm_name filesep subname 'Betas.mat'];
        
        % make design and evaluate
        command = 'limo_batch_design_matrix(files_in)';
        pipeline(subject).design.command = command;
        pipeline(subject).design.files_in = pipeline(subject).import.files_out;
        pipeline(subject).design.files_out = [root filesep glm_name filesep subname 'Yr.mat'];
        
        % run GLM
        command = 'limo_eeg(4,files_in)';
        pipeline(subject).glm.command = command;
        pipeline(subject).glm.files_in = pipeline(subject).import.files_out;
        pipeline(subject).glm.files_out = [root filesep glm_name filesep subname 'Betas.mat'];
    end
    
end

if strcmp(option,'contrast only') || strcmp(option,'both')
    
    if ~exist('model','var')
        model.defaults.bootstrap = 0;
        model.defaults.tfce      = 0;
    end
    
    for subject = 1:length(batch_contrast.LIMO_files)
        command = 'limo_batch_contrast(files_in,opt.C)';
        pipeline(subject).n_contrast.command = command;
        pipeline(subject).n_contrast.files_in = batch_contrast.LIMO_files{subject};
        if iscell(batch_contrast.mat)
            pipeline(subject).n_contrast.opt.C = cell2mat(batch_contrast.mat);
        else
            pipeline(subject).n_contrast.opt.C = batch_contrast.mat;
        end
        
        if exist(batch_contrast.LIMO_files{subject},'file')
            sub_LIMO = load(batch_contrast.LIMO_files{subject});
            if ~isfield(sub_LIMO.LIMO,'contrast')
                start = 0;
            else
                start = length(sub_LIMO.LIMO.contrast);
            end
        else
            start = 0;
        end
        
        subname = STUDY.datasetinfo(subject).subject;
        for c=1:size(batch_contrast.mat,1)
            name{c} = [fileparts(batch_contrast.LIMO_files{subject}) filesep subname '_desc-con_' num2str(c+start) '.mat'];
        end
        pipeline(subject).n_contrast.files_out = name; % name{1};
        LIMO_files.con{subject} = name;
    end
end

%% -------------------------------------
%% run the analyses
%% -------------------------------------

% run pipelines and report
if strcmp(option,'model specification') || strcmp(option,'both')
    N               = size(model.set_files,1);
    LIMO_files.mat  = LIMO_files.mat';
    LIMO_files.Beta = LIMO_files.Beta';
    remove_limo     = zeros(1,N);
else
    N               = length(batch_contrast.LIMO_files);
end
procstatus = zeros(1,N);

if isfield(LIMO_files,'con')
    LIMO_files.con = LIMO_files.con';
    remove_con     = zeros(1,N);
else
    remove_con     = 0;
end


% ----------------------
%% Save pipeline
% useful to re-run, simply calling psom_run_pipeline
if ~exist('glm_name','var') && strcmp(option,'contrast only')
    [~,glm_name]=fileparts(fileparts(pipeline(1).n_contrast.files_in));
end

if strcmp(option,'contrast only')
    save([LIMO_files.LIMO filesep 'limo_con_pipeline_' glm_name '.mat'],'pipeline')
else
    save([LIMO_files.LIMO filesep 'limo_pipeline_' glm_name '.mat'],'pipeline')
end

% allocate names
for subject = 1:N
    limopt{subject} = opt;
    limopt{subject}.path_logs = [LIMO_files.LIMO filesep 'limo_batch_report' filesep glm_name filesep 'subject' num2str(subject)];
end

limo_settings_script;
if model.defaults.bootstrap ~= 0 || ~limo_settings.psom % debugging mode, serial analysis
    
    for subject = 1:N 
        disp('--------------------------------')
        fprintf('processing model %g/%g \n',subject,N)
        disp('--------------------------------')
        
        psom_pipeline_debug(pipeline(subject));
        if strcmp(option,'contrast only')
            name = fileparts(batch_contrast.LIMO_files{subject}); 
        else
            [~,name]=fileparts(model.set_files{subject}); 
        end
        sub = min(strfind(name,'sub-'));
        ses = min(strfind(name,'ses-'));
        und = strfind(name,'_');
        if ~isempty(sub) && ~isempty(ses) && ~isempty(und)
            try
                sub_und = und(und>sub); ses_und = und(und>ses);
                report{subject} = ['subject ' name(sub+4:sub+min(abs(sub_und-sub))-1) ' session ' name(ses+4:ses+min(abs(ses_und-ses))-1) ' processed'];
            catch
                report{subject} = ['subject ' num2str(subject) ' processed'];
            end
        else
            report{subject} = ['subject ' num2str(subject) ' processed'];
        end
        procstatus(subject) = 1;
    end
    
else % parallel call to the pipeline , the usual way

    limo_check_ppool
    parfor subject = 1:N 
        disp('--------------------------------')
        fprintf('processing model %g/%g \n',subject,N)
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
            % limo_eeg(4,fileparts(pipeline(subject).glm.files_in))
            % limo_batch_contrast(pipeline(subject).n_contrast.files_in,pipeline(subject).n_contrast.opt.C)
            
            if strcmp(option,'contrast only')
                name = fileparts(batch_contrast.LIMO_files{subject}); %#ok<PFBNS>
            else
                [~,name]=fileparts(model.set_files{subject}); %#ok<PFBNS>
            end
            sub = min(strfind(name,'sub-'));
            ses = min(strfind(name,'ses-'));
            und = strfind(name,'_');
            
            if ~isempty(sub) && ~isempty(ses) && ~isempty(und)
                try
                    sub_und = und(und>sub); ses_und = und(und>ses);
                    if strcmp(option,'contrast only')
                        report{subject} = ['subject ' name(sub:sub+min(abs(sub_und-sub))-1) ' processed'];
                    else
                        report{subject} = ['subject ' name(sub+4:sub+min(abs(sub_und-sub))-1) ' session ' name(ses+4:ses+min(abs(ses_und-ses))-1) ' processed'];
                    end
                catch
                    report{subject} = ['subject ' num2str(subject) ' processed'];
                end
            else
                report{subject} = ['subject ' num2str(subject) ' processed'];
            end
            procstatus(subject) = 1;
        catch ME
            report{subject} = sprintf('subject %g failed: %s',subject,ME.message');
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
    try
        poolobj = gcp('nocreate'); 
        delete(poolobj); % close parallel pool;
    catch closepool
        disp('no parpool to close %s\n',closepool.message)
    end
end

%% Save txt files
cd(LIMO_files.LIMO)
if strcmp(option,'model specification') || strcmp(option,'both')
    % check that the dimensions of the model are identical otherwise tell the
    % user -- condition might be missing from design or from preprocessing
    good_subjects = find(remove_limo==0);
    for s = 1:length(good_subjects)
        LIMO = load(LIMO_files.mat{good_subjects(s)});
        design_dim(s) = size(LIMO.LIMO.design.X,2);
        names{s} = limo_get_subname(LIMO.LIMO.dir);
    end
    
    if exist("design_dim","var")
        if length(unique(design_dim)) > 1
            limo_warndlg(sprintf('Some subjects have different design dimensions.\nThis might be a feature of experiment (like encoding errors but some subjects have none) but is often not'))
            for s = 1:length(good_subjects)
                warning('subject %s design with %g regressors', names{s},design_dim(s));
            end
        end
    end

    % save as txt file the list of .set, Betas, LIMO and con
    % these lists can then be used in second level analyses
    if ~all(remove_limo)
        cell2csv([LIMO_files.LIMO filesep 'LIMO_files_' glm_name '.txt'], LIMO_files.mat(find(~remove_limo),:))
        cell2csv([LIMO_files.LIMO filesep 'Beta_files_' glm_name '.txt'], LIMO_files.Beta(find(~remove_limo),:))
    end
end

if strcmp(option,'contrast only') || strcmp(option,'both')
    for c=1:size(batch_contrast.mat,1)
        index = 1; clear name
        for subject = 1:N
            if strcmp(option,'contrast only')
                LIMO = load([fileparts(pipeline(subject).n_contrast.files_in) filesep 'LIMO.mat']); LIMO = LIMO.LIMO;
                if isfield(LIMO,'contrast')
                    con_num = max(find(cellfun(@(x) isequal(x.C,limo_contrast_checking(LIMO.dir,LIMO.design.X,batch_contrast.mat(c,:))),LIMO.contrast))); %#ok<*MXFND> % if several identical contrasts, take max
                else
                    con_num = c;
                end
                path = fileparts(pipeline(subject).n_contrast.files_in);
                name{index} = fullfile(path,[limo_get_subname(path) '_desc-con_' num2str(con_num) '.mat']);
            else
                con_num = c;
                path = fileparts(pipeline(subject).glm.files_out);
                name{index} = fullfile(path,[limo_get_subname(path) '_desc-con_' num2str(con_num) '.mat']);
            end
            index = index + 1;
        end
        name = name';
        
        if ~all(remove_con)
            cell2csv([LIMO_files.LIMO filesep 'con_' num2str(con_num) '_files_' glm_name '.txt'], name(find(~remove_con),:)); %#ok<*FNDSB>
        end
    end
end

% save the report from psom
cell2csv([LIMO_files.LIMO filesep 'limo_batch_report' filesep 'batch_report_' glm_name '.txt'], report')
failed = zeros(1,N);
for subject=1:N
    if strfind(report{subject},'failed')
        failed(subject) = 1;
    end
end

if sum(failed) == 0
    disp('LIMO batch processing finished succesfully')
else
    if sum(failed) == N % all subjects
        warning('LIMO batch done but all subjects failed. This can be a psom/disk access issue, try setting psom to false in limo_settings_script.m')
    else
        failed = find(failed);
        for s = 1:length(failed)
            warning('LIMO batch done, some errors where detected\ncheck limo batch report %s',limo_get_subname(fileparts(pipeline(failed(s)).n_contrast.files_in)));
        end
    end
end

% if EEGLAB STUDY check for groups and sessions 
% and further export txt files
if exist('STUDY','var')
    try
        if isfield(model, 'set_files')
            cell2csv([LIMO_files.LIMO filesep 'EEGLAB_set_' glm_name '.txt'],model.set_files)
        end
        
        if ~isempty(STUDY.datasetinfo(subject).session)
            sesvalues = unique(arrayfun(@(x) x.session, STUDY.datasetinfo));
        else
            sesvalues = 1;
        end
        
        % split txt files if more than 1 group or session
        if length(STUDY.group) > 1 || length(sesvalues)>1
            for s=1:length(sesvalues)
                for g= 1:length(STUDY.group)
                    if length(STUDY.group) > 1
                        subset = arrayfun(@(x)(strcmpi(x.group,STUDY.group{g})), STUDY.datasetinfo);
                    end
                    
                    if length(sesvalues) > 1
                        sesset = arrayfun(@(x) x.session==s, STUDY.datasetinfo);
                    end
                    
                    if isfield(LIMO_files,'mat') && isfield(LIMO_files,'Beta')
                        if length(STUDY.group) > 1 && length(sesvalues)==1 %#ok<*ISCL> % only groups
                            if any(subset)
                                cell2csv(fullfile(LIMO_files.LIMO, ['LIMO_files_Gp-' STUDY.group{g} '_' glm_name '.txt']), LIMO_files.mat(subset));
                                cell2csv(fullfile(LIMO_files.LIMO, ['Beta_files_Gp-' STUDY.group{g} '_' glm_name '.txt']), LIMO_files.Beta(subset));
                            end
                        elseif length(STUDY.group) == 1 && length(sesvalues) > 1 % only sessions
                            if any(sesset)
                                cell2csv(fullfile(LIMO_files.LIMO, ['LIMO_files_ses-' num2str(s) '_' glm_name '.txt']), LIMO_files.mat(sesset));
                                cell2csv(fullfile(LIMO_files.LIMO, ['Beta_files_ses-' num2str(s) '_' glm_name '.txt']), LIMO_files.Beta(sesset));
                            end
                        else % groups and sessions
                            if any(subset.*sesset)
                                cell2csv(fullfile(LIMO_files.LIMO, ['LIMO_files_ses-' num2str(s) '_Gp-' STUDY.group{g} '_' glm_name '.txt']), LIMO_files.mat(logical(subset.*sesset)));
                                cell2csv(fullfile(LIMO_files.LIMO, ['Beta_files_ses-' num2str(s) '_Gp-' STUDY.group{g} '_' glm_name '.txt']), LIMO_files.Beta(logical(subset.*sesset)));
                            end
                        end
                    end
                    
                    if isfield(LIMO_files,'con')
                        if length(STUDY.group) > 1 && length(sesvalues)==1 % only groups
                            tmpcell = LIMO_files.con(subset);
                            if ~isempty(tmpcell{1})
                                for c=1:length(tmpcell{1})
                                    [~,con_name,~] = fileparts(LIMO_files.con{1}{c});
                                    cell2csv(fullfile(LIMO_files.LIMO, [con_name '_files_Gp-' STUDY.group{g} '_' glm_name '.txt']),cellfun(@(x) x(c), tmpcell));
                                end
                            end
                        elseif length(STUDY.group) == 1 && length(sesvalues) > 1 % only sessions
                            tmpcell = LIMO_files.con(sesset);
                            if ~isempty(tmpcell{1})
                                for c=1:length(tmpcell{1})
                                    [~,con_name,~] = fileparts(LIMO_files.con{1}{c});
                                    cell2csv(fullfile(LIMO_files.LIMO, [con_name '_files_ses-' num2str(s) '_' glm_name '.txt']),cellfun(@(x) x(c), tmpcell));
                                end
                            end
                        else
                            tmpcell = LIMO_files.con(logical(subset.*sesset));
                            if ~isempty(tmpcell)
                                for c=1:length(tmpcell{1})
                                    [~,con_name,~] = fileparts(LIMO_files.con{1}{c});
                                    cell2csv(fullfile(LIMO_files.LIMO, [con_name '_files_ses-' num2str(s) '_Gp-' STUDY.group{g} '_' glm_name '.txt']),cellfun(@(x) x(c), tmpcell));
                                end
                            end
                        end
                    end
                end
            end
        end
    catch writtingerr
        if sum(failed) == 0
            warning(writtingerr.identifier,'all LIMO files created but failing to write some metadata txt files ''%s''\n ',writtingerr.message);
        else
            warning(writtingerr.identifier,'also failing to write some metadata txt files ''%s''\n ',writtingerr.message);
        end
    end
end

disp('LIMO batch works thanks to PSOM by Bellec et al. (2012)')
disp('The Pipeline System for Octave and Matlab. Front. Neuroinform. 6:7')
