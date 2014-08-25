function limo_batch(varargin)

% interactive function to run several 1st level analyses
% select directories and files - possibly enter contrasts of
% interests and let it run. The batch relies on PSOM (see Ref)
% see opt.mode for parallel computing on grid using qsub or msub
% <https://code.google.com/p/psom/wiki/ConfigurationPsom>
%
% FORMAT limo_batch
% limo_batch(option,model,contrast)
%
% INPUT if empty uses GUI
% option should be 'model specification' 'contrast only' or 'both'
% model is a structure that specifiy information to build a model
% model.set_files: a cell array of EEG.set (full path) for the different subjects
% model.cat_files: a cell array of categorial variable files
% model.cont_files: a cell array of continuous variable files
% model.defaults: specifiy the parameters to use for each subject
% model.defaults.analysis 'Time' 'Frequency' or 'Time-Frequency'
% model.defaults.fullfactorial 0/1
% model.defaults.zscore 0/1
% model.defaults.start starting time in ms
% model.defaults.end ending time in ms
% model.defaults.lowf starting point in Hz
% model.defaults.highf ending point in Hz
% model.defaults.bootstrap 0/1
% model.defaults.tfce 0/1
% model.defaults.channloc common channel locations (necessary if bootstrap = 1)
% contrast is a structure that specify which contrasts to run for which subject
% contrast.LIMO_files: a list of LIMO.mat (full path) for the different subjects
% this is optional if option 'both' is selected
% contrast.mat: a matrix of contrasts to run (assumes the same for all subjects)
%
% OUTPUT none - generate a directory per subject with GLM results in it
%
% see also limo_eeg limo_import_t limo_import_f limo_import_tf and psom in external folder
%
% Reference for pipeline engine
% Bellec P, Lavoie-Courchesne S, Dickinson P, Lerch JP, Zijdenbos AP and Evans AC (2012)
% The pipeline system for Octave and Matlab (PSOM): a lightweight scripting framework and
% execution engine for scientific workflows. Front. Neuroinform. 6:7.
% doi: 10.3389/fninf.2012.00007
% -----------------------------
% Copyright (C) LIMO Team 2014

% Cyril Pernet and Nicolas Chauveau 2012
% CP 24-06-2013 updated to be even more automatic + fix for new designs
% Cyril Pernet May 2014 - redesigned it using psom


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
    end
    % contrast
    if strcmp(option,'both')
        [FileName,PathName,FilterIndex]=uigetfile({'*.mat','MAT-files (*.mat)'; ...
            '*.txt','Text (*.txt)'}, 'Pick a matrix of contrasts');
        if strcmp(FileName(end-3:end),'.txt')
            contrast.mat = importdata(FileName);
        elseif strcmp(FileName(end-3:end),'.mat')
            FileName = load([PathName FileName]);
            contrast.mat = getfield(FileName,cell2mat(fieldnames(FileName)));
        else
            disp('limo batch aborded')
        end
        % update paths
        for f=1:size(model.set_files,1)
            [root,~,~] = fileparts(model.set_files{f});
            folder = ['GLM_' model.defaults.analysis];
            contrast.LIMO_files{f} = [root filesep folder filesep 'LIMO.mat'];
        end
    end
    
    if strcmp(option,'contrast only')
        
        % get paths
        [FileName,PathName,FilterIndex]=uigetfile({'*.txt','Text (*.txt)'; ...
            '*.mat','MAT-files (*.mat)'}, 'Pick a list of LIMO.mat files');
        if strcmp(FileName(end-3:end),'.txt')
            contrast.LIMO_files = importdata(FileName);
        elseif strcmp(FileName(end-3:end),'.mat')
            FileName = load([PathName FileName]);
            contrast.LIMO_files = getfield(FileName,cell2mat(fieldnames(FileName)));
        else
            disp('limo batch aborded')
        end
        
        % get the constrasts
        [FileName,PathName,FilterIndex]=uigetfile({'*.mat','MAT-files (*.mat)'; ...
            '*.txt','Text (*.txt)'}, 'Pick a matrix of contrasts');
        if strcmp(FileName(end-3:end),'.txt')
            contrast.mat = importdata(FileName);
        elseif strcmp(FileName(end-3:end),'.mat')
            FileName = load([PathName FileName]);
            contrast.mat = getfield(FileName,cell2mat(fieldnames(FileName)));
        else
            disp('limo batch aborded')
        end
    end
else
    option = varargin{1};
    if strcmp(option,'model specification') || strcmp(option,'both')
        model = varargin{2};
    end
    if strcmp(option,'contrast only') || strcmp(option,'both')
        contrast = varargin{3};
        if strcmp(option,'both') && ~isfield(contrast,'LIMO_files')
            for f=1:size(model.set_files,1)
                [root,~,~] = fileparts(model.set_files{f});
                folder = ['GLM_' model.defaults.analysis];
                contrast.LIMO_files = [root filesep folder filesep 'LIMO.mat'];
            end
        end
        if ~isfield(contrast,'mat')
            errordlg('the field contrast.mat is missing'); return
        end
    end
end

current =pwd;
mkdir('limo_batch_report')

%% -------------------------------------
%% build pipelines
%% -------------------------------------

if strcmp(option,'model specification') || strcmp(option,'both')
    % quick check
    if ~isempty(size(model.cat_files,1))
        if size(model.cat_files,1) ~= size(model.set_files,1)
            error('the number of set and cat files disagree')
        end
    end
    if ~isempty(size(model.cont_files,1))
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
        [root,~,~] = fileparts(model.set_files{subject});
        pipeline(subject).import.files_out = [root filesep 'GLM_' model.defaults.analysis filesep 'LIMO.mat'];
        pipeline(subject).import.opt.cat = model.cat_files{subject};
        pipeline(subject).import.opt.cont = model.cont_files{subject};
        pipeline(subject).import.opt.defaults = model.defaults;
        
        % make design and evaluate
        command = 'limo_batch_design_matrix(files_in)';
        pipeline(subject).design.command = command;
        pipeline(subject).design.files_in = pipeline(subject).import.files_out;
        pipeline(subject).design.files_out = [root filesep 'GLM_' model.defaults.analysis filesep 'Yr.mat'];
        
        % run GLM
        if strcmp(model.defaults.analysis,'Time') || strcmp(model.defaults.analysis,'Frequency');
            command = 'cd(fileparts(files_in)), limo_eeg(4)';
        else strcmp(model.defaults.analysis,'Time-Frequency');
            command = 'cd(fileparts(files_in)), limo_eeg_tf(4)';
        end
        pipeline(subject).glm.command = command;
        pipeline(subject).glm.files_in = pipeline(subject).import.files_out;
        pipeline(subject).glm.files_out = [root filesep 'GLM_' model.defaults.analysis filesep 'Betas.mat'];
    end
    
elseif strcmp(option,'contrast only') || strcmp(option,'both')
    
    for subject = 1:size(contrast.LIMO_files,1)
        command = 'limo_batch_contrast(files_in,opt.C)';
        pipeline(subject).contrast.command = command;
        pipeline(subject).contrast.files_in = contrast.LIMO_files{subject};
        pipeline(subject).contrast.opt.C = contrast.mat;
    end
end

%% -------------------------------------
%% run the analyses
%% -------------------------------------

% run pipelines and report
try
    N = size(model.set_files,1);
catch
    N = size(contrast.LIMO_files,1);
end

for subject = 1:N
    try
        opt.path_logs = [current filesep 'limo_batch_report' filesep 'subject' num2str(subject)];
        psom_run_pipeline(pipeline(subject),opt)
        report{subject} = ['subject ' num2str(subject) ' processed'];
    catch ME
        report{subject} = ['subject ' num2str(subject) ' failed'];
    end
end
cd([current filesep 'limo_batch_report'])
cell2csv('batch_report', report')
cd(current); disp('LIMO batch processing done');

