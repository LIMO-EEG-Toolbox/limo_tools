function LIMOPath = limo_random_select(stattest,expected_chanlocs,varargin)

% This function is used to combine parameters computed at the 1st level
% using LIMO_glm. Whereas in LIMO_glm observations are assumed independent
% (i.e. N-way ANOVA/ANCOVA or Regression), LIMO_random_effect distinguishes
% independents and non-independent (repeated) measures.
%
% No statistical test is done in LIMO_random_select, only the grouping
% and organization of the data - also creating the LIMO.mat structure.
% Once data are re-created, call is made to LIMO_random_robust.
%
% FORMAT: limo_random_select(stattest,expected_chanlocs)
%         limo_random_select(stattest,expected_chanlocs,options)
%
% INPUTS: stattest defines the statitiscal analysis 'one sample t-test'
%                                                   'two-samples t-test'
%                                                   'paired t-test'
%                                                   'regression'
%                                                   'N-Ways ANOVA'
%                                                   'ANCOVA'
%                                                   'Repeated Measures ANOVA'
%
%         expected_chanlocs: the EEGlab structure defining all channels
%                           (this file is created by STUDY or LIMO_tools)
%
%         Optional inputs (appears in LIMO.mat):
%                'LIMOfiles'  Cell array with the full paths to the subject file or file
%                             who contains the list of path to a group of sets. The
%                             dimensions of the cell corresponds to group (rows), factor
%                             and level respectively (columns).
%                'parameters' Cell array of parameters to be tested, relative to Beta files.
%                            e.g. {[1 3],[2 4]} or {[1 3],[2 4];[1 3],[2 4]} in case of 2 groups.
%                            use ones for con files, e.g. {[1 1],[1 1]}
%                            Add nested cells for more repetition levels.
%       --> for LIMOfiles and parameters the rule is groups in rows, repeated measures in columns
%                (at the expection of paired t-test where group applies ie use rows)
%                'regressor_file' a file or matrix of data to regress when stattest = 4
%                'analysis_type' is 'Full scalp analysis' or '1 channel/component only'
%                'channel' Index of the electrode(s) to use if '1 channel/component only'
%                            is selected in analysis_type
%                'type' is 'Channels' or 'Component'
%                'method': 'robust','weighted','mean'; 
%                'nboot' is the number of bootstrap to do (default = 1000)
%                'tfce' 0/1 indicates to computes tfce or not (default = 0)
%                'zscore' for regression design, default [] will ask user
%                         to zscore or not, set to 'Yes' or 'No'
%                'skip design check' for Regression, ANOVA, ANCOVA the
%                                    design will pop-up and user is asked
%                                    to start the analysis, set to 'Yes' to
%                                    skip this step
%
% OUTPUT filepath is the Path to the LIMO structure.
%        all files are created where the function is called (ie no directory as argument)
%
% Examples
% - repeated measure ANOVA with command line using Betas
% LIMOPath = limo_random_select('Repeated Measures ANOVA',chanlocs,'LIMOfiles',...
%     {'F:\WakemanHenson_Faces\eeg\derivatives\LIMO_Face_detection\Beta_files_FaceRepAll_GLM_Channels_Time_WLS.txt'},...
%     'analysis type','Full scalp analysis','parameters',{[1 2 3],[4 5 6],[7 8 9]},...
%     'factor names',{'face','repetition'},'type','Channels','nboot',0,'tfce',0);
% - t-test with command line using con files
%     for N=length(STUDY.subject):-1:1
%         data{1,N} = con1_files{N}(1);
%         data{2,N} = con2_files{N}(2); 
%     end
%     LIMOPath = limo_random_select('paired t-test',STUDY.limo.chanloc,...
%         'LIMOfiles',data,'analysis_type','Full scalp analysis', 'type','Channels','nboot',101,'tfce',1);
%
% Cyril Pernet, Ramon Martinez-Cancino, Arnaud Delorme
%
% ----------------------------------------
%  Copyright (C) LIMO Team 2020

LIMOPath = [];

%% check inputs

if nargin == 0
    help LIMO_random_select
    return
elseif nargin < 2
    error('LIMO_random_select: not enough arguments');
end

tests = {'one sample t-test', 'two-samples t-test', 'paired t-test', ...
    'regression', 'N-Ways ANOVA', 'ANCOVA', 'Repeated Measures ANOVA'};
if sum(strcmpi(stattest,tests)) == 0
    error('input argument error, stat test unknown')
end

try
    if ischar(expected_chanlocs)
        LIMO.data = load(expected_chanlocs);
        if isfield(LIMO.data,'expected_chanlocs')
            LIMO.data.chanlocs = LIMO.data.expected_chanlocs;
        end
        if isfield(LIMO.data,'channeighbstructmat')
            LIMO.data = renameStructField(LIMO.data, 'channeighbstructmat', 'neighbouring_matrix');
        end
    else
        LIMO.data.chanlocs            = expected_chanlocs.expected_chanlocs;
        LIMO.data.expected_chanlocs   = expected_chanlocs.expected_chanlocs;
        LIMO.data.neighbouring_matrix = expected_chanlocs.channeighbstructmat;
    end
catch chanloc_err
    error('expected_chanloc file data not recognized \n%s',chanloc_err.message)
end

% build LIMO.mat from inputs
LIMO.dir               = pwd;
LIMO.Level             = 2;
LIMO.Analysis          = [];
LIMO.Type              = [];
LIMO.data.data         = [];
LIMO.design.bootstrap  = 0;
LIMO.design.tfce       = 0;
LIMO.design.electrode  = [];
LIMO.design.component  = [];
LIMO.design.parameters = [];
LIMO.design.method     = [];
regressor_file         = [];
analysis_type          = [];
zopt                   = [];
skip_design_check      = 'No';
warning on

for in = 1:2:(nargin-2)
    if strcmpi(varargin{in},'LIMOfiles')
        if ~iscell(varargin{in+1})
            LIMO.data.data = {varargin{in+1}}; %#ok<CCAT1>
        else
            LIMO.data.data = varargin{in+1};
        end
    elseif strcmpi(varargin{in},'analysis type') || strcmpi(varargin{in},'analysis_type')
        if any(strcmpi(varargin{in+1},{'Full scalp analysis','1 channel/component only'}))
            analysis_type = varargin{in+1};
        else
            error('analysis type argument unrecognized')
        end
    elseif contains(varargin{in},'regressor','IgnoreCase',true)
        regressor_file = varargin{in+1};
    elseif strcmpi(varargin{in},'zscore')
        zopt = varargin{in+1};
    elseif strcmpi(varargin{in},'skip design check') || strcmpi(varargin{in},'skip_design_check')
        skip_design_check = varargin{in+1};
    elseif strcmpi(varargin{in},'channel')
        LIMO.design.electrode = varargin{in+1};
    elseif contains(varargin{in},'parameter','IgnoreCase',true)
        if iscell(varargin{in+1})
            LIMO.design.parameters = varargin{in+1};
        else
            LIMO.design.parameters = {varargin{in+1}};
        end
    elseif contains(varargin{in},'factor')
        LIMO.design.factor_names = varargin{in+1};
    elseif strcmpi(varargin{in},'method')
        if any(strcmpi(varargin{in+1},{'robust','weighted','mean'}))
            LIMO.design.method = varargin{in+1};
        else
            error('unrecognized method selected')
        end
    elseif strcmpi(varargin{in},'type')
        LIMO.Type = varargin{in+1};
    elseif strcmpi(varargin{in},'nboot')
        LIMO.design.bootstrap = varargin{in+1};
    elseif strcmpi(varargin{in},'tfce')
        LIMO.design.tfce = varargin{in+1};
    end
end

if isempty(analysis_type)
    analysis_type = limo_questdlg('Do you want to run a full analysis or a single channel/component analysis?','type of analysis?','1 channel/component only','Full scalp analysis','Full scalp analysis');
    if isempty(analysis_type)
        return
    end
end

if evalin( 'base', 'exist(''STUDY'',''var'') == 1' )
    global STUDY %#ok<*GVMIS,TLEV> 
end

% ----------------------------------
%%  One sample t-test and regression
% ----------------------------------
if strcmpi(stattest,'one sample t-test') || strcmpi(stattest,'regression')

    % get files
    % ---------
    if isempty(LIMO.data.data)
        [Names,Paths,LIMO.data.data] = limo_get_files;
    else
        if ischar(LIMO.data.data{1}) && length(LIMO.data.data) == 1 % Case for path to the files
            [Names,Paths,LIMO.data.data] = limo_get_files([],[],[],LIMO.data.data{1});
        else % Case when all paths are provided
            if size(LIMO.data.data,1) == 1
                LIMO.data.data = LIMO.data.data';
            end
            [Names,Paths,LIMO.data.data] = breaklimofiles(LIMO.data.data);
        end
    end
    LIMO.data.data_dir = Paths;

    if isempty(Names)
        limo_warndlg('no files selected')
        return
    end

    % check type of files and returns which beta param to test
    % -------------------------------------------------------
    if ~isfield(LIMO.design,'parameters') || isempty(LIMO.design.parameters)
        parameters = check_files(Paths, Names,1);
    else
        parameters = check_files(Paths, Names,1,LIMO.design.parameters{1});
    end

    if isempty(parameters)
        limo_errordlg('file selection failed or canceled, only Beta and Con files are supported','Selection error'); 
        return
    end

    % match frames, update LIMO
    % --------------------------
    [first_frame,last_frame,subj_chanlocs,~,LIMO] = match_frames(Paths,LIMO);

    % match channels, update LIMO
    % -----------------------------
    if strcmpi(stattest,'one sample t-test')
        LIMO = match_channels(1,analysis_type,LIMO);
    else
        LIMO = match_channels(4,analysis_type,LIMO);
    end

    % get data for all parameters
    % -----------------------------
    [data,removed] = getdata(1,analysis_type,first_frame,last_frame,subj_chanlocs,LIMO);
    if isempty(data)
        limo_errordlg('no data were retreived - check inputs and data files','limo_random_select');
        return
    end

    % if regression get regressor(s)
    % -------------------------------
    if strcmpi(stattest,'regression')
        if isempty(regressor_file)

            FilterIndex = 0;
            if ~isempty(STUDY)
                indvars = pop_listfactors(STUDY, 'gui', 'off', 'level', 'two', 'vartype', 'continuous');
                if ~isempty(indvars)
                    % get variable from study, DOES NOT HANDLE multiple sessions
                    uiList = { { 'style' 'text' 'string' 'Select subject specific variable(s) from the EEGLAB study' } ...
                               { 'style' 'listbox' 'string' { indvars.label } 'max' 2} ...
                               { 'style' 'text' 'string' 'These variables will be saved in the current folder as "regression_vars.txt"' } ...' ...
                               { 'style' 'text' 'string' 'Alternatively, press browse to load a text file with values to regress on'} };
                    res = inputgui('uilist', uiList, 'geometry', { [1] [1] [1] [1]}, 'geomvert', [1 3 1 1], 'cancel', 'Browse');
                    if isempty(res)
                        [FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','select regressor file');
                    else
                        FilterIndex = 1;
                        PathName = pwd;
                        FileName = 'regression_vars.txt';
                        std_saveindvar(STUDY, { indvars(res{1}).label }, fullfile(PathName, FileName));
                    end
                end
            end
            if FilterIndex == 0
                limo_warndlg( [ 'No EEGLAB STUDY present or not relevant variable found in the STUDY.' 10 ...
                    'Next a browsing window will ask you to select a text file to regress on' 10 ...
                    'with one row per subject and as many columns as continuous regressors.' 10 ...
                    'If you want to use continuous regressors and categorical variables ' 10 ...
                    '(age and gender for example), use ANCOVA instead.' ], 'LIMO regressor file');
                [FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','select regressor file');
            end
            if FilterIndex == 0
                return
            end
        elseif ischar(regressor_file) && exist(regressor_file,'file')
            [PathName,nametmp,exttmp] = fileparts(regressor_file);
            FileName = [nametmp exttmp]; clear nametmp exttmp;
        elseif isnumeric(regressor_file)
            X = regressor_file;
        else
            error('LIMO_random_select error: Provide a valid regressor file');
        end

        if ~exist('X','var')
            if strcmp(FileName(end-3:end),'.txt')
                X = load(fullfile(PathName,FileName));
            elseif strcmp(FileName(end-3:end),'.mat')
                X = load(fullfile(PathName,FileName));
                X = X.(cell2mat(fieldnames(X)));
            end
        end
        disp('Regressor(s) loaded');

        % check size and orientation
        N = size(data);
        N = N(end);
        if size(X,2) == N || size(X,2) == size(Paths,2)
            disp('X has been transposed'); X = X';
        end

        % adjust covariate(s) or data
        if size(X,1) > N 
            if sum(removed) ~=0
                try
                    index = 0;
                    for i=1:size(Paths,2)
                        index = index + 1;
                        if removed(i) == 1
                            X(index,:) = [];
                        end
                    end
                    warning('covariate adjusted for delete subjects');
                catch ME
                    limo_errordlg(sprintf('the number of regression value %g differs from the number of subjects %g',size(X,1),N),'Covariate error');
                    fprintf('%s',ME.message); return
                end
            end
            limo_errordlg(sprintf('the number of regression value %g differs from the number of subjects %g',size(X,1),N),'Covariate error');
        elseif sum(isnan(X(:))) ~= 0
            if ~skip_design_check
                if sum(sum(isnan(X),2)) == 1
                    limo_warndlg('loaded regressor(s) include NaN(s) - corresponding subjects are removed')
                end
            else
                fprintf(2, 'loaded regressor(s) include NaN(s) - corresponding subjects are removed\n');
            end
            
            sub_toremove = find(sum(isnan(X),2));
            X(sub_toremove,:) = [];
            if numel(size(data)) == 5 
                data(:,:,:,:,sub_toremove) = [];%<--- dim 4 = parameters
            else
                data(:,:,:,sub_toremove) = []; %<--- dim 3 = parameters
            end
        end
    end

    % load the right parameter and compute
    % ------------------------------------
    root = LIMO.dir;
    LIMOPath = cell(1,length(parameters));
    for i=1:length(parameters)
        if length(parameters) ~= 1 % make subfolders
            mkdir(fullfile(root,['parameter_' num2str(parameters(i))]));
            LIMO.dir = fullfile(root,['parameter_' num2str(parameters(i))]);
        end

        if strcmp(analysis_type,'1 channel/component only') && size(data,1) == 1
            if all(contains(Names,'con'))
                tmp_data = data;
            else
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    tmp               = squeeze(data(:,:,:,parameters(i),:));
                    tmp_data          = NaN(1,size(tmp,1),size(tmp,2),size(tmp,3)); % add dim 1 = 1 channel
                    tmp_data(1,:,:,:) = tmp; clear tmp;
                else
                    tmp               = squeeze(data(:,:,parameters(i),:));
                    tmp_data          = NaN(1,size(tmp,1),size(tmp,2));
                    tmp_data(1,:,:)   = tmp; clear tmp;
                end
            end
        else
            if strcmp(LIMO.Analysis,'Time-Frequency')
                tmp_data          = squeeze(data(:,:,:,parameters(i),:));
            else
                tmp_data          = squeeze(data(:,:,parameters(i),:));
            end
            if size(data,1) == 1
                tmp_data = reshape(tmp_data, [ 1 size(tmp_data)]);
            end
        end

        if strcmp(LIMO.Analysis,'Time-Frequency')
            LIMO.data.size3D = [size(tmp_data,1) size(tmp_data,2)*size(tmp_data,3) size(tmp_data,4)];
            LIMO.data.size4D = [size(tmp_data,1) size(tmp_data,2) size(tmp_data,3) size(tmp_data,4)];
        end

        % one-sample t-test
        % -----------------
        if strcmpi(stattest,'one sample t-test')
            Yr                 = tmp_data; clear tmp_data
            LIMO.design.name   = 'Robust one sample t-test';
            LIMO.design.X      = ones(size(data,4),1);
            if strcmpi(LIMO.design.method,'robust')
                LIMO.design.method = 'Trimmed mean';
            elseif strcmpi(LIMO.design.method,'weighted')
                LIMO.design.method = 'Weighted mean';
                [LIMO.design.weight.global,LIMO.design.weight.local] = ...
                    limo_group_outliers(Beta_files,LIMO.data.expected_chanlocs,LIMO.design.X);
            else
                LIMO.design.method = 'Mean';
            end
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
            save(fullfile(LIMO.dir,'Yr.mat'),'Yr','-v7.3');
            tmpname = limo_random_robust(1,fullfile(LIMO.dir,'Yr.mat'),...
                parameters(i),LIMO);
            if nargout ~= 0
                LIMOPath{i} = tmpname;
            end
            cd(root)

            % ------------
            % regression
            % -------------
        elseif strcmpi(stattest,'regression')
            if strcmp(LIMO.Analysis,'Time-Frequency')
                Yr = tmp_data(:,:,:,~isnan(sum(X,2)));
            else
                Yr = tmp_data(:,:,~isnan(sum(X,2)));
            end
            clear tmp_data
            LIMO.design.name = 'Regression';
            save(fullfile(LIMO.dir,'Yr.mat'),'Yr','-v7.3');
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
            if isempty(zopt)
                tmpname = limo_random_robust(4,Yr,X,...
                parameters(i),LIMO,'go',skip_design_check);
            else
                tmpname = limo_random_robust(4,Yr,X,...
                parameters(i),LIMO,'zscore',zopt,'go',skip_design_check);
            end

            if nargout ~= 0
                LIMOPath{i} = tmpname;
            end
        end
    end

    % ------------------------------
    %%  Two samples t-test
    % -----------------------------
elseif strcmpi(stattest,'two-samples t-test')

    LIMO.design.X = [];
    if isempty(LIMO.data.data)
        [Names{1},Paths{1},LIMO.data.data{1}] = limo_get_files(1);
        LIMO.data.data_dir{1}                 = Paths;
    else
        % input must be a 2*N cell array of cells
        [gp,list]=size(LIMO.data.data);
        if gp ==1 && list >=2
            LIMO.data.data = LIMO.data.data'; % fix user input
            [gp,list]=size(LIMO.data.data);
        end

        if gp>2
            msg = 'input must be a cell array of dimension 2 (gps) * N (list)';
            limo_error(sprintf('%s\n observed input is %g * %g',msg,size(LIMO.data.data)))
        end

        % now read
        if list == 1 && ischar(LIMO.data.data{1}) % Case for path to the files
            [Names{1},Paths{1},LIMO.data.data{1}] = limo_get_files([],[],[],LIMO.data.data{1});
            LIMO.data.data_dir{1}                 = Paths{1};
        else % Case when all paths are provided
            [Names{1},Paths{1}]   = breaklimofiles(LIMO.data.data(1,:));
            LIMO.data.data_dir{1} = Paths{1};
        end
    end

    if isempty(Names{1})
        limo_warndlg('Could not parse files names - function aborded')
        return
    end

    % check type of files and returns which beta param to test
    % -------------------------------------------------------
    if isfield(LIMO.design,'parameters')
        parameters =  cell2mat(LIMO.design.parameters);
    else
        LIMO.design.parameters = [];
        parameters             = [];
    end

    if isempty(parameters)
        parameters = check_files(Paths{1}, Names{1},1,[],'selectone');
        if isempty(parameters)
            return
        elseif length(parameters) > 1
            warning('only 1 parameter at a time for two-samples analysis, restricting to the 1st indicated')
            parameters = parameters(1);
        end
    else
        parameters(1) = check_files(Paths{1}, Names{1},1,parameters(1));
    end

    if length(LIMO.data.data) == 1
        [Names{2},Paths{2},LIMO.data.data{2}] = limo_get_files(2);
        LIMO.data.data_dir{2}                 = Paths;
    else
        if list == 1 && ischar(LIMO.data.data{2}) % Case for path to the files
            [Names{2},Paths{2},LIMO.data.data{2}] = limo_get_files([],[],[],LIMO.data.data{2});
            LIMO.data.data_dir{2}                 = Paths{2};
            LIMO.data.data_dir                    = LIMO.data.data_dir';
        else % Case when all paths are provided
            [Names{2},Paths{2}]   = breaklimofiles(LIMO.data.data(2,:));
            LIMO.data = rmfield(LIMO.data,'data');
            for sub = length(Paths{1}):-1:1; LIMO.data.data_dir{1}{sub} = Paths{1}{sub}; end
            for sub = length(Paths{1}):-1:1; LIMO.data.data{1}{sub} = fullfile(Paths{1}{sub},Names{1}{sub}); end
            for sub = length(Paths{2}):-1:1; LIMO.data.data_dir{2}{sub} = Paths{2}{sub}; end
            for sub = length(Paths{2}):-1:1; LIMO.data.data{2}{sub} = fullfile(Paths{2}{sub},Names{2}{sub}); end
            LIMO.data.data_dir = LIMO.data.data_dir';
            LIMO.data.data     = LIMO.data.data';
        end
    end

    if isempty(Names{2})
       limo_warndlg('Could not parse files names - function aborded')
       return
    end

    % check type of files and returns which beta param to test
    % -------------------------------------------------------
    if length(parameters) == 1
        parameters(2) = check_files(Paths{2}, Names{2},1,[],'selectone');
    else
        parameters(2) = check_files(Paths{2}, Names{2},1,parameters(2));
    end

    if ~isfield(LIMO.design,'parameters')
        LIMO.design.parameters = parameters;
    end

    % check type of files and returns which beta param to test
    % -------------------------------------------------------
    if ~all(contains(Names{1},'betas')) && ~all(contains(Names{1},'Betas')) && ~all(contains(Names{1},'con'))
        % do this only if betas - for con paramters = [1 1]
        for gp = 1:2
            for sub=1:size(LIMO.data.data_dir{gp},2)
                sub_LIMO = load(cell2mat(fullfile(LIMO.data.data_dir{gp}{1}(sub),'LIMO.mat')));
                if parameters(gp) > size(sub_LIMO.LIMO.design.X,2)-1
                    error('invalid parameter %g - design subject %s inconsistent',parameters(gp),cell2mat(fullfile(LIMO.data.data_dir{gp}{1}{sub})));
                end
            end
        end
    end

    % match frames, update LIMO
    [first_frame,last_frame,subj_chanlocs,~,LIMO] = match_frames(Paths,LIMO);

    % match channels, update LIMO
    LIMO = match_channels(2,analysis_type,LIMO);

    % get data for all parameters dim [channel, frame, param, nb, subjects]
    % -----------------------------------------------------------------
    data = getdata(2,analysis_type,first_frame,last_frame,subj_chanlocs,LIMO);

    % compute
    % --------
    % free some memory
    clear Names Paths channeighbstructmat expected_chanlocs subj_chanlocs

    if strcmp(LIMO.Analysis,'Time-Frequency')
        if strcmp(analysis_type,'1 channel/component only')
            tmp                = squeeze(data{1}(:,:,:,parameters(1),:));
            tmp_data1          = ones(1,size(tmp,1),size(tmp,2),size(tmp,3)); % add dim 1 = 1 channel
            tmp_data1(1,:,:,:) = tmp; clear tmp
            tmp                = squeeze(data{2}(:,:,:,parameters(2),:));
            tmp_data2          = ones(1,size(tmp,1),size(tmp,2),size(tmp,3));
            tmp_data2(1,:,:,:) = tmp; clear tmp
        else
            tmp_data1          = squeeze(data{1}(:,:,:,parameters(1),:));
            tmp_data2          = squeeze(data{2}(:,:,:,parameters(2),:));
        end

        if size(tmp_data1,1) ~= size(tmp_data2,1) || size(tmp_data1,2) ~= size(tmp_data2,2) || size(tmp_data1,3) ~= size(tmp_data2,3)
                limo_errordlg('file selection is corrupted, data sizes don''t match');
            return
        end

    else
        if strcmp(analysis_type,'1 channel/component only')
            tmp              = squeeze(data{1}(:,:,parameters(1),:));
            tmp_data1        = ones(1,size(tmp,1),size(tmp,2)); % add dim 1 = 1 channel
            tmp_data1(1,:,:) = tmp; clear tmp
            tmp              = squeeze(data{2}(:,:,parameters(2),:));
            tmp_data2        = ones(1,size(tmp,1),size(tmp,2));
            tmp_data2(1,:,:) = tmp; clear tmp
        else
            tmp_data1        = squeeze(data{1}(:,:,parameters(1),:));
            tmp_data2        = squeeze(data{2}(:,:,parameters(2),:));
        end

        if size(tmp_data1,1) ~= size(tmp_data2,1) || size(tmp_data1,2) ~= size(tmp_data2,2)
            limo_errordlg('file selection is corrupted, data sizes don''t match');
            return
        end
    end
    clear data

    if strcmp(LIMO.Analysis,'Time-Frequency')
        LIMO.data.size3D = [size(tmp_data1,1) size(tmp_data1,2)*size(tmp_data1,3) 5];
        LIMO.data.size4D = [size(tmp_data1,1) size(tmp_data1,2) size(tmp_data1,3) 5];
    end

    Y1r = tmp_data1; save Y1r Y1r, clear Y1r
    Y2r = tmp_data2; save Y2r Y2r, clear Y2r
    LIMO.design.method = 'Yuen t-test (Trimmed means)'; save LIMO LIMO
    tmpname = limo_random_robust(2,tmp_data1,tmp_data2,parameters,LIMO);
    if nargout ~= 0
        LIMOPath = tmpname;
    end
    if exist('data.mat','file')
        try delete data.mat; end %#ok<TRYNC>
    end

    % ------------------------------
    %%  Paired t-test
    % -----------------------------
elseif strcmpi(stattest,'paired t-test')

    LIMO.design.X = [];
    if isempty(LIMO.data.data)
        [Names{1},Paths{1},LIMO.data.data{1}] = limo_get_files;
        LIMO.data.data_dir                    = Paths{1};
    else
        % input must be a N*2 cell array of cells
        [list,pair]=size(LIMO.data.data);
        if list == 2 && pair >= 2
            LIMO.data.data = LIMO.data.data'; % fix user input
            [list,pair]=size(LIMO.data.data);
        end

        if pair>2
            msg = 'input must be a cell array of dimension N (list) *2 (pairs)';
            limo_error(sprintf('%s\n observed input is %g * %g',msg,size(LIMO.data.data)))
            return
        end

        % now read
        if list == 1 && ischar(LIMO.data.data{1}) % Case for path to the files
            [Names{1},Paths{1},LIMO.data.data{1}] = limo_get_files([],[],[],LIMO.data.data{1});
            LIMO.data.data_dir{1}                 = Paths{1};
        else % Case when all paths are provided
            [Names{1},Paths{1}]   = breaklimofiles(LIMO.data.data(:,1));
            LIMO.data.data_dir{1} = Paths{1};
        end
    end

    if isempty(Names{1})
        limo_warndlg('Could not parse files names - function aborded')
        return
    end

    % check type of files and returns which beta param to test
    % -------------------------------------------------------
    if isfield(LIMO.design,'parameters')
        parameters = cell2mat(LIMO.design.parameters);
    else
        LIMO.design.parameters = [];
        parameters             = [];
    end

    if isempty(parameters)
        parameters = check_files(Paths{1}, Names{1},1, [], 'selecttwo');
    else
        parameters = check_files(Paths{1}, Names{1},1,parameters);
    end

    if size(parameters,2) ~= 2
        if isempty(LIMO.design.parameters) && contains(Names{1}{1}, 'Betas')
            % hack only availbale if beta files NOT command line argument
            limo_errordlg('Error: for paired t-test, please select 2 parameters','Paired t-test');
            return
        end
    end

    if size(parameters,2) == 1 % either con file, or command line beta with one parameter
        % do 2nd pair of data
        if size(LIMO.data.data,2) == 1 % only 1st group of files
            [Names{2},Paths{2},LIMO.data.data{2}] = limo_get_files([],[],'select paired file');
            LIMO.data.data_dir{2} = Paths{2};
        else
            if list == 1 && ischar(LIMO.data.data{2})
                [Names{2},Paths{2},LIMO.data.data{2}] = limo_get_files([],[],[],LIMO.data.data{2});
                LIMO.data.data_dir{2} = Paths{2};
            else % Case when all paths are provided
                [Names{2},Paths{2}] = breaklimofiles(LIMO.data.data(:,2));
                LIMO.data = rmfield(LIMO.data,'data');
                for sub = length(Paths{1}):-1:1; LIMO.data.data_dir{1}{sub} = Paths{1}{sub}; end
                for sub = length(Paths{1}):-1:1; LIMO.data.data{1}{sub} = fullfile(Paths{1}{sub},Names{1}{sub}); end
                for sub = length(Paths{2}):-1:1; LIMO.data.data_dir{2}{sub} = Paths{2}{sub}; end
                for sub = length(Paths{2}):-1:1; LIMO.data.data{2}{sub} = fullfile(Paths{2}{sub},Names{2}{sub}); end
            end
        end

        if isempty(Names{2})
            limo_warndlg('Could not parse files names - function aborded')
            return
        else
            if ~isempty(LIMO.design.parameters)
                % hack only availbale if beta files and command line argument // not allowed otherwise because it's a paired design
                parameters(2) = check_files(Paths{2}, Names{2},1,parameters(2));
            else
                newparameters = check_files(Paths{2}, Names{2},1);
                if newparameters ~= 1
                    error('paired t-test second set must also be con files')
                else
                    parameters = 1;
                end
            end
            con_parameters    = [str2double(unique(cellfun(@(x) x(5:end-4),Names{1}))) ...
                str2double(unique(cellfun(@(x) x(5:end-4),Names{2})))];
            if all(isnan(con_parameters))
                clear con_parameters % was betas from command line
            end
        end

        if size(Names{1},2) ~= size(Names{2},2)
            limo_errordlg('the nb of files differs between pairs 1 and 2','Paired t-test error'); 
            return
        end
    elseif size(parameters,2) ~=2 % if it was beta file one needs a pair of parameters
            limo_errordlg('2 parameters must be selected for beta files','Paired t-test error'); 
            return
    else % check Betas match design
        for s = 1:length(Paths{1})
            sub_LIMO = load(fullfile(Paths{1}{s},'LIMO.mat'));
            if max(parameters) > size(sub_LIMO.LIMO.design.X,2)
                limo_errordlg('invalid parameter(s)','Paired t-test error'); 
                return
            end
        end
        cd(LIMO.dir);
    end

    % match frames, update LIMO
    [first_frame,last_frame,subj_chanlocs,~,LIMO] = match_frames(Paths,LIMO);

    % match channels, update LIMO
    LIMO = match_channels(3,analysis_type,LIMO);

    % get data for all parameters dim [channel, frame, param, nb, subjects]
    % -----------------------------------------------------------------
    if any(size(LIMO.data.data) == 2)
        data = getdata(2,analysis_type,first_frame,last_frame,subj_chanlocs,LIMO);
    else
        data = getdata(1,analysis_type,first_frame,last_frame,subj_chanlocs,LIMO);
    end

    % compute
    % --------
    if strcmpi(LIMO.Analysis,'Time-Frequency')
        if strcmpi(analysis_type,'1 channel/component only')
            if size(parameters,2) == 2 % beta files
                tmp                = squeeze(data(:,:,:,parameters(1),:));
                tmp_data1          = ones(1,size(tmp,1),size(tmp,2),size(tmp,3)); % add dim 1 = 1 channel
                tmp_data1(1,:,:,:) = tmp; clear tmp
                tmp                = squeeze(data(:,:,:,parameters(2),:));
                tmp_data2          = ones(1,size(tmp,1),size(tmp,2),size(tmp,3));
                tmp_data2(1,:,:,:) = tmp; clear tmp
            else % con files
                tmp                = squeeze(data{1}(:,:,:,:,:));
                tmp_data1          = ones(1,size(tmp,1),size(tmp,2),size(tmp,3)); % add dim 1 = 1 channel
                tmp_data1(1,:,:,:) = tmp; clear tmp
                tmp                = squeeze(data{2}(:,:,:,:,:));
                tmp_data2          = ones(1,size(tmp,1),size(tmp,2),size(tmp,3));
                tmp_data2(1,:,:,:) = tmp; clear tmp
            end
        else
            if size(parameters,2) == 2 % beta files
                if iscell(data)
                    tmp_data1 = squeeze(data{1}(:,:,:,parameters(1),:));
                    tmp_data2 = squeeze(data{2}(:,:,:,parameters(2),:));
                else
                    tmp_data1 = squeeze(data(:,:,:,parameters(1),:));
                    tmp_data2 = squeeze(data(:,:,:,parameters(2),:));
                end
            else % con files
                tmp_data1 = squeeze(data{1}(:,:,:,:));
                tmp_data2 = squeeze(data{2}(:,:,:,:));
            end
        end
    else
        if strcmpi(analysis_type,'1 channel/component only')
            if size(parameters,2) == 2 % beta files
                tmp              = squeeze(data(:,:,parameters(1),:));
                tmp_data1        = ones(1,size(tmp,1),size(tmp,2)); % add dim 1 = 1 channel
                tmp_data1(1,:,:) = tmp; clear tmp
                tmp              = squeeze(data(:,:,parameters(2),:));
                tmp_data2        = ones(1,size(tmp,1),size(tmp,2));
                tmp_data2(1,:,:) = tmp; clear tmp
            else % con files
                tmp              = squeeze(data{1}(:,:,:,:));
                tmp_data1        = ones(1,size(tmp,1),size(tmp,2)); % add dim 1 = 1 channel
                tmp_data1(1,:,:) = tmp; clear tmp
                tmp              = squeeze(data{2}(:,:,:,:));
                tmp_data2        = ones(1,size(tmp,1),size(tmp,2));
                tmp_data2(1,:,:) = tmp; clear tmp
            end
        else
            if size(parameters,2) == 2 % beta files
                if iscell(data)
                    tmp_data1 = squeeze(data{1}(:,:,parameters(1),:));
                    tmp_data2 = squeeze(data{2}(:,:,parameters(2),:));
                else
                    tmp_data1 = squeeze(data(:,:,parameters(1),:));
                    tmp_data2 = squeeze(data(:,:,parameters(2),:));
                end
            else % con files
                tmp_data1 = squeeze(data{1}(:,:,:));
                tmp_data2 = squeeze(data{2}(:,:,:));
            end
        end
    end

    if strcmpi(LIMO.design.method,'Robust')
        LIMO.design.method = 'Yuen t-test (Trimmed means)';
    elseif strcmpi(LIMO.design.method,'weighted')
        LIMO.design.method = 'Weighted mean';
        [LIMO.design.weight.global,LIMO.design.weight.local] = ...
            limo_group_outliers(Beta_files,LIMO.data.expected_chanlocs,LIMO.design.X);
    else
        LIMO.design.method = 'Mean';
    end

    if strcmp(LIMO.Analysis,'Time-Frequency')
        LIMO.data.size3D = [size(tmp_data1,1) size(tmp_data1,2)*size(tmp_data1,3) 5];
        LIMO.data.size4D = [size(tmp_data1,1) size(tmp_data1,2) size(tmp_data1,3) 5];
    end
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')

    Y1r = tmp_data1; save Y1r Y1r, clear Y1r
    Y2r = tmp_data2; save Y2r Y2r, clear Y2r
    if exist('con_parameters','var')
        parameters = con_parameters; % substitute param of data by con values for consistent naming
    end
    LIMO.design.parameters = parameters; % update in any cases
    tmpname = limo_random_robust(3,tmp_data1,tmp_data2,parameters,LIMO);
    if nargout ~= 0
        LIMOPath = tmpname;
    end

    % -----------------------------------
    %%  Various sorts of ANOVAs/ANCOVAs
    % -----------------------------------
elseif strcmpi(stattest,'N-Ways ANOVA') || strcmpi(stattest,'ANCOVA')

    % ---------------------------------------------------------------------
    %              One Way ANOVA / ANCOVA
    % ---------------------------------------------------------------------

    % Ask for Gp
    % -------------
    if ~isempty(LIMO.data.data)
        [gp_nb,maxsub] = size(LIMO.data.data);
        % to compute we must have more data than groups
        if gp_nb == 1 || gp_nb > maxsub % groups are always in row
            LIMO.data.data = LIMO.data.data';
            [gp_nb,maxsub] = size(LIMO.data.data);
            limo_warndlg('same number of groups and files (%g) in - be sure groups are in rows',gp_nb)
        elseif gp_nb == maxsub
            warning('input transposed %g groups',gp_nb)
        elseif gp_nb == 1 && maxsub == 1
            limo_errordlg('only one group detected - wrong input or test'); 
            return
        end
    else
        gp_nb = cell2mat(limo_inputdlg('How many independent groups? e.g. 3 or [3 2] for nested gps','Groups'));
        if isempty(gp_nb)
            return
        elseif sum(str2double(gp_nb) <= 2) && strcmpi(stattest,'N-Ways ANOVA')
            limo_errordlg('at least 3 groups are expected for a N-ways ANOVA')
            return
        elseif sum(str2double(gp_nb) <= 1) && strcmpi(stattest,'ANCOVA')
            limo_errordlg('at least 2 groups are expected for an ANCOVA')
            return
        else
            gp_nb          = str2double(gp_nb);
            Names          = cell(gp_nb,1);
            Paths          = cell(gp_nb,1);
            LIMO.data.data = cell(gp_nb,1);
        end
    end

    % select data per gp / conditions
    % ---------------------------------
    for i=1:gp_nb
        if isempty(LIMO.data.data{i})
            [Names{i},Paths{i},LIMO.data.data{i}] = limo_get_files([' beta or con file gp ',num2str(i)]);
        else
            if maxsub == 1 && ischar(LIMO.data.data{i}) % Case for path to the files
                [Names{i},Paths{i},LIMO.data.data{i}] = limo_get_files([],[],[],LIMO.data.data{i});
            else % Case when all paths are provided
                [Names{i},Paths{i},LIMO.data.data{i}] = breaklimofiles(LIMO.data.data(i,:));
                if i == gp_nb
                    LIMO.data.data = LIMO.data.data(:,1); % remove trailing names
                end
            end
        end
        if isempty(Names{i}); return; end
    end

    if isempty(LIMO.design.parameters)
        % if all con - no need to ask
        if ~isempty(STUDY) && ~any(contains(Names{1}, 'con_'))
            [param,~] = get_beta_indices('selectmulti', fullfile(Paths{1}{1}, Names{1}{1}));
        else
            param = check_files(Paths, Names,gp_nb);
            if param ~= 1
                param = cell2mat(limo_inputdlg('which parameters to test e.g [2]','parameters option'));
            else
                param = num2str(param); % pretend it's limo_inputdlg output
            end
        end

        if isempty(param)
            disp('selection aborded'); return
        else
            if ischar(param)
                parameters = eval(['[' param ']']);
            else
                parameters = param;
            end
        end

        if contains(stattest,'ANCOVA','IgnoreCase',true) && ...
            length(parameters) > 1 % this is only group * cov without repeated measures
            limo_errordlg(sprintf('The ANCOVA model doesn''t deal with repeated measures,\n you could use contrasts per subject to create an effect to covary on'))
            return
        end

        if length(param) == 1
            parameters = repmat(param,1,gp_nb);
        end
    else
        if iscell(LIMO.design.parameters)
            LIMO.design.parameters = cell2mat(LIMO.design.parameters);
        end
        if length(LIMO.design.parameters) == 1
            parameters = repmat(LIMO.design.parameters,1,gp_nb);
            LIMO.design.parameters = parameters;
        else
            parameters = LIMO.design.parameters;
        end
    end

    if length(parameters)>gp_nb
        limo_errordlg('%g parameter selected, longer than the number of groups %g',length(parameters),gp_nb)
        return
    else
        parameters = check_files(Paths, Names,gp_nb,parameters);
    end

    if isempty(parameters)
        limo_errordlg('selected files are not Beta or con files');
        return
    end
    LIMO.data.data_dir     = Paths';
    LIMO.design.parameters = parameters;

    % organize the data in such a way that we can easily compute stuff
    % ---------------------------------------------------------------------
    % match frames, update LIMO
    [first_frame,last_frame,subj_chanlocs,~,LIMO] = match_frames(limo_concatcells(Paths,'char'),LIMO);

    % match channels, update LIMO
    LIMO = match_channels(5,analysis_type,LIMO);

    % get data for all parameters dim [channel, frame, param, nb, subjects]
    % -----------------------------------------------------------------
    [data,removed] = getdata(2,analysis_type,first_frame,last_frame,subj_chanlocs,LIMO);
    nb_subjects    = cellfun(@(x) size(x,ndims(x)), data);

    % get some nice comment in LIMO.mat
    if strcmpi(analysis_type,'Full scalp analysis')
        if strcmpi(LIMO.Type,'Components')
            if contains(stattest,'N-Ways','IgnoreCase',true)
                LIMO.design.name = 'N-ways ANOVA all components';
            else
                LIMO.design.name = 'ANCOVA all components';
            end
        else
            if contains(stattest,'N-Ways','IgnoreCase',true)
                LIMO.design.name = 'N-ways ANOVA all channels';
            else
                LIMO.design.name = 'ANCOVA all channels';
            end
        end

    elseif strcmpi(analysis_type,'1 channel/component only')
        if strcmpi(LIMO.Type,'Components')
            if contains(stattest,'N-Ways','IgnoreCase',true)
                LIMO.design.name = 'N-ways ANOVA one component';
            else
                LIMO.design.name = 'ANCOVA one component';
            end
        else
            if contains(stattest,'N-Ways','IgnoreCase',true)
                LIMO.design.name = 'N-ways ANOVA one channel';
            else
                LIMO.design.name = 'ANCOVA one channel';
            end
        end
    end

    % now load covariates and check it matches data
    if strcmpi(stattest,'ANCOVA')
        if isempty(regressor_file)
            if ~isempty(STUDY)
                indvars = pop_listfactors(STUDY, 'gui', 'off', 'level', 'two', 'vartype', 'continuous');
                if ~isempty(indvars)
                    % get variable from study, DOES NOT HANDLE multiple sessions
                    uiList = { { 'style' 'text' 'string' 'Select subject specific variable(s) from the EEGLAB study' } ...
                               { 'style' 'listbox' 'string' { indvars.label } 'max' 2} ...
                               { 'style' 'checkbox' 'string' 'Compute interaction with sessions (will show as an additional covariate)' } ...
                               { 'style' 'checkbox' 'string' 'Compute interaction with groups (will show as an additional covariate)' } ...
                               { 'style' 'text' 'string' '' } ...
                               { 'style' 'text' 'string' 'These variables will be saved in the current folder as "covariate_vars.txt"' } ...' ...
                               { 'style' 'text' 'string' 'Alternatively, press browse to load a text file with values to regress on'} };
                    if length(STUDY.group)   < 2, uiList(4) = []; end
                    if length(STUDY.session) < 2, uiList(3) = []; end
                    res = inputgui('uilist', uiList, 'geometry', mattocell(ones(1,length(uiList))), 'geomvert', [1 3 ones(1,length(uiList)-2)], 'cancel', 'Browse');
                    if isempty(res)
                        [FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','select regressor file');
                    else
                        interaction_sess  = false;
                        interaction_group = false;
                        FilterIndex = 1;
                        PathName = pwd;
                        FileName = 'covariate_vars.txt';
                        if length(res) > 1
                            if length(STUDY.session) >=2, interaction_sess  = res{2}; end
                            if length(STUDY.group)   >=2, interaction_group = res{end}; end
                        end
                        std_saveindvar(STUDY, { indvars(res{1}).label }, fullfile(PathName, FileName), STUDY.session, STUDY.group, interaction_sess, interaction_group); % note that here it should depend on what the user select (group/session)
                    end
                end
            else
                [FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','select covariate file');
            end
            if FilterIndex == 0
                return
            else
                regressor_file = fullfile(PathName,FileName);
            end
        end

        if ischar(regressor_file)
            X = load(regressor_file);
            if ~isnumeric(X)
                X = X.(cell2mat(fieldnames(X)));
            end
        else
            X = regressor_file;
        end
        disp('Regressor(s) loaded');

        try
            if size(X,2) == sum(nb_subjects)
                disp('regressor transposed');
                X = X';
            end

            if size(X,1) ~= sum(nb_subjects)
                try
                    index = 0;
                    for h=1:gp_nb
                        for i=1:size(Paths{h},2)
                            index = index + 1;
                            if removed{h}(i) == 1
                                X(index,:) = [];
                            end
                        end
                    end
                    disp('covariate adjusted for delete subjects');
                catch ME
                    if size(X,1) ~= sum(nb_subjects)
                        if exist('errordlg2','file')
                            errordlg2(sprintf('the number of regression value %g differs from the number of subjects %g',size(X,1),sum(nb_subjects))); return
                        else
                            errordlg(sprintf('the number of regression value %g differs from the number of subjects %g',size(X,1),sum(nb_subjects))); return
                        end
                    else
                        if exist('errordlg2','file')
                            errordlg2(sprintf('log error:%s',ME.message),'fail adjusting covarate(s)'); return
                        else
                            errordlg(sprintf('log error:%s',ME.message),'fail adjusting covarate(s)'); return
                        end
                    end
                end
            end
        catch ME
            if exist('errordlg2','file')
                errordlg2('data format or length of the covariate does not match');
            else
                errordlg('data format or length of the covariate does not match');
            end
            error('%s',ME.message)
        end
        Cont = X;
        clear X
    else
        Cont = [];
    end

    % 4th send relevant info to LIMO_random_robust
    % ----------------------------------------------
    % for N ways ANOVA we pass the data and X, so we simply stack
    % the data on top of each other per gp - this allows using betas
    % and still compare possibly different parameters
    index = 1;
    for i=1:gp_nb
        current_param  = parameters(i); % select only relevant parameters (usually the same)
        if strcmp(LIMO.Analysis,'Time-Frequency')
            if i==1
                sz = [ size(data{i}) 1 1 1];
                tmp_data = NaN(sz([1 2 3 5]));
            end
            tmp_data(:,:,:,index:(sum(nb_subjects(1:i))))  = squeeze(data{i}(:,:,:,current_param,:));
        else
            if i==1
                sz = [ size(data{i}) 1 1 1];
                tmp_data = NaN(sz([1 2 4]));
            end
            tmp_data(:,:,index:(sum(nb_subjects(1:i)))) = squeeze(data{i}(:,:,current_param,:));
        end
        index = index + nb_subjects(i);
    end

    % make categorical variable
    % ------------------------
    index = 1;
    Cat = zeros(sum(nb_subjects),1);
    for n=1:length(nb_subjects)
        Cat(index:nb_subjects(n)+index-1) = n;
        index = index + nb_subjects(n);
    end

    if strcmp(LIMO.Analysis,'Time-Frequency')
        LIMO.data.size3D = [size(tmp_data,1) size(tmp_data,2)*size(tmp_data,3) size(tmp_data,4)];
        LIMO.data.size4D = [size(tmp_data,1) size(tmp_data,2) size(tmp_data,3) size(tmp_data,4)];
    end
    save LIMO LIMO

    % clear some memory
    clear Names Paths data channeighbstructmat expected_chanlocs subj_chanlocs

    % do the analysis
    Yr = tmp_data; clear tmp_data; save Yr Yr
    tmpname = limo_random_robust(5,Yr,Cat,Cont,LIMO,'go',skip_design_check);
    if nargout ~= 0
        LIMOPath = tmpname;
    end

elseif strcmpi(stattest,'Repeated measures ANOVA')
    % ---------------------------------------------------------------------
    %              Repeated measures ANOVA
    % ---------------------------------------------------------------------

    % get some comment in LIMO.mat
    if strcmp(analysis_type,'Full scalp analysis')
        if strcmp(LIMO.Type,'Components')
            LIMO.design.name = 'Repeated measures ANOVA all components';
        else
            LIMO.design.name = 'Repeated measures ANOVA all channels';
        end
    elseif strcmp(analysis_type,'1 channel/component only')
        if strcmp(LIMO.Type,'Components')
            LIMO.design.name = 'Repeated measures ANOVA one component';
        else
            LIMO.design.name = 'Repeated measures ANOVA one channel';
        end
    end

    % Ask for Gp
    % ----------
    gp_nb = [];
    if ~isempty(LIMO.data.data)
        gp_nb = size(LIMO.data.data,1);
    else
        if exist('STUDY','var')
            gp_val = length(STUDY.group);
        else
            gp_val = 1;
        end
        gp_nb = cell2mat(limo_inputdlg('How many independent groups of subjects or session per subject?','Groups', 1, {gp_val}));
    end

    if isempty(gp_nb)
        return
    elseif length(gp_nb) > 1
        limo_errordlg('only 1 independent factor (with n groups) is handled by LIMO')
        return
    else
        if ischar(gp_nb); gp_nb = str2double(gp_nb); end
        if gp_nb == 0; gp_nb = 1; end
    end

    % Ask for Repeated Measures
    % --------------------------
    if ~isempty(LIMO.design.parameters) % infer factors from parameters
        if size(LIMO.design.parameters,1) == 1 && gp_nb > 1
            LIMO.design.parameters = repmat(LIMO.design.parameters,gp_nb,1);
        end

        for g=gp_nb:-1:1 % strings like GUI
            if all(size(LIMO.design.parameters(g,:))==1)
                factor_nb{g} = num2str(length(LIMO.design.parameters{g}));
            else
                factor_nb{g} = num2str(getlevels(LIMO.design.parameters(g,:)));
            end
        end

        % check that groups have the same number of factors
        if ~cellfun(@(x) length(x)==length(LIMO.design.parameters{1}),LIMO.design.parameters)
            error('parameters input sizes different between groups, the number of factors to infer must be identical')
        else
            factor_nb = factor_nb{1};
        end
    else
        uiList = { { 'style' 'text' 'string' 'Enter repeated factors level' 'fontweight' 'bold'} ...
                     { 'style' 'text' 'string' '(for more than 3 factors use command line, see help limo_random_select)' } ...
                 { 'style' 'text' 'string' '' } ...
                   { 'style' 'text' 'string' 'Name' } ...
                   { 'style' 'text' 'string' 'Number of measures' } ...
                   { 'style' 'text' 'string' 'Factor 1' } ...
                   { 'style' 'edit' 'string' 'Factor 1' } ...
                   {} { 'style' 'edit' 'string' '3' } {} ...
                   { 'style' 'text' 'string' 'Factor 2 (if any)' } ...
                   { 'style' 'edit' 'string' '' } ...
                   {} { 'style' 'edit' 'string' '' } {}...
                   { 'style' 'text' 'string' 'Factor 3 (if any)' } ...
                   { 'style' 'edit' 'string' '' } ...
                   {} { 'style' 'edit' 'string' '' } {} };
        uiGeom = { [1] [1] [1 1 1] [1 1 0.3 0.3 0.3] [1 1 0.3 0.3 0.3] [1 1 0.3 0.3 0.3] };
        res = inputgui('uilist', uiList, 'geometry', uiGeom);
        if isempty(res)
            return;
        end
        
        factor_nb = [ res{2} ' ' res{4} ' ' res{6} ];
        factor_names = { res{1} res{3} res{5} };
    end

    % Cases of wrong input
    try
        factor_nb = eval(['[' factor_nb ']']);
        if ~exist('factor_names','var')
            factor_names = LIMO.design.factor_names;
        end
        factor_names = factor_names(1:length(factor_nb));
    catch ME
        limo_errordlg(sprintf('log error: %s',ME.message),'could not evaluate factors')
        return
    end

    if isempty(factor_nb) || length(factor_nb)==1 && factor_nb == 0
        limo_warndlg('no factor entered, Rep. ANOVA aborded');
    end

    % 2nd select data per gp / conditions
    % ---------------------------------------------------------------
    if ~isempty(LIMO.data.data)
        for gp = gp_nb:-1:1
            isbeta(gp) = mean(cellfun(@(x) contains(x,'Beta'),LIMO.data.data(gp,:)));
            iscon(gp)  = mean(cellfun(@(x) contains(x,'con'),LIMO.data.data(gp,:)));
        end

        if any(isbeta(:)) && any(iscon(:))
            error('input data mix Beta and con files - not supported')
        elseif sum(isbeta(:)) == 0 && sum(iscon(:)) == 0 % maybe it's a custom file name
            for gp = gp_nb:-1:1
                all_files  = limo_get_files([],[],[],LIMO.data.data{gp});
                isbeta(gp) = mean(cellfun(@(x) contains(x,'Beta'),all_files));
                iscon(gp)  = mean(cellfun(@(x) contains(x,'con'),all_files));
            end
            if any(isbeta) && any(iscon)
                error('input data mix Beta and con files - not supported')
            elseif sum(isbeta) == 0 && sum(iscon) == 0
                error('not all betas or all con?? issue reading files')
            else
                if all(isbeta); a = 'beta';
                elseif all(iscon); a = 'con'; end
            end
        else
            if all(isbeta); a = 'beta';
            elseif all(iscon); a = 'con'; end
        end
    else
        a = limo_questdlg('Do you want to load several contrast files per subject or a single beta file per group?','ANOVA loading files','con','beta','beta');
    end

    % beta files
    if strcmp(a,'beta')
        for i=1:gp_nb
            if length(LIMO.data.data) < i
                [Names{i},Paths{i},LIMO.data.data{i}] = limo_get_files([' beta file gp ',num2str(i)]);
            else
                if ischar(LIMO.data.data{i}) && size(LIMO.data.data,2) == 1 % Case for path to the files
                    [Names{i},Paths{i},LIMO.data.data{i}] = limo_get_files([],[],[],LIMO.data.data{i});
                else % Case when all paths are provided
                    [Names{i},Paths{i},LIMO.data.data{i}] = breaklimofiles(LIMO.data.data{i});
                end
            end

            if isempty(Names{i})
                warndlg2('files names not loaded'); return
            end
        
            if isfield(LIMO.design,'parameters') && ~isempty(LIMO.design.parameters)
                if length(factor_nb) <=2
                    parameters(i,:) = check_files(Paths{i}, Names{i},1,cell2mat(LIMO.design.parameters(i,:)));
                else % cell of cells
                    all_param = LIMO.design.parameters;
                    while any(cellfun(@iscell,all_param))
                        all_param = [all_param{cellfun(@iscell,all_param)} all_param(~cellfun(@iscell,all_param))];
                    end
                    parameters(i,:) = check_files(Paths{i}, Names{i},1,cell2mat(all_param));
                end
            else % GUI
                if i > 1
                    parameters(i,:) = parameters(1,:); % copy design for other groups
                else
                    [parameters(i,:),betas] = check_files(Paths{i}, Names{i},1,[], '');
                    if length(factor_nb) > 1
                        paramTmp  = reorganize_params(parameters(i,:), betas,factor_names, factor_nb);
                        if ~isempty(paramTmp)
                            parameters(i,:) = paramTmp;
                        else
                            return
                        end
                    end
                end
            end
        end

        if size(parameters,2) ~= prod(factor_nb)
            limo_warndlg(['the number of beta parameter chosen (',num2str(size(parameters,2)), ...
                ') does not match the total number of factor levels (',num2str(prod(factor_nb)),')'])
            return
        end

        % check if variables names are accessible
        if exist(fullfile(Paths{1}{1}, 'LIMO.mat'),'file')
            LIMOtmp = load('-mat', fullfile(Paths{1}{1}, 'LIMO.mat'));
            if isfield(LIMOtmp.LIMO.design, 'labels')
                paramLinear = parameters(i,:);
                if iscell(paramLinear) paramLinear = [ paramLinear{:} ]; end
                LIMO.design.labels = LIMOtmp.LIMO.design.labels(paramLinear);
            end
        end

        if size(LIMO.data.data,2) == gp_nb
            LIMO.data.data = LIMO.data.data'; % gps in rows
        end
        LIMO.data.data_dir = Paths';

    else  % multiple con files
        for i=1:gp_nb
            if isempty(LIMO.data.data) % GUI
                for j=1:length(factor_nb)
                    for k=1:factor_nb(j)
                        [Names{i,k},Paths{i,k},LIMO.data.data{i,k}] = limo_get_files([' gp ',num2str(i),' factor ',num2str(j),' level ',num2str(k)]);

                        if isempty(Names{i,k})
                            warning('no files found - selection aborded');
                            return
                        else
                            LIMO.data.data_dir{i,k} = Paths{i,k};
                        end
                    end
                end
            else
                for j=1:length(LIMO.data.data(i,:))
                    % Case for path to the files
                    if all(size(LIMO.data.data(i,j)))
                        [Names{i,j},Paths{i,j},LIMO.data.data{i,j}] = limo_get_files([],[],[],LIMO.data.data{i,j});
                        % Case when all paths are provided
                    elseif size(LIMO.data.data{num2str(i),num2str(j)},1) > 1
                        [Names{i,j},Paths{i,j},LIMO.data.data{i,j}] = breaklimofiles(LIMO.data.data{i}{i,j});
                    end

                    if isempty(Names{i,j})
                        limo_warndlg('no files found - selection aborded');
                        return;
                    else
                        LIMO.data.data_dir{i,j} = Paths{i,j};
                    end
                end
            end
            gp_param(i) = check_files(Paths(i,:), Names(i,:),j);
        end

        if all(gp_param) % comfirms it's all con files
            parameters = repmat(1:length(LIMO.data.data),gp_nb,1);
        else
            limo_errordlg('Inconsistent Names between groups - all con files expected of 1 parameter to extract');
            return
        end
    end
    LIMO.design.parameters = parameters; % move from factor to which param in file is picked

    % 3rd organize data
    % ---------------------------------------------------------------------
    % match frames, update LIMO
    [first_frame,last_frame,subj_chanlocs,~,LIMO] = match_frames(Paths(:)',LIMO);

    % match channels, update LIMO
    LIMO = match_channels(6,analysis_type,LIMO);

    % get data for all parameters
    % -----------------------------
    subject_index = 1;
    matrix_index  = 1;
    for h = 1:gp_nb % each group
        nb_subjects(h) = 0;
        for i=1:size(Paths{h},2)
            if all(contains(LIMO.data.data{h},'Betas')) % set of beta files
                tmp = load(cell2mat(LIMO.data.data{h}(i)));
                tmp = tmp.(cell2mat(fieldnames(tmp)));
            else % iterate through con files making tmp like a beta files
                for c = size(LIMO.data.data,2):-1:1
                    con = load(cell2mat(LIMO.data.data{h,c}(i)));
                    con = con.(cell2mat(fieldnames(con)));
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                         tmp(:,:,:,c) = con(:,:,:,1);
                    else
                        tmp(:,:,c) = con(:,:,1);
                    end
                end
            end

            % get indices to trim data
            if strcmp(LIMO.Analysis,'Time-Frequency')
                begins_at = fliplr((max(first_frame) - first_frame(subject_index,:) + 1)); % returns time/freq/or freq-time
                ends_at(1) = size(tmp,2) - (last_frame(subject_index,2) - min(last_frame(:,2)));
                ends_at(2) = size(tmp,3) - (last_frame(subject_index,1) - min(last_frame(:,1)));
            else
                begins_at = max(first_frame) - first_frame(subject_index) + 1;
                ends_at = size(tmp,2) - (last_frame(subject_index) - min(last_frame));
            end

            % data are of dim size(expected_chanlocs,2), latter start/earlier stop across subjects, parameters, nb of subjects
            if strcmp(analysis_type,'Full scalp analysis') %&& size(subj_chanlocs(subject_index).chanlocs,2) == size(tmp,1)

                if strcmpi(LIMO.Type,'Channels') && size(subj_chanlocs(subject_index).chanlocs,1) == size(tmp,1) || ...
                        strcmpi(LIMO.Type,'Channels') && size(subj_chanlocs(subject_index).chanlocs,2) == size(tmp,1)
                    matched_data = limo_match_elec(subj_chanlocs(subject_index).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                elseif  strcmpi(LIMO.Type,'Components')
                    matched_data = tmp(:,begins_at:ends_at,:);
                else
                    if ~isfield(LIMO,'Type')
                        error('LIMO.Type not found - can''t match data across subjects')
                    elseif strcmpi(LIMO.Type,'Channels')
                        error('error while matching data sizes subject %g group %g has %g channels vs, %g in the data',i,h,length(subj_chanlocs(subject_index).chanlocs),size(tmp,1))
                    end
                end

                if matrix_index == 1
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        data(:,:,:,:,matrix_index) = matched_data;
                    else
                        data(:,:,:,matrix_index) = matched_data;
                    end
                else
                    if strcmp(LIMO.Analysis,'Time-Frequency') && ...
                            all(size(matched_data) == size(squeeze(data(:,:,:,:,1))))
                        data(:,:,:,:,matrix_index) = matched_data;
                    elseif ~strcmp(LIMO.Analysis,'Time-Frequency') && ...
                            all(size(matched_data) == size(squeeze(data(:,:,:,1))))
                        data(:,:,:,matrix_index) = matched_data;
                    else
                        error('The data from subject %g have a different size than previous subjects?',i)
                    end
                end
                matrix_index  = matrix_index+1;
                removed{h}(i) = 0; nb_subjects(h) = nb_subjects(h)+1;

                % Use single channel
            elseif strcmp(analysis_type,'1 channel/component only') %&& size(subj_chanlocs(subject_index).chanlocs,2) == size(tmp,1)
                if strcmpi(LIMO.Type,'Channels') && length(subj_chanlocs(subject_index).chanlocs) == size(tmp,1)
                    if length(LIMO.design.electrode) == 1
                        matched_data = limo_match_elec(subj_chanlocs(subject_index).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                    else
                        out = limo_match_elec(subj_chanlocs(subject_index).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                        matched_data = out(i,:,:); % matches the expected chanloc of the subject
                    end
                elseif strcmpi(LIMO.Type,'Components')
                    if size(LIMO.design.component,2) == 1
                        matched_data = tmp(LIMO.design.component,begins_at:ends_at,:); % all param for beta, if con, adjust dim
                    else
                        matched_data = tmp(LIMO.design.component(subject_index),begins_at:ends_at,:); % matches the expected chanloc of the subject
                    end
                end

                if matrix_index == 1
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        data(1,:,:,:,matrix_index) = matched_data;
                    else
                        data(1,:,:,matrix_index) = matched_data;
                    end
                else
                    if strcmp(LIMO.Analysis,'Time-Frequency') && ...
                            all(size(squeeze(matched_data)) == size(squeeze(data(:,:,:,:,1))))
                        data(1,:,:,:,matrix_index) = matched_data;
                    elseif ~strcmp(LIMO.Analysis,'Time-Frequency') && ...
                            all(size(squeeze(matched_data)) == size(squeeze(data(:,:,:,1))))
                        data(1,:,:,matrix_index) = matched_data;
                    else
                        error('The data from subject %g have a different size than previous subjects?',i)
                    end
                end
                matrix_index = matrix_index+1;
                removed{h}(i) = 0; nb_subjects(h) = nb_subjects(h)+1;

            else
                fprintf('subject %g gp %g discarded, channel description and data size don''t match',i,h); disp(' ')
                removed{h}(i) = 1;
            end
            clear tmp Betas con matched_data
            subject_index = subject_index+1;
        end
    end

    % last re-check dimensions
    if sum(single(isnan(data(:)))) == numel(data)
        if exist('errordlg2','file')
            errordlg2('the data matrix is empty! either betas.mat/con.mat files are empty or there is a bug'); return
        else
            errordlg('the data matrix is empty! either betas.mat files are empty or there is a bug'); return
        end
    end
    
    if gp_nb ==1 && size(data,numel(size(data))-1) <= 2
        warning('the concatenated data have %g repeated measures, consider using a t-test',size(data,numel(size(data))-1))
        return
    elseif gp_nb >1 && size(data,numel(size(data))-1) == 1
        warning('the concatenated data have %g groups but only 1 measure, consider using a N-way ANOVA',gp_nb)
        return
    end

    % 4th send relevant info to LIMO_random_robust
    % ----------------------------------------------
    for h=1:gp_nb
        for i=1:size(Paths{h},2)
            if removed{h}(i) == 1
                LIMO.data.data{h}(i) = [];
                LIMO.data.data_dir{h}(i) = []; % somehow to indicate this subject is removed
            end
        end
    end

    % final check
    if sum(nb_subjects) < prod(factor_nb)
       error('there are more variables than observations, some factors can''t be estimated')
    end

    % data dim [channel * frames * all param * subjects]
    % the expected dim in LIMO_rep_anova are [channel * frames * subjects * conditions]
    % + one vector describing the group belonging

    if strcmp(LIMO.Analysis,'Time-Frequency')
        tmp_data = NaN(size(data,1),size(data,2),size(data,3),sum(nb_subjects), prod(factor_nb));
    else
        tmp_data = NaN(size(data,1),size(data,2),sum(nb_subjects), prod(factor_nb));
    end
    gp = NaN(sum(nb_subjects),1);

    for i=1:gp_nb
        % select only relevant parameters (could be different for different groups)
        current_param = parameters(i,:);
        % pick up all subjects of a group
        if i==1
            from = 1;
            to = nb_subjects(i);
        else
            from = from+nb_subjects(i-1);
            to = to+nb_subjects(i);
        end
        gp(from:to) = i;
        % iterate
        for j=1:prod(factor_nb)
            if current_param(j)>size(data,3)
                warning('The parameter %g requested (gp %g) is not valid, beta max=%g ',current_param(j),i,size(data,3))
                return
            end
            if strcmp(LIMO.Analysis,'Time-Frequency')
                tmp_data(:,:,:,from:to,j) = squeeze(data(:,:,:,current_param(j),from:to));
            else
                tmp_data(:,:,from:to,j) = squeeze(data(:,:,current_param(j),from:to));
            end
        end
    end
    clear data

    % finally, since all seems to match up, ask for factor names
    % ---------------------------------------------------------
    if ~isfield(LIMO.design, 'factor_names')
        LIMO.design.factor_names = factor_names;
    end

    if length(LIMO.design.factor_names) ~= length(factor_nb)
        warning('factor names (n=%g) removed as it does not match number of factors entered (n=%g)',length(LIMO.design.factor_names),length(factor_nb));
        LIMO.design = rmfield(LIMO.design,'factor_names');
    end

    clear Betas channeighbstructmat expected_chanlocs Names Paths subj_chanlocs
    cd(LIMO.dir);
    if strcmp(LIMO.Analysis,'Time-Frequency')
        LIMO.data.size4D = size(tmp_data);
        LIMO.data.size3D = [LIMO.data.size4D(1) LIMO.data.size4D(2)*LIMO.data.size4D(3) LIMO.data.size4D(4)];
    end
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
    Yr = tmp_data; save(fullfile(LIMO.dir,'Yr.mat'),'Yr','-v7.3');
    clear tmp_data

    % compute
    % --------
    tmpname = limo_random_robust(6,Yr,gp,factor_nb,LIMO,'go',skip_design_check);
    if nargout ~= 0
        LIMOPath = tmpname;
    end
end

end % closes the function


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                                   ROUTINES
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% fileparts for multiple cell entries in row (ie per group)
function [Names,Paths,Files] = breaklimofiles(cellfiles)
for ifiles = length(cellfiles):-1:1
    if ~isempty(cellfiles{ifiles})
        [Paths{ifiles}, filename, ext] = fileparts(cellfiles{ifiles});
        Names{ifiles}                  = [filename ext];
        Files{ifiles}                  = fullfile(Paths{ifiles},[filename ext]);
    end
end
end

%% frame matching
function [first_frame,last_frame,subj_chanlocs,channeighbstructmat,LIMO] = match_frames(Paths,LIMO)

% once we have all the files, we need to collect information to match the frames across subjects
% OUTPUT: first_frame and last _frame returns the beginning and end in terms of indices
%                                     for time-frequency these are vectors with time then frequency
%         subj_chanlocs the chanlocs per subjects
%         channeighbstructmat the neighbourg matrices
%
% the LIMO structure is also updated the reflect the smallest interval(s) across subjects,
% which is used for the second leve analysis

channeighbstructmat = [];
disp('match frames between subjects ...')
% check Paths format
if iscell(Paths{1})
    tmp = Paths; clear Paths
    index = 1;
    for gp=1:size(tmp,2)
        for s=1:size(tmp{gp},2)
            Paths{index} = tmp{gp}(s);
            index = index + 1;
        end
    end
end

% now loop loading the LIMO.mat for each subject to collect information
% ---------------------------------------------------------------------

ME = [];
for i=size(Paths,2):-1:1
    try
        if iscell(Paths{i})
            cd (cell2mat(Paths{i}));
        else
            cd (Paths{i});
        end
        limo = load('LIMO.mat');
        limo = limo.LIMO;
    catch
        if iscell(Paths{i})
            p=cell2mat(Paths{i});
        else
            p=Paths{i};
        end
        error('cannot load/find the LIMO.mat in %s',p)
    end

    if i==size(Paths,2)
        Analysis = limo.Analysis;
    else
        if ~strcmp(limo.Analysis,Analysis)
            error('Looks like different type of analyses (Time/Freq/Time-Freq) are mixed up')
        end
    end

    sampling_rate(i) = limo.data.sampling_rate;
    if strcmpi(LIMO.Type,'Channels')
        subj_chanlocs(i).chanlocs = limo.data.chanlocs;
        if isfield(limo.data,'channeighbstructmat')
            channeighbstructmat = limo.data.channeighbstructmat;
        end
    else
        subj_chanlocs(i).chanlocs = [];
    end

    if strcmp(Analysis,'Time-Frequency')
        first_frame(i,1)            = limo.data.trim1;
        last_frame(i,1)             = limo.data.trim2;
        start(i,1)                  = limo.data.start;
        stop(i,1)                   = limo.data.end;

        first_frame(i,2)            = limo.data.trim_lowf;
        last_frame(i,2)             = limo.data.trim_highf;
        start(i,2)                  = limo.data.tf_freqs(1);
        stop(i,2)                   = limo.data.tf_freqs(end);

        tf_times{i}(1,:)            = limo.data.tf_times;
        tf_freqs{i}(1,:)            = limo.data.tf_freqs;
    else
        first_frame(i)              = limo.data.trim1;
        last_frame(i)               = limo.data.trim2;
        start(i)                    = limo.data.start;
        stop(i)                     = limo.data.end;

        if strcmp(Analysis,'Frequency')
            freqlist{i}(1,:)        = limo.data.freqlist;
        end
    end
end
clear limo

% quick check things are ok
if  strcmpi(LIMO.Type,'Channels') && ~isempty(ME) && isempty(LIMO.data.neighbouring_matrix)
    error('some subject(s) have a different channel structure \nplease load an expected chanloc when choosing a test');
end

if length(unique(sampling_rate)) ~= 1
    error('data have different sampling rates')
end

% match and return into LIMO
LIMO.Analysis           = Analysis;
LIMO.data.sampling_rate = sampling_rate(1);

% we need 1) to find the highest start in time and freq 2) the lowest end
% in time and freq and 3) match that on freqlist or tf_times/tf_freqs

[v,c] = max(first_frame);
if strcmp(Analysis,'Time-Frequency')
    LIMO.data.trim1 = v(1);
    LIMO.data.start = start(c(1),1);
    LIMO.data.trim_lowf = v(2);
    LIMO.data.lowf = start(c(2),2);
    % loop to adjust each subject list
    for i=1:size(Paths,2)
        [~,ind] = min(abs(tf_times{i}-LIMO.data.start));
        tf_times{i} = tf_times{i}(ind:end);
        [~,ind] = min(abs(tf_freqs{i}-LIMO.data.lowf));
        tf_freqs{i} = tf_freqs{i}(ind:end);
    end
else
    LIMO.data.trim1 = v;
    LIMO.data.start = start(c);
    if strcmp(Analysis,'Frequency')
        % loop to adjust each subject list
        for i=1:size(Paths,2)
            [~,ind] = min(abs(freqlist{i}-LIMO.data.start));
            freqlist{i} = freqlist{i}(ind:end);
        end
    end
end

[v,c] = min(last_frame);
if strcmp(Analysis,'Time-Frequency')
    LIMO.data.trim2 = v(1);
    LIMO.data.end = stop(c(1),1);
    LIMO.data.trim_highf = v(2);
    LIMO.data.highf = stop(c(2),2);
    % loop to adjust each subject list
    for i=1:size(Paths,2)
        [~,ind] = min(abs(tf_times{i}-LIMO.data.end));
        tf_times{i} = tf_times{i}(1:ind);
        [~,ind] = min(abs(tf_freqs{i}-LIMO.data.highf));
        tf_freqs{i} = tf_freqs{i}(1:ind);
    end
else
    LIMO.data.trim2 = v;
    LIMO.data.end = stop(c);
    if strcmp(Analysis,'Frequency')
        % loop to adjust each subject list
        for i=1:size(Paths,2)
            [~,ind] = min(abs(freqlist{i}-LIMO.data.end));
            freqlist{i} = freqlist{i}(1:ind);
        end
    end
end

% finally match everything
if strcmp(Analysis,'Frequency')
    % check all lists match ; if sampled the same across subject, all same sizes
    try
        freqlist = cell2mat(freqlist');
    catch list_issue
        assignin('base','freqlist',freqlist)
        error('the resolution of frequency lists doesn''t match between subjects \n%s',list_iisue.message)
    end
    LIMO.data.freqlist = mean(freqlist,1);
    LIMO.data.start    = LIMO.data.freqlist(1);
    LIMO.data.end      = LIMO.data.freqlist(end);

elseif strcmp(Analysis,'Time-Frequency')
    % check all lists match
    try
        tf_times = cell2mat(tf_times');
        tf_freqs = cell2mat(tf_freqs');
    catch list_issue
        error('the resolution of time/frequency lists doesn''t match between subjects \n%s',list_issue.message)
    end
    LIMO.data.tf_times = mean(tf_times,1);
    LIMO.data.tf_freqs = mean(tf_freqs,1);
    LIMO.data.start    = LIMO.data.tf_times(1);
    LIMO.data.lowf     = LIMO.data.tf_freqs(1);
    LIMO.data.end      = LIMO.data.tf_times(end);
    LIMO.data.highf    = LIMO.data.tf_freqs(end);
end
cd(LIMO.dir)
end

%% match channels and update LIMO
function LIMO = match_channels(stattest,analysis_type,LIMO)

% note channel can be a singleton or a vector (channel/component optimized analysis)
if strcmpi(analysis_type,'1 channel/component only')
    if isempty(LIMO.design.electrode) && strcmp(LIMO.Type,'Channels')
        channel = limo_inputdlg('which electrode to analyse [?]','Electrode option');
    elseif isempty(LIMO.design.electrode) && strcmp(LIMO.Type,'Components')
        channel = limo_inputdlg('which component to analyse [?]','Component option'); % can be 1 nb or a vector of channels (channel optimized analysis)
    else
        if ischar(LIMO.design.electrode)
            LIMO.design.electrode = load(LIMO.design.electrode);
            LIMO.design.electrode = LIMO.design.electrode.(cell2mat(fieldnames(LIMO.design.electrode)));
        end
        channel = {num2str(LIMO.design.electrode)}; % reformat temporarilly as if from limo_inputdlg
    end

    if isempty(cell2mat(channel))
        [file,filepath,index] = uigetfile('*.mat',['select a ' LIMO.Type ' file']);
        if isempty(file) || index == 0
            return
        else
            % check the vector has the same length as the number of files
            channel_vector = load(fullfile(filepath,file));
            tmpname        = fieldnames(channel_vector);
            channel_vector = getfield(channel_vector,tmpname{1}); %#ok<GFLD>
            clear tmpname
            if length(channel_vector) ~= length(Paths)
                if exist('errordlg2','file')
                    errordlg2(['the nb of ' LIMO.Type ' does not match the number of subjects'],'Error');
                else
                    errordlg(['the nb of ' LIMO.Type ' does not match the number of subjects'],'Error');
                end
                return
            end

            % add the name to LIMO if absent
            if ~isfield(LIMO.design,'name')
                if stattest == 1
                    LIMO.design.name = ['one sample t-test one ' LIMO.Type(1:end-1)];
                elseif stattest == 2
                    LIMO.design.name = ['two samples t-test one ' LIMO.Type(1:end-1)];
                elseif stattest == 3
                    LIMO.design.name = ['paired t-test one ' LIMO.Type(1:end-1)];
                elseif stattest == 4
                    LIMO.design.name = ['regression analysis one ' LIMO.Type(1:end-1)];
                elseif stattest == 5
                    LIMO.design.name = ['AN(C)OVA analysis one ' LIMO.Type(1:end-1)];
                end
            end

            % restric the channels
            if strcmp(LIMO.Type,'Channels')
                LIMO.data.chanlocs          = LIMO.data.expected_chanlocs;
                LIMO.data.expected_chanlocs = LIMO.data.expected_chanlocs(channel_vector);
                LIMO.design.electrode       = channel_vector;
            else
                LIMO.data.chanlocs    = [];
                LIMO.design.component = channel_vector;
            end
        end

    elseif numel(str2num(channel{1})) || ...
            max(size(cell2mat(channel))) == numel(LIMO.data.data) || ...
            max(size(cell2mat(channel))) == sum(cellfun(@numel,LIMO.data.data))

        if ~isfield(LIMO.design,'name')
            if stattest == 1
                LIMO.design.name = ['one sample t-test one ' LIMO.Type(1:end-1)];
            elseif stattest == 2
                LIMO.design.name = ['two samples t-test one ' LIMO.Type(1:end-1)];
            elseif stattest == 3
                LIMO.design.name = ['paired t-test one ' LIMO.Type(1:end-1)];
            elseif stattest == 4
                LIMO.design.name = ['regression analysis one ' LIMO.Type(1:end-1)];
            else
                LIMO.design.name = ['AN(C)OVA analysis one ' LIMO.Type(1:end-1)];
            end
        end

        if strcmp(LIMO.Type,'Channels')
            LIMO.design.electrode       = str2num(cell2mat(channel));
            LIMO.data.chanlocs          = LIMO.data.expected_chanlocs;
            LIMO.data.expected_chanlocs = LIMO.data.expected_chanlocs(LIMO.design.electrode);
        else
            LIMO.design.component = eval(cell2mat(channel));
            LIMO.data.chanlocs    = [];
        end
    else
        error(['the nb of ' LIMO.Type ' does not match the number of subjects'],[LIMO.Type(1:end-1) ' error']);
    end

    % ---------------
else % Full scalp
    % ---------------

    if ~isfield(LIMO.design,'name')
        if stattest == 1
            LIMO.design.name = ['One sample t-test all ' LIMO.Type];
        elseif stattest == 2
            LIMO.design.name = ['Two samples t-test all ' LIMO.Type];
        elseif stattest == 3
            LIMO.design.name = ['Paired t-test all ' LIMO.Type];
        elseif stattest == 4
            LIMO.design.name = ['Regression analysis all ' LIMO.Type];
        else
            LIMO.design.name = ['AN(C)OVA analysis all ' LIMO.Type];
        end
    end

    if isfield(LIMO.data,'expected_chanlocs')
        LIMO.data.chanlocs = LIMO.data.expected_chanlocs;
    end

    if strcmpi(LIMO.Type,'Components')
        LIMO.design.component = [];
    else
        LIMO.design.electrode = [];
    end
end
end

%% assemble the data matrix
function [data,removed] = getdata(stattest,analysis_type,first_frame,last_frame,subj_chanlocs,LIMO)

data = [];
removed = [];
disp('gathering data ...');
if stattest == 1 % one sample
    index = 1;
    if all(size(LIMO.data.data)==[1 1]) % cell of cell
        LIMO.data.data = LIMO.data.data{1};
    end

    for i=1:size(LIMO.data.data,2) % for each subject
        tmp = load(LIMO.data.data{i});

        % get indices to trim data
        if strcmp(LIMO.Analysis,'Time-Frequency')
            if contains(LIMO.data.data{i},'Betas')
                tmp = tmp.Betas;
            else
                tmp = tmp.con(:,:,:,1);
            end
            begins_at = fliplr((max(first_frame) - first_frame(i,:) + 1)); % returns time/freq/or freq-time
            ends_at(1) = size(tmp,2) - (last_frame(i,2) - min(last_frame(:,2)));
            ends_at(2) = size(tmp,3) - (last_frame(i,1) - min(last_frame(:,1)));
        else
            if contains(LIMO.data.data{i},'Betas')
                tmp = tmp.Betas;
            else
                tmp = tmp.con(:,:,1);
            end
            begins_at = max(first_frame) - first_frame(i) + 1;
            ends_at = size(tmp,2) - (last_frame(i) - min(last_frame));
        end

        % data dim [channel, freq/time, param, nb subjects]
        if isempty(LIMO.Type)
            LIMO.Type = limo_questdlg('Is the analysis on','Type is empty','Channels','Components','Channels');
        end

        if strcmp(analysis_type,'Full scalp analysis')
            if strcmpi(LIMO.Type,'Channels') && length(subj_chanlocs(i).chanlocs) == size(tmp,1)
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    data(:,:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                else
                    data(:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                end
            elseif strcmpi(LIMO.Type,'Components')
                try
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        data(:,:,:,:,index) = tmp(:,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                    else
                        data(:,:,:,index) = tmp(:,begins_at:ends_at,:);
                    end
                catch dim_error
                    if strcmp(dim_error,'Subscripted assignment dimension mismatch.')
                        disp('you are trying to match matrices of ICs of different size')
                        if isempty(expected_chanlocs)
                            disp('either cluster data are run 1st level batch, or input cluster file')
                        end
                    end
                end
            end
            index = index + 1;
            removed(i) = 0;

        elseif strcmp(analysis_type,'1 channel/component only') %&& size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)

            % Use single channel
            if  length(LIMO.design.electrode) == 1
                if strcmpi(LIMO.Type,'Channels') && length(subj_chanlocs(i).chanlocs) == size(tmp,1)
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta
                        else
                            data(1,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                        end
                    else
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                        else
                            data(1,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                        end
                    end
                elseif strcmpi(LIMO.Type,'Components')
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,:,index) = tmp(LIMO.design.component,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                        else
                            data(1,:,:,index) = tmp(LIMO.design.component,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                        end
                    else
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,index) = tmp(LIMO.design.component,begins_at:ends_at,:);
                        else
                            data(1,:,index) = tmp(LIMO.design.component,begins_at:ends_at,:);
                        end
                    end
                end
                index = index + 1;
                removed(i) = 0;

            else  % Use multiple single channels
                if strcmpi(LIMO.Type,'Channels')
                    out = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,:,index) = out(i,:,:,:);
                        else
                            data(1,:,:,index) = out(i,:,:,:);
                        end
                    else
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,index) = out(i,:,:); % matches the expected chanloc of the subject
                        else
                            data(1,:,index) = out(i,:,:); % matches the expected chanloc of the subject
                        end
                    end
                elseif strcmpi(LIMO.Type,'Components')
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,:,index) = tmp(LIMO.design.component(i),begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                        else
                            data(1,:,:,index) = tmp(LIMO.design.component(i),begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                        end
                    else
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,index) = tmp(LIMO.design.component(i),begins_at:ends_at,:);
                        else
                            data(1,:,index) = tmp(LIMO.design.component(i),begins_at:ends_at,:);
                        end
                    end
                end
                index = index +1;
                removed(i) = 0;
            end
        else
            fprintf('subject %g discarded, channel description and data size don''t match \n',i);
            removed(i) = 1;
        end
        clear tmp
    end

elseif stattest == 2 % several samples
    subject_nb = 1;
    for igp = 1:length(LIMO.data.data)
        index = 1;
        for i=1:size(LIMO.data.data{igp},2) % for each subject per group
            tmp = load(cell2mat(LIMO.data.data{igp}(i)));

            % get indices to trim data
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                    tmp = tmp.Betas;
                else
                    tmp = tmp.con(:,:,:,1);
                end
                begins_at = fliplr((max(first_frame) - first_frame(subject_nb,:) + 1)); % returns time/freq/or freq-time
                ends_at(1) = size(tmp,2) - (last_frame(subject_nb,2) - min(last_frame(:,2)));
                ends_at(2) = size(tmp,3) - (last_frame(subject_nb,1) - min(last_frame(:,1)));
            else
                if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                    tmp = tmp.Betas;
                else
                    tmp = tmp.con(:,:,1);
                end
                begins_at = max(first_frame) - first_frame(subject_nb) + 1;
                ends_at = size(tmp,2) - (last_frame(subject_nb) - min(last_frame));
            end

            if strcmpi(analysis_type,'Full scalp analysis') %&& size(subj_chanlocs(subject_nb).chanlocs,2) == size(tmp,1)
                if strcmpi(LIMO.Type,'Channels') && length(subj_chanlocs(subject_nb).chanlocs) == size(tmp,1)
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        tmp_data(:,:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                    else
                        tmp_data(:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                    end
                elseif strcmpi(LIMO.Type,'Components')
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        tmp_data(:,:,:,:,index) = tmp(:,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                    else
                        tmp_data(:,:,:,index) = tmp(:,begins_at:ends_at,:);
                    end

                end
                removed(igp,i) = 0;
                index = index + 1;

            elseif strcmpi(analysis_type,'1 channel/component only') %&& size(subj_chanlocs(subject_nb).chanlocs,2) == size(tmp,1)

                % Use single channel
                if length(LIMO.design.electrode) == 1
                    if strcmpi(LIMO.Type,'Channels') && length(subj_chanlocs(subject_nb).chanlocs) == size(tmp,1)
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                            else
                                tmp_data(1,:,:,1,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                            end
                        else
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                            else
                                tmp_data(1,:,1,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                            end
                        end
                    elseif strcmpi(LIMO.Type,'Components')
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,:,index) = tmp(LIMO.design.electrode,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:); % all param for beta, if con, adjust dim
                            else
                                tmp_data(1,:,:,1,index) = tmp(LIMO.design.electrode,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:); % all param for beta, if con, adjust dim
                            end
                        else
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,index) = tmp(LIMO.design.electrode,begins_at:ends_at,:); % all param for beta, if con, adjust dim
                            else
                                tmp_data(1,:,1,index) = tmp(LIMO.design.electrode,begins_at:ends_at,:); % all param for beta, if con, adjust dim
                            end
                        end
                    end
                    removed(igp,i) = 0;
                    index = index + 1;

                else % Use multiple single channels
                    if strcmpi(LIMO.Type,'Channels')
                        out = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,:,index) = out(subject_nb,:,:,:); % matches the expected chanloc of the subject
                            else
                                tmp_data(1,:,:,1,index) = out(subject_nb,:,:,:); % matches the expected chanloc of the subject
                            end
                        else
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,index) = out(subject_nb,:,:);     % matches the expected chanloc of the subject
                            else
                                tmp_data(1,:,1,index) = out(subject_nb,:,:);     % matches the expected chanloc of the subject
                            end
                        end
                    elseif strcmpi(LIMO.Type,'Components')
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,:,index) = tmp(LIMO.design.electrode(subject_nb),begins_at(1):ends_at(1),begins_at(2):ends_at(2),:); % matches the expected chanloc of the subject
                            else
                                tmp_data(1,:,:,1,index) = tmp(LIMO.design.electrode(subject_nb),begins_at(1):ends_at(1),begins_at(2):ends_at(2),:); % matches the expected chanloc of the subject
                            end
                        else
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,index) = tmp(LIMO.design.electrode(subject_nb),begins_at:ends_at,:);     % matches the expected chanloc of the subject
                            else
                                tmp_data(1,:,1,index) = tmp(LIMO.design.electrode(subject_nb),begins_at:ends_at,:);     % matches the expected chanloc of the subject
                            end
                        end
                    end
                    removed(igp,i) = 0;
                    index = index +1;
                end
            else
                fprintf('subject %g of group %g discarded, channel description and data size don''t match \n',i, igp);
                removed(igp,i) = 1;
            end
            clear tmp
            subject_nb = subject_nb + 1;
        end
        data{igp} = tmp_data;
        clear tmp tmp_data
    end
end
end

%% repeated measure levels
function levels = getlevels(params)
% drill-down the cell array of repeated measures parameters
if ~iscell(params)
    levels = length(params);
    if levels == 1
        levels = [];
    end
else
    levels = [length(params) getlevels(params{1}) ];
end
end

% get beta indices from study
% ---------------------------
function [param,betas] = get_beta_indices(selectmode, betaFile,factorname,factorn)

param = [];
betas = {};
if nargin < 3
    factorname = {};
end

if nargin > 1
    try
        limoFile = strrep(betaFile, 'Betas.mat', 'LIMO.mat');
        LIMO = load('-mat', limoFile);
        betas = { LIMO.LIMO.design.labels.description };
    catch
        disp('Warning: could not find associated LIMO file');
    end
end

if ~isempty(betas)
    % simple selection
    for iBeta = 1:length(betas)
        betas{iBeta} = [ int2str(iBeta) ' - ' betas{iBeta}];
    end
    if strcmpi(selectmode, 'selectone')
        strBeta = 'Pick a single beta parameters below';
        strInds = 'Or ignore selection above and enter beta index';
        maxval = 1;
    elseif strcmpi(selectmode, 'selecttwo')
        strBeta = 'Pick TWO beta parameters below';
        strInds = 'Or ignore selection above and enter beta index';
        maxval = 2;
    else
        strBeta = 'Pick one or more beta parameters below';
        strInds = 'Or ignore selection above and enter beta indices';
        maxval = 2;
    end
    uiList = { {'style' 'text' 'string' strBeta } ...
        { 'style' 'listbox' 'string' betas 'max' maxval } ...
        {'style' 'text' 'string' strInds } ...
        {'style' 'edit' 'string' '' } };
    res = inputgui('uilist', uiList, 'geometry', { [1] [1] [3 1] }, 'geomvert', [1 length(betas)/2+1 1]);
    if isempty(res), return; end
    if ~isempty(res{2})
        param =  eval( [ '[' res{2} ']' ] );
    else
        param =  res{1};
    end
else
    param = eval( [ '[' cell2mat(limo_inputdlg('which parameters to test e.g [1:3]','parameters option')) ']' ]);
end
end

%% reorder parameters
function parameters = reorganize_params(parameters, betas, factorname, factorn)

% simple selection
str = ['Reorganize values as needed: ' factorname{1} '=rows and ' factorname{2} '=columns'];
uiList = { {'style' 'text' 'string' str 'fontweight' 'bold'} };
uiGeom = { [1] };
uiVert = 1;
for iRow = 1:length(parameters)
    uiList = [ uiList(:)' ...
        {{ 'style' 'popupmenu' 'string' betas(parameters) 'value' iRow }} ];
    if mod(iRow, factorn(1)) == 0
        uiGeom = [ uiGeom(:)' {ones(1,factorn(1))} ];
    end
end
res = inputgui('uilist', uiList, 'geometry', uiGeom, 'minwidth', 800);

if isempty(res), parameters = []; return; end
parameters = parameters(cell2mat(res));
end

%% file checking
function [parameters,betas] = check_files(Paths,Names,gp,parameters,selectmode)
% after selecting file, check they are all the same type (betas or con)
% return parameters that match with files (eg 1 for con, or whatever value 
% for the beta file)

betas = {};
if nargin < 4
    parameters = [];
end
if nargin < 5 || isempty(selectmode)
    selectmode = 'multi';
end
if nargin < 6
    factorname = '';
    factorn    = [];
end

if gp == 1
    if iscell(Names{gp})
        Names = Names{gp};
        Paths = Paths{gp};
    end

    % one sample case
    % ---------------
    is_beta = []; is_con = [];
    for i=size(Names,2):-1:1
        if strfind(Names{i},'Betas')
            is_beta(i) = 1;
        elseif strfind(Names{i},'con')
            is_con(i) = 1;
        end
    end

    if (isempty(is_beta)) == 0 && sum(is_beta) ~= size(Names,2) || (isempty(is_con)) == 0 && sum(is_con) ~= size(Names,2)
        errordlg2('file selection failed, only Beta or Con files are supported'); return
    elseif (isempty(is_beta)) == 0 && sum(is_beta) == size(Names,2) && nargout ~= 0
        if isempty(parameters)
            if isempty(factorname) || length(factorname) > 2
                [parameters,betas] = get_beta_indices(selectmode, fullfile(Paths{1}, Names{1}));
            end
            %parameters = { [1 2 3] [4 5 6] [7 8 9] };
            if isempty(parameters)
                return
            end
        end
    elseif (isempty(is_con)) == 0 && sum(is_con) == size(Names,2)
        parameters = 1;
    end

elseif gp > 1

    % several samples case
    % -------------------
    for g = gp:-1:1
        is_beta = []; is_con = [];
        for i=1:size(Names{g},2)
            if contains(Names{g}(i),'Betas')
                is_beta(i) = 1;
            elseif strfind(Names{g}{i},'con')
                is_con(i) = 1;
            end
        end

        if ~isempty(is_beta)
            test{g} = sum(is_beta) == size(Names{g},2);
        elseif ~isempty(is_con)
            test{g} = sum(is_con) == size(Names{g},2);
        end
    end

    if sum(cell2mat(test)) ~= length(Names)
        error('file selection failed, only sets of Beta or sets of Con files are supported');
    elseif ~isempty(is_beta) && sum(cell2mat(test)) == length(Names) && nargout ~= 0
        if isempty(parameters)
            parameters = eval(cell2mat(limo_inputdlg('which parameter(s) to test e.g 1','parameters option')));
        elseif ~isempty(parameters) && size(parameters,2) ~=1 && size(parameters,2) ~=gp
            warndlg2('A valid parameter value must be provided - selection aborded');
            return
        end
        if isempty(parameters) || size(parameters,2) ~=1 && size(parameters,2) ~=gp
            warndlg2('A valid parameter value must be provided - selection aborded');
            return
        end
    elseif ~isempty(is_con) && sum(cell2mat(test)) == length(Names)
        if isempty(parameters)
            parameters = 1;
        else
            parameters = ones(1,length(parameters));
        end
    end
end

if isempty(betas) && isempty(parameters)
    error('LIMO could not match betas across subjects, the array is empty')
end

end
