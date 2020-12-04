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
% FORMAT: LIMO_random_select(stattest,expected_chanlocs)
%         LIMO_random_select(stattest,expected_chanlocs,options)
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
%                             dimensions of the cell correspond with group, factor and
%                             level respectively. (no default)
%                'parameters' Cell array of parameters to be tested, relative to LIMOfiles.
%                            ie. {[1 2]} or {[1 2];[1 2]} in case of 2
%                            groups. Add nested cells for more repetition levels.
%       --> for LIMOfiles and parameters the rule is groups in rows, repeated measures in columns
%                'regressor_file' a file or matrix of data to regress when stattest = 4
%                'analysis_type' is 'Full scalp analysis' or '1 channel/component only'
%                'channel' Index of the electrode(s) to use if '1 channel/component only'
%                            is selected in analysis_type
%                'type' is 'Channels' or 'Component'
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
%
% Example for a repeated measure ANOVA with command line
% LIMOPath = limo_random_select('Repeated Measures ANOVA',chanlocs,'LIMOfiles',...
%     {'F:\WakemanHenson_Faces\eeg\derivatives\LIMO_Face_detection\Beta_files_FaceRepAll_GLM_Channels_Time_WLS.txt'},...
%     'analysis type','Full scalp analysis','parameters',{[1 2 3],[4 5 6],[7 8 9]},...
%     'factor names',{'face','repetition'},'type','Channels','nboot',0,'tfce',0);
%
% Cyril Pernet - The University of Edinburgh
% Ramon Martinez-Cancino - UC San Diego
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
regressor_file         = [];
analysis_type          = [];
zopt                   = [];
skip_design_check      = 'No';

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
        LIMO.design.parameters = varargin{in+1};
    elseif strcmpi(varargin{in},'factor names')
        LIMO.design.factor_names = varargin{in+1};
    elseif strcmpi(varargin{in},'type')
        LIMO.Type = varargin{in+1};
    elseif strcmpi(varargin{in},'nboot')
        LIMO.design.bootstrap = varargin{in+1};
    elseif strcmpi(varargin{in},'tfce')
        LIMO.design.tfce = varargin{in+1};
    end
end

if isempty(analysis_type)
    analysis_type   = questdlg('Rdx option','type of analysis?','Full scalp analysis','1 channel/component only','Full scalp analysis');
    if isempty(analysis_type)
        return
    end
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
        if ischar(LIMO.data.data{1}) % Case for path to the files
            [Names,Paths,LIMO.data.data] = limo_get_files([],[],[],LIMO.data.data{1});
        else % Case when all paths are provided
            [Names,Paths,LIMO.data.data] = breaklimofiles(LIMO.data.data);
        end
    end
    LIMO.data.data_dir = Paths;
    
    if isempty(Names)
        disp('no files selected')
        return
    end
    
    % check type of files and returns which beta param to test
    % -------------------------------------------------------
    if ~isfield(LIMO.design,'parameters')
        parameters = check_files(Names,1);
    else
        parameters = check_files(Names,1,LIMO.design.parameters);
    end
    
    if isempty(parameters)
        errordlg2('file selection failed, only Beta and Con files are supported','Selection error'); return
    end
    
    % match frames, update LIMO
    % --------------------------
    [first_frame,last_frame,subj_chanlocs,~,LIMO] = match_frames(Paths,LIMO);
    
    % match channels, update LIMO
    % -----------------------------
    LIMO = match_channels(1,analysis_type,LIMO);
    
    % get data for all parameters
    % -----------------------------
    [data,removed] = getdata(1,analysis_type,first_frame,last_frame,subj_chanlocs,LIMO);
    
    % if regression get regressor(s)
    % -------------------------------
    if strcmpi(stattest,'regression')
        if isempty(regressor_file)
            [FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','select regressor file');
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
        
        if strcmp(FileName(end-3:end),'.txt')
            X = load(fullfile(PathName,FileName));
        elseif strcmp(FileName(end-3:end),'.mat')
            X = load(fullfile(PathName,FileName));
            X = X.(cell2mat(fieldnames(X)));
        end
        disp('Regressor(s) loaded');
        
        % check size and orientation
        if strcmp(LIMO.Analysis,'Time-Frequency')
            N = size(data,5);
        else
            N = size(data,4);
        end
        
        if size(X,2) == N || size(X,2) == size(Paths,2)
            disp('X has been transposed'); X = X';
        end
        
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
                    disp('covariate adjusted for delete subjects');
                catch ME
                    errordlg2(sprintf('the number of regression value %g differs from the number of subjects %g',size(X,1),N),'Covariate error');
                    fprintf('%s',ME.message); return
                end
            end
            errordlg2(sprintf('the number of regression value %g differs from the number of subjects %g',size(X,1),N),'Covariate error');
        end
        
        if size(X,2)==1 && LIMO.design.bootstrap < 599
            if LIMO.design.bootstrap ~= 0
                LIMO.design.bootstrap = 599;
                disp('nb of bootstrap adjusted to 599 for a simple regression');
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
            if strcmp(LIMO.Analysis,'Time-Frequency')
                tmp               = squeeze(data(:,:,:,parameters(i),:));
                tmp_data          = NaN(1,size(tmp,1),size(tmp,2),size(tmp,3)); % add dim 1 = 1 channel
                tmp_data(1,:,:,:) = tmp; clear tmp;
            else
                tmp               = squeeze(data(:,:,parameters(i),:));
                tmp_data          = NaN(1,size(tmp,1),size(tmp,2));
                tmp_data(1,:,:)   = tmp; clear tmp;
            end
        else
            if strcmp(LIMO.Analysis,'Time-Frequency')
                tmp_data          = squeeze(data(:,:,:,parameters(i),:));
            else
                tmp_data          = squeeze(data(:,:,parameters(i),:));
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
            LIMO.design.method = 'Trimmed mean';
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
            save(fullfile(LIMO.dir,'Yr.mat'),'Yr','-v7.3');
            tmpname = limo_random_robust(1,fullfile(LIMO.dir,'Yr.mat'),...
                parameters(i),LIMO); % ,'zscore',zopt,'go',skip_design_check);
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
            LIMO.design.name = 'Robust regression';
            save(fullfile(LIMO.dir,'Yr.mat'),'Yr','-v7.3');
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
            if isempty(zopt)
                tmpname = limo_random_robust(4,Yr,X(~isnan(sum(X,2)),:),...
                parameters(i),LIMO,'go',skip_design_check);
            else
                tmpname = limo_random_robust(4,Yr,X(~isnan(sum(X,2)),:),...
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
        [Names,Paths,LIMO.data.data{1}] = limo_get_files(1);
        LIMO.data.data_dir{1}           = Paths;
    else
        if ischar(LIMO.data.data{1}) % Case for path to the files
            [Names,Paths,LIMO.data.data{1}] = limo_get_files([],[],[],LIMO.data.data{1});
            LIMO.data.data_dir{1}           = Paths;
        else % Case when all paths are provided
            [Names,Paths]         = breaklimofiles(LIMO.data.data{1});
            LIMO.data.data_dir{1} = Paths;
        end
    end
    if isempty(Names)
        return
    end
    % N = size(Names,2);
    
    % check type of files and returns which beta param to test
    % -------------------------------------------------------
    if isfield(LIMO.design,'parameters')
        parameters =  cell2mat(LIMO.design.parameters);
    else
        LIMO.design.parameters = [];
        parameters             = [];
    end
    
    if isempty(parameters)
        parameters = check_files(Names,1);
        if isempty(parameters)
            return
        elseif length(parameters) > 1
            warning on
            warning('only 1 parameter at a time for two-samples analysis, restricting to the 1st indicated')
            parameters = parameters(1);
        end
    else
        parameters(1) = check_files(Names,1,parameters(1));
    end
    
    if length(LIMO.data.data) == 1
        [Names,Paths,LIMO.data.data{2}] = limo_get_files(2);
        LIMO.data.data_dir{2}           = Paths;
    else
        if ischar(LIMO.data.data{2}) % Case for path to the files
            [Names,Paths,LIMO.data.data{2}] = limo_get_files([],[],[],LIMO.data.data{2});
            LIMO.data.data_dir{2}           = Paths;
        else % Case when all paths are provided
            [Names,Paths]         = breaklimofiles(LIMO.data.data{2});
            LIMO.data.data_dir{2} = Paths;
        end
    end
    if isempty(Names)
        return
    end
    % N = N+size(Names,2);
    
    % check type of files and returns which beta param to test
    % -------------------------------------------------------
    if length(parameters) == 1
        parameters(2) = check_files(Names,1);
    else
        parameters(2) = check_files(Names,1,parameters(2));
    end
    
    if ~isfield(LIMO.design,'parameters')
        LIMO.design.parameters = parameters;
    end
    
    % check type of files and returns which beta param to test
    % -------------------------------------------------------
    for gp = 1:2
        for sub=1:size(LIMO.data.data_dir{gp},2)
            sub_LIMO = load(cell2mat(fullfile(LIMO.data.data_dir{gp}(sub),'LIMO.mat')));
            if parameters(gp) > size(sub_LIMO.LIMO.design.X,2)-1
                error('invalid parameter %g - design subject %s inconsistent',parameters(gp),cell2mat(fullfile(LIMO.data.data_dir{gp}(sub))));
            end
        end
    end
    
    % match frames, update LIMO
    [first_frame,last_frame,subj_chanlocs,~,LIMO] = match_frames([LIMO.data.data_dir{1}' ; LIMO.data.data_dir{2}']',LIMO);
    
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
            errordlg('file selection is corrupted, data sizes don''t match');
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
            errordlg('file selection is corrupted, data sizes don''t match');
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
        [Names,Paths,LIMO.data.data] = limo_get_files;
        LIMO.data.data_dir           = Paths;
    else
        if ischar(LIMO.data.data{1}) % Case for path to the files
            [Names,Paths,LIMO.data.data{1}] = limo_get_files([],[],[],LIMO.data.data{1});
            LIMO.data.data_dir{1}           = Paths;
        else % Case when all paths are provided
            [Names,Paths]         = breaklimofiles(LIMO.data.data{1});
            LIMO.data.data_dir{1} = Paths;
        end
    end
    if isempty(Names)
        return
    end
    % N = size(Names,2);
    
    % check type of files and returns which beta param to test
    % -------------------------------------------------------
    if isfield(LIMO.design,'parameters')
        parameters = cell2mat(LIMO.design.parameters);
    else
        LIMO.design.parameters = [];
        parameters             = [];
    end
    
    if isempty(parameters)
        parameters = check_files(Names,1);
    else
        parameters = check_files(Names,1,parameters(1));
    end
    
    if size(parameters,2) == 1 % either con file, or command line beta with one parameter
        n = Names; clear Names; Names{1} = n; clear n;
        p = Paths; clear Paths; Paths{1} = p; clear p;
        if ~iscell(LIMO.data.data{1}) % size(LIMO.data.data,2) == 1
            d = LIMO.data.data; LIMO.data = rmfield(LIMO.data,'data'); LIMO.data.data{1} = d; clear d;
            d = LIMO.data.data_dir; LIMO.data = rmfield(LIMO.data,'data_dir'); LIMO.data.data_dir{1} = d; clear d;
        end
        
        if size(LIMO.data.data,1) == 1 % only 1st group of files in {1}
            [Names{2},Paths{2},LIMO.data.data{2}] = limo_get_files([],[],'select paired file');
        else
            if ischar(LIMO.data.data{2})% Case for path to the files
                [Names{2},Paths{2},LIMO.data.data{2}] = limo_get_files([],[],[],LIMO.data.data{2});
            else % Case when all paths are provided
                [Names{2},Paths{2}] = breaklimofiles(LIMO.data.data{2});
            end
        end
        
        if isempty(Names{2})
            return
        else
            if ~isempty(LIMO.design.parameters)
                % hack only availbale if beta files and command line argument // not allowed otherwise because it's a paired design
                parameters(2) = check_files(Names{2},1,parameters(2));
            else
                newparameters = check_files(Names{2},1);
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
        LIMO.data.data_dir{2} = Paths{2};
        % N = N + size(Names{2},2);
        if size(Names{1},2) ~= size(Names{2},2)
            errordlg('the nb of files differs between pairs 1 and 2','Paired t-test error'); return
        end
    elseif size(parameters,2) ~=2 % if it was beta file one needs a pair of parameters
        errordlg('2 parameters must be selected for beta files','Paired t-test error'); return
    else % leave it as is, but check this is valid
        for s = 1:size(Paths,2)
            sub_LIMO = load(fullfile(cell2mat(Paths(s)),'LIMO.mat'));
            if max(parameters) > size(sub_LIMO.LIMO.design.X,2)
                errordlg('invalid parameter(s)','Paired t-test error'); return
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
                tmp_data1(:,:,:,1,:) = squeeze(data{1}(:,:,:,:));
                tmp_data2(:,:,:,1,:) = squeeze(data{2}(:,:,:,:));
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
                tmp_data1(:,:,1,:) = squeeze(data{1}(:,:,:));
                tmp_data2(:,:,1,:) = squeeze(data{2}(:,:,:));
            end
        end
    end
    
    LIMO.design.method = 'Yuen t-test (Trimmed means)';
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
    tmpname = limo_random_robust(3,squeeze(tmp_data1),squeeze(tmp_data2),parameters,LIMO);
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
        gp_nb = size(LIMO.data.data,1);
    else
        gp_nb = cell2mat(inputdlg('How many independent groups? e.g. 3 or [3 2] for nested gps','Groups'));
        if isempty(gp_nb)
            return
        elseif sum(str2double(gp_nb) <= 2) && strcmpi(stattest,'N-Ways ANOVA')
            errordlg('at least 3 groups are expected for a N-ways ANOVA')
            return
        elseif sum(str2double(gp_nb) <= 1) && strcmpi(stattest,'ANCOVA')
            errordlg('at least 2 groups are expected for an ANCOVA')
            return
        else
            gp_nb          = str2double(gp_nb);
            a              = questdlg('load con files or beta file','ANOVA loading files','con','beta','beta');
            Names          = cell(gp_nb,1);
            Paths          = cell(gp_nb,1);
            LIMO.data.data = cell(gp_nb,1);
        end
    end
    
    % select data per gp / conditions
    % ---------------------------------
    for i=1:gp_nb
        if isempty(LIMO.data.data{i})
            if strcmp(a,'beta') % beta files
                [Names{i},Paths{i},LIMO.data.data{i}] = limo_get_files([' beta file gp ',num2str(i)]);
            else
                [Names{i},Paths{i},LIMO.data.data{i}] = limo_get_files([' con file gp ',num2str(i)]);
            end
        else
            if ischar(LIMO.data.data{i}) % Case for path to the files
                [Names{i},Paths{i},LIMO.data.data{i}] = limo_get_files([],[],[],LIMO.data.data{i});
            else % Case when all paths are provided
                [Names{i},Paths{i},LIMO.data.data{i}] = breaklimofiles(LIMO.data.data{i});
            end
        end        
        if isempty(Names{i}); return; end
    end
    
    if isempty(LIMO.design.parameters)
        param = cell2mat(inputdlg('which parameters to test e.g [1 3 3]','parameters option'));
        if isempty(param)
            disp('selection aborded'); return
        else
            if contains(param,'[') && contains(param,']')
                parameters = eval(param);
            else
                parameters = eval(['[' param ']']);
            end
        end
        
        if length(parameters) == 1
            parameters = repmat(parameters,1,gp_nb);
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
        error('%g parameter selected, longer than the number of groups %g',length(parameters),gp_nb)
    else
        parameters = check_files(Names,size(Names,2),parameters);
    end
    
    if isempty(parameters)
        error('selected files are not Beta or con files');
    end
    LIMO.data.data_dir     = Paths;
    LIMO.design.parameters = parameters;
    
    % organize the data in such a way that we can easily compute stuff
    % ---------------------------------------------------------------------
    % match frames, update LIMO
    [first_frame,last_frame,subj_chanlocs,~,LIMO] = match_frames(Paths,LIMO);
    
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
        [FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','select covariate file');
        if FilterIndex == 0
            return
        else
            X = load(fullfile(PathName,FileName));
            if strcmp(FileName(end-3:end),'.mat')
                X = X.(cell2mat(fieldnames(X)));
            end
        end
        disp('Regressor(s) loaded');
        
        try
            if size(X,2) == sum(nb_subjects)
                disp('regressor transposed');
                X = X';
            end
            
            if size(X,1) ~= size(data,ndims(data))
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
                    if size(X,1) ~= size(data,4)
                        errordlg('the number of regression value differs from the number of subjects'); return
                    else
                        errordlg(sprintf('log error:%s',ME.message),'fail adjusting covarate(s)'); retun
                    end
                end
            end
        catch ME
            errordlg('data format or length of the covariate does not match');
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
                tmp_data = NaN(size(data{i},[1 2 3 5]));
            end
            tmp_data(:,:,:,index:(sum(nb_subjects(1:i))))  = squeeze(data{i}(:,:,:,current_param,:));
        else
            if i==1
                tmp_data = NaN(size(data{i},[1 2 4]));
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
    % -------------
    if ~isempty(LIMO.data.data)
        gp_nb = size(LIMO.data.data,1);
    else
        gp_nb = cell2mat(inputdlg('How many independent groups? e.g. 2','Groups'));
    end
    
    if isempty(gp_nb)
        return
    elseif length(gp_nb) > 1
        errordlg2('only 1 independent factor (with n groups) is handled by LIMO')
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
            
        for g=gp_nb:-1:1
            if all(size(LIMO.design.parameters(g,:))==1)
                factor_nb{g} = num2str(length(LIMO.design.parameters{g}));
            else
                factor_nb{g} = num2str(getlevels(LIMO.design.parameters(g,:)));
            end
        end
        
        % check all factor numbers match between groups and reduce
        if ~all(cellfun(@(x) strcmpi(x,factor_nb{1}),factor_nb))
            error('parameters input sizes different between groups, the number of factors to infer must be identical')
        else
            factor_nb = factor_nb{1};
        end
    else
        factor_nb = cell2mat(inputdlg('Enter repeated factors level? e.g. [2 3] for 2 levels F1 and 3 levels F2','Factors'));
    end
    
    if isempty(factor_nb)
        disp('selection aborded'); return
    else
        if contains(factor_nb,'[') && contains(factor_nb,']')
            factor_nb = eval(factor_nb);
        else
            try
                factor_nb = eval(['[' factor_nb ']']);
            catch ME
                errordlg2(sprintf('log error: %s',ME.message),'could not evaluate factors')
                return
            end
        end
    end
    
    % Cases of wrong input
    if isempty(factor_nb) || length(factor_nb)==1 && factor_nb == 0
        disp('no factor entered, Rep. ANOVA aborded');
        return
    end
    
    % 2nd select data per gp / conditions
    % ---------------------------------------------------------------
    if ~isempty(LIMO.data.data)
        if ischar(LIMO.data.data{1}) % && ischar(LIMO.data.data{2})
            isbeta = cellfun(@(x) contains(x,'Beta'),LIMO.data.data);
            iscon  = cellfun(@(x) contains(x,'con'),LIMO.data.data);
        else
            for gp = 1:gp_nb
                isbeta = mean(cellfun(@(x) contains(x,'Beta'),LIMO.data.data{gp}));
                iscon  = mean(cellfun(@(x) contains(x,'con'),LIMO.data.data{gp}));
            end
        end
        
        if any(isbeta) && any(iscon)
            error('input data mix Beta and con files - not supported')
        elseif sum(isbeta) == 0 && sum(iscon) == 0
            all_files = limo_get_files([],[],[],LIMO.data.data);
            isbeta    = cellfun(@(x) contains(x,'Beta'),all_files);
            iscon     = cellfun(@(x) contains(x,'con'),all_files);
            if any(isbeta) && any(iscon)
                error('input data mix Beta and con files - not supported')
            elseif sum(isbeta) == 0 && sum(iscon) == 0
                error('not all betas or all con?? issue reading files')
            else
                if all(isbeta)
                    a = 'beta';
                elseif all(iscon)
                    a = 'con';
                end
            end
        else
            if all(isbeta)
                a = 'beta';
            elseif all(iscon)
                a = 'con';
            end
        end
    else
        a = questdlg('load several con files per subject or one beta file','ANOVA loading files','con','beta','beta');
    end
    
    N = 0; cell_nb = 1;
    % beta files
    if strcmp(a,'beta')
        for i=1:gp_nb 
            if length(LIMO.data.data) < i
                [Names{cell_nb},Paths{cell_nb},LIMO.data.data{cell_nb}] = limo_get_files([' beta file gp ',num2str(i)]);
            else
                if ischar(LIMO.data.data{i}) % Case for path to the files
                    [Names{cell_nb},Paths{cell_nb},LIMO.data.data{cell_nb}] = limo_get_files([],[],[],LIMO.data.data{i});
                else % Case when all paths are provided
                    [Names{cell_nb},Paths{cell_nb},LIMO.data.data{cell_nb}] = breaklimofiles(LIMO.data.data{i});
                end
            end
            
            if isempty(Names{cell_nb})
                return
            end
            
            if isfield(LIMO.design,'parameters')
                if ~isempty(LIMO.design.parameters)
                    if length(factor_nb) <=2
                        parameters(:,i) = check_files(Names,1,cell2mat(LIMO.design.parameters(i,:)));
                    else % cell of cells
                        all_param = LIMO.design.parameters;
                        while any(cellfun(@iscell,all_param))
                            all_param = [all_param{cellfun(@iscell,all_param)} all_param(~cellfun(@iscell,all_param))];
                        end
                        parameters(:,i) = check_files(Names,1,cell2mat(all_param));
                    end
                else
                    parameters(:,i) = check_files(Names,1);
                end
            else
                parameters(:,i) = check_files(Names,1);
            end
            
            N = N + size(Names{cell_nb},2);
            cell_nb = cell_nb +1;
        end
        
        if size(parameters,1) ~= prod(factor_nb)
            error(['the number of parameter chosen (',num2str(length(parameters)), ...
                ') does not match the total number of levels (',num2str(prod(factor_nb)),')'])
        end
        LIMO.data.data_dir = Paths;
        
    else  % multiple con files
        for i=1:gp_nb
            if isempty(LIMO.data.data) % GUI
                for j=1:length(factor_nb)
                    for k=1:factor_nb(j)
                        [names{k},paths{k},full_names{k}] = limo_get_files([' gp ',num2str(i),' factor ',num2str(j),' level ',num2str(k)]);
                        if isempty(names{k}); disp('no files found - selection aborded'); return; end
                    end
                end
            else
                for j=1:prod(factor_nb)
                    % Case for path to the files
                    if size(LIMO.data.data{i,j},1) == 1
                        [names{j},paths{j},full_names{j}] = limo_get_files([],[],[],LIMO.data.data{i,j});
                        % Case when all paths are provided
                    elseif size(LIMO.data.data{num2str(i),num2str(j)},1) > 1
                        [names{j},paths{j},full_names{j}] = breaklimofiles(LIMO.data.data{i}{i,j});
                    end
                    if isempty(names{j}); disp('no files found - selection aborded'); return; end
                end
            end
            N = N + size(names{cell_nb},2);
            
            Names{cell_nb}          = names{1};
            Paths{cell_nb}          = paths{1};
            LIMO.data.data{cell_nb} = full_names{1};
            for l=2:size(names,2)
                Names{cell_nb}          = [Names{cell_nb} names{l}];
                Paths{cell_nb}          = [Paths{cell_nb} paths{l}];
                LIMO.data.data{cell_nb} = [LIMO.data.data{cell_nb} full_names{l}];
            end
            LIMO.data.data_dir{cell_nb} = Paths{cell_nb};
            check_files(Names,size(Names,2));
            cell_nb         = cell_nb+1;
            parameters(:,i) = 1:prod(factor_nb);
        end
    end
    LIMO.design.parameters = parameters;
    
    % 3rd organize data
    % ---------------------------------------------------------------------
    % match frames, update LIMO
    [first_frame,last_frame,subj_chanlocs,~,LIMO] = match_frames(Paths,LIMO);
    
    % match channels, update LIMO
    LIMO = match_channels(6,analysis_type,LIMO);
    
    % get data for all parameters
    % -----------------------------
    subject_index = 1;
    matrix_index  = 1;
    for h = 1:gp_nb % each group
        nb_subjects(h) = 0;
        for i=1:size(Paths{h},2)
            tmp = load(cell2mat(LIMO.data.data{h}(i)));
            tmp = tmp.(cell2mat(fieldnames(tmp)));
            
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
                
                if strcmpi(LIMO.Type,'Channels') && length(subj_chanlocs(subject_index).chanlocs) == size(tmp,1)
                    matched_data = limo_match_elec(subj_chanlocs(subject_index).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                elseif  strcmpi(LIMO.Type,'Components')
                    matched_data = tmp(:,begins_at:ends_at,:);
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
                    if size(LIMO.design.electrode,2) == 1
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
                        data(:,:,:,:,matrix_index) = matched_data;
                    else
                        data(:,:,:,matrix_index) = matched_data;
                    end
                else
                    if strcmp(LIMO.Analysis,'Time-Frequency') && ...
                            all(size(squeeze(matched_data)) == size(squeeze(data(:,:,:,:,1))))
                        data(1,:,:,:,matrix_index) = matched_data;
                    elseif ~strcmp(LIMO.Analysis,'Time-Frequency') && ...
                            all(size(squeeze(matched_data)) == size(squeeze(data(:,:,:,:,1))))
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
    
    % if con files, re-stack on parameter dimension
    if strcmp(a,'con')
        nb_subjects = single(nb_subjects./prod(factor_nb));
        tmp_data = squeeze(data); clear data; % parameter dim = 1
        for c = prod(factor_nb):-1:1
            sub_idx1 = nb_subjects*c; % since we loaded data in reverse order
            sub_idx2 = sub_idx1-nb_subjects+1;
            if strcmp(LIMO.Analysis,'Time-Frequency')
                data(:,:,:,c,1:nb_subjects) = tmp_data(:,:,:,sub_idx2:sub_idx1);
            else
                data(:,:,c,1:nb_subjects) = tmp_data(:,:,sub_idx2:sub_idx1);
            end
        end
    end
    clear tmp_data
    
    % 4th send relevant info to LIMO_random_robust
    % ----------------------------------------------
    for h=1:gp_nb
        for i=1:size(Paths{h},2)
            if removed{h}(i) == 1
                LIMO.data.data{h}(i) = []; LIMO.data.data_dir{h}(i) = []; % somehow to indicate this subject is removed
            end
        end
    end
    
    % final check
    if sum(nb_subjects) < prod(factor_nb)
       error('there are more variables than observations, some factors can''t be estimated') 
    end

    % the expected dim in LIMO_rep_anova are [frames, subjects,
    % conditions] so we send to LIMO_random_robust data with dim
    % [channels, frames, subjects, conditions] + one vector
    % describing the group belonging
    
    if strcmp(LIMO.Analysis,'Time-Frequency')
        tmp_data = NaN(size(data,1),size(data,2),size(data,3),sum(nb_subjects), prod(factor_nb));
    else
        tmp_data = NaN(size(data,1),size(data,2),sum(nb_subjects), prod(factor_nb));
    end
    gp = NaN(sum(nb_subjects),1);
    
    for i=1:gp_nb
        % select only relevant parameters (could be different from different groups)
        current_param = parameters(:,i);
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
                error('The parameter %g requested (gp %g) is not valid, beta max=%g ',current_param(j),i,size(data,3))
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
        LIMO.design.factor_names = cell(1,length(factor_nb));
        for i=1:length(factor_nb)
            LIMO.design.factor_names{i} = cell2mat(inputdlg(['name factor ' num2str(i) ': ' num2str(factor_nb(i)) ' levels'],'Rep. Measures Names',1,{''},'on'));
            if isempty(LIMO.design.factor_names{i})
                LIMO.design.factor_names{i} = ['Factor_' num2str(i)];
            end
        end
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
N     = size(cellfiles,2);
Paths = cell(1,N);
Names = cell(1,N);
Files = cell(1,N);
for ifiles = 1:N
    [Paths{ifiles}, filename, ext] = fileparts(cellfiles{ifiles});
    Names{ifiles}                  = [filename ext];
    Files{ifiles}                  = fullfile(Paths{ifiles},[filename ext]);
end
end

%% file checking
function parameters = check_files(Names,gp,parameters)
% after selecting file, check they are all the same type (betas or con)
% return parameters that match with files (eg 1 for con, or whatever value
% for the beta file up to the number of regressors in the design matrix

if nargin < 3
    parameters = [];
end
if gp == 1
    if iscell(Names{gp})
        Names = Names{gp};
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
        error('file selection failed, only Beta or Con files are supported')
    elseif (isempty(is_beta)) == 0 && sum(is_beta) == size(Names,2) && nargout ~= 0
        if isempty(parameters)
            param = cell2mat(inputdlg('which parameters to test e.g [1:3]','parameters option'));
            if isempty(param)
                disp('selection aborded'); return
            else
                if contains(param,'[') && contains(param,']')
                    parameters = eval(param);
                else
                    parameters = eval(['[' param ']']);
                end
            end
        end
        
    elseif (isempty(is_con)) == 0 && sum(is_con) == size(Names,2)
        parameters = 1;
    end
    
elseif gp > 1
    
    % several sample cases
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
    
    if (isempty(is_beta)) == 0 && sum(cell2mat(test)) ~= size(Names,2) || (isempty(is_con)) == 0 && sum(cell2mat(test)) ~= size(Names,2)
        error('file selection failed, only Beta or Con files are supported');
    elseif (isempty(is_beta)) == 0 && sum(cell2mat(test)) == size(Names,2) && nargout ~= 0
        if isempty(parameters)
            parameters = eval(cell2mat(inputdlg('which parameter(s) to test e.g 1','parameters option')));
        elseif ~isempty(parameters) && size(parameters,2) ~=1 && size(parameters,2) ~=gp
            fprintf(2,'A valid parameter value must be provided \n');
            return
        end
        if isempty(parameters) || size(parameters,2) ~=1 && size(parameters,2) ~=gp
            return
        end
    elseif (isempty(is_con)) == 0 && sum(cell2mat(test)) == size(Names,2)
        parameters = 1;
    end
end
end

%% frame matching
function [first_frame,last_frame,subj_chanlocs,channeighbstructmat,LIMO] = match_frames(Paths,LIMO)

% once we have all the files, we need to collect information to match the
% frames across subjects
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
        channel = inputdlg('which electrode to analyse [?]','Electrode option');
    elseif isempty(LIMO.design.electrode) && strcmp(LIMO.Type,'Components')
        channel = inputdlg('which component to analyse [?]','Component option'); % can be 1 nb or a vector of channels (channel optimized analysis)
    else
        channel = {num2str(LIMO.design.electrode)};
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
                errordlg(['the nb of ' LIMO.Type ' does not match the number of subjects'],'Error');
                return;
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
        
    elseif size(eval(cell2mat(channel)),2) == 1 || size(eval(cell2mat(channel)),2) == size(Names,2)
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
            LIMO.design.electrode       = eval(cell2mat(channel));
            LIMO.data.chanlocs          = LIMO.data.expected_chanlocs;
            LIMO.data.expected_chanlocs = LIMO.data.expected_chanlocs(LIMO.design.electrode);
        else
            LIMO.design.component = eval(cell2mat(channel));
            LIMO.data.chanlocs    = [];
        end
    else
        errordlg(['the nb of ' LIMO.Type ' does not match the number of subjects'],[LIMO.Type(1:end-1) ' error']);
        return
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

disp('gathering data ...');
if stattest == 1 || stattest == 4
    index = 1;
    for i=size(LIMO.data.data,2):-1:1 % for each subject
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
                        data(1,:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                    else
                        data(1,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                    end
                elseif strcmpi(LIMO.Type,'Components')
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        data(:,:,:,:,index) = tmp(LIMO.design.component,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                    else
                        data(:,:,:,index) = tmp(LIMO.design.component,begins_at:ends_at,:);
                    end
                end
                index = index + 1;
                removed(i) = 0;
                
                % Use multiple single channels
            else
                if strcmpi(LIMO.Type,'Channels')
                    out = LIMO_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        data(1,:,:,:,index) = out(i,:,:,:);
                    else
                        data(1,:,:,index) = out(i,:,:); % matches the expected chanloc of the subject
                    end
                elseif strcmpi(LIMO.Type,'Components')
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        data(1,:,:,:,index) = tmp(LIMO.design.component(i),begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                    else
                        data(1,:,:,index) = tmp(LIMO.design.component(i),begins_at:ends_at,:);
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
    
elseif stattest == 2 % t-tests and gp ANOVA
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
                if size(LIMO.design.electrode,2) == 1
                    if strcmpi(LIMO.Type,'Channels') && length(subj_chanlocs(subject_nb).chanlocs) == size(tmp,1)
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            tmp_data(1,:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                        else
                            tmp_data(1,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                        end
                    elseif strcmpi(LIMO.Type,'Components')
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            tmp_data(1,:,:,:,index) = tmp(LIMO.design.electrode,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:); % all param for beta, if con, adjust dim
                        else
                            tmp_data(1,:,:,index) = tmp(LIMO.design.electrode,begins_at:ends_at,:); % all param for beta, if con, adjust dim
                        end
                    end
                    removed(igp,i) = 0;
                    index = index + 1;
                    
                    % Use multiple single channels
                else
                    if strcmpi(LIMO.Type,'Channels')
                        out = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            tmp_data(1,:,:,:,index) = out(subject_nb,:,:,:); % matches the expected chanloc of the subject
                        else
                            tmp_data(1,:,:,index) = out(subject_nb,:,:);     % matches the expected chanloc of the subject
                        end
                    elseif strcmpi(LIMO.Type,'Components')
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            tmp_data(1,:,:,:,index) = tmp(LIMO.design.electrode(subject_nb),begins_at(1):ends_at(1),begins_at(2):ends_at(2),:); % matches the expected chanloc of the subject
                        else
                            tmp_data(1,:,:,index) = tmp(LIMO.design.electrode(subject_nb),begins_at:ends_at,:);     % matches the expected chanloc of the subject
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
