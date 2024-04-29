function result=limo_central_tendency_and_ci(varargin)

% The function computes estimates of central tendency (mean, trimmed mean,
% Harell-Davis 0.5 decile, median) with 95% Bayesian Highest Density Intervals.
% INPUTS are either a data matrix, con files or LIMO.mat files. 
% If you input LIMO files, estimates of the raw data for the 
% categorical variables will be performed (it makes no sense to summarize continuous
% variables). Non overlap of 95% HDI shows univariate and 'non-corrected' significant
% differences. This can also be assessed directly using limo_plot_difference
%
% FORMAT
% limo_central_tendency_and_ci(varargin)
% result = limo_central_tendency_and_ci(varargin)
%
% INPUTS
% limo_central_tendency_and_ci(expected_chan_loc)
%                expected_chan_loc is the name of the channel structure from EEGLAB
%                this option calls the GUI
%
% limo_central_tendency_and_ci(data, 'Analysis_type',selected_channels,savename)
%                data is a [channel * [freq/time] frames * trials/subjects] matrix
%                Analysis_type should be 'Mean', 'Trimmed mean', 'HD' or 'Median'
%                selected_channels can be [] for all brain or 1 or many channels (1 per trial/subject)
%                savename (optional) name for saving files
%
% limo_central_tendency_and_ci(Files, parameters, expected_chan_loc, 'Estimator1', 'Estimator2', selected_channels,savename)
%                Files are the full names (with paths) of LIMO.mat files or con files
%                parameters are which part of the raw data to analyse based on the design matrix, e.g. [1 2]; 
%                           it can also be 'con_X' (x being the contrast number) meaning that the columns of the design matrix spanned by 
%                           computed contrasts will be used (useful when the design change among subjects)
%                           if Files are contrast files, parameter must be 1              
%                expected_chan_loc is the channel structure from EEGLAB but for the group of subjects
%                Estimator should be 'Mean', 'Weighted mean', 'Trimmed mean', 'HD' or 'Median' (doesn't matter for con files)
%                Estimator 1 is applied to trials within-subjects
%                Estimator 2 is applied across subjects
%                selected_channels can be [] for all brain or 1 or many channels (=nb files) 
%                savename (optional) name for saving files
%
% OUTPUTS result = limo_central_tendency_and_ci()
%         result is a structure with the fields 'subject' and 'central'
%         if not called, the equivalent of the results fields are saved on the drive
%         result.subjects returns the estimator1 computed per subject DIM [channel freq/time parameter subject]
%         result.estimator2 returns the estimator 2 computed across subjects DIM [channel freq/time parameter 3]
%                        the last dim is 3 for low CI bound, estimator value, high CI bound
%                estimator2 van be Median, Harrell_Davis, trimmed_mean, or mean
%
%         if empty files are created on the drive (typically when called via GUI)
%         for single subject the fomat is channel, freq/time, parameters, subjects
%         for the group it's a structure data.estimator2 name and data.limo
%
% Examples Call the GUI
%          ------------
%          limo_central_tendency_and_ci('limo_gp_level_chanlocs.mat')
%
%          Trimmed mean of Beta parameters
%          -------------------------------
%          data = load('Yr.mat'); for instance betas from a rep- measure ANOVAs
%          % being in the ANOVA directory, it will also use LIMO.mat for extra info
%          for condition = 1:9
%              tmp = squeeze(data.Yr(:,:,:,condition));
%              limo_central_tendency_and_ci(tmp, 'Trimmed mean',[],['Condition' num2str(condition)])
%          end
%
%          Weighed means ERP per subject + group level trimmed mean and 95% HDI
%          ---------------------------------------------------------------------
%          Files = fullfile('..derivatives/LIMO_studyname','LIMO_files_face_detection_all_Face_time_GLM_Channels_Time_WLS.txt')
%          expected_chan_loc = fullfile('.../derivatives','limo_gp_level_chanlocs.mat')
%          limo_central_tendency_and_ci(Files, [1 2 3], expected_chan_loc, 'Weighted mean', 'Trimmed mean', [], 'ERPs')
%
%          Trimmed mean ERP + 95% HDI for a condition for a single subject from EEGLAB .daterp
%          ------------------------------------------------------------------------------------
%          data = load('sub-002.daterp','-mat'); % read single trials
%          index = arrayfun(@(x) contains(x.type,'famous','IgnoreCase',true), data.trialinfo); % get condition of intetest
%          FN = fieldnames(data); all_channels = find(contains(FN,'chan'));
%          for channel = length(all_channels):-1:1
%              data_matrix(channel,:,:) = data.(FN{channel})(:,index); % make a data matrix
%          end 
%          limo_central_tendency_and_ci(data_matrix, 'Trimmed mean',[],'Famous_trimmed_mean')
%
% Guillaume Rousselet provided the initial code to do the stats
% Cyril Pernet made the interface, organize to suite EEG data etc - version 1. 18 May 2010
% June/July 2013 - Fixed some bugs CP / thx to Benedikt Ehinger
% Novembre 2013 - fixed further issues related to parameter selection CP / thx to Matt Craddock
% version 2 September 2015 - included within subject weighted mean + update for time frequency
% version 3 February 2016 - CP/GAR updated for Bayesian HDI
% version 4+ CP maintenance of inputs/arguments/beta or con/etc .. see gitlog
%
% see also limo_central_estimator.m limo_add_plots.m limo_plot_difference.m
% -------------------------------------------------------------------------
%  Copyright (C) LIMO Team 2021


%% file selection and checkings
% -----------------------------
current_dir = pwd; warning off
result = []; % the output if requested
data   = []; % the matrix of data to compute summary stats on

if nargin == 3 || nargin == 4
    % ------------------------
    
    data = varargin{1};
    if ischar(data)
        if ~exist(data,'file')
            limo_errordlg('%s does not exist',data); 
            return
        else
            try
                data = load(data);
                data = data.(cell2mat(fieldnames(data)));
            catch
                limo_errordlg('%s is not a matrix',data); 
                return
            end
        end
    end
    
    if ndims(data)<3 || ndims(data) >4 %#ok<*ISMAT>
        if ndims(data) == 2
            disp('for 2D data, try using limo_central_estimator.m');
        end
        limo_errordlg('data in must be 3 or 4 dimensional: [1/all channels], [freq/time] frames, subjects'); 
        return
    elseif ndims(data) == 4
        limo.Analysis = 'Time-Frequency';
         if exist('LIMO.mat','file')
            disp('updating data structure with local LIMO.mat')
            LIMO                    = load('LIMO.mat');
            limo.Level              = LIMO.LIMO.Level;
            limo.Type               = LIMO.LIMO.Type;
            limo.data.sampling_rate = LIMO.LIMO.data.sampling_rate;
            limo.data.trim1         = LIMO.LIMO.data.trim1;
            limo.data.trim2         = LIMO.LIMO.data.trim2;
            limo.data.start         = LIMO.LIMO.data.start;
            limo.data.end           = LIMO.LIMO.data.end;
            limo.data.trim_lowf     = LIMO.LIMO.data.trim_lowf;
            limo.data.trim_highf    = LIMO.LIMO.data.trim_highf;
            limo.data.lowf          = LIMO.LIMO.data.lowf;
            limo.data.highf         = LIMO.LIMO.data.highf;
            limo.data.tf_times      = LIMO.LIMO.data.tf_times;
            limo.data.tf_freqs      = LIMO.LIMO.data.tf_freqs;
            if isfield(LIMO.LIMO.data, 'neighbouring_matrix')
                limo.data.neighbouring_matrix = LIMO.LIMO.data.neighbouring_matrix;
            end
            if isfield(LIMO.LIMO.data, 'expected_chanlocs')
                limo.data.expected_chanlocs = LIMO.LIMO.data.expected_chanlocs;
            end
            if isfield(LIMO.LIMO.data, 'chanlocs')
                limo.data.expected_chanlocs = LIMO.LIMO.data.chanlocs;
            end
         end
    else
        if exist('LIMO.mat','file')
            disp('updating data structure with local LIMO.mat')
            LIMO                    = load('LIMO.mat');
            limo.Level              = LIMO.LIMO.Level;
            limo.Analysis           = LIMO.LIMO.Analysis;
            limo.Type               = LIMO.LIMO.Type;
            limo.data.sampling_rate = LIMO.LIMO.data.sampling_rate;
            limo.data.trim1         = LIMO.LIMO.data.trim1;
            limo.data.trim2         = LIMO.LIMO.data.trim2;
            limo.data.start         = LIMO.LIMO.data.start;
            limo.data.end           = LIMO.LIMO.data.end;
            if isfield(LIMO.LIMO.data, 'timevect')
                limo.data.timevect = LIMO.LIMO.data.timevect;
            end
            if isfield(LIMO.LIMO.data, 'freqlist')
                limo.data.expected_chanlocs = LIMO.LIMO.data.freqlist;
            end
            if isfield(LIMO.LIMO.data, 'neighbouring_matrix')
                limo.data.neighbouring_matrix = LIMO.LIMO.data.neighbouring_matrix;
            end
            if isfield(LIMO.LIMO.data, 'expected_chanlocs')
                limo.data.expected_chanlocs = LIMO.LIMO.data.expected_chanlocs;
            end
            if isfield(LIMO.LIMO.data, 'chanlocs')
                limo.data.expected_chanlocs = LIMO.LIMO.data.chanlocs;
            end
        else
            limo.Analysis = 'Time or Frequency';
        end
    end
    
    Estimator2 = varargin{2};
    if strcmpi(Estimator2,'Trimmed mean') || strcmpi(Estimator2,'HD') ...
            || strcmpi(Estimator2,'Median') || strcmpi(Estimator2,'Mean') ...
            || strcmpi(Estimator2,'All')
        parameters = 1; %#ok<NASGU>
    else
        limo_errordlg('type of estimator not recognized'); 
        return
    end
    
    selected_channels = varargin{3};
    if ~isempty(selected_channels)
        Analysis_type = '1 channel only';
        if strcmpi(limo.Analysis,'Time-Frequency')
            data = data(selected_channels,:,:,:);
        else
            data = data(selected_channels,:,:);
        end
    else
        Analysis_type = 'Full brain analysis';
    end
          
    if nargin == 4
        savename = varargin{4};
        [p,f,ext]=fileparts(savename);
        if strcmp(ext,'.mat')
            savename=fullfile(p,f);
        end
    end
    
elseif nargin == 6 || nargin == 7
    % ---------------------------
    
    if exist(varargin{1},'file')
        Files = varargin{1};
        if size(Files,1) == 1 % select a txt file listing all files
            [Names,Paths,Files] = limo_get_files([],[],[],Files);
        end
    else
        limo_errordlg('input file not found'); 
        return
    end
    
    parameters = varargin{2};
    is_limo    = zeros(1,size(Names,2));
    is_con     = zeros(1,size(Names,2));
    for i=size(Names,2):-1:1
        if strfind(Names{i},'LIMO'); is_limo(i) = 1;
        elseif strfind(Names{i},'con'); is_con(i) = 1; end
    end
    if all(is_con) && parameters ~=1; parameters = 1; 
        warning on; warning('all con files in, parameter set to 1'); warning off
    end

    if exist(fullfile(pwd,'LIMO.mat'),'file')
        disp('updating data structure with local LIMO.mat')
        LIMO                    = load('LIMO.mat');
        limo.Level              = LIMO.LIMO.Level;
        limo.Analysis           = LIMO.LIMO.Analysis;
        limo.data.sampling_rate = LIMO.LIMO.data.sampling_rate;
        limo.data.trim1         = LIMO.LIMO.data.trim1;
        limo.data.trim2         = LIMO.LIMO.data.trim2;
        limo.data.start         = LIMO.LIMO.data.start;
        limo.data.end           = LIMO.LIMO.data.end;
        if isfield(LIMO.LIMO.data, 'timevect')
            limo.data.timevect = LIMO.LIMO.data.timevect;
        end
        if isfield(LIMO.LIMO.data, 'freqlist')
            limo.data.expected_chanlocs = LIMO.LIMO.data.freqlist;
        end
        if isfield(LIMO.LIMO.data, 'neighbouring_matrix')
            limo.data.neighbouring_matrix = LIMO.LIMO.data.neighbouring_matrix;
        end
        if isfield(LIMO.LIMO.data, 'expected_chanlocs')
            limo.data.expected_chanlocs = LIMO.LIMO.data.expected_chanlocs;
        end
        if isfield(LIMO.LIMO.data, 'chanlocs')
            limo.data.expected_chanlocs = LIMO.LIMO.data.chanlocs;
        end
    else
        limo.Analysis = 'Time or Frequency';
    end

    expected_chanlocs   = varargin{3};
    if ischar(expected_chanlocs)
        expected_chanlocs             = load(expected_chanlocs);
        limo.data.neighbouring_matrix = expected_chanlocs.channeighbstructmat;
        limo.data.expected_chanlocs   = expected_chanlocs.expected_chanlocs;
        expected_chanlocs             = limo.data.expected_chanlocs;
    end
    Estimator1          = varargin{4};
    Estimator2          = varargin{5};
    selected_channels   = varargin{6};
    if isempty(selected_channels)
        Analysis_type = 'Full brain analysis';
    else
        Analysis_type = '1 channel only';
        expected_chanlocs = expected_chanlocs(selected_channels);
    end
    
    % match frames
    % -------------
    [first_frame,last_frame,subj_chanlocs,limo] = limo_match_frames(Paths,limo);
    
    % get data for all parameters dim [channel, frame, param, nb subjects
    % ---------------------------------------------------------------------
    disp('gathering data ...'); 
    for i=size(Paths,2):-1:1 % for each subject
        fprintf('processing subject %g\n',i);
        LIMO         = load(fullfile(Paths{i},'LIMO.mat')); 
        LIMO         = LIMO.LIMO;
        limo.Type{i} = LIMO.Type;

        if all(is_limo)
            Yr = load(fullfile(Paths{i},'Yr.mat'));   
        elseif all(is_con)
            Yr = load(Files{i}); 
        end
        Yr = Yr.(cell2mat(fieldnames(Yr)));

        if strcmpi(LIMO.Analysis,'Time-Frequency')
            begins_at = fliplr((max(first_frame) - first_frame(i,:) + 1)); % returns time/freq/or freq-time
            ends_at(1) = size(Yr,2) - (last_frame(i,2) - min(last_frame(:,2)));
            ends_at(2) = size(Yr,3) - (last_frame(i,1) - min(last_frame(:,1)));
        else
            begins_at = max(first_frame) - first_frame(i) + 1;
            ends_at = size(Yr,2) - (last_frame(i) - min(last_frame));
        end
        
        if max(parameters) <= sum(LIMO.design.nb_conditions+LIMO.design.nb_interactions) || ...
                max(parameters) == size(LIMO.design.X,2) || ...% any categorial or the constant
                (~isnumeric(parameters) && contains(parameters,'con'))
            if all(is_limo)
                if isnumeric(parameters)
                    index = logical(sum(LIMO.design.X(:,parameters)==1,2));
                else
                    if contains(parameters,'con')
                        tmp=find(LIMO.contrast{str2double(parameters(5:end))}.C);
                        index = logical(sum(LIMO.design.X(:,tmp)==1,2)); clear tmp
                    else
                        limo_error('unrecognized input parameter')
                    end
                end

                for channel=size(Yr,1):-1:1
                    if strcmpi(Estimator1,'Weighted Mean')
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            for f=size(Yr,2):-1:1
                                fw(1,f,:,:)  = squeeze(Yr(channel,f,:,index)).*repmat(squeeze(LIMO.design.weights(channel,f,index))',size(Yr,3),1);
                            end
                            tmp(channel,:,:) = limo_tf_4d_reshape(fw,LIMO.data.size3D);
                            clear fw;
                        else
                            tmp(channel,:,:) = squeeze(Yr(channel,:,index)).*repmat(LIMO.design.weights(channel,index),size(Yr,2),1);
                        end
                    else
                        tmp(channel,:,index) = squeeze(Yr(channel,:,index));
                    end
                end
                
                % 1st level analysis
                % --------------------
                if strcmpi(Estimator1,'Trimmed mean') % trim raw data @ 20%
                    tmp = limo_trimmed_mean(tmp,20);
                elseif strcmpi(Estimator1,'Median') % median raw data
                    tmp = nanmedian(tmp,3);
                elseif strcmpi(Estimator1,'HD') % mid-decile Harrell-Davis of raw data
                    tmp = limo_harrell_davis(tmp,0.5);
                elseif strcmpi(Estimator1,'Mean') || strcmpi(Estimator1,'Weighted Mean') % mean of raw data
                    tmp = nanmean(tmp,3);
                end
            else
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    tmp = limo_tf_4d_reshape(Yr,LIMO.data.size3D);
                    tmp = squeeze(tmp(:,:,1));
                else
                    tmp = Yr;
                    tmp = squeeze(tmp(:,:,1));
                end
            end
            clear Yr
            
            if strcmpi(Analysis_type,'Full brain analysis') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
                 if strcmpi(LIMO.Analysis,'Time-Frequency')
                     data(:,:,:,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,reshape(tmp,LIMO.data.size4D(1:3)));
                 else
                     data(:,:,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                 end
            elseif strcmpi(Analysis_type,'1 channel only') && length(subj_chanlocs(i).chanlocs) == size(tmp,1)
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    if size(selected_channels,2) == 1
                        data(1,:,:,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,reshape(tmp,LIMO.data.size4D(1:3)));
                    else
                        out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,reshape(tmp,LIMO.data.size4D(1:3))); 
                        data(1,:,:,i) = out(i,:,:); 
                    end
                else
                    if size(selected_channels,2) == 1
                        data(1,:,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                    else
                        out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                        data(1,:,i) = out(i,:,:); % matches the expected chanloc of the subject
                    end
                end
            end
            clear tmp
        else
            if max(parameters) > size(LIMO.design.X,2)
                warning('subject %g, parameter %g not computed: \n the design only includes %g regressors plus the constant',Paths{i},size(LIMO.design.X,2));
            else
                warning('subject %g, \n parameter %g not computed - continuous regressor',Paths{i},max(parameters));
            end
        end
    end
    
    if nargout == 1
        result.subjects = data;
    end
    
    if nargin == 7
        savename = varargin{7};
        [p,f,ext]=fileparts(savename);
        if strcmp(ext,'.mat')
            savename=fullfile(p,f);
        end
    end

    if all(cellfun(@(x) strcmpi(x,limo.Type{1}), limo.Type))
        limo.Type = limo.Type{1};
    else
        limo_errordlg('despite successful data aggregation, LIMO.Type differ?? channels/compomnents/sources - check your data')
        return
    end
    
elseif nargin == 1
    % ---------------------------
    
    % Expected_chanlocs
    expected_chanlocs = load(varargin{1});
    
    % check if Betas/Con 
    option = limo_questdlg('type of analysis','what data to analyse?','Raw Data','Betas','Con','Betas');
    if isempty(option)
        return
    end
    
    % -----------------------------
    % ANALYSIS ON BETAS PARAMETERS
    % -----------------------------
    if strcmpi(option,'Betas') || strcmpi(option,'Con')
        
        Estimator1 = option;
        Estimator2 = limo_questdlg('Estimation option','which estimator?','Mean','Trimmed mean','HD/Median','Trimmed mean');
        if strcmpi(Estimator2,'HD/Median')
            Estimator2 = 'HD';
        end
        
        % get the data
        % ------------
        Names = {}; %#ok<NASGU>
        [Names,Paths,Files] = limo_get_files([],{'*.mat;*.txt','matlab or text'},sprintf('Select %s files',option)); %#ok<ASGLU>
        if isempty(Names)
            return
        elseif size(Names,2) < 3
            limo_errordlg('LIMO cannot do group bootrap estimates - too few subjects')
            return
        end
        
        % check type of files and returns which beta param to test
        % -------------------------------------------------------
        is_betas = [];
        is_con   = [];
        for i=size(Names,2):-1:1
            if strfind(Names{i},'Betas')
                is_betas(i) = 1;
            elseif strfind(Names{i},'con')
                is_con(i) = 1;
            end
        end
        
        if (isempty(is_betas)) == 0 && sum(is_betas) == size(Names,2)
            if strcmpi(Estimator1,'Con')
                limo_warndlg('you indicated computation for contrasts, but all files are beta parameters - still computing though',...
                    'selection warning');
                Estimator1 = 'Betas';
            end
            parameters = limo_inputdlg('which parameters to test e.g [1:3]','parameters option');
            if isempty(parameters)
                return
            else
                parameters = cell2mat(parameters);
                if ~strcmp(parameters(1),'[') && ~strcmp(parameters(end),']')
                    parameters = ['[' parameters ']'];
                end
                parameters = eval(parameters);
            end
        elseif (isempty(is_con)) == 0 && sum(is_con) == size(Names,2)
            if strcmpi(Estimator1,'Betas')
                limo_warndlg('you indicated computation for Betas, but all files are contrasts - still computing though',...
                    'selection warning')
                Estimator1 = 'Con';
            end
            parameters = 1; 
        else
            limo_errordlg('file selection failed, only Betas.mat files are supported'); 
            return
        end
        
        % match frames
        % ------------
        limo.data.neighbouring_matrix               = expected_chanlocs.channeighbstructmat;
        limo.data.expected_chanlocs                 = expected_chanlocs.expected_chanlocs;
        [first_frame,last_frame,subj_chanlocs,limo] = limo_match_frames(Paths,limo);
        limo.Level                                  = 2;
        
        % match channels
        % --------------
        Analysis_type   = limo_questdlg('Rdx option','type of analysis?','Full brain analysis','1 channel only','Full brain analysis');
        if isempty(Analysis_type)
            return
        end
        
        if strcmpi(Analysis_type,'1 channel only')
            channel = limo_inputdlg('which channel to analyse [?]','channel option'); % can be 1 nb or a vector of channels (channel optimized analysis)
            if isempty(cell2mat(channel))
                [file,dirf,index] = uigetfile('*.mat','select your channel file');
                if index == 0
                    return
                else
                    channel_vector = load(fullfile(dirf,file));
                    channel_vector = channel_vector.cell2mat(fieldname(channel_vector));
                    % check the vector has the same length as the number of files
                    if length(channel_vector) ~= size(Names,2)
                        errordlg('the nb of channels does not match the number of subjects','channel error'); return;
                    end
                    % restric the channels
                    expected_chanlocs = limo.data.expected_chanlocs(channel_vector);
                end
            elseif size(eval(cell2mat(channel)),2) == 1 || size(eval(cell2mat(channel)),2) == size(Names,2)
                selected_channels = eval(cell2mat(channel));
                expected_chanlocs = limo.data.expected_chanlocs(selected_channels);
            else
                limo_errordlg('the nb of channels does not match the number of subjects','channel error'); 
                return
            end
        else
            expected_chanlocs = limo.data.expected_chanlocs;
        end
        
        % make one large matrix
        disp('gathering data ...'); index = 1;
        for i=size(Paths,2):-1:1 % for each subject
            fprintf('processing subject %g\n',i);
            % load file and store contend
            LIMO = load([Paths{i} filesep 'LIMO.mat']); LIMO = LIMO.LIMO;
            Yr   = load([Paths{i} filesep Names{i}]);
            Yr   = Yr.(cell2mat(fieldnames(Yr)));
            
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                begins_at  = fliplr((max(first_frame) - first_frame(i,:) + 1)); % returns time/freq/or freq-time
                ends_at(1) = size(Yr,2) - (last_frame(i,2) - min(last_frame(:,2)));
                ends_at(2) = size(Yr,3) - (last_frame(i,1) - min(last_frame(:,1)));
            else
                begins_at = max(first_frame) - first_frame(i) + 1;
                ends_at   = size(Yr,2) - (last_frame(i) - min(last_frame));
            end
            
            if strcmpi(Analysis_type,'Full brain analysis')
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    data(:,:,:,:,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,squeeze(Yr(:,:,:,parameters)));
                else
                    data(:,:,i)     = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,squeeze(Yr(:,:,parameters)));
                end
            elseif strcmpi(Analysis_type,'1 channel only')
                if size(selected_channels,2) == 1
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        data(1,:,:,1:length(parameters),i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,squeeze(Yr(:,:,:,parameters)));
                    else
                        data(1,:,1:length(parameters),i)   = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,squeeze(Yr(:,:,parameters)));
                    end
                else % optimized channel
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        out             = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,squeeze(Yr(:,:,:,parameters))); % out is for all expected chanlocs, i.e. across subjects
                        data(1,:,:,:,i) = out(i,:,:,:); % matches the expected chanloc of the subject
                    else
                        out           = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,squeeze(Yr(:,:,parameters))); % out is for all expected chanlocs, i.e. across subjects
                        data(1,:,:,i) = out(i,:,:); % matches the expected chanloc of the subject
                    end
                end
            end
            clear tmp
        end
        
    else
        % -------------------
        % ANALYSIS ON ERPs
        % -------------------
        
        % select data
        % -----------
        [Names,Paths,Files] = limo_get_files([],{'*.mat;*.txt','matlab or text'},'Select LIMO files'); %#ok<ASGLU>
        if isempty(Names)
            return
        elseif size(Names,2) < 3
            limo_errordlg('LIMO cannot do group bootrap estimates - too few subjects')
            return
        end
        
        % check it's LIMO.mat files and which param to test
        % --------------------------------------------------
        is_limo = [];
        for i=size(Names,2):-1:1
            if strcmpi(Names{i},'LIMO.mat')
                is_limo(i) = 1;
            end
        end
        
        if (isempty(is_limo)) == 0 && sum(is_limo) == size(Names,2)
            Q = limo_questdlg('Type of merging','Options','Evaluate single conditions','Pool Conditions','Evaluate single conditions');
            if strcmpi(Q,'Evaluate single conditions')
                parameters = limo_inputdlg('which parameters to test e.g [1:3]','parameters option');
            else
                parameters = limo_inputdlg('which parameters to pool e.g [1 3 5]','parameters option');
            end
            
            if isempty(parameters)
                return
            else
                parameters = eval(cell2mat(parameters));
                if isnan(parameters)
                    parameters = str2double(cell2mat(parameters));
                end
            end
        
        else
            limo_errordlg('file selection failed, only LIMO.mat files are supported'); 
            return
        end
        
        % check what type of analysis
        % ---------------------------
        Analysis_type   = limo_questdlg('Rdx option','type of analysis?','Full brain analysis','1 channel only','Full brain analysis');
        if isempty(Analysis_type)
            return;
        else
            limo.Type  = 'Channel';
        end
        
        limo.data.neighbouring_matrix  = expected_chanlocs.channeighbstructmat;
        if strcmpi(Analysis_type,'1 channel only')
           channel = limo_inputdlg('which channel to analyse [?]','channel option'); % can be 1 nb or a vector of channels (channel optimized analysis)
            if isempty(cell2mat(channel))
                [file,dir,index] = uigetfile('*.mat','select your channel file');
                if isempty(file)
                    return
                else
                    cd(dir); 
                    channel_vector = load(file);
                    channel_vector = channel_vector.getfield(channel_vector);
                    % check the vector has the same length as the number of files
                    if length(channel_vector) ~= length(Paths)
                        errordlg('the nb of channels does not match the number of subjects','channel error'); return;
                    end
                    selected_channels = channel_vector;
                    expected_chanlocs = expected_chanlocs.expected_chanlocs(selected_channels);
                end
            elseif size(eval(cell2mat(channel)),2) == 1 || size(eval(cell2mat(channel)),2) == size(Names,2)
                selected_channels = eval(cell2mat(channel));
                expected_chanlocs = expected_chanlocs.expected_chanlocs(selected_channels);
            else
                limo_errordlg('the nb of channels does not match the number of subjects','channel error'); 
                return;
            end
        else
            selected_channels = [];
            expected_chanlocs = expected_chanlocs.expected_chanlocs;
        end
        limo.data.expected_chanlocs = expected_chanlocs;
        
        % select method
        % -------------
        [Estimator1,Estimator2]  = limo_central_tendency_questdlg;
        if isempty(Estimator1) && isempty(Estimator2)
            return
        end
        
        if strcmpi(Estimator1,'All') || strcmpi(Estimator1,'Mean')
            weighted_mean = limo_questdlg('do you want to use weights to compute means?','saving option','yes','no','yes');
        end
        
        % match frames
        % -------------
        [first_frame,last_frame,subj_chanlocs,limo] = limo_match_frames(Paths,limo);
        limo.Level = 2;
        
        % get data for all parameters dim [channel, frame, param, nb subjects
        % ---------------------------------------------------------------------
        disp('gathering data ...'); 
        for i=size(Paths,2):-1:1 % for each subject
            fprintf('processing subject %g',i); disp(' ')
            LIMO = load(fullfile(Paths{i},'LIMO.mat')); LIMO = LIMO.LIMO;
            Yr   = load(fullfile(Paths{i},'Yr.mat'));   Yr = Yr.Yr;
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                begins_at = fliplr((max(first_frame) - first_frame(i,:) + 1)); % returns time/freq/or freq-time
                ends_at(1) = size(Yr,2) - (last_frame(i,2) - min(last_frame(:,2)));
                ends_at(2) = size(Yr,3) - (last_frame(i,1) - min(last_frame(:,1)));
            else
                
                begins_at = max(first_frame) - first_frame(i) + 1;
                ends_at = size(Yr,2) - (last_frame(i) - min(last_frame));
            end
            
            if strcmpi(Q,'Evaluate single conditions')
                for j=length(parameters):-1:1
                    if parameters(j) <= sum(LIMO.design.nb_conditions+LIMO.design.nb_interactions) || ...
                            parameters(j) == size(LIMO.design.X,2)
                        index = LIMO.design.X(:,parameters(j))==1;
                        if strcmpi(weighted_mean,'yes')
                            for channel=1:size(Yr,1)
                                if strcmpi(LIMO.Analysis,'Time-Frequency')
                                    for f=size(Yr,2):-1:1
                                        fw(1,f,:,:) = squeeze(Yr(channel,f,:,index)).*repmat(squeeze(LIMO.design.weights(channel,f,index))',size(Yr,3),1);
                                    end
                                    tmp(channel,:,:) = limo_tf_4d_reshape(fw,LIMO.data.size3D);
                                    clear fw;
                                else
                                    tmp(channel,:,:) = squeeze(Yr(channel,:,index)).*repmat(LIMO.design.weights(channel,index),size(Yr,2),1);
                                end
                            end
                        else
                            tmp = squeeze(Yr(:,:,index)); % retain those trials only
                        end
                        
                        % 1st level analysis
                        % --------------------
                        if strcmpi(Estimator1,'Trimmed mean') % trim raw data @ 20%
                            tmp = limo_trimmed_mean(tmp,20);
                        elseif strcmpi(Estimator1,'Median') % median raw data
                            tmp = nanmedian(tmp,3);
                        elseif strcmpi(Estimator1,'HD') % mid-decile Harrell-Davis of raw data
                            tmp = limo_harrell_davis(tmp,0.5);
                        elseif strcmpi(Estimator1,'Mean') % mean of raw or weighted data
                            tmp = nanmean(tmp,3);
                        end
                        
                        if strcmpi(Analysis_type,'Full brain analysis') && length(subj_chanlocs(i).chanlocs) == size(tmp,1)
                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                data(:,:,:,j,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,reshape(tmp,LIMO.data.size4D(1:3)));
                            else
                                data(:,:,j,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                            end
                        elseif strcmpi(Analysis_type,'1 channel only') && length(subj_chanlocs(i).chanlocs) == size(tmp,1)
                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                if size(selected_channels,2) == 1
                                    data(1,:,:,j,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,reshape(tmp,LIMO.data.size4D(1:3)));
                                else
                                    out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,reshape(tmp,LIMO.data.size4D(1:3))); 
                                    data(1,:,:,j,i) = out(i,:,:); 
                                end
                            else
                                if size(selected_channels,2) == 1
                                    data(1,:,j,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                                else
                                    out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                                    data(1,:,j,i) = out(i,:,:); % matches the expected chanloc of the subject
                                end
                            end
                        end
                        clear tmp
                    else
                        if max(j) > size(LIMO.design.X,2)
                            warning('subject %g, parameter %g not computed: \n the design only includes %g regressors plus the constant',Paths{i},size(LIMO.design.X,2));
                        else
                            warning('subject %g, \n parameter %g not computed - continuous regressor',Paths{i},j);
                        end
                    end
                end
            elseif strcmpi(Q,'Pool Conditions')
                if max(parameters) <= sum(LIMO.design.nb_conditions)+sum(LIMO.design.nb_interactions) || ...
                        max(parameters) == size(LIMO.design.X)
                    index = find(sum(LIMO.design.X(:,parameters)==1,2)); % find all trials from selected columns
                    if strcmpi(weighted_mean,'yes')
                        for channel=size(Yr,1):-1:1
                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                for f=size(Yr,2):-1:1
                                     fw(1,f,:,:) = squeeze(Yr(channel,f,:,index)).*repmat(squeeze(LIMO.design.weights(channel,f,index))',size(Yr,3),1);
                                end
                                tmp(channel,:,:) = limo_tf_4d_reshape(fw,LIMO.data.size3D);
                                clear fw;
                            else
                                tmp(channel,:,:) = squeeze(Yr(channel,:,index)).*repmat(LIMO.design.weights(channel,index),size(Yr,2),1);
                            end
                        end
                    else
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            tmp =  limo_tf_4d_reshape(squeeze(Yr(:,:,:,index)),LIMO.data.size3D);
                        else
                            tmp =  squeeze(Yr(:,:,index)); % retain those trials only
                        end
                    end
                    
                    % 1st level analysis
                    % --------------------
                    if strcmpi(Estimator1,'Trimmed mean') % trim raw data @ 20%
                        tmp=limo_trimmed_mean(tmp,20);
                    elseif strcmpi(Estimator1,'Median') % median raw data
                        tmp = nanmedian(tmp,3);
                    elseif strcmpi(Estimator1,'HD') % mid-decile Harrell-Davis of raw data
                        tmp = limo_harrell_davis(tmp,0.5);
                    elseif strcmpi(Estimator1,'Mean') % mean of raw data on which we do across subjects TM, HD and Median
                        tmp = nanmean(tmp,3);
                    end
                    
                    if strcmpi(Analysis_type,'Full brain analysis') && length(subj_chanlocs(i).chanlocs) == size(tmp,1)
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            data(:,:,:,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,reshape(tmp,LIMO.data.size4D(1:3)));
                        else
                            data(:,:,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                        end
                    elseif strcmpi(Analysis_type,'1 channel only') && length(subj_chanlocs(i).chanlocs) == size(tmp,1)
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            if size(selected_channels,2) == 1
                                data(1,:,:,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,reshape(tmp,LIMO.data.size4D(1:3)));
                            else
                                out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,reshape(tmp,LIMO.data.size4D(1:3))); % out is for all expected chanlocs, ie across subjects
                                data(1,:,:,i) = out(i,:,:,:); % matches the expected chanloc of the subject
                            end
                        else
                            if size(selected_channels,2) == 1
                                data(1,:,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                            else
                                out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                                data(1,:,i) = out(i,:,:); % matches the expected chanloc of the subject
                            end
                        end
                    end
                    clear tmp
                else
                   fprintf('pooling not computed - one or more continuous regressor selected \n');
                end
            end
            clear Yr
        end
        % update estimator1 name
        if strcmpi(weighted_mean,'yes')
            Estimator1 = 'Weighted mean';
        end
    end
    
else
    limo_errordlg('nb of arguments incorrect');
    return
end % closes varargin


%% Analysis part
% --------------
cd(current_dir)

if ~isempty(data)
    % Data is either [channel, frame, trials/subject] or [channel,
    % frame, conditions (from parameters), subjects] but we always want 4D
    % or 5D data with 1 or more conditions
    if ~strcmpi(limo.Analysis,'Time-Frequency') && ndims(data) == 3
        tmp = data; clear data
        for i=size(tmp,1):-1:1
            for j=size(tmp,2):-1:1
                data(i,j,1,:) = tmp(i,j,:); % now data is 4D
            end
        end
    elseif strcmpi(limo.Analysis,'Time-Frequency') && ndims(data) == 4
        tmp = data; clear data
        for i=size(tmp,1):-1:1
            for j=size(tmp,2):-1:1
                for k=size(tmp,3):-1:1
                    data(i,j,k,1,:) = tmp(i,j,k,:); % now data is 5D
                end
            end
        end
    end
    
    n = size(data,ndims(data)); % number of subjects always last
    if ndims(data) < 4 
        limo_errordlg('an unexpected issue occured, the number of dimensions is too low, likely caused by selected only 1 subject')
        return
    elseif n < 3
        limo_errordlg('LIMO cannot do group bootrap estimates - too few subjects')
        return
    end
    
    if n<=10 && strcmpi(Estimator2,'HD')
        msgbox('CI of the Harell Davis estimates cannot be computed for less than 11 observations - switched to median','Computation info');
        Estimator2 = 'Median';
    end
    
    % save as
    if nargout ==0
        if exist('savename','var')
            name = savename;
        else
            name = cell2mat(limo_inputdlg('save as [?]','name option'));
            if isempty(name)
                disp('no name selected - aborded'); return
            end
        end
        
        if exist('Estimator1','var')
            newname = sprintf('%s_single_subjects_%s',name,Estimator1);
            if ~strcmpi(limo.Analysis,'Time-Frequency') 
                Data.data = data; Data.limo = limo;
                save (newname,'Data'); clear Data
            elseif strcmpi(limo.Analysis,'Time-Frequency')
                Data.data = data; Data.limo = limo;
                save (newname,'Data'); clear Data
            end
        end
        
    else
        result.subjects = data;
        if exist('limo','var')
            result.limo = limo;
        end
    end
    
    disp('processing data across subjects ..')
    % --------------------------------------------------------------
    if nargout == 1 && exist('limo','var')
        result.limo = limo;
    end
        
    if strcmpi(Estimator2,'Mean') || strcmpi(Estimator2,'All')
        disp('Compute the Mean estimator and 95% CI ...')
        index = 1; h = waitbar(0,'computing','name','% done');
        if strcmpi(limo.Analysis,'Time-Frequency')
            M = NaN(size(data,1),size(data,2),size(data,3),size(data,4),3);
            for k = 1:size(data,4)
                for channel =1:size(data,1)
                    waitbar(index/(size(data,4)*size(data,1)));
                    index              = index+1;
                    if  strcmpi(Analysis_type,'1 channel only')
                        for f=size(data,2):-1:1
                            tmp(1,f,:,:) = data(1,f,:,k,:);
                        end
                        tmp            = limo_tf_4d_reshape(tmp,...
                            [size(data,1) size(data,2)*size(data,3) size(data,5)]);
                    else
                        tmp            = limo_tf_4d_reshape(squeeze(data(:,:,:,k,:)),...
                            [size(data,1) size(data,2)*size(data,3) size(data,5)]);
                    end
                    tmp                = squeeze(tmp(channel,:,:));
                    Y                  = tmp(:,~isnan(tmp(1,:)));
                    [est,ci]           = limo_central_estimator(Y,'mean');
                    M(channel,:,:,k,1) = reshape(ci(1,:),size(data,2),size(data,3));
                    M(channel,:,:,k,2) = reshape(est,size(data,2),size(data,3));
                    M(channel,:,:,k,3) = reshape(ci(2,:),size(data,2),size(data,3));
                end
            end
        else
            M = NaN(size(data,1),size(data,2),size(data,3),3);
            for k = 1:size(data,3)
                for channel =1:size(data,1)
                    waitbar(index/(size(data,3)*size(data,1)));
                    index            = index+1;
                    tmp              = squeeze(data(channel,:,k,:));
                    Y                = tmp(:,~isnan(tmp(1,:)));
                    [est,ci]         = limo_central_estimator(Y,'mean');
                    M(channel,:,k,1) = ci(1,:);
                    M(channel,:,k,2) = est;
                    M(channel,:,k,3) = ci(2,:);
                end
            end
        end
        close(h);
        
        if nargout ==0
            if nargin == 3 || nargin == 4
                newname = sprintf('%s_Mean',name);
            else
                newname = sprintf('%s_Mean_of_%s',name,Estimator1);
            end
            Data.mean = M;
            if exist('limo','var')
                Data.limo = limo;
            end
            save (newname,'Data');
        else
            result.mean = M;
        end
    end
    
    % --------------------------------------------------------------
    if strcmpi(Estimator2,'Trimmed mean') || strcmpi(Estimator2,'All')
        disp('Compute 20% Trimmed Mean estimator and 95% CI ...')
        index = 1; h = waitbar(0,'computing','name','% done');
        if strcmpi(limo.Analysis,'Time-Frequency')
            TM = NaN(size(data,1),size(data,2),size(data,3),size(data,4),3);
            for k = 1:size(data,4)
                for channel =1:size(data,1)
                    waitbar(index/(size(data,4)*size(data,1)));
                    index              = index+1;
                    if  strcmpi(Analysis_type,'1 channel only')
                        for f=size(data,2):-1:1
                            tmp(1,f,:,:) = data(1,f,:,k,:);
                        end
                        tmp            = limo_tf_4d_reshape(tmp,...
                            [size(data,1) size(data,2)*size(data,3) size(data,5)]);
                    else
                        tmp            = limo_tf_4d_reshape(squeeze(data(:,:,:,k,:)),...
                            [size(data,1) size(data,2)*size(data,3) size(data,5)]);
                    end
                    tmp                = squeeze(tmp(channel,:,:));
                    Y                   = tmp(:,~isnan(tmp(1,:)));
                    [est,ci]            = limo_central_estimator(Y,'trimmed mean');
                    TM(channel,:,:,k,1) = reshape(ci(1,:),size(data,2),size(data,3));
                    TM(channel,:,:,k,2) = reshape(est,size(data,2),size(data,3));
                    TM(channel,:,:,k,3) = reshape(ci(2,:),size(data,2),size(data,3));
                end
            end
        else
            TM = NaN(size(data,1),size(data,2),size(data,3),3);
            for k=1:size(data,3) % for each parameter
                for channel =1:size(data,1)
                    waitbar(index/(size(data,3)*size(data,1)));
                    index = index+1;
                    tmp               = squeeze(data(channel,:,k,:));
                    Y                 = tmp(:,~isnan(tmp(1,:)));
                    [est,ci]          = limo_central_estimator(Y,'trimmed mean');
                    TM(channel,:,k,1) = ci(1,:);
                    TM(channel,:,k,2) = est;
                    TM(channel,:,k,3) = ci(2,:);
                end
            end
        end
        close(h);
        
        if nargout ==0
            if nargin == 3 || nargin == 4
                newname = sprintf('%s_Trimmed_mean',name);
            else
                newname = sprintf('%s_Trimmed_mean_of_%s',name,Estimator1);
            end
            Data.trimmed_mean = TM;
            
            if exist('limo','var')
                Data.limo = limo;
            end
            save (newname,'Data');
        else
            result.trimmed_mean= TM;
        end
    end
    
    % -----------------------------------------------------
    if strcmpi(Estimator2,'HD') || strcmpi(Estimator2,'All')
        disp('Compute Harrell-Davis estimator and 95% CI ...')
        if strcmpi(limo.Analysis,'Time-Frequency')
            HD = NaN(size(data,1),size(data,2),size(data,3),size(data,4),3);
            for k = 1:size(data,4)
                for channel =1:size(data,1)
                    waitbar(index/(size(data,4)*size(data,1)));
                    index              = index+1;
                    if  strcmpi(Analysis_type,'1 channel only')
                        for f=size(data,2):-1:1
                            tmp(1,f,:,:) = data(1,f,:,k,:);
                        end
                        tmp            = limo_tf_4d_reshape(tmp,...
                            [size(data,1) size(data,2)*size(data,3) size(data,5)]);
                    else
                        tmp            = limo_tf_4d_reshape(squeeze(data(:,:,:,k,:)),...
                            [size(data,1) size(data,2)*size(data,3) size(data,5)]);
                    end
                    tmp                = squeeze(tmp(channel,:,:));
                    Y                   = tmp(:,~isnan(tmp(1,:)));
                    [est,ci]            = limo_central_estimator(Y,'HD');
                    HD(channel,:,:,k,1) = reshape(ci(1,:),size(data,2),size(data,3));
                    HD(channel,:,:,k,2) = reshape(est,size(data,2),size(data,3));
                    HD(channel,:,:,k,3) = reshape(ci(2,:),size(data,2),size(data,3));
                end
            end
        else
            HD = NaN(size(data,1),size(data,2),size(data,3),3);
            index = 1; h = waitbar(0,'computing','name','% done');
            for k=1:size(data,3)
                for channel =1:size(data,1)
                    waitbar(index/(size(data,3)*size(data,1)));
                    index             = index+1;
                    tmp               = squeeze(data(channel,:,k,:));
                    Y                 = tmp(:,~isnan(tmp(1,:)));
                    [est,ci]          = limo_central_estimator(Y,'HD');
                    HD(channel,:,k,1) = ci(1,:);
                    HD(channel,:,k,2) = est;
                    HD(channel,:,k,3) = ci(2,:);
                end
            end
        end
        close(h);
        
        if nargout ==0
            if nargin == 3 || nargin == 4
                newname = sprintf('%s_Mid_decile_Harrell_David',name);
            else
                newname = sprintf('%s_Mid_decile_Harrell_David_of_%s',name,Estimator1);
            end
            Data.Harrell_Davis = HD;
            if exist('limo','var')
                Data.limo = limo;
            end
            save (newname,'Data');
        else
            result.Harrell_Davis = HD;
        end
    end
    
    % -------------------------------------------------
    if strcmpi(Estimator2,'Median') || strcmpi(Estimator2,'All')
        disp('Compute Median estimator and 95% CI ...')
        index = 1; h = waitbar(0,'computing','name','% done');
        if strcmpi(limo.Analysis,'Time-Frequency')
            Med = NaN(size(data,1),size(data,2),size(data,3),size(data,4),3);
            for k = 1:size(data,4)
                for channel =1:size(data,1)
                    waitbar(index/(size(data,4)*size(data,1)));
                    index              = index+1;
                    if  strcmpi(Analysis_type,'1 channel only')
                        for f=size(data,2):-1:1
                            tmp(1,f,:,:) = data(1,f,:,k,:);
                        end
                        tmp            = limo_tf_4d_reshape(tmp,...
                            [size(data,1) size(data,2)*size(data,3) size(data,5)]);
                    else
                        tmp            = limo_tf_4d_reshape(squeeze(data(:,:,:,k,:)),...
                            [size(data,1) size(data,2)*size(data,3) size(data,5)]);
                    end
                    tmp                = squeeze(tmp(channel,:,:));
                    Y                    = tmp(:,~isnan(tmp(1,:)));
                    [est,ci]             = limo_central_estimator(Y,'median');
                    Med(channel,:,:,k,1) = reshape(ci(1,:),size(data,2),size(data,3));
                    Med(channel,:,:,k,2) = reshape(est,size(data,2),size(data,3));
                    Med(channel,:,:,k,3) = reshape(ci(2,:),size(data,2),size(data,3));
                end
            end
        else
            Med = NaN(size(data,1),size(data,2),size(data,3),3);
            for k=1:size(data,3)
                for channel =1:size(data,1)
                    waitbar(index/(size(data,3)*size(data,1)));
                    index = index+1;
                    tmp                = squeeze(data(channel,:,k,:));
                    Y                  = tmp(:,~isnan(tmp(1,:)));
                    [est,ci]           = limo_central_estimator(Y,'median');
                    Med(channel,:,k,1) = ci(1,:);
                    Med(channel,:,k,2) = est;
                    Med(channel,:,k,3) = ci(2,:);
                end
            end
        end
        close(h);
        
        if nargout ==0
            if nargin == 3 || nargin == 4
                newname = sprintf('%s_Median',name);
            else
                newname = sprintf('%s_Median_of_%s',name,Estimator1);
            end
            Data.median = Med;
            if exist('limo','var')
                Data.limo = limo;
            end
            save (newname,'Data');
        else
            result.Median = Med;
        end
    end
elseif isempty(data(:))
    limo_errordlg('computed central tendency is empty - nothing obtained')
    return
end

% ------------------------------
if isempty(result)
    result = 'computation done';
else
    disp('computation done');
end
