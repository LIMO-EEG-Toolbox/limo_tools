function result=limo_central_tendency_and_ci(varargin)

% The function computes robust estimates of central tendency
% (trimmed mean, Harell-Davis 0.5 decile, median) with 95% confidence intervals.
% If you input EEG data, all trials/sujects will be taken into account.
% If you input LIMO files, estimates of the 
% raw data for the categorical variables will be performed (continuous variables
% are not supported). Non overlap of 95% CI shows univariate and non
% corrected significant differences.
% If you input Betas files, estimates and 95% CI are computed for
% categorical and continuous variables of your choice.
%
% FORMAT
% limo_central_tendency_and_ci(varargin)
% result=limo_central_tendency_and_ci(varargin)
%
% if no output, ask name and save on disk
% if output, returns the 2nd level estimator
%
% INPUTS
% limo_central_tendency_and_ci(expected_chan_loc)
%                expected_chan_loc is the name of the electrode structure from EEGLAB
%                this option calls the GUI
%
% limo_central_tendency_and_ci(data, 'Analysis_type',selected_electrodes)
%                data is a [electrode * frame * trials/subjects] matrix
%                Analysis_type should be 'Mean', 'Trimmed mean', 'HD' or 'Median'
%                selected_electrodes can be [] for all brain or 1 or many channels
%                Note that this entry always does full scalp analysis
%
% limo_central_tendency_and_ci(Files, parameters, expected_chan_loc, 'Estimator1', 'Estimator2', selected_electrodes)
%                Files are the full names (with paths) of LIMO.mat files
%                parameters are which part of the raw data to analyse based on the design matrix, e.g. [1 2]
%                expected_chan_loc is the electrode structure from EEGLAB
%                Estimator should be 'Mean', 'Trimmed mean', 'HD' or 'Median'
%                Estimator 1 is applied to trials within-subjects
%                Estimator 2 is applied across subjects
%                selected_electrodes can be [] for all brain or 1 or many (=nb files) channels
%
% Guillaume Rousselet provided the initial code to do the stats
% Cyril Pernet made the interface, organize to suite EEG data etc ..
% version 1. 18 May 2010
% 26-08-2010 - GAR: figure output
% June/July 2013 - Fixed some bugs CP / thx to Benedikt Ehinger
% Novembre 2013 - fixed further issues related to parameter selection CP / thx to Matt Craddock
% version 2 - included within subjectg weighted mean + update for time frequency
% -----------------------------
%  Copyright (C) LIMO Team 2015


%% file selection and checkings
% -----------------------------
current_dir = pwd; warning off
answer = [];
result = [];
data =[];

if nargin == 3
    % ------------------------
    
    data = varargin{1};
    if length(size(data)) == 2
        tmp = NaN(1,size(data,1),size(data,2));
        tmp(1,:,:) = data;
        data = tmp; clear tmp;
    end
    
    Estimator2 = varargin{2};
    if strcmp(Estimator2,'Trimmed mean') || strcmp(Estimator2,'HD') ...
            || strcmp(Estimator2,'Median') || strcmp(Estimator2,'Mean') || strcmp(Estimator2,'All')
        Analysis_type = 'Full scalp analysis';
        parameters = 1;
    else
        error('type of estimator not recognized') ;
    end
    selected_electrodes = varargin{3};
    if ~isempty(selected_electrodes)
        data = data(selected_electrodes,:,:);
    end
    
    
elseif nargin == 6
    % ---------------------------
    
    Files = varargin{1};
    for i=1:size(Files,2)
        Paths{i} = Files{i}(1:end-9);
    end
    
    parameters          = varargin{2};
    expected_chanlocs   = varargin{3};
    Estimator1          = varargin{4};
    Estimator2          = varargin{5};
    selected_electrodes = varargin{6};
    if isempty(selected_electrodes)
        Analysis_type = 'Full scalp analysis';
    else
        Analysis_type = '1 electrode only';
        expected_chanlocs = expected_chanlocs(selected_electrodes);
    end
    
    if nargout == 0
        answer = questdlg('do you want to save individual estimates too?','saving option','yes','no','yes');
    end

    % match frames
    % -------------
    limo.data.neighbouring_matrix = expected_chanlocs;
    [first_frame,last_frame,subj_chanlocs,~,~] = limo_match_frames(Paths,limo);
    
%     disp('matching frames across subjects')
%     first_frame = NaN(1,size(Paths,2));last_frame = first_frame;
%     start = last_frame; stop = start; sampling_rate = stop;
%     for i=1:size(Paths,2)
%         cd (Paths{i});
%         load LIMO;
%         sampling_rate(i)          = LIMO.data.sampling_rate;
%         first_frame(i)            = LIMO.data.trim1;
%         last_frame(i)             = LIMO.data.trim2;
%         start(i)                  = LIMO.data.start;
%         stop(i)                   = LIMO.data.end;
%         subj_chanlocs(i).chanlocs = LIMO.data.chanlocs;
%         clear LIMO
%     end
%     
%     if (sum(sampling_rate == sampling_rate(1))) ~= length(sampling_rate)
%         errordlg('data have different sampling rates'); return
%     end
    
    
    % get data for all parameters dim [electrode, frame, param, nb subjects
    % ---------------------------------------------------------------------
    disp('gathering data ...'); index = 1;
    for i=1:size(Paths,2) % for each subject
        fprintf('processing subject %g',i); disp(' ')
        cd(Paths{i});
        load LIMO ; load Yr;
        if strcmp(LIMO.Analysis,'Time-Frequency')
            begins_at = fliplr((max(first_frame) - first_frame(i,:) + 1)); % returns time/freq/or freq-time
            ends_at(1) = size(Yr,2) - (last_frame(i,2) - min(last_frame(:,2)));
            ends_at(2) = size(Yr,3) - (last_frame(i,1) - min(last_frame(:,1)));
        else
            begins_at = max(first_frame) - first_frame(i) + 1;
            ends_at = size(Yr,2) - (last_frame(i) - min(last_frame));
        end

%         begins_at = max(first_frame) - first_frame(i) + 1;
%         ends_at = size(Yr,2) - (last_frame(i) - min(last_frame));
               
        for j=length(parameters)
            if parameters(j) < sum(LIMO.design.nb_conditions+LIMO.design.nb_interactions)
                for k=1:size(Yr,1)
                    for l=1:size(Yr,2)
                        tmp(k,l,:) =  squeeze(Yr(k,l,LIMO.design.X(:,parameters(j))==1));
                    end
                end
                
                % 1st level analysis
                % --------------------
                if strcmp(Estimator1,'Trimmed mean') % trim raw data @ 20%
                    tmp=limo_trimmed_mean(tmp,20);
                elseif strcmp(Estimator1,'Median') % median raw data
                    tmp = nanmedian(tmp,3);
                elseif strcmp(Estimator1,'HD') % mid-decile Harrell-Davis of raw data
                    tmp = limo_harrell_davis(tmp,0.5);
                elseif strcmp(Estimator1,'All') || strcmp(Estimator1,'Mean') % mean of raw data on which we do across subjects TM, HD and Median
                    % apply or not weights 
                    tmp = mean(tmp,3);
                end
                
                if strcmp(Analysis_type,'Full scalp analysis') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
                    data(:,:,j,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                elseif strcmp(Analysis_type,'1 electrode only') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
                    if size(selected_electrodes,2) == 1;
                        data(1,:,j,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                    else
                        out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                        data(1,:,j,i) = out(i,:,:); % matches the expected chanloc of the subject
                    end
                end
                clear tmp
            else
                fprintf('parameter %g not computed - continuous regressor',j);disp(' ')
            end
        end
    end
    
    
elseif nargin == 1
    % ---------------------------
    
    % Expected_chanlocs
    load(varargin{1});
    
    % check if Betas or ERPs
    option = questdlg('type of analysis','what data to analyse?','ERPs','Betas','ERPs');
    
    % -----------------------------
    % ANALYSIS ON BETAS PARAMETERS
    % -----------------------------
    if strcmp(option,'Betas')
        
        Estimator1 = 'Betas';
        Estimator2 = questdlg('Estimation option','which estimator?','Mean','Trimmed mean','HD/Median','Trimmed mean');
        if strcmp(Estimator2,'HD/Median')
            Estimator2 = 'HD';
        end
        
        % get the data
        % ------------
        go = 1; index = 1; Names = {};
        [Names,Paths,Files] = limo_get_files;
        if isempty(Names)
            return
        end
        
        
        % check type of files and returns which beta param to test
        % -------------------------------------------------------
        is_betas = [];
        for i=1:size(Names,2)
            if strcmp(Names{i},'Betas.mat')
                is_betas(i) = 1;
            end
        end
        
        if (isempty(is_betas)) == 0 && sum(is_betas) == size(Names,2)
            parameters = eval(cell2mat(inputdlg('which parameters to test e.g [1:3]','parameters option')));
            if isempty(parameters)
                return
            end
        else
            errordlg('file selection failed, only Betas.mat files are supported'); return
        end
        
        % match frames
        % ------------
     limo.data.neighbouring_matrix = expected_chanlocs;
    [first_frame,last_frame,subj_chanlocs,~,~] = limo_match_frames(Paths,limo);

%     disp('matching frames across subjects')
%         first_frame = NaN(1,size(Paths,2));last_frame = first_frame;
%         start = last_frame; stop = start; sampling_rate = stop;
%         for i=1:size(Paths,2)
%             cd (Paths{i});
%             load LIMO;
%             sampling_rate(i)          = LIMO.data.sampling_rate;
%             first_frame(i)            = LIMO.data.trim1;
%             last_frame(i)             = LIMO.data.trim2;
%             start(i)                  = LIMO.data.start;
%             stop(i)                   = LIMO.data.end;
%             subj_chanlocs(i).chanlocs = LIMO.data.chanlocs;
%             clear LIMO
%         end
%         
%         if (sum(sampling_rate == sampling_rate(1))) ~= length(sampling_rate)
%             errordlg('data have different sampling rates'); return
%         end
        
        
        % match electrodes
        % --------------
        Analysis_type   = questdlg('Rdx option','type of analysis?','Full scalp analysis','1 electrode only','Full scalp analysis');
        if isempty(Analysis_type)
            return;
        end
        
        if strcmp(Analysis_type,'1 electrode only')
            electrode = inputdlg('which electrode to analyse [?]','Electrode option'); % can be 1 nb or a vector of electrodes (electrode optimized analysis)
            if isempty(cell2mat(electrode))
                [file,dir,index] = uigetfile('*.mat','select your electrode file');
                if isempty(file)
                    return
                else
                    cd(dir); load(file); 
                    % check the vector has the same length as the number of files
                    if length(electrode_vector) ~= size(Names,2)
                        errordlg('the nb of electrodes does not match the number of subjects','Electrode error'); return;
                    end
                    % restric the channels
                    expected_chanlocs = expected_chanlocs(electrode_vector);
                end
            elseif size(eval(cell2mat(electrode)),2) == 1 || size(eval(cell2mat(electrode)),2) == size(Names,2);
                selected_electrodes = eval(cell2mat(electrode));
                expected_chanlocs = expected_chanlocs(selected_electrodes);
            else
                errordlg('the nb of electrodes does not match the number of subjects','Electrode error'); return;
            end
        end
        
        % make one large matrix
        disp('gathering data ...'); index = 1;
        for i=1:size(Paths,2) % for each subject
            fprintf('processing subject %g',i); disp(' ')
            cd(Paths{i});
            load LIMO ; load Betas;
        if strcmp(LIMO.Analysis,'Time-Frequency')
            begins_at = fliplr((max(first_frame) - first_frame(i,:) + 1)); % returns time/freq/or freq-time
            ends_at(1) = size(Yr,2) - (last_frame(i,2) - min(last_frame(:,2)));
            ends_at(2) = size(Yr,3) - (last_frame(i,1) - min(last_frame(:,1)));
        else
            begins_at = max(first_frame) - first_frame(i) + 1;
            ends_at = size(Yr,2) - (last_frame(i) - min(last_frame));
        end
        
%         begins_at = max(first_frame) - first_frame(i) + 1;
%             ends_at = size(Betas,2) - (last_frame(i) - min(last_frame));
            
            if strcmp(Analysis_type,'Full scalp analysis')
                data(:,:,:,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,squeeze(Betas(:,:,parameters)));
            elseif strcmp(Analysis_type,'1 electrode only')
                if size(selected_electrodes,2) == 1;
                    data(1,:,1:length(parameters),i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,squeeze(Betas(:,:,parameters)));
                else
                    out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,squeeze(Betas(:,:,parameters))); % out is for all expected chanlocs, i.e. across subjects
                    data(1,:,:,i) = out(i,:,:); % matches the expected chanloc of the subject
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
        [Names,Paths,Files] = limo_get_files(1);

%         go = 1; index = 1;
%         while go == 1
%             [name,path] = uigetfile('LIMO.mat; *.txt',['select LIMO file or list ',num2str(index)]);
%             if name == 0
%                 go = 0;
%             elseif name(end-2:end)=='mat'
%                 Names{index} = name;
%                 Paths{index} = path;
%                 Files{index} = sprintf('%s\%s',path,name);
%                 cd(path); cd ..
%                 index = index + 1;
%             else
%                 group_files=textread([path name],'%s','delimiter','');
%                 for f=1:length(groupe_files)
%                     temp=groupe_files{f};
%                     t=strfind(temp,'\');
%                     Names{f}=temp(t(end)+1:end);
%                     Paths{f}=temp(1:t(end));
%                     Files{f}=temp;
%                 end
%                 index = f;
%                 go = 0;
%             end
%         end
        
        % check it's LIMO.mat files and which param to test
        % --------------------------------------------------
        is_limo = [];
        for i=1:size(Names,2)
            if strcmp(Names{i},'LIMO.mat')
                is_limo(i) = 1;
            end
        end
        
        if (isempty(is_limo)) == 0 && sum(is_limo) == size(Names,2)
            Q = questdlg('Type of merging','Options','Evaluate single conditions','Pool Conditions','Evaluate single conditions');
            if strcmp(Q,'Evaluate single conditions')
                parameters = eval(cell2mat(inputdlg('which parameters to test e.g [1:3]','parameters option')));
                if isempty(parameters)
                    return
                end
            else
                parameters = eval(cell2mat(inputdlg('which parameters to pool e.g [1 3 5]','parameters option')));
            end
        else
            errordlg('file selection failed, only LIMO.mat files are supported'); return
        end
        
        % check what type of analysis
        % ---------------------------
        Analysis_type   = questdlg('Rdx option','type of analysis?','Full scalp analysis','1 electrode only','Full scalp analysis');
        if isempty(Analysis_type)
            return;
        end
        
        if strcmp(Analysis_type,'1 electrode only')
            electrode = inputdlg('which electrode to analyse [?]','Electrode option'); % can be 1 nb or a vector of electrodes (electrode optimized analysis)
            if isempty(cell2mat(electrode))
                [file,dir,index] = uigetfile('*.mat','select your electrode file');
                if isempty(file)
                    return
                else
                    cd(dir); load(file);
                    % check the vector has the same length as the number of files
                    if length(electrode_vector) ~= length(Paths)
                        errordlg('the nb of electrodes does not match the number of subjects','Electrode error'); return;
                    end
                    selected_electrodes = electrode_vector;
                    expected_chanlocs = expected_chanlocs(selected_electrodes);
                end
            elseif size(eval(cell2mat(electrode)),2) == 1 || size(eval(cell2mat(electrode)),2) == size(Names,2);
                selected_electrodes = eval(cell2mat(electrode));
                expected_chanlocs = expected_chanlocs(selected_electrodes);
            else
                errordlg('the nb of electrodes does not match the number of subjects','Electrode error'); return;
            end
        else
            selected_electrodes = [];
        end
        
        % select method
        % -------------
        [Estimator1,Estimator2]  = limo_central_tendency_questdlg(varargin);
        if isempty(Estimator1) && isempty(Estimator2)
            return
        end
        
        if strcmp(Estimator1,'All') || strcmp(Estimator1,'Mean')
            weighted_mean = questdlg('do you want to use weights to compute means?','saving option','yes','no','yes');
        end
        
        answer = questdlg('do you want to save individual estimates too?','saving option','yes','no','yes');
        
        % match frames
        % -------------
     limo.data.neighbouring_matrix = expected_chanlocs;
    [first_frame,last_frame,subj_chanlocs,~,~] = limo_match_frames(Paths,limo);

    %         disp('matching frames across subjects')
%         first_frame = NaN(1,size(Paths,2));last_frame = first_frame;
%         start = last_frame; stop = start; sampling_rate = stop;
%         for i=1:size(Paths,2)
%             cd (Paths{i});
%             load LIMO;
%             sampling_rate(i)          = LIMO.data.sampling_rate;
%             first_frame(i)            = LIMO.data.trim1;
%             last_frame(i)             = LIMO.data.trim2;
%             start(i)                  = LIMO.data.start;
%             stop(i)                   = LIMO.data.end;
%             subj_chanlocs(i).chanlocs = LIMO.data.chanlocs;
%             clear LIMO
%         end
%         
%         if (sum(sampling_rate == sampling_rate(1))) ~= length(sampling_rate)
%             errordlg('data have different sampling rates'); return
%         end
        
        
        % get data for all parameters dim [electrode, frame, param, nb subjects
        % ---------------------------------------------------------------------
        disp('gathering data ...'); index = 1;
        for i=1:size(Paths,2) % for each subject
            fprintf('processing subject %g',i); disp(' ')
            cd(Paths{i});
            load LIMO ; load Yr;
            if strcmp(LIMO.Analysis,'Time-Frequency')
                begins_at = fliplr((max(first_frame) - first_frame(i,:) + 1)); % returns time/freq/or freq-time
                ends_at(1) = size(Yr,2) - (last_frame(i,2) - min(last_frame(:,2)));
                ends_at(2) = size(Yr,3) - (last_frame(i,1) - min(last_frame(:,1)));
            else
                begins_at = max(first_frame) - first_frame(i) + 1;
                ends_at = size(Yr,2) - (last_frame(i) - min(last_frame));
            end
        
        
%         begins_at = max(first_frame) - first_frame(i) + 1;
%             ends_at = size(Yr,2) - (last_frame(i) - min(last_frame));
            
            if strcmp(Q,'Evaluate single conditions')
                for j=length(parameters)
                    if parameters(j) <= sum(LIMO.design.nb_conditions+LIMO.design.nb_interactions)
                        index = LIMO.design.X(:,parameters(j))==1;
                        if strcmp(weighted_mean,'yes')
                            for electrode=1:size(Yr,1)
                                tmp(electrode,:,:) =  squeeze(Yr(electrode,:,index)).*repmat(LIMO.design.weights(electrode,index),size(Yr,2),1);
                            end
                        else
                            tmp =  squeeze(Yr(:,:,index)); % retain those trials only
                        end
                    
                        % 1st level analysis
                        % --------------------
                        if strcmp(Estimator1,'Trimmed mean') % trim raw data @ 20%
                            tmp=limo_trimmed_mean(tmp,20);
                        elseif strcmp(Estimator1,'Median') % median raw data
                            tmp = nanmedian(tmp,3);
                        elseif strcmp(Estimator1,'HD') % mid-decile Harrell-Davis of raw data
                            tmp = limo_harrell_davis(tmp,0.5);
                        elseif strcmp(Estimator1,'All') || strcmp(Estimator1,'Mean') % mean of raw data on which we do across subjects TM, HD and Median
                            tmp = mean(tmp,3);
                        end
                        
                        if strcmp(Analysis_type,'Full scalp analysis') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
                            data(:,:,j,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                        elseif strcmp(Analysis_type,'1 electrode only') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
                            if size(selected_electrodes,2) == 1;
                                data(1,:,j,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                            else
                                out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                                data(1,:,j,i) = out(i,:,:); % matches the expected chanloc of the subject
                            end
                        end
                        clear tmp
                    else
                        fprintf('parameter %g not computed - continuous regressor',j);disp(' ')
                    end
                end
            elseif strcmp(Q,'Pool Conditions')
                if max(parameters) <= sum(LIMO.design.nb_conditions)+sum(LIMO.design.nb_interactions)
                    index = find(sum(LIMO.design.X(:,parameters)==1,2)); % find all trials from selected columns 
                    if strcmp(weighted_mean,'yes')
                        for electrode=1:size(Yr,1)
                            tmp(electrode,:,:) =  squeeze(Yr(electrode,:,index)).*repmat(LIMO.design.weights(electrode,index),size(Yr,2),1);
                        end
                    else
                        tmp =  squeeze(Yr(:,:,index)); % retain those trials only
                    end
                    
                    % 1st level analysis
                    % --------------------
                    if strcmp(Estimator1,'Trimmed mean') % trim raw data @ 20%
                        tmp=limo_trimmed_mean(tmp,20);
                    elseif strcmp(Estimator1,'Median') % median raw data
                        tmp = nanmedian(tmp,3);
                    elseif strcmp(Estimator1,'HD') % mid-decile Harrell-Davis of raw data
                        tmp = limo_harrell_davis(tmp,0.5);
                    elseif strcmp(Estimator1,'All') || strcmp(Estimator1,'Mean') % mean of raw data on which we do across subjects TM, HD and Median
                        tmp = mean(tmp,3);
                    end
                    
                    if strcmp(Analysis_type,'Full scalp analysis') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
                        data(:,:,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                    elseif strcmp(Analysis_type,'1 electrode only') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
                        if size(selected_electrodes,2) == 1;
                            data(1,:,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                        else
                            out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                            data(1,:,i) = out(i,:,:); % matches the expected chanloc of the subject
                        end
                    end
                    clear tmp
                else
                    fprintf('pooling not computed - one or more continuous regressor selected \n');
                end
            end
        end
    end
    
else
    error('nb of arguments incorrect')
end % closes varargin
  
    
%% Analysis part
cd(current_dir)

if ~isempty(data)
    % Data is either [electrode, frame, trials/subject] or [electrode,
    % frame, conditions (from parameters), subjects]
    if length(size(data)) == 3
        tmp = data; clear data
        for i=1:size(tmp,1)
            for j=1:size(tmp,2)
                data(i,j,1,:) = tmp(i,j,:); % now data is 4D
            end
        end
    end
    
    n = size(data,4);
    if n<=10 && strcmp(Estimator2,'HD')
        msgbox('CI of the Harell Davis estimates cannot be computed for less than 11 observations - switched to median','Computation info');
        Estimator2 = 'Median';
    end
    
    % save as
    if nargout ==0
        name = cell2mat(inputdlg('save as [?]','name option'));
        if ~isempty(name)
            newname = sprintf('%s_single_subjects_%s',name,Estimator1);
            save (newname,'data');
        end
    end
    
    disp('processing data across subjects ..')
    % --------------------------------------------------------------
    if strcmp(Estimator2,'Mean') || strcmp(Estimator2,'All')
        myestimator='Mean';
        disp('Compute the Mean estimator and 95% CI ...')
        M = NaN(size(data,1),size(data,2),size(data,3),3);
        for k = 1:size(data,3)
            for electrode =1:size(data,1)
                tmp = squeeze(data(electrode,:,k,:));
                Y = tmp(:,find(~isnan(tmp(1,:))));
                [m,dfe,ci,sd,n,t,p] = limo_ttest(1,Y,0,0.05);
                M(electrode,:,k,2) = m;
                M(electrode,:,k,1) = ci(:,1);
                M(electrode,:,k,3) = ci(:,2);
            end
        end
        
        if nargout ==0
            if nargin == 3
                newname = sprintf('%s_Mean',name);
            else
                newname = sprintf('%s_Mean_of_%s',name,Estimator1);
            end
            save (newname,'M');
        else
            result = M;
        end
    end
    
    % --------------------------------------------------------------
    if strcmp(Estimator2,'Trimmed mean') || strcmp(Estimator2,'All')
        myestimator='20% trimmed mean';
        disp('Compute 20% Trimmed Mean estimator and 95% CI ...')
        TM = NaN(size(data,1),size(data,2),size(data,3),3);
        for k=1:size(data,3)
            if size(data,1) == 1;
                tmp_data = NaN(1,size(data,2),size(data,4));
                tmp_data(1,:,:) = squeeze(data(:,:,k,:));
            else
                tmp_data = squeeze(data(:,:,k,:));
            end
            TM(:,:,k,:)=limo_trimmed_mean(tmp_data,20,5/100);
        end
        
        if nargout ==0
            if nargin == 3
                newname = sprintf('%s_Trimmed_mean',name);
            else
                newname = sprintf('%s_Trimmed_mean_of_%s',name,Estimator1);
            end
            save ([newname],'TM');
        else
            result = TM;
        end
    end
    
    % -----------------------------------------------------
    if strcmp(Estimator2,'HD') || strcmp(Estimator2,'All')
        myestimator='Harrell-Davis';
        disp('Compute Harrell-Davis estimator and 95% CI ...')
        HD = NaN(size(data,1),size(data,2),size(data,3),3);
        for k=1:size(data,3)
            if size(data,1) == 1;
                tmp_data = NaN(1,size(data,2),size(data,4));
                tmp_data(1,:,:) = squeeze(data(:,:,k,:));
            else
                tmp_data = squeeze(data(:,:,k,:));
            end
            HD(:,:,k,:)=limo_harrell_davis(tmp_data,.5,100);
        end
        
        if nargout ==0
            if nargin == 3
                newname = sprintf('%s_Mid_decile_Harrell_David',name);
            else
                newname = sprintf('%s_Mid_decile_Harrell_David_of_%s',name,Estimator1);
            end
            save ([newname],'HD');
        else
            result = HD;
        end
    end
    
    % -------------------------------------------------
    if strcmp(Estimator2,'Median') || strcmp(Estimator2,'All')
        myestimator='Median';
        disp('Compute Median estimator and 95% CI ...')
        Med = NaN(size(data,1),size(data,2),size(data,3),3);
        for k=1:size(data,3)
            if size(data,1) == 1;
                tmp_data = NaN(1,size(data,2),size(data,4));
                tmp_data(1,:,:) = squeeze(data(:,:,k,:));
            else
                tmp_data = squeeze(data(:,:,k,:));
            end
            Med(:,:,k,:)=limo_median(tmp_data,100);
        end
        
        if nargout ==0
            if nargin == 3
                newname = sprintf('%s_Median',name);
            else
                newname = sprintf('%s_Median_of_%s',name,Estimator1);
            end
            save ([newname],'Med');
        else
            result = Med;
        end
    end
elseif isempty(data(:))    
    warndlg('computed central tendency is empty','nothing obtained')
end
disp('computation done')

