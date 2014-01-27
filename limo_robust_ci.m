function limo_robust_ci(varargin)

% the function creates robust estimates (trimmed mean, Harell-Davis -mid 
% decile, Median) with 95% CI. It takes either EEG data and all trials/sujects
% will be taken into account or LIMO files and selective estimates from the
% raw data for categorical variable will be performed (continuous variables
% are not supported). Non overlap of 95% CI shows (univariate non corrected)
% significant differences
%
% INPUTS
% limo_robust_ci(expected_chan_loc)
%                expected_chan_loc is the name of the electrode structure from EEGlab
%                this option calls the GUI
%
% limo_robust_ci(data, 'Analysis_type',selected_electrodes)
%                data is a [electrode * frame * trials/subjects] matrix              
%                Analysis_type should be 'Trimmed mean' or 'HD' or 'Median'
%                selected_electrodes can be [] for all brain or 1 or many channels
%                Note that this entry always does full scalp analysis
%
% limo_robust_ci(Files, parameters, expected_chan_loc, 'Estimator1', 'Estimator2', selected_electrodes)
%                Files are the full names (with paths) of LIMO.mat files
%                parameters are which part of the raw data to analyse based
%                on the design matrix eg. [1 2]
%                expected_chan_loc is the electrode structure from EEGlab
%                Estimator should be 'Trimmed mean' or 'HD' or 'Median'
%                Estimator 1 is applied to trials subject-wise, Estimator 2 across subjects
%                selected_electrodes can be [] for all brain or 1 or many (=nb files) channels
%
% Guillaume Rousselet provided the initial code to do the stats
% Cyril Pernet made the interface, organize to suite EEG data etc .. 
% version 1. 18 May 2010
% -----------------------------
%  Copyright (C) LIMO Team 2010


%% file selection and checkings
% -----------------------------
current_dir = pwd; warning off

if length(varargin) == 3
    
    data = varargin{1};
    if length(size(data)) == 2  % if only one electrode the size is adjusted
        tmp = NaN(1,size(data,1),size(data,2));
        clear Data; data = tmp; clear tmp
    end
    Estimator2 = varargin{2};
    if strcmp(Estimator2,'Trimmed mean') || strcmp(Estimator,'HD') ...
            || strcmp(Estimator2,'Median') || strcmp(Estimator,'All')
        Analysis_type = 'Full scalp analysis';
        parameters = 1;
    else
        error('type of estimator not recognized') ;
    end
    selected_electrodes = varargin{3};
    if ~isempty(selected_electrodes)
        data = data(selected_electrodes,:,:);
    end 

elseif length(varargin) == 6 || length(varargin) == 1
    
    if length(varargin) == 6
        
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
         
    elseif length(varargin) == 1
        
        % Expected_chanlocs
        load(varargin{1});
        
        % select data
        % -----------
        go = 1; index = 1;
        while go == 1
            [name,path] = uigetfile('LIMO.mat',['select a LIMO file',num2str(index)]);
            if name == 0
                go = 0;
            else
                Names{index} = name;
                Paths{index} = path;
                Files{index} = sprintf('%s\%s',path,name);
                cd(path); cd ..
                index = index + 1;
            end
        end
           
        % check it's LIMO.mat files and which param to test
        % --------------------------------------------------
        is_limo = [];
        for i=1:size(Names,2)
            if strcmp(Names{i},'LIMO.mat')
                is_limo(i) = 1;
            end
        end
        
        if (isempty(is_limo)) == 0 && sum(is_limo) == size(Names,2)
            parameters = eval(cell2mat(inputdlg('which parameters to test e.g [1:3]','parameters option')));
            if isempty(parameters)
                return
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
                return
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
        Estimator1   = questdlg('type of robust estimator within subject?','Estimator option','HD', 'Median','Trimmed mean','HD');
        if isempty(Estimator1)== 1
            return;
        end
        
        Estimator2   = questdlg('type of robust estimator across subjects?','Estimator option','HD', 'Median','Trimmed mean','HD');
        if isempty(Estimator2)== 1
            return;
        end
        
    end % end variable loading
    
    % match frames
    % -------------
    disp('matching frames across subjects')
    first_frame = NaN(1,size(Paths,2));last_frame = first_frame;
    start = last_frame; stop = start; sampling_rate = stop;
    for i=1:size(Paths,2)
        cd (Paths{i});
        load LIMO;
        sampling_rate(i)          = LIMO.data.sampling_rate;
        first_frame(i)            = LIMO.data.trim1;
        last_frame(i)             = LIMO.data.trim2;
        start(i)                  = LIMO.data.start;
        stop(i)                   = LIMO.data.end;
        subj_chanlocs(i).chanlocs = LIMO.data.chanlocs;
        clear LIMO
    end
    
    if (sum(sampling_rate == sampling_rate(1))) ~= length(sampling_rate)
        errordlg('data have different sampling rates'); return
    end
    
    
    % get data for all parameters dim [electrode, frame, param, nb subjects
    % ---------------------------------------------------------------------
    disp('gathering data ...'); index = 1;
    for i=1:size(Paths,2) % for each subject
        fprintf('processing subject %g',i); disp(' ')
        cd(Paths{i});
        load LIMO ; load Yr;
        begins_at = max(first_frame) - first_frame(i) + 1;
        ends_at = size(Yr,2) - (last_frame(i) - min(last_frame));
        
        for j=parameters
            if j <= LIMO.design.nb_conditions
                for k=1:size(Yr,1)
                    for l=1:size(Yr,2)
                        tmp(k,l,:) =  squeeze(Yr(k,l,LIMO.design.X(:,j)==1));
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
                elseif strcmp(Estimator1,'All') % mean of raw data on which we do across subjects TM, HD and Median
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

else
    error('nb of arguments incorrect')
    
end % closes varargin

    
%% Analysis part
cd(current_dir)
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

% --------------------------------------------------------------
if strcmp(Estimator2,'Trimmed mean') || strcmp(Estimator2,'All')
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
    save 20percent_Trimmed_mean TM
end
    
% -----------------------------------------------------
if strcmp(Estimator2,'HD') || strcmp(Estimator2,'All')
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
    save Mid_decile_Harrell_David HD
end

% -------------------------------------------------
if strcmp(Estimator2,'Median') || strcmp(Estimator2,'All')
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
    save Median Med
end

disp('computation done')

