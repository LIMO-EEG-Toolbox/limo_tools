function limo_random_select(varargin)

% This function is used to combine parameters computed at the 1st level
% using limo_glm1. Whereas in limo_glm1 observations are assumed independent
% (i.e. N-way ANOVA/ANCOVA or Regression), limo_random_effect distinguishes
% independents and non-independent (repeated) measures.
%
% In the case of repeated measures ANOVA sphericity is accounted for using 
% a multivariate approach (effectively computing a Hotelling T2 test).
% Note that no statistical test is done in limo_random_select, only the grouping
% and organization of the data (hence the name select) - once data are
% selected and re-organized they are send to limo_random_robust which deals
% with the data structures and call the stat functions
%
% FORMAT
% limo_random_select(type,expected_chanlocs)
% limo_random_select(type,expected_chanlocs,nboot,tfce)
%
% INPUT
% type = 1 for a one sample t-test
% type = 2 for a two-samples t-test
% type = 3 for a paired t-test
% type = 4 for a regression
% type = 5 for an ANOVA
% expected_chanlocs: the EEGlab structure defining all electrodes
% (this file can be created with via limo_tools)
% nboot the number of bootstrap to do (default = 0)
% tfce 0/1 indicates to computes tfce or not (default = 0)
%
% see also limo_random_robust, limo_expected_chanlocs
% Cyril Pernet - Guillaume Rousselet v1 18-05-2009
% Cyril Pernet v3 20-05-2010
% Cyril Pernet 22-04-2011 fixed the repeated ANOVA (thx to Nicolas) 
% and added tfce
% Mariane Latinus 2013 - Update for tfce + some fix
% Cyril Pernet changed Regression / ANOVA to get structure handled within
% limo_random_effect May 2013
% Add fix from Marlene Poncet for N by N repeated neasures - July 2013
% Update for 'Frequency' and 'Time-Freuqncy' Cyril Pernet May 2014
% ---------------------------------------------------------
%  Copyright (C) LIMO Team 2014


%% take the inputs and load some parameters
nboot = 0;
tfce  = 0;

if nargin < 1
    error('not enough arguments');
elseif nargin <= 4
    type              = varargin{1};
    expected_chanlocs = varargin{2};
    if nargin == 4
        nboot         = varargin{3};
        tfce          = varargin{4};
    end
else
   error('too many arguments');    
end


Analysis_type   = questdlg('Rdx option','type of analysis?','Full scalp analysis','1 electrode only','Full scalp analysis');
if isempty(Analysis_type)
    return
end

% check chanlocs and nboot
global limo
limo.dir = pwd;

if ~isempty (expected_chanlocs)
    chan_name = expected_chanlocs;
    load(chan_name);
    % expected_chanlocs = eval(chan_name(find(chan_name == '/',1,'last')+1:end-4));
    limo.data.chanlocs = expected_chanlocs;
    limo.data.neighbouring_matrix = channeighbstructmat;
else
    limo.data.chanlocs = [];
    limo.data.neighbouring_matrix = [];
end

limo.design.bootstrap = nboot;
limo.design.tfce = tfce;
limo.Level = 2;

% ----------------------------------
%%  One sample t-test and regression
% ----------------------------------
if type == 1 || type == 4

    % get files
    % ---------
    [Names,Paths,limo.data.data] = limo_get_files(' ');
    if isempty(Names)
        disp('no files selected')
        return
    end
    limo.data.data_dir = Paths;
    
    % check type of files and returns which beta param to test
    % -------------------------------------------------------
    parameters = check_files(Names,1);
    if isempty(parameters)
        errordlg('file selection failed, only Beta and Con files are supported','Selection error'); return
    end
    
    % match frames
    % ------------
    [first_frame,last_frame,subj_chanlocs,channeighbstructmat] = match_frames(Paths);
    
    % update chanloc if needed
    % ------------------------
    if isempty (expected_chanlocs)
        limo.data.chanlocs = subj_chanlocs(1).chanlocs;
        limo.data.neighbouring_matrix = channeighbstructmat;
    end
    
    % match electrodes
    % --------------
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
                % add the name to LIMO
                if type == 1
                    limo.design.name = 'one sample t-test one electrode';
                else
                    limo.design.name = 'regression analysis one electrode';
                end
                % restric the channels
                limo.design.electrode = electrode_vector;
                expected_chanlocs = expected_chanlocs(electrode_vector);
                limo.data.chanlocs = expected_chanlocs;
            end
            
        elseif size(eval(cell2mat(electrode)),2) == 1 || size(eval(cell2mat(electrode)),2) == size(Names,2);
            if type == 1
                limo.design.name = 'one sample t-test one electrode';
            else
                limo.design.name = 'regression analysis one electrode';
            end
            % restric the channels
            limo.design.electrode = eval(cell2mat(electrode));
            expected_chanlocs = expected_chanlocs(limo.design.electrode);
            limo.data.chanlocs = expected_chanlocs;
        else
            errordlg('the nb of electrodes does not match the number of subjects','Electrode error'); return;
        end
    else
        if type == 1
            limo.design.name = 'one sample t-test all electrodes';
            limo.data.chanlocs = expected_chanlocs;
        else
            limo.design.name = 'regression analysis all electrodes';
            limo.data.chanlocs = expected_chanlocs;
        end
        limo.design.electrode = [];
    end
    
    % get data for all parameters dim [electrode, frame, param, nb subjects
    % -----------------------------------------------------------------
    disp('gathering data ...'); index = 1;
    for i=1:size(Paths,2) % for each subject
        load(limo.data.data{i});
        try
            tmp = eval(str2mat(Names{1}(1:end-4)));
        catch
            tmp = eval(str2mat(Names{1}(1:end-6)));
            tmp = squeeze(tmp(:,:,1));
        end
        begins_at = max(first_frame) - first_frame(i) + 1;
        ends_at = size(tmp,2) - (last_frame(i) - min(last_frame));
        if strcmp(Analysis_type,'Full scalp analysis') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
            data(:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
            index = index + 1; removed(i) = 0;
        elseif strcmp(Analysis_type,'1 electrode only') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
            if  length(limo.design.electrode) == 1;
                data(1,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                index = index + 1; removed(i) = 0;
            else
                out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                data(1,:,:,index) = out(i,:,:); % matches the expected chanloc of the subject
                index = index +1; removed(i) = 0;
            end
        else
            fprintf('subject %g discarded, channel description and data size don''t match',i);
            removed(i) = 1; disp(' ')
        end
        clear tmp
    end
    
    if strcmp(Analysis_type,'1 electrode only') && size(data,2) == 1
        tmp_data = squeeze(data); clear data
        data(1,1:size(tmp_data,1),1,1:size(tmp_data,2)) = tmp_data;
    end
    
    
    % one-sample t-test
    % -----------------
    if type == 1
        limo.design.X = [];
        LIMO = limo; 
        
        % clear some memory
        clear Betas Names Paths channeighbstructmat expected_chanlocs limo subj_chanlocs
        if parameters == 0; parameters = [1:size(data,3)]; end
        
        for i=parameters
            cd(LIMO.dir);
            if length(parameters) > 1
                dir_name = sprintf('parameter_%g',i);
                mkdir(dir_name); cd(dir_name);
            end
            LIMO.dir = pwd;
            
            if strcmp(Analysis_type,'1 electrode only') && size(data,1) == 1
                tmp = squeeze(data(:,:,i,:));
                tmp_data = ones(1,size(tmp,1),size(tmp,2)); % add dim 1 = 1 electrode
                tmp_data(1,:,:) = tmp; clear tmp;
            else
                tmp_data = squeeze(data(:,:,i,:));
            end
            save LIMO LIMO
            Yr = tmp_data; save Yr Yr, clear Yr % just to be consistent with name
            limo_random_robust(type,tmp_data,i,nboot,tfce)
        end
        
        
        % regression
        % -------------
        
    elseif type == 4
        cd(limo.dir)
        [FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','select regressor file');
        if FilterIndex == 0
            return
        elseif strcmp(FileName(end-3:end),'.txt')
            cd(PathName)
            X = load(FileName);
        elseif strcmp(FileName(end-3:end),'.mat')
            cd(PathName)
            load(FileName);
            X = eval(FileName(1:end-4));
        end
        disp('Regressor(s) loaded');
        
        % check size and orientation
        try
            if size(X,2) == size(data,4) || size(X,2) == size(Paths,2); disp('X has been transposed'); X = X'; end
            if size(X,1) > size(data,4);
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
                    errordlg('the number of regression value differs from the number of subjects','Covariate error'); return
                end
            end
                       
            if size(X,2)==2 && nboot ~= 599; 
                limo.design.bootstrap = 599;
                disp('nb of bootstrap adjusted to 599 for a simple regression'); 
            end
                        
        catch ME
            errordlg('covariate not loaded - make sure data are in lines or columns','Covariate error'); return
        end
        
        LIMO = limo;
        clear Betas Names Paths channeighbstructmat expected_chanlocs limo subj_chanlocs
        if parameters == 0; parameters = [1:size(data,3)]; end
        
        for i=parameters
            cd(LIMO.dir);
            if length(parameters) > 1
                dir_name = sprintf('parameter_%g',i);
                mkdir(dir_name); cd(dir_name);
            end
            LIMO.dir = pwd;
            
            if strcmp(Analysis_type,'1 electrode only')
                tmp = squeeze(data(:,:,i,:));
                tmp_data = ones(1,size(tmp,1),size(tmp,2)); % add dim 1 = 1 electrode
                tmp_data(1,:,:) = tmp; clear tmp;
            else
                tmp_data = squeeze(data(:,:,i,:));
            end
            
            % compute
            save LIMO LIMO; clear LIMO ;
            limo_random_robust(type,tmp_data,X,i,nboot,tfce)
        end
    end

    % ------------------------------
    %%  Two samples t-test
    % -----------------------------
elseif type == 2

    limo.design.X = [];
    if strcmp(Analysis_type,'Full scalp analysis') || strcmp(Analysis_type,'1 electrode only')

        N = 0;
        for gp = 1:2
            [Names{gp},Paths{gp},limo.data.data{gp}] = limo_get_files([' gp' num2str(gp)]);
            if isempty(Names{gp})
                return
            end
            limo.data.data_dir{gp} = Paths{gp};
            N = N + size(Names{gp},2);
            parameters(:,gp) = check_files(Names{gp},1);
        end

        % check type of files and returns which beta param to test
        % -------------------------------------------------------
        cd (cell2mat(Paths{1}(1)))
        load LIMO
        if max(parameters) > size(LIMO.design.X,2)
            errordlg('invalid parameter(s)','Parameters error'); return
        end

        % match frames
        % ------------
        [first_frame,last_frame,subj_chanlocs] = match_frames(Paths);


        % match electrodes
        % --------------
        if strcmp(Analysis_type,'1 electrode only')
            electrode = inputdlg('which electrode to analyse [?]','Electrode option'); % can be 1 nb or a vector of electrodes (electrode optimized analysis)
            if isempty(cell2mat(electrode))
                [file,dir,index] = uigetfile('*.mat','select your electrode file');
                if isempty(file)
                    return
                else
                    cd(dir); load(file); 
                    % check the vector has the same length as the number of files
                    if length(electrode_vector) ~= N
                        errordlg('the nb of electrodes does not match the number of subjects','Electrode error'); return;
                    end
                    % add the name to LIMO
                    limo.design.name = 'two samples t-test one electrode';
                    % restric the channels
                    limo.design.electrode = electrode_vector;
                    expected_chanlocs = expected_chanlocs(electrode_vector);
                    limo.data.chanlocs = expected_chanlocs;
                end
            elseif size(eval(cell2mat(electrode)),2) == 1 || size(eval(cell2mat(electrode)),2) == N;
                limo.design.name = 'two samples t-test one electrode';
                limo.design.electrode = eval(cell2mat(electrode));
                expected_chanlocs = expected_chanlocs(limo.design.electrode);
                limo.data.chanlocs = expected_chanlocs;
            else
                errordlg('the nb of electrodes does not match the number of subjects','Electrode error'); return;
            end
        else
            limo.design.name = 'two samples t-test all electrodes';
            limo.data.chanlocs = expected_chanlocs;
            limo.design.electrode = [];
        end


        % get data for all parameters dim [electrode, frame, param, nb subjects
        % -----------------------------------------------------------------
        disp('gathering data ...');
        subject_nb = 1;
        for g = 1:gp
            index = 1;
            for i=1:size(Paths{g},2) % for each subject per group
                load(cell2mat(limo.data.data{g}(i)));
                name = str2mat(cell2mat(Names{g}(i)));
                try
                    tmp = eval(name(1:end-4));
                catch
                    tmp = eval(name(1:end-6));
                    tmp = squeeze(tmp(:,:,1));
                end
                begins_at = max(first_frame) - first_frame(subject_nb) + 1;
                ends_at = size(tmp,2) - (last_frame(subject_nb) - min(last_frame));
                if strcmp(Analysis_type,'Full scalp analysis') && size(subj_chanlocs(subject_nb).chanlocs,2) == size(tmp,1)
                    tmp_data(:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                    index = index + 1;
                elseif strcmp(Analysis_type,'1 electrode only') && size(subj_chanlocs(subject_nb).chanlocs,2) == size(tmp,1)
                    if size(limo.design.electrode,2) == 1;
                        tmp_data(1,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                        index = index + 1;
                    else
                        out = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                        tmp_data(1,:,:,index) = out(subject_nb,:,:); % matches the expected chanloc of the subject
                        index = index +1;
                    end
                else
                    fprintf('subject %g of group %g discarded, channel description and data size don''t match',i, g); disp(' ')
                end
                clear tmp
                subject_nb = subject_nb + 1;
            end

            if strcmp(Analysis_type,'1 electrode only') && size(tmp_data,2) == 1
                tmp_data2 = squeeze(tmp_data); clear tmp_data
                tmp_data(1,1:size(tmp_data2,1),1,1:size(tmp_data2,2)) = tmp_data2; clear tmp_data2
            end

            data{g} = tmp_data;
            clear tmp tmp_data
        end
    end


    % compute
    % --------
    LIMO = limo; cd(limo.dir); 
    % clear some memory
    clear Betas Names Paths channeighbstructmat expected_chanlocs limo subj_chanlocs
    for i=parameters  % note that it allows to perfom a series of two-sample t-testss
        if strcmp(Analysis_type,'1 electrode only')
            tmp = squeeze(data{1}(:,:,i,:));
            tmp_data1 = ones(1,size(tmp,1),size(tmp,2)); % add dim 1 = 1 electrode
            tmp_data1(1,:,:) = tmp; clear tmp
            tmp = squeeze(data{2}(:,:,i,:));
            tmp_data2 = ones(1,size(tmp,1),size(tmp,2));
            tmp_data2(1,:,:) = tmp; clear tmp
        else
            tmp_data1 = squeeze(data{1}(:,:,i,:));
            tmp_data2 = squeeze(data{2}(:,:,i,:));
        end
        
        if size(tmp_data1) ~= size(tmp_data2);
            errordlg('file selection is corrupted, data sizes don''t match');
            return
        end
        
        LIMO.dir = pwd;
        save LIMO LIMO; 
        Y1r = tmp_data1; save Y1r Y1r, clear Y1r 
        Y2r = tmp_data2; save Y2r Y2r, clear Y2r 
        limo_random_robust(type,tmp_data1,tmp_data2,i,nboot,tfce)
    end
    delete data.mat


    % ------------------------------
    %%  Paired t-test
    % -----------------------------
elseif type == 3

    limo.design.X = [];
    if strcmp(Analysis_type,'Full scalp analysis') || strcmp(Analysis_type,'1 electrode only')

        [Names,Paths,limo.data.data] = limo_get_files(' ');
        if isempty(Names)
            return
        end
        limo.data.data_dir = Paths;
        N = size(Names,2);


        % check type of files and returns which beta param to test
        % -------------------------------------------------------
        parameters = check_files(Names,1);
        if size(parameters,2) == 1  % was not beta files, ie was con files
            n = Names; clear Names; Names{1} = n; clear n;
            p = Paths; clear Paths; Paths{1} = p; clear p;
            d = limo.data.data; clear limo.data.data; limo.data.data{1} = d; clear d;
            d = limo.data.data_dir; clear limo.data.data_dir; limo.data.data_dir{1} = d; clear d;
            [Names{2},Paths{2},limo.data.data{2}] = limo_get_files([' gp2']);
            if isempty(Names{2})
                return
            end
            limo.data.data_dir{2} = Paths{2};
            N = N + size(Names{2},2);
            if size(Names{1},2) ~= size(Names{2},2)
                errordlg('the nb of files differs between pairs 1 and 2','Paired t-test error'); return
            end
        elseif size(parameters,2) ~=2 % if it was beta file one needs a pair of parameters
            errordlg('2 parameters must be selected for beta files','Paired t-test error'); return
        else
            cd (cell2mat(Paths(1))); load LIMO
            if max(parameters) > size(LIMO.design.X,2)
                errordlg('invalid parameter(s)','Paired t-test error'); return
            end
        end


        % match frames
        % ------------
        [first_frame,last_frame,subj_chanlocs] = match_frames(Paths);


        % match electrodes
        % --------------
        if strcmp(Analysis_type,'1 electrode only')
            electrode = inputdlg('which electrode to analyse [?]','Electrode option'); % can be 1 nb or a vector of electrodes (electrode optimized analysis)
            if isempty(cell2mat(electrode))
                [file,dir,index] = uigetfile('*.mat','select your electrode file');
                if isempty(file)
                    return
                else
                    cd(dir); load(file); 
                    % check the vector has the same length as the number of files
                    if length(electrode_vector) ~= N
                        errordlg('the nb of electrodes does not match the number of subjects','Electrode error'); return;
                    end
                    % add the name to LIMO
                    limo.design.name = 'paired t-test one electrode';
                    % restric the channels
                    limo.design.electrode = electrode_vector;
                    expected_chanlocs = expected_chanlocs(electrode_vector);
                    limo.data.chanlocs = expected_chanlocs;
                end
            elseif size(eval(cell2mat(electrode)),2) == 1 || size(eval(cell2mat(electrode)),2) == N
                limo.design.name = 'paired t-test one electrode';
                limo.design.electrode = eval(cell2mat(electrode));
                expected_chanlocs = expected_chanlocs(limo.design.electrode);
                limo.data.chanlocs = expected_chanlocs;
            else
                errordlg('the nb of electrodes does not match the number of subjects','Electrode error'); return
            end
        else
            limo.design.name = 'paired t-test all electrodes';
            limo.data.chanlocs = expected_chanlocs;
            limo.design.electrode = [];
        end


        % get data
        % ----------
        disp('gathering data ...')

        if size(parameters,2) == 1 % select twice con files ; works as two-sample t-test

            subject_nb = 1;
            for g = 1:2
                index = 1;
                for i=1:size(Paths{g},2) % for each subject per group
                    load(cell2mat(limo.data.data{g}(i)));
                    name = str2mat(cell2mat(Names{g}(i)));
                    tmp = eval(name(1:end-6));
                    tmp = squeeze(tmp(:,:,1));
                    begins_at = max(first_frame) - first_frame(subject_nb) + 1;
                    ends_at = size(tmp,2) - (last_frame(subject_nb) - min(last_frame));
                    if strcmp(Analysis_type,'Full scalp analysis') && size(subj_chanlocs(subject_nb).chanlocs,2) == size(tmp,1)
                        tmp_data(:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                        index = index + 1;
                    elseif strcmp(Analysis_type,'1 electrode only') && size(subj_chanlocs(subject_nb).chanlocs,2) == size(tmp,1)
                        if size(limo.design.electrode,2) == 1;
                            tmp_data(1,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                            index = index + 1;
                        else
                            out = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                            tmp_data(1,:,:,index) = out(subject_nb,:,:); % matches the expected chanloc of the subject
                            index = index +1;
                        end
                    else
                        fprintf('subject %g of group %g discarded, channel description and data size don''t match',i, g); disp(' ')
                    end
                    clear tmp
                    subject_nb = subject_nb + 1;
                end

                if strcmp(Analysis_type,'1 electrode only') && size(tmp_data,2) == 1
                    tmp_data2 = squeeze(tmp_data); clear tmp_data
                    tmp_data(1,1:size(tmp_data2,1),1,1:size(tmp_data2,2)) = tmp_data2; clear tmp_data2
                end

                data{g} = tmp_data;
                clear tmp tmp_data
            end

            % final check
            if size(data{1},4) ~= size(data{2},4)
                errordlg('sample sizes don''t match, maybe some files were discarded - check the command window for messages'); return
            end

        else  % was betas files, works as one-sample t-test

            index = 1;
            for i=1:size(Paths,2) % for each subject
                load(limo.data.data{i});
                tmp = eval(str2mat(Names{1}(1:end-4)));
                begins_at = max(first_frame) - first_frame(i) + 1;
                ends_at = size(tmp,2) - (last_frame(i) - min(last_frame));
                if strcmp(Analysis_type,'Full scalp analysis') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
                    data(:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                    index = index + 1;
                elseif strcmp(Analysis_type,'1 electrode only') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
                    if size(limo.design.electrode,2) == 1;
                        data(1,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                        index = index + 1;
                    else
                        out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                        data(1,:,:,index) = out(i,:,:); % matches the expected chanloc of the subject
                        index = index +1;
                    end
                else
                    fprintf('subject %g discarded, channel description and data size don''t match',i); disp(' ')
                end
                clear tmp
            end

            if strcmp(Analysis_type,'1 electrode only') && size(data,2) == 1
                tmp_data = squeeze(data); clear data
                data(1,1:size(tmp_data,1),1,1:size(tmp_data,2)) = tmp_data;
            end
        end
    end


    % compute
    % --------
    LIMO = limo; cd(limo.dir); save LIMO LIMO;

    if strcmp(Analysis_type,'1 electrode only')
        if size(parameters,2) == 2 % beta files
            tmp = squeeze(data(:,:,parameters(1),:));
            tmp_data1 = ones(1,size(tmp,1),size(tmp,2)); % add dim 1 = 1 electrode
            tmp_data1(1,:,:) = tmp; clear tmp
            tmp = squeeze(data(:,:,parameters(2),:));
            tmp_data2 = ones(1,size(tmp,1),size(tmp,2));
            tmp_data2(1,:,:) = tmp; clear tmp
        else % con files
            tmp = squeeze(data{1}(:,:,:,:));
            tmp_data1 = ones(1,size(tmp,1),size(tmp,2)); % add dim 1 = 1 electrode
            tmp_data1(1,:,:) = tmp; clear tmp
            tmp = squeeze(data{2}(:,:,:,:));
            tmp_data2 = ones(1,size(tmp,1),size(tmp,2));
            tmp_data2(1,:,:) = tmp; clear tmp
        end
    else
        if size(parameters,2) == 2 % beta files
            tmp_data1 = squeeze(data(:,:,parameters(1),:));
            tmp_data2 = squeeze(data(:,:,parameters(2),:));
        else % con files
            tmp_data1(:,:,1,:) = squeeze(data{1}(:,:,:));
            tmp_data2(:,:,1,:) = squeeze(data{2}(:,:,:));
        end
    end
    
    Y1r = tmp_data1; save Y1r Y1r, clear Y1r 
    Y2r = tmp_data2; save Y2r Y2r, clear Y2r 
    limo_random_robust(type,tmp_data1,tmp_data2,parameters,nboot,tfce)

    % -----------------------------------
    %%  Various sorts of ANOVAs/ANCOVAs
    % -----------------------------------
elseif type == 5

    % 1st know what design it is
    % --------------------------
    
    answer = questdlg('What ANOVA model do you want to run?', 'Model selection', 'Repeated Measures', ...
        'N-Ways','ANCOVA','Repeated Measures');
    
    % get some nice comment in LIMO.mat
    if strcmp(Analysis_type,'Full scalp analysis')
        if strcmp(answer,'Repeated Measures')
            limo.design.name = 'Repeated measures ANOVA all electrodes';
        elseif strcmp(answer,'N-Ways')
            limo.design.name = 'N-ways ANOVA all electrodes';
        else
            limo.design.name = 'ANCOVA all electrodes';
        end

    elseif strcmp(Analysis_type,'1 electrode only')
        if strcmp(answer,'Repeated Measures')
            limo.design.name = 'Repeated measures ANOVA one electrode';
        elseif strcmp(answer,'N-Ways')
            limo.design.name = 'N-ways ANOVA one electrode';
        else
            limo.design.name = 'ANCOVA one electrode';
        end
    end
    
   
    % ---------------------------------------------------------------------
    %              One Way ANOVA / ANCOVA
    % ---------------------------------------------------------------------
    if strcmp(answer,'N-Ways') || strcmp(answer,'ANCOVA')
        
        % Ask for Gp 
        % -------------
        gp_nb = eval(cell2mat(inputdlg('How many independent groups? e.g. [3 2] for 3x2 ANOVA','Groups')));
        if isempty(gp_nb)
            return;
        elseif gp_nb <=1
            errordlg('at least 2 groups are expected for an ANOVA')
            return
        end
        
        % select data per gp / conditions
        % ---------------------------------
        a = questdlg('load con files or beta file','ANOVA loading files','con','beta','beta');
        N = 0; cell_nb = 1;
        if strcmp(a,'beta') % beta files
            for i=1:prod(gp_nb)
                [Names{cell_nb},Paths{cell_nb},limo.data.data{cell_nb}] = limo_get_files([' beta file gp ',num2str(i)]);
                if isempty(Names{cell_nb}); return; end
                parameters(:,i) = check_files(Names,1);
                if size(parameters,1) > 1
                    errordlg('no 2 parameters from a given subject should be used in a N-way ANOVA'); return
                end
                % create categorical file to be passed to design matrix
                N = N + size(Names{cell_nb},2);
                cell_nb = cell_nb +1;
            end
            limo.data.data_dir = Paths;
        else  % multiple con files
            for i=1:prod(gp_nb)
                [Names{cell_nb},Paths{cell_nb},limo.data.data{cell_nb}] = limo_get_files([' con file gp ',num2str(i)]);
                if isempty(Names{cell_nb}); return; end
                parameters(i) = 1;
                limo.data.data_dir{cell_nb} = Paths{cell_nb};
                N = N + size(Names{cell_nb},2);
                cell_nb=cell_nb+1;
            end
            check_files(Names,size(Names,2));
        end

        % organize the data in such a way that we can easily compute stuff
        % ---------------------------------------------------------------------
        % match frames
        [first_frame,last_frame,subj_chanlocs] = match_frames(Paths);

        % match electrodes
        if strcmp(Analysis_type,'1 electrode only')
            electrode = inputdlg('which electrode to analyse [?]','Electrode option'); % can be 1 nb or a vector of electrodes (electrode optimized analysis)
            if isempty(cell2mat(electrode))
                [file,dir,index] = uigetfile('*.mat','select your electrode file');
                if isempty(file)
                    return
                else
                    cd(dir); load(file); 
                    % check the vector has the same length as the number of files
                    if length(electrode_vector) ~= N
                        errordlg('the nb of electrodes does not match the number of subjects','Electrode error'); return;
                    end
                    % restric the channels
                    limo.design.electrode = electrode_vector;
                    expected_chanlocs = expected_chanlocs(electrode_vector);
                    limo.data.chanlocs = expected_chanlocs;
                end
            elseif size(eval(cell2mat(electrode)),2) == 1 || size(eval(cell2mat(electrode)),2) == N;
                limo.design.electrode = eval(cell2mat(electrode));
                expected_chanlocs = expected_chanlocs(limo.design.electrode);
                limo.data.chanlocs = expected_chanlocs;
            else
                error('the nb of electrodes does not match the number of subjects')
            end
        else
            limo.data.chanlocs = expected_chanlocs;
            limo.design.electrode = [];
        end

        % get all the data [electrode, frame, param] per gp
        disp('gathering data ...'); subject_index = 1; matrix_index = 1;
        for h = 1:prod(gp_nb) % each group
            nb_subjects(h) = 0;
            for i=1:size(Paths{h},2) % for each subject of the gp h
                load(cell2mat(limo.data.data{h}(i)));
                try tmp = Betas; catch tmp = con; end
                begins_at = max(first_frame) - first_frame(subject_index) + 1;
                ends_at = size(tmp,2) - (last_frame(subject_index) - min(last_frame));

                % data are of dim size(expected_chanlocs,2), latter start/earlier stop across subjects, parameters, nb of subjects
                if strcmp(Analysis_type,'Full scalp analysis') && size(subj_chanlocs(subject_index).chanlocs,2) == size(tmp,1)
                    data(:,:,:,matrix_index) = limo_match_elec(subj_chanlocs(subject_index).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                    matrix_index = matrix_index+1; removed{h}(i) = 0; nb_subjects(h) = nb_subjects(h)+1;
                elseif strcmp(Analysis_type,'1 electrode only') && size(subj_chanlocs(subject_index).chanlocs,2) == size(tmp,1)
                    if size(limo.design.electrode,2) == 1;
                        data(1,:,:,matrix_index) = limo_match_elec(subj_chanlocs(subject_index).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                    else
                        out = limo_match_elec(subj_chanlocs(subject_index).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, i.e. across subjects
                        data(1,:,:,matrix_index) = out(i,:,:); % matches the expected chanloc of the subject
                    end
                    matrix_index = matrix_index+1; removed{h}(i) = 0; nb_subjects(h) = nb_subjects(h)+1;
                else
                    fprintf('subject %g gp %g discarded, channel description and data size don''t match',i,h); disp(' ')
                    removed{h}(i) = 1;
                end
                clear tmp
                subject_index = subject_index+1;
            end
        end

        % now load covariates and check it matches data
        if strcmp(answer,'ANCOVA')
            Cont = []; index = 0;
            [FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','select covariate file');
            if FilterIndex == 0
                return
            elseif strcmp(FileName(end-3:end),'.txt')
                cd(PathName)
                X = load(FileName);
            elseif strcmp(FileName(end-3:end),'.mat')
                cd(PathName)
                load(FileName);
                X = eval(FileName(1:end-4));
            end
            disp('Regressor(s) loaded');

            try
                if size(X,2) == N; disp('regressor transposed'); X = X'; end
                if size(X,1) ~= size(data,4);
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
                        if size(X,1) ~= size(data,4);
                            errordlg('the number of regression value differs from the number of subjects'); return
                        end
                    end
                end
            catch ME
                errordlg('data format or length of the covariate does not match'); return
            end
            cov_nb = size(X,2); Cont = X; clear X
        else
            Cont = 0;
        end


        % 4th send relevant info to limo_random_robust
        % ----------------------------------------------
        for h=1:prod(gp_nb)
            for i=1:size(Paths{h},2)
                if removed{h}(i) == 1
                    limo.data.data{h}(i) = []; limo.data.data_dir{h}(i) = []; % somehow to indicate this subject is removed
                end
            end
        end

        tmp_data = NaN(size(data,1),size(data,2),size(data,4));
        % for N ways ANOVA we pass the data and X, so we simply stack
        % the data on top of each other per gp
        for i=1:prod(gp_nb)
            if i==1; from = 1; to = nb_subjects(i);
            else from = from+nb_subjects(i); to = to+nb_subjects(i); end
            current_param = parameters(i); % select only relevant parameters
            tmp_data(:,:,from:to) = squeeze(data(:,:,current_param,from:to));
        end

        % make categorical variable
        % ------------------------
        clear index
        % which numbers to use
        for i=1:length(gp_nb)
            index{i} = [1:gp_nb(i)];
            m(i) = size(index{i},2);
        end
        
        % make columns of levels per factor
        n = prod(gp_nb);
        for i=1:length(gp_nb)
            nb_values = length(index{i});
            C =  kron(eye(nb_values),ones(n/nb_values,1));
            c{i} = sum(C.*repmat(index{i},length(C),1),2); % make columns
            n = n/length(index{i});
        end
        
        % replicate columns to match size
        n = prod(gp_nb);
        for i=1:length(c)
            template(:,i) = repmat(c{i},n/length(c{i}),1);
        end
        
        % make Cat variable
        Cat = [];
        for g = 1:length(nb_subjects)
            Cat = [Cat ; repmat(template(g,:),nb_subjects(g),1)];
        end
        
            % clear some memory
            cd(limo.dir); LIMO = limo;
            save LIMO LIMO; clear LIMO limo
            clear Names Paths data channeighbstructmat expected_chanlocs subj_chanlocs
            
            % do the analysis
            limo_random_robust(type,tmp_data,Cat,Cont,nboot,tfce)
   
            
    else
        % ---------------------------------------------------------------------
        %              Repeated measure ANOVA
        % ---------------------------------------------------------------------
        
        % Ask for Gp 
        % -------------
        gp_nb = eval(cell2mat(inputdlg('How many independent groups? e.g. 2','Groups')));
        if isempty(gp_nb)
            return;
        elseif length(gp_nb) > 1
            errordlg('only 1 independent factor (with n groups) is handled by LIMO')
            return
        elseif gp_nb == 0;
            gp_nb = 1;
        end
        
        % Ask for Repeated Measures
        % --------------------------
        factor_nb = eval(cell2mat(inputdlg('How many repeated factors? e.g. [2 3] for 2 levels F1 and 3 levels F2','Factors')));
        % in case the inpuit is [] or [0]
        if isempty(factor_nb)
            return;
        end
        if length(factor_nb)==1 && factor_nb == 0
            return
        end

        % 2nd select data per gp / conditions
        % ---------------------------------------------------------------
        a = questdlg('load several con files per subject or one beta file','ANOVA loading files','con','beta','beta');
        N = 0; cell_nb = 1;
        % beta files
        if strcmp(a,'beta')
            for i=1:gp_nb
                [Names{cell_nb},Paths{cell_nb},limo.data.data{cell_nb}] = limo_get_files([' beta file gp ',num2str(i)]);
                if isempty(Names{cell_nb}); return; end
                parameters(:,i) = check_files(Names,1);
                if length(parameters(:,i)) ~= prod(factor_nb)
                    error(['the number of parameter chosen (',num2str(length(parameters)), ...
                        ') does not match the total number of levels (',num2str(prod(factor_nb)),')'])
                end
                N = N + size(Names{cell_nb},2);
                cell_nb = cell_nb +1;
            end
            limo.data.data_dir = Paths;

        else  % multiple con files
            for i=1:gp_nb
                for j=1:length(factor_nb)
                    for k=1:factor_nb(j)
                        [names{k},paths{k},full_names{k}] = limo_get_files([' gp ',num2str(i),' factor ',num2str(j),' level ',num2str(k)]);
                        if isempty(names{k}); return; end
                        N = N + size(names{cell_nb},2);
                    end
                end
                
                Names{cell_nb} = names{1};
                Paths{cell_nb} = paths{1};
                limo.data.data{cell_nb} = full_names{1};
                for l=2:size(names,2)
                    Names{cell_nb} = [Names{cell_nb} names{l}];
                    Paths{cell_nb} = [Paths{cell_nb} paths{l}];
                    limo.data.data{cell_nb} = [limo.data.data{cell_nb} full_names{l}];
                end
                limo.data.data_dir{cell_nb} = Paths{cell_nb};
                check_files(Names,size(Names,2));
                cell_nb=cell_nb+1;
                parameters(i) = 1;
            end
        end
      
        
        % 3rd organize data
        % ---------------------------------------------------------------------
        % match frames
        [first_frame,last_frame,subj_chanlocs] = match_frames(Paths);

        % match electrodes
        if strcmp(Analysis_type,'1 electrode only')
            electrode = inputdlg('which electrode to analyse [?]','Electrode option'); % can be 1 nb or a vector of electrodes (electrode optimized analysis)
            if isempty(cell2mat(electrode))
                [file,dir,index] = uigetfile('*.mat','select your electrode file');
                if isempty(file)
                    return
                else
                    cd(dir); load(file); 
                    % check the vector has the same length as the number of files
                    if length(electrode_vector) ~= N
                        errordlg('the nb of electrodes does not match the number of subjects','Electrode error'); return;
                    end
                    % restric the channels
                    limo.design.electrode = electrode_vector;
                    expected_chanlocs = expected_chanlocs(electrode_vector);
                    limo.data.chanlocs = expected_chanlocs;
                end
            elseif size(eval(cell2mat(electrode)),2) == 1 || size(eval(cell2mat(electrode)),2) == N;
                limo.design.electrode = eval(cell2mat(electrode));
                expected_chanlocs = expected_chanlocs(limo.design.electrode);
                limo.data.chanlocs = expected_chanlocs;
            else
                error('the nb of electrodes does not match the number of subjects')
            end
        else
            limo.data.chanlocs = expected_chanlocs;
            limo.design.electrode = [];
        end

        % get data for all parameters dim [electrode, frame, param]
        disp('gathering data ...');        
        subject_index = 1;
        matrix_index = 1;
        for h = 1:gp_nb % each group
            nb_subjects(h) = 0;
            for i=1:size(Paths{h},2) % for each subject of the gp h % i=1:size(Paths{h*sum(factor_levels)},2)
                load(cell2mat(limo.data.data{h}(i)));
                try
                    tmp = Betas;
                catch
                    tmp = squeeze(con(:,:,1));
                end
                begins_at = max(first_frame) - first_frame(subject_index) + 1;
                ends_at = size(tmp,2) - (last_frame(subject_index) - min(last_frame));

                % data are of dim size(expected_chanlocs,2), latter start/earlier stop across subjects, parameters, nb of subjects
                if strcmp(Analysis_type,'Full scalp analysis') && size(subj_chanlocs(subject_index).chanlocs,2) == size(tmp,1)
                    data(:,:,:,matrix_index) = limo_match_elec(subj_chanlocs(subject_index).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                    matrix_index = matrix_index+1; removed{h}(i) = 0; nb_subjects(h) = nb_subjects(h)+1;
                elseif strcmp(Analysis_type,'1 electrode only') && size(subj_chanlocs(subject_index).chanlocs,2) == size(tmp,1)
                    if size(limo.design.electrode,2) == 1;
                        data(1,:,:,matrix_index) = limo_match_elec(subj_chanlocs(subject_index).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                    else
                        out = limo_match_elec(subj_chanlocs(subject_index).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                        data(1,:,:,matrix_index) = out(i,:,:); % matches the expected chanloc of the subject
                    end
                    matrix_index = matrix_index+1; removed{h}(i) = 0; nb_subjects(h) = nb_subjects(h)+1;
                else
                    fprintf('subject %g gp %g discarded, channel description and data size don''t match',i,h); disp(' ')
                    removed{h}(i) = 1;
                end
                clear tmp
                subject_index = subject_index+1;
            end
        end

        % 4th send relevant info to limo_random_robust
        % ----------------------------------------------
        for h=1:gp_nb
            for i=1:size(Paths{h},2)
                if removed{h}(i) == 1
                    limo.data.data{h}(i) = []; limo.data.data_dir{h}(i) = []; % somehow to indicate this subject is removed
                end
            end
        end

        % the expected dim in limo_rep_anova are [frames, subjects,
        % conditions] so we send to limo_random_robust data with dim
        % [electrodes, frames, subjects, conditions] + one vector
        % describing the group belonging
        
        if strcmp(a,'con')
            nb_subjects = nb_subjects / prod(factor_nb);
            tmp_data = NaN(size(data,1),size(data,2),sum(nb_subjects), prod(factor_nb));
            gp = NaN(sum(nb_subjects),1);
            
            for i=1:gp_nb
                if i==1; from = 1; index1 = from ; to = nb_subjects(i); index2 = to;
                else
                    from = from+nb_subjects(i-1); to = to+nb_subjects(i);
                    index2 = index1+nb_subjects(i)-1;
                end
                gp(from:to) = i;
                
                for j=1:prod(factor_nb)
                    tmp_data(:,:,from:to,j) = squeeze(data(:,:,1,index1:index2));
                    index1 = index1+nb_subjects(i); index2 = index2+nb_subjects(i);
                end
            end
            
        else
            tmp_data = NaN(size(data,1),size(data,2),sum(nb_subjects), prod(factor_nb));
            gp = NaN(sum(nb_subjects),1);
            
            for i=1:gp_nb
                current_param = parameters(:,i); % select only relevant parameters (could be different from different groups)
                
                if i==1; from = 1; to = nb_subjects(i);
                else from = from+nb_subjects(i-1); to = to+nb_subjects(i); end
                gp(from:to) = i;
                
                for j=1:prod(factor_nb)
                    tmp_data(:,:,from:to,j) = squeeze(data(:,:,current_param(j),from:to));
                end
            end
        end
        
        LIMO = limo; cd(limo.dir); save LIMO LIMO
        Yr = tmp_data; save Yr Yr;
        clear Betas LIMO Yr channeighbstructmat data expected_chanlocs Names Paths limo subj_chanlocs
        limo_random_robust(type+1,tmp_data,gp,factor_nb,nboot,tfce)
    end
end % closes type


end % closes the function


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                                   ROUTINES
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% file checking subfunction
function parameters = check_files(Names,gp)
% after selecting file, check they are all the same type

parameters = [];
if gp == 1

    % one sample case
    % ---------------
    is_beta = []; is_con = [];
    for i=1:size(Names,2)
        if strcmp(Names{i},'Betas.mat')
            is_beta(i) = 1;
        elseif strcmp(Names{i}(1:3),'con')
            is_con(i) = 1;
        end
    end

    if (isempty(is_beta)) == 0 && sum(is_beta) ~= size(Names,2) || (isempty(is_con)) == 0 && sum(is_con) ~= size(Names,2)
        error('file selection failed, only Beta or Con files are supported')
    elseif (isempty(is_beta)) == 0 && sum(is_beta) == size(Names,2)
        parameters = eval(cell2mat(inputdlg('which parameters to test e.g [1:3]','parameters option')));
        if isempty(parameters)
            return
        end
    elseif (isempty(is_con)) == 0 && sum(is_con) == size(Names,2)
        parameters = 1;
    end

elseif gp > 1

    % several sample cases
    % -------------------
    for g = 1:gp
        is_beta = []; is_con = [];
        for i=1:size(Names{g},2)
            if strcmp(Names{g}(i),'Betas.mat')
                is_beta(i) = 1;
            elseif strcmp(Names{g}{i}(1:3),'con')
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
    elseif (isempty(is_beta)) == 0 && sum(cell2mat(test)) == size(Names,2)
        parameters = eval(cell2mat(inputdlg('which parameter(s) to test e.g 1','parameters option')));
        if isempty(parameters) || size(parameters,2) ~=1 && size(parameters,2) ~=gp
            return
        end
    elseif (isempty(is_con)) == 0 && sum(cell2mat(test)) == size(Names,2)
        parameters = 1;
    end

end
end


%% frame matching subfunction
function [first_frame,last_frame,subj_chanlocs,channeighbstructmat] = match_frames(Paths)

% once we have all the files, we need to collect
% information to match the frames across subjects

global limo
channeighbstructmat = []; ME = [];

disp('match frames between subjects ...')
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

first_frame = NaN(1,size(Paths,2));last_frame = first_frame;
start = last_frame; stop = start; sampling_rate = stop;
for i=1:size(Paths,2)
    try
        cd (Paths{i});
    catch
        cd (cell2mat(Paths{i}))
    end
    load LIMO;
    
    if i==1
        Analysis = LIMO.Analysis;
    else
        if ~strcmp(LIMO.Analysis,Analysis)
            error('Looks like different type of analyses (Time/Freq/Time-Freq) are mixed up')
        end
    end
    sampling_rate(i)          = LIMO.data.sampling_rate;
    first_frame(i)            = LIMO.data.trim1;
    last_frame(i)             = LIMO.data.trim2;
    start(i)                  = LIMO.data.start;
    stop(i)                   = LIMO.data.end;
    subj_chanlocs(i).chanlocs = LIMO.data.chanlocs;
    try
       channeighbstructmat = LIMO.data.channeighbstructmat;
    catch ME
    end
end

if ~isempty(ME) && isempty(limo.data.neighbouring_matrix)
    error(sprintf('some subject(s) have a different channel structure \nplease load an expected chanloc when choosing a test'));                          
end

if (sum(sampling_rate == sampling_rate(1))) ~= length(sampling_rate)
    error('data have different sampling rates')
end

limo.Analysis = Analysis;
limo.data.sampling_rate = sampling_rate(1);
[v,c] = max(first_frame);
limo.data.trim1 = v;
limo.data.start = start(c);
[v,c] = min(last_frame);
limo.data.trim2 = v;
limo.data.end = stop(c);
if strcmp(Analysis,'Frequency')
   a = find((LIMO.data.freqlist-limo.data.start)==0);
   b = find((LIMO.data.freqlist-limo.data.end)==0);
   limo.data.freqlist = LIMO.data.freqlist(a:b);  % assumes all subject match this
end
end



