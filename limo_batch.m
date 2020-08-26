function limo_batch

% interactive function to run several 1st level analyses one after the
% others - select directories and files - possibly enter contrasts of
% interests and let it run.
%
% Cyril Pernet and Nicolas Chauveau
% CP 24-06-2013 updated to be even more automatic + fix for new designs
% -----------------------------
% Copyright (C) LIMO Team 2010


%% Defaults
bootstrap = 0;
tfce      = 0;

%% -------------------------------------
% 1 - specify the folder of each subject
% and for each subject specify data set
% categorical vector, continuous vector,
% start and end of analysis and if the
% continuous regressor should not be
% z-scored (by default it is)

current_dir = pwd;
default = [];

option = questdlg('batch mode','option','model specification','contrast only','model specification');
N = inputdlg('how many subjects do you want to batch?', 'Define N');
if isempty(N) || str2double(cell2mat(N)) == 0
    disp('exiting batch mode'); limo_gui; return
else
    N = str2double(cell2mat(N));
end


% -----------------------------------------
% load the expected channel locations
% -----------------------------------------
if bootstrap == 1
    [chan_file,chan_path,whatsup]=uigetfile('expected_chanlocs.mat','Select channel location file');
    if whatsup == 1
        load (sprintf('%s%s',chan_path,chan_file))
        test = eval(chan_file(1:end-4));
        if isstruct(test) && ~isempty(test(1).labels) && ~isempty(test(1).theta) && ~isempty(test(1).radius) ...
                && ~isempty(test(1).X) && ~isempty(test(1).Y) && ~isempty(test(1).Z) && ~isempty(test(1).sph_theta) ...
                && ~isempty(test(1).sph_phi) && ~isempty(test(1).sph_radius) && ~isempty(test(1).urchan)
            disp('channel location loaded');
        else
            warndlg('this file is not recognize as a channel location file or informations are missing','file error')
        end
    else
        disp('exiting batch mode'); limo_gui; return
    end
end

if strcmp(option,'model specification')
    
    % --------------------------------------------------------------
    % do subject 1 and ask if the same operation has to be repeated
    % --------------------------------------------------------------
    
    dirname = inputdlg('default directory name', 'Name');
    
    % directory
    % ---------
    subject_dir = uigetdir(current_dir,'Select data directory for subject 1');
    if subject_dir == 0
        disp('exiting batch mode'); limo_gui; return
    else
        cd(subject_dir)
        if ~strcmp(dirname,'')
            mkdir(cell2mat(dirname))
            cd(cell2mat(dirname))
        end
        batch_input{1}.directory = pwd;
        cd(subject_dir)
    end
    
    % data
    % -----
    name = dir('*.set');
    try % if .set in subject dir just load it
        EEG=pop_loadset(name.name);
        batch_input{1}.data_dir = pwd;
        batch_input{1}.data     = name.name;
        batch_input{1}.chanlocs = EEG.chanlocs;
        batch_input{1}.start    = EEG.xmin;
        batch_input{1}.end      = EEG.xmax;
        batch_input{1}.rate     = EEG.srate;
        fprintf('Data set %s loaded',name.name); disp(' ')
    catch
        [FileName,PathName,FilterIndex]=uigetfile('*.set','EEGLAB EEG epoch data - Subject 1');
        if FilterIndex == 0
            disp('exiting batch mode'); limo_gui; return
        else
            current_dir = pwd;
            cd(PathName)
            try
                EEG=pop_loadset(FileName);
                batch_input{1}.data_dir = PathName;
                batch_input{1}.data     = FileName;
                batch_input{1}.chanlocs = EEG.chanlocs;
                batch_input{1}.start    = EEG.xmin;
                batch_input{1}.end      = EEG.xmax;
                batch_input{1}.rate     = EEG.srate;
                fprintf('Data set %s loaded',FileName); disp(' ')
            catch
                errordlg('pop_loadset eeglab function not found','error');
            end
            cd(batch_input{1}.data_dir);%Nicolas 18/08/2011
        end
    end
    
    % categorical vector
    % ------------------
    batch_input{1}.Cat = [];
    batch_input{1}.design_type = 0;
    [FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','LIMO categorical data - Subject 1');
    if FilterIndex == 1
        cd(PathName);
        if strcmp(FileName(end-3:end),'.txt')
            batch_input{1}.Cat = load(FileName);
        else
            load(FileName)
            batch_input{1}.Cat = eval(FileName(1:end-4));
        end
        disp('Categorical data loaded');
        
        if size(batch_input{1}.Cat,2) > 1
            button = questdlg('do a full factorial?','multi-factors analysis','YES','NO','NO');
            if strcmp(button,'NO')
                batch_input{1}.design_type = 0;
            else
                batch_input{1}.design_type = 1;
            end
        else
            batch_input{1}.design_type = 0;
        end
        CatName = FileName;
    end
    cd(batch_input{1}.data_dir); %Nicolas 18/08/2011
    
    % continuous vector
    % ------------------
    batch_input{1}.Cont = [];
    [FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','LIMO continuous data - Subject 1');
    if FilterIndex == 1
        cd(PathName);
        if strcmp(FileName(end-3:end),'.txt')
            batch_input{1}.Cont = load(FileName);
        else
            load(FileName)
            batch_input{1}.Cont = eval(FileName(1:end-4));
        end
        disp('Continuous data loaded');
        ContName = FileName;
    end
    cd(batch_input{1}.directory);
    
    % quick check
    if isempty(batch_input{1}.Cat) && isempty(batch_input{1}.Cont)
        errordlg('no variables selected - batch aborded')
        return
    end
    
    % zscoring
    % ----------
    if ~isempty(batch_input{1}.Cont)
        button = questdlg('zscore continuous regressor(s) [default]?','ZSCORE','YES','NO','YES');
        if strcmp(button,'YES')
            batch_input{1}.zscore = 1;
            disp('zscoring on');
        else
            batch_input{1}.zscore = 0;
            disp('zscoring off');
        end
    end
    
    if N>1 && ~isempty(batch_input{1}.Cont)
        button = questdlg('apply for same zscoring option for all subjects?','Default','YES','NO','YES');
        if strcmp(button,'YES')
            default.zscore = 1;
            default.zscore_value = batch_input{1}.zscore;
        else
            default.zscore = 0;
        end
    else
        batch_input{1}.zscore = 0;
        default.zscore = 0;
    end
    
    % start and end
    % -------------
    prompt = {'Starting time of the analysis in sec:','Ending time of the analysis in sec:'};
    dlg_title = 'Start/Stop timinig';
    num_lines = 1;
    def = {num2str(EEG.xmin), num2str(EEG.xmax)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    % start
    v = str2num(answer{1})*1000; % change to ms
    difference = rem(v,(1/EEG.srate*1000));
    if difference ~=0
        v = v+difference;
        [value,position]=min(abs(EEG.times - v));
        v = EEG.times(position);
        h = warndlg(sprintf('adjusting to sampling rate start at %g ms',v),'adjusting stating point'); uiwait(h)
    end
    
    start = v/1000;  % back in sec
    if start < EEG.xmin
        errordlg('error in the starting point input')
    else
        batch_input{1}.start    = start;
        batch_input{1}.trim1    = find(EEG.times == v); % gives the 1st column to start the analysis
    end
    
    % stop
    v = str2num(answer{2})*1000;
    difference = rem(v,(1/EEG.srate*1000));
    if difference ~=0
        v = v+difference;
        [value,position]=min(abs(EEG.times - v));
        v = EEG.times(position);
        h = warndlg(sprintf('adjusting to sampling rate stops at %g ms',v),'adjusting stating point'); uiwait(h)
    end
    
    ending = v/1000;  % back in sec
    if ending > EEG.xmax
        errordlg('error in the ending point input')
    else
        batch_input{1}.end     = ending;
        batch_input{1}.trim2   = find(EEG.times == v); % gives the last column to end the analysis
    end
    
    % default start - stop
    if N>1
        button = questdlg(['use same start ' num2str(start) ' and end ' num2str(ending) ' points for all subjects?'],'Default timing','YES','NO','YES');
        if strcmp(button,'YES')
            default.timing = 1;
            default.start    = start;
            default.trim1    = batch_input{1}.trim1; % gives the 1st column to start the analysis
            default.end     = ending;
            default.trim2   = batch_input{1}.trim2; % gives the last column to end the analysis
        else
            default.timing = 0;
        end
    end
    
    % create LIMO.mat
    % ---------------
    LIMO.dir                    = batch_input{1}.directory;
    LIMO.data.data_dir          = batch_input{1}.data_dir;
    LIMO.data.data              = batch_input{1}.data;
    LIMO.data.chanlocs          = batch_input{1}.chanlocs;
    LIMO.data.start             = batch_input{1}.start;
    LIMO.data.end               = batch_input{1}.end ;
    LIMO.data.sampling_rate     = batch_input{1}.rate;
    LIMO.data.Cat               = batch_input{1}.Cat;
    LIMO.data.Cont              = batch_input{1}.Cont;
    LIMO.design.zscore          = batch_input{1}.zscore;
    LIMO.data.trim1             = batch_input{1}.trim1;
    LIMO.data.trim2             = batch_input{1}.trim2;
    LIMO.design.method           = 'OLS';
    LIMO.design.type_of_analysis = 'Mass-univariate';
    LIMO.design.fullfactorial    = batch_input{1}.design_type;
    LIMO.design.bootstrap        = bootstrap;
    LIMO.design.tfce             = tfce;
    LIMO.Level                   = 1;
    cd (LIMO.dir);
    save LIMO LIMO
    
    % ---------------------------------------------------
    % loop to create the LIMO.mat for all subjects
    % ---------------------------------------------------
    
    
    S = 1;
    while S<N
        S = S+1; % which subjects
        
        % directory
        % ---------
        subject_dir = uigetdir(current_dir,['Select data directory for subject ' num2str(S)]);
        if subject_dir == 0;
            disp('batching aborded ..')
            return
        end
        
        cd(subject_dir);
        if ~strcmp(dirname,'')
            mkdir(cell2mat(dirname)) ;
            cd(cell2mat(dirname)) ;
        end
        batch_input{S}.directory = pwd;
        cd(subject_dir);
        
        % data
        % -----
        name = dir('*.set');
        try % if .set in subject dir just load it
            EEG=pop_loadset(name.name);
            batch_input{S}.data_dir = pwd;
            batch_input{S}.data     = name.name;
            batch_input{S}.chanlocs = EEG.chanlocs;
            batch_input{S}.start    = EEG.xmin;
            batch_input{S}.end      = EEG.xmax;
            batch_input{S}.rate     = EEG.srate;
            fprintf('Data set %s loaded',name.name); disp(' ')
        catch
            [FileName,PathName,FilterIndex]=uigetfile('*.set',['EEGLAB EEG epoch data - Subject ' num2str(S)]);
            if FilterIndex == 0
                S = S-1; % allows user to cancel and restrat at subject S-1
            else
                current_dir = pwd;
                cd(PathName)
                try
                    EEG=pop_loadset(FileName);
                    batch_input{S}.data_dir = PathName;
                    batch_input{S}.data     = FileName;
                    batch_input{S}.chanlocs = EEG.chanlocs;
                    batch_input{S}.start    = EEG.xmin;
                    batch_input{S}.end      = EEG.xmax;
                    batch_input{S}.rate     = EEG.srate;
                    fprintf('Data set %s loaded',FileName); disp(' ')
                catch ME
                    error(['pop_loadset eeglab function not found subject' num2str(S)],'error');
                end
            end
        end
        cd(batch_input{S}.directory);
        
        
        % categorical vector
        % ------------------
        batch_input{S}.Cat = [];
        try
            if strcmp(CatName(end-3:end),'.txt')
                batch_input{S}.Cat = load(CatName);
            else
                load(CatName)
                batch_input{S}.Cat = eval(CatName(1:end-4));
            end
            disp('Categorical data loaded');
        catch
            [FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat',['LIMO categorical data - Subject ' num2str(S)]);
            if FilterIndex == 1
                cd(PathName);
                if strcmp(FileName(end-3:end),'.txt')
                    batch_input{S}.Cat = load(FileName);
                else
                    load(FileName)
                    batch_input{S}.Cat = eval(FileName(1:end-4));
                end
                disp('Categorical data loaded');
            end
        end
        cd(batch_input{S}.directory);
        
        % continuous vector
        % ------------------
        if ~isempty( batch_input{1}.Cont )
            batch_input{S}.Cont = [];
            try
                if strcmp(ContName(end-3:end),'.txt')
                    batch_input{S}.Cont = load(ContName);
                else
                    load(ContName)
                    batch_input{S}.Cont = eval(ContName(1:end-4));
                end
                disp('Continuous data loaded');
            catch
                [FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat',['LIMO continuous data - Subject ' num2str(S)]);
                if FilterIndex == 1
                    cd(PathName);
                    if strcmp(FileName(end-3:end),'.txt')
                        batch_input{S}.Cont = load(FileName);
                    else
                        load(FileName)
                        batch_input{S}.Cont = eval(FileName(1:end-4));
                    end
                    disp('Continuous data loaded');
                end
            end
        else
            batch_input{S}.Cont =[];
        end
        cd(batch_input{1}.directory);
        
        % zscoring
        % ----------
        if ~isempty( batch_input{1}.Cont )
            if default.zscore == 0
                button = questdlg('zscore continuous regressor(s)?','ZSCORE','YES','NO','NO');
                if strcmp(button,'YES')
                    batch_input{S}.zscore = 1;
                    disp('zscoring on');
                else
                    batch_input{S}.zscore = default.zscore_value;
                    disp('zscoring off');
                end
            else
                batch_input{S}.zscore = default.zscore_value;
            end
        else
            batch_input{S}.zscore = default.zscore;
        end
        
        % start and end
        % -------------
        if default.timing == 0
            
            prompt = {'Starting time of the analysis in sec:','Ending time of the analysis in sec:'};
            dlg_title = 'Start/Stop timinig';
            num_lines = 1;
            def = {num2str(EEG.xmin), num2str(EEG.xmax)};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            
            % start
            v = str2num(answer{1})*1000; % change to ms
            difference = rem(v,(1/EEG.srate*1000));
            if difference ~=0
                v = v+difference;
                [value,position]=min(abs(EEG.times - v));
                v = EEG.times(position);
                h=warndlg(sprintf('adjusting to sampling rate start at %g ms',v),'adjusting stating point'); uiwait(h)
            end
            
            start = v/1000;  % back in sec
            if start < EEG.xmin
                errordlg('error in the starting point input')
            else
                batch_input{S}.start    = start;
                batch_input{S}.trim1    = find(EEG.times == v); % gives the 1st column to start the analysis
            end
            
            % stop
            v = str2num(answer{2})*1000;
            difference = rem(v,(1/EEG.srate*1000));
            if difference ~=0
                v = v+difference;
                [value,position]=min(abs(EEG.times - v));
                v = EEG.times(position);
                h=warndlg(sprintf('adjusting to sampling rate start at %g ms',v),'adjusting stating point'); uiwait(h)
            end
            
            ending = v/1000;  % back in sec
            if ending > EEG.xmax
                errordlg('error in the ending point input')
            else
                batch_input{S}.end     = ending;
                batch_input{S}.trim2   = find(EEG.times == v); % gives the last column to end the analysis
            end
            
        else
            batch_input{S}.start   = default.start;
            batch_input{S}.trim1   = default.trim1;
            batch_input{S}.end     = default.end;
            batch_input{S}.trim2   = default.trim2;
        end
        
        % create LIMO.mat
        % ---------------
        LIMO.dir                     = batch_input{S}.directory;
        LIMO.data.data_dir           = batch_input{S}.data_dir;
        LIMO.data.data               = batch_input{S}.data;
        LIMO.data.chanlocs           = batch_input{S}.chanlocs;
        LIMO.data.start              = batch_input{S}.start;
        LIMO.data.end                = batch_input{S}.end ;
        LIMO.data.sampling_rate      = batch_input{S}.rate;
        LIMO.data.Cat                = batch_input{S}.Cat;
        LIMO.data.Cont               = batch_input{S}.Cont;
        LIMO.design.zscore           = batch_input{S}.zscore;
        LIMO.data.trim1              = batch_input{S}.trim1;
        LIMO.data.trim2              = batch_input{S}.trim2;
        LIMO.design.method           = 'OLS';
        LIMO.design.type_of_analysis = 'Mass-univariate';
        LIMO.design.fullfactorial    = batch_input{1}.design_type;
        LIMO.design.bootstrap        = bootstrap;
        LIMO.design.tfce             = tfce;
        LIMO.Level                   = 1;
        cd (LIMO.dir);
        save LIMO LIMO
    end
end


% check chanlocs for bootstrap
% ----------------------------
if bootstrap == 1
    neighbourdist = eval(cell2mat(inputdlg('enter neighbourhood distance','neighbourhood distance'))); % 0.37 for biosemi 128;
    for S = 1:N
        cd(batch_input{S}.directory); load LIMO
        channeighbstructmat = limo_expected_chanlocs(LIMO.data.data, LIMO.data.data_dir,neighbourdist);
        LIMO.data.neighbouring_matrix = channeighbstructmat;
        save LIMO LIMO
    end
end


disp('all data loaded - ready to start analyzing ...')


%% ----------------------------------------
% Contrast option

if strcmp(option,'model specification') % if model spec, option to add contrasts
    option = questdlg('add contrast','option','add contrasts','no contrasts','add contrasts');
end

if strcmp(option,'contrast only') || strcmp(option,'add contrasts')
    contrast_kind = questdlg('T or F contrast for all subjects?','Kind of contrast', 'T','F','T');
    
    if strcmp(option,'contrast only') % if contrast only, can be different for each subjects
        button = questdlg('Making the same contrasts for all subjects?','Contrasts dlg','YES','NO','YES');
    else
        button = 'YES';
    end
    
    if strcmp(button,'YES')
        go = 1; index = 1;
        while go ~=0
            C1 = inputdlg(['Enter contrast ' num2str(index)], 'Contrasts');
            if isempty(C1)
                go = 0; break
            else
                C1 = str2num(cell2mat(C1));
                default.contrast{index,1}= C1;
                index=index+1;
            end
        end
    end
end

%% ------------------------------
% loop through subjects to create the design matrix, evaluate
% and possibly do contrasts

if strcmp(option,'contrast only')
    
    % directory
    % ---------
    for S=1:N
        subject_dir = uigetdir(current_dir,['Select LIMO directory for subject ' num2str(S)]);
        batch_input{S}.directory = subject_dir;
    end
    
    % contrasts
    % ---------
    if isfield(default,'contrast')
        for S=1:N
            cd(batch_input{S}.directory);
            fprintf('Analyzing subject %g \n in %s \n',S, pwd)
            
            % correct to add 0 for subject 1
            % then assume same matrix for all
            % other subjects
            try load LIMO
                if S == 1
                    disp('checking contrast validity ..')
                    for j=1:length(default.contrast)
                        C(j,:) = limo_contrast_checking(LIMO.dir, LIMO.design.X, default.contrast{j});
                        default.contrast{j} = C(j,:);
                    end
                    
                    % check overall validity
                    go = limo_contrast_checking(C,LIMO.design.X);
                    if go == 0
                        error('at least one contrast is invalid')
                    end
                end
                
                % Perform the analysis
                if isfield(LIMO,'contrast')
                    previous_con = size(LIMO.contrast,2);
                else
                    previous_con = 0;
                end
                
                load Yr; load Betas;
                
                for i=1:size(C,1)  % for each new contrast
                    
                    % update LIMO.mat with contrasts
                    switch contrast_kind
                        case'T'
                            LIMO.contrast{previous_con+i}.C = C(i,:); % if T
                            result = limo_contrast(Yr, Betas, LIMO, 0,1);
                        case'F'
                            LIMO.contrast{previous_con+i}.C = diag(C(i,:));% if F
                            LIMO.contrast{previous_con+i}.V = contrast_kind;
                            result = limo_contrast(Yr, Betas, LIMO, 1,1);
                    end
                    save LIMO LIMO
                end
                
            catch ME
                sprintf('Data not found! skipped subject %g', num2str(S))
            end
        end
        
    else % different contrast for each subject
        
        % input a contrast for each subject
        for S=1:N
            cd(batch_input{S}.directory);
            fprintf('Subject %g \n in %s \n',S, pwd)
            
            try load LIMO
                limo_review(LIMO)
                c = inputdlg(['Enter contrast subject ' num2str(S)], 'Contrasts');
                if isempty(c)
                    go = 0; return
                else
                    c = str2num(cell2mat(c));
                    C(S,:) = limo_contrast_checking(LIMO.dir, LIMO.design.X, c);
                    go = limo_contrast_checking(C(S,:),LIMO.design.X);
                    if go == 0
                        error('invalid contrast')
                    else
                        close('Review Design');
                    end
                end
            catch ME
                batch_input{S} = [];
                sprintf('Data not found! skipped subject %g', S)
            end
        end
        
        % run the contrasts
        for S=1:N
            try
                cd(batch_input{S}.directory);
                fprintf('running contrast subject %g \n in %s \n',S, pwd)
                
                % Perform the analysis
                if isfield(LIMO,'contrast')
                    previous_con = size(LIMO.contrast,2);
                else
                    previous_con = 0;
                end
                
                load Yr; load Betas;
                
                % update LIMO.mat with contrasts and create con/ess file
                switch contrast_kind
                    case'T'
                        LIMO.contrast{previous_con+1}.C = C(S,:); % if T
                        LIMO.contrast{previous_con+1}.V = contrast_kind;
                        result = limo_contrast(Yr, Betas, LIMO, 0,1);
                    case'F'
                        LIMO.contrast{previous_con+1}.C = diag(C(S,:));% if F
                        LIMO.contrast{previous_con+1}.V = contrast_kind;
                        result = limo_contrast(Yr, Betas, LIMO, 1,1);
                end
                
                save LIMO LIMO
            catch ME
                disp(' '); % for bad subjects - empty dir
            end
        end
    end
    
else  % same as above for contrasts but in the context of the model spec and estimation first
    
    for S=1:N
        cd(batch_input{S}.directory);
        fprintf('Analyzing subject %g \n',S)
        
        % 3 - make the design matrix (limo_design_matrix)
        % -----------------------------------------------
        load LIMO
        cd(LIMO.data.data_dir);
        try
            EEG=pop_loadset(LIMO.data.data);
            Y = EEG.data(:,LIMO.data.trim1:LIMO.data.trim2,:);
            clear EEG ALLCOM ALLEEG CURRENTSET CURRENTSTUDY LASTCOM STUDY
            cd (LIMO.dir)
            
            % make the design matrix
            disp('computing design matrix');
            [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,LIMO.design.nb_continuous] = limo_design_matrix(Y, LIMO,0);
            
            % update LIMO.mat
            LIMO.design.name  = 'batch processing';
            disp('design matrix done ...')
            LIMO.design.status = 'to do';
            save LIMO LIMO
            clear Y Cat Cont
            
            % 4 - Analyze data
            % -----------------
            limo_eeg(4)
            
            % 5 - run constrasts
            % ---------------------
            if isfield(default,'contrast')
                
                % correct to add 0 for subject 1
                % then assume same matrix for all
                % other subjects
                disp('checking contrast validity ..')
                if S == 1
                    for j=1:length(default.contrast)
                        C(j,:) = limo_contrast_checking(LIMO.dir, LIMO.design.X, default.contrast{j});
                        default.contrast{j} = C(j,:);
                    end
                    
                    % check overall validity
                    go = limo_contrast_checking(C,LIMO.design.X);
                    if go == 0
                        error('at least one contrast is invalid')
                    end
                    
                end
                
                % Perform the analysis
                
                if isfield(LIMO,'contrast')
                    previous_con = size(LIMO.contrast,2);
                else
                    previous_con = 0;
                end
                
                load Yr; load Betas;
                
                for i=1:size(C,1)  % for each new contrast
                    fprintf('computing contrast %g',i); disp(' ');
                    
                    % update LIMO.mat with contrasts
                    switch contrast_kind
                        case'T'
                            LIMO.contrast{previous_con+i}.C = C(i,:); % if T
                            result = limo_contrast(Yr, Betas, LIMO, 0,1);
                        case'F'
                            LIMO.contrast{previous_con+i}.C = diag(C(i,:));% if F
                            LIMO.contrast{previous_con+i}.V = contrast_kind;
                            result = limo_contrast(Yr, Betas, LIMO, 1,1);
                    end
                    save LIMO LIMO
                end
            end
            
        catch ME
            fprintf('An error occured during processing of subject %g \n %s \n', S,  ME.message)
        end
    end
end
