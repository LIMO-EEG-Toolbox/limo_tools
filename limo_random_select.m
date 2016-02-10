function filepath = limo_random_select(type,expected_chanlocs,varargin)

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
% limo_random_select(type,expected_chanlocs,'nboot',100,'tfce',1)
% limo_random_select(type,expected_chanlocs,'nboot',nbootval,'tfce',tfceval,'analysis_type','singlechan','electrode',2,'parameters',{[1:3]});
%
% INPUT
% type = 1 for a one sample t-test
% type = 2 for a two-samples t-test
% type = 3 for a paired t-test
% type = 4 for a regression
% type = 5 for an ANOVA
% expected_chanlocs: the EEGlab structure defining all electrodes
% (this file can be created with via limo_tools)
%
% Optional inputs:
%  'type'           - 'ica' or 'chan'
%  'nboot'          - the number of bootstrap to do (default = 0)
%  'tfce'           - 0/1 indicates to computes tfce or not (default = 0)
%  'analysis_type'  - |'fullchan'|'singlechan'| Type of analysis to perform. (no default)
%                     Select all channels ('fullchan') or single channel ('singlechan')
%  'limofiles'      - Cell array with the full paths to the subject file or file
%                     who contains the list of path to a group of sets. The
%                     dimensions of the cell correspond with group, factor and
%                     level respectively. (no default)
%  'electrode'      - Index of the electrode to use if 'singlechan' is
%                     selected in analysis_type.(no default)
%  'parameters'     - Cell array of parameters to be tested. 
%                     ie. {[1 2]} or {[1 2],[1 2]} in case of 2 groups. (no default) 
%  'regfile'        - Full path to Regressor file in case type = 4 (no default)
%  'folderprefix'   - Prefix for the results folder.
%  'folderpath'     - Path to save the results. Default is current
%                     directory
% Note: If the values of the parameters without default values are not
%       provided, a window will pop asking for the value.
%
% OUTPUT
% filepath - Path to the contrast result file. Mainly for EEGALB functionality to
%            allow loading test directly.
%
% Cyril Pernet - The University of Edinburgh
% Ramon Martinez-Cancino - UC San Diego
% ---------------------------------------------------------
% Copyright (C) LIMO Team 2015

%% take the inputs and load some parameters

if nargin < 2
    fprintf(2,'limo_random_select error: not enough arguments \n');
    help limo_random_select;
    return;
end

chanfile     = load(expected_chanlocs);
maxchan_indx = length(chanfile.expected_chanlocs); % Getting number of electrodes

g = finputcheck(varargin, { 'nboot'          'integer'  []                             0  ;     % Bootstrap
                            'tfce'           'integer'  [0 1]                          0  ;     % tfce
                            'analysis_type'  'string'   {'fullchan','singlechan',''}   '' ;     % Analysis Type (Full scalp or single electrode)
                            'limofiles'      'cell'     {}                             {} ;     % Path to subject file or group file Cell array with dimensions {group,,level}
                            'electrode'      'integer'  [1:maxchan_indx]               [] ;     % Electrode index
                            'regfile'        'string'   ''                             '' ;     % Path to regressor files
                            'folderprefix'   'string'   ''                             '' ;     % Prefix for folder to save 
                            'folderpath'     'string'   ''                             '' ;     % Path to folder to save
                            'type'           'string'   {'chan','ica'}             'chan' ;     % Type of measure ['ica', 'chan']
                            'parameters'     'cell'     {}                             {} ;});  % Parameters to analyze (one cell p/group)
if isstr(g), error(g); end; 
clear  chanfile maxchan_indx;

% Check Analysis Type
if isempty(g.analysis_type)
    g.analysis_type   = questdlg('Rdx option','type of analysis?','Full scalp analysis','1 electrode only','Full scalp analysis');
    if isempty(g.analysis_type), return; end
else
    analysis_type_argument = {'Full scalp analysis','1 electrode only'};
    tmpindx                = find(strcmp(g.analysis_type,{'fullchan','singlechan'}));
    g.analysis_type        = analysis_type_argument(tmpindx); clear tmpindx;
end

% check chanlocs and g.nboot
global limo
 if ~isempty(g.folderpath)
     cd(g.folderpath);
 end
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

limo.design.bootstrap = g.nboot;
limo.design.tfce = g.tfce;
limo.Level = 2;

% ----------------------------------
%%  One sample t-test and regression
% ----------------------------------
if type == 1 || type == 4

    % get files
    % ---------
    if isempty(g.limofiles)
        [Names,Paths,limo.data.data] = limo_get_files;
    % Case for path to the files
    elseif size(g.limofiles{1},1) == 1 
        [Names,Paths,limo.data.data] = limo_get_files([],[],[],g.limofiles{1});
    % Case when all paths are provided
    elseif size(g.limofiles{1},1) > 1   
        [Names,Paths,limo.data.data] = breaklimofiles(g.limofiles{1}); 
    end
    
    if isempty(Names)
        disp('no files selected')
        return
    end
    limo.data.data_dir = Paths;
    
    % check type of files and returns which beta param to test
    % -------------------------------------------------------
    if isempty(g.parameters)
        parameters = check_files(Names,1);
    else
        parameters = check_files(Names,1,g.parameters{1});
    end
    if isempty(parameters)
        errordlg('file selection failed, only Beta and Con files are supported','Selection error'); return
    end
    
    % match frames
    % ------------
    [first_frame,last_frame,subj_chanlocs,channeighbstructmat] = match_frames(Paths);
    
   
    % match electrodes
    % --------------
    if strcmp(g.analysis_type,'1 electrode only')
        if isempty(g.electrode)
            electrode = inputdlg('which electrode to analyse [?]','Electrode option'); % can be 1 nb or a vector of electrodes (electrode optimized analysis)
        else
            electrode = {num2str(g.electrode)};
        end
        
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
    else % Full scalp
        if type == 1
            limo.design.name = 'one sample t-test all electrodes';
            limo.data.chanlocs = expected_chanlocs;
        else
            limo.design.name = 'regression analysis all electrodes';
            limo.data.chanlocs = expected_chanlocs;
        end
        limo.design.electrode = [];
    end
    
    % get data for all parameters 
    % -----------------------------
    disp('gathering data ...'); index = 1;
    for i=1:size(Paths,2) % for each subject
        load(limo.data.data{i});
        try
            tmp = eval(str2mat(Names{1}(1:end-4))); % load Betas values
        catch
            tmp = eval(str2mat(Names{1}(1:end-6)));
            tmp = squeeze(tmp(:,:,1));
        end
        
        % get indices to trim data
        if strcmp(limo.Analysis,'Time-Frequency')
            begins_at = fliplr((max(first_frame) - first_frame(i,:) + 1)); % returns time/freq/or freq-time
            ends_at(1) = size(tmp,2) - (last_frame(i,2) - min(last_frame(:,2)));
            ends_at(2) = size(tmp,3) - (last_frame(i,1) - min(last_frame(:,1)));
        else
            begins_at = max(first_frame) - first_frame(i) + 1;
            ends_at = size(tmp,2) - (last_frame(i) - min(last_frame));
        end
        
        % data dim [electrode, freq/time, param, nb subjects]
        if strcmp(g.analysis_type,'Full scalp analysis') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
            if strcmp(limo.Analysis,'Time-Frequency')
                data(:,:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
            else
                data(:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
            end
            index = index + 1; removed(i) = 0;
        
        elseif strcmp(g.analysis_type,'1 electrode only') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
            if  length(limo.design.electrode) == 1;
                if strcmp(limo.Analysis,'Time-Frequency')
                    data(1,:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                else
                    data(1,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                end
                index = index + 1; removed(i) = 0;
            else % use multiple single electrodes
                out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                if strcmp(limo.Analysis,'Time-Frequency')
                    data(1,:,:,:,index) = out(i,:,:,:);
                else
                    data(1,:,:,index) = out(i,:,:); % matches the expected chanloc of the subject
                end
                index = index +1; removed(i) = 0;
            end
        else
            fprintf('subject %g discarded, channel description and data size don''t match',i);
            removed(i) = 1; disp(' ')
        end
        clear tmp
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
            if length(parameters) > 1 && i == parameters(1)
                foldername = 'parameter_%g';
                if ~isempty(g.folderprefix), foldername = [g.folderprefix foldername]; end
                dir_name = sprintf(foldername,i);
                mkdir(dir_name); cd(dir_name);
            elseif length(parameters) > 1 && i ~= parameters(1)
                cd ..
                foldername = 'parameter_%g';
                if ~isempty(g.folderprefix), foldername = [g.folderprefix foldername]; end
                dir_name = sprintf(foldername,i);
                mkdir(dir_name); cd(dir_name);
            end
            LIMO.dir = pwd;
            
            if strcmp(g.analysis_type,'1 electrode only') && size(data,1) == 1
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    tmp = squeeze(data(:,:,:,i,:));
                    tmp_data = NaN(1,size(tmp,1),size(tmp,2),size(tmp,3)); % add dim 1 = 1 electrode
                    tmp_data(1,:,:,:) = tmp; clear tmp;
                else
                    tmp = squeeze(data(:,:,i,:));
                    tmp_data = NaN(1,size(tmp,1),size(tmp,2)); % add dim 1 = 1 electrode
                    tmp_data(1,:,:) = tmp; clear tmp;
                end
            else
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    tmp_data = squeeze(data(:,:,:,i,:));
                else
                    tmp_data = squeeze(data(:,:,i,:));
                end
            end
            
            if strcmp(LIMO.Analysis,'Time-Frequency')
                LIMO.data.size3D = [size(tmp_data,1) size(tmp_data,2)*size(tmp_data,3) size(tmp_data,4)];
                LIMO.data.size4D = [size(tmp_data,1) size(tmp_data,2) size(tmp_data,3) size(tmp_data,4)];
            end
            
            LIMO.design.method = 'Trimmed means'; save LIMO LIMO
            Yr = tmp_data; save Yr Yr, clear Yr % just to be consistent with name
            tmpname = limo_random_robust(type,tmp_data,i,g.nboot,g.tfce);
            if nargout ~= 0, filepath{i} = tmpname; end
        end
        
        
        % regression
        % -------------
        
    elseif type == 4
        cd(limo.dir)
        if isempty(g.regfile)
            [FileName,PathName,FilterIndex]=uigetfile('*.txt;*.mat','select regressor file');
        elseif ~isempty(g.regfile) && exist(g.regfile,'file')
            [PathName,nametmp,exttmp] = fileparts(g.regfile);
            FileName = [nametmp exttmp]; 
            clear nametmp exttmp;
        else
            error('limo_random_select error: Provide a valid regressor file');
        end
        if FilterIndex == 0
            return
        elseif strcmp(FileName(end-3:end),'.txt')
            cd(PathName)
            X = load(FileName);
        elseif strcmp(FileName(end-3:end),'.mat')
            cd(PathName)
            % load(FileName);
            % X = eval(FileName(1:end-4));
            X = load(FileName); 
            X = getfield(X,cell2mat(fieldnames(X)));
        end
        disp('Regressor(s) loaded');
        
        % check size and orientation
        try
            if strcmp(limo.Analysis,'Time-Frequency')
                N = size(data,5);
            else
                N = size(data,4);
            end
            
            if size(X,2) == N || size(X,2) == size(Paths,2)
                disp('X has been transposed'); X = X';
            end

            if size(X,1) > N;
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
                       
            if size(X,2)==1 && g.nboot < 599; 
                limo.design.bootstrap = 599;
                disp('nb of bootstrap adjusted to 599 for a simple regression'); 
            end
                        
        catch ME
            errordlg('covariate error - make sure data are in lines or columns','Covariate error'); 
            return
        end
        
        LIMO = limo;
        clear Betas Names Paths channeighbstructmat expected_chanlocs limo subj_chanlocs
        if parameters == 0; parameters = [1:size(data,3)]; end
        
        for i=parameters
            cd(LIMO.dir);
            if length(parameters) > 1 && i == parameters(1)
                foldername = 'parameter_%g';
                if ~isempty(g.folderprefix), foldername = [g.folderprefix foldername]; end
                dir_name = sprintf(foldername,i);
                mkdir(dir_name); cd(dir_name);
            elseif length(parameters) > 1 && i ~= parameters(1)
                cd ..
                foldername = 'parameter_%g';
                if ~isempty(g.folderprefix), foldername = [g.folderprefix foldername]; end
                dir_name = sprintf(foldername,i);
                mkdir(dir_name); cd(dir_name);
            end
            
            if strcmp(g.analysis_type,'1 electrode only')
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    tmp = squeeze(data(:,:,:,i,:));
                    tmp_data = NaN(1,size(tmp,1),size(tmp,2),size(tmp,3)); % add dim 1 = 1 electrode
                    tmp_data(1,:,:,:) = tmp; clear tmp;
                else
                    tmp = squeeze(data(:,:,i,:));
                    tmp_data = NaN(1,size(tmp,1),size(tmp,2)); % add dim 1 = 1 electrode
                    tmp_data(1,:,:) = tmp; clear tmp;
                end
            else
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    tmp_data = squeeze(data(:,:,:,i,:));
                else
                    tmp_data = squeeze(data(:,:,i,:));
                end
            end
            
            if strcmp(LIMO.Analysis,'Time-Frequency')
                LIMO.data.size3D = [size(tmp_data,1) size(tmp_data,2)*size(tmp_data,3) size(tmp_data,4)];
                LIMO.data.size4D = [size(tmp_data,1) size(tmp_data,2) size(tmp_data,3) size(tmp_data,4)];
            end

            % compute
            save LIMO LIMO; clear LIMO ;
            tmpname = limo_random_robust(type,tmp_data,X,i,g.nboot,g.tfce);
            if nargout ~= 0, filepath{i} = tmpname; end
        end
    end

    % ------------------------------
    %%  Two samples t-test
    % -----------------------------
elseif type == 2

    limo.design.X = [];
    if strcmp(g.analysis_type,'Full scalp analysis') || strcmp(g.analysis_type,'1 electrode only')

        N = 0;
        for gp = 1:2
            if isempty(g.limofiles)
                [Names{gp},Paths{gp},limo.data.data{gp}] = limo_get_files([' gp' num2str(gp)]);
            % Case for path to the files
            elseif size(g.limofiles{gp},1) == 1
                [Names{gp},Paths{gp},limo.data.data{gp}] = limo_get_files([],[],[],g.limofiles{gp});
            % Case when all paths are provided
            elseif size(g.limofiles{gp},1) > 1
                [Names{gp},Paths{gp},limo.data.data{gp}] = breaklimofiles(g.limofiles{gp});
            end
            if isempty(Names{gp})
                return
            end
            limo.data.data_dir{gp} = Paths{gp};
            N = N + size(Names{gp},2);
            if isempty(g.parameters)
                parameters(:,gp) = check_files(Names{gp},1);
            else
                parameters(:,gp) = check_files(Names{gp},1,g.parameters{gp});
            end
        end

        if ~exist('parameters','var')
            disp('selection aborded')
            return
        end
        
        % check type of files and returns which beta param to test
        % -------------------------------------------------------
        for a=1:size(Paths,2)
            for b=1:size(Paths{a},2)
                cd (cell2mat(Paths{a}(b))); load LIMO
                if max(parameters(a)) > size(LIMO.design.X,2)
                    errordlg('invalid parameter(s)','Parameters error'); return
                end
            end
        end
        cd(limo.dir);
        
        % match frames
        % ------------
        [first_frame,last_frame,subj_chanlocs] = match_frames(Paths);


        % match electrodes
        % --------------
        if strcmp(g.analysis_type,'1 electrode only')
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
        for igp = 1:gp
            index = 1;
            for i=1:size(Paths{igp},2) % for each subject per group
                load(cell2mat(limo.data.data{igp}(i)));
                name = str2mat(cell2mat(Names{igp}(i)));
                try
                    tmp = eval(name(1:end-4));
                catch
                    tmp = eval(name(1:end-6));
                    tmp = squeeze(tmp(:,:,1));
                end
                
                % get indices to trim data
                if strcmp(limo.Analysis,'Time-Frequency')
                    begins_at = fliplr((max(first_frame) - first_frame(subject_nb,:) + 1)); % returns time/freq/or freq-time
                    ends_at(1) = size(tmp,2) - (last_frame(subject_nb,2) - min(last_frame(:,2)));
                    ends_at(2) = size(tmp,3) - (last_frame(subject_nb,1) - min(last_frame(:,1)));
                else
                    begins_at = max(first_frame) - first_frame(subject_nb) + 1;
                    ends_at = size(tmp,2) - (last_frame(subject_nb) - min(last_frame));
                end
                
                if strcmp(g.analysis_type,'Full scalp analysis') && size(subj_chanlocs(subject_nb).chanlocs,2) == size(tmp,1)
                    if strcmp(limo.Analysis,'Time-Frequency')
                        tmp_data(:,:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                        index = index + 1;
                    else
                        tmp_data(:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                        index = index + 1;
                    end
                    
                elseif strcmp(g.analysis_type,'1 electrode only') && size(subj_chanlocs(subject_nb).chanlocs,2) == size(tmp,1)
                    if size(limo.design.electrode,2) == 1;
                        if strcmp(limo.Analysis,'Time-Frequency')
                            tmp_data(1,:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                            index = index + 1;
                        else
                            tmp_data(1,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                            index = index + 1;
                        end
                    else
                        out = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                        if strcmp(limo.Analysis,'Time-Frequency')
                            tmp_data(1,:,:,:,index) = out(subject_nb,:,:,:); % matches the expected chanloc of the subject
                            index = index +1;
                        else
                            tmp_data(1,:,:,index) = out(subject_nb,:,:); % matches the expected chanloc of the subject
                            index = index +1;
                        end
                    end
                else
                    fprintf('subject %g of group %g discarded, channel description and data size don''t match',i, igp); disp(' ')
                end
                
                clear tmp
                subject_nb = subject_nb + 1;
            end
            
            data{igp} = tmp_data;
            clear tmp tmp_data
        end
    end


    % compute
    % --------
    LIMO = limo; cd(limo.dir);  i=parameters;
    % free some memory
    clear Betas Names Paths channeighbstructmat expected_chanlocs limo subj_chanlocs
     
    if strcmp(LIMO.Analysis,'Time-Frequency')
        if strcmp(g.analysis_type,'1 electrode only')
            tmp = squeeze(data{1}(:,:,:,i(1),:));
            tmp_data1 = ones(1,size(tmp,1),size(tmp,2),size(tmp,3)); % add dim 1 = 1 electrode
            tmp_data1(1,:,:,:) = tmp; clear tmp
            tmp = squeeze(data{2}(:,:,:,i,:));
            tmp_data2 = ones(1,size(tmp,1),size(tmp,2),size(tmp,3));
            tmp_data2(1,:,:,:) = tmp; clear tmp
        else
            tmp_data1 = squeeze(data{1}(:,:,:,i(1),:));
            tmp_data2 = squeeze(data{2}(:,:,:,i(2),:));
        end
        
        if size(tmp_data1,1) ~= size(tmp_data2,1) || size(tmp_data1,2) ~= size(tmp_data2,2) || size(tmp_data1,3) ~= size(tmp_data2,3)
            errordlg('file selection is corrupted, data sizes don''t match');
            return
        end
        
    else
        if strcmp(g.analysis_type,'1 electrode only')
            tmp = squeeze(data{1}(:,:,i,:));
            tmp_data1 = ones(1,size(tmp,1),size(tmp,2)); % add dim 1 = 1 electrode
            tmp_data1(1,:,:) = tmp; clear tmp
            tmp = squeeze(data{2}(:,:,i,:));
            tmp_data2 = ones(1,size(tmp,1),size(tmp,2));
            tmp_data2(1,:,:) = tmp; clear tmp
        else
            tmp_data1 = squeeze(data{1}(:,:,i(1),:));
            tmp_data2 = squeeze(data{2}(:,:,i(2),:));
        end
        
        if size(tmp_data1,1) ~= size(tmp_data2,1) || size(tmp_data1,2) ~= size(tmp_data2,2)
            errordlg('file selection is corrupted, data sizes don''t match');
            return
        end
    end
    
    if strcmp(LIMO.Analysis,'Time-Frequency')
        LIMO.data.size3D = [size(tmp_data1,1) size(tmp_data1,2)*size(tmp_data1,3) 5];
        LIMO.data.size4D = [size(tmp_data1,1) size(tmp_data1,2) size(tmp_data1,3) 5];
    end
    
    Y1r = tmp_data1; save Y1r Y1r, clear Y1r
    Y2r = tmp_data2; save Y2r Y2r, clear Y2r
    LIMO.design.method = 'Yuen t-test (trimmed means)'; save LIMO LIMO
    tmpname = limo_random_robust(type,tmp_data1,tmp_data2,i,g.nboot,g.tfce);
    if nargout ~= 0, filepath = tmpname; end
    delete data.mat


    % ------------------------------
    %%  Paired t-test
    % -----------------------------
elseif type == 3

    limo.design.X = [];
    if strcmp(g.analysis_type,'Full scalp analysis') || strcmp(g.analysis_type,'1 electrode only')
        if isempty(g.limofiles)
            [Names,Paths,limo.data.data] = limo_get_files;
        % Case for path to the files
        elseif size(g.limofiles{1},1) == 1
            [Names,Paths,limo.data.data] = limo_get_files([],[],[],g.limofiles{1});
        % Case when all paths are provided
        elseif size(g.limofiles{1},1) > 1
            [Names,Paths,limo.data.data] = breaklimofiles(g.limofiles{1});
        end
        if isempty(Names)
            return
        end
        limo.data.data_dir = Paths;
        N = size(Names,2);


        % check type of files and returns which beta param to test
        % -------------------------------------------------------
        if isempty(g.parameters)
            parameters = check_files(Names,1);
        else
        parameters = check_files(Names,1,g.parameters{1});
        end
        if size(parameters,2) == 1  % was not beta files, ie was con files
            n = Names; clear Names; Names{1} = n; clear n;
            p = Paths; clear Paths; Paths{1} = p; clear p;
            d = limo.data.data; clear limo.data.data; limo.data.data{1} = d; clear d;
            d = limo.data.data_dir; clear limo.data.data_dir; limo.data.data_dir{1} = d; clear d;
            
            if isempty(g.limofiles)
                [Names{2},Paths{2},limo.data.data{2}] = limo_get_files([' gp2']);
            % Case for path to the files
            elseif size(g.limofiles{2},1) == 1
                [Names{2},Paths{2},limo.data.data{2}] = limo_get_files([],[],[],g.limofiles{2});
            % Case when all paths are provided
            elseif size(g.limofiles{2},1) > 1
                [Names{2},Paths{2},limo.data.data{2}] = breaklimofiles(g.limofiles{2});
            end
            
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
            for s = 1:size(Paths,2)
                cd (cell2mat(Paths(s))); load LIMO
                if max(parameters) > size(LIMO.design.X,2)
                    errordlg('invalid parameter(s)','Paired t-test error'); return
                end
            end
            cd(limo.dir);
        end

        
        % match frames
        % ------------
        [first_frame,last_frame,subj_chanlocs] = match_frames(Paths);


        % match electrodes
        % --------------
        if strcmp(g.analysis_type,'1 electrode only')
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

        if size(parameters,2) == 1 % groups and cons

            subject_nb = 1;
            for g = 1:size(Paths,2)
                index = 1;
                for i=1:size(Paths{g},2) % for each subject per group
                    load(cell2mat(limo.data.data{g}(i)));
                    name = str2mat(cell2mat(Names{g}(i)));
                    if strcmp(name,'Betas.mat')
                        tmp = eval(name(1:end-4));
                    else
                        tmp = eval(name(1:end-6));
                    end
                    
                    % get indices to trim data
                    if strcmp(limo.Analysis,'Time-Frequency')
                        tmp = squeeze(tmp(:,:,:,1));
                        begins_at = fliplr((max(first_frame) - first_frame(subject_nb,:) + 1)); % returns time/freq/or freq-time
                        ends_at(1) = size(tmp,2) - (last_frame(subject_nb,2) - min(last_frame(:,2)));
                        ends_at(2) = size(tmp,3) - (last_frame(subject_nb,1) - min(last_frame(:,1)));
                    else
                        tmp = squeeze(tmp(:,:,1));
                        begins_at = max(first_frame) - first_frame(subject_nb) + 1;
                        ends_at = size(tmp,2) - (last_frame(subject_nb) - min(last_frame));
                    end
                    
                    if strcmp(g.analysis_type,'Full scalp analysis') && size(subj_chanlocs(subject_nb).chanlocs,2) == size(tmp,1)
                        if strcmp(limo.Analysis,'Time-Frequency')
                            tmp_data(:,:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                            index = index + 1;
                        else
                            tmp_data(:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                            index = index + 1;
                        end
                    elseif strcmp(g.analysis_type,'1 electrode only') && size(subj_chanlocs(subject_nb).chanlocs,2) == size(tmp,1)
                        if strcmp(limo.Analysis,'Time-Frequency')
                            if size(limo.design.electrode,2) == 1;
                                tmp_data(1,:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                                index = index + 1;
                            else
                                out = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                                tmp_data(1,:,:,:,index) = out(subject_nb,:,:); % matches the expected chanloc of the subject
                                index = index +1;
                            end
                        else
                            if size(limo.design.electrode,2) == 1;
                                tmp_data(1,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                                index = index + 1;
                            else
                                out = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                                tmp_data(1,:,:,index) = out(subject_nb,:,:); % matches the expected chanloc of the subject
                                index = index +1;
                            end
                        end
                    else
                        fprintf('subject %g of group %g discarded, channel description and data size don''t match',i, g); disp(' ')
                    end
                    clear tmp
                    subject_nb = subject_nb + 1;
                end

                if strcmp(g.analysis_type,'1 electrode only') && size(tmp_data,2) == 1
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
                
                % get indices to trim data
                if strcmp(limo.Analysis,'Time-Frequency')
                    begins_at = fliplr((max(first_frame) - first_frame(i,:) + 1)); % returns time/freq/or freq-time
                    ends_at(1) = size(tmp,2) - (last_frame(i,2) - min(last_frame(:,2)));
                    ends_at(2) = size(tmp,3) - (last_frame(i,1) - min(last_frame(:,1)));
                else
                    begins_at = max(first_frame) - first_frame(i) + 1;
                    ends_at = size(tmp,2) - (last_frame(i) - min(last_frame));
                end
                
                if strcmp(g.analysis_type,'Full scalp analysis') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
                    if strcmp(limo.Analysis,'Time-Frequency')
                        data(:,:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                        index = index + 1;
                    else
                        data(:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                        index = index + 1;
                    end
                elseif strcmp(g.analysis_type,'1 electrode only') && size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)
                    if size(limo.design.electrode,2) == 1;
                        if strcmp(limo.Analysis,'Time-Frequency')
                            data(1,:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                            index = index + 1;
                        else
                            data(1,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                            index = index + 1;
                        end
                    else
                        out = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                        if strcmp(limo.Analysis,'Time-Frequency')
                            data(1,:,:,:,index) = out(i,:,:,:); % matches the expected chanloc of the subject
                            index = index +1;
                        else
                            data(1,:,:,index) = out(i,:,:); % matches the expected chanloc of the subject
                            index = index +1;
                        end
                    end
                else
                    fprintf('subject %g discarded, channel description and data size don''t match',i); disp(' ')
                end
                clear tmp
            end
        end
    end


    % compute
    % --------
    if strcmp(limo.Analysis,'Time-Frequency')
        if strcmp(g.analysis_type,'1 electrode only')
            if size(parameters,2) == 2 % beta files
                tmp = squeeze(data(:,:,:,parameters(1),:));
                tmp_data1 = ones(1,size(tmp,1),size(tmp,2),size(tmp,3)); % add dim 1 = 1 electrode
                tmp_data1(1,:,:,:) = tmp; clear tmp
                tmp = squeeze(data(:,:,:,parameters(2),:));
                tmp_data2 = ones(1,size(tmp,1),size(tmp,2),size(tmp,3));
                tmp_data2(1,:,:,:) = tmp; clear tmp
            else % con files
                tmp = squeeze(data{1}(:,:,:,:,:));
                tmp_data1 = ones(1,size(tmp,1),size(tmp,2),size(tmp,3)); % add dim 1 = 1 electrode
                tmp_data1(1,:,:,:) = tmp; clear tmp
                tmp = squeeze(data{2}(:,:,:,:,:));
                tmp_data2 = ones(1,size(tmp,1),size(tmp,2),size(tmp,3));
                tmp_data2(1,:,:,:) = tmp; clear tmp
            end
        else
            if size(parameters,2) == 2 % beta files
                tmp_data1 = squeeze(data(:,:,:,parameters(1),:));
                tmp_data2 = squeeze(data(:,:,:,parameters(2),:));
            else % con files
                tmp_data1(:,:,:,1,:) = squeeze(data{1}(:,:,:,:));
                tmp_data2(:,:,:,1,:) = squeeze(data{2}(:,:,:,:));
            end
        end
    else
        if strcmp(g.analysis_type,'1 electrode only')
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
    end
    
    LIMO = limo; cd(limo.dir); 
    LIMO.design.method = 'Yuen t-test (Trimmed means)'; 
    if strcmp(LIMO.Analysis,'Time-Frequency')
        LIMO.data.size3D = [size(tmp_data1,1) size(tmp_data1,2)*size(tmp_data1,3) 5];
        LIMO.data.size4D = [size(tmp_data1,1) size(tmp_data1,2) size(tmp_data1,3) 5];
    end
    save LIMO LIMO
    
    Y1r = tmp_data1; save Y1r Y1r, clear Y1r
    Y2r = tmp_data2; save Y2r Y2r, clear Y2r
    tmpname = limo_random_robust(type,tmp_data1,tmp_data2,parameters,g.nboot,g.tfce);
    if nargout ~= 0, filepath = tmpname; end
    
    % -----------------------------------
    %%  Various sorts of ANOVAs/ANCOVAs
    % -----------------------------------
elseif type == 5

    % 1st know what design it is
    % --------------------------
    
    answer = questdlg('What ANOVA model do you want to run?', 'Model selection', 'Repeated Measures', ...
        'N-Ways','ANCOVA','Repeated Measures');
    
    % get some nice comment in LIMO.mat
    if strcmp(g.analysis_type,'Full scalp analysis')
        if strcmp(answer,'Repeated Measures')
            limo.design.name = 'Repeated measures ANOVA all electrodes';
        elseif strcmp(answer,'N-Ways')
            limo.design.name = 'N-ways ANOVA all electrodes';
        else
            limo.design.name = 'ANCOVA all electrodes';
        end

    elseif strcmp(g.analysis_type,'1 electrode only')
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
        elseif sum(gp_nb <= 1)
            errordlg('at least 2 groups are expected for an ANOVA')
            return
        end
        
        % select data per gp / conditions
        % ---------------------------------
        a = questdlg('load con files or beta file','ANOVA loading files','con','beta','beta');
        N = 0; cell_nb = 1;
        if strcmp(a,'beta') % beta files
            for i=1:prod(gp_nb)
                if isempty(g.limofiles)
                    [Names{cell_nb},Paths{cell_nb},limo.data.data{cell_nb}] = limo_get_files([' beta file gp ',num2str(i)]);
                % Case for path to the files
                elseif size(g.limofiles{i},1) == 1
                    [Names{cell_nb},Paths{cell_nb},limo.data.data{cell_nb}] = limo_get_files([],[],[],g.limofiles{i});
                % Case when all paths are provided
                elseif size(g.limofiles{i},1) > 1
                    [Names{cell_nb},Paths{cell_nb},limo.data.data{cell_nb}] = breaklimofiles(g.limofiles{i});
                end
                if isempty(Names{cell_nb}); return; end
                if isempty(g.parameters)
                    parameters(:,i) = check_files(Names,1);
                else
                    parameters(:,i) = check_files(Names,1,g.parameters{i});
                end
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
                if isempty(g.limofiles)
                    [Names{cell_nb},Paths{cell_nb},limo.data.data{cell_nb}] = limo_get_files([' con file gp ',num2str(i)]);
                % Case for path to the files
                elseif size(g.limofiles{i},1) == 1
                    [Names{cell_nb},Paths{cell_nb},limo.data.data{cell_nb}] = limo_get_files([],[],[],g.limofiles{i});
                % Case when all paths are provided
                elseif size(g.limofiles{i},1) > 1
                    [Names{cell_nb},Paths{cell_nb},limo.data.data{cell_nb}] = breaklimofiles(g.limofiles{i});
                end
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
        if strcmp(g.analysis_type,'1 electrode only')
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
                
                % get indices to trim data
                if strcmp(limo.Analysis,'Time-Frequency')
                    begins_at = fliplr((max(first_frame) - first_frame(subject_index,:) + 1)); % returns time/freq/or freq-time
                    ends_at(1) = size(tmp,2) - (last_frame(subject_index,2) - min(last_frame(:,2)));
                    ends_at(2) = size(tmp,3) - (last_frame(subject_index,1) - min(last_frame(:,1)));
                else
                    begins_at = max(first_frame) - first_frame(subject_index) + 1;
                    ends_at = size(tmp,2) - (last_frame(subject_index) - min(last_frame));
                end
                
                % data are of dim size(expected_chanlocs,2), latter start/earlier stop across subjects, parameters, nb of subjects
                if strcmp(g.analysis_type,'Full scalp analysis') && size(subj_chanlocs(subject_index).chanlocs,2) == size(tmp,1)
                    if strcmp(limo.Analysis,'Time-Frequency')
                        data(:,:,:,:,matrix_index) = limo_match_elec(subj_chanlocs(subject_index).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                        matrix_index = matrix_index+1; removed{h}(i) = 0; nb_subjects(h) = nb_subjects(h)+1;
                    else
                        data(:,:,:,matrix_index) = limo_match_elec(subj_chanlocs(subject_index).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                        matrix_index = matrix_index+1; removed{h}(i) = 0; nb_subjects(h) = nb_subjects(h)+1;
                    end
                elseif strcmp(g.analysis_type,'1 electrode only') && size(subj_chanlocs(subject_index).chanlocs,2) == size(tmp,1)
                    if strcmp(limo.Analysis,'Time-Frequency')
                        if size(limo.design.electrode,2) == 1;
                            data(1,:,:,:,matrix_index) = limo_match_elec(subj_chanlocs(subject_index).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                        else
                            out = limo_match_elec(subj_chanlocs(subject_index).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, i.e. across subjects
                            data(1,:,:,:,matrix_index) = out(i,:,:); % matches the expected chanloc of the subject
                        end
                        matrix_index = matrix_index+1; removed{h}(i) = 0; nb_subjects(h) = nb_subjects(h)+1;
                    else
                        if size(limo.design.electrode,2) == 1;
                            data(1,:,:,matrix_index) = limo_match_elec(subj_chanlocs(subject_index).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                        else
                            out = limo_match_elec(subj_chanlocs(subject_index).chanlocs,expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, i.e. across subjects
                            data(1,:,:,matrix_index) = out(i,:,:); % matches the expected chanloc of the subject
                        end
                        matrix_index = matrix_index+1; removed{h}(i) = 0; nb_subjects(h) = nb_subjects(h)+1;
                    end
                else
                    fprintf('subject %g gp %g discarded, channel description and data size don''t match',i,h); disp(' ')
                    removed{h}(i) = 1;
                end
                clear tmp
                subject_index = subject_index+1;
            end
        end
        cd(limo.dir);
        
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
                if size(X,1) ~= size(data,ndims(data));
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
            Cont = [];
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

        
        % for N ways ANOVA we pass the data and X, so we simply stack
        % the data on top of each other per gp - this allows using betas
        % and still compare possibly different parameters
        if strcmp(limo.Analysis,'Time-Frequency')
            tmp_data = NaN(size(data,1),size(data,2),size(data,3),size(data,5));
            for i=1:prod(gp_nb)
                if i==1; from = 1; to = nb_subjects(i);
                else from = from+nb_subjects(i); to = to+nb_subjects(i); end
                current_param = parameters(i); % select only relevant parameters
                tmp_data(:,:,:,from:to) = squeeze(data(:,:,:,current_param,from:to));
            end
        else
            tmp_data = NaN(size(data,1),size(data,2),size(data,4));
            for i=1:prod(gp_nb)
                if i==1; from = 1; to = nb_subjects(i);
                else from = from+nb_subjects(i); to = to+nb_subjects(i); end
                current_param = parameters(i); % select only relevant parameters
                tmp_data(:,:,from:to) = squeeze(data(:,:,current_param,from:to));
            end
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
        
        LIMO = limo; cd(limo.dir);
        if size(Cat,2) == 1 && isempty(Cont)
            LIMO.design.method = 'Robust ANOVA (Trimmed means)';
        else
            LIMO.design.method = 'Robust ANOVA (Iterative Reweighted Least Square)';
        end
                
        if strcmp(LIMO.Analysis,'Time-Frequency')
            LIMO.data.size3D = [size(tmp_data,1) size(tmp_data,2)*size(tmp_data,3) size(tmp_data,4)];
            LIMO.data.size4D = [size(tmp_data,1) size(tmp_data,2) size(tmp_data,3) size(tmp_data,4)];
        end
        save LIMO LIMO
        
        % clear some memory
        clear LIMO limo Names Paths data channeighbstructmat expected_chanlocs subj_chanlocs
        
        % do the analysis
        Yr = tmp_data; clear tmp_data; save Yr Yr
        if isempty(Cont); Cont = 0; end
        tmpname = limo_random_robust(type,Yr,Cat,Cont,g.nboot,g.tfce);
        if nargout ~= 0, filepath = tmpname; end
        
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
                if isempty(g.limofiles)
                    [Names{cell_nb},Paths{cell_nb},limo.data.data{cell_nb}] = limo_get_files([' beta file gp ',num2str(i)]);
                % Case for path to the files
                elseif size(g.limofiles{i},1) == 1
                    [Names{cell_nb},Paths{cell_nb},limo.data.data{cell_nb}] = limo_get_files([],[],[],g.limofiles{i});
                % Case when all paths are provided
                elseif size(g.limofiles{i},1) > 1
                    [Names{cell_nb},Paths{cell_nb},limo.data.data{cell_nb}] = breaklimofiles(g.limofiles{i});
                end
                if isempty(Names{cell_nb}); return; end
                if isempty(g.parameters)
                    parameters(:,i) = check_files(Names,1);
                else
                    parameters(:,i) = check_files(Names,1,g.parameters{i});
                end
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
                        if isempty(g.limofiles)
                            [names{k},paths{k},full_names{k}] = limo_get_files([' gp ',num2str(i),' factor ',num2str(j),' level ',num2str(k)]);
                        % Case for path to the files
                        elseif size(g.limofiles{num2str(i),num2str(j),num2str(k)},1) == 1
                            [names{k},paths{k},full_names{k}] = limo_get_files([],[],[],g.limofiles{num2str(i),num2str(j),num2str(k)});
                        % Case when all paths are provided
                        elseif size(g.limofiles{num2str(i),num2str(j),num2str(k)},1) > 1
                            [names{k},paths{k},full_names{k}] = breaklimofiles(g.limofiles{num2str(i),num2str(j),num2str(k)});
                        end
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
        if strcmp(g.analysis_type,'1 electrode only')
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
                if strcmp(g.analysis_type,'Full scalp analysis') && size(subj_chanlocs(subject_index).chanlocs,2) == size(tmp,1)
                    data(:,:,:,matrix_index) = limo_match_elec(subj_chanlocs(subject_index).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
                    matrix_index = matrix_index+1; removed{h}(i) = 0; nb_subjects(h) = nb_subjects(h)+1;
                elseif strcmp(g.analysis_type,'1 electrode only') && size(subj_chanlocs(subject_index).chanlocs,2) == size(tmp,1)
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
        tmpname = limo_random_robust(type+1,tmp_data,gp,factor_nb,g.nboot,g.tfce);
        if nargout ~= 0, filepath = tmpname; end
    end
end % closes type


end % closes the function


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                                   ROUTINES
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% file checking subfunction
function parameters = check_files(Names,gp,parameters)
% after selecting file, check they are all the same type
if nargin < 3
    parameters = [];
end

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
    elseif (isempty(is_beta)) == 0 && sum(is_beta) == size(Names,2) && nargout ~= 0
        if isempty(parameters)
            parameters = eval(cell2mat(inputdlg('which parameters to test e.g [1:3]','parameters option')));
        end
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


%% frame matching subfunction
function [first_frame,last_frame,subj_chanlocs,channeighbstructmat] = match_frames(Paths)

% once we have all the files, we need to collect
% information to match the frames across subjects
% first_frame last _frame returns the beginning and
% end in terms of indices ; for time-frequency these
% are vectors with time st then frequency
% the limo structure is also updated the reflect the
% smallest interval(s) across subjects, which is used
% fore the second leve analysis
 

global limo
channeighbstructmat = []; ME = [];

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
    subj_chanlocs(i).chanlocs = LIMO.data.chanlocs;
    try
       channeighbstructmat = LIMO.data.channeighbstructmat;
    catch ME
    end
       
    if strcmp(Analysis,'Time-Frequency')
        first_frame(i,1)            = LIMO.data.trim1;
        last_frame(i,1)             = LIMO.data.trim2;
        start(i,1)                  = LIMO.data.start;
        stop(i,1)                   = LIMO.data.end;
        
        first_frame(i,2)            = LIMO.data.trim_low_f;
        last_frame(i,2)             = LIMO.data.trim_high_f;
        start(i,2)                  = LIMO.data.tf_freqs(1);
        stop(i,2)                   = LIMO.data.tf_freqs(end);
        
        tf_times{i}(1,:)            = LIMO.data.tf_times;
        tf_freqs{i}(1,:)            = LIMO.data.tf_freqs;
    else
        first_frame(i)              = LIMO.data.trim1;
        last_frame(i)               = LIMO.data.trim2;
        start(i)                    = LIMO.data.start;
        stop(i)                     = LIMO.data.end;
        
        if strcmp(Analysis,'Frequency')
            freqlist{i}(1,:)        = LIMO.data.freqlist;
        end
    end
end

% quick check things are ok
if ~isempty(ME) && isempty(limo.data.neighbouring_matrix)
    error(sprintf('some subject(s) have a different channel structure \nplease load an expected chanloc when choosing a test'));                          
end

if (sum(sampling_rate == sampling_rate(1))) ~= length(sampling_rate)
    error('data have different sampling rates')
end

% match and return into limo - the temp structure passed as global
limo.Analysis = Analysis;
limo.data.sampling_rate = sampling_rate(1);

% we need 1) to find the highest start in time and freq 2) the lowest end
% in time and freq and 3) match that on freqlist or tf_times/tf_freqs

[v,c] = max(first_frame);
if strcmp(Analysis,'Time-Frequency')
    limo.data.trim1 = v(1);
    limo.data.start = start(c(1),1);
    limo.data.trim_low_f = v(2);
    limo.data.low_f = start(c(2),2);
    % loop to adjust each subject list
    for i=1:size(Paths,2)
        [~,ind] = min(abs(tf_times{i}-limo.data.start));
        tf_times{i} = tf_times{i}(ind:end);
        [~,ind] = min(abs(tf_freqs{i}-limo.data.low_f));
        tf_freqs{i} = tf_freqs{i}(ind:end);
    end
else
    limo.data.trim1 = v;
    limo.data.start = start(c);
    if strcmp(Analysis,'Frequency')
        % loop to adjust each subject list
        for i=1:size(Paths,2)
            [~,ind] = min(abs(freqlist{i}-limo.data.start));
            freqlist{i} = freqlist{i}(ind:end);
        end
    end
end

[v,c] = min(last_frame);
if strcmp(Analysis,'Time-Frequency')
    limo.data.trim2 = v(1);
    limo.data.end = stop(c(1),1);
    limo.data.trim_high_f = v(2);    
    limo.data.high_f = stop(c(2),2);
    % loop to adjust each subject list
    for i=1:size(Paths,2)
        [~,ind] = min(abs(tf_times{i}-limo.data.end));
        tf_times{i} = tf_times{i}(1:ind);
        [~,ind] = min(abs(tf_freqs{i}-limo.data.high_f));
        tf_freqs{i} = tf_freqs{i}(1:ind);
    end
else
    limo.data.trim2 = v;
    limo.data.end = stop(c);
    if strcmp(Analysis,'Frequency')
        % loop to adjust each subject list
        for i=1:size(Paths,2)
            [~,ind] = min(abs(freqlist{i}-limo.data.end));
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
        error('the resolution of frequency lists doesn''t match between subjects')
    end
    limo.data.freqlist = mean(freqlist,1);
    limo.data.start    = limo.data.freqlist(1);
    limo.data.end      = limo.data.freqlist(end);
    
elseif strcmp(Analysis,'Time-Frequency')
    % check all lists match
    try
        tf_times = cell2mat(tf_times');
        tf_freqs = cell2mat(tf_freqs');
    catch list_issue
        error('the resolution of time/frequency lists doesn''t match between subjects')
    end
    limo.data.tf_times = mean(tf_times,1);
    limo.data.tf_freqs = mean(tf_freqs,1);
    limo.data.start    = limo.data.tf_times(1);
    limo.data.low_f    = limo.data.tf_freqs(1);
    limo.data.end      = limo.data.tf_times(end);
    limo.data.high_f   = limo.data.tf_freqs(end);
end
end
function [Names,Paths,Files] = breaklimofiles(cellfiles)
for ifiles = 1:size(cellfiles,1)
    [Paths{ifiles} filename ext] = fileparts(cellfiles{ifiles});
    Names{ifiles} = [filename ext];
    Files{ifiles} = fullfile(Paths{ifiles},[filename ext]);
end
end


