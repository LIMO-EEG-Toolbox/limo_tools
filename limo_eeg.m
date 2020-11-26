function limo_eeg(varargin)

% LIMO_EEGLIMO - start up and master function of the LIMO_EEGLIMO toolbox
% Calling this function brings up different GUIs.
% Each time an option is used it calls subroutines.
% In this function is also implemented the call to the GLM, creating files
% etc .. see input
%
% LIMO_EEGLIMO is designed to perform a hierarchical LInear MOdeling of EEGLIMO data
% All analyses can be performed with this toolbox but the visualization
% relies heavily on EEGLIMOlab functions http://sccn.ucsd.edu/eeglab/
% In addition, the data format is the one used by EEGLIMOlab.
%
% INPUT limo_eeg(value,option)
%                 1 - load the GUI
%                 2,X - call limo_import (time X=1 or freuqency X=2), creating LIMO.mat file and call limo_egg(3)
%                 3 - call limo_design_matrix and populate LIMO.design
%                 4,fullfile - call limo_glm (mass univariate) or limo_glm2 (multivariate) 
%                 5 - shortcut to limo_results, look at possible results and print a report
%                 6,C - shortcut to limo_contrast for the current directory,
%                 ask for a list of contrasts if not given as 2nd argument) and run them all
%                 e.g. C = [1 1 -1 -1; 1 -1 1 -1]; limo_eeg(6,C) would do those
%                 two contrasts for the data located in the current dir
%
% Cyril Pernet & Andrew Stewart v6 21/01/2014
% Cyril Pernet & Ramon Martinez-Cancino 23-10-2014 updates for components (ICA)
%
% ------------------------------
%  Copyright (C) LIMO Team 2020

% make sure paths are ok
root = fileparts(which('limo_eeg'));
pathCell = regexp(path, pathsep, 'split');
onPath = all([sum(strcmp([root filesep 'help'],pathCell))~=0,...
    sum(strcmp([root filesep 'limo_cluster_functions'],pathCell))~=0,...
    sum(strcmp([root filesep 'external' filesep 'psom'],pathCell))~=0,...
    sum(strcmp([root filesep 'deprecated'], pathCell))~=0]);
if onPath == 0
    addpath([root filesep 'limo_cluster_functions'])
    addpath([root filesep 'external'])
    addpath([root filesep 'external' filesep 'psom'])
    addpath([root filesep 'help'])
    addpath([root filesep 'deprecated'])
end


% in case data are already there
if isempty(varargin)
    varargin={1};
end

% start
switch varargin{1}
    
    %------
    case {1}
        
        % ------------------------------------------------------------------------
        %                       GUI
        % ------------------------------------------------------------------------
        % if not called via the eeglab menu but via the matlab command window
        % show the GUI
        
        disp(' ')
        disp('LIMO_EEG was primarily designed by Cyril Pernet and Guillaume Rousselet,');
        disp('with the contributon of Andrew Stewart, Nicolas Chauveau, Carl Gaspar,');
        disp('Luisa Frei, Ignacio Suay Mas and Marianne Latinus, Ramon Martinez-Cancino,');
        disp('and Arnaud Delorme. These authors are thereafter referred as the LIMO Team');
        disp(' ')
        disp('LIMO_EEG  Copyright (C) 2015  LIMO TEAM');
        disp('This program comes with ABSOLUTELY NO WARRANTY.');
        disp('This is free software, and you are welcome to redistribute');
        disp('it under certain conditions - type help limo_eeg for details');
        disp(' ');
        disp('Please use our boilerplate Citation and Reporting:')
        disp('https://github.com/LIMO-EEG-Toolbox/limo_tools/wiki/Reporting-methods-and-results')
        disp('References are in the citations.nbid file')
        disp(' ')
        limo_gui
        
        %------
    case {2}
        
        % ------------------------------------------------------------------------
        %                       IMPORT
        % ------------------------------------------------------------------------
        % the EEGLIMO data are not imported but path / name is saved in LIMO.mat
        % Cat and Cont are imported manually from a txt or mat file
        % Other informations are i) the starting time point (sec), ii) the method to
        % use (if multivariate stats have to be computed) and iii) the working
        % directory where all informations will be saved
        
        clc; 
        if varargin{2} == 1
           out = limo_import_t;  % Data from electrodes over time in each trial
        elseif varargin{2} == 2
            out = limo_import_f;  % Data from electrodes spectral power in each trial
        elseif varargin{2} == 3
            out = limo_import_tf; % Data from electrodes spectral power over time in each trial
        end
        
        % if bootstrap with tfce - get the neighbourghing matrix now so
        % the estimation and results can be all computed without any other
        % input from user (see limo_eeg(5))
        % if bootstrap do TFCE
        if ~strcmp(out,'LIMO import aborded')
            try
                load LIMO
                if LIMO.design.bootstrap == 1
                    if ~isfield(LIMO.data,'neighbouring_matrix')
                        answer = questdlg('load or compute neighbouring matrix?','channel neighbouring definition','Load','Compute','Compute');
                        if strcmp(answer,'Load')
                            [file,newpath,whatsup] = uigetfile('*.mat','select neighbourghing matrix (or expected chanloc file)');
                            if whatsup == 0
                                disp('selection aborded');
                                return
                            else
                                channeighbstructmat = load(sprintf('%s%s',chan_path,chan_file));
                                fn = fieldnames(channeighbstructmat);
                                index = find(ismember(fn,'channeighbstructmat'));
                                channeighbstructmat = getfield(channeighbstructmat,fn{index});
                                cd(LIMO.dir);
                            end
                        else
                            channeighbstructmat = limo_expected_chanlocs(LIMO.data.data, LIMO.data.data_dir);
                        end
                        LIMO.data.neighbouring_matrix = channeighbstructmat;
                        save LIMO LIMO
                    end
                end
                disp('import done');
                
            catch
                disp('errors related to bootstrap ?? ');
                return
            end
            
            % now estimate the design matrix
            limo_eeg(3)
        else
            disp('import aborded')
        end
        
        
        %------
    case {3}
        
        % ------------------------------------------------------------------------
        %                       DESIGN MATRIX
        % ------------------------------------------------------------------------
        % returns the design matrix and some info about the matrix
        % some files are also created to be filled during the model computation
        
        % get the LIMO.mat
        try
            load LIMO
        catch
            [file,dir_newpath] = uigetfile('LIMO.mat','select a LIMO.mat file');
            if file ==0
                return
            else
                cd (dir_newpath); load LIMO.mat;
            end
        end
        cd (LIMO.dir);
        
        
        % Check data where specified and load
        if exist('EEGLIMO','var')
            if strcmp([LIMO.data.data_dir filesep LIMO.data.data],[EEGLIMO.filepath filesep EEGLIMO.filename])
                disp('Using Global variable EEGLIMO')
            else
                cd (LIMO.data.data_dir);
                disp('reloading data ..');
                EEGLIMO=pop_loadset(LIMO.data.data);
            end
        else
            cd (LIMO.data.data_dir);
            disp('reloading data ..');
            EEGLIMO=pop_loadset(LIMO.data.data);
        end
        
        % Load either elec voltage over time, elec power over frequency, or
        % electrode time-frequency -  depending on declared analysis
        
        if strcmp(LIMO.Analysis,'Time')
            if strcmp(LIMO.Type,'Components')
                if isfield(EEGLIMO.etc.datafiles,'icaerp')
                    if ~iscell(EEGLIMO.etc.datafiles.icaerp) && strcmp(EEGLIMO.etc.datafiles.icaerp(end-3:end),'.mat')
                        Y = load(EEGLIMO.etc.datafiles.icaerp);
                        if isstruct(Y)
                            Y = getfield(Y,cell2mat(fieldnames(Y)));
                        end
                    else
                        try
                            signal = load('-mat',EEGLIMO.etc.datafiles.icaerp);
                            if isstruct(signal); signal  = limo_struct2mat(signal); end
                        catch
                            for d=1:length(EEGLIMO.etc.datafiles.icaerp)
                                signal{d} = load('-mat',cell2mat(EEGLIMO.etc.datafiles.icaerp(d)));
                                if isstruct(signal{d}); signal{d}  = limo_struct2mat(signal{d}); end
                            end
                        end
                        signal = limo_concatcells(signal);
                    end
                else
                    signal = eeg_getdatact(EEGLIMO,'component',[1:size(EEGLIMO.icawinv,2)]);
                end
                Y = signal(:,LIMO.data.trim1:LIMO.data.trim2,:); clear signal

            else % channels
                if isfield(EEGLIMO.etc,'datafiles.daterp')
                    if ~iscell(EEGLIMO.etc.datafiles.daterp) && strcmp(EEGLIMO.etc.datafiles.daterp(end-3:end),'.mat')
                        Y = load(EEGLIMO.etc.datafiles.daterp);
                        if isstruct(Y)
                            Y = getfield(Y,cell2mat(fieldnames(Y)));
                        end
                    else
                        for d=1:length(EEGLIMO.etc.datafiles.daterp)
                            Y{d} = load('-mat',cell2mat(EEGLIMO.etc.datafiles.daterp(d)));
                            if isstruct(Y{d}); Y{d}  = limo_struct2mat(Y{d}); end
                        end
                        Y = limo_concatcells(Y);
                    end
                else
                    disp('the field EEG.etc.datspec pointing to the data is missing - using a hack')
                    try
                        cd(LIMO.data.data_dir); spec = dir('*.datspec');
                        for d=1:length(spec)
                            Y{d} = load('-mat',spec(d).name);
                            if isstruct(Y{d}); Y{d}  = limo_struct2mat(Y{d}); end
                        end
                        Y = limo_concatcells(Y);
                    catch
                        Y = EEGLIMO.data;
                    end
                end
                Y = Y(:,LIMO.data.trim1:LIMO.data.trim2,:);
            end
            clear EEGLIMO
            
            
        elseif strcmp(LIMO.Analysis,'Frequency')
            if strcmp(LIMO.Type,'Components')
                if isfield(EEGLIMO.etc.datafiles,'icaspec')
                    if ~iscell(EEGLIMO.etc.datafiles.icaspec) && strcmp(EEGLIMO.etc.datafiles.icaspec(end-3:end),'.mat')
                        Y = load(EEGLIMO.etc.datafiles.icaspec);
                        if isstruct(Y)
                            Y = getfield(Y,cell2mat(fieldnames(Y)));
                        end
                    else
                         try
                            signal = load('-mat',EEGLIMO.etc.datafiles.icaspec);
                            if isstruct(signal); signal  = limo_struct2mat(signal); end
                        catch
                            for d=1:length(EEGLIMO.etc.datafiles.icaspec)
                                signal{d} = load('-mat',cell2mat(EEGLIMO.etc.datafiles.icaspec(d)));
                                if isstruct(signal{d}); signal{d}  = limo_struct2mat(signal{d}); end
                            end
                        end
                        signal = limo_concatcells(signal);
                    end
                else
                    signal = eeg_getdatact(EEGLIMO,'component',[1:length(EEGLIMO.icawinv)]);
                end
                Y = signal(:,LIMO.data.trim1:LIMO.data.trim2,:); clear signal
                
            else % channels
                if isfield(EEGLIMO.etc.datafiles,'datspec')
                    if ~iscell(EEGLIMO.etc.datafiles.datspec) && strcmp(EEGLIMO.etc.datafiles.datspec(end-3:end),'.mat')
                        Y = load(EEGLIMO.etc.datafiles.datspec);
                        if isstruct(Y)
                            Y = getfield(Y,cell2mat(fieldnames(Y)));
                        end
                    else
                        for d=1:length(EEGLIMO.etc.datafiles.datspec)
                            Y{d} = load('-mat',cell2mat(EEGLIMO.etc.datafiles.datspec(d)));
                            if isstruct(Y{d}); Y{d}  = limo_struct2mat(Y{d}); end
                        end
                        Y = limo_concatcells(Y); clear EEGLIMO
                    end
                else
                    disp('the field EEG.etc.datspec pointing to the data is missing - using a hack')
                    try
                        cd(LIMO.data.data_dir); spec = dir('*.datspec');
                        for d=1:length(spec)
                            Y{d} = load('-mat',spec(d).name);
                            if isstruct(Y{d}); Y{d}  = limo_struct2mat(Y{d}); end
                        end
                        Y = limo_concatcells(Y);
                    catch
                        Y = EEGLIMO.data;
                    end
                end
                Y = Y(:,LIMO.data.trim1:LIMO.data.trim2,:);
            end
            clear EEGLIMO
            
        elseif strcmp(LIMO.Analysis,'Time-Frequency')
            disp('Time-Frequency implementation - loading tf data...');
            
            if strcmp(LIMO.Type,'Components')
                if isfield(EEGLIMO.etc.datafiles,'icatimef')
                    try
                        signal = load('-mat',EEGLIMO.etc.datafiles.icatimef);
                        if isstruct(signal); signal  = limo_struct2mat(signal); end
                    catch
                        for d=1:length(EEGLIMO.etc.datafiles.icaspec)
                            signal{d} = load('-mat',cell2mat(EEGLIMO.etc.datafiles.icatimef(d)));
                            if isstruct(signal{d}); signal{d}  = limo_struct2mat(signal{d}); end
                        end
                    end
                    signal = limo_concatcells(signal);
                else
                    signal = eeg_getdatact(EEGLIMO,'component',[1:length(EEGLIMO.icawinv)]);
                end
                Y = signal(:,LIMO.data.trim_low_f:LIMO.data.trim_high_f,LIMO.data.trim1:LIMO.data.trim2,:); clear signal

            else % channels
                if isfield(EEGLIMO.etc.datafiles,'dattimef')
                    if exist(EEGLIMO.etc.datafiles.dattimef,'file')
                        for d=1:length(EEGLIMO.etc.datafiles.dattimef)
                        Y{d} = load('-mat',cell2mat(EEGLIMO.etc.datafiles.dattimef(d)));
                        if isstruct(Y{d}); Y{d}  = limo_struct2mat(Y{d}); end
                    end
                    Y = limo_concatcells(Y);
                    clear EEGLIMO
                    else
                       error('The file %s \n in EEG.etc.datafile is not found', EEGLIMO.etc.datafiles.dattimef)
                    end
                elseif EEGLIMO.etc.datafiles.datersp % .mat file
                    if exist(EEGLIMO.etc.datafiles.datersp,'file')
                        Y = load(EEGLIMO.etc.datafiles.datersp);
                        if isstruct(Y)
                            Y = getfield(Y,cell2mat(fieldnames(Y)));
                        end
                    else
                        [~,name,ext]=fileparts(EEGLIMO.etc.datafiles.datersp);
                        try
                            Y = load([LIMO.dir filesep name ext]);
                            if isstruct(Y)
                                Y = getfield(Y,cell2mat(fieldnames(Y)));
                            end
                        catch
                            Y = load('-mat',[LIMO.dir filesep name '.datersp']);
                        end
                    end
                else
                    disp('no data found, the field EEG.etc.dattimef or EEGLIMO.etc.datersp pointing to the data is missing - using a hack')
                    try
                        cd(LIMO.data.data_dir); ersp = dir('*.dattimef');
                        for d=1:length(ersp)
                            Y{d} = load('-mat',ersp(d).name);
                            if isstruct(Y{d}); Y{d}  = limo_struct2mat(Y{d}); end
                        end
                        Y = limo_concatcells(Y);
                    catch
                        Y = EEGLIMO.data;
                    end
                end
                Y = Y(:,LIMO.data.trim_low_f:LIMO.data.trim_high_f,LIMO.data.trim1:LIMO.data.trim2,:);
            end
            LIMO.data.size4D= size(Y);
            LIMO.data.size3D= [LIMO.data.size4D(1) LIMO.data.size4D(2)*LIMO.data.size4D(3) LIMO.data.size4D(4)];
        end
        
        clear ALLCOM ALLEEGLIMO CURRENTSET CURRENTSTUDY LASTCOM STUDY 
        cd (LIMO.dir) ; save LIMO LIMO

        % make the design matrix
        disp('computing design matrix');
        if strcmp(LIMO.Analysis,'Time-Frequency') % use limo_design_matrix_tf
            [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
                LIMO.design.nb_continuous] = limo_design_matrix_tf(Y, LIMO,1);
        else  % for time or power use limo_design_matrix
            [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
                LIMO.design.nb_continuous] = limo_design_matrix(Y, LIMO,1);
        end
       
        % update LIMO.mat
        if prod(LIMO.design.nb_conditions) > 0 && LIMO.design.nb_continuous == 0
            if length(LIMO.design.nb_conditions) == 1
                if LIMO.design.nb_conditions == 2
                    LIMO.design.name  = sprintf('GLM Categorical: T-test i.e. %g conditions',LIMO.design.nb_conditions);
                else
                    LIMO.design.name  = sprintf('GLM Categorical: 1 way ANOVA with %g conditions',LIMO.design.nb_conditions);
                end
            else
                LIMO.design.name  = sprintf('GLM Categorical: N way ANOVA with %g factors',length(LIMO.design.nb_conditions));
            end
            
        elseif prod(LIMO.design.nb_conditions) == 0 && LIMO.design.nb_continuous > 0
            if LIMO.design.nb_continuous == 1
                LIMO.design.name  = sprintf('GLM Continuous: Simple Regression');
            else
                LIMO.design.name  = sprintf('GLM Continuous: Multiple Regression with %g continuous variables',LIMO.design.nb_continuous);
            end
            
        elseif prod(LIMO.design.nb_conditions) > 0 && LIMO.design.nb_continuous > 0
            if length(LIMO.design.nb_conditions) == 1
                LIMO.design.name  = sprintf('GLM AnCOVA with %g conditions and %g continuous variable(s)',LIMO.design.nb_conditions,LIMO.design.nb_continuous);
            else
                LIMO.design.name  = sprintf('GLM AnCOVA with %g factors and %g continuous variable(s)',length(LIMO.design.nb_conditions),LIMO.design.nb_continuous);
            end
        else
            LIMO.design.name = 'Mean';
        end
        
        % if you run several subjects in a row with
        % the GUI and use contrasts - a new subject can have a contrast field
        if isfield(LIMO,'contrast')
            LIMO = removefields(LIMO, 'contrast');
        end
        disp('design matrix done ...')
        
        % ---------------
        LIMO.design.status = 'to do';
        save LIMO LIMO; clear Y 
        
        a = questdlg('run the analysis?','Start GLM analysis','Yes','No','Yes');
        if strcmp(a,'Yes')
            limo_eeg(4); clear LIMO; limo_gui
        else
            return
        end
        
        %% ------------------------------------------------------------------------
        %                       ANALYZE
        % ------------------------------------------------------------------------
        % estimates the model specified in (2)
        % save all info onto disk
        
    case{4}
        
        % NBOOT (updated if specified in LIMO.design)
        % ------------------------------------------
        
        % get the LIMO.mat
        if nargin == 2
            if ischar(varargin{2})
                LIMO = load(varargin{2});
                if ~isfield(LIMO,'LIMO')
                    error('input file not recognized as a LIMO.mat structure')
                end
            else
                LIMO = varargin{2};
            end
        elseif exist(fullfile(pwd,'LIMO.mat'),'file')
            % maybe just here in the current directory
            LIMO = load('LIMO.mat'); 
        else % ask user
            [file,dir_newpath,ind] = uigetfile('LIMO.mat','select a LIMO.mat file');
            if ind ==0
                return
            else
                if strcmpi(file,'LIMO')
                    LIMO = load(fullfile(dir_newpath,'LIMO.mat'));
                else
                    error('not a LIMO.mat file')
                end
            end
        end
        
        if isfield(LIMO,'LIMO')
            LIMO = LIMO.LIMO;
        end
        
        % ---------------- univariate analysis ------------------
        % --------------------------------------------------------
        if strcmp(LIMO.design.type_of_analysis,'Mass-univariate')
            
            limo_glm_handling(LIMO)
          
            % ----------------------------------------------------------
            %% ---------------- multivariate analysis ------------------
            % --------------------------------------------------------
       
        elseif strcmp(LIMO.design.type_of_analysis,'Multivariate')
           
            % to do limo_glm_handling(LIMO)

            update = 1;
            
            % --------- load files created by limo_design_matrix ------------------
            load Yr; load Yhat; load Res; load Betas;
            
            % ------------- prepare weight matrice  -------------------------------------
            if strcmp(LIMO.design.method,'WLS') || strcmp(LIMO.design.method,'OLS')
                W = ones(size(Yr,1),size(Yr,3));
            elseif strcmp(LIMO.design.method,'IRLS')
                W = ones(size(Yr));
            end
            
            % -------------- loop the analysis time frames per time frames
            
            if strcmp(LIMO.design.status,'to do')
                
                % 1st get weights based on time
                if strcmp(LIMO.design.method,'WLS')
                    fprintf('getting trial weights \n')
                    array = find(~isnan(Yr(:,1,1))); % skip empty electrodes
                    for e = 1:size(Yr,1)
                        electrode = array(e); [Betas,W(e,:)] = limo_WLS(LIMO.design.X,squeeze(Yr(electrode,:,:))');
                    end
                    LIMO.design.weights = W;
                end
                
                % 2nd run the multivariate analysis over electrodes
                for t = 1:size(Yr,2)
                    fprintf('analysing time frame %g/%g \n',t,size(Yr,2));
                    model = limo_mglm(squeeze(Yr(:,t,:))',LIMO); warning off;
                    
                    % update the LIMO.mat
                    if update == 1
                        if LIMO.design.nb_conditions ~=0
                            LIMO.model.conditions_df  = [model.conditions.Roy.df'  model.conditions.Roy.dfe'  model.conditions.Pillai.df'  model.conditions.Pillai.dfe'];
                        end
                        if LIMO.design.nb_interactions ~=0
                            LIMO.model.interactions_df  = [model.interactions.Roy.df' model.interactions.Roy.dfe' model.interactions.Pillai.df' model.interactions.Pillai.dfe' ];
                        end
                        if LIMO.design.nb_continuous ~=0
                            LIMO.model.continuous_df  = [model.continuous.Roy.df model.continuous.Roy.dfe];
                        end
                        update = 0;
                    end
                    
                    % update the files to be stored on the disk
                    fitted_data = LIMO.design.X*model.betas;
                    Yhat(:,t,:) = fitted_data';
                    Res(:,t,:)  = squeeze(Yr(:,t,:)) - fitted_data'; clear fitted_data
                    R2{t}       = model.R2;
                    Betas(:,t,:) = model.betas';
                    
                    if prod(LIMO.design.nb_conditions) ~=0
                        if length(LIMO.design.nb_conditions) == 1
                            tmp_Condition_effect{t} = model.conditions;
                        else
                            for i=1:length(LIMO.design.nb_conditions)
                                tmp_Condition_effect{t}(i).EV = model.conditions.EV(i,:);
                                tmp_Condition_effect{t}(i).Roy.F = model.conditions.Roy.F(i);
                                tmp_Condition_effect{t}(i).Roy.p = model.conditions.Roy.p(i);
                                tmp_Condition_effect{t}(i).Pillai.F = model.conditions.Pillai.F(i);
                                tmp_Condition_effect{t}(i).Pillai.p = model.conditions.Pillai.p(i);
                            end
                        end
                    end
                    
                    if LIMO.design.fullfactorial == 1
                        if length(LIMO.design.nb_interactions) == 1
                            tmp_Interaction_effect{t} = model.interactions;
                        else
                            for i=1:length(LIMO.design.nb_interactions)
                                tmp_Interaction_effect{t}(i).EV = model.conditions.EV(i,:);
                                tmp_Interaction_effect{t}(i).Roy.F = model.conditions.Roy.F(i);
                                tmp_Interaction_effect{t}(i).Roy.p = model.conditions.Roy.p(i);
                                tmp_Interaction_effect{t}(i).Pillai.F = model.conditions.Pillai.F(i);
                                tmp_Interaction_effect{t}(i).Pillai.p = model.conditions.Pillai.p(i);
                            end
                        end
                    end
                    
                    if LIMO.design.nb_continuous ~=0
                        if LIMO.design.nb_continuous == 1
                            tmp_Covariate_effect{t} = model.continuous;
                        else
                            for i=1:LIMO.design.nb_continuous
                                tmp_Covariate_effect{t}(i).EV = model.conditions.EV(i,:);
                                tmp_Covariate_effect{t}(i).Roy.F = model.conditions.Roy.F(i);
                                tmp_Covariate_effect{t}(i).Roy.p = model.conditions.Roy.p(i);
                                tmp_Covariate_effect{t}(i).Pillai.F = model.conditions.Pillai.F(i);
                                tmp_Covariate_effect{t}(i).Pillai.p = model.conditions.Pillai.p(i);
                            end
                        end
                    end
                end
                
                % save data on the disk and clean out
                LIMO.design.weights = W;
                LIMO.design.status = 'done';
                save LIMO LIMO; save Yhat Yhat;
                save Res Res; save Betas Betas;
                clear Yhat Res Betas
                
                % R2 data
                name = sprintf('R2_EV',i); R2_EV = NaN(size(Yr,1),size(Yr,2));
                for t=1:size(Yr,2); R2_EV(:,t) = real(R2{t}.EV); end
                save(name,'R2_EV','-v7.3')
                name = sprintf('R2'); tmp = NaN(size(Yr,2),5);
                for t=1:size(Yr,2); tmp(t,:) = [R2{t}.V R2{t}.Roy.F R2{t}.Roy.p R2{t}.Pillai.F R2{t}.Pillai.p]; end
                R2 = tmp; save(name,'R2','-v7.3')
                
                % condition effects
                if prod(LIMO.design.nb_conditions) ~=0
                    for i=1:length(LIMO.design.nb_conditions)
                        name = sprintf('Condition_effect_%g_EV',i);
                        if length(LIMO.design.nb_conditions) == 1
                            for t=1:size(Yr,2); Condition_effect_EV(:,t) = real(tmp_Condition_effect{t}.EV); end
                            save(name,'Condition_effect_EV','-v7.3')
                            name = sprintf('Condition_effect_%g',i);
                            for t=1:size(Yr,2); Condition_effect(t,:) = [tmp_Condition_effect{t}.Roy.F tmp_Condition_effect{t}.Roy.p tmp_Condition_effect{t}.Pillai.F tmp_Condition_effect{t}.Pillai.p]; end
                            save(name,'Condition_effect','-v7.3')
                        else
                            for t=1:size(Yr,2); Condition_effect_EV(:,t) = real(tmp_Condition_effect{t}(i).EV); end
                            save(name,'Condition_effect_EV','-v7.3')
                            name = sprintf('Condition_effect_%g',i);
                            for t=1:size(Yr,2); Condition_effect(t,:) = [tmp_Condition_effect{t}(i).Roy.F tmp_Condition_effect{t}(i).Roy.p tmp_Condition_effect{t}(i).Pillai.F tmp_Condition_effect{t}(i).Pillai.p]; end
                            save(name,'Condition_effect','-v7.3')
                        end
                    end
                    clear Condition_effect Condition_effect_EV tmp_Condition_effect
                end
                
                % interaction effects
                if LIMO.design.fullfactorial == 1
                    for i=1:length(LIMO.design.nb_interactions)
                        name = sprintf('Interaction_effect_%g_EV',i);
                        if length(LIMO.design.nb_interactions) == 1
                            for t=1:size(Yr,2); Interaction_effect_EV(:,t) = real(tmp_Interaction_effect{t}.EV); end
                            save(name,'Interaction_effect_EV','-v7.3')
                            name = sprintf('Interaction_effect_%g',i);
                            for t=1:size(Yr,2); Interaction_effect(t,:) = [tmp_Interaction_effect{t}.Roy.F tmp_Interaction_effect{t}.Roy.p tmp_Interaction_effect{t}.Pillai.F tmp_Interaction_effect{t}.Pillai.p]; end
                            save(name,'Interaction_effect','-v7.3')
                        else
                            for t=1:size(Yr,2); Interaction_effect_EV(:,t) = real(tmp_Interaction_effect{t}(i).EV); end
                            save(name,'Interaction_effect_EV','-v7.3')
                            name = sprintf('Interaction_effect_%g',i);
                            for t=1:size(Yr,2); Interaction_effect(t,:) = [tmp_Interaction_effect{t}(i).Roy.F tmp_Interaction_effect{t}(i).Roy.p tmp_Interaction_effect{t}(i).Pillai.F tmp_Interaction_effect{t}(i).Pillai.p]; end
                            save(name,'Interaction_effectV','-v7.3')
                        end
                    end
                    clear Interaction_effect Interaction_effect_EV tmp_Interaction_effect
                end
                
                if LIMO.design.nb_continuous ~=0
                    for i=1:LIMO.design.nb_continuous
                        name = sprintf('Covariate_effect_%g_EV',i);
                        if LIMO.design.nb_continuous == 1
                            for t=1:size(Yr,2); Covariate_effect_EV(:,t) = real(tmp_Covariate_effect{t}.EV); end
                            save(name,'Covariate_effect_EV','-v7.3')
                            name = sprintf('Covariate_effect_%g',i);
                            for t=1:size(Yr,2); Covariate_effect(t,:) = [tmp_Covariate_effect{t}.Roy.F tmp_Covariate_effect{t}.Roy.p tmp_Covariate_effect{t}.Pillai.F tmp_Covariate_effect{t}.Pillai.p]; end
                            save(name,'Covariate_effect','-v7.3')
                        else
                            for t=1:size(Yr,2); Covariate_effect_EV(:,t) = real(tmp_Covariate_effect{t}(i).EV); end
                            save(name,'Covariate_effect_EV','-v7.3')
                            name = sprintf('Covariate_effect_%g',i);
                            for t=1:size(Yr,2); Covariate_effect(t,:) = [tmp_Covariate_effect{t}(i).Roy.F tmp_Covariate_effect{t}(i).Roy.p tmp_Covariate_effect{t}(i).Pillai.F tmp_Covariate_effect{t}(i).Pillai.p]; end
                            save(name,'Covariate_effect','-v7.3')
                        end
                    end
                    clear Covariate_effect Covariate_effect_EV tmp_Covariate_effect
                end
                clear file electrode filename model reg dir i W
            end
            
            
            % if bootsrrap
            if LIMO.design.bootstrap == 1
            end
            
            % TFCE if requested
            if LIMO.design.tfce == 1
            end
            
        end
        warning on;
        
    case{5}
        
        
        %% ------------------------------------------------------------------------
        %                       Results
        % ------------------------------------------------------------------------
        
        % short cut to limo_results
        % check which files are there
        % -------------------------
        try
            LIMO = load('LIMO.mat');
            LIMO = LIMO.LIMO;
        catch
            [file,dir_newpath] = uigetfile('LIMO.mat','select a LIMO.mat file');
            if file ==0
                return
            else
                cd(dir_newpath); 
                LIMO = load('LIMO.mat');
                LIMO = LIMO.LIMO;
            end
        end
        cd(LIMO.dir);
        
        % R2
        % ---
        if exist('R2.mat','file')
            if LIMO.design.bootstrap ~=0
                if LIMO.design.tfce == 1
                    limo_display_results(1,'R2.mat',pwd,0.05,3,LIMO,0);
                else
                    limo_display_results(1,'R2.mat',pwd,0.05,2,LIMO,0);
                end
            else
                limo_display_results(1,'R2.mat',pwd,0.05,1,LIMO,0);
            end
            saveas(gcf, 'R2.fig','fig'); close(gcf)
            clear R2.mat
        end
        
        % conditions
        if isfield(LIMO.design,'nb_conditions') ...
                && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true)
            if prod(LIMO.design.nb_conditions) ~=0
                for i=1:length(LIMO.design.nb_conditions)
                    name = sprintf('Condition_effect_%g.mat',i);
                    if LIMO.design.bootstrap ~=0
                        if LIMO.design.tfce == 1
                            limo_display_results(1,name,pwd,0.05,3,LIMO,0);
                        else
                            limo_display_results(1,name,pwd,0.05,2,LIMO,0);
                        end
                    else
                        limo_display_results(1,name,pwd,0.05,1,LIMO,0);
                    end
                    savename = sprintf('Condition_effect_%g.fig',i);
                    saveas(gcf, savename,'fig'); close(gcf)
                end
            end
        end
        
        % interactions
        if isfield(LIMO.design,'nb_interactions') ...
                && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true)
            if LIMO.design.fullfactorial == 1
                for i=1:length(LIMO.design.nb_interactions)
                    name = sprintf('Interaction_effect_%g.mat',i);
                    if LIMO.design.bootstrap ~=0
                        if LIMO.design.tfce == 1
                            limo_display_results(1,name,pwd,0.05,3,LIMO,0);
                        else
                            limo_display_results(1,name,pwd,0.05,2,LIMO,0);
                        end
                    else
                        limo_display_results(1,name,pwd,0.05,1,LIMO,0);
                    end
                    savename = sprintf('Interaction_effect_%g.fig',i);
                    saveas(gcf, savename,'fig'); close(gcf)
                end
            end
        end
        
        % covariates / continuous regressors
        if isfield(LIMO.design,'nb_continuous')
            if LIMO.design.nb_continuous ~=0
                for i=1:LIMO.design.nb_continuous
                    name = sprintf('Covariate_effect_%g.mat',i);
                    if LIMO.design.bootstrap ~=0
                        if LIMO.design.tfce == 1
                            limo_display_results(1,name,pwd,0.05,3,LIMO,0);
                        else
                            limo_display_results(1,name,pwd,0.05,2,LIMO,0);
                        end
                    else
                        limo_display_results(1,name,pwd,0.05,1,LIMO,0);
                    end
                    savename = sprintf('Covariate_effect_%g.fig',i);
                    saveas(gcf, savename,'fig'); close(gcf)
                end
            end
        end
         
        check_semi = dir('semi_partial_coef*.mat');
        if ~isempty(check_semi)
            for i=1:size(check_semi,2)
                name = sprintf('semi_partial_coef_%g.mat',i);
                if LIMO.design.bootstrap ~=0
                    if LIMO.design.tfce == 1
                        limo_display_results(1,name,pwd,0.05,3,LIMO,0);
                    else
                        limo_display_results(1,name,pwd,0.05,2,LIMO,0);
                    end
                else
                    limo_display_results(1,name,pwd,0.05,1,LIMO,0);
                end
                savename = sprintf('semi_partial_coef_%g.fig',i);
                saveas(gcf, savename,'fig'); close(gcf)
            end
        end
        
        other_names = {'*ttest*.mat','Rep_ANOVA*.mat','con*.mat','ess*.mat'};
        for check = 1:length(other_names)
            check_test = dir(cell2mat(other_names(check)));
            if ~isempty(check_test)
                for file = 1:size(check_test,1)
                    name = check_test(file).name;
                    if LIMO.design.bootstrap ~=0
                        if LIMO.design.tfce == 1
                            limo_display_results(1,name,pwd,0.05,3,LIMO,0);
                        else
                            limo_display_results(1,name,pwd,0.05,2,LIMO,0);
                        end
                    else
                        limo_display_results(1,name,pwd,0.05,1,LIMO,0);
                    end
                    savename = sprintf('%s.fig',name(1:end-4));
                    saveas(gcf, savename,'fig'); close(gcf)
                end
            end
        end
        
    case{6}
        
        
        %% ------------------------------------------------------------------------
        %                       Contrast
        % ------------------------------------------------------------------------
        
        
        % from the result GUI call the contrast manager; here we load a
        % series contrast -- this could be commented and put the contrast
        % right away. IMPORTANT by using limo_eeg(6) the .mat for the
        % contrast must be called C (also the name used in the contrast
        % manager. This bit is usuful for batching (replicate some part of
        % code of the contrast manager)
        
        % load LIMO and C
        if exist('LIMO.mat','file')
            load LIMO; F = getfield(LIMO,'contrast');
            for f=1:length(F)
                C(f,:) = F{f}.C;
            end
        else
            [LIMO_file,LIMO_dir] = uigetfile('.mat','select a LIMO.mat file');
            cd (LIMO_dir); load LIMO.mat;
        end
        
        previous_con = 0;
        if ~exist('C','var')
            [contrast_file,contrast_dir] = uigetfile({'*.mat';'*.txt'},'select your contrast file');
            cd (contrast_dir); load(contrast_file);  % problm here it has to be named C
            if strcmp(FileName(end-3:end),'.txt')
                C = importdata(contrast_file);
            elseif strcmp(FileName(end-3:end),'.mat')
                contrast_file = load(contrast_file);
                C = getfield(contrast_file,cell2mat(fieldnames(contrast_file)));
            end
            cd (LIMO.dir);
            if isfield(LIMO,'contrast')
                previous_con = size(LIMO.contrast,2);
            end
        end
        
        % Check dimensions
        C = limo_contrast_checking(LIMO.dir, LIMO.design.X, C);
        
        % Perform the analysis
        load Yr; load Betas;
        
        for i=1:size(C,1)  % for each contrast
            
            % check validity
            go = limo_contrast_checking(C(i,:),LIMO.design.X);
            if go == 0
                fprintf('the contrast %g is not valid',i)
                error('error line 281 in limo_eeg')
            end
            
            % update LIMO.mat
            LIMO.contrast{previous_con+i}.C = C(i,:);
            
            % create con file
            con = zeros(size(Yr,1),size(Yr,2),3); % dim 3 =F/t/p
            filename = sprintf('con_%g.mat',(i+previous_con));
            save ([filename], 'con'); clear con;
            
            % update con file
            fprintf('computing contrast %g',i); disp(' ');
            result = limo_contrast(Yr, Betas, LIMO, 0,1);
            
            % update multivariate results
            if strfind(LIMO.design.type_of_analysis,'multivariate')
                LIMO.contrast{i}.multivariate = result;
            end
            save LIMO LIMO
        end
        
        clear Yr LIMO_dir LIMO_file contrast_dir contrast_file electrode filename previous_con result C;
        
        
    case{7}
        
        % ------------------------------------------------------------------------
        %                       Gp Effects
        % ------------------------------------------------------------------------

        limo_random_effect
end

