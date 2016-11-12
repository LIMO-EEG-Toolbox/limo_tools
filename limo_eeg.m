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
%                 4 - call limo_glm1 (mass univariate) or limo_glm2 (multivariate) 
%                 5 - shortcut to limo_results, look at possible results and print a report
%                 6,C - shortcut to limo_contrast for the current directory,
%                 ask for a list of contrasts if not given as 2nd argument) and run them all
%                 e.g. C = [1 1 -1 -1; 1 -1 1 -1]; limo_eeg(6,C) would do those
%                 two contrasts for the data located in the current dir
%
% LIMO_EEGLIMO was primarily designed by Cyril Pernet and Guillaume Rousselet,
% with the contributon of Andrew Stewart, Nicolas Chauveau, Carl Gaspar, 
% Luisa Frei, Ignacio Suay Mas and Marianne Latinus. These authors are thereafter
% referred as the LIMO Team
%
% THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
% APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS
% AND/OR OTHER PARTIES PROVIDE THE PROGRAM “AS IS�? WITHOUT WARRANTY OF ANY KIND,
% EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
% THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.
% SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
% REPAIR OR CORRECTION.
%
% IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY
% COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE PROGRAM AS
% PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL
% OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING
% BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
% THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH
% HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
%
% Cyril Pernet & Andrew Stewart v6 21/01/2014
% Cyril Pernet & Ramon Martinez-Cancino 23-10-2014 updates for components (ICA)
%
% ------------------------------------------
% Copyright (C) LIMO Team 2016

% make sure paths are ok
local_path = which('limo_eeg');
root = fileparts(local_path);
pathCell = regexp(path, pathsep, 'split');
onPath = any(strcmp([root filesep 'help'], pathCell));
if onPath == 0
    addpath([root filesep 'limo_cluster_functions'])
    addpath([root filesep 'external'])
    addpath([root filesep 'external' filesep 'psom'])
    addpath([root filesep 'help'])
end


% in case data are already there
if isempty(varargin);
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
        disp('LIMO EEG Ref:')
        disp('Pernet, C.R., Chauveau, N., Gaspar, C., Rousselet, G.A. (2011).')
        disp('LIMO EEGLIMO: a toolbox for hierarchical LInear MOdeling of ElectroEncephaloGraphic data.')
        disp('Computational Intelligence and Neuroscience, Volume 2011')
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
                                cd(newpath);
                                channeighbstructmat = load(file);
                                channeighbstructmat = getfield(channeighbstructmat,cell2mat(fieldnames(channeighbstructmat)));
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
                    disp('the field EEGLIMO.etc.datafiles.daterp or .icaerp pointing to the data is missing - using EEGLIMO.data')
                    Y = EEGLIMO.data(:,LIMO.data.trim1:LIMO.data.trim2,:);
                end
                clear EEGLIMO
            end
            
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
                    error('the field EEGLIMO.etc.datspec pointing to the data is missing')
                end
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
                    for d=1:length(EEGLIMO.etc.datafiles.dattimef)
                        Y{d} = load('-mat',cell2mat(EEGLIMO.etc.datafiles.dattimef(d)));
                        if isstruct(Y{d}); Y{d}  = limo_struct2mat(Y{d}); end
                    end
                    Y = limo_concatcells(Y);
                    clear EEGLIMO
                elseif EEGLIMO.etc.datafiles.datersp % .mat file
                    Y = load(EEGLIMO.etc.datafiles.datersp);
                    if isstruct(Y)
                        Y = getfield(Y,cell2mat(fieldnames(Y)));
                    end
                else
                    error('no data found, the field EEGLIMO.etc.dattimef or EEGLIMO.etc.datersp pointing to the data is missing')
                end
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
                    LIMO.design.name  = sprintf('Categorical: T-test i.e. %g conditions',LIMO.design.nb_conditions);
                else
                    LIMO.design.name  = sprintf('Categorical: 1 way ANOVA with %g conditions',LIMO.design.nb_conditions);
                end
            else
                LIMO.design.name  = sprintf('Categorical: N way ANOVA with %g factors',length(LIMO.design.nb_conditions));
            end
            
        elseif prod(LIMO.design.nb_conditions) == 0 && LIMO.design.nb_continuous > 0
            if LIMO.design.nb_continuous == 1
                LIMO.design.name  = sprintf('Continuous: Simple Regression');
            else
                LIMO.design.name  = sprintf('Continuous: Multiple Regression with %g continuous variables',LIMO.design.nb_continuous);
            end
            
        elseif prod(LIMO.design.nb_conditions) > 0 && LIMO.design.nb_continuous > 0
            if length(LIMO.design.nb_conditions) == 1
                LIMO.design.name      = sprintf('AnCOVA with %g conditions and %g continuous variable(s)',LIMO.design.nb_conditions,LIMO.design.nb_continuous);
            else
                LIMO.design.name      = sprintf('AnCOVA with %g factors and %g continuous variable(s)',length(LIMO.design.nb_conditions),LIMO.design.nb_continuous);
            end
        end
        
        disp('design matrix done ...')
        
        
        % fix a bug which occurs if you run several subjects in a row with
        % the GUI and use contrasts - a new subject will have a contrast field
        % must be a way to solve this properly ??
        tofix = isfield(LIMO,'contrast');
        if tofix == 1
            LIMO.contrast = [];
        end
        
        % ---------------
        LIMO.design.status = 'to do';
        save LIMO LIMO; clear Y 
        
        a = questdlg('run the analysis?','Start GLM analysis','Yes','No','Yes');
        if strcmp(a,'Yes')
            if strcmp(LIMO.Analysis,'Time-Frequency')  
                limo_eeg_tf(4);
            else
                limo_eeg(4);
            end
            clear LIMO
            limo_gui
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
        nboot =  800;
        % ----------
        
        % get the LIMO.mat
        try
            load('LIMO.mat');
        catch
            [file,dir_newpath,ind] = uigetfile('LIMO.mat','select a LIMO.mat file');
            if ind ==0
                return
            else
                cd (dir_newpath); load LIMO.mat;
            end
        end
        cd (LIMO.dir);
        
        
        % ---------------- univariate analysis ------------------
        % --------------------------------------------------------
        if strcmp(LIMO.design.type_of_analysis,'Mass-univariate')
            
            % --------- load files created by limo_design_matrix ------------------
            if strcmp(LIMO.design.status,'to do')
                load Yr; load Yhat; load Res; load R2; load Betas;
                
                % ------------ prepare condition/covariates -------------------
                if LIMO.design.nb_conditions ~=0
                    tmp_Condition_effect = NaN(size(Yr,1),size(Yr,2),length(LIMO.design.nb_conditions),2);
                end
                
                if LIMO.design.nb_interactions ~=0
                    tmp_Interaction_effect = NaN(size(Yr,1),size(Yr,2),length(LIMO.design.nb_interactions),2);
                end
                
                if LIMO.design.nb_continuous ~=0
                    tmp_Covariate_effect = NaN(size(Yr,1),size(Yr,2),LIMO.design.nb_continuous,2);
                end
                
                % -------------- loop the analysis electrode per electrode
                if size(Yr,1) == 1
                    array = 1;
                else
                    array = find(~isnan(Yr(:,1,1))); % skip empty electrodes
                end
                
                % ------------- prepare weight matrix  -------------------------------------
                try 
                    [dist,out] = limo_pcout(squeeze(Yr(array(1),:,:))');
                catch pcout_error
                    if strcmp(pcout_error,'Principal Component Projection cannot be computed, more observations than variables are needed')
                        LIMO.design.method = 'OLS';
                        disp('Cannot use WLS, not enough observations - switching to OLS')
                    end
                end
                
                if strcmp(LIMO.design.method,'WLS') || strcmp(LIMO.design.method,'OLS')
                    W = ones(size(Yr,1),size(Yr,3));
                elseif strcmp(LIMO.design.method,'IRLS')
                    W = zeros(size(Yr));
                end

                % ------------ run limo_glm1 per electrodes ---------------------------
                update = 1;
                X = LIMO.design.X;
                for e = 1:size(array,1)
                    electrode = array(e); warning off;
                    if LIMO.Level == 2
                        fprintf('analyzing channel %g/%g \n',e,size(array,1));
                        Y = squeeze(Yr(electrode,:,:));
                        index = find(~isnan(Y(1,:)));
                        Y = Y(:,index);
                        LIMO.design.X = X(index,:);
                        model = limo_glm1(Y',LIMO); warning on;
                        if isempty(index)
                            index = [1:size(Y,2)];
                        end
                    else % level 1 we should not have any NaNs
                        if strcmp(LIMO.Type,'Channels')
                            fprintf('analyzing channel %g/%g \n',e,size(array,1));
                        else
                            fprintf('analyzing component %g/%g \n',e,size(array,1));
                        end
                        index = [1:size(Yr,3)];
                        model = limo_glm1(squeeze(Yr(electrode,:,:))',LIMO);
                    end
                    
                    % update the LIMO.mat (do it only once)
                    if update == 1
                        LIMO.model.model_df = model.df;
                        if LIMO.design.nb_conditions ~=0
                            LIMO.model.conditions_df  = model.conditions.df;
                        end
                        if LIMO.design.nb_interactions ~=0
                            LIMO.model.interactions_df  = model.interactions.df;
                        end
                        if LIMO.design.nb_continuous ~=0
                            LIMO.model.continuous_df  = model.continuous.df;
                        end
                        update = 0;
                    end
                    
                    % update the files to be stored on the disk
                    if strcmp(LIMO.design.method,'IRLS')
                        W(electrode,:,index) = model.W;
                    elseif strcmp(LIMO.design.method,'WLS')
                        W(electrode,index) = model.W;
                    end
                    fitted_data = LIMO.design.X*model.betas;
                    Yhat(electrode,:,index) = fitted_data';
                    Res(electrode,:,index)  = squeeze(Yr(electrode,:,index)) - fitted_data'; 
                    clear fitted_data
                    R2(electrode,:,1) = model.R2_univariate;
                    R2(electrode,:,2) = model.F;
                    R2(electrode,:,3) = model.p;
                    Betas(electrode,:,:,1) = model.betas';
                    
                    if prod(LIMO.design.nb_conditions) ~=0
                        if length(LIMO.design.nb_conditions) == 1
                            tmp_Condition_effect(electrode,:,1,1) = model.conditions.F;
                            tmp_Condition_effect(electrode,:,1,2) = model.conditions.p;
                        else
                            for i=1:length(LIMO.design.nb_conditions)
                                tmp_Condition_effect(electrode,:,i,1) = model.conditions.F(i,:);
                                tmp_Condition_effect(electrode,:,i,2) = model.conditions.p(i,:);
                            end
                        end
                    end
                    
                    if LIMO.design.fullfactorial == 1
                        if length(LIMO.design.nb_interactions) == 1
                            tmp_Interaction_effect(electrode,:,1,1) = model.interactions.F;
                            tmp_Interaction_effect(electrode,:,1,2) = model.interactions.p;
                        else
                            for i=1:length(LIMO.design.nb_interactions)
                                tmp_Interaction_effect(electrode,:,i,1) = model.interactions.F(:,i);
                                tmp_Interaction_effect(electrode,:,i,2) = model.interactions.p(:,i);
                            end
                        end
                    end
                    
                    if LIMO.design.nb_continuous ~=0
                        if LIMO.design.nb_continuous == 1
                            tmp_Covariate_effect(electrode,:,1,1) = model.continuous.F;
                            tmp_Covariate_effect(electrode,:,1,2) = model.continuous.p;
                        else
                            for i=1:LIMO.design.nb_continuous
                                tmp_Covariate_effect(electrode,:,i,1) = model.continuous.F(:,i);
                                tmp_Covariate_effect(electrode,:,i,2) = model.continuous.p(:,i);
                            end
                        end
                    end
                end
                
                % save data on the disk and clean out
                disp('saving data to disk')
                LIMO.design.X       = X;
                LIMO.design.weights = W;
                LIMO.design.status = 'done';
                save LIMO LIMO; save Yhat Yhat -v7.3;
                save Res Res; save Betas Betas -v7.3;
                save R2 R2 -v7.3; clear Yhat Res Betas R2
                
                if prod(LIMO.design.nb_conditions) ~=0
                    for i=1:length(LIMO.design.nb_conditions)
                        name = sprintf('Condition_effect_%g',i);
                        if size(tmp_Condition_effect,1) == 1
                            tmp = squeeze(tmp_Condition_effect(1,:,i,:));
                            Condition_effect = NaN(1,size(tmp_Condition_effect,2),2);
                            Condition_effect(1,:,:) = tmp;
                        else
                            Condition_effect = squeeze(tmp_Condition_effect(:,:,i,:));
                        end
                        save(name,'Condition_effect','-v7.3')
                    end
                    clear Condition_effect tmp_Condition_effect
                end
                
                if LIMO.design.fullfactorial == 1
                    for i=1:length(LIMO.design.nb_interactions)
                        name = sprintf('Interaction_effect_%g',i);
                        if size(tmp_Interaction_effect,1) == 1
                            tmp = squeeze(tmp_Interaction_effect(1,:,i,:));
                            Interaction_effect = NaN(1,size(tmp_Interaction_effect,2),2);
                            Interaction_effect(1,:,:) = tmp;
                        else
                            Interaction_effect = squeeze(tmp_Interaction_effect(:,:,i,:));
                        end
                        save(name,'Interaction_effect','-v7.3')
                    end
                    clear Interaction_effect tmp_Interaction_effect
                end
                
                if LIMO.design.nb_continuous ~=0
                    for i=1:LIMO.design.nb_continuous
                        name = sprintf('Covariate_effect_%g',i);
                        if size(tmp_Covariate_effect,1) == 1
                            tmp = squeeze(tmp_Covariate_effect(1,:,i,:));
                            Covariate_effect = NaN(1,size(tmp_Covariate_effect,2),2);
                            Covariate_effect(1,:,:) = tmp;
                        else
                            Covariate_effect = squeeze(tmp_Covariate_effect(:,:,i,:));
                        end
                        save(name,'Covariate_effect','-v7.3')
                    end
                    clear Covariate_effect tmp_Covariate_effect
                end
                clear file electrode filename model reg dir i W
            end
            
            
            % as above for bootstrap under H0
            % -------------------------------
            boot_go = 0;
            if LIMO.design.bootstrap ~=0
                % avoid overwriting / recomputing H0 if done
                % (limo_eeg(4) called via the results interface)
                if ~exist('H0','dir')
                    boot_go = 1;
                else
                    ow = questdlg('overwrie H0?','limo check','yes','no','yes');
                    if strcmp(ow,'yes')
                        boot_go = 1;
                    end
                end
                
                if ~isfield(LIMO.data,'neighbouring_matrix')
                    answer = questdlg('load or compute neighbouring matrix?','channel neighbouring definition','Load','Compute','Compute');
                    if strcmp(answer,'Load')
                        [file,newpath,whatsup] = uigetfile('*.mat','select neighbourghing matrix (or expected chanloc file)');
                        if whatsup == 0
                            disp('selection aborded');
                            return
                        else
                            cd(newpath); load(file); cd(LIMO.dir);
                        end
                    else
                        channeighbstructmat = limo_expected_chanlocs(LIMO.data.data, LIMO.data.data_dir);
                    end
                    LIMO.data.neighbouring_matrix = channeighbstructmat;
                    save LIMO LIMO
                end
            end
            
            if boot_go == 1
                try
                    fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%% \n Bootstrapping data with the GLM can take a while, be patient .. \n %%%%%%%%%%%%%%%%%%%%%%%% \n')
                    mkdir H0; load Yr;
                    if size(Yr,1) == 1
                        array = 1;
                    else
                        array = find(~isnan(Yr(:,1,1))); % skip empty electrodes
                    end
                
                    if LIMO.design.bootstrap > 800
                        nboot = LIMO.design.bootstrap;
                    end
                    
                    if LIMO.Level == 2
                        boot_table = limo_create_boot_table(Yr,nboot);
                    else
                        boot_table = randi(size(Yr,3),size(Yr,3),nboot);
                    end
                    H0_Betas = NaN(size(Yr,1), size(Yr,2), size(LIMO.design.X,2), nboot);
                    H0_R2 = NaN(size(Yr,1), size(Yr,2), 3, nboot); % stores R, F and p values for each boot
                    
                    if LIMO.design.nb_conditions ~= 0
                        tmp_H0_Conditions = NaN(size(Yr,1), size(Yr,2), length(LIMO.design.nb_continuous), 2, nboot);
                    end
                    
                    if LIMO.design.nb_interactions ~=0
                        tmp_H0_Interaction_effect = NaN(size(Yr,1),size(Yr,2),length(LIMO.design.nb_interactions), 2, nboot);
                    end
                    
                    if LIMO.design.nb_continuous ~= 0
                        tmp_H0_Covariates = NaN(size(Yr,1), size(Yr,2), LIMO.design.nb_continuous, 2, nboot);
                    end
                    
                    warning off;
                    X = LIMO.design.X;
                    h = waitbar(0,'bootstraping data','name','% done');
                    for e = 1:size(array,1)
                        electrode = array(e);
                        waitbar(e/size(array,1))
                        fprintf('bootstrapping electrode %g \n',electrode);
                        if LIMO.Level == 2
                            Y = squeeze(Yr(electrode,:,:));
                            index = find(~isnan(Y(1,:)));
                            model = limo_glm1_boot(Y(:,index)',X(index,:),LIMO.design.nb_conditions,LIMO.design.nb_interactions,LIMO.design.nb_continuous,LIMO.design.zscore,LIMO.design.method,LIMO.Analysis,[],[],boot_table{electrode});
                        else
                            % index = [1:size(Yr,3)];
                            model = limo_glm1_boot(squeeze(Yr(electrode,:,:))',LIMO,boot_table);
                        end
                        
                        % update the files to be stored on the disk
                        H0_Betas(electrode,:,:,:) = model.Betas;
                        
                        for B = 1:nboot % now loop because we use cells
                            H0_R2(electrode,:,1,B) = model.R2{B};
                            H0_R2(electrode,:,2,B) = model.F{B};
                            H0_R2(electrode,:,3,B) = model.p{B};
                            
                            if prod(LIMO.design.nb_conditions) ~=0
                                if length(LIMO.design.nb_conditions) == 1
                                    tmp_H0_Conditions(electrode,:,1,1,B) = model.conditions.F{B};
                                    tmp_H0_Conditions(electrode,:,1,2,B) = model.conditions.p{B};
                                else
                                    for i=1:length(LIMO.design.nb_conditions)
                                        tmp_H0_Conditions(electrode,:,i,1,B) = model.conditions.F{B}(i,:);
                                        tmp_H0_Conditions(electrode,:,i,2,B) = model.conditions.p{B}(i,:);
                                    end
                                end
                            end
                            
                            if LIMO.design.fullfactorial == 1
                                if length(LIMO.design.nb_interactions) == 1
                                    tmp_H0_Interaction_effect(electrode,:,1,1,B) = model.interactions.F{B};
                                    tmp_H0_Interaction_effect(electrode,:,1,2,B) = model.interactions.p{B};
                                else
                                    for i=1:length(LIMO.design.nb_interactions)
                                        tmp_H0_Interaction_effect(electrode,:,i,1,B) = model.interactions.F{B}(i,:);
                                        tmp_H0_Interaction_effect(electrode,:,i,2,B) = model.interactions.p{B}(i,:);
                                    end
                                end
                            end
                            
                            if LIMO.design.nb_continuous ~=0
                                if LIMO.design.nb_continuous == 1
                                    tmp_H0_Covariates(electrode,:,1,1,B) = model.continuous.F{B};
                                    tmp_H0_Covariates(electrode,:,1,2,B) = model.continuous.p{B};
                                else
                                    for i=1:LIMO.design.nb_continuous
                                        tmp_H0_Covariates(electrode,:,i,1,B) = model.continuous.F{B}(:,i);
                                        tmp_H0_Covariates(electrode,:,i,2,B) = model.continuous.p{B}(:,i);
                                    end
                                end
                            end
                        end
                    end
                    close(h)
                    warning on;
                    
                    % save data on the disk and clear out
                    cd H0
                    save boot_table boot_table
                    save H0_Betas H0_Betas -v7.3
                    save H0_R2 H0_R2 -v7.3
                    
                    
                    if prod(LIMO.design.nb_conditions) ~=0
                        for i=1:length(LIMO.design.nb_conditions)
                            name = sprintf('H0_Condition_effect_%g',i);
                            H0_Condition_effect = squeeze(tmp_H0_Conditions(:,:,i,:,:));
                            save(name,'H0_Condition_effect','-v7.3');
                            clear H0_Condition_effect
                        end
                        clear tmp_H0_Conditions
                    end
                    
                    if LIMO.design.fullfactorial == 1
                        for i=1:length(LIMO.design.nb_interactions)
                            name = sprintf('H0_Interaction_effect_%g',i);
                            H0_Interaction_effect = squeeze(tmp_H0_Interaction_effect(:,:,i,:,:));
                            save(name,'H0_Interaction_effect','-v7.3');
                            clear H0_Interaction_effect
                        end
                        clear tmp_H0_Interaction_effect
                    end
                    
                    if LIMO.design.nb_continuous ~=0
                        for i=1:LIMO.design.nb_continuous
                            name = sprintf('H0_Covariate_effect_%g',i);
                            H0_Covariate_effect = squeeze(tmp_H0_Covariates(:,:,i,:,:));
                            save(name,'H0_Covariate_effect','-v7.3');
                            clear H0_Covariate_effect
                        end
                        clear tmp_H0_Covariates
                    end
                    
                    clear electrode model H0_R2; cd ..
                    disp(' ');
                    
                catch boot_error
                    disp('an error occured while attempting to bootstrap the data')
                    fprintf('%s \n',boot_error.message); return
                end
            end
            
            
            % TFCE if requested
            % --------------
            if LIMO.design.tfce == 1
                % load Yr;
                if isfield(LIMO.data,'neighbouring_matrix') && LIMO.design.bootstrap ~=0
                    % clear Yr;
                    if exist('TFCE','dir')
                        if strcmp(questdlg('TFCE directory detected, overwrite?','data check','Yes','No','No'),'No');
                            return
                        end
                    end
                    
                    fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%% \n Computing TFCE for GLM takes a while, be patient .. \n %%%%%%%%%%%%%%%%%%%%%%%% \n')
                    mkdir TFCE; PCT_test = ver('distcomp');
                    
                    % R2
                    load R2.mat; fprintf('Creating R2 TFCE scores \n'); cd('TFCE');
                    if size(R2,1) == 1
                        tfce_score(1,:) = limo_tfce(1, squeeze(R2(:,:,2)),LIMO.data.neighbouring_matrix);
                    else
                        tfce_score = limo_tfce(2, squeeze(R2(:,:,2)),LIMO.data.neighbouring_matrix);
                    end
                    save('tfce_R2','tfce_score'); clear R2; cd ..;
                    
                    cd('H0'); fprintf('Thresholding H0_R2 using TFCE \n'); load H0_R2;
                    if size(H0_R2,1) == 1
                        if ~isempty(PCT_test)
                            tfce_H0_score = NaN(1,size(H0_R2,2),LIMO.design.bootstrap);
                            parfor b=1:nboot
                                tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_R2(:,:,2,b)),LIMO.data.neighbouring_matrix,0);
                            end
                        else
                            tfce_H0_score(1,:,:) = limo_tfce(1, squeeze(H0_R2(:,2,:)),LIMO.data.neighbouring_matrix);
                        end
                    else
                        if ~isempty(PCT_test)
                            tfce_H0_score = NaN(size(H0_R2,1),size(H0_R2,2),LIMO.design.bootstrap);
                            parfor b=1:nboot
                                tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_R2(:,:,2,b)),LIMO.data.neighbouring_matrix,0);
                            end
                        else
                            tfce_H0_score = limo_tfce(2, squeeze(H0_R2(:,:,2,:)),LIMO.data.neighbouring_matrix);
                        end
                    end
                    save('tfce_H0_R2','tfce_H0_score'); clear H0_R2; cd ..;
                    
                    % conditions
                    if prod(LIMO.design.nb_conditions) ~=0
                        for i=1:length(LIMO.design.nb_conditions)
                            name = sprintf('Condition_effect_%g.mat',i); load(name);
                            cd('TFCE'); fprintf('Creating Condition %g TFCE scores \n',i)
                            if size(Condition_effect,1) == 1
                                tfce_score(1,:) = limo_tfce(1, squeeze(Condition_effect(:,:,1)),LIMO.data.neighbouring_matrix);
                                tfce_score = limo_tfce(2, squeeze(Condition_effect(:,:,1)),LIMO.data.neighbouring_matrix);
                            end
                            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_score');
                            clear Condition_effect tfce_score; cd ..
                        end
                        
                        cd('H0'); fprintf('Creating H0 Condition(s) TFCE scores \n');
                        for i=1:length(LIMO.design.nb_conditions)
                            name = sprintf('H0_Condition_effect_%g.mat',i); load(name);
                            if size(H0_Condition_effect,1) == 1
                                if exist('parfor','file')
                                    tfce_H0_score = NaN(1,size(H0_Condition_effect,2),LIMO.design.bootstrap);
                                    parfor b=1:nboot
                                        tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Condition_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                                    end
                                else
                                    tfce_H0_score(1,:,:) = limo_tfce(1,squeeze(H0_Condition_effect(:,:,1,:)),LIMO.data.neighbouring_matrix);
                                end
                            else
                                if exist('parfor','file')
                                    tfce_H0_score = NaN(size(H0_Condition_effect,1),size(H0_R2,2),LIMO.design.bootstrap);
                                    parfor b=1:nboot
                                        tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Condition_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                                    end
                                else
                                    tfce_H0_score = limo_tfce(2,squeeze(H0_Condition_effect(:,:,1,:)),LIMO.data.neighbouring_matrix);
                                end
                            end
                            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_H0_score');
                            clear H0_Condition_effect tfce_H0_score;
                        end
                        cd ..
                    end
                    
                    % interactions
                    if LIMO.design.fullfactorial == 1
                        for i=1:length(LIMO.design.fullfactorial)
                            name = sprintf('Interaction_effect_%g.mat',i); load(name);
                            cd('TFCE'); fprintf('Creating Interaction %g TFCE scores \n',i)
                            if size(Interaction_effect,1) == 1
                                tfce_score(1,:) = limo_tfce(1,squeeze(Interaction_effect(:,:,1)),LIMO.data.neighbouring_matrix);
                            else
                                tfce_score = limo_tfce(2,squeeze(Interaction_effect(:,:,1)),LIMO.data.neighbouring_matrix);
                            end
                            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_score');
                            clear Interaction_effect tfce_score; cd ..
                        end
                        
                        cd('H0'); fprintf('Creating H0 Interaction(s) TFCE scores \n');
                        for i=1:length(LIMO.design.fullfactorial)
                            name = sprintf('H0_Interaction_effect_%g.mat',i); load(name);
                            if size(H0_Interaction_effect,1) == 1
                                if exist('parfor','file')
                                    tfce_H0_score = NaN(1,size(H0_Interaction_effect,2),LIMO.design.bootstrap);
                                    parfor b=1:nboot
                                        tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Interaction_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                                    end
                                else
                                    tfce_H0_score(1,:,:) = limo_tfce(1,squeeze(H0_Interaction_effect(:,:,1,:)),LIMO.data.neighbouring_matrix);
                                end
                            else
                                if exist('parfor','file')
                                    tfce_H0_score = NaN(size(H0_Interaction_effect,1),size(H0_Interaction_effect,2),LIMO.design.bootstrap);
                                    parfor b=1:nboot
                                        tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Interaction_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                                    end
                                else
                                    tfce_H0_score = limo_tfce(2,squeeze(H0_Interaction_effect(:,:,1,:)),LIMO.data.neighbouring_matrix);
                                end
                            end
                            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_H0_score');
                            clear H0_Interaction_effect tfce_H0_score;
                        end
                        cd ..
                    end
                    
                    % covariates / continuous regressors
                    if LIMO.design.nb_continuous ~=0
                        for i=1:LIMO.design.nb_continuous
                            name = sprintf('Covariate_effect_%g.mat',i); load(name);
                            cd('TFCE'); fprintf('Creating Covariate %g TFCE scores \n',i);
                            if size(Covariate_effect,1) == 1
                                tfce_score(1,:) = limo_tfce(1,squeeze(Covariate_effect(:,:,1)),LIMO.data.neighbouring_matrix);
                            else
                                tfce_score = limo_tfce(2,squeeze(Covariate_effect(:,:,1)),LIMO.data.neighbouring_matrix);
                            end
                            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_score');
                            clear Covariate_effect tfce_score; cd ..
                        end
                        
                        cd('H0'); fprintf('Creating H0 Covariate(s) TFCE scores \n');
                        for i=1:LIMO.design.nb_continuous
                            name = sprintf('H0_Covariate_effect_%g.mat',i); load(name);
                            if size(H0_Covariate_effect,1) == 1
                                if ~isempty(PCT_test)
                                    tfce_H0_score = NaN(1,size(H0_Covariate_effect,2),LIMO.design.bootstrap);
                                    parfor b=1:nboot
                                        tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Covariate_effect(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                                    end
                                else
                                    tfce_H0_score(1,:,:) = limo_tfce(1,squeeze(H0_Covariate_effect(:,:,1,:)),LIMO.data.neighbouring_matrix);
                                end
                            else
                                if ~isempty(PCT_test)
                                    tfce_H0_score = NaN(size(H0_Covariate_effect,1),size(H0_Covariate_effect,2),LIMO.design.bootstrap);
                                    parfor b=1:nboot
                                        tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Covariate_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                                    end
                                else
                                    tfce_H0_score = limo_tfce(2,squeeze(H0_Covariate_effect(:,:,1,:)),LIMO.data.neighbouring_matrix);
                                end
                            end
                            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_H0_score');
                            clear H0_Covariate_effect tfce_H0_score
                        end
                        cd ..
                    end
                elseif ~isfield(LIMO.data,'neighbouring_matrix')
                    disp('No TFCE performed, neighbourhood matrix missing')
                elseif  LIMO.design.bootstrap ==0
                    disp('No TFCE performed, since there was no bootstraps computed')
                end
            end
            
            
          
            % ----------------------------------------------------------
            %% ---------------- multivariate analysis ------------------
            % --------------------------------------------------------
       
        
        elseif strcmp(LIMO.design.type_of_analysis,'Multivariate')
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
            load('LIMO.mat');
        catch
            [file,dir_newpath] = uigetfile('LIMO.mat','select a LIMO.mat file');
            if file ==0
                return
            else
                cd (dir_newpath); load LIMO.mat;
            end
        end
        cd (LIMO.dir);
        
        % R2
        % ---
        if LIMO.design.bootstrap == 1
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
        
        % conditions
        if prod(LIMO.design.nb_conditions) ~=0
            for i=1:length(LIMO.design.nb_conditions)
                name = sprintf('Condition_effect_%g.mat',i);
                if LIMO.design.bootstrap == 1
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
        
        % interactions
        if LIMO.design.fullfactorial == 1
            for i=1:length(LIMO.design.nb_interactions)
                name = sprintf('Interaction_effect_%g.mat',i);
                if LIMO.design.bootstrap == 1
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
        
        % covariates / continuous regressors
        if LIMO.design.nb_continuous ~=0
            for i=1:LIMO.design.nb_continuous
                name = sprintf('Covariate_effect_%g.mat',i);
                if LIMO.design.bootstrap == 1
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
        if ~exist(LIMO,'var')
            load LIMO; cd (LIMO.dir);
        else
            [LIMO_file,LIMO_dir] = uigetfile('.mat','select a LIMO.mat file');
            cd (LIMO_dir); load LIMO.mat;
        end
        
        if ~exist(C,'var')
            [contrast_file,contrast_dir] = uigetfile({'*.mat';'*.txt'},'select your contrast file');
            cd (contrast_dir); load(contrast_file);  % problm here it has to be named C
            if strcmp(FileName(end-3:end),'.txt')
                C = importdata(contrast_file);
            elseif strcmp(FileName(end-3:end),'.mat')
                contrast_file = load(contrast_file);
                C = getfield(contrast_file,cell2mat(fieldnames(contrast_file)));
            end
            cd (LIMO.dir);
        end
        
        % Check dimensions
        C = limo_contrast_checking(LIMO.dir, LIMO.design.X, C);
        
        % Perform the analysis
        try
            previous_con = size(LIMO.contrast,2);
        catch
            previous_con = 0;
        end
        
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
            fprintf('compute contrast %g',i); disp(' ');
            % loop for each electrodes
            for electrode = 1:size(Yr,1)
                fprintf('electrode %g',electrode); disp(' ');
                result = limo_contrast(squeeze(Yr(electrode,:,:))', squeeze(Betas(electrode,:,:))', electrode, LIMO);
                
                % update multivariate results
                if LIMO.Method == 2
                    LIMO.contrast{i}.multivariate{electrode} = result;
                end
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

