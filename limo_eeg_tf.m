function limo_eeg_tf(varargin)

% Forked from limo_eeg to run time-frequncy analyses
% This mostly runs parallel to limo_eeg{3} and +, but adapted for 4D tf
% data - the LIMO.mat is created and updated via limo_eeg.m - if data are
% 4D then limo_eeg_tf is called.
%
% IMPUT limo_eeg_tf(4) to run the GLM with/without bootstrap and tfce
%       limo_eeg_tf(5) make figures
%       limo_eeg(6) contrasts
%
% Andrew X Stewart, Nov 2013
% Cyril Pernet February 2014
% ---------------------------
% Copyright (C) LIMO Team 2014

global LIMO
if ~isnumeric(cell2mat(varargin))
    error('case number expected - see help')
end

switch varargin{1}
    
    % case{1} - see limo_eeg = GUI
    % case{2} - see limo_eeg = IMPORT
    % case{3} - see limo_eeg = DESIGN MATRIX
    
    case{4}
        
        
        % NBOOT (updated if specified in LIMO.design)
        % ------------------------------------------
        nboot = 1000;
        % ----------
        
        % get the LIMO.mat
        try
            load('LIMO.mat');
        catch
            [file,dir_path,ind] = uigetfile('LIMO.mat','select a LIMO.mat file');
            if ind ==0
                return
            else
                cd (dir_path); load LIMO.mat;
            end
        end
        cd(LIMO.dir);
        
        % ---------------- univariate analysis ------------------
        % --------------------------------------------------------
        if strcmp(LIMO.design.type_of_analysis,'Mass-univariate')
            
            % --------- load files created by limo_design_matrix ------------------
            load Yr; % load a 4D Yr with tf data,
            if sum(size(Yr) ~= LIMO.data.size4D)~=0; % then check it is so
                errordlg('Is 4D data given to limo_design_matrix_tf?','LIMO.data.size4D'); return
            end
            
            Yr = limo_tf_4d_reshape(Yr); % reshape to 3D
            if sum(size(Yr) ~= LIMO.data.size3D)~=0; % confirm now shaped 3D
                errordlg('4D data are not reshaped correctly!','LIMO.data.size3D'); return
            end
            
            Yhat = zeros(LIMO.data.size3D);
            Res = zeros(LIMO.data.size3D);
            Betas = zeros(LIMO.data.size3D(1),LIMO.data.size3D(2),size(LIMO.design.X,2));
            R2 = zeros(LIMO.data.size3D(1),LIMO.data.size3D(2),3);
            
            % ------------- prepare weight matrix -------------------------------------
            if strcmp(LIMO.design.method,'WLS') || strcmp(LIMO.design.method,'OLS')
                W = ones(LIMO.data.size3D(1),LIMO.data.size4D(2),LIMO.data.size3D(3));
            elseif strcmp(LIMO.design.method,'IRLS')
                W = zeros(size(Yr));
            end
            
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
            
            if strcmp(LIMO.design.status,'to do')
                update = 1;
                X = LIMO.design.X;
                for e = 1:size(array,1)
                    electrode = array(e); warning off;
                    fprintf('analyzing channel %g/%g \n',electrode,size(Yr,1));
                    if LIMO.Level == 2
                        Y = squeeze(Yr(electrode,:,:));
                        index = find(~isnan(Y(1,:)));
                        Y = Y(:,index);
                        LIMO.design.X = X(index,:);
                        if size(LIMO.design.X,1) <= size(LIMO.design.X,2)
                            fprintf('skipping channel %g not enough data \n',electrode);
                        else
                            model = limo_glm1(Y',LIMO); warning on;
                        end
                        if isempty(index)
                            index = [1:size(Y,2)];
                        end
                    else % level 1 we should not have any NaNs
                        index = [1:size(Yr,3)];
                        model = limo_glm1(squeeze(Yr(electrode,:,:))',LIMO);
                    end
                    
                    % update the LIMO.mat (do it only once)
                    if update == 1
                        LIMO.model.model_df = model.df;
                        if LIMO.design.nb_conditions ~=0
                            LIMO.model.conditions_df = model.conditions.df;
                        end
                        if LIMO.design.nb_interactions ~=0
                            LIMO.model.interactions_df = model.interactions.df;
                        end
                        if LIMO.design.nb_continuous ~=0
                            LIMO.model.continuous_df = model.continuous.df;
                        end
                        update = 0;
                    end
                    
                    % update the files to be stored on the disk
                    if strcmp(LIMO.design.method,'IRLS')
                        W(electrode,:,:,index) = model.W';
                    else
                        W(electrode,:,index) = model.W';
                    end
                    fitted_data = LIMO.design.X*model.betas;
                    Yhat(electrode,:,index) = fitted_data';
                    Res(electrode,:,index) = squeeze(Yr(electrode,:,index)) - fitted_data'; clear fitted_data
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
                                tmp_Interaction_effect(electrode,:,i,1) = model.interactions.F(i,:);
                                tmp_Interaction_effect(electrode,:,i,2) = model.interactions.p(i,:);
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
                LIMO.design.X = X;
                LIMO.design.weights = W;
                LIMO.design.status = 'done';
                
                
                % Save all tf data as 4D elec x freqs x times x trials
                disp('saving results to drive ..'); save LIMO LIMO;
                Yhat = limo_tf_4d_reshape(Yhat); save Yhat Yhat -v7.3;
                Res = limo_tf_4d_reshape(Res); save Res Res -v7.3;
                Betas = limo_tf_4d_reshape(Betas); save Betas Betas -v7.3;
                R2 = limo_tf_4d_reshape(R2); save R2 R2 -v7.3;
                clear Yhat Res Betas R2
                
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
                        Condition_effect = limo_tf_4d_reshape(Condition_effect);
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
                        Interaction_effect = limo_tf_4d_reshape(Interaction_effect);
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
                        Covariate_effect = limo_tf_4d_reshape(Covariate_effect);
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
                if exist('H0','dir')
                    if strcmp(questdlg('H0 directory detected, overwrite?','data check','Yes','No','No'),'No');
                        if LIMO.design.tfce == 1
                            errordlg2('bootstrap skipped - attempting to continue with tfce');
                        else
                            return
                        end
                    else
                        boot_go = 1;
                    end
                else
                    boot_go = 1;
                end
            else
                clear Yr
            end
            
            if boot_go == 1
                try
                    fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%% \n Bootstrapping data with the GLM can take a while, be patient .. \n %%%%%%%%%%%%%%%%%%%%%%%% \n')
                    mkdir H0;
                    
                    if LIMO.design.bootstrap > 599 || ~exist(nboot,'var')
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
                            [rows,cols] = size(X(index,:));
                            if (rows >= cols)
                                model = limo_glm1_boot(Y(:,index)',X(index,:),LIMO.design.nb_conditions,LIMO.design.nb_interactions,LIMO.design.nb_continuous,LIMO.design.zscore,LIMO.design.method,boot_table{electrode});
                            end
                        else
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
                                    tmp_H0_Interaction_effect(electrode,:,1,1,:) = model.interactions.F{B};
                                    tmp_H0_Interaction_effect(electrode,:,1,2,:) = model.interactions.p{B};
                                else
                                    for i=1:length(LIMO.design.nb_interactions)
                                        tmp_H0_Interaction_effect(electrode,:,i,1,:) = model.interactions.F{B}(i,:);
                                        tmp_H0_Interaction_effect(electrode,:,i,2,:) = model.interactions.p{B}(i,:);
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
                    close(h); warning on;
                    clear electrode model;
                    
                    % save data on the disk and clear out
                    disp('saving H0 data to disk ... ');
                    save(['H0' filesep 'boot_table'],'boot_table')
                    H0_Betas = limo_tf_5d_reshape(H0_Betas);
                    save(['H0' filesep 'H0_Betas'],'H0_Betas','-v7.3')
                    clear H0_Betas
                    H0_R2 = limo_tf_5d_reshape(H0_R2);
                    save(['H0' filesep 'H0_R2'], 'H0_R2', '-v7.3');
                    clear H0_R2
                    
                    if prod(LIMO.design.nb_conditions) ~=0
                        for i=1:length(LIMO.design.nb_conditions)
                            name = sprintf('H0_Condition_effect_%g',i);
                            H0_Condition_effect = squeeze(tmp_H0_Conditions(:,:,i,:,:));
                            H0_Condition_effect = limo_tf_5d_reshape(H0_Condition_effect);
                            save(['H0' filesep name],'H0_Condition_effect','-v7.3');
                            clear H0_Condition_effect
                        end
                        clear tmp_H0_Conditions
                    end
                    
                    if LIMO.design.fullfactorial == 1
                        for i=1:length(LIMO.design.nb_interactions)
                            name = sprintf('H0_Interaction_effect_%g',i);
                            H0_Interaction_effect = squeeze(tmp_H0_Interaction_effect(:,:,i,:,:));
                            H0_Interaction_effect = limo_tf_5d_reshape(H0_Interaction_effect);
                            save(['H0' filesep name],'H0_Interaction_effect','-v7.3');
                            clear H0_Interaction_effect
                        end
                        clear tmp_H0_Interaction_effect
                    end
                    
                    if LIMO.design.nb_continuous ~=0
                        for i=1:LIMO.design.nb_continuous
                            name = sprintf('H0_Covariate_effect_%g',i);
                            H0_Covariate_effect = squeeze(tmp_H0_Covariates(:,:,i,:,:));
                            H0_Covariate_effect = limo_tf_5d_reshape(H0_Covariate_effect);
                            save(['H0' filesep name],'H0_Covariate_effect','-v7.3');
                            clear H0_Covariate_effect
                        end
                        clear tmp_H0_Covariates
                    end
                    disp(' ');
                    
                catch boot_error
                    disp('an error occured while attempting to bootstrap or save the data')
                    fprintf('%s \n',boot_error.message); return
                end
                cd LIMO.dir
            end
            
            
            % TFCE if requested
            % --------------
            if LIMO.design.tfce == 1
                if isfield(LIMO.data,'neighbouring_matrix') && LIMO.design.bootstrap ~=0
                    if exist('TFCE','dir')
                        if strcmp(questdlg('TFCE directory detected, overwrite?','data check','Yes','No','No'),'No');
                            return
                        end
                    end
                    
                    fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%% \n Computing TFCE for GLM takes a while, be patient .. \n %%%%%%%%%%%%%%%%%%%%%%%% \n')
                    mkdir TFCE;
                    
                    % R2
                    load R2.mat; fprintf('Creating R2 TFCE scores \n');
                    if size(R2,1) == 1
                        tfce_score(1,:,:) = limo_tfce(2, squeeze(R2(:,:,:,2)),[]);
                    else
                        tfce_score = limo_tfce(3, squeeze(R2(:,:,:,2)),LIMO.data.neighbouring_matrix);
                    end
                    save(['H0' filesep 'tfce_R2'],'tfce_score'); clear R2;
                    
                    cd('H0'); fprintf('Thresholding H0_R2 using TFCE \n'); load H0_R2;
                    if size(H0_R2,1) == 1
                        if exist('parfor','file')
                            tfce_H0_score = NaN(1,size(H0_R2,2),size(H0_R2,3),LIMO.design.bootstrap);
                            parfor b=1:nboot
                                tfce_H0_score(1,:,:,b) = limo_tfce(2,squeeze(H0_R2(:,:,:,2,b)),LIMO.data.neighbouring_matrix,0);
                            end
                        else
                            tfce_H0_score(1,:,:,:) = limo_tfce(2, squeeze(H0_R2(:,:,:,2,:)),[]);
                        end
                    else
                        if exist('parfor','file')
                            tfce_H0_score = NaN(size(H0_R2,1),size(H0_R2,2),size(H0_R2,3),LIMO.design.bootstrap);
                            parfor b=1:nboot
                                tfce_H0_score(:,:,:,b) = limo_tfce(3,squeeze(H0_R2(:,:,:,2,b)),LIMO.data.neighbouring_matrix,0);
                            end
                        else
                            tfce_H0_score = limo_tfce(3, squeeze(H0_R2(:,:,:,2,:)),LIMO.data.neighbouring_matrix);
                        end
                    end
                    save('tfce_H0_R2','tfce_H0_score'); clear H0_R2; cd ..;
                    
                    % conditions
                    if prod(LIMO.design.nb_conditions) ~=0
                        for i=1:length(LIMO.design.nb_conditions)
                            name = sprintf('Condition_effect_%g.mat',i); load(name);
                            cd('TFCE'); fprintf('Creating Condition %g TFCE scores \n',i)
                            if size(Condition_effect,1) == 1
                                tfce_score(1,:,:) = limo_tfce(2, squeeze(Condition_effect(:,:,:,1)),[]);
                            else
                                tfce_score = limo_tfce(3, squeeze(Condition_effect(:,:,:,1)),LIMO.data.neighbouring_matrix);
                            end
                            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_score');
                            clear Condition_effect tfce_score; cd ..
                        end
                        
                        cd('H0'); fprintf('Creating H0 Condition(s) TFCE scores \n');
                        for i=1:length(LIMO.design.nb_conditions)
                            name = sprintf('H0_Condition_effect_%g.mat',i); load(name);
                            if size(H0_Condition_effect,1) == 1
                                if exist('parfor','file')
                                    tfce_H0_score = NaN(1,size(H0_Condition_effect,2),size(H0_Condition_effect,3),LIMO.design.bootstrap);
                                    parfor b=1:nboot
                                        tfce_H0_score(:,:,:,b) = limo_tfce(2,squeeze(H0_Condition_effect(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                                    end
                                else
                                    tfce_H0_score(2,:,:,:) = limo_tfce(2,squeeze(H0_Condition_effect(:,:,:,1,:)),[]);
                                end
                            else
                                if exist('parfor','file')
                                    tfce_H0_score = NaN(size(H0_Condition_effect,1),size(H0_Condition_effect,2),size(H0_Condition_effect,3),LIMO.design.bootstrap);
                                    parfor b=1:nboot
                                        tfce_H0_score(:,:,:,b) = limo_tfce(3,squeeze(H0_Condition_effect(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                                    end
                                else
                                    tfce_H0_score = limo_tfce(3,squeeze(H0_Condition_effect(:,:,:,1,:)),LIMO.data.neighbouring_matrix);
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
                                tfce_score(1,:) = limo_tfce(2,squeeze(Interaction_effect(:,:,:,1)),[]);
                            else
                                tfce_score = limo_tfce(3,squeeze(Interaction_effect(:,:,:,1)),LIMO.data.neighbouring_matrix);
                            end
                            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_score');
                            clear Interaction_effect tfce_score; cd ..
                        end
                        
                        cd('H0'); fprintf('Creating H0 Interaction(s) TFCE scores \n');
                        for i=1:length(LIMO.design.fullfactorial)
                            name = sprintf('H0_Interaction_effect_%g.mat',i); load(name);
                            if size(H0_Interaction_effect,1) == 1
                                if exist('parfor','file')
                                    tfce_H0_score = NaN(1,size(H0_Interaction_effect,2),size(H0_Interaction_effect,3),LIMO.design.bootstrap);
                                    parfor b=1:nboot
                                        tfce_H0_score(:,:,:,b) = limo_tfce(2,squeeze(H0_Interaction_effect(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                                    end
                                else
                                    tfce_H0_score(1,:,:) = limo_tfce(2,squeeze(H0_Interaction_effect(:,:,:,1,:)),[]);
                                end
                            else
                                if exist('parfor','file')
                                    tfce_H0_score = NaN(size(H0_Interaction_effect,1),size(H0_Interaction_effect,2),size(H0_Interaction_effect,3),LIMO.design.bootstrap);
                                    parfor b=1:nboot
                                        tfce_H0_score(:,:,:,b) = limo_tfce(3,squeeze(H0_Interaction_effect(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                                    end
                                else
                                    tfce_H0_score = limo_tfce(3,squeeze(H0_Interaction_effect(:,:,:,1,:)),LIMO.data.neighbouring_matrix);
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
                                tfce_score(1,:,:) = limo_tfce(2,squeeze(Covariate_effect(:,:,:,1)),[]);
                            else
                                tfce_score = limo_tfce(3,squeeze(Covariate_effect(:,:,:,1)),LIMO.data.neighbouring_matrix);
                            end
                            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_score');
                            clear Covariate_effect tfce_score; cd ..
                        end
                        
                        cd('H0'); fprintf('Creating H0 Covariate(s) TFCE scores \n');
                        for i=1:LIMO.design.nb_continuous
                            name = sprintf('H0_Covariate_effect_%g.mat',i); load(name);
                            if size(H0_Covariate_effect,1) ==1
                                if exist('parfor','file')
                                    tfce_H0_score = NaN(1,size(H0_Covariate_effect,2),size(H0_Covariate_effect,3),LIMO.design.bootstrap);
                                    parfor b=1:nboot
                                        tfce_H0_score(:,:,:,b) = limo_tfce(2,squeeze(H0_Covariate_effect(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                                    end
                                else
                                    tfce_H0_score(1,:,:,:) = limo_tfce(2,squeeze(H0_Covariate_effect(:,:,:,1,:)),[]);
                                end
                            else
                                if exist('parfor','file')
                                    tfce_H0_score = NaN(size(H0_Covariate_effect,1),size(H0_Covariate_effect,2),size(H0_Covariate_effect,3),LIMO.design.bootstrap);
                                    parfor b=1:nboot
                                        tfce_H0_score(:,:,:,b) = limo_tfce(3,squeeze(H0_Covariate_effect(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                                    end
                                else
                                    tfce_H0_score = limo_tfce(3,squeeze(H0_Covariate_effect(:,:,:,1,:)),LIMO.data.neighbouring_matrix);
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
            
            % ---------------- multivariate analysis ------------------
            % --------------------------------------------------------
        elseif strcmp(LIMO.design.type_of_analysis,'Multivariate')
            disp('multivariate in time*freq space not supported')
            % update = 1;
            %
            % % --------- load files created by limo_design_matrix ------------------
            % load Yr; load Yhat; load Res; load Betas;
            %
            % % ------------- prepare weight matrice -------------------------------------
            % if strcmp(LIMO.design.method,'WLS') || strcmp(LIMO.design.method,'OLS')
            % W = ones(size(Yr,1),size(Yr,3));
            % elseif strcmp(LIMO.design.method,'IRLS')
            % W = ones(size(Yr));
            % end
            %
            % % -------------- loop the analysis time frames per time frames
            %
            % if strcmp(LIMO.design.status,'to do')
            %
            % % 1st get weights based on time
            % if strcmp(LIMO.design.method,'WLS')
            % fprintf('getting trial weights \n')
            % array = find(~isnan(Yr(:,1,1))); % skip empty electrodes
            % for e = 1:size(Yr,1)
            % electrode = array(e); [Betas,W(e,:)] = limo_WLS(LIMO.design.X,squeeze(Yr(electrode,:,:))');
            % end
            % LIMO.design.weights = W;
            % end
            %
            % % 2nd run the multivative analysis
            % for t = 1:size(Yr,2)
            % fprintf('analyzing time frame %g/%g \n',t,size(Yr,2));
            % model = limo_glmm1(squeeze(Yr(:,t,:))',LIMO); warning off;
            %
            % % update the LIMO.mat
            % if update == 1
            % if LIMO.design.nb_conditions ~=0
            % LIMO.model.conditions_df = [model.conditions.Roy.df' model.conditions.Roy.dfe' model.conditions.Pillai.df' model.conditions.Pillai.dfe'];
            % end
            % if LIMO.design.nb_interactions ~=0
            % LIMO.model.interactions_df = [model.interactions.Roy.df' model.interactions.Roy.dfe' model.interactions.Pillai.df' model.interactions.Pillai.dfe' ];
            % end
            % if LIMO.design.nb_continuous ~=0
            % LIMO.model.continuous_df = [model.continuous.Roy.df model.continuous.Roy.dfe];
            % end
            % update = 0;
            % end
            %
            % % update the files to be stored on the disk
            % fitted_data = LIMO.design.X*model.betas;
            % Yhat(:,t,:) = fitted_data';
            % Res(:,t,:) = squeeze(Yr(:,t,:)) - fitted_data'; clear fitted_data
            % R2{t} = model.R2;
            % Betas(:,t,:) = model.betas';
            %
            % if prod(LIMO.design.nb_conditions) ~=0
            % if length(LIMO.design.nb_conditions) == 1
            % tmp_Condition_effect{t} = model.conditions;
            % else
            % for i=1:length(LIMO.design.nb_conditions)
            % tmp_Condition_effect{t}(i).EV = model.conditions.EV(i,:);
            % tmp_Condition_effect{t}(i).Roy.F = model.conditions.Roy.F(i);
            % tmp_Condition_effect{t}(i).Roy.p = model.conditions.Roy.p(i);
            % tmp_Condition_effect{t}(i).Pillai.F = model.conditions.Pillai.F(i);
            % tmp_Condition_effect{t}(i).Pillai.p = model.conditions.Pillai.p(i);
            % end
            % end
            % end
            %
            % if LIMO.design.fullfactorial == 1
            % if length(LIMO.design.nb_interactions) == 1
            % tmp_Interaction_effect{t} = model.interactions;
            % else
            % for i=1:length(LIMO.design.nb_interactions)
            % tmp_Interaction_effect{t}(i).EV = model.conditions.EV(i,:);
            % tmp_Interaction_effect{t}(i).Roy.F = model.conditions.Roy.F(i);
            % tmp_Interaction_effect{t}(i).Roy.p = model.conditions.Roy.p(i);
            % tmp_Interaction_effect{t}(i).Pillai.F = model.conditions.Pillai.F(i);
            % tmp_Interaction_effect{t}(i).Pillai.p = model.conditions.Pillai.p(i);
            % end
            % end
            % end
            %
            % if LIMO.design.nb_continuous ~=0
            % if LIMO.design.nb_continuous == 1
            % tmp_Covariate_effect{t} = model.continuous;
            % else
            % for i=1:LIMO.design.nb_continuous
            % tmp_Covariate_effect{t}(i).EV = model.conditions.EV(i,:);
            % tmp_Covariate_effect{t}(i).Roy.F = model.conditions.Roy.F(i);
            % tmp_Covariate_effect{t}(i).Roy.p = model.conditions.Roy.p(i);
            % tmp_Covariate_effect{t}(i).Pillai.F = model.conditions.Pillai.F(i);
            % tmp_Covariate_effect{t}(i).Pillai.p = model.conditions.Pillai.p(i);
            % end
            % end
            % end
            % end
            %
            % % save data on the disk and clean out
            % LIMO.design.weights = W;
            % LIMO.design.status = 'done';
            % save LIMO LIMO; save Yhat Yhat;
            % save Res Res; save Betas Betas;
            % clear Yhat Res Betas
            %
            % % R2 data
            % name = sprintf('R2_EV',i); R2_EV = NaN(size(Yr,1),size(Yr,2));
            % for t=1:size(Yr,2); R2_EV(:,t) = real(R2{t}.EV); end
            % save(name,'R2_EV','-v7.3')
            % name = sprintf('R2'); tmp = NaN(size(Yr,2),5);
            % for t=1:size(Yr,2); tmp(t,:) = [R2{t}.V R2{t}.Roy.F R2{t}.Roy.p R2{t}.Pillai.F R2{t}.Pillai.p]; end
            % R2 = tmp; save(name,'R2','-v7.3')
            %
            % % condition effects
            % if prod(LIMO.design.nb_conditions) ~=0
            % for i=1:length(LIMO.design.nb_conditions)
            % name = sprintf('Condition_effect_%g_EV',i);
            % if length(LIMO.design.nb_conditions) == 1
            % for t=1:size(Yr,2); Condition_effect_EV(:,t) = real(tmp_Condition_effect{t}.EV); end
            % save(name,'Condition_effect_EV','-v7.3')
            % name = sprintf('Condition_effect_%g',i);
            % for t=1:size(Yr,2); Condition_effect(t,:) = [tmp_Condition_effect{t}.Roy.F tmp_Condition_effect{t}.Roy.p tmp_Condition_effect{t}.Pillai.F tmp_Condition_effect{t}.Pillai.p]; end
            % save(name,'Condition_effect','-v7.3')
            % else
            % for t=1:size(Yr,2); Condition_effect_EV(:,t) = real(tmp_Condition_effect{t}(i).EV); end
            % save(name,'Condition_effect_EV','-v7.3')
            % name = sprintf('Condition_effect_%g',i);
            % for t=1:size(Yr,2); Condition_effect(t,:) = [tmp_Condition_effect{t}(i).Roy.F tmp_Condition_effect{t}(i).Roy.p tmp_Condition_effect{t}(i).Pillai.F tmp_Condition_effect{t}(i).Pillai.p]; end
            % save(name,'Condition_effect','-v7.3')
            % end
            % end
            % clear Condition_effect Condition_effect_EV tmp_Condition_effect
            % end
            %
            % % interaction effects
            % if LIMO.design.fullfactorial == 1
            % for i=1:length(LIMO.design.nb_interactions)
            % name = sprintf('Interaction_effect_%g_EV',i);
            % if length(LIMO.design.nb_interactions) == 1
            % for t=1:size(Yr,2); Interaction_effect_EV(:,t) = real(tmp_Interaction_effect{t}.EV); end
            % save(name,'Interaction_effect_EV','-v7.3')
            % name = sprintf('Interaction_effect_%g',i);
            % for t=1:size(Yr,2); Interaction_effect(t,:) = [tmp_Interaction_effect{t}.Roy.F tmp_Interaction_effect{t}.Roy.p tmp_Interaction_effect{t}.Pillai.F tmp_Interaction_effect{t}.Pillai.p]; end
            % save(name,'Interaction_effect','-v7.3')
            % else
            % for t=1:size(Yr,2); Interaction_effect_EV(:,t) = real(tmp_Interaction_effect{t}(i).EV); end
            % save(name,'Interaction_effect_EV','-v7.3')
            % name = sprintf('Interaction_effect_%g',i);
            % for t=1:size(Yr,2); Interaction_effect(t,:) = [tmp_Interaction_effect{t}(i).Roy.F tmp_Interaction_effect{t}(i).Roy.p tmp_Interaction_effect{t}(i).Pillai.F tmp_Interaction_effect{t}(i).Pillai.p]; end
            % save(name,'Interaction_effectV','-v7.3')
            % end
            % end
            % clear Interaction_effect Interaction_effect_EV tmp_Interaction_effect
            % end
            %
            % if LIMO.design.nb_continuous ~=0
            % for i=1:LIMO.design.nb_continuous
            % name = sprintf('Covariate_effect_%g_EV',i);
            % if LIMO.design.nb_continuous == 1
            % for t=1:size(Yr,2); Covariate_effect_EV(:,t) = real(tmp_Covariate_effect{t}.EV); end
            % save(name,'Covariate_effect_EV','-v7.3')
            % name = sprintf('Covariate_effect_%g',i);
            % for t=1:size(Yr,2); Covariate_effect(t,:) = [tmp_Covariate_effect{t}.Roy.F tmp_Covariate_effect{t}.Roy.p tmp_Covariate_effect{t}.Pillai.F tmp_Covariate_effect{t}.Pillai.p]; end
            % save(name,'Covariate_effect','-v7.3')
            % else
            % for t=1:size(Yr,2); Covariate_effect_EV(:,t) = real(tmp_Covariate_effect{t}(i).EV); end
            % save(name,'Covariate_effect_EV','-v7.3')
            % name = sprintf('Covariate_effect_%g',i);
            % for t=1:size(Yr,2); Covariate_effect(t,:) = [tmp_Covariate_effect{t}(i).Roy.F tmp_Covariate_effect{t}(i).Roy.p tmp_Covariate_effect{t}(i).Pillai.F tmp_Covariate_effect{t}(i).Pillai.p]; end
            % save(name,'Covariate_effect','-v7.3')
            % end
            % end
            % clear Covariate_effect Covariate_effect_EV tmp_Covariate_effect
            % end
            % clear file electrode filename model reg dir i W
            % end
            %
            %
            % % if bootsrrap
            % if LIMO.design.bootstrap == 1
            %
            % end
            %
            % % TFCE if requested
            % if LIMO.design.tfce == 1
            % end
            
        end
        warning on;
        
        
end

