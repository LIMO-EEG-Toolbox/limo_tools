function limo_glm_handling(LIMO)

% data handling function for GLM
% this function calls limo_glm, limo_glm_boot to get the analysis done, and
% organize all the files around that - externalized from limo_eeg(4)
%
% FORMAT limo_glm_handling(varargin)
%
% INPUTS
%
% OUTPUTS
%
% Cyril Pernet
% ------------------------------------------------------------------
%  Copyright (C) LIMO Team 2019

cd (LIMO.dir);

%% Compute GLM and save stats files

if strcmp(LIMO.design.status,'to do')
    Yr    = load('Yr');    Yr    = Yr.(cell2mat(fieldnames(Yr)));  
    Yhat  = load('Yhat');  Yhat  = Yhat.(cell2mat(fieldnames(Yhat)));
    Res   = load('Res');   Res   = Res.(cell2mat(fieldnames(Res)));
    R2    = load('R2');    R2    = R2.(cell2mat(fieldnames(R2)));
    Betas = load('Betas'); Betas = Betas.(cell2mat(fieldnames(Betas)));
    
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
        limo_pcout(squeeze(Yr(array(1),:,:))');
    catch pcout_error
        if strcmp(pcout_error,'Principal Component Projection cannot be computed, more observations than variables are needed')
            LIMO.design.method = 'OLS'; disp('Cannot use WLS, not enough observations - switching to OLS')
        end
    end
    
    if strcmp(LIMO.design.method,'WLS') || strcmp(LIMO.design.method,'OLS')
        W = ones(size(Yr,1),size(Yr,3));
    elseif strcmp(LIMO.design.method,'IRLS')
        W = zeros(size(Yr));
    end
    
    % ------------ run limo_glm per electrodes ---------------------------
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
            model = limo_glm(Y',LIMO); warning on;
            if isempty(index)
                index = 1:size(Y,2);
            end
        else % level 1 we should not have any NaNs
            if strcmp(LIMO.Type,'Channels')
                fprintf('analyzing channel %g/%g \n',e,size(array,1));
            else
                fprintf('analyzing component %g/%g \n',e,size(array,1));
            end
            index = 1:size(Yr,3);
            model = limo_glm(squeeze(Yr(electrode,:,:))',LIMO);
        end
        
        % update the LIMO.mat (do it only once)
        if update == 1 && ~strcmpi(LIMO.design.method,'IRLS')
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
        elseif update == 1
            LIMO.model.model_df(electrode,:) = model.df;
            if LIMO.design.nb_conditions ~=0
                LIMO.model.conditions_df(electrode,:,:)  = model.conditions.df;
            end
            if LIMO.design.nb_interactions ~=0
                LIMO.model.interactions_df(electrode,:,:)  = model.interactions.df;
            end
            if LIMO.design.nb_continuous ~=0
                LIMO.model.continuous_df(electrode,:,:)  = model.continuous.df;
            end
        end
        
        % update the files to be stored on the disk
        if strcmp(LIMO.design.method,'IRLS')
            W(electrode,:,index) = model.W';
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
            for i=1:length(LIMO.design.nb_interactions)
                tmp_Interaction_effect(electrode,:,i,1) = model.interactions.F(i,:);
                tmp_Interaction_effect(electrode,:,i,2) = model.interactions.p(i,:);
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
    LIMO.design.name   = 'GLM';
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


%% Bootstrap under H0
% -------------------------------
boot_go = 0;
if LIMO.design.bootstrap ~=0
    % avoid overwriting / recomputing H0 if done
    % (limo_eeg(4) called via the results interface)
    if ~exist('H0','dir')
        boot_go = 1;
    else
        ow = questdlg('overwrite H0?','limo check','yes','no','yes');
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
                channeighbstructmat = load([newpath filesep file]); 
                channeighbstructmat = channeighbstructmat.channeighbstructmat;
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
        mkdir H0; fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%% \n Bootstrapping GLM, ... \n %%%%%%%%%%%%%%%%%%%%%%%% \n')
        
        if ~exist('Yr','var')
            Yr = load('Yr');
            Yr = Yr.(cell2mat(fieldnames(Yr)));
        end
        
        if size(Yr,1) == 1
            array = 1;
        else
            array = find(~isnan(Yr(:,1,1))); % skip empty electrodes
        end
        
        if LIMO.design.bootstrap <= 800
            fprintf('setting bootstrap to the minimum required, i.e. 800 instead of %g\n',LIMO.design.bootstrap)
            LIMO.design.bootstrap = 800;
        end
        nboot = LIMO.design.bootstrap;
        
        if LIMO.Level == 2
            boot_table = limo_create_boot_table(Yr,nboot);
        else
            boot_table = randi(size(Yr,3),size(Yr,3),nboot);
        end
        H0_R2    = NaN(size(Yr,1), size(Yr,2), 3, nboot); % stores R, F and p values for each boot
        H0_Betas = NaN(size(Yr,1), size(Yr,2), size(LIMO.design.X,2), nboot);
        
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
                index = find(~isnan(Y(1,:))); % because across subjects, we can have missing data
                model = limo_glm_boot(Y(:,index)',X(index,:),LIMO.design.nb_conditions,LIMO.design.nb_interactions,LIMO.design.nb_continuous,LIMO.design.method,LIMO.Analysis,[],[],boot_table{electrode});
            else
                model = limo_glm_boot(squeeze(Yr(electrode,:,:))',LIMO,boot_table);
            end
            
            % update the files to be stored on the disk
            for B = 1:nboot % now loop because we use cells
                H0_Betas(electrode,:,:,B) = model.betas{B}';
                H0_R2(electrode,:,1,B) = model.R2_univariate{B};
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
                else
                    tfce_score      = limo_tfce(2, squeeze(Condition_effect(:,:,1)),LIMO.data.neighbouring_matrix);
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
                        tfce_H0_score = NaN(size(H0_Condition_effect,1),size(H0_Condition_effect,2),LIMO.design.bootstrap);
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



